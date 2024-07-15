! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_spectra
! unused at the moment
#if defined(_USEPOST)
  use mpi
  use cudecomp
  use mod_precision
  use mod_common_mpi, only:myid
  implicit none
  type(cudecompGridDesc) :: gd_spect_r,gd_spect_c
  integer, dimension(2) :: arrplan
  private
  public init_spectra,cmpt_spectra
  contains
  subroutine init_spectra(ng)
    !
    ! to be called under main, before calling init_wspace_arrays(),
    ! since this subroutine will fetch `wsize_fft` and increase it
    ! in case more memory is needed
    !
    use cufft
    use mod_common_cudecomp, only: handle,gd_poi
    use mod_fft            , only: wsize_fft
    use mod_utils          , only: f_sizeof
    implicit none
#if defined(_SINGLE_PRECISION)
    integer, parameter :: CUFFT_R2C_TYPE = CUFFT_R2C
    integer, parameter :: CUFFT_C2C_TYPE = CUFFT_C2C
#else
    integer, parameter :: CUFFT_R2C_TYPE = CUFFT_D2Z
    integer, parameter :: CUFFT_C2C_TYPE = CUFFT_Z2Z
#endif
    integer, intent(in ), dimension(3) :: ng
    integer :: cufft_plan_r2c_x, cufft_plan_c2c_y
    integer, dimension(3) :: n_x,nh_y
    integer :: nn,batch
    integer(int_ptr_kind()) :: wsize,max_wsize
    type(cudecompGridDescConfig) :: conf_poi,conf_spect
    type(cudecompPencilInfo) :: ap_r_x,ap_c_y
    integer :: istat
    max_wsize = 0
    !
    ! initialize grid descriptors for spectra calculation, create cuFFT plans, and
    ! determine workspace size
    !
    istat = cudecompGetGridDescConfig(handle,gd_poi,conf_poi)
    istat = cudecompGridDescConfigSetDefaults(conf_spect)
    conf_spect%transpose_comm_backend    = conf_poi%transpose_comm_backend
    conf_spect%transpose_axis_contiguous = conf_poi%transpose_axis_contiguous
    conf_spect%gdims(:) = [2*(ng(1)/2+1),ng(2),ng(3)]
    conf_spect%pdims(:) = conf_poi%pdims(1:2)
    conf_spect%gdims_dist(:) = ng(:)
    istat = cudecompGridDescCreate(handle,gd_spect_r,conf_spect)
    conf_spect%gdims(:) = [ng(1)/2+1,ng(2),ng(3)]
    conf_spect%gdims_dist(:) = 0
    istat = cudecompGridDescCreate(handle,gd_spect_c,conf_spect)
    istat = cudecompGetPencilInfo(handle,gd_spect_r,ap_r_x,1)
    istat = cudecompGetPencilInfo(handle,gd_spect_c,ap_c_y,2)
    n_x(:)  = ap_r_x%shape(:)
    nh_y(:) = ap_c_y%shape(:)
    !
    ! x-axis transforms (axis contiguous layout)
    !
    batch = product(n_x(2:3))
    nn = ng(1)
    istat = cufftCreate(cufft_plan_r2c_x)
    istat = cufftSetAutoAllocation(cufft_plan_r2c_x,0)
    istat = cufftMakePlanMany(cufft_plan_r2c_x,1,nn,null(),1,nn+2,null(),1,nn/2+1,CUFFT_R2C_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)
    !
    ! y-axis transforms (axis contiguous layout)
    !
    batch = product(nh_y(2:3))
    nn = ng(2)
    istat = cufftCreate(cufft_plan_c2c_y)
    istat = cufftSetAutoAllocation(cufft_plan_c2c_y,0)
    istat = cufftMakePlanMany(cufft_plan_c2c_y,1,nn,null(),1,nn,null(),1,nn,CUFFT_C2C_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)/f_sizeof(1._rp)
    !
    istat = cudecompGetTransposeWorkspaceSize(handle,gd_spect_r,wsize); max_wsize = max(  wsize,max_wsize)
    istat = cudecompGetTransposeWorkspaceSize(handle,gd_spect_c,wsize); max_wsize = max(2*wsize,max_wsize)
    !
    wsize_fft = max(wsize_fft,max_wsize)
    arrplan(:) = [cufft_plan_r2c_x,cufft_plan_c2c_y]
  end subroutine init_spectra
  !
  subroutine cmpt_spectra(fname,n,nh,ng,z_g,is_staggered_z,p)
    use cufft
    use cudafor
    use mod_workspaces     , only: set_cufft_wspace
    use mod_common_cudecomp, only: handle,work,solver_buf_0
    use mod_common_mpi     , only: myid,ierr
    use mod_param          , only: visc,velf,bcvel,lz
    implicit none
    integer, parameter :: nplanes = 8
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n(3),nh,ng(3)
    real(rp), intent(in), dimension(0:) :: z_g
    logical , intent(in) :: is_staggered_z
    real(rp), intent(in), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer :: istat,dtype_r,dtype_c,iunit,rlen
    integer, dimension(3) :: n_x,n_y,n_z,nh_x,nh_y,nh_z,loh_y,hih_y
    integer, dimension(nplanes) :: kplanes
    real(rp) :: ssum
    integer :: i,j,k,ii,jj,kk,jh,kplane,iim,iip,ks
    integer :: npts
    real(rp), allocatable, dimension(:,:)  , save :: spect1d_x,spect1d_x_0,spect1d_x_sum
    real(rp), allocatable, dimension(:,:)  , save :: spect1d_y,spect1d_y_0,spect1d_y_sum
    real(rp), allocatable, dimension(:,:,:), save :: spect2d
    real(rp)   , device, pointer, dimension(:,:,:) :: p_x  ,p_y  ,p_z
    complex(rp), device, pointer, dimension(:,:,:) :: p_x_c,p_y_c
    type(cudecompPencilInfo) :: ap_r_x,ap_r_y,ap_r_z,ap_c_x,ap_c_y
    real(rp) :: ub,h,reb,retau
    integer(int_ptr_kind()) :: bufsize
    !
#if defined(_SINGLE_PRECISION)
    dtype_r = CUDECOMP_FLOAT
    dtype_c = CUDECOMP_FLOAT_COMPLEX
#else
    dtype_r = CUDECOMP_DOUBLE
    dtype_c = CUDECOMP_DOUBLE_COMPLEX
#endif
    !
    ! fetch pencil configurations
    !
    istat = cudecompGetPencilInfo(handle,gd_spect_r,ap_r_x,1)
    istat = cudecompGetPencilInfo(handle,gd_spect_r,ap_r_y,2)
    istat = cudecompGetPencilInfo(handle,gd_spect_r,ap_r_z,3)
    istat = cudecompGetPencilInfo(handle,gd_spect_c,ap_c_x,1)
    istat = cudecompGetPencilInfo(handle,gd_spect_c,ap_c_y,2)
    n_x(:)  = ap_r_x%shape(:)
    n_y(:)  = ap_r_y%shape(:)
    n_z(:)  = ap_r_z%shape(:)
    nh_x(:) = ap_c_x%shape(:)
    nh_y(:) = ap_c_y%shape(:)
    loh_y(ap_c_y%order(:)) = ap_c_y%lo(:)
    hih_y(ap_c_y%order(:)) = ap_c_y%hi(:)
    !
    !$acc wait
    if(.not. allocated(spect2d)) then
      allocate(spect2d(ng(1)/2+1,ng(2)/2+1,nplanes))
      allocate(spect1d_x_0(ng(1)/2+1,ng(3)),spect1d_x_sum(ng(1)/2+1,ng(3)))
      allocate(spect1d_y_0(ng(2)/2+1,ng(3)),spect1d_y_sum(ng(2)/2+1,ng(3)))
      !$acc enter data create(spect2d,spect1d_x_0,spect1d_x_sum,spect1d_y_0,spect1d_y_sum)
      call set_cufft_wspace(arrplan)
    end if
    !$acc enter data copyin(ng,nh_y,hih_y,loh_y)
    !
    bufsize = max(ap_r_x%size,ap_r_y%size,ap_r_z%size,ap_c_x%size,ap_c_y%size)
    if(size(solver_buf_0) < bufsize) then
      !$acc exit data delete(solver_buf_0)
      deallocate(solver_buf_0)
      allocate(solver_buf_0(bufsize))
      !$acc enter data create(solver_buf_0)
    end if
    !
    ! have p_x/y/z_r and p_x/y/z_c pointing to the same memory space
    !
    !
    !$acc host_data use_device(solver_buf_0)
    call c_f_pointer(c_devloc(solver_buf_0),p_x  ,n_x )
    call c_f_pointer(c_devloc(solver_buf_0),p_y  ,n_y )
    call c_f_pointer(c_devloc(solver_buf_0),p_z  ,n_z )
    call c_f_pointer(c_devloc(solver_buf_0),p_x_c,nh_x)
    call c_f_pointer(c_devloc(solver_buf_0),p_y_c,nh_y)
    !$acc end host_data
    !
#if defined(_DECOMP_X)
    !$acc parallel loop collapse(3) default(present)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p_x(i,j,k) = p(i,j,k)
        end do
      end do
    end do
#elif defined(_DECOMP_Y)
    !$acc parallel loop collapse(3) default(present)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p_y(i,j,k) = p(i,j,k)
        end do
      end do
    end do
#elif defined(_DECOMP_Z)
    !$acc parallel loop collapse(3) default(present)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p_z(i,j,k) = p(i,j,k)
        end do
      end do
    end do
#endif
    !
    ! all these operations are done in place
    !
#if defined(_DECOMP_Z)
    !$acc host_data use_device(work)
    istat = cudecompTransposeZtoY(handle,gd_spect_r,p_z,p_y,work,dtype_r)
    istat = cudecompTransposeYtoX(handle,gd_spect_r,p_y,p_x,work,dtype_r)
    !$acc end host_data
#elif defined(_DECOMP_Y)
    !$acc host_data use_device(work)
    istat = cudecompTransposeYtoX(handle,gd_spect_r,p_y,p_x,work,dtype_r)
    !$acc end host_data
#endif
    !$acc host_data use_device(work)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecR2C(arrplan(1),p_x,p_x_c)
#else
    istat = cufftExecD2Z(arrplan(1),p_x,p_x_c)
#endif
    istat = cudecompTransposeXtoY(handle,gd_spect_c,p_x_c,p_y_c,work,dtype_c)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecC2C(arrplan(2),p_y_c,p_y_c,CUFFT_FORWARD)
#else
    istat = cufftExecZ2Z(arrplan(2),p_y_c,p_y_c,CUFFT_FORWARD)
#endif
    !$acc end host_data
    !
    ! note: local array ordering is: ap_cmplx_y%order(:) = [2,3,1];
    !
    !$acc parallel loop collapse(3) default(present)
    do k=1,nh_y(3)
      do j=1,nh_y(2)
        do i=1,nh_y(1)
          p_y_c(i,j,k) = p_y_c(i,j,k)/(1._rp*ng(1)*ng(2))
        end do
      end do
    end do
    !
    ! extract some z-planes of interest
    !
    ub    = velf(1)-0.5*(bcvel(0,3,1)+bcvel(1,3,1)) ! just in case we solve on a convective reference frame
    reb   = ub*lz/visc
    retau = 0.09_rp*reb**0.88_rp
    h     = lz/2.   ! needs to be changed for open channel
    !
    kplanes(1) = 1  ! first off-wall cell
    do k=1,ng(3)
      if ( z_g(k)/h*retau > 15 ) exit ! z^+ ~= 15
    end do
    kplanes(2) = k
    do k=1,ng(3)
     if (z_g(k).gt.0.2) exit ! z = 0.2
    end do
    kplanes(3) = k
    do k=1,ng(3)
     if (z_g(k).gt.0.5) exit ! z = 0.5
    end do
    kplanes(4) = k
    ks = 0; if(is_staggered_z) ks = 1
    do k=5,8 ! add mirror planes too
      kplanes(k) = ng(3)-kplanes(8-k+1)-ks+1
    end do
    !
    if (myid.eq.0) then
      open(10,file=fname//'_kplanes.out',form='formatted')
      do kplane=1,nplanes
        write(10,*) kplanes(kplane),z_g(kplanes(kplane))
      end do
      close(10)
    end if
    !
    ! compute 2D spectra
    !
    !$acc kernels
    spect2d(:,:,:) = 0._rp
    !$acc end kernels
    !
    ! remove mode (kx,ky) = (0,0), i.e.: the mean flow
    !
    if(loh_y(1) == 1 .and. loh_y(2) == 1) then
      !$acc parallel loop default(present)
      do k=1,hih_y(3)-loh_y(3)+1
        p_y_c(1,k,1) = (0._rp,0._rp)
      end do
    end if
    do kplane=1,nplanes
      k = kplanes(kplane)
      if(k >= loh_y(3) .and. k <= hih_y(3)) then
        !$acc parallel loop collapse(2) default(present) private(ii,jj,kk,jh)
        do i=loh_y(1),hih_y(1)
          do j=loh_y(2),hih_y(2)
            ii=j-loh_y(2)+1; jj=k-loh_y(3)+1; kk=i-loh_y(1)+1
            !
            ! n.b.: if important for performance to avoid the atomic updates below, 
            !       one can unfold the j loop into two (see jh)
            !
            jh = min(j,ng(2)-j+2)! jh = (1, 2, ..., ny/2, ny/2+1, ny/2, ..., 2)
            !$acc atomic update
            spect2d(i,jh,kplane) = spect2d(i,jh,kplane) + abs(p_y_c(ii,jj,kk))**2
            !$acc end atomic
          end do
        end do
      end if
    end do
    npts = (ng(1)/2+1)*(ng(2)/2+1)*nplanes
    !$acc update self(spect2d)
    call mpi_allreduce(MPI_IN_PLACE,spect2d(1,1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) then
      inquire(iolength=rlen) 1._rp
      open(newunit=iunit,file=fname//'_2d.bin',access='direct',recl=npts*rlen)
      write(iunit,rec=1) spect2d
      close(iunit)
    end if
    !
    ! 1D x spectra
    !
    !$acc kernels
    spect1d_x_sum(:,:) = 0._rp
    spect1d_x_0(:,:)   = 0._rp ! ky = 0 corresponds to the 1D x spectrum averaged along y
    !$acc end kernels
    !
    !$acc parallel loop collapse(2) default(present) private(ssum,ii,jj,kk)
    do i=loh_y(1),hih_y(1)
      do k=loh_y(3),hih_y(3)
        ssum = 0._rp
        jj=k-loh_y(3)+1; kk=i-loh_y(1)+1
        !$acc loop seq
        do j=loh_y(2),hih_y(2)
          ii=j-loh_y(2)+1
          ssum = ssum + abs(p_y_c(ii,jj,kk))**2
        end do
        spect1d_x_sum(i,k) = ssum
        ii=1-loh_y(2)+1
        if(loh_y(2) == 1) &
          spect1d_x_0(i,k) = abs(p_y_c(ii,jj,kk))**2
      end do
    end do
    npts = (ng(1)/2+1)*ng(3)
    !$acc update self(spect1d_x_sum,spect1d_x_0)
    call mpi_allreduce(MPI_IN_PLACE,spect1d_x_sum(1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,spect1d_x_0(  1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) then
      inquire(iolength=rlen) 1._rp
      open(newunit=iunit,file=fname//'_1d_x.bin',access='direct',recl=npts*rlen)
      write(iunit,rec=1) spect1d_x_0
      write(iunit,rec=2) spect1d_x_sum
      close(iunit)
    end if
    !
    ! 1D y spectra
    !
    !$acc kernels
    spect1d_y_sum(:,:) = 0._rp
    spect1d_y_0(:,:)   = 0._rp ! kx = 0 corresponds to the 1D y spectrum averaged along x
    !$acc end kernels
    !
    !$acc parallel loop collapse(2) default(present) private(ssum,ii,jj,kk,jh)
    do k=loh_y(3),hih_y(3)
      do j=loh_y(2),hih_y(2)
        ssum = 0._rp
        ii=j-loh_y(2)+1; jj=k-loh_y(3)+1
        !
        ! n.b.: if important for performance to avoid the atomic updates below, 
        !       one can unfold the j loop into two (see jh)
        !
        jh = min(j,ng(2)-j+2) ! jh = (1, 2, ..., ny/2, ny/2+1, ny/2, ..., 2)
        !$acc loop seq
        do i=loh_y(1),hih_y(1)
          kk=i-loh_y(1)+1
          ssum = ssum + abs(p_y_c(ii,jj,kk))**2
        end do
        !$acc atomic update
        spect1d_y_sum(jh,k) = spect1d_y_sum(jh,k) + ssum
        !$acc end atomic
        kk=1-loh_y(1)+1
        !$acc atomic update
        if(loh_y(1) == 1) &
          spect1d_y_0(jh,k) = spect1d_y_0(jh,k) + abs(p_y_c(ii,jj,kk))**2
        !$acc end atomic
      end do
    end do
    npts = (ng(2)/2+1)*ng(3)
    !$acc update self(spect1d_y_sum,spect1d_y_0)
    call mpi_allreduce(MPI_IN_PLACE,spect1d_y_sum(1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,spect1d_y_0(  1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) then
      inquire(iolength=rlen) 1._rp
      open(newunit=iunit,file=fname//'_1d_y.bin',access='direct',recl=npts*rlen)
      write(iunit,rec=1) spect1d_y_0
      write(iunit,rec=2) spect1d_y_sum
      close(iunit)
    end if
    !$acc exit data delete(ng,nh_y,hih_y,loh_y)
    !$acc wait
  end subroutine cmpt_spectra
#endif
end module mod_spectra
