! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_sgs
  use mpi
  use mod_common_mpi, only: ierr
  use mod_precision
  use mod_param, only: c_smag,big
  use mod_post, only: strain_rate
  use mod_typedef, only: cond_bound
  use mod_bound, only: boundp,bounduvw
  implicit none
  private
  public cmpt_sgs
  contains
    !
  subroutine cmpt_sgs(sgstype,n,ng,lo,hi,cbcvel,cbcpre,bcu,bcv,bcw,bcp,nb,is_bound,lwm,l,dl, &
                      zc,zf,dzc,dzf,visc,h,ind,u,v,w,dw,dw_plus,s0,uc,vc,wc,uf,vf,wf, &
                      sij,lij,mij,bcuf,bcvf,bcwf,visct)
    !
    ! compute subgrid viscosity at cell centers
    ! the current implementation of the dynamcic version is two times the 
    ! cost of the static one. Acceleration can be further achieved by
    ! call only one boundp to do the data exchange of multiple variables.
    ! Note that it does not save time by simply masking wall bc's in boundp
    ! 
    ! The dynamic version yields quite good results for Re_tau=395,550,1000,
    ! with <=5% error in the friction coeffficient.
    ! 
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    integer , intent(in ), dimension(3) :: n,ng,lo,hi
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    type(cond_bound), intent(in )                :: bcu,bcv,bcw,bcp
    type(cond_bound), intent(inout)              :: bcuf,bcvf,bcwf
    integer , intent(in ), dimension(0:1,3)      :: nb,lwm,ind
    logical , intent(in ), dimension(0:1,3)      :: is_bound
    real(rp), intent(in ), dimension(3)          :: l,dl
    real(rp), intent(in ), dimension(0:)         :: zc,zf,dzc,dzf
    real(rp), intent(in )                        :: visc,h
    real(rp), intent(in ), dimension(0:,0:,0:)   :: u,v,w,dw
    real(rp), intent(out), dimension(0:,0:,0:)   :: dw_plus,s0,uc,vc,wc,uf,vf,wf,visct
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: sij,lij,mij
    real(rp), dimension(3)        :: dli
    real(rp), dimension(0:n(3)+1) :: dzci,dzfi,alph2
    integer :: i,j,k,m
    !
    dli(:)  = dl( :)**(-1)
    dzci(:) = dzc(:)**(-1)
    dzfi(:) = dzf(:)**(-1)
    !
    select case(trim(sgstype))
    case('none')
      visct(:,:,:) = 0._rp
    case('smag')
      call strain_rate(n,dli,dzci,dzfi,is_bound,lwm,u,v,w,s0,sij)
      call cmpt_dw_plus(cbcvel,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw,dw_plus)
      call sgs_smag(n,dl,dzf,dw_plus,s0,visct)
    case('dsmag')
#if !defined(_FILTER_2D)
      alph2(:) = 4.00_rp
      alph2(1) = 2.52_rp
      alph2(n(3)) = 2.52_rp
#else
      alph2(:) = 2.52_rp
#endif
      !
      call strain_rate(n,dli,dzci,dzfi,is_bound,lwm,u,v,w,s0,sij)
      !
      visct = s0
      !
      ! Lij
      !
      call interpolate(n,u,v,w,uc,vc,wc)
      ! only periodic/patched bc's are used, since filtering is not performed
      ! in the wall-normal direction for the first off-wall layer of cells
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,uc)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,vc)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,wc)
#if !defined(_FILTER_2D)
      call filter(uc*uc,lij(:,:,:,1),is_fil2d_wall=.true.)
      call filter(vc*vc,lij(:,:,:,2),is_fil2d_wall=.true.)
      call filter(wc*wc,lij(:,:,:,3),is_fil2d_wall=.true.)
      call filter(uc*vc,lij(:,:,:,4),is_fil2d_wall=.true.)
      call filter(uc*wc,lij(:,:,:,5),is_fil2d_wall=.true.)
      call filter(vc*wc,lij(:,:,:,6),is_fil2d_wall=.true.)
      call filter(uc,uf,is_fil2d_wall=.true.)
      call filter(vc,vf,is_fil2d_wall=.true.)
      call filter(wc,wf,is_fil2d_wall=.true.)
#else
      call filter2d(uc*uc,lij(:,:,:,1))
      call filter2d(vc*vc,lij(:,:,:,2))
      call filter2d(wc*wc,lij(:,:,:,3))
      call filter2d(uc*vc,lij(:,:,:,4))
      call filter2d(uc*wc,lij(:,:,:,5))
      call filter2d(vc*wc,lij(:,:,:,6))
      call filter2d(uc,uf)
      call filter2d(vc,vf)
      call filter2d(wc,wf)
#endif
      lij(:,:,:,1) = lij(:,:,:,1) - uf*uf
      lij(:,:,:,2) = lij(:,:,:,2) - vf*vf
      lij(:,:,:,3) = lij(:,:,:,3) - wf*wf
      lij(:,:,:,4) = lij(:,:,:,4) - uf*vf
      lij(:,:,:,5) = lij(:,:,:,5) - uf*wf
      lij(:,:,:,6) = lij(:,:,:,6) - vf*wf
      !
      ! Mij
      !
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,s0)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,1))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,2))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,3))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,4))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,5))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,6))
#if !defined(_FILTER_2D)
      call filter(s0*sij(:,:,:,1),mij(:,:,:,1),is_fil2d_wall=.true.)
      call filter(s0*sij(:,:,:,2),mij(:,:,:,2),is_fil2d_wall=.true.)
      call filter(s0*sij(:,:,:,3),mij(:,:,:,3),is_fil2d_wall=.true.)
      call filter(s0*sij(:,:,:,4),mij(:,:,:,4),is_fil2d_wall=.true.)
      call filter(s0*sij(:,:,:,5),mij(:,:,:,5),is_fil2d_wall=.true.)
      call filter(s0*sij(:,:,:,6),mij(:,:,:,6),is_fil2d_wall=.true.)
      call filter(u,uf,is_fil2d_wall=.true. )
      call filter(v,vf,is_fil2d_wall=.true. )
      call filter(w,wf,is_fil2d_wall=.false.)
#else
      call filter2d(s0*sij(:,:,:,1),mij(:,:,:,1))
      call filter2d(s0*sij(:,:,:,2),mij(:,:,:,2))
      call filter2d(s0*sij(:,:,:,3),mij(:,:,:,3))
      call filter2d(s0*sij(:,:,:,4),mij(:,:,:,4))
      call filter2d(s0*sij(:,:,:,5),mij(:,:,:,5))
      call filter2d(s0*sij(:,:,:,6),mij(:,:,:,6))
      call filter2d(u,uf)
      call filter2d(v,vf)
      call filter2d(w,wf)
#endif
      ! all bc's are used here
      call bounduvw(cbcvel,n,bcuf,bcvf,bcwf,nb,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,h,ind, &
                    .true.,.false.,uf,vf,wf)
      call strain_rate(n,dli,dzci,dzfi,is_bound,lwm,uf,vf,wf,s0,sij)
      do m = 1,6
        do k = 1,n(3)
          mij(:,:,k,m) = 2._rp*(mij(:,:,k,m)-alph2(k)*s0(:,:,k)*sij(:,:,k,m))
        end do
      end do
      !
      ! cs = c_smag^2*del**2
      !
      s0(:,:,:) = mij(:,:,:,1)*lij(:,:,:,1) + &
                  mij(:,:,:,2)*lij(:,:,:,2) + &
                  mij(:,:,:,3)*lij(:,:,:,3) + &
                 (mij(:,:,:,4)*lij(:,:,:,4) + &
                  mij(:,:,:,5)*lij(:,:,:,5) + &
                  mij(:,:,:,6)*lij(:,:,:,6))*2._rp
      call ave2d(ng,lo,hi,3,l,dl,dzf,s0)
      visct = visct*s0
      s0(:,:,:) = mij(:,:,:,1)*mij(:,:,:,1) + &
                  mij(:,:,:,2)*mij(:,:,:,2) + &
                  mij(:,:,:,3)*mij(:,:,:,3) + &
                 (mij(:,:,:,4)*mij(:,:,:,4) + &
                  mij(:,:,:,5)*mij(:,:,:,5) + &
                  mij(:,:,:,6)*mij(:,:,:,6))*2._rp
      call ave2d(ng,lo,hi,3,l,dl,dzf,s0)
      visct = visct/s0
    case('amd')
      print*, 'ERROR: AMD model not yet implemented'
    case default
      print*, 'ERROR: unknown SGS model'
    end select
  end subroutine cmpt_sgs
  !
  subroutine ave2d(ng,lo,hi,idir,l,dl,dz,p)
    !
    ! average a variable over two domain directions
    !
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,l  -> uniform grid spacing and length arrays
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D scalar field
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp) :: p1d(ng(idir))
    real(dp) :: grid_area_ratio,p1d_s
    integer :: i,j,k
    !
    !$acc enter data create(p1d)
    !$acc kernels default(present)
    p1d(:) = 0._rp
    !$acc end kernels
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc parallel loop gang default(present) private(p1d_s)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*grid_area_ratio
          end do
        end do
        p1d(k) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      do k=lo(3),hi(3)
        p(:,:,k) = p1d(k)
      end do
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(j) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      do j=lo(2),hi(2)
        p(:,j,:) = p1d(j)
      end do
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(i) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      do i=lo(1),hi(1)
        p(i,:,:) = p1d(i)
      end do
    end select
  end subroutine ave2d
  !
  subroutine filter(p,pf,is_fil2d_wall)
    !
    ! top-hat filter 3D, second-order trapezoidal rule
    ! we do not define temporary variables to store array elements, since
    ! each element is used only once, which is different from mom_xyz_ad.
    ! The current implementation's cost is ~50% of that of computing eight
    ! cubes as done in Bae's code and Davidson's book.
    !
    implicit none
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    logical , intent(in ) :: is_fil2d_wall
    integer :: n(3),i,j,k
    !
    n = shape(p)-2
    pf = 0._rp
    !
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = 8._rp*(p(i,j,k)) + &
                      4._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i,j,k-1) + &
                             p(i+1,j,k) + p(i,j+1,k) + p(i,j,k+1)) + &
                      2._rp*(p(i,j+1,k+1) + p(i+1,j,k+1) + p(i+1,j+1,k) + &
                             p(i,j-1,k-1) + p(i-1,j,k-1) + p(i-1,j-1,k) + &
                             p(i,j-1,k+1) + p(i-1,j,k+1) + p(i-1,j+1,k) + &
                             p(i,j+1,k-1) + p(i+1,j,k-1) + p(i+1,j-1,k)) + &
                      1._rp*(p(i-1,j-1,k-1) + p(i+1,j-1,k-1) + p(i-1,j+1,k-1) + p(i+1,j+1,k-1) + &
                             p(i-1,j-1,k+1) + p(i+1,j-1,k+1) + p(i-1,j+1,k+1) + p(i+1,j+1,k+1))
          pf(i,j,k) = pf(i,j,k)/64._rp
        end do
      end do
    end do
    ! first off-wall layer, zc(1) and zc(n(3))
    if(is_fil2d_wall) then
      ! bottom wall
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,1) = 4._rp*(p(i,j,1)) + &
                      2._rp*(p(i-1,j,1) + p(i,j-1,1) + p(i+1,j,1) + p(i,j+1,1)) + &
                      1._rp*(p(i-1,j-1,1) + p(i+1,j-1,1) + p(i-1,j+1,1) + p(i+1,j+1,1))
          pf(i,j,1) = pf(i,j,1)/16._rp
        end do
      end do
      ! top wall
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,n(3)) = 4._rp*(p(i,j,n(3))) + &
                         2._rp*(p(i-1,j,n(3)) + p(i,j-1,n(3)) + p(i+1,j,n(3)) + p(i,j+1,n(3))) + &
                         1._rp*(p(i-1,j-1,n(3)) + p(i+1,j-1,n(3)) + p(i-1,j+1,n(3)) + p(i+1,j+1,n(3)))
          pf(i,j,n(3)) = pf(i,j,n(3))/16._rp
        end do
      end do
    end if
  end subroutine filter
  !
  subroutine filter2d(p,pf)
    !
    ! top-hat filter 2D, second-order trapezoidal rule
    ! filtering along only the homogeneous directions might be more suitable,
    ! so cs and del can be taken out of the second filtering operation,
    ! as indicated by R. Agrawal (2022)
    !
    implicit none
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer :: n(3),i,j,k
    !
    n = shape(p)-2
    pf = 0._rp
    !
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = 4._rp*(p(i,j,k)) + &
                      2._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i+1,j,k) + p(i,j+1,k)) + &
                      1._rp*(p(i-1,j-1,k) + p(i+1,j-1,k) + p(i-1,j+1,k) + p(i+1,j+1,k))
          pf(i,j,k) = pf(i,j,k)/16._rp
        end do
      end do
    end do
  end subroutine filter2d
  !
  subroutine filter_old(p,pf,is_fil2d_wall)
    !
    ! top-hat filter, second-order trapezoidal rule
    !
    implicit none
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    logical , intent(in ) :: is_fil2d_wall
    real(rp) :: p_s(8)
    integer :: n(3),i,j,k,ii,jj,kk
    !
    n = shape(p)-2
    pf = 0._rp
    !
    do kk = 1,n(3)
      do jj = 1,n(2)
        do ii = 1,n(1)
          p_s = 0._rp
          ! (1,1,1)
          do k = 0,1
            do j = 0,1
              do i = 0,1
                p_s(1) = p_s(1) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (-1,1,1)
          do k = 0,1
            do j = 0,1
              do i = -1,0
                p_s(2) = p_s(2) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (1,-1,1)
          do k = 0,1
            do j = -1,0
              do i = 0,1
                p_s(3) = p_s(3) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (1,1,-1)
          do k = -1,0
            do j = 0,1
              do i = 0,1
                p_s(4) = p_s(4) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (1,-1,-1)
          do k = -1,0
            do j = -1,0
              do i = 0,1
                p_s(5) = p_s(5) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (-1,1,-1)
          do k = -1,0
            do j = 0,1
              do i = -1,0
                p_s(6) = p_s(6) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (-1,-1,1)
          do k = 0,1
            do j = -1,0
              do i = -1,0
                p_s(7) = p_s(7) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          ! (-1,-1,-1)
          do k = -1,0
            do j = -1,0
              do i = -1,0
                p_s(8) = p_s(8) + p(ii+i,jj+j,kk+k)
              end do
            end do
          end do
          !
          pf(ii,jj,kk) = sum(p_s(1:8))/64._rp
          !
        end do
      end do
    end do
    ! quantities stored at zc(1) and zc(n(3))
    if(is_fil2d_wall) then
      ! bottom wall
      do jj = 1,n(2)
        do ii = 1,n(1)
          p_s = 0._rp
          ! (1,1)
          do j = 0,1
            do i = 0,1
              p_s(1) = p_s(1) + p(ii+i,jj+j,1)
            end do
          end do
          ! (-1,1)
          do j = 0,1
            do i = -1,0
              p_s(2) = p_s(2) + p(ii+i,jj+j,1)
            end do
          end do
          ! (1,-1)
          do j = -1,0
            do i = 0,1
              p_s(3) = p_s(3) + p(ii+i,jj+j,1)
            end do
          end do
          ! (-1,-1)
          do j = -1,0
            do i = -1,0
              p_s(4) = p_s(4) + p(ii+i,jj+j,1)
            end do
          end do
          !
          pf(ii,jj,1) = sum(p_s(1:4))/16._rp
          !
        end do
      end do
      ! top wall
      do jj = 1,n(2)
        do ii = 1,n(1)
          p_s = 0._rp
          ! (1,1)
          do j = 0,1
            do i = 0,1
              p_s(1) = p_s(1) + p(ii+i,jj+j,n(3))
            end do
          end do
          ! (-1,1)
          do j = 0,1
            do i = -1,0
              p_s(2) = p_s(2) + p(ii+i,jj+j,n(3))
            end do
          end do
          ! (1,-1)
          do j = -1,0
            do i = 0,1
              p_s(3) = p_s(3) + p(ii+i,jj+j,n(3))
            end do
          end do
          ! (-1,-1)
          do j = -1,0
            do i = -1,0
              p_s(4) = p_s(4) + p(ii+i,jj+j,n(3))
            end do
          end do
          !
          pf(ii,jj,n(3)) = sum(p_s(1:4))/16._rp
          !
        end do
      end do
    end if
  end subroutine filter_old
  !
  subroutine interpolate(n,u,v,w,uc,vc,wc)
    !
    ! interpolate velocity to cell centers,
    ! equivalent to reconstruction (FV)
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: uc,vc,wc
    integer :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc(i,j,k) = 0.5_rp*(u(i,j,k)+u(i-1,j,k))
          vc(i,j,k) = 0.5_rp*(v(i,j,k)+v(i,j-1,k))
          wc(i,j,k) = 0.5_rp*(w(i,j,k)+w(i,j,k-1))
        end do
      end do
    end do
  end subroutine
  !
  subroutine sgs_smag(n,dl,dzf,dw_plus,s0,visct)
    !
    ! classical Smagorinsky model with van Driest damping
    ! 
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dl
    real(rp), intent(in ), dimension(0:)       :: dzf
    real(rp), intent(in ), dimension(0:,0:,0:) :: dw_plus,s0
    real(rp), intent(out), dimension(0:,0:,0:) :: visct
    real(rp) :: del,fd
    integer :: i,j,k
    !
    do k=1,n(3)
      del = (dl(1)*dl(2)*dzf(k))**(1./3.)
      do j=1,n(2)
        do i=1,n(1)
          fd = 1._rp-exp(-dw_plus(i,j,k)/25._rp)
          visct(i,j,k) = (c_smag*del*fd)**2*s0(i,j,k)
        end do
      end do
    end do
  end subroutine sgs_smag
  !
  subroutine cmpt_dw_plus(cbc,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw,dw_plus)
    !
    ! inner-scaled distance to the nearest wall. We assume that a wall only
    ! affects its neighboring block, which requires that block to have enough
    ! off-wall height. Perfect partitioning has <= 2 blocks between two
    ! opposite walls. dw_plus is calculated based on minimum distance dw,
    ! instead of dw_plus, so the implementation ensures the same dw_plus
    ! under different partitionings.
    !
    ! distance to the nearest wall. Identification of walls is based on
    ! non-penetration boundary conditions, which applied to both the no-slip
    ! and free-slip (wall model) cases.
    !
    ! it is unacceptable to assume zero velocity at the wall. For no-slip walls,
    ! the velocity at the wall is zero. When a wall model is applied, tauw must 
    ! be computed using the first off-wall and ghost cells. It is incorrect to
    ! assume non-slip wall, which can lead to large errors.
    ! 
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(0:)       :: zc,dzc
    real(rp), intent(in )                      :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in ), dimension(0:,0:,0:) :: dw
    real(rp), intent(out), dimension(0:,0:,0:) :: dw_plus
    real(rp) :: tauw(2),tauw_tot,this_dw,visci
    integer :: i,j,k
    !
    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: dw1
    !
    visci = 1._rp/visc
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do i = 1,n(1)
            this_dw = dl(1)*(i-0.5)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(n(1),j,k)-v(n(1)+1,j,k)+v(n(1),j-1,k)-v(n(1)+1,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(n(1),j,k)-w(n(1)+1,j,k)+w(n(1),j,k-1)-w(n(1)+1,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do i = 1,n(1)
            this_dw = dl(1)*(n(1)-i+0.5)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do j = 1,n(2)
            this_dw = dl(2)*(j-0.5)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,n(2),k)-u(i,n(2)+1,k)+u(i-1,n(2),k)-u(i-1,n(2)+1,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,n(2),k)-w(i,n(2)+1,k)+w(i,n(2),k-1)-w(i,n(2)+1,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do j = 1,n(2)
            this_dw = dl(2)*(n(2)-j+0.5)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,1)-u(i,j,0)+u(i-1,j,1)-u(i-1,j,0))/dzc(0)
          tauw(2) = visc*0.5_rp*(v(i,j,1)-v(i,j,0)+v(i,j-1,1)-v(i,j-1,0))/dzc(0)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do k = 1,n(3)
            this_dw = zc(k)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,n(3))-u(i,j,n(3)+1)+u(i-1,j,n(3))-u(i-1,j,n(3)+1))/dzc(n(3))
          tauw(2) = visc*0.5_rp*(v(i,j,n(3))-v(i,j,n(3)+1)+v(i,j-1,n(3))-v(i,j-1,n(3)+1))/dzc(n(3))
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          do k = 1,n(3)
            this_dw = l(3)-zc(k)
            if(this_dw == dw(i,j,k)) then
              dw1(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    ! check if dw and dw1 are the same. Remove this in the future
    if (.not.all(dw(1:n(1),1:n(2),1:n(3)) == dw1(1:n(1),1:n(2),1:n(3)))) then
      print*, "dw and dw1 are different"
      call MPI_FINALIZE(ierr)
      error stop
    end if
  end subroutine cmpt_dw_plus
end module mod_sgs