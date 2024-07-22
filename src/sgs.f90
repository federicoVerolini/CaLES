! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#define _CHANNEL
module mod_sgs
  use mpi
  use mod_precision
  use mod_common_mpi, only: ierr
  use mod_param, only: c_smag,big
  use mod_typedef, only: bound
  use mod_bound, only: boundp,bounduvw
  implicit none
  private
  public cmpt_sgs
  contains
  !
  subroutine cmpt_sgs(sgstype,n,ng,lo,hi,cbcvel,cbcsgs,bcs,nb,is_bound,lwm,l,dl,dli,zc,zf,dzc,dzf, &
                      dzci,dzfi,visc,h,ind,u,v,w,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,visct)
    !
    ! compute subgrid viscosity at cell centers
    ! the LES with the dynamic model is ~2 times the cost of the LES with the static one.
    ! Acceleration can be further achieved by calling only one boundp to do the data
    ! exchange of multiple variables. Note that it does not save time to simply mask
    ! wall bc's in boundp. The dynamic version yields quite good results for Re_tau=
    ! 395,550,1000, with <=5% errors in the friction coefficient. Clipping is necessary
    ! to avoid negative values of averaged eddy viscosity (duct flow).
    !
    ! 2D filter is used in the first off-wall layer (Bae, Orlandi, LESGO and Balaras, 1995).
    ! It is difficult to do 3D filtering of Sij in the first layer, though feasible for
    ! velocity.
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    integer , intent(in ), dimension(3) :: n,ng,lo,hi
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcsgs
    type(bound), intent(in   ) :: bcs
    type(bound), intent(inout) :: bcuf,bcvf,bcwf
    type(bound), intent(in   ) :: bcu_mag,bcv_mag,bcw_mag
    integer , intent(in ), dimension(0:1,3)      :: nb,lwm,ind
    logical , intent(in ), dimension(0:1,3)      :: is_bound
    real(rp), intent(in ), dimension(3)          :: l,dl,dli
    real(rp), intent(in ), dimension(0:)         :: zc,zf,dzc,dzf,dzci,dzfi
    real(rp), intent(in )                        :: visc,h
    real(rp), intent(in ), dimension(0:,0:,0:)   :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:)   :: visct
    real(rp), allocatable, dimension(:,:,:)  , save :: dw_plus,s0,uc,vc,wc,uf,vf,wf, &
                                                       wk,wk1,wk2,wk3,alph2
    real(rp), allocatable, dimension(:,:,:,:), save :: sij,lij,mij
    logical, save :: is_first = .true.
    integer :: i,j,k,m
    !
    select case(trim(sgstype))
    case('none')
      !$acc kernels default(present) async(1)
      visct(:,:,:) = 0._rp
      !$acc end kernels
    case('smag')
      if(is_first) then
        is_first = .false.
        allocate(wk     (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk1    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk2    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk3    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 s0     (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 dw_plus(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 sij    (0:n(1)+1,0:n(2)+1,0:n(3)+1,6))
        !$acc enter data create(wk,wk1,wk2,wk3) async(1)
        !$acc enter data create(s0,dw_plus,sij) async(1)
      end if
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      call cmpt_dw_plus(cbcvel,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw_plus)
      call sgs_smag(n,dl,dzf,dw_plus,s0,visct)
    case('dsmag')
      if(is_first) then
        is_first = .false.
        allocate(wk (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk1(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk2(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk3(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 uc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 uf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 s0 (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 sij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 lij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 mij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6))
        allocate(alph2(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        !$acc enter data create(wk,wk1,wk2,wk3) async(1)
        !$acc enter data create(uc,vc,wc,uf,vf,wf) async(1)
        !$acc enter data create(s0,sij,lij,mij) async(1)
        !$acc enter data create(alph2) async(1)
        call cmpt_alph2(n,is_bound,cbcvel,alph2)
      end if
      !
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      !
      !$acc kernels default(present) async(1)
      visct = s0
      !$acc end kernels
      !
      ! Lij
      !
      call interpolate(n,u,v,w,uc,vc,wc)
      ! periodic/patched bc's are handled, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,uc)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,vc)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,wc)
#if !defined(_FILTER_2D)
      !$acc kernels default(present) async(1)
      lij(:,:,:,1) = uc*uc
      lij(:,:,:,2) = vc*vc
      lij(:,:,:,3) = wc*wc
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,lij(:,:,:,1),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,lij(:,:,:,2),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,lij(:,:,:,3),wk3,iface=0,cbc=cbcvel)
      call filter3d(n,wk1,lij(:,:,:,1))
      call filter3d(n,wk2,lij(:,:,:,2))
      call filter3d(n,wk3,lij(:,:,:,3))
      !$acc kernels default(present) async(1)
      lij(:,:,:,4) = uc*vc
      lij(:,:,:,5) = uc*wc
      lij(:,:,:,6) = vc*wc
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,lij(:,:,:,4),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,lij(:,:,:,5),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,lij(:,:,:,6),wk3,iface=0,cbc=cbcvel)
      call filter3d(n,wk1,lij(:,:,:,4))
      call filter3d(n,wk2,lij(:,:,:,5))
      call filter3d(n,wk3,lij(:,:,:,6))
      call extrapolate(n,is_bound,dzci,uc,wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,vc,wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wc,wk3,iface=0,cbc=cbcvel)
      call filter3d(n,wk1,uf)
      call filter3d(n,wk2,vf)
      call filter3d(n,wk3,wf)
#else
      !$acc kernels default(present) async(1)
      wk1 = uc*uc
      wk2 = vc*vc
      wk3 = wc*wc
      !$acc end kernels
      call filter2d(n,wk1,lij(:,:,:,1))
      call filter2d(n,wk2,lij(:,:,:,2))
      call filter2d(n,wk3,lij(:,:,:,3))
      !$acc kernels default(present) async(1)
      wk1 = uc*vc
      wk2 = uc*wc
      wk3 = vc*wc
      !$acc end kernels
      call filter2d(n,wk1,lij(:,:,:,4))
      call filter2d(n,wk2,lij(:,:,:,5))
      call filter2d(n,wk3,lij(:,:,:,6))
      call filter2d(n,uc,uf)
      call filter2d(n,vc,vf)
      call filter2d(n,wc,wf)
#endif
      !$acc kernels default(present) async(1)
      lij(:,:,:,1) = lij(:,:,:,1) - uf*uf
      lij(:,:,:,2) = lij(:,:,:,2) - vf*vf
      lij(:,:,:,3) = lij(:,:,:,3) - wf*wf
      lij(:,:,:,4) = lij(:,:,:,4) - uf*vf
      lij(:,:,:,5) = lij(:,:,:,5) - uf*wf
      lij(:,:,:,6) = lij(:,:,:,6) - vf*wf
      !$acc end kernels
      !
      ! Mij
      !
      ! periodic/patched bc's are handled, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,s0)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,1))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,2))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,3))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,4))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,5))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,6))
#if !defined(_FILTER_2D)
      !$acc kernels default(present) async(1)
      mij(:,:,:,1) = s0*sij(:,:,:,1)
      mij(:,:,:,2) = s0*sij(:,:,:,2)
      mij(:,:,:,3) = s0*sij(:,:,:,3)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,mij(:,:,:,1),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,mij(:,:,:,2),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,mij(:,:,:,3),wk3,iface=0,cbc=cbcvel)
      call filter3d(n,wk1,mij(:,:,:,1))
      call filter3d(n,wk2,mij(:,:,:,2))
      call filter3d(n,wk3,mij(:,:,:,3))
      !$acc kernels default(present) async(1)
      mij(:,:,:,4) = s0*sij(:,:,:,4)
      mij(:,:,:,5) = s0*sij(:,:,:,5)
      mij(:,:,:,6) = s0*sij(:,:,:,6)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,mij(:,:,:,4),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,mij(:,:,:,5),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,mij(:,:,:,6),wk3,iface=0,cbc=cbcvel)
      call filter3d(n,wk1,mij(:,:,:,4))
      call filter3d(n,wk2,mij(:,:,:,5))
      call filter3d(n,wk3,mij(:,:,:,6))
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,cbc=cbcvel)
      call filter3d(n,wk1,uf)
      call filter3d(n,wk2,vf)
      call filter3d(n,wk3,wf)
#else
      !$acc kernels default(present) async(1)
      wk1 = s0*sij(:,:,:,1)
      wk2 = s0*sij(:,:,:,2)
      wk3 = s0*sij(:,:,:,3)
      !$acc end kernels
      call filter2d(n,wk1,mij(:,:,:,1))
      call filter2d(n,wk2,mij(:,:,:,2))
      call filter2d(n,wk3,mij(:,:,:,3))
      !$acc kernels default(present) async(1)
      wk1 = s0*sij(:,:,:,4)
      wk2 = s0*sij(:,:,:,5)
      wk3 = s0*sij(:,:,:,6)
      !$acc end kernels
      call filter2d(n,wk1,mij(:,:,:,4))
      call filter2d(n,wk2,mij(:,:,:,5))
      call filter2d(n,wk3,mij(:,:,:,6))
      call filter2d(n,u,uf)
      call filter2d(n,v,vf)
      call filter2d(n,w,wf)
#endif
      ! when not using a wall model, all bc's are used for computing strain rate, and
      ! bcv/v/wf are equal to bcu/v/w=0 for no-slip walls. When using a wall model,
      ! wall bc's for wall-parallel velocity are not used due to extrapolation, but
      ! the wall bc's for the wall-normal velicity are used to impose no-penetration
      ! on the filtered velocity. Also, periodic/patched bc's are used. Hence,
      ! there is no need to update wall model bc's in bounduvw. Filtered velocity
      ! should satisfy the wall bc's, evidenced by the fact that each mode satisfies
      ! the no-slip/no-penetration bc's if a spectral filter is applied.
      call bounduvw(cbcvel,n,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,nb,is_bound,lwm, &
                    l,dl,zc,zf,dzc,dzf,visc,h,ind,.false.,.false.,uf,vf,wf)
      call extrapolate(n,is_bound,dzci,uf,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,vf,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wf,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      !$acc kernels default(present) async(1)
      do m = 1,6
        mij(:,:,:,m) = 2._rp*(mij(:,:,:,m)-alph2*s0*sij(:,:,:,m))
      end do
      !$acc end kernels
      !
      ! cs = c_smag^2*del**2
      !
      !$acc kernels default(present) async(1)
      s0(:,:,:) = mij(:,:,:,1)*lij(:,:,:,1) + &
                  mij(:,:,:,2)*lij(:,:,:,2) + &
                  mij(:,:,:,3)*lij(:,:,:,3) + &
                 (mij(:,:,:,4)*lij(:,:,:,4) + &
                  mij(:,:,:,5)*lij(:,:,:,5) + &
                  mij(:,:,:,6)*lij(:,:,:,6))*2._rp
      !$acc end kernels
#if defined(_CHANNEL)
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,s0)
#elif defined(_DUCT)
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,s0)
#elif defined(_CAVITY)
      !
#endif
      !$acc kernels default(present) async(1)
      visct = visct*s0
      !$acc end kernels
      !
      !$acc kernels default(present) async(1)
      s0(:,:,:) = mij(:,:,:,1)*mij(:,:,:,1) + &
                  mij(:,:,:,2)*mij(:,:,:,2) + &
                  mij(:,:,:,3)*mij(:,:,:,3) + &
                 (mij(:,:,:,4)*mij(:,:,:,4) + &
                  mij(:,:,:,5)*mij(:,:,:,5) + &
                  mij(:,:,:,6)*mij(:,:,:,6))*2._rp
      !$acc end kernels
#if defined(_CHANNEL)
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,s0)
#elif defined(_DUCT)
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,s0)
#elif defined(_CAVITY)
      !
#endif
      !$acc kernels default(present) async(1)
      visct = max(visct/s0,0._rp)
      !$acc end kernels
    case('amd')
      print*, 'ERROR: AMD model not yet implemented'
    case default
      print*, 'ERROR: unknown SGS model'
    end select
  end subroutine cmpt_sgs
  !
  subroutine ave1d_channel(ng,lo,hi,idir,l,dl,dz,p)
    !
    ! average a variable over two domain directions
    ! adapted from out1d, channel
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
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc data copyout(p1d) async(1)
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
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
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do k=lo(3),hi(3)
        p(:,:,k) = p1d(k)
      end do
      !$acc end kernels
      !$acc end data
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
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
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do j=lo(2),hi(2)
        p(:,j,:) = p1d(j)
      end do
      !$acc end kernels
      !$acc end data
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
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
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do i=lo(1),hi(1)
        p(i,:,:) = p1d(i)
      end do
      !$acc end kernels
      !$acc end data
    end select
  end subroutine ave1d_channel
  !
  subroutine ave2d_duct(ng,lo,hi,idir,l,dl,dz,p)
    !
    ! average a variable over one domain direction,
    ! adapted from out2d_duct
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp), allocatable, dimension(:,:) :: p2d
    real(dp) :: grid_area_ratio,p2d_s
    integer :: i,j,k
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      allocate(p2d(ng(1),ng(3)))
      !$acc data copyout(p2d) async(1)
      !$acc parallel loop gang collapse(2) default(present) private(p2d_s) async(1)
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          p2d_s = 0._rp
          !$acc loop reduction(+:p2d_s)
          do j=lo(2),hi(2)
            p2d_s = p2d_s + p(i,j,k)
          end do
          p2d(i,k) = p2d_s
        end do
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(1)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p2d) async(1)
      !$acc kernels default(present) async(1)
      do j=lo(2),hi(2)
        p(lo(1):hi(1),j,lo(3):hi(3)) = p2d(lo(1):hi(1),lo(3):hi(3))*grid_area_ratio
      end do
      !$acc end kernels
      !$acc end data
    case(1)
      grid_area_ratio = dl(1)/l(1)
      allocate(p2d(ng(2),ng(3)))
      !$acc data copyout(p2d) async(1)
      !$acc parallel loop gang collapse(2) default(present) private(p2d_s) async(1)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          p2d_s = 0._rp
          !$acc loop reduction(+:p2d_s)
          do i=lo(1),hi(1)
            p2d_s = p2d_s + p(i,j,k)
          end do
          p2d(j,k) = p2d_s
        end do
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p2d) async(1)
      !$acc kernels default(present) async(1)
      do i=lo(1),hi(1)
        p(i,lo(2):hi(2),lo(3):hi(3)) = p2d(lo(2):hi(2),lo(3):hi(3))*grid_area_ratio
      end do
      !$acc end kernels
      !$acc end data
    end select
  end subroutine ave2d_duct
  !
  subroutine filter3d(n,p,pf)
    !
    ! 3D top-hat filter, second-order trapezoidal rule
    ! it is useless to define temporary variables to touch array elements, since
    ! each element is used only once. The current implementation is ~50% the cost of
    ! computing eight cubes (Bae's code and Davidson's book).
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = (8._rp*(p(i,j,k)) + &
                       4._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i,j,k-1) + &
                              p(i+1,j,k) + p(i,j+1,k) + p(i,j,k+1)) + &
                       2._rp*(p(i,j+1,k+1) + p(i+1,j,k+1) + p(i+1,j+1,k) + &
                              p(i,j-1,k-1) + p(i-1,j,k-1) + p(i-1,j-1,k) + &
                              p(i,j-1,k+1) + p(i-1,j,k+1) + p(i-1,j+1,k) + &
                              p(i,j+1,k-1) + p(i+1,j,k-1) + p(i+1,j-1,k)) + &
                       1._rp*(p(i-1,j-1,k-1) + p(i+1,j-1,k-1) + p(i-1,j+1,k-1) + p(i+1,j+1,k-1) + &
                              p(i-1,j-1,k+1) + p(i+1,j-1,k+1) + p(i-1,j+1,k+1) + p(i+1,j+1,k+1)))/64._rp
        end do
      end do
    end do
  end subroutine filter3d
  !
  subroutine extrapolate(n,is_bound,dzci,p,wk,iface,cbc,lwm)
    !
    ! linear extrapolation of wall-parallel velocity
    !
    ! called before filter/strain_rate for 2D filter/one-sided difference. The
    ! extrapolation requires >= 2 off-wall layers. Wall-normal velocity is
    ! extrapolated when stored at cell centers (iface=0), but is not
    ! when stored at cell faces (iface=1,2,3). i.e., no-penetration bc.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(0:)       :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: wk
    integer , intent(in ) :: iface
    character(len=1), intent(in), dimension(0:1,3,3), optional :: cbc
    integer, intent(in), dimension(0:1,3), optional :: lwm
    logical, dimension(0:1,3) :: is_extra
    real(rp) :: dzc(0:n(3)+1),factor0,factor1
    integer :: i,j,k
    !
    dzc = 1._rp/dzci
    !
    if(present(cbc)) then
      factor0 = 1._rp
      factor1 = 1._rp
      is_extra(:,1) = (is_bound(:,1).and.cbc(:,1,1)=='D'.and.iface/=1)
      is_extra(:,2) = (is_bound(:,2).and.cbc(:,2,2)=='D'.and.iface/=2)
      is_extra(:,3) = (is_bound(:,3).and.cbc(:,3,3)=='D'.and.iface/=3)
    elseif(present(lwm)) then
      factor0 = dzc(0)*dzci(1)
      factor1 = dzc(n(3))*dzci(n(3)-1)
      is_extra(:,1) = (is_bound(:,1).and.lwm(:,1)/=0.and.iface/=1)
      is_extra(:,2) = (is_bound(:,2).and.lwm(:,2)/=0.and.iface/=2)
      is_extra(:,3) = (is_bound(:,3).and.lwm(:,3)/=0.and.iface/=3)
    end if
    !$acc kernels default(present) async(1)
    wk(:,:,:) = p(:,:,:)
    !$acc end kernels
    if(is_extra(0,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          wk(0,j,k) = 2._rp*wk(1,j,k) - wk(2,j,k)
        end do
      end do
    end if
    if(is_extra(1,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          wk(n(1)+1,j,k) = 2._rp*wk(n(1),j,k) - wk(n(1)-1,j,k)
        end do
      end do
    end if
    if(is_extra(0,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do i = 0,n(1)+1
          wk(i,0,k) = 2._rp*wk(i,1,k) - wk(i,2,k)
        end do
      end do
    end if
    if(is_extra(1,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do i = 0,n(1)+1
          wk(i,n(2)+1,k) = 2._rp*wk(i,n(2),k) - wk(i,n(2)-1,k)
        end do
      end do
    end if
    if(is_extra(0,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n(2)+1
        do i = 0,n(1)+1
          wk(i,j,0) = (1._rp+factor0)*wk(i,j,1) - factor0*wk(i,j,2)
        end do
      end do
    end if
    if(is_extra(1,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n(2)+1
        do i = 0,n(1)+1
          wk(i,j,n(3)+1) = (1._rp+factor1)*wk(i,j,n(3)) - factor1*wk(i,j,n(3)-1)
        end do
      end do
    end if
  end subroutine extrapolate
  !
  subroutine cmpt_alph2(n,is_bound,cbc,alph2)
    !
    ! compute filter ratio used in the dynamic Smagorinsky model
    ! set alph2=2.52 in the first off-wall layer to yield more accurate results than 4.00
    ! in practical simulations. Specifically, the near-wall velocity profile is more
    ! accurate. The effect also depends on the grid aspect ratio, with more obvious
    ! effects at AR=1 than AR=2. The choice has negligible influence on WRLES.
    !
    implicit none
    integer, intent(in), dimension(3)        :: n
    logical, intent(in), dimension(0:1,3)    :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3), optional :: cbc
    real(rp), intent(out), dimension(0:,0:,0:) :: alph2
    !
#if !defined(_FILTER_2D)
    !$acc kernels default(present) async(1)
    alph2 = 4.00_rp
    !$acc end kernels
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      !$acc kernels default(present) async(1)
      alph2(1,:,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      !$acc kernels default(present) async(1)
      alph2(n(1),:,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,1,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,n(2),:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,:,1) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,:,n(3)) = 2.52_rp
      !$acc end kernels
    end if
#else
    !$acc kernels default(present) async(1)
    alph2 = 2.52_rp
    !$acc end kernels
#endif
  end subroutine cmpt_alph2
  !
  subroutine filter2d(n,p,pf)
    !
    ! 2D top-hat filter, second-order trapezoidal rule
    ! filtering along only the homogeneous directions might be more suitable,
    ! so cs and del can be taken out of the second filtering operation,
    ! as indicated by R. Agrawal (2022)
    ! this subroutine only applies to channel flow at the moment
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = (4._rp*(p(i,j,k)) + &
                       2._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i+1,j,k) + p(i,j+1,k)) + &
                       1._rp*(p(i-1,j-1,k) + p(i+1,j-1,k) + p(i-1,j+1,k) + p(i+1,j+1,k)))/16._rp
        end do
      end do
    end do
  end subroutine filter2d
  !
  subroutine interpolate(n,u,v,w,uc,vc,wc)
    !
    ! interpolate velocity to cell centers, equivalent to reconstruction (FV)
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: uc,vc,wc
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc(i,j,k) = 0.5_rp*(u(i,j,k)+u(i-1,j,k))
          vc(i,j,k) = 0.5_rp*(v(i,j,k)+v(i,j-1,k))
          wc(i,j,k) = 0.5_rp*(w(i,j,k)+w(i,j,k-1))
        end do
      end do
    end do
  end subroutine interpolate
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
    !$acc parallel loop gang default(present) private(del,fd) async(1)
    do k=1,n(3)
      del = (dl(1)*dl(2)*dzf(k))**(1._rp/3._rp)
      !$acc loop vector collapse(2)
      do j=1,n(2)
        do i=1,n(1)
          fd = 1._rp-exp(-dw_plus(i,j,k)/25._rp)
          visct(i,j,k) = (c_smag*del*fd)**2*s0(i,j,k)
        end do
      end do
    end do
  end subroutine sgs_smag
  !
  subroutine cmpt_dw_plus(cbc,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw_plus)
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
    real(rp), intent(out), dimension(0:,0:,0:) :: dw_plus
    real(rp), allocatable, dimension( :, :, :), save :: dw
    logical, save :: is_first = .true.
    real(rp) :: tauw(2),tauw_tot,this_dw,visci
    integer :: i,j,k
    !
    visci = 1._rp/visc
    !
    if(is_first) then
      is_first = .false.
      allocate(dw(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      !$acc enter data create(dw) async(1)
    end if
    !$acc kernels default(present) async(1)
    dw(:,:,:) = big
    !$acc end kernels
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do i = 1,n(1)
            this_dw = dl(1)*(i-0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(n(1),j,k)-v(n(1)+1,j,k)+v(n(1),j-1,k)-v(n(1)+1,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(n(1),j,k)-w(n(1)+1,j,k)+w(n(1),j,k-1)-w(n(1)+1,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do i = 1,n(1)
            this_dw = dl(1)*(n(1)-i+0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do j = 1,n(2)
            this_dw = dl(2)*(j-0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,n(2),k)-u(i,n(2)+1,k)+u(i-1,n(2),k)-u(i-1,n(2)+1,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,n(2),k)-w(i,n(2)+1,k)+w(i,n(2),k-1)-w(i,n(2)+1,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do j = 1,n(2)
            this_dw = dl(2)*(n(2)-j+0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,1)-u(i,j,0)+u(i-1,j,1)-u(i-1,j,0))/dzc(0)
          tauw(2) = visc*0.5_rp*(v(i,j,1)-v(i,j,0)+v(i,j-1,1)-v(i,j-1,0))/dzc(0)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do k = 1,n(3)
            this_dw = zc(k)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    !
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      !$acc parallel loop gang collapse(2) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,n(3))-u(i,j,n(3)+1)+u(i-1,j,n(3))-u(i-1,j,n(3)+1))/dzc(n(3))
          tauw(2) = visc*0.5_rp*(v(i,j,n(3))-v(i,j,n(3)+1)+v(i,j-1,n(3))-v(i,j-1,n(3)+1))/dzc(n(3))
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
          !$acc loop vector
          do k = 1,n(3)
            this_dw = l(3)-zc(k)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
  end subroutine cmpt_dw_plus
  !
  subroutine strain_rate(n,dli,dzci,dzfi,ux,uy,uz,s0,sij)
    !
    ! compute the strain rate field
    !
    ! Sij is first computed at (or averaged to) cell center, then s0=sqrt(2SijSij)
    ! at cell center (Bae and Orlandi's codes). The implementation is efficient, since
    ! it avoids repetitive computation of derivatives. Costa first averages SijSij to
    ! cell center, then computes s0, which leads to larger s0, especially when Sij
    ! at the cell edges have opposite signs. The second and third loops cannot combine.
    ! please modify this back to Costa's version for memory saving
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(0:,0:,0:) :: s0
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: sij
    real(rp) :: dxi,dyi
    integer :: i,j,k
    !
    dxi = dli(1)
    dyi = dli(2)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 0,n(3)
      do j = 0,n(2)
        do i = 0,n(1)
          ! compute at cell edge
          sij(i,j,k,1) = 0.5_rp*((ux(i,j+1,k)-ux(i,j,k))*dyi     + (uy(i+1,j,k)-uy(i,j,k))*dxi)
          sij(i,j,k,2) = 0.5_rp*((ux(i,j,k+1)-ux(i,j,k))*dzci(k) + (uz(i+1,j,k)-uz(i,j,k))*dxi)
          sij(i,j,k,3) = 0.5_rp*((uy(i,j,k+1)-uy(i,j,k))*dzci(k) + (uz(i,j+1,k)-uz(i,j,k))*dyi)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          ! average to cell center
          sij(i,j,k,4) = 0.25_rp*(sij(i,j,k,1)+sij(i-1,j,k,1)+sij(i,j-1,k,1)+sij(i-1,j-1,k,1))
          sij(i,j,k,5) = 0.25_rp*(sij(i,j,k,2)+sij(i-1,j,k,2)+sij(i,j,k-1,2)+sij(i-1,j,k-1,2))
          sij(i,j,k,6) = 0.25_rp*(sij(i,j,k,3)+sij(i,j-1,k,3)+sij(i,j,k-1,3)+sij(i,j-1,k-1,3))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          ! compute at cell center
          sij(i,j,k,1) = (ux(i,j,k)-ux(i-1,j,k))*dxi
          sij(i,j,k,2) = (uy(i,j,k)-uy(i,j-1,k))*dyi
          sij(i,j,k,3) = (uz(i,j,k)-uz(i,j,k-1))*dzfi(k)
        end do
      end do
    end do
    !$acc kernels default(present) async(1)
    s0 = sij(:,:,:,1)**2 + sij(:,:,:,2)**2 + sij(:,:,:,3)**2 + &
        (sij(:,:,:,4)**2 + sij(:,:,:,5)**2 + sij(:,:,:,6)**2)*2._rp
    s0 = sqrt(2._rp*s0)
    !$acc end kernels
  end subroutine strain_rate
end module mod_sgs