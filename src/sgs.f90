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
  use mod_typedef, only: cond_bound
  use mod_bound, only: boundp,bounduvw
  use mod_post, only: strain_rate
  implicit none
  private
  public cmpt_sgs
  contains
  !
  subroutine cmpt_sgs(sgstype,n,ng,lo,hi,cbcvel,cbcpre,bcp,nb,is_bound,lwm,l,dl,zc,zf,dzc,dzf, &
                      visc,h,ind,u,v,w,dw,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,visct)
    !
    ! compute subgrid viscosity at cell centers
    ! the current implementation of the dynamcic version is two times the cost of
    ! the static one. Acceleration can be further achieved by calling only one boundp to do
    ! the data exchange of multiple variables. Note that it does not save time to simply mask
    ! wall bc's in boundp
    ! 
    ! The dynamic version yields quite good results for Re_tau=395,550,1000,
    ! with <=5% error in the friction coefficient. Clipping is necessary to avoid negative
    ! eddy viscosity after averaging.
    !
    ! 2D filter is used for the first off-wall layer, as done by Bae, Orlandi, LESGO,
    ! and Balaras (1995), Finite-Difference Computations of High Reynolds Number
    ! Flows Using the Dynamic Subgrid-Scale Model.
    ! It is difficult to do 3D filtering of Sij for the first layer, though
    ! feasible for the velocity.
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    integer , intent(in ), dimension(3) :: n,ng,lo,hi
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    type(cond_bound), intent(in )                :: bcp
    type(cond_bound), intent(inout)              :: bcuf,bcvf,bcwf
    type(cond_bound), intent(in )                :: bcu_mag,bcv_mag,bcw_mag
    integer , intent(in ), dimension(0:1,3)      :: nb,lwm,ind
    logical , intent(in ), dimension(0:1,3)      :: is_bound
    real(rp), intent(in ), dimension(3)          :: l,dl
    real(rp), intent(in ), dimension(0:)         :: zc,zf,dzc,dzf
    real(rp), intent(in )                        :: visc,h
    real(rp), intent(in ), dimension(0:,0:,0:)   :: u,v,w,dw
    real(rp), intent(out), dimension(0:,0:,0:)   :: visct
    !
    real(rp), allocatable, dimension(:,:,:)  , save :: dw_plus,s0,uc,vc,wc,uf,vf,wf, &
                                                       wk,wk1,wk2,wk3,alph2
    real(rp), allocatable, dimension(:,:,:,:), save :: sij,lij,mij
    real(rp), dimension(3)        :: dli
    real(rp), dimension(0:n(3)+1) :: dzci,dzfi
    logical, save :: is_first = .true.
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
      if(is_first) then
        is_first = .false.
        allocate(wk     (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk1    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk2    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk3    (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 s0     (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 dw_plus(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 sij    (0:n(1)+1,0:n(2)+1,0:n(3)+1,6))
      end if
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      call cmpt_dw_plus(cbcvel,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw,dw_plus) ! need modification on tauw?
      call sgs_smag(n,dl,dzf,dw_plus,s0,visct)
    case('dsmag')
      if(is_first) then
        is_first = .false.
        allocate(wk (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk1(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk2(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk3(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 s0 (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 uc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 uf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 sij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 lij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 mij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6))
        allocate(alph2(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        call cmpt_alph2(n,is_bound,cbcvel,alph2)
      end if
      !
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      !
      visct = s0
      !
      ! Lij
      !
      call interpolate(n,u,v,w,uc,vc,wc)
      ! periodic and patched bc's are assigned, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,uc)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,vc)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,wc)
#if !defined(_FILTER_2D)
      call extrapolate(n,is_bound,dzci,uc*uc,wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,vc*vc,wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wc*wc,wk3,iface=0,cbc=cbcvel)
      call filter(wk1,lij(:,:,:,1))
      call filter(wk2,lij(:,:,:,2))
      call filter(wk3,lij(:,:,:,3))
      call extrapolate(n,is_bound,dzci,uc*vc,wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,uc*wc,wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,vc*wc,wk3,iface=0,cbc=cbcvel)
      call filter(wk1,lij(:,:,:,4))
      call filter(wk2,lij(:,:,:,5))
      call filter(wk3,lij(:,:,:,6))
      call extrapolate(n,is_bound,dzci,uc,wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,vc,wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wc,wk3,iface=0,cbc=cbcvel)
      call filter(wk1,uf)
      call filter(wk2,vf)
      call filter(wk3,wf)
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
      ! periodic and patched bc's are assigned, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,s0)
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,1))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,2))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,3))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,4))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,5))
      call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,sij(:,:,:,6))
#if !defined(_FILTER_2D)
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,1),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,2),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,3),wk3,iface=0,cbc=cbcvel)
      call filter(wk1,mij(:,:,:,1))
      call filter(wk2,mij(:,:,:,2))
      call filter(wk3,mij(:,:,:,3))
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,4),wk1,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,5),wk2,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,s0*sij(:,:,:,6),wk3,iface=0,cbc=cbcvel)
      call filter(wk1,mij(:,:,:,4))
      call filter(wk2,mij(:,:,:,5))
      call filter(wk3,mij(:,:,:,6))
      call extrapolate(n,is_bound,dzci,u,wk1,iface=1,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,v,wk2,iface=2,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,w,wk3,iface=3,cbc=cbcvel)
      call filter(wk1,uf)
      call filter(wk2,vf)
      call filter(wk3,wf)
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
      ! when using no wall model, all bc's are used for computing strain rate.
      ! When using a wall model, wall bc's for wall-parallel valocity are not used due
      ! to one-sided finite difference. bcuf, bcvf and bcwf are equal to bcu/v/w=0 for
      ! no-slip walls. When using a wall model, only the wall bc's for the wall-normal
      ! velicity are essentially used. Filtered velocity should satisfy the wall bc's,
      ! evidenced by that each mode satisfies the no-slip/no-penetration bc's
      ! when a spectral filter is applied.
      !
      call bounduvw(cbcvel,n,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,nb,is_bound,lwm,l,dl,zc,zf,dzc,dzf, &
                    visc,h,ind,.true.,.false.,uf,vf,wf)
      call extrapolate(n,is_bound,dzci,uf,wk1,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,vf,wk2,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wf,wk3,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk1,wk2,wk3,s0,sij)
      do m = 1,6
        mij(:,:,:,m) = 2._rp*(mij(:,:,:,m)-alph2(:,:,:)*s0(:,:,:)*sij(:,:,:,m))
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
#if defined(_CHANNEL)
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,s0)
#elif defined(_DUCT)
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,s0)
#elif defined(_CAVITY)
      !
#endif
      visct = visct*s0
      s0(:,:,:) = mij(:,:,:,1)*mij(:,:,:,1) + &
                  mij(:,:,:,2)*mij(:,:,:,2) + &
                  mij(:,:,:,3)*mij(:,:,:,3) + &
                 (mij(:,:,:,4)*mij(:,:,:,4) + &
                  mij(:,:,:,5)*mij(:,:,:,5) + &
                  mij(:,:,:,6)*mij(:,:,:,6))*2._rp
#if defined(_CHANNEL)
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,s0)
#elif defined(_DUCT)
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,s0)
#elif defined(_CAVITY)
      !
#endif
      visct = visct/s0
      visct = max(visct,0._rp)
    case('amd')
      print*, 'ERROR: AMD model not yet implemented'
    case default
      print*, 'ERROR: unknown SGS model'
    end select
  end subroutine cmpt_sgs
  !
  subroutine strain_rate_old(n,dli,dzci,dzfi,is_bound,lwm,ux,uy,uz,s0,sij)
    !
    ! Sij should be first computed at (or averaged to) cell center, then s0=sqrt(2SijSij)
    ! at cell center. This implementation is also adopted by Bae and Orlandi. Costa averages
    ! SijSij to cell center first, then computes s0, which always leads to larger s0,
    ! especially when Sij at the cell edges have opposite signs. The current implementation
    ! avoids repetitive computation of derivatives, so it is much more efficient. Note that
    ! three seperate loops are required; the second and third loops (special treatment does
    ! not count) cannot be combined.
    !
    ! when a wall model is applied, the first layer of cells is large that discontinuity
    ! appears near the wall, one-sided average (difference) should be used; averaging
    ! should not be done across the discontinuity, since the wall-normal derivatives can be
    ! as different as several orders of magnitude. Using second-order one-sided approximation
    ! does not bring any benefit.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    integer , intent(in ), dimension(0:1,3)    :: lwm
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(0:,0:,0:) :: s0
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: sij
    real(rp) :: dxi,dyi
    integer :: i,j,k
    !
    real(rp) :: d12,d13,f1,f2,f3,dfdz,dzc(0:n(3)+1)
    dzc = 1./dzci
    !
    dxi = dli(1)
    dyi = dli(2)
    !
    ! compute s0 = sqrt(2*sij*sij), where sij = (1/2)(du_i/dx_j + du_j/dx_i)
    !
    !$acc parallel loop collapse(3) default(present) private(s11,s12,s13,s22,s23,s33)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(s11,s12,s13,s22,s23,s33)
    ! compute at cell edge
    do k = 0,n(3)
      do j = 0,n(2)
        do i = 0,n(1)
          sij(i,j,k,1) = 0.5_rp*((ux(i,j+1,k)-ux(i,j,k))*dyi     + (uy(i+1,j,k)-uy(i,j,k))*dxi)
          sij(i,j,k,2) = 0.5_rp*((ux(i,j,k+1)-ux(i,j,k))*dzci(k) + (uz(i+1,j,k)-uz(i,j,k))*dxi)
          sij(i,j,k,3) = 0.5_rp*((uy(i,j,k+1)-uy(i,j,k))*dzci(k) + (uz(i,j+1,k)-uz(i,j,k))*dyi)
        end do
      end do
    end do
    ! extrapolate to walls
    ! better to extrapolate the wall-parallel velocity components to the ghost cells,
    ! rather than to assign Sij at the walls. This makes the strain rate calculation clean
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      do k = 0,n(3)
        do j = 0,n(2)
          sij(0,j,k,1) = sij(1,j,k,1)
          sij(0,j,k,2) = sij(1,j,k,2)
        end do
      end do
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      do k = 0,n(3)
        do j = 0,n(2)
          sij(n(1),j,k,1) = sij(n(1)-1,j,k,1)
          sij(n(1),j,k,2) = sij(n(1)-1,j,k,2)
        end do
      end do
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      do k = 0,n(3)
        do i = 0,n(1)
          sij(i,0,k,1) = sij(i,1,k,1)
          sij(i,0,k,3) = sij(i,1,k,3)
        end do
      end do
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      do k = 0,n(3)
        do i = 0,n(1)
          sij(i,n(2),k,1) = sij(i,n(2)-1,k,1)
          sij(i,n(2),k,3) = sij(i,n(2)-1,k,3)
        end do
      end do
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      do j = 0,n(2)
        do i = 0,n(1)
          sij(i,j,0,2) = sij(i,j,1,2)
          sij(i,j,0,3) = sij(i,j,1,3)
        end do
      end do
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      do j = 0,n(2)
        do i = 0,n(1)
          sij(i,j,n(3),2) = sij(i,j,n(3)-1,2)
          sij(i,j,n(3),3) = sij(i,j,n(3)-1,3)
        end do
      end do
    end if
    ! average to cell center
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          sij(i,j,k,4) = 0.25_rp*(sij(i,j,k,1)+sij(i-1,j,k,1)+sij(i,j-1,k,1)+sij(i-1,j-1,k,1))
          sij(i,j,k,5) = 0.25_rp*(sij(i,j,k,2)+sij(i-1,j,k,2)+sij(i,j,k-1,2)+sij(i-1,j,k-1,2))
          sij(i,j,k,6) = 0.25_rp*(sij(i,j,k,3)+sij(i,j-1,k,3)+sij(i,j,k-1,3)+sij(i,j-1,k-1,3))
        end do
      end do
    end do
    ! compute at cell center
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          sij(i,j,k,1) = (ux(i,j,k)-ux(i-1,j,k))*dxi
          sij(i,j,k,2) = (uy(i,j,k)-uy(i,j-1,k))*dyi
          sij(i,j,k,3) = (uz(i,j,k)-uz(i,j,k-1))*dzfi(k)
        end do
      end do
    end do
    !
    s0 = sij(:,:,:,1)**2 + sij(:,:,:,2)**2 + sij(:,:,:,3)**2 + &
        (sij(:,:,:,4)**2 + sij(:,:,:,5)**2 + sij(:,:,:,6)**2)*2._rp
    s0 = sqrt(2._rp*s0)
  end subroutine strain_rate_old
  !
  subroutine extrapolate_old(n,dli,dzci,dzfi,is_bound,lwm,ux,uy,uz,wk1,wk2,wk3)
    !
    ! linear extrapolation of wall-parallel velocity to ghost cells
    ! for a bottom wall, u(0:n(1),1:n(2),0) and v(1:n(1),0:n(2),0) must be assigned.
    !
    ! when a wall model is applied, the first-layer wall-parallel velocity is large
    ! that discontinuity appears near the wall, one-sided difference should be used.
    ! Averaging should not be done across the discontinuity, because the wall-normal
    ! derivatives can be as different as several orders of magnitude. Application of
    ! second-order one-sided difference does not bring benefits.
    !
    ! check channel, duct and cavity on nonuniform grid
    ! try to remove this subroutine
    ! confirm full index get the same results at used points.
    ! then, can be combined.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    integer , intent(in ), dimension(0:1,3)    :: lwm
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(0:,0:,0:) :: wk1,wk2,wk3
    integer  :: i,j,k
    real(rp) :: dzc(0:n(3)+1),factor
    !
    dzc = 1._rp/dzci
    wk1 = ux
    wk2 = uy
    wk3 = uz
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      factor = 1._rp
      wk2(0,:,:) = wk2(1,:,:) - factor*(wk2(2,:,:)-wk2(1,:,:))
      wk3(0,:,:) = wk3(1,:,:) - factor*(wk3(2,:,:)-wk3(1,:,:))
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      factor = 1._rp
      wk2(n(1)+1,:,:) = wk2(n(1),:,:) + factor*(wk2(n(1),:,:)-wk2(n(1)-1,:,:))
      wk3(n(1)+1,:,:) = wk3(n(1),:,:) + factor*(wk3(n(1),:,:)-wk3(n(1)-1,:,:))
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      factor = 1._rp
      wk1(:,0,:) = wk1(:,1,:) - factor*(wk1(:,2,:)-wk1(:,1,:))
      wk3(:,0,:) = wk3(:,1,:) - factor*(wk3(:,2,:)-wk3(:,1,:))
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      factor = 1._rp
      wk1(:,n(2)+1,:) = wk1(:,n(2),:) + factor*(wk1(:,n(2),:)-wk1(:,n(2)-1,:))
      wk3(:,n(2)+1,:) = wk3(:,n(2),:) + factor*(wk3(:,n(2),:)-wk3(:,n(2)-1,:))
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      factor = dzci(1)*dzc(0)
      wk1(:,:,0) = wk1(:,:,1) - factor*(wk1(:,:,2)-wk1(:,:,1))
      wk2(:,:,0) = wk2(:,:,1) - factor*(wk2(:,:,2)-wk2(:,:,1))
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      factor = dzci(n(3)-1)*dzc(n(3))
      wk1(:,:,n(3)+1) = wk1(:,:,n(3)) + factor*(wk1(:,:,n(3))-wk1(:,:,n(3)-1))
      wk2(:,:,n(3)+1) = wk2(:,:,n(3)) + factor*(wk2(:,:,n(3))-wk2(:,:,n(3)-1))
    end if
  end subroutine extrapolate_old
  !
  subroutine extrapolate_old_old(n,dli,dzci,dzfi,is_bound,lwm,ux,uy,uz,wk1,wk2,wk3)
    !
    ! linear extrapolation of wall-parallel velocity to ghost cells
    ! for a bottom wall, u(0:n(1),1:n(2),0) and v(1:n(1),0:n(2),0) must be assigned.
    !
    ! when a wall model is applied, the first-layer wall-parallel velocity is large
    ! that discontinuity appears near the wall, one-sided difference should be used.
    ! Averaging should not be done across the discontinuity, because the wall-normal
    ! derivatives can be as different as several orders of magnitude. Application of
    ! second-order one-sided difference does not bring benefits.
    !
    ! check channel, duct and cavity on nonuniform grid
    ! try to remove this subroutine
    ! confirm full index get the same results at used points.
    ! then, can be combined.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    integer , intent(in ), dimension(0:1,3)    :: lwm
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(0:,0:,0:) :: wk1,wk2,wk3
    integer  :: i,j,k
    real(rp) :: dzc(0:n(3)+1),factor
    !
    dzc = 1._rp/dzci
    wk1 = ux
    wk2 = uy
    wk3 = uz
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      factor = 1._rp
      do k = 1,n(3)
        do j = 0,n(2)
          wk2(0,j,k) = wk2(1,j,k) - factor*(wk2(2,j,k)-wk2(1,j,k))
        end do
      end do
      do k = 0,n(3)
        do j = 1,n(2)
          wk3(0,j,k) = wk3(1,j,k) - factor*(wk3(2,j,k)-wk3(1,j,k))
        end do
      end do
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      factor = 1._rp
      do k = 1,n(3)
        do j = 0,n(2)
          wk2(n(1)+1,j,k) = wk2(n(1),j,k) + factor*(wk2(n(1),j,k)-wk2(n(1)-1,j,k))
        end do
      end do
      do k = 0,n(3)
        do j = 1,n(2)
          wk3(n(1)+1,j,k) = wk3(n(1),j,k) + factor*(wk3(n(1),j,k)-wk3(n(1)-1,j,k))
        end do
      end do
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      do k = 1,n(3)
        do i = 0,n(1)
          wk1(i,0,k) = wk1(i,1,k) - factor*(wk1(i,2,k)-wk1(i,1,k))
        end do
      end do
      do k = 0,n(3)
        do i = 1,n(1)
          wk3(i,0,k) = wk3(i,1,k) - factor*(wk3(i,2,k)-wk3(i,1,k))
        end do
      end do
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      do k = 1,n(3)
        do i = 0,n(1)
          wk1(i,n(2)+1,k) = wk1(i,n(2),k) + factor*(wk1(i,n(2),k)-wk1(i,n(2)-1,k))
        end do
      end do
      do k = 0,n(3)
        do i = 1,n(1)
          wk3(i,n(2)+1,k) = wk3(i,n(2),k) + factor*(wk3(i,n(2),k)-wk3(i,n(2)-1,k))
        end do
      end do
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      factor = dzci(1)*dzc(0)
      do j = 1,n(2)
        do i = 0,n(1)
          wk1(i,j,0) = wk1(i,j,1) - factor*(wk1(i,j,2)-wk1(i,j,1))
        end do
      end do
      do j = 0,n(2)
        do i = 1,n(1)
          wk2(i,j,0) = wk2(i,j,1) - factor*(wk2(i,j,2)-wk2(i,j,1))
        end do
      end do
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      factor = dzci(n(3)-1)*dzc(n(3))
      do j = 1,n(2)
        do i = 0,n(1)
          wk1(i,j,n(3)+1) = wk1(i,j,n(3)) + factor*(wk1(i,j,n(3))-wk1(i,j,n(3)-1))
        end do
      end do
      do j = 0,n(2)
        do i = 1,n(1)
          wk2(i,j,n(3)+1) = wk2(i,j,n(3)) + factor*(wk2(i,j,n(3))-wk2(i,j,n(3)-1))
        end do
      end do
    end if
  end subroutine extrapolate_old_old
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
    real(dp) :: grid_area_ratio
    integer :: i,j,k
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      allocate(p2d(ng(1),ng(3)))
      !
      p2d(:,:) = 0._rp
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          do j=lo(2),hi(2)
            p2d(i,k) = p2d(i,k) + p(i,j,k)
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(1)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      p2d(:,:) = p2d(:,:)*grid_area_ratio
      do j=lo(2),hi(2)
        p(lo(1):hi(1),j,lo(3):hi(3)) = p2d(lo(1):hi(1),lo(3):hi(3))
      end do
    case(1)
      grid_area_ratio = dl(1)/l(1)
      allocate(p2d(ng(2),ng(3)))
      !
      p2d(:,:) = 0._rp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p2d(j,k) = p2d(j,k) + p(i,j,k)
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      p2d(:,:) = p2d(:,:)*grid_area_ratio
      do i=lo(1),hi(1)
        p(i,lo(2):hi(2),lo(3):hi(3)) = p2d(lo(2):hi(2),lo(3):hi(3))
      end do
    end select
  end subroutine ave2d_duct
  !
  subroutine filter(p,pf)
    !
    ! 3D top-hat filter, second-order trapezoidal rule
    ! it is useless to define temporary variables to touch array elements, since
    ! each element is used only once. The current implementation is ~50% the cost of
    ! computing eight cubes (Bae's code and Davidson's book).
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
          pf(i,j,k) = 8._rp*(p(i,j,k)) + &
                      4._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i,j,k-1) + &
                             p(i+1,j,k) + p(i,j+1,k) + p(i,j,k+1)) + &
                      2._rp*(p(i,j+1,k+1) + p(i+1,j,k+1) + p(i+1,j+1,k) + &
                             p(i,j-1,k-1) + p(i-1,j,k-1) + p(i-1,j-1,k) + &
                             p(i,j-1,k+1) + p(i-1,j,k+1) + p(i-1,j+1,k) + &
                             p(i,j+1,k-1) + p(i+1,j,k-1) + p(i+1,j-1,k)) + &
                      1._rp*(p(i-1,j-1,k-1) + p(i+1,j-1,k-1) + p(i-1,j+1,k-1) + p(i+1,j+1,k-1) + &
                             p(i-1,j-1,k+1) + p(i+1,j-1,k+1) + p(i-1,j+1,k+1) + p(i+1,j+1,k+1))
        end do
      end do
    end do
    pf = pf/64._rp
  end subroutine filter
    !
  subroutine filter_old(cbc,is_bound,p,pf,iface)
    !
    ! 3D top-hat filter, second-order trapezoidal rule
    ! we do not define temporary variables to touch array elements, since each element is
    ! used only once. The current implementation's cost is ~50% of that of computing eight
    ! cubes as done in Bae's code and Davidson's book. 2D filtering is always used for the
    ! first off-wall layer, and it is achieved via linear extrapolation to the ghost
    ! points from the interior, so the filtering operation requires >= 2 off-wall layers.
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    logical , intent(in), dimension(0:1,3)     :: is_bound
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer , intent(in), optional :: iface
    integer :: n(3),i,j,k
    !
    real(rp), allocatable, dimension(:,:,:), save :: wk
    logical, save :: is_first = .true.
    !
    n = shape(p)-2
    !
    if(is_first) then
      is_first = .false.
      allocate(wk(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    end if
    wk = p
    pf = 0._rp
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D'.and.iface/=1) then
      wk(0,:,:) = 2._rp*wk(1,:,:) - wk(2,:,:)
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D'.and.iface/=1) then
      wk(n(1)+1,:,:) = 2._rp*wk(n(1),:,:) - wk(n(1)-1,:,:)
    end if
    if(is_bound(0,2).and.cbc(0,2,2)=='D'.and.iface/=2) then
      wk(:,0,:) = 2._rp*wk(:,1,:) - wk(:,2,:)
    end if
    if(is_bound(1,2).and.cbc(1,2,2)=='D'.and.iface/=2) then
      wk(:,n(2)+1,:) = 2._rp*wk(:,n(2),:) - wk(:,n(2)-1,:)
    end if
    if(is_bound(0,3).and.cbc(0,3,3)=='D'.and.iface/=3) then
      wk(:,:,0) = 2._rp*wk(:,:,1) - wk(:,:,2)
    end if
    if(is_bound(1,3).and.cbc(1,3,3)=='D'.and.iface/=3) then
      wk(:,:,n(3)+1) = 2._rp*wk(:,:,n(3)) - wk(:,:,n(3)-1)
    end if
    !
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = 8._rp*(wk(i,j,k)) + &
                      4._rp*(wk(i-1,j,k) + wk(i,j-1,k) + wk(i,j,k-1) + &
                             wk(i+1,j,k) + wk(i,j+1,k) + wk(i,j,k+1)) + &
                      2._rp*(wk(i,j+1,k+1) + wk(i+1,j,k+1) + wk(i+1,j+1,k) + &
                             wk(i,j-1,k-1) + wk(i-1,j,k-1) + wk(i-1,j-1,k) + &
                             wk(i,j-1,k+1) + wk(i-1,j,k+1) + wk(i-1,j+1,k) + &
                             wk(i,j+1,k-1) + wk(i+1,j,k-1) + wk(i+1,j-1,k)) + &
                      1._rp*(wk(i-1,j-1,k-1) + wk(i+1,j-1,k-1) + wk(i-1,j+1,k-1) + wk(i+1,j+1,k-1) + &
                             wk(i-1,j-1,k+1) + wk(i+1,j-1,k+1) + wk(i-1,j+1,k+1) + wk(i+1,j+1,k+1))
          pf(i,j,k) = pf(i,j,k)/64._rp
        end do
      end do
    end do
    ! ! first off-wall layer, zc(1) and zc(n(3))
    ! if(is_fil2d_wall) then
    !   ! bottom wall
    !   do j = 1,n(2)
    !     do i = 1,n(1)
    !       pf(i,j,1) = 4._rp*(p(i,j,1)) + &
    !                   2._rp*(p(i-1,j,1) + p(i,j-1,1) + p(i+1,j,1) + p(i,j+1,1)) + &
    !                   1._rp*(p(i-1,j-1,1) + p(i+1,j-1,1) + p(i-1,j+1,1) + p(i+1,j+1,1))
    !       pf(i,j,1) = pf(i,j,1)/16._rp
    !     end do
    !   end do
    !   ! top wall
    !   do j = 1,n(2)
    !     do i = 1,n(1)
    !       pf(i,j,n(3)) = 4._rp*(p(i,j,n(3))) + &
    !                      2._rp*(p(i-1,j,n(3)) + p(i,j-1,n(3)) + p(i+1,j,n(3)) + p(i,j+1,n(3))) + &
    !                      1._rp*(p(i-1,j-1,n(3)) + p(i+1,j-1,n(3)) + p(i-1,j+1,n(3)) + p(i+1,j+1,n(3)))
    !       pf(i,j,n(3)) = pf(i,j,n(3))/16._rp
    !     end do
    !   end do
    ! end if
  end subroutine filter_old
  !
  subroutine extrapolate(n,is_bound,dzci,p,wk,iface,cbc,lwm)
    !
    ! linear extrapolation of wall-parallel velocity
    !
    ! called before filter/strain_rate for 2D filter/one-sided difference. The
    ! extrapolation requires >= 2 off-wall layers. Wall-normal velocity is
    ! extrapolated when it is stored at cell centers (iface=0), but is not
    ! when it is stored at cell faces (iface=1,2,3) and keeps no-penetration bc.
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
    logical, dimension(0:1,3) :: is_ext
    real(rp) :: dzc(0:n(3)+1),factor(0:1)
    !
    dzc = 1._rp/dzci
    wk = p
    !
    if(present(cbc)) then
      factor(0) = 1._rp
      factor(1) = 1._rp
      is_ext(:,1) = is_bound(:,1).and.cbc(:,1,1)=='D'.and.iface/=1
      is_ext(:,2) = is_bound(:,2).and.cbc(:,2,2)=='D'.and.iface/=2
      is_ext(:,3) = is_bound(:,3).and.cbc(:,3,3)=='D'.and.iface/=3
    elseif(present(lwm)) then
      factor(0) = dzc(0)*dzci(1)
      factor(1) = dzc(n(3))*dzci(n(3)-1)
      is_ext(:,1) = is_bound(:,1).and.lwm(:,1)/=0.and.iface/=1
      is_ext(:,2) = is_bound(:,2).and.lwm(:,2)/=0.and.iface/=2
      is_ext(:,3) = is_bound(:,3).and.lwm(:,3)/=0.and.iface/=3
    end if
    !
    if(is_ext(0,1)) wk(0     ,:     ,:     ) =             2._rp*wk(1   ,:   ,:   ) -           wk(2     ,:     ,:     )
    if(is_ext(1,1)) wk(n(1)+1,:     ,:     ) =             2._rp*wk(n(1),:   ,:   ) -           wk(n(1)-1,:     ,:     )
    if(is_ext(0,2)) wk(:     ,0     ,:     ) =             2._rp*wk(:   ,1   ,:   ) -           wk(:     ,2     ,:     )
    if(is_ext(1,2)) wk(:     ,n(2)+1,:     ) =             2._rp*wk(:   ,n(2),:   ) -           wk(:     ,n(2)-1,:     )
    if(is_ext(0,3)) wk(:     ,:     ,0     ) = (1._rp+factor(0))*wk(:   ,:   ,1   ) - factor(0)*wk(:     ,:     ,2     )
    if(is_ext(1,3)) wk(:     ,:     ,n(3)+1) = (1._rp+factor(1))*wk(:   ,:   ,n(3)) - factor(1)*wk(:     ,:     ,n(3)-1)
  end subroutine extrapolate
  !
  subroutine cmpt_alph2(n,is_bound,cbc,alph2)
    !
    ! compute filter ratio, alph2, used in the dynamic Smagorinsky model
    ! using alph2=2.52 in the first off-wall layer yields more accurate results than 4.00
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
    alph2 = 4.00_rp
    if(is_bound(0,1).and.cbc(0,1,1)=='D') alph2(1   ,:,:) = 2.52_rp
    if(is_bound(1,1).and.cbc(1,1,1)=='D') alph2(n(1),:,:) = 2.52_rp
    if(is_bound(0,2).and.cbc(0,2,2)=='D') alph2(:,1   ,:) = 2.52_rp
    if(is_bound(1,2).and.cbc(1,2,2)=='D') alph2(:,n(2),:) = 2.52_rp
    if(is_bound(0,3).and.cbc(0,3,3)=='D') alph2(:,:,1   ) = 2.52_rp
    if(is_bound(1,3).and.cbc(1,3,3)=='D') alph2(:,:,n(3)) = 2.52_rp
#else
    alph2 = 2.52_rp
#endif
  end subroutine cmpt_alph2
  !
  subroutine filter2d(p,pf)
    !
    ! 2D top-hat filter, second-order trapezoidal rule
    ! filtering along only the homogeneous directions might be more suitable,
    ! so cs and del can be taken out of the second filtering operation,
    ! as indicated by R. Agrawal (2022)
    ! this subroutine only applies to channel flow at the moment
    !
    implicit none
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer :: n(3),i,j,k
    !
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
  subroutine filter_old_old(p,pf,is_fil2d_wall)
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
  end subroutine filter_old_old
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