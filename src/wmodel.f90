! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_wmodel    !!!!!rename as wallstress
  use mod_precision
  use mod_typedef, only: cond_bound
  use mod_param, only: kap_log,b_log
  implicit none
  private
  public cmpt_bcuvw,cmpt_bcp
  contains
    !
  subroutine cmpt_bcuvw(n,is_bound,lwm,lo,l,dl,zc,visc,h,u,v,w,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw, determined via bcvel or wall model
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: lwm
    integer , intent(in), dimension(3) :: lo
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: visc,h
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(cond_bound), intent(inout) :: bcu,bcv,bcw
    real(rp) :: wei,uh,vh,wh,u1,u2,v1,v2,w1,w2,tauw(2)
    integer  :: i,j,k,j1,j2,k1,k2
    !
    ! to be simplified
    ! modify bc 1, consistent with bc 0, using h, rather than l-h
    if(is_bound(0,3).and.lwm(0,3)>0) then  !to be improved by modifying input.nml straightly
      !lower wall
      k = 1
      do while(zc(k) < h)
        k = k + 1
      end do
      if(k > n(3)) print *, 'error with wall model'
      k2 = k
      k1 = k - 1
      wei= (h-zc(k1))/(zc(k2)-zc(k1))
      do j = 1,n(2)
        do i = 1,n(1)
          u1 = u(i,j,k1)
          u2 = u(i,j,k2)
          v1 = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          v2 = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          uh = (1._rp - wei)*u1 + wei*u2
          vh = (1._rp - wei)*v1 + wei*v2
          call wallmodel(uh,vh,h,visc,tauw)
          bcu%z(i,j,0) = 1._rp/visc*tauw(1) !tauw, x direction
          !
          u1 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          u2 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          v1 = v(i,j,k1)
          v2 = v(i,j,k2)
          uh = (1._rp - wei)*u1 + wei*u2
          vh = (1._rp - wei)*v1 + wei*v2
          call wallmodel(uh,vh,h,visc,tauw)
          bcv%z(i,j,0) = 1._rp/visc*tauw(2) !tauw, y direction
        end do
      end do
    end if
    if(is_bound(1,3).and.lwm(1,3)>0) then
      !upper wall
      k = n(3)
      do while(zc(k) > l(3)-h)
        k = k - 1
      end do
      if(k < 1) print *, 'error with wall model'
      k2 = k
      k1 = k + 1
      wei= ((l(3)-h)-zc(k1))/(zc(k2)-zc(k1))
      do j = 1,n(2)
        do i = 1,n(1)
          u1 = u(i,j,k1)
          u2 = u(i,j,k2)
          v1 = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          v2 = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          uh = (1._rp - wei)*u1 + wei*u2
          vh = (1._rp - wei)*v1 + wei*v2
          call wallmodel(uh,vh,h,visc,tauw)
          bcu%z(i,j,1) = -1._rp/visc*tauw(1) !tauw, x direction
          !
          u1 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          u2 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          v1 = v(i,j,k1)
          v2 = v(i,j,k2)
          uh = (1._rp - wei)*u1 + wei*u2
          vh = (1._rp - wei)*v1 + wei*v2
          call wallmodel(uh,vh,h,visc,tauw)
          bcv%z(i,j,1) = -1._rp/visc*tauw(2) !tauw, y direction
        end do
      end do
    end if
    !used in square duct (four walls)
    !bug here, should consider non-uniform spacings in the z direction
    if(is_bound(0,2).and.lwm(0,2)>0) then
      !lower wall
      j = 1
      do while((j-0.5)*dl(2) < h)
        j = j + 1
      end do
      if(j > n(2)) print *, 'error with wall model'
      j2 = j
      j1 = j - 1
      wei= (h-(j1-0.5)*dl(2))/((j2-j1)*dl(2))
      do k = 1,n(3)
        do i = 1,n(1)
          u1 = u(i,j1,k)
          u2 = u(i,j2,k)
          w1 = 0.25_rp*(w(i,j1,k) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k))
          w2 = 0.25_rp*(w(i,j2,k) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k))
          uh = (1._rp - wei)*u1 + wei*u2
          wh = (1._rp - wei)*w1 + wei*w2
          call wallmodel(uh,wh,h,visc,tauw)
          bcu%y(i,k,0) = 1._rp/visc*tauw(1) !tauw, x direction
          !
          u1 = 0.25_rp*(u(i-1,j1,k) + u(i,j1,k) + u(i-1,j1,k+1) + u(i,j1,k+1))
          u2 = 0.25_rp*(u(i-1,j2,k) + u(i,j2,k) + u(i-1,j2,k+1) + u(i,j2,k+1))
          w1 = w(i,j1,k)
          w2 = w(i,j2,k)
          uh = (1._rp - wei)*u1 + wei*u2
          wh = (1._rp - wei)*w1 + wei*w2
          call wallmodel(uh,wh,h,visc,tauw)
          bcw%y(i,k,0) = 1._rp/visc*tauw(2) !tauw, z direction
        end do
      end do
    end if
    if(is_bound(1,2).and.lwm(1,2)>0) then
      !upper wall
      j = n(2)
      do while(((j-1+lo(2))-0.5)*dl(2) > l(2)-h)
        j = j - 1
      end do
      if(j < 1) print *, 'error with wall model'
      j2 = j
      j1 = j + 1
      wei= ((l(2)-h)-((j1-1+lo(2))-0.5)*dl(2))/((j2-j1)*dl(2))
      do k = 1,n(3)
        do i = 1,n(1)
          u1 = u(i,j1,k)
          u2 = u(i,j2,k)
          w1 = 0.25_rp*(w(i,j1,k  ) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k))
          w2 = 0.25_rp*(w(i,j2,k  ) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k))
          uh = (1._rp - wei)*u1 + wei*u2
          wh = (1._rp - wei)*w1 + wei*w2
          call wallmodel(uh,wh,h,visc,tauw)
          bcu%y(i,k,1) = -1._rp/visc*tauw(1) !tauw, x direction
          !
          u1 = 0.25_rp*(u(i-1,j1,k  ) + u(i,j1,k  ) + u(i-1,j1,k+1) + u(i,j1,k+1))
          u2 = 0.25_rp*(u(i-1,j2,k  ) + u(i,j2,k  ) + u(i-1,j2,k+1) + u(i,j2,k+1))
          w1 = w(i,j1,k)
          w2 = w(i,j2,k)
          uh = (1._rp - wei)*u1 + wei*u2
          wh = (1._rp - wei)*w1 + wei*w2
          call wallmodel(uh,wh,h,visc,tauw)
          bcw%y(i,k,1) = -1._rp/visc*tauw(2) !tauw, z direction
        end do
      end do
    end if
    !
  end subroutine cmpt_bcuvw
  !
  subroutine cmpt_bcp(cbc,n,bc,is_bound,is_wm,p,bcp)
    !
    ! bcp, determined via bcpre
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3) :: bc
    logical , intent(in), dimension(0:1,3) :: is_bound
    logical , intent(in) :: is_wm
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    type(cond_bound), intent(inout) :: bcp
    !
    if(is_bound(0,1).and.cbc(0,1)/='P') bcp%x(:,:,0) = bc(0,1)
    if(is_bound(1,1).and.cbc(1,1)/='P') bcp%x(:,:,1) = bc(1,1)
    if(is_bound(0,2).and.cbc(0,2)/='P') bcp%y(:,:,0) = bc(0,2)
    if(is_bound(1,2).and.cbc(1,2)/='P') bcp%y(:,:,1) = bc(1,2)
    if(is_bound(0,3).and.cbc(0,3)/='P') bcp%z(:,:,0) = bc(0,3)
    if(is_bound(1,3).and.cbc(1,3)/='P') bcp%z(:,:,1) = bc(1,3)
    !
  end subroutine cmpt_bcp
  !
  ! subroutine cmpt_bctau(cbc,n,bc,is_bound,zc,visc,kap,b,h,u,v,w,bctau1,bctau2)
  !   !
  !   ! compute bctau at cell centers via wall models
  !   !
  !   implicit none
  !   character(len=1), intent(in), dimension(0:1,3,3) :: cbc
  !   integer , intent(in), dimension(3) :: n
  !   real(rp), intent(in), dimension(0:1,3,3) :: bc
  !   logical , intent(in), dimension(0:1,3) :: is_bound
  !   real(rp), intent(in), dimension(0:) :: zc
  !   real(rp), intent(in) :: visc,kap,b,h
  !   real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
  !   type(cond_bound), intent(inout) :: bctau1,bctau2
  !   real(rp) :: wei,uc,vc,tauw(2)
  !   integer  :: i,j,k,k1,k2
  !   !lower wall
  !   k = 1
  !   do while(zc(k) < h)
  !     k = k + 1
  !   end do
  !   k2 = k
  !   k1 = k - 1
  !   wei= (h-zc(k1))/(zc(k2)-zc(k1))
  !   do j = 1,n(2)+1
  !     do i = 1,n(1)+1
  !       uc = (1._rp - wei)*(u(i-1,j  ,k1) + u(i,j,k1)) &
  !                   + wei *(u(i-1,j  ,k2) + u(i,j,k2))
  !       vc = (1._rp - wei)*(v(i  ,j-1,k1) + v(i,j,k1)) &
  !                   + wei *(v(i  ,j-1,k2) + v(i,j,k2))
  !       uc = uc*0.5_rp
  !       vc = vc*0.5_rp
  !       call wallmodel(uc,vc,h,visc,kap,b,tauw)
  !       bctau1%z(i,j,0) = tauw(1)
  !       bctau2%z(i,j,0) = tauw(2)
  !     end do
  !   end do
  !   !upper wall
  !   k = n(3)
  !   do while(zc(k) > 1._rp-h)
  !     k = k - 1
  !   end do
  !   k2 = k
  !   k1 = k + 1
  !   wei= ((1._rp-h)-zc(k1))/(zc(k2)-zc(k1))
  !   do j = 1,n(2)+1
  !     do i = 1,n(1)+1
  !       uc = (1._rp - wei)*(u(i-1,j  ,k1) + u(i,j,k1)) &
  !                   + wei *(u(i-1,j  ,k2) + u(i,j,k2))
  !       vc = (1._rp - wei)*(v(i  ,j-1,k1) + v(i,j,k1)) &
  !                   + wei *(v(i  ,j-1,k2) + v(i,j,k2))
  !       uc = uc*0.5_rp
  !       vc = vc*0.5_rp
  !       call wallmodel(uc,vc,h,visc,kap,b,tauw)
  !       bctau1%z(i,j,1) = tauw(1)
  !       bctau2%z(i,j,1) = tauw(2)
  !     end do
  !   end do
  ! end subroutine cmpt_bctau
  !
  subroutine wallmodel(uh,vh,h,visc,tauw)
    implicit none
    real(rp), intent(in)  :: uh,vh,h,visc
    real(rp), intent(out), dimension(2) :: tauw
    real(rp) :: upar,utau,f,fp,conv,tauw_old,tauw_tot
    real(rp) :: umax
    !
    if(.true.) then
      !log-law profile
      upar = sqrt(uh*uh+vh*vh)
      utau = sqrt(upar/h*visc)
      conv = 1._rp
      do while(conv > 1.e-4)
        tauw_old = utau*utau
        f  = upar/utau - 1._rp/kap_log*log(h*utau/visc) - b_log
        fp = -1._rp/utau*(upar/utau + 1._rp/kap_log)
        utau = abs(utau - f/fp)  !to be improved
        tauw_tot = utau*utau
        conv = abs(tauw_tot-tauw_old)/tauw_old
      end do
      tauw(1)  = tauw_tot*uh/upar
      tauw(2)  = tauw_tot*vh/upar
    else
      !parabolic profile
      upar = sqrt(uh*uh+vh*vh)
      umax = upar/(1._rp-((h-0.5_rp)/0.5_rp)**2)
      tauw_tot = 4._rp*umax*visc
      tauw(1)  = tauw_tot*uh/upar
      tauw(2)  = tauw_tot*vh/upar
    end if
  end subroutine wallmodel
end module mod_wmodel