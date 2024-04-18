! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_wmodel
  use mod_precision
  use mod_typedef, only: cond_bound
  use mod_param, only: kap_log,b_log,eps
  implicit none
  private
  public cmpt_bcuvw
  contains
    !
  subroutine cmpt_bcuvw(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,h,ind,u,v,w,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw determined via wall model
    !
    ! index 0 must be calculated for the right/front/top walls, but not necessary 
    ! for the opposite walls. However, index 0 for the left/back/bottom walls
    ! is necessary when eddy viscosity is introduced.
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: lwm,ind
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,h
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(cond_bound), intent(inout) :: bcu,bcv,bcw
    real(rp) :: visci,wei,uh,vh,wh,u1,u2,v1,v2,w1,w2,tauw(2),coef
    integer  :: i,j,k,i1,i2,j1,j2,k1,k2
    !
    visci = 1._rp/visc
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      i2 = ind(0,1)
      i1 = ind(0,1) - 1
      coef = (h-(i1-0.5)*dl(1))/dl(1)
      do k = 1,n(3)
        do j = 0,n(2)
          v1 = v(i1,j,k)
          v2 = v(i2,j,k)
          w1 = 0.25_rp*(w(i1,j,k) + w(i1,j+1,k) + w(i1,j,k-1) + w(i1,j+1,k-1))
          w2 = 0.25_rp*(w(i2,j,k) + w(i2,j+1,k) + w(i2,j,k-1) + w(i2,j+1,k-1))
          vh = (1._rp-coef)*v1 + coef*v2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(0,1),vh,wh,h,l(1),visc,tauw)
          bcv%x(j,k,0) = visci*tauw(1)
        end do
      end do
      do k = 0,n(3)
        do j = 1,n(2)
          wei = (zf(k)-zc(k))/dzc(k)
          v1 = (1._rp-wei)*(v(i1,j-1,k) + v(i1,j,k)) + wei*(v(i1,j-1,k+1) + v(i1,j,k+1))
          v2 = (1._rp-wei)*(v(i2,j-1,k) + v(i2,j,k)) + wei*(v(i2,j-1,k+1) + v(i2,j,k+1))
          v1 = 0.5_rp*v1
          v2 = 0.5_rp*v2
          w1 = w(i1,j,k)
          w2 = w(i2,j,k)
          vh = (1._rp-coef)*v1 + coef*v2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(0,1),vh,wh,h,l(1),visc,tauw)
          bcw%x(j,k,0) = visci*tauw(2)
        end do
      end do
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      i2 = ind(1,1)
      i1 = ind(1,1) + 1
      coef = (h-(n(1)-i1+0.5)*dl(1))/dl(1)
      do k = 1,n(3)
        do j = 0,n(2)
          v1 = v(i1,j,k)
          v2 = v(i2,j,k)
          w1 = 0.25_rp*(w(i1,j,k) + w(i1,j+1,k) + w(i1,j,k-1) + w(i1,j+1,k-1))
          w2 = 0.25_rp*(w(i2,j,k) + w(i2,j+1,k) + w(i2,j,k-1) + w(i2,j+1,k-1))
          vh = (1._rp-coef)*v1 + coef*v2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(1,1),vh,wh,h,l(1),visc,tauw)
          bcv%x(j,k,1) = -visci*tauw(1)
        end do
      end do
      do k = 0,n(3)
        do j = 1,n(2)
          wei = (zf(k)-zc(k))/dzc(k)
          v1 = (1._rp-wei)*(v(i1,j-1,k) + v(i1,j,k)) + wei*(v(i1,j-1,k+1) + v(i1,j,k+1))
          v2 = (1._rp-wei)*(v(i2,j-1,k) + v(i2,j,k)) + wei*(v(i2,j-1,k+1) + v(i2,j,k+1))
          v1 = 0.5_rp*v1
          v2 = 0.5_rp*v2
          w1 = w(i1,j,k)
          w2 = w(i2,j,k)
          vh = (1._rp-coef)*v1 + coef*v2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(1,1),vh,wh,h,l(1),visc,tauw)
          bcw%x(j,k,1) = -visci*tauw(2)
        end do
      end do
    end if
    !
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      j2 = ind(0,2)
      j1 = ind(0,2) - 1
      coef = (h-(j1-0.5)*dl(2))/dl(2)
      do k = 1,n(3)
        do i = 0,n(1)
          u1 = u(i,j1,k)
          u2 = u(i,j2,k)
          w1 = 0.25_rp*(w(i,j1,k) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k-1))
          w2 = 0.25_rp*(w(i,j2,k) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k-1))
          uh = (1._rp-coef)*u1 + coef*u2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(0,2),uh,wh,h,l(2),visc,tauw)
          bcu%y(i,k,0) = visci*tauw(1)
        end do
      end do
      do k = 0,n(3)
        do i = 1,n(1)
          wei = (zf(k)-zc(k))/dzc(k)
          u1 = (1._rp-wei)*(u(i-1,j1,k) + u(i,j1,k)) + wei*(u(i-1,j1,k+1) + u(i,j1,k+1))
          u2 = (1._rp-wei)*(u(i-1,j2,k) + u(i,j2,k)) + wei*(u(i-1,j2,k+1) + u(i,j2,k+1))
          u1 = 0.5_rp*u1
          u2 = 0.5_rp*u2
          w1 = w(i,j1,k)
          w2 = w(i,j2,k)
          uh = (1._rp-coef)*u1 + coef*u2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(0,2),uh,wh,h,l(2),visc,tauw)
          bcw%y(i,k,0) = visci*tauw(2)
        end do
      end do
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      j2 = ind(1,2)
      j1 = ind(1,2) + 1
      coef = (h-(n(2)-j1+0.5)*dl(2))/dl(2)
      do k = 1,n(3)
        do i = 0,n(1)
          u1 = u(i,j1,k)
          u2 = u(i,j2,k)
          w1 = 0.25_rp*(w(i,j1,k) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k-1))
          w2 = 0.25_rp*(w(i,j2,k) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k-1))
          uh = (1._rp-coef)*u1 + coef*u2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(1,2),uh,wh,h,l(2),visc,tauw)
          bcu%y(i,k,1) = -visci*tauw(1)
        end do
      end do
      do k = 0,n(3)
        do i = 1,n(1)
          wei = (zf(k)-zc(k))/dzc(k)
          u1 = (1._rp-wei)*(u(i-1,j1,k) + u(i,j1,k)) + wei*(u(i-1,j1,k+1) + u(i,j1,k+1))
          u2 = (1._rp-wei)*(u(i-1,j2,k) + u(i,j2,k)) + wei*(u(i-1,j2,k+1) + u(i,j2,k+1))
          u1 = 0.5_rp*u1
          u2 = 0.5_rp*u2
          w1 = w(i,j1,k)
          w2 = w(i,j2,k)
          uh = (1._rp-coef)*u1 + coef*u2
          wh = (1._rp-coef)*w1 + coef*w2
          call wallmodel(lwm(1,2),uh,wh,h,l(2),visc,tauw)
          bcw%y(i,k,1) = -visci*tauw(2)
        end do
      end do
    end if
    !
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      k2 = ind(0,3)
      k1 = ind(0,3) - 1
      coef = (h-zc(k1))/dzc(k1)
      do j = 1,n(2)
        do i = 0,n(1)
          u1 = u(i,j,k1)
          u2 = u(i,j,k2)
          v1 = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          v2 = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          uh = (1._rp-coef)*u1 + coef*u2
          vh = (1._rp-coef)*v1 + coef*v2
          call wallmodel(lwm(0,3),uh,vh,h,l(3),visc,tauw)
          bcu%z(i,j,0) = visci*tauw(1)
        end do
      end do
      do j = 0,n(2)
        do i = 1,n(1)
          u1 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          u2 = 0.25_rp*(u(i-1,j,k2) + u(i,j,k2) + u(i-1,j+1,k2) + u(i,j+1,k2))
          v1 = v(i,j,k1)
          v2 = v(i,j,k2)
          uh = (1._rp-coef)*u1 + coef*u2
          vh = (1._rp-coef)*v1 + coef*v2
          call wallmodel(lwm(0,3),uh,vh,h,l(3),visc,tauw)
          bcv%z(i,j,0) = visci*tauw(2)
        end do
      end do
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      k2 = ind(1,3)
      k1 = ind(1,3) + 1
      coef = (h-(l(3)-zc(k1)))/(dzc(k2))
      do j = 1,n(2)
        do i = 0,n(1)
          u1 = u(i,j,k1)
          u2 = u(i,j,k2)
          v1 = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          v2 = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          uh = (1._rp-coef)*u1 + coef*u2
          vh = (1._rp-coef)*v1 + coef*v2
          call wallmodel(lwm(1,3),uh,vh,h,l(3),visc,tauw)
          bcu%z(i,j,1) = -visci*tauw(1)
        end do
      end do
      do j = 0,n(2)
        do i = 1,n(1)
          u1 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          u2 = 0.25_rp*(u(i-1,j,k2) + u(i,j,k2) + u(i-1,j+1,k2) + u(i,j+1,k2))
          v1 = v(i,j,k1)
          v2 = v(i,j,k2)
          uh = (1._rp-coef)*u1 + coef*u2
          vh = (1._rp-coef)*v1 + coef*v2
          call wallmodel(lwm(1,3),uh,vh,h,l(3),visc,tauw)
          bcv%z(i,j,1) = -visci*tauw(2)
        end do
      end do
    end if
  end subroutine cmpt_bcuvw
  !
  subroutine wallmodel(mtype,uh,vh,h,l1d,visc,tauw)
    implicit none
    integer, parameter :: WM_LAM = -1, &
                          WM_LOG =  1
    integer, intent(in)  :: mtype
    real(rp), intent(in) :: uh,vh,h,l1d,visc
    real(rp), intent(out), dimension(2) :: tauw
    real(rp) :: upar,utau,f,fp,conv,tauw_old,tauw_tot
    real(rp) :: umax,del
    !
    select case(mtype)
    case(WM_LOG)
      upar = sqrt(uh*uh+vh*vh)
      utau = sqrt(upar/h*visc)
      conv = 1._rp
      do while(conv > 1.e-4)
        ! Newton-Raphson, ~6 iterations
        tauw_old = utau*utau
        f  = upar/utau - 1._rp/kap_log*log(h*utau/visc) - b_log
        fp = -1._rp/utau*(upar/utau + 1._rp/kap_log)
        utau = abs(utau - f/fp) ! robust
        tauw_tot = utau*utau
        conv = abs(tauw_tot-tauw_old)/tauw_old
      end do
      tauw(1)  = tauw_tot*uh/upar
      tauw(2)  = tauw_tot*vh/upar
    case(WM_LAM)
      upar = sqrt(uh*uh+vh*vh)
      del  = 0.5_rp*l1d
      umax = upar/(h/del*(2._rp-h/del))
      tauw_tot = 2._rp/del*umax*visc
      tauw(1)  = tauw_tot*uh/upar
      tauw(2)  = tauw_tot*vh/upar
    end select
    !
  end subroutine wallmodel
end module mod_wmodel