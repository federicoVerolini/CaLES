! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_wmodel
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan,is_finite => ieee_is_finite
  use mpi
  use mod_precision
  use mod_typedef, only: cond_bound
  use mod_param, only: kap_log,b_log,eps
  implicit none
  private
  public updt_wallmodelbc
  contains
    !
  subroutine updt_wallmodelbc(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,h,ind,u,v,w,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw determined via wall model
    ! wall-parallel velocity at ghost cells is used only for computing the viscous terms,
    ! including its left- and right-hand sides. It is not used for computing the convective
    ! terms, or in the correction procedure. When a wall model is implemented, we can
    ! (1) use Neumann bc, by modifying the ghost cell wall-parallel velocity
    ! (2) use no-slip bc, but replace the viscous stress on the wall
    ! For implicit schemes, when the second approach is used, one must care about the extra
    ! term added to the right-hand-side for the first-layer of cells. The extra term must be
    ! computed using velocity gradient (neumann bc), rather than no-slip bc. In contrast,
    ! the first method completely lets the boundary condition serve for correctly
    ! computing the viscous term, including both the left- and right-hand sides, at the first
    ! off-wall layer of cells.
    ! In WMLES, the wall should be considered as a no-slip wall with modified (more accurate)
    ! wall stress. When the wall stress is required, such as viscous stress, it should be
    ! regarded as a neumann bc. When the wall velocity is required, such as work,
    ! it should be regarded as a no-slip wall, so zero work is done at the wall. When filtering
    ! is required, the wall is a slip wall, with the velocity extrapolated from the interior.
    ! Hence, a ghost point can have three possible values for different purposes.
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
    real(rp) :: wei,coef,uh,vh,wh,u1,u2,v1,v2,w1,w2,tauw(2),visci
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
  end subroutine updt_wallmodelbc
  !
  subroutine wallmodel(mtype,uh,vh,h,l1d,visc,tauw)
    !
    ! Newton-Raphson, 3~7 iters when initialized with an assumed linear profile
    ! between the wall and the wall model height. It is a bad idea to initialize
    ! with an assumed linear profile below the first cell center, because its
    ! slope is quite uncertain. The uncertainty is from the erroneous velocity
    ! at the first cell center, which could make the computed gradient very
    ! different from the converged value. We do not save tauw_tot from previous
    ! step, so those values are not directly available. At zero velocity, utau
    ! is visc/h*exp(-kap_log*b_log). It is necessary to use abs to avoid negative
    ! values of utau. Note that arrays (u/v/w) are initialized as zero (ghost points),
    ! which makes a reasonable guess at end points.
    ! The convergence criterion is |tauw/taw_old-1.0| < 1.0e-4, corresponding
    ! approximately to |utau/utau_old-1.0| < 0.5e-4.
    ! For a channel, it is good to use the computed utau available at the nearest
    ! grid point, which helps reduce ~50% of the iterations. However, several "if"
    ! statements have to be introduced for special points (square duct/cavity),
    ! which is cubersome, considering that the current purpose is to develop
    ! machine learning wall models, which do not involve iterations, and that
    ! the simple log-law wall model accounts less than 1% of the total time.
    !
    implicit none
    integer, parameter :: WM_LAM = -1, &
                          WM_LOG =  1
    integer, intent(in)  :: mtype
    real(rp), intent(in) :: uh,vh,h,l1d,visc
    real(rp), intent(out), dimension(2) :: tauw
    real(rp) :: upar,utau,f,fp,conv,utau_old,tauw_tot
    real(rp) :: umax,del
    !
    integer :: i,ierr
    !
    select case(mtype)
    case(WM_LOG)
      conv = 1._rp
      upar = sqrt(uh*uh+vh*vh)
      utau = max(sqrt(upar/h*visc),visc/h*exp(-kap_log*b_log))
      do while(conv>0.5e-4_rp)
        utau_old = utau
        f  = upar/utau-1._rp/kap_log*log(h*utau/visc)-b_log
        fp = -1._rp/utau*(upar/utau+1._rp/kap_log)
        utau = abs(utau-f/fp)
        conv = abs(utau/utau_old-1._rp)
      end do
      tauw_tot = utau*utau
      tauw(1) = tauw_tot*uh/(upar+eps)
      tauw(2) = tauw_tot*vh/(upar+eps)
    case(WM_LAM)
      upar = sqrt(uh*uh+vh*vh)
      del  = 0.5_rp*l1d
      umax = upar/(h/del*(2._rp-h/del))
      tauw_tot = 2._rp/del*umax*visc
      tauw(1) = tauw_tot*uh/(upar+eps)
      tauw(2) = tauw_tot*vh/(upar+eps)
    end select
    !
  end subroutine wallmodel
  !
  ! subroutine correc_1st_point(n,is_bound,lwm,dl,dzc,dzf,u,v,w)
  !   !
  !   ! extrapolate the first off-wall point from the second and third off-wall points
  !   ! This corrects only the inner points (1:n), ghost points are handled by bounduvw
  !   ! as normally done. The wall-normal velocity is not corrected, since it does not
  !   ! contain discontinuities. The correction procedure is done for correct computation
  !   ! of the viscous terms (and advective terms), so it should be called before computing
  !   ! them. Also, the wall-normal gradient is represented using the first off-wall and
  !   ! ghost points, so the correction procedure should be called before bounduvw.
  !   ! The correction is done for both the prediction and projection steps, because we assume
  !   ! that the first off-wall point is always incorrect, so should be always replaced by 
  !   ! the extrapolation. The correction may not be necessary for the intermidiate velocity.
  !   ! The correction makes the velocity near the first off-wall point is not divergence-free.
  !   !
  !   !
  !   implicit none
  !   integer , intent(in), dimension(3) :: n
  !   logical , intent(in), dimension(0:1,3) :: is_bound
  !   integer , intent(in), dimension(0:1,3) :: lwm
  !   real(rp), intent(in), dimension(3)  :: dl
  !   real(rp), intent(in), dimension(0:) :: dzc,dzf
  !   real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
  !   real(rp), dimension(0:n(3)+1) :: dzci,dzfi
  !   real(rp) :: d(2:4),dd
  !   integer  :: i,j,k
  !   !
  !   if(is_bound(0,3).and.lwm(0,3)/=0) then
  !     d(2) = 0.5_rp*dzf(2)
  !     d(3) = d(2) + dzc(2)
  !     d(4) = d(3) + dzc(3)
  !     do j = 1,n(2)
  !       do i = 1,n(1)
  !         dd = deriv_1st_ord(u(i,j,2:3),d(2:3),ibound=0)
  !         ! dd = deriv_2nd_ord(u(i,j,2:4),d(2:4),ibound=0)
  !         u(i,j,1) = u(i,j,2) - dd*dzc(1)
  !         dd = deriv_1st_ord(v(i,j,2:3),d(2:3),ibound=0)
  !         ! dd = deriv_2nd_ord(v(i,j,2:4),d(2:4),ibound=0)
  !         v(i,j,1) = v(i,j,2) - dd*dzc(1)
  !       end do
  !     end do
  !   end if
  !   !
  !   if(is_bound(1,3).and.lwm(1,3)/=0) then
  !     d(2) = 0.5_rp*dzf(n(3)-1)
  !     d(3) = d(2) + dzc(n(3)-2)
  !     d(4) = d(3) + dzc(n(3)-3)
  !     do j = 1,n(2)
  !       do i = 1,n(1)
  !         dd = deriv_1st_ord(u(i,j,n(3)-1:n(3)-2:-1),d(2:3),ibound=1)
  !         ! dd = deriv_2nd_ord(u(i,j,n(3)-1:n(3)-3:-1),d(2:4),ibound=1)
  !         u(i,j,n(3)) = u(i,j,n(3)-1) + dd*dzc(n(3)-1)
  !         dd = deriv_1st_ord(v(i,j,n(3)-1:n(3)-2:-1),d(2:3),ibound=1)
  !         ! dd = deriv_2nd_ord(v(i,j,n(3)-1:n(3)-3:-1),d(2:4),ibound=1)
  !         v(i,j,n(3)) = v(i,j,n(3)-1) + dd*dzc(n(3)-1)
  !       end do
  !     end do
  !   end if
  ! end subroutine correc_1st_point
end module mod_wmodel