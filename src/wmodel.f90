! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_wmodel
  use mod_precision
  use mod_typedef, only: cond_bound
  implicit none
  private
  public cmpt_bcuvw,cmpt_bcp,cmpt_bctau
  contains
    !
  subroutine cmpt_bcuvw(cbc,n,bc,is_bound,is_wm,visc,u,v,w,bctau1,bctau2,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw, determined via bcvel or wall model
    !
    implicit none
    character(len=1), intent(inout), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in), dimension(0:1,3) :: is_bound
    logical , intent(in) :: is_wm
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(cond_bound), intent(in) :: bctau1,bctau2
    type(cond_bound), intent(inout) :: bcu,bcv,bcw
    !
    if(.not.is_wm) then
      !bcu,bcv,bcw, determined via bcvel
      if(is_bound(0,1)) then
        if(cbc(0,1,1)/='P') bcu%x(:,:,0) = bc(0,1,1)
        if(cbc(0,1,2)/='P') bcv%x(:,:,0) = bc(0,1,2)
        if(cbc(0,1,3)/='P') bcw%x(:,:,0) = bc(0,1,3)
      end if
      if(is_bound(1,1)) then
        if(cbc(1,1,1)/='P') bcu%x(:,:,1) = bc(1,1,1)
        if(cbc(1,1,2)/='P') bcv%x(:,:,1) = bc(1,1,2)
        if(cbc(1,1,3)/='P') bcw%x(:,:,1) = bc(1,1,3)
      end if
      if(is_bound(0,2)) then
        if(cbc(0,2,1)/='P') bcu%y(:,:,0) = bc(0,2,1)
        if(cbc(0,2,2)/='P') bcv%y(:,:,0) = bc(0,2,2)
        if(cbc(0,2,3)/='P') bcw%y(:,:,0) = bc(0,2,3)
      end if
      if(is_bound(1,2)) then
        if(cbc(1,2,1)/='P') bcu%y(:,:,1) = bc(1,2,1)
        if(cbc(1,2,2)/='P') bcv%y(:,:,1) = bc(1,2,2)
        if(cbc(1,2,3)/='P') bcw%y(:,:,1) = bc(1,2,3)
      end if
      if(is_bound(0,3)) then
        if(cbc(0,3,1)/='P') bcu%z(:,:,0) = bc(0,3,1)
        if(cbc(0,3,2)/='P') bcv%z(:,:,0) = bc(0,3,2)
        if(cbc(0,3,3)/='P') bcw%z(:,:,0) = bc(0,3,3)
      end if
      if(is_bound(1,3)) then
        if(cbc(1,3,1)/='P') bcu%z(:,:,1) = bc(1,3,1)
        if(cbc(1,3,2)/='P') bcv%z(:,:,1) = bc(1,3,2)
        if(cbc(1,3,3)/='P') bcw%z(:,:,1) = bc(1,3,3)
      end if
    else
      if(is_bound(0,3)) then
        cbc(0,3,1) = 'N'
        cbc(0,3,2) = 'N'
        if(cbc(0,3,1)/='P') bcu%z(1:n(1),1:n(2),0) =  0.5_rp/visc*(bctau1%z(1:n(1),1:n(2),0) + bctau1%z(2:n(1)+1,1:n(2)  ,0))
        if(cbc(0,3,2)/='P') bcv%z(1:n(1),1:n(2),0) =  0.5_rp/visc*(bctau2%z(1:n(1),1:n(2),0) + bctau2%z(1:n(1)  ,2:n(2)+1,0))
      end if
      if(is_bound(1,3)) then
        cbc(1,3,1) = 'N'
        cbc(1,3,2) = 'N'
        if(cbc(1,3,1)/='P') bcu%z(1:n(1),1:n(2),1) = -0.5_rp/visc*(bctau1%z(1:n(1),1:n(2),1) + bctau1%z(2:n(1)+1,1:n(2)  ,1))
        if(cbc(1,3,2)/='P') bcv%z(1:n(1),1:n(2),1) = -0.5_rp/visc*(bctau2%z(1:n(1),1:n(2),1) + bctau2%z(1:n(1)  ,2:n(2)+1,1))
      end if
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
  subroutine cmpt_bctau(cbc,n,bc,is_bound,zc,visc,kap,b,h,u,v,w,bctau1,bctau2)
    !
    ! compute bctau at cell centers via wall models
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: visc,kap,b,h
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(cond_bound), intent(inout) :: bctau1,bctau2
    real(rp) :: wei,uc,vc,tauw(2)
    integer  :: i,j,k,k1,k2
    !lower wall
    k = 1
    do while(zc(k) < h)
      k = k + 1
    end do
    k2 = k
    k1 = k - 1
    wei= (h-zc(k1))/(zc(k2)-zc(k1))
    do j = 1,n(2)+1
      do i = 1,n(1)+1
        uc = (1._rp - wei)*(u(i-1,j  ,k1) + u(i,j,k1)) &
                    + wei *(u(i-1,j  ,k2) + u(i,j,k2))
        vc = (1._rp - wei)*(v(i  ,j-1,k1) + v(i,j,k1)) &
                    + wei *(v(i  ,j-1,k2) + v(i,j,k2))
        uc = uc*0.5_rp
        vc = vc*0.5_rp
        call wallmodel(uc,vc,h,visc,kap,b,tauw)
        bctau1%z(i,j,0) = tauw(1)
        bctau2%z(i,j,0) = tauw(2)
      end do
    end do
    !upper wall
    k = n(3)
    do while(zc(k) > 1._rp-h)
      k = k - 1
    end do
    k2 = k
    k1 = k + 1
    wei= ((1._rp-h)-zc(k1))/(zc(k2)-zc(k1))
    do j = 1,n(2)+1
      do i = 1,n(1)+1
        uc = (1._rp - wei)*(u(i-1,j  ,k1) + u(i,j,k1)) &
                    + wei *(u(i-1,j  ,k2) + u(i,j,k2))
        vc = (1._rp - wei)*(v(i  ,j-1,k1) + v(i,j,k1)) &
                    + wei *(v(i  ,j-1,k2) + v(i,j,k2))
        uc = uc*0.5_rp
        vc = vc*0.5_rp
        call wallmodel(uc,vc,h,visc,kap,b,tauw)
        bctau1%z(i,j,1) = tauw(1)
        bctau2%z(i,j,1) = tauw(2)
      end do
    end do
  end subroutine cmpt_bctau
  !
  subroutine wallmodel(uh,vh,h,visc,kap,b,tauw)
    real(rp), intent(in)  :: uh,vh,h,visc,kap,b
    real(rp), intent(out), dimension(2) :: tauw
    real(rp) :: upar,utau,f,fp,conv,tauw_old,tauw_tot
    real(rp) :: umax
    !
    if(.false.) then
      !log-law profile
      upar = sqrt(uh*uh+vh*vh)
      utau = sqrt(upar/h*visc)
      conv = 1._rp
      do while(conv > 1.e-4)
        tauw_old = utau*utau
        f  = upar/utau - 1._rp/kap*log(h*utau/visc) - b
        fp = -1._rp/utau*(upar/utau + 1._rp/kap)
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