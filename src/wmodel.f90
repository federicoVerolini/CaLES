! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_wm
  use mod_const
  use mod_typedef, only: cond_bound
  implicit none
  private
  public comput_bcuvw,comput_bcp
  contains
      !
  subroutine comput_bcuvw(cbc,n,bc,is_bound,u,v,w,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw, determined via bcvel or wall model
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(cond_bound), intent(inout) :: bcu,bcv,bcw
    !
    if(.false.) then
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
        if(cbc(0,3,1)/='P') bcu%z(:,:,0) = -1.*u(:,:,1   )
        if(cbc(0,3,2)/='P') bcv%z(:,:,0) = -v(:,:,1   )
        if(cbc(0,3,3)/='P') bcw%z(:,:,0) =  w(:,:,1   )
      end if
      if(is_bound(1,3)) then
        if(cbc(1,3,1)/='P') bcu%z(:,:,1) = -1.*u(:,:,n(3))
        if(cbc(1,3,2)/='P') bcv%z(:,:,1) = -v(:,:,n(3))
        if(cbc(1,3,3)/='P') bcw%z(:,:,1) =  w(:,:,n(3))
      end if
    end if
    !
  end subroutine comput_bcuvw
  !
  subroutine comput_bcp(cbc,n,bc,is_bound,bcp)
    !
    ! bcp, determined via bcpre
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3) :: bc
    logical , intent(in), dimension(0:1,3) :: is_bound
    type(cond_bound), intent(inout) :: bcp
    !
    if(is_bound(0,1).and.cbc(0,1)/='P') bcp%x(:,:,0) = bc(0,1)
    if(is_bound(1,1).and.cbc(1,1)/='P') bcp%x(:,:,1) = bc(1,1)
    if(is_bound(0,2).and.cbc(0,2)/='P') bcp%y(:,:,0) = bc(0,2)
    if(is_bound(1,2).and.cbc(1,2)/='P') bcp%y(:,:,1) = bc(1,2)
    if(is_bound(0,3).and.cbc(0,3)/='P') bcp%z(:,:,0) = bc(0,3)
    if(is_bound(1,3).and.cbc(1,3)/='P') bcp%z(:,:,1) = bc(1,3)
    !
  end subroutine comput_bcp
end module mod_wm