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
  subroutine comput_bcuvw(cbc,n,bc,bcu,bcv,bcw)
    !
    ! bcu,bcv,bcw, determined via bcvel or wall model
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    real(rp), intent(inout), dimension(0:,0:,0:) :: bcu,bcv,bcw
    !
    if(.true.) then
      bcu(:,:,0) = bc(0,3,1)
      bcv(:,:,0) = bc(0,3,2)
      bcw(:,:,0) = bc(0,3,3)
      bcu(:,:,1) = bc(1,3,1)
      bcv(:,:,1) = bc(1,3,2)
      bcw(:,:,1) = bc(1,3,3)
    else !%to be replaced by wall model
      !
    end if
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
    type(cond_bound), intent(out) :: bcp
    !
    if(is_bound(0,1)) bcp%x(:,:,0) = bc(0,1)
    if(is_bound(1,1)) bcp%x(:,:,1) = bc(1,1)
    if(is_bound(0,2)) bcp%y(:,:,0) = bc(0,2)
    if(is_bound(1,2)) bcp%y(:,:,1) = bc(1,2)
    if(is_bound(0,3)) bcp%z(:,:,0) = bc(0,3)
    if(is_bound(1,3)) bcp%z(:,:,1) = bc(1,3)
    !
  end subroutine comput_bcp
end module mod_wm


!   subroutine bounduvw(cbc,n,bc,nb,is_bound,is_correc,dl,dzc,dzf,u,v,w)
!     !
!     ! imposes velocity boundary conditions
!     !
!     implicit none
!     character(len=1), intent(in), dimension(0:1,3,3) :: cbc
!     integer , intent(in), dimension(3) :: n
!     real(rp), intent(in), dimension(0:1,3,3) :: bc
!     integer , intent(in), dimension(0:1,3  ) :: nb
!     logical , intent(in), dimension(0:1,3  ) :: is_bound
!     logical , intent(in)                     :: is_correc
!     real(rp), intent(in), dimension(3 ) :: dl
!     real(rp), intent(in), dimension(0:) :: dzc,dzf
!     real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
!     logical :: impose_norm_bc
!     integer :: idir,nh
!     !
!     nh = 1 !number of ghost points
!     !
! #if !defined(_OPENACC)
!     do idir = 1,3
!       call updthalo(nh,halo(idir),nb(:,idir),idir,u) !info from neighboring block (MPI data exchange)
!       call updthalo(nh,halo(idir),nb(:,idir),idir,v)
!       call updthalo(nh,halo(idir),nb(:,idir),idir,w)
!     end do
! #else
!     call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u)
!     call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v)
!     call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w)
! #endif
!     !
!     impose_norm_bc = (.not.is_correc).or.(cbc(0,1,1)//cbc(1,1,1) == 'PP')
!     if(is_bound(0,1)) then
!       if(impose_norm_bc) call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),dl(1),u) !periodic for channel, unchanged
!                          call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),dl(1),v) !periodic for channel, unchanged
!                          call set_bc(cbc(0,1,3),0,1,nh,.true. ,bc(0,1,3),dl(1),w) !periodic for channel, unchanged
!     end if
!     if(is_bound(1,1)) then
!       if(impose_norm_bc) call set_bc(cbc(1,1,1),1,1,nh,.false.,bc(1,1,1),dl(1),u) !periodic for channel, unchanged
!                          call set_bc(cbc(1,1,2),1,1,nh,.true. ,bc(1,1,2),dl(1),v) !periodic for channel, unchanged
!                          call set_bc(cbc(1,1,3),1,1,nh,.true. ,bc(1,1,3),dl(1),w) !periodic for channel, unchanged
!     end if
!     impose_norm_bc = (.not.is_correc).or.(cbc(0,2,2)//cbc(1,2,2) == 'PP')
!     if(is_bound(0,2)) then                                                        
!                          call set_bc(cbc(0,2,1),0,2,nh,.true. ,bc(0,2,1),dl(2),u) !unused in channel
!       if(impose_norm_bc) call set_bc(cbc(0,2,2),0,2,nh,.false.,bc(0,2,2),dl(2),v) !unused in channel
!                          call set_bc(cbc(0,2,3),0,2,nh,.true. ,bc(0,2,3),dl(2),w) !unused in channel
!      end if
!     if(is_bound(1,2)) then                                                        
!                          call set_bc(cbc(1,2,1),1,2,nh,.true. ,bc(1,2,1),dl(2),u) !unused in channel
!       if(impose_norm_bc) call set_bc(cbc(1,2,2),1,2,nh,.false.,bc(1,2,2),dl(2),v) !unused in channel
!                          call set_bc(cbc(1,2,3),1,2,nh,.true. ,bc(1,2,3),dl(2),w) !unused in channel
!     end if
!     impose_norm_bc = (.not.is_correc).or.(cbc(0,3,3)//cbc(1,3,3) == 'PP')
!     if(is_bound(0,3)) then                                                            !dzc, between centers; dzf, between faces
!                          call set_bc(cbc(0,3,1),0,3,nh,.true. ,bc(0,3,1),dzc(0)   ,u) !Dirichlet, u(0) = 2*0-u(1); Neumann, u(0) = -dzc*0 + u(1)
!                          call set_bc(cbc(0,3,2),0,3,nh,.true. ,bc(0,3,2),dzc(0)   ,v) !Dirichlet, v(0) = 2*0-v(1); Neumann, v(0) = -dzc*0 + v(1)
!       if(impose_norm_bc) call set_bc(cbc(0,3,3),0,3,nh,.false.,bc(0,3,3),dzf(0)   ,w) !Dirichlet, w(0) = 0, no w(-1);
!     end if
!     if(is_bound(1,3)) then
!                          call set_bc(cbc(1,3,1),1,3,nh,.true. ,bc(1,3,1),dzc(n(3)),u) !Dirichlet, u(n+1) = 2*0-u(n); Neumann, u(n+1) = dzc*0 + u(n)
!                          call set_bc(cbc(1,3,2),1,3,nh,.true. ,bc(1,3,2),dzc(n(3)),v) !Dirichlet, v(n+1) = 2*0-v(n); Neumann, u(n+1) = dzc*0 + u(n)
!       if(impose_norm_bc) call set_bc(cbc(1,3,3),1,3,nh,.false.,bc(1,3,3),dzf(n(3)),w) !Dirichlet, w(n) = 0, w(n+1) = w(n-1);
!     end if
!   end subroutine bounduvw
