! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_precision
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dl,dzci,dzfi,visc,visct,u,v,w,dtmax)
    !
    ! computes maximum allowed time step, refer to Pieter Wesseling (P200, 2009)
    ! for the stability conditions of the advective and diffusion terms
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: visct,u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti
    real(rp) :: dtidx,dtidy,dtidz,dtid,dl2i
    real(rp) :: viscx,viscy,viscz
    integer :: i,j,k
    !
    dxi = 1._rp/dl(1)
    dyi = 1._rp/dl(2)
    dzi = 1._rp/dl(3)
    !
    dti = 0._rp
    dtid = 0._rp
    !$acc data copy(dti) async(1)
    !$acc parallel loop collapse(3) default(present) private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) reduction(max:dti) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) REDUCTION(max:dti)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti  = max(dti,dtix,dtiy,dtiz)
          !
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
          ! nothing to do
#else
          viscx = visc + 0.5_rp*(visct(i,j,k)+visct(i+1,j  ,k  ))
          viscy = visc + 0.5_rp*(visct(i,j,k)+visct(i  ,j+1,k  ))
          viscz = visc + 0.5_rp*(visct(i,j,k)+visct(i  ,j  ,k+1))
          dl2i  = dxi*dxi+dyi*dyi
          dtidx = dl2i
          dtidy = dl2i
          dtidz = dl2i
#if !defined(_IMPDIFF_1D)
          dtidx = dtidx + dzfi(k)*dzfi(k)
          dtidy = dtidy + dzfi(k)*dzfi(k)
          dtidz = dtidz + dzci(k)*dzci(k)
#endif
          dtidx = viscx*dtidx
          dtidy = viscy*dtidy
          dtidz = viscz*dtidz
          dtid  = max(dtid,dtidx,dtidy,dtidz)
#endif
        end do
      end do
    end do
    if(dti == 0._rp) dti = 1._rp
    #if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
        dtmax = sqrt(3._rp)/dti
    #else
        dtmax = min(0.4125_rp/dtid,sqrt(3._rp)/dti)
        !the viscous CFL might be as high as 1.5
    #endif
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,dtmax,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  end subroutine chkdt
end module mod_chkdt
