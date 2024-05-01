! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_sgs
  use mod_precision
  use mod_param, only: c_smag,big
  use mod_post, only: strain_rate
  implicit none
  private
  public cmpt_sgs
  contains
    !
  subroutine cmpt_sgs(sgstype,cbc,n,is_bound,l,dl,zc,dzc,dzf,visc,u,v,w,visct)
    !
    ! compute subgrid viscosity at cell centers
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(0:)       :: zc,dzc,dzf
    real(rp), intent(in )                      :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: visct
    real(rp), allocatable, dimension(:,:,:)    :: str,dw_plus
    real(rp), dimension(3)        :: dli
    real(rp), dimension(0:n(3)+1) :: dzci,dzfi
    !
    select case(trim(sgstype))
    case('none')
      visct(:,:,:) = 0._rp
    case('smag')
      dli(:)  = dl( :)**(-1)
      dzci(:) = dzc(:)**(-1)
      dzfi(:) = dzf(:)**(-1)
      allocate(str    (1:n(1),1:n(2),1:n(3)), &
               dw_plus(1:n(1),1:n(2),1:n(3)))
      call strain_rate(n,dli,dzci,dzfi,u,v,w,str)
      call cmpt_dw_plus(cbc,n,is_bound,l,dl,zc,dzc,visc,u,v,w,dw_plus)
      call sgs_smag(n,dl,dzf,str,dw_plus,visct)
      deallocate(str, &
                 dw_plus)
    case('amd')
      print*, 'ERROR: AMD model not implemented yet'
    case default
      print*, 'ERROR: unknown SGS model'
    end select
  end subroutine cmpt_sgs
  !
  subroutine sgs_smag(n,dl,dzf,str,dw_plus,visct)
    !
    ! classical Smagorinsky model with van Driest damping
    ! 
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dl
    real(rp), intent(in ), dimension(0:)       :: dzf
    real(rp), intent(in ), dimension(1:,1:,1:) :: str,dw_plus
    real(rp), intent(out), dimension(0:,0:,0:) :: visct
    real(rp) :: del,fd
    integer :: i,j,k
    !
    do k=1,n(3)
      del = (dl(1)*dl(2)*dzf(k))**(1./3.)
      do j=1,n(2)
        do i=1,n(1)
          fd = 1._rp-exp(-dw_plus(i,j,k)/25._rp)
          visct(i,j,k) = (c_smag*del*fd)**2*sqrt(2._rp*str(i,j,k))
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
    ! identification of walls is based on boundary conditions, which might
    ! be problematic in some cases.
    !
    ! It is inappropriate to assume no-slip wall bc. We must use the ghost
    ! cells to compute wall friction.
    ! 
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(0:)       :: zc,dzc
    real(rp), intent(in )                      :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(1:,1:,1:) :: dw_plus
    real(rp), dimension(:,:,:), allocatable    :: dw
    real(rp) :: tauw(2),tauw_tot,this_dw_plus,this_dw,visci
    real(rp) :: u_cci,u_cco,u_mci,u_mco,v_cci,v_cco,v_cmi,v_cmo, &
                u_cic,u_coc,u_mic,u_moc,w_cic,w_coc,w_cim,w_com, &
                v_icc,v_occ,v_imc,v_omc,w_icc,w_occ,w_icm,w_ocm
    integer :: i,j,k
    !
    visci = 1._rp/visc
    !
    allocate(dw(1:n(1),1:n(2),1:n(3)))
    dw = big
    dw_plus = big
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
      do k=1,n(3)
        do j=1,n(2)
          tauw(1) = visc*0.5_rp*(v(n(1),j,k)-v(n(1)+1,j,k)+v(n(1),j-1,k)-v(n(1)+1,j-1,k))/dl(1)
          tauw(2) = visc*0.5_rp*(w(n(1),j,k)-w(n(1)+1,j,k)+w(n(1),j,k-1)-w(n(1)+1,j,k-1))/dl(1)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
      do k=1,n(3)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,n(2),k)-u(i,n(2)+1,k)+u(i-1,n(2),k)-u(i-1,n(2)+1,k))/dl(2)
          tauw(2) = visc*0.5_rp*(w(i,n(2),k)-w(i,n(2)+1,k)+w(i,n(2),k-1)-w(i,n(2)+1,k-1))/dl(2)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,1)-u(i,j,0)+u(i-1,j,1)-u(i-1,j,0))/dzc(0)
          tauw(2) = visc*0.5_rp*(v(i,j,1)-v(i,j,0)+v(i,j-1,1)-v(i,j-1,0))/dzc(0)
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
      do j=1,n(2)
        do i=1,n(1)
          tauw(1) = visc*0.5_rp*(u(i,j,n(3))-u(i,j,n(3)+1)+u(i-1,j,n(3))-u(i-1,j,n(3)+1))/dzc(n(3))
          tauw(2) = visc*0.5_rp*(v(i,j,n(3))-v(i,j,n(3)+1)+v(i,j-1,n(3))-v(i,j-1,n(3)+1))/dzc(n(3))
          tauw_tot= sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))
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
end module mod_sgs