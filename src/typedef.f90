! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_typedef
  use mod_precision
  type cond_bound
  real(rp), allocatable, dimension(:,:,:) :: x
  real(rp), allocatable, dimension(:,:,:) :: y
  real(rp), allocatable, dimension(:,:,:) :: z
  end type cond_bound
end module mod_typedef
