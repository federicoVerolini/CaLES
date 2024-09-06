! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_typedef
  use mod_precision
  type bound
  real(rp), allocatable, dimension(:,:,:) :: x
  real(rp), allocatable, dimension(:,:,:) :: y
  real(rp), allocatable, dimension(:,:,:) :: z
  end type bound
end module mod_typedef