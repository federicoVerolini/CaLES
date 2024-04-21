program main
  use mod_precision
  use mod_wmodel    , only: wallmodel
  implicit none
  real(rp) :: tauw(2)
  integer :: mtype
  real(rp) :: reb, uh, vh, h, l1d, visci, visc, tauw_ini, kap, b

  kap = 0.41_rp
  b = 5.2_rp

  reb = 250000._rp

  mtype = 1
  uh = 0._rp
  vh = 0._rp
  h = 0.1_rp
  l1d = 2._rp
  visci = (reb/2._rp)
  visc = 1._rp/visci

  tauw_ini = (visc/h*exp(-kap*b))**2
  print*, tauw_ini
  tauw_ini = epsilon(1._rp)
  tauw_ini = 0._rp
  call wallmodel(mtype,uh,vh,h,l1d,visc,tauw_ini,tauw)

  ! better to initialize with small values

end program main