module mod_sgs
  implicit none
  private
  public :: mySubroutine
  contains
    subroutine mySubroutine(c)
      integer, dimension(0:), intent(in) :: c
      print*, 1
    end subroutine mySubroutine
  end module mod_sgs

  program main
  use mod_sgs

  integer, allocatable :: c(:)
  integer :: b(5),d(5,10)

  ! allocate(c(0:10))
  ! c = 1
  call mySubroutine(c)

  b = 0
  d = 0

  print*, b*d

end program main