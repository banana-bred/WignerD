! ================================================================================================================================ !
module wignerd__characters
  !! Contains procedures related to characters and character arrays,
  !! such as converting integers to characters

  use wignerd__types, only: rp

  implicit none

  private

  ! -- procedures
  public :: ndigits
  public :: int2char

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  recursive function int2char(i) result(output)
    !! Convert the integer i to a character (array) "i"

    implicit none

    integer, intent(in) :: i
    character(:), allocatable :: output
      !! The integer "i" as a character array, e.g.
      !! int2char(2)  \(\to\) "2"
      !! int2char(16) \(\to\) "16"

    integer :: n

    n = ndigits(i)

    ! -- allocate characters so that they're wide enough to write to
    allocate(character(n) :: output)

    write(output, '(I0)') i

  end function int2char

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function ndigits(n) result(num)
    !! Returns number of characters an integer will occupy

    use wignerd__constants, only: one

    implicit none

    integer, intent(in) :: n
    integer :: num
      !! The number of characters it takes to represent a number, e.g.
      !! ndigits(7)    \(to\) 1
      !! ndigits(-7)   \(to\) 2
      !! ndigits(38)   \(to\) 2
      !! ndigits(3877) \(to\) 7

    num = 1

    if(n .eq. 0) return

    num = floor(log10(abs(n) * one)) + 1

    ! -- account for minus sign
    if(n.lt.1) num = num + 1

  end function ndigits

! ================================================================================================================================ !
end module wignerd__characters
! ================================================================================================================================ !
