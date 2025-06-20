! WignerD: calculate the Wigner matrices d and D
! Copyright (C) 2025 Josh Forer <j.forer@posteo.net>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
