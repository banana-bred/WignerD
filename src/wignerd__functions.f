! WignerD: calculate the Wigner d and D matrices
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
module wignerd__functions
  !! Contains various functions used in the code

  implicit none

  private

  public delta
  public factorial

  interface factorial
    !! Returns \(n!\) for integer or double precision real \(n\)
    module procedure :: factorial_int
    module procedure :: factorial_real
  end interface factorial

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function delta(m, n) result(res)
    !! Return the Kronecker delta function \(δ_{m,n}\)
    use wignerd__types,     only: rp
    use wignerd__constants, only: zero, one
    implicit none
    integer, intent(in) :: m, n
    real(rp) :: res
    res = zero ; if(m .eq. n) res = one
  end function delta

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function factorial_int(n) result(res)
    !! Returns \(n!\) for integer \(n\)
    use wignerd__types, only: rp
    implicit none
    integer, intent(in) :: n
    real(rp) :: res
    res = gamma(real(n + 1, kind = rp))
  end function factorial_int

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function factorial_real(n) result(res)
    !! Returns \(n!\) for double precision real \(n\)
    use wignerd__types,     only: rp
    use wignerd__constants, only: one
    implicit none
    real(rp), intent(in) :: n
    real(rp) :: res
    res = gamma(n + one)
  end function factorial_real

! ================================================================================================================================ !
end module wignerd__functions
! ================================================================================================================================ !
