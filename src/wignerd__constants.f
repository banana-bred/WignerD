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
module wignerd__constants
  !! Contains defined constants

  use wignerd__types, only: rp

  implicit none

  private

  ! -- numbers
  real(rp), parameter, public :: ahalf = 0.5_rp
  real(rp), parameter, public :: zero  = 0._rp
  real(rp), parameter, public :: one   = 1._rp
  real(rp), parameter, public :: two   = 2._rp
  real(rp), parameter, public :: three = 3._rp
  real(rp), parameter, public :: four  = 4._rp
  real(rp), parameter, public :: pi    = atan(1._rp) * 4._rp
    !! π
  complex(rp), parameter, public :: im = (zero, one)
    !! the square root of -1

! ================================================================================================================================ !
end module wignerd__constants
! ================================================================================================================================ !
