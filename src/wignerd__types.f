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
module wignerd__types
  !! Module containing defining the precision for real and complex types

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private

  ! -- types
  integer, parameter, public :: rp = real64
    !! The precision for real and complex types.
  real(rp), parameter, public :: rp_huge = huge(1.0_rp)
    !! The largest representable real number with precision `real(rp)`

! ================================================================================================================================ !
end module wignerd__types
! ================================================================================================================================ !
