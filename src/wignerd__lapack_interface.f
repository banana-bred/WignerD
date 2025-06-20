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
module wignerd__lapack_interface
  !! Select at build time which lapack implementation to use

  implicit none

  private

  public :: zhbev

  interface zhbev
    !! This interface exists in case the user wants to use their own
    !! implementation of the LAPACK routine `ZHBEV` and not the one
    !! found in the standard library.
    subroutine zhbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info)
      use, intrinsic  :: iso_fortran_env, only: real64
      character       :: jobz
      character       :: uplo
      integer         :: n
      integer         :: kd
      complex(real64) :: ab(ldab, *)
      integer         :: ldab
      real(real64)    :: w(*)
      complex(real64) :: z(ldz, *)
      integer         :: ldz
      complex(real64) :: work(*)
      real(real64)    :: rwork(*)
      integer         :: info
    end subroutine
  end interface zhbev

! ================================================================================================================================ !
end module wignerd__lapack_interface
! ================================================================================================================================ !
