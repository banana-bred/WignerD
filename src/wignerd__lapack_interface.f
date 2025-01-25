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
