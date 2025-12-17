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
