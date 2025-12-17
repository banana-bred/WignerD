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
