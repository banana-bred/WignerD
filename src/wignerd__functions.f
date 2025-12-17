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
