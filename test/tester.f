! ================================================================================================================================ !
module test_wignerd
  use wignerd,            only: wigner_d, wigner_little_d, wigner_big_D
  use testdrive,          only: new_unittest, unittest_type, error_type, check
  use wignerd__types,     only: rp
  use wignerd__constants, only: zero, one, two, three, four, pi

  implicit none

  private

  public :: collect_wignerd

  real(rp), parameter :: atol = 1e-14_rp

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine collect_wignerd(testsuites)
    !! Collect all exported unit tests

    implicit none

    type(unittest_type), allocatable, intent(out) :: testsuites(:)
      !! collection of tests

    testsuites = [ &
                new_unittest("d_analytic special values", test_little_d_special_values_analytic) &
              , new_unittest("d_diag special values", test_little_d_special_values_diag) &
      ]
  end subroutine collect_wignerd

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_little_d_special_values_analytic(error)
    !! Test known analytic values of $d^j_{m',m}(\beta)$ for small j using the analytic d-matrix routine

    implicit none

    type(error_type), allocatable, intent(out) :: error

    integer  :: ibeta
    real(rp) :: j
    real(rp), target :: beta_array(8)
    real(rp), allocatable :: d(:,:)
    real(rp), pointer :: beta

    ! -- the four quadrants and the four axes
    beta_array = [(pi * ibeta / four, ibeta = 0, 7 )]

    ! -- j = 0 --------------------------------------------------- !
    j    = zero
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      ! -- integer j
      d = wigner_little_d(nint(j), beta, use_analytic = .true.)

      call check(error, d(1,1), one, thr = atol)
      if(allocated(error)) return

      ! -- real j
      d = wigner_little_d(j, beta, use_analytic = .true.)

      call check(error, d(1,1), one)
      if(allocated(error)) return

    enddo

    ! -- j = 1/2 ------------------------------------------------- !
    j    = one / two
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      d    = wigner_little_d(j, beta, use_analytic = .true.)

      ! -- m' = 1/2, m = -1/2
      call check(error, d(2,1), -sin(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m =  1/2
      call check(error, d(2,2), cos(beta / 2), thr = atol)
      if(allocated(error)) return

    enddo

    ! -- j = 1 --------------------------------------------------- !
    j = one
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      ! -- integer j
      d    = wigner_little_d(nint(j), beta, use_analytic = .true.)

      ! -- m' = 0, m = 0
      call check(error, d(2,2), cos(beta), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = -1
      call check(error, d(3,1), (one - cos(beta)) / two, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 0
      call check(error, d(3,2), -sin(beta) / sqrt(two), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 1
      call check(error, d(3,3), (one + cos(beta)) / two, thr = atol)
      if(allocated(error)) return

      ! -- real j
      d    = wigner_little_d(j, beta, use_analytic = .true.)

      ! -- m' = 0, m = 0
      call check(error, d(2,2), cos(beta), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = -1
      call check(error, d(3,1), (one - cos(beta)) / two, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 0
      call check(error, d(3,2), -sin(beta) / sqrt(two), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 1
      call check(error, d(3,3), (one + cos(beta)) / two, thr = atol)
      if(allocated(error)) return

    enddo

    ! -- j = 3/2 ------------------------------------------------- !
    j    = three / two
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      d    = wigner_little_d(j, beta, use_analytic = .true.)

      ! -- m' = 3/2, m = 3/2
      call check(error, d(4,4), cos(beta / 2) * (1 + cos(beta)) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = 1/2
      call check(error, d(4,3), -sqrt(three)/2 * (1 + cos(beta)) * sin(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = -1/2
      call check(error, d(4,2), sqrt(three)/2 * (1 - cos(beta)) * cos(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = -3/2
      call check(error, d(4,1), -(1 - cos(beta)) * sin(beta / 2) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m = 1/2
      call check(error, d(3,3), (3*cos(beta) - 1) * cos(beta / 2) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m = -1/2
      call check(error, d(3,2), -(3*cos(beta) + 1) * sin(beta / 2) / 2, thr = atol)
      if(allocated(error)) return

    enddo

    nullify(beta)


  end subroutine test_little_d_special_values_analytic

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_little_d_special_values_diag(error)
    !! Test known analytic values of $d^j_{m',m}(\beta)$ for small j using the diagonalization d-matrix routine

    type(error_type), allocatable, intent(out) :: error

    integer  :: ibeta
    real(rp) :: j
    real(rp), target :: beta_array(8)
    real(rp), allocatable :: d(:,:)
    real(rp), pointer :: beta

    ! -- the four quadrants and the four axes
    beta_array = [(pi * ibeta / four, ibeta = 0, 7 )]

    ! -- j = 0 --------------------------------------------------- !
    j = zero
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      ! -- integer j
      d    = wigner_little_d(nint(j), beta)

      call check(error, d(1,1), one, thr = atol)
      if(allocated(error)) return

      ! -- real j
      d = wigner_little_d(j, beta)

      call check(error, d(1,1), one, thr = atol)
      if(allocated(error)) return

    enddo

    ! -- j = 1/2 ------------------------------------------------- !
    j    = one / two
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      d    = wigner_little_d(j, beta)

      ! -- m' = 1/2, m = -1/2
      call check(error, d(2,1), -sin(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m =  1/2
      call check(error, d(2,2), cos(beta / 2), thr = atol)
      if(allocated(error)) return

    enddo

    ! -- j = 1 --------------------------------------------------- !
    j = one
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      ! -- integer j
      d    = wigner_little_d(nint(j), beta)

      ! -- m' = 0, m = 0
      call check(error, d(2,2), cos(beta), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = -1
      call check(error, d(3,1), (one - cos(beta)) / two, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 0
      call check(error, d(3,2), -sin(beta) / sqrt(two), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 1
      call check(error, d(3,3), (one + cos(beta)) / two, thr = atol)
      if(allocated(error)) return

      ! -- real j
      d    = wigner_little_d(j, beta)

      ! -- m' = 0, m = 0
      call check(error, d(2,2), cos(beta), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = -1
      call check(error, d(3,1), (one - cos(beta)) / two, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 0
      call check(error, d(3,2), -sin(beta) / sqrt(two), thr = atol)
      if(allocated(error)) return
      ! -- m' = 1, m = 1
      call check(error, d(3,3), (one + cos(beta)) / two, thr = atol)
      if(allocated(error)) return

    enddo


    ! -- j = 3/2 ------------------------------------------------- !
    j    = three / two
    do ibeta = 1, size(beta_array, 1)

      beta => beta_array(ibeta)

      d    = wigner_little_d(j, beta)

      ! -- m' = 3/2, m = 3/2
      call check(error, d(4,4), cos(beta / 2) * (1 + cos(beta)) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = 1/2
      call check(error, d(4,3), -sqrt(three)/2 * (1 + cos(beta)) * sin(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = -1/2
      call check(error, d(4,2), sqrt(three)/2 * (1 - cos(beta)) * cos(beta / 2), thr = atol)
      if(allocated(error)) return
      ! -- m' = 3/2, m = -3/2
      call check(error, d(4,1), -(1 - cos(beta)) * sin(beta / 2) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m = 1/2
      call check(error, d(3,3), (3*cos(beta) - 1) * cos(beta / 2) / 2, thr = atol)
      if(allocated(error)) return
      ! -- m' = 1/2, m = -1/2
      call check(error, d(3,2), -(3*cos(beta) + 1) * sin(beta / 2) / 2, thr = atol)
      if(allocated(error)) return

    enddo

    nullify(beta)

  end subroutine test_little_d_special_values_diag

! ================================================================================================================================ !
end module test_wignerd
! ================================================================================================================================ !

! ================================================================================================================================ !
program tester

  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_wignerd, only: collect_wignerd

  implicit none

  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("WignerD", collect_wignerd) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat .gt. 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

! ================================================================================================================================ !
end program tester
! ================================================================================================================================ !
