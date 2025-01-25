! ================================================================================================================================ !
module wignerd
  !! Module for computing the Wigner D-matrix and d-matrix (big D and little d; Fortran is case insensitive).
  !! The d-matrix can be calculated via the analytic expression, but this can overflow for large j due to
  !! the ratios of large factoials. The d-matrix can instead be calculated via sparse banded matrix
  !! diagonalization. See DOI: 10.1103PhysRevE.92.043307 for more details.

  implicit none

  private

  public :: wigner_d
  public :: wigner_big_D
  public :: wigner_little_d

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  interface wigner_d
    !! Interface to access the routines to calculate the matrix D(α,β,γ) or d(β)
    !! based on the number of arguments
    module procedure :: wigner_big_D_jint
    module procedure :: wigner_big_D_jreal
    module procedure :: wigner_little_d_jint
    module procedure :: wigner_little_d_jreal
  end interface wigner_d

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  interface wigner_big_D
    !! Interface to access the routines to calculate D(α,β,γ) via matrix diagonalization
    !! or the analytic expression
    module procedure :: wigner_big_D_jint
    module procedure :: wigner_big_D_jreal
  end interface wigner_big_D

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  interface wigner_little_d
    !! Interface to access the routines to calculate d(β) via matrix diagonalization
    !! or the analytic expression
    module procedure :: wigner_little_d_jint
    module procedure :: wigner_little_d_jreal
  end interface wigner_little_d

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_big_D_jint(j, euler_alpha, euler_beta, euler_gamma, use_analytic) result(D)
    !! Return the wigner D-matrix \(D^j_{m',m}(\alpha, \beta, \gamma)\) for integer j

    use wignerd__types,     only: rp
    use wignerd__constants, only: im

    implicit none

    integer,  intent(in) :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_alpha
      !! The Euler angle α
    real(rp), intent(in) :: euler_beta
      !! The Euler angle β
    real(rp), intent(in) :: euler_gamma
      !! The Euler angle γ
    logical,  intent(in), optional :: use_analytic
      !! Force the use of the analytic expression to obtain \(d^j(\beta)\) ?
    complex(rp), allocatable :: D(:,:)
      !! The Wigner D-matrix \(D^j(\alpha,\beta,\gamma)\)

    logical :: use_analytic_local

    use_analytic_local = .false. ; if(present(use_analytic)) use_analytic_local = use_analytic

    D = wigner_big_D_jreal(real(j, kind = rp), euler_alpha, euler_beta, euler_gamma, use_analytic_local)

  end function wigner_big_D_jint
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_big_D_jreal(j, euler_alpha, euler_beta, euler_gamma, use_analytic) result(D)
    !! Return the wigner D-matrix \(D^j_{m',m}(\alpha, \beta, \gamma)\) for real j

    use wignerd__types,     only: rp
    use wignerd__constants, only: one, im

    implicit none

    real(rp), intent(in)  :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_alpha
      !! The Euler angle α
    real(rp), intent(in) :: euler_beta
      !! The Euler angle β
    real(rp), intent(in) :: euler_gamma
      !! The Euler angle γ
    logical,  intent(in), optional :: use_analytic
      !! Force the use of the analytic expression to obtain \(d^j(\beta)\) ?
    complex(rp), allocatable :: D(:,:)
      !! The Wigner D-matrix \(D^j(\alpha,\beta,\gamma)\)

    logical :: use_analytic_local

    integer :: n
    integer :: i, k
    real(rp) :: m, mp
    real(rp), allocatable :: little_d(:,:)

    n = nint(2*j) + 1

    use_analytic_local = .false. ; if(present(use_analytic)) use_analytic_local = use_analytic

    ! -- allocate D
    D = wigner_little_d_jreal(j, euler_beta, use_analytic)

    little_d = D

    do concurrent (i = 1 : n, k = 1 : n)
      mp = -j + i - 1
      m  = -j + k - 1
      D(i, k) = exp(-im * mp * euler_alpha) * little_d(i, k) * exp(-im * m * euler_gamma)
    enddo

  end function wigner_big_D_jreal

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_little_d_jint(j, euler_beta, use_analytic) result(d)
    !! Return the Wigner d-matrix \(d^j_{m'm}(\beta)\) for integer j

    use wignerd__types, only: rp

    implicit none

    integer, intent(in) :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_beta
      !! The euler angle β
    logical, intent(in), optional :: use_analytic
      !! Force the use of the analytic expression to obtain \(d^j(\beta)\) ?
    real(rp), allocatable :: d(:,:)
      !! The matrix \(d^j(\beta)\)

    logical :: use_analytic_Local

    use_analytic_local = .false. ; if(present(use_analytic)) use_analytic_local = use_analytic

    d = wigner_little_d_jreal(real(j, kind = rp), euler_beta, use_analytic_local)

  end function wigner_little_d_jint
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_little_d_jreal(j, euler_beta, use_analytic) result(d)
    !! Return the Wigner d-matrix \(d^j_{m'm}(\beta)\) for real j

    use wignerd__types,     only: rp
    use wignerd__constants, only: one

    implicit none

    real(rp), intent(in) :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_beta
      !! The euler angle β
    logical, intent(in), optional :: use_analytic
      !! Force the use of the analytic expression to obtain \(d^j(\beta)\) ?
    real(rp), allocatable :: d(:,:)
      !! The matrix \(d^j(\beta)\)

    logical :: use_analytic_local

    use_analytic_local = .false. ; if(present(use_analytic)) use_analytic_local = use_analytic

    if(use_analytic_local .eqv. .true.) then
      ! -- use the analytic expression
      d = wigner_little_d_analytic(j, euler_beta)
      return
    endif

    ! -- use the matrix diagonalization
    d = wigner_little_d_diag(j, euler_beta)

  end function wigner_little_d_jreal

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_little_d_analytic(j, euler_beta) result(little_d)
    !! Return the Wigner d-matrix \(d^j_{m',m}(\beta)\) via the analytic expression for real \(j\)
    use wignerd__types,     only: rp
    use wignerd__system,    only: die
    use wignerd__constants, only: zero, two
    use wignerd__functions, only: fact => factorial

    implicit none

    real(rp), intent(in) :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_beta
      !! The euler angle β
    real(rp), allocatable :: little_d(:,:)
      !! The matrix \(d^j(\beta)\)

    integer :: i, k
    integer :: n
    integer :: two_s

    real(rp) :: d
    real(rp) :: m, mp
    real(rp) :: s
    real(rp) :: smin, smax
    real(rp) :: numer, denom

    if(nint(2*j) .lt. 0) call die("The angular momentum j cannot be negative")

    n = nint(2*j) + 1

    allocate(little_d(n, n))

    do concurrent(i = 1 : n, k = 1 : n)

      mp = -j + i - 1
      m  = -j + k - 1

      smin = max(zero, m - mp)
      smax = min(j + m, j - mp)

      d = 0
      do  two_s = nint(2*smin), nint(2*smax), 2
        s = two_s / two
        numer =  (-1) ** (mp - m + s) &
          * cos(euler_beta / 2) ** (2*j + m - mp - 2*s) &
          * sin(euler_beta / 2) ** (mp - m + 2*s)
        denom = fact(j + m - s) * fact(s) * fact(mp - m + s) * fact(j - mp - s)
        d = d + numer / denom
      enddo

      little_d(i, k) = d * sqrt( fact(j + mp) * fact(j - mp) * fact(j + m) * fact(j - m) )

    enddo

  end function wigner_little_d_analytic

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function wigner_little_d_diag(j, euler_beta) result(little_d)
    !! Return the Wigner d-matrix \(d^j_{m',m}(\beta)\) via matrix diagonalization for real j
    !! see DOI: 10.1103PhysRevE.92.043307 for more details

    use wignerd__types,      only: rp
    use wignerd__system,     only: die
    use wignerd__constants,  only: zero
    use wignerd__characters, only: int2char

    implicit none

    real(rp), intent(in) :: j
      !! The angular momentum \(j\)
    real(rp), intent(in) :: euler_beta
      !! The Euler angle β
    complex(rp), allocatable :: Jy(:,:)
      !! The operator \(J_y\) in matrix form
    real(rp), allocatable :: little_d(:,:)
      !! The matrix \(d^j(\beta)\)

    integer :: n

    integer  :: row, col
    real(rp) :: mrow, mcol

    real(rp), allocatable :: eigvals(:)
    complex(rp), allocatable :: eigvecs(:,:)

    if(nint(2*j) .lt. 0) call die("The angular momentum j cannot be negative")

    n = nint(2*j) + 1

    allocate(Jy(n, n)) ; Jy = zero

    ! -- iterate over the rows and columns to populate Jy(:,:).
    !    The projections m of j for the rows and columns go from -j to j in steps of 1
    do concurrent( row = 1 : n, col = 1 : n )
      mrow = -j + (row - 1)
      mcol = -j + (col - 1)
      Jy(row, col) = Jy_element(j, mrow, mcol)
    enddo

    ! -- diagonalize the Jy matrix
    diag: block

      use wignerd__characters,             only: int2char
#ifdef USE_EXTERNAL_LAPACK
      use wignerd__lapack_interface, only: zhbev
#else
      use stdlib_linalg_lapack, only: zhbev => stdlib_zhbev
#endif

      integer :: i, k
      integer :: kd
      integer :: ldab, ldz
      integer :: info
      integer :: nrwork
      real(rp), allocatable :: rwork(:)
      complex(rp), allocatable :: work(:)
      complex(rp), allocatable :: AB(:,:)
      character :: jobz, uplo

      jobz   = "V" ! -- get eigenvectors too
      uplo   = "U" ! -- store upper triangle of Jy
      kd     = 1   ! -- one superdiagonal
      ldab   = kd + 1
      ldz    = n
      nrwork = max(1, 3*n - 2)

      allocate(eigvals(n))
      allocate(work(n))
      allocate(rwork(nrwork))
      allocate(eigvecs(ldz, n))
      allocate(AB(ldab, n))

      ! -- build the band matrix AB
      AB = zero
      do k = 1, n
        do i = max(1, k - kd), k
          AB(kd + 1 + i - k, k) = Jy(i, k)
        enddo
      enddo

      call zhbev(jobz, uplo, n, kd, AB, ldab, eigvals, eigvecs, ldz, work, rwork, info)

      if(info .ne. 0) call die("Procedure ZHBEV returned with INFO = " // int2char(info))

    end block diag

    build_little_d: block

      use wignerd__constants, only: im

      integer :: i, k
      integer :: imu
      real(rp) :: cos_term(n), sin_term(n)

      allocate(little_d(n, n))

      ! -- fill the wigner d-matrix with its elements
      do concurrent (i = 1 : n, k = 1 : n)

        ! -- the inner product $$<i|\mu>y$$ just picks out the element i of the vector $$|\mu>_y$$ because
        !    the basis Jz is just the unit vectors. The product $$e^{-i\mu\beta} \left<i|\mu\right>_y{}_y\left<\mu|k\right>$$ is purely real, so
        !    so we only compute the real part
        associate(mu => eigvals)
          cos_term = cos(mu(:) * euler_beta) * (eigvecs(i, :) % re * eigvecs(k, :) % re + eigvecs(i, :) % im * eigvecs(k, :) % im)
          sin_term = sin(mu(:) * euler_beta) * (eigvecs(i, :) % im * eigvecs(k, :) % re - eigvecs(i, :) % re * eigvecs(k, :) % im)
          little_d(i, k) = sum(cos_term + sin_term)
        end associate

      enddo

    end block build_little_d

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !

    ! -- !
    pure elemental function Jy_element(j, m, n) result(output)
      !! Returns the matrix element \(<jm'|J_y|jm>\) from the above reference:
      !!   $$<jm | J_y | jn> = [ X_{-n} δ_{m,n+1} - X_{n}δ_{m,n-1} ] / (2i)$$
      !! where \(X_m = sqrt( (j+m)(j-m+1) )\)

      use wignerd__types,     only: rp
      use wignerd__functions, only: delta
      use wignerd__constants, only: im, zero

      implicit none

      real(rp), intent(in) :: j
        !! The angular momentum \(j\)
      real(rp), intent(in) :: m, n
        !! The projections of \(j\)

      integer :: double_m, double_n
        !! for passing integer arguments to delta() in case of half integer projections
      complex(rp) :: output

      double_m = nint(2*m)
      double_n = nint(2*n)
      output = ( X(j, -n) * delta(double_m, double_n + 2) - X(j,n) * delta(double_m, double_n - 2) ) / (2 * im)

    end function Jy_element

    ! -- !
    pure elemental function X(j, m) result(output)
      !! Returns the quantity \(X_m = sqrt( (j + m) (j - m + 1) )\)

      implicit none

      real(rp), intent(in) :: j, m
      real(rp) :: arg
      real(rp) :: output

      output = sqrt( (j + m) * (j - m + 1) )

    end function X

  end function wigner_little_d_diag

! ================================================================================================================================ !
end module wignerd
! ================================================================================================================================ !
