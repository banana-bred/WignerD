---
project: WignerD
summary: Provides routines to return the Wigner matrices D(αβγ) and d(β).
         Public interfaces for accessing these routines are made available in the `wignerd` module.
extensions: f
fixed_extensions:
max_frontpage_items: 6
display: public
---
Provides procedures to calculate the Wigner matrices
$$
\begin{equation}
D^j_{m'm}(\alpha,\beta,\gamma) = e^{-im'\alpha} d^j_{m'm}(\beta) e^{-im\gamma}
\end{equation}
$$
and
$$
\begin{equation}
\begin{aligned}
    d^j_{m'm}(\beta) &=
    \sqrt{ (j+m')!(j-m')!(j+m)!(j-m)! }
    \\
    &\quad \times
    \sum\limits_{s = s_\text{min}}^{s_\text{max}}
    \left[
        \frac{
            (-1)^{m'-m+s} \left(\cos \frac{\beta}{2}\right)^{2j+m-m'-2s} \left(\sin \frac{\beta}{2}\right)^{m'-m+2s}
        }{
            (j+m-s)!s!(m'-m+s)!(j-m'-s)!
        }
    \right]
\end{aligned}
\end{equation}
$$
in their real forms.
The \(d\) matrix can be obtained via its analytic expression (2) or via matrix diagonalization<sup>[[1]](#1)</sup>.
**The matrix diagonalization method is used by default**, but the analytic method can be forced with an optional parameter in all interfaces.
The factorial terms in (2) can easily overflow double precision for approximately the following values of j:

<center>

| # bits | critical \(j\) |
| ------ | ------------ |
| 32  | 17  |
| 64  | 86  |
| 128 | 878 |

</center>

This problem is avoided when using the matrix diagonalization method of Feng <i>et al.</i><sup>[[1]](#1)</sup>.
For larger \(j\), the matrix diagonalization method quickly becomes faster than the analytic version.

### Dependencies
- LAPACK's `ZHBEV` routine, which is made available in the [Fortran standard library](https://stdlib.fortran-lang.org/).
- [*optional*] The [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm) for easy building.

### Building with [fpm](https://github.com/fortran-lang/fpm)
In the package directory, run

    $ fpm build --profile release

The archive file `libWignerD.a` and several `.mod` files will be placed in the generated `build` subdirectory.
If you'd rather use your local version of LAPACK, use the flag

    $ fpm build --profile release --flag="-DUSE_EXTERNAL_LAPACK"

### Building without [fpm](https://github.com/fortran-lang/fpm)
Assuming you have a local installation of LAPACK and that your linker program knows where to find it, just run the provided compile script:

    $ ./compile

The default compiler is `gfortran`.
The archive file `libwignerd.a` and several `.mod` files will be placed in the generated `build/lib` and `build/mod` subdirectories.
These will be needed for reference by another program.

### Testing
A few tests are included for explicit values of the \(d^j_{m'm}(\beta)\) for several values of the angle \(\beta\).
Just run

    $ fpm test

## Usage

To use this project within your [fpm](https://github.com/fortran-lang/fpm) project, add the following to your `fpm.toml` file:

    [dependencies]
    wignerd = { git = "https://github.com/banana-bred/WignerD" }

Otherwise, you will just need the generated archive and mod files mentioned above.
Don't forget to tell your compiler where they are.

### Available routines/interfaces/procedures
The module `wignerd` contains the following public interfaces, which can be accessed via the `use` statement :

| interface               | description |
| ----------------------- | ----------- |
| `wigner_d(...)`         | This is an interface wrapper for `wigner_little_d` and `wigner_big_d`, depending on the number of arguments that you provide.|
| `wigner_little_d(...)`  | Returns \(d^j_{m'm}(\beta)\) via the analytic expression or matrix diagonalization. See docs for details.|
| `wigner_big_D(...)`     | Returns \(D^j_{m'm}(\alpha,\beta,\gamma)\) via the analytic expression or matrix diagonalization. See docs for details.|

More info on input/output types throughout the docs.

### Example
Calculate \(D^{1/2}(\alpha,\beta,\gamma)\) and \(d^{1/2}(\beta)\):

    program D

        use, intrinsic :: iso_fortran_env, only: wp => real64
        use wignerd,                       only: wigner_d, wigner_big_D, wigner_little_d
        use wignerd__constants,            only: one, two, pi

        real(wp), allocatable :: little_d1(:,:)
        real(wp), allocatable :: little_d2(:,:)
        real(wp), allocatable :: big_D1(:,:)
        real(wp), allocatable :: big_D2(:,:)
        real(wp) :: j, euler_alpha, euler_beta, euler_gamma

        j = one / two

        euler_alpha = pi / 6
        euler_beta  = pi / 2
        euler_gamma = pi

        little_d1 = wigner_d(j, euler_beta)                            ! -- calculate d^j(β)
        little_d2 = wigner_little_d(j, euler_beta)                     ! -- calculate d^j(β)
        big_D1 = wigner_d(j, euler_alpha, euler_beta, euler_gamma)     ! -- calculate D^j(α,β,γ)
        big_D2 = wigner_big_D(j, euler_alpha, euler_beta, euler_gamma) ! -- calculate D^j(α,β,γ)

    end program D

In the above, arrays `big_D1` and `big_D2` will hold the same information because they end up calling the same routines (the same goes for `little_d1` and `little_d2`).
The above calls default to the diagonalization method to obtain \(D^{1/2}(\alpha,\beta,\gamma)\) and \(d^{1/2}(\beta)\).
We can force the use of the analytic expression via the optional input parameter `use_analytic`:

    little_d1 = wigner_little_d(j, euler_beta) ! <--------------------------------- matrix diagonalization
    little_d2 = wigner_little_d(j, euler_beta, use_analytic = .true.) ! <---------- analytic
    big_D1    = wigner_big_D(j, euler_alpha, euler_beta, euler_beta) ! <----------- matrix diagonalization
    big_D2    = wigner_big_D(j,  euler_alpha, euler_beta, euler_beta, .true.) ! <-- analytic

#### Happy rotating !

---

## Reference(s)

<a id="1">[1]</a>
X. M. Feng, P. Wang, W. Yang, and G. R. Jin,
*High-precision evaluation of Wigner's \(d\) matrix by exact diagonalization*,
Phys. Rev. E. 2015, 92, 043307,
URL: [https://doi.org/10.1103/PhysRevE.92.043307](https://doi.org/10.1103/PhysRevE.92.043307)
