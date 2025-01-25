# WignerD

Provides routines to return the Wigner D-matrix

$$
\begin{equation}
D^j_{m'm}(\alpha,\beta,\gamma) = e^{-im'\alpha} d^j_{m'm}(\beta) e^{-im\gamma}
\end{equation}
$$

and the d-matrix

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

for a given value of the angular momentum j.
By default, these are calculated via matrix diagonalization using the method of Feng <i>et al.</i>[[1]](1), but the use of the analytic expression for the d-matrix can be forced as well.

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
A few tests are included for explicit values of the d-matrix for several values of the angle β.
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
| `wigner_little_d(...)`  | Returns $d^j_{m'm}(\beta)$ via the analytic expression or matrix diagonalization. See [docs](https://banana-bred.github.io/WignerD/) for details.|
| `wigner_big_D(...)`     | Returns $D^j_{m'm}(\alpha,\beta,\gamma)$ via the analytic expression or matrix diagonalization. See [docs](https://banana-bred.github.io/WignerD/) for details.|

More info on input/output types throughout the [docs](https://banana-bred.github.io/WignerD/).


[[1]](1)

---

## Reference(s)

<a id="1">[1]</a>
X. M. Feng, P. Wang, W. Yang, and G. R. Jin,
*High-precision evaluation of Wigner's \(d\) matrix by exact diagonalization*,
Phys. Rev. E. 2015, 92, 043307,
URL: [https://doi.org/10.1103/PhysRevE.92.043307](https://doi.org/10.1103/PhysRevE.92.043307)
