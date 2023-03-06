# TVproximal
This code implements Newton-like algorithms for computing the [proximal operator](https://en.wikipedia.org/wiki/Proximal_operator) of the [anisotropic total-variation](https://en.wikipedia.org/wiki/Total_variation_denoising) map, i.e., solving the optimization problem:

$$ \mathrm{prox}(y) = \underset{x\in\mathbb{R}^n}{\text{argmin}}~\frac{1}{2}\Vert x - y\Vert_2^2 + \Vert x\Vert_{TV,p}$$

where $\Vert\cdot\Vert_{TV,p}$ denotes the $
\ell-p$ [total-variation map](https://en.wikipedia.org/wiki/Total_variation) in one or two variables for $p\in\{0,1\}$.  More explicitly, for one dimensional variables $x\in\mathbb{R}^n$:

$$ \Vert x\Vert_{TV,p} = \Vert Dx\Vert_p = \left(\sum_{i=1}^{n-1} \vert x_{i+1}-x_i\vert^p\right)^{\frac{1}{p}}.$$

or if $x\in\mathbb{R}^{m\times n}$, the anisotropic total variation becomes:

$$ \Vert x\Vert_{TV,p} = \sum_{i=1}^{n-1}\sum_{j=1}^{m-1}\Vert x_{:,j}\Vert_{TV,p} + \Vert x_{i,:}\Vert_{TV,p}.$$
The code solves these optimization problesm for a given input $x$, effectively solving the "total variation denoising" problem.  The algorithm uses a second-order Newton method to rapidly solve the optimization problem, as Newton algorithms exhibit (roughly) linear convergence rates.

# Dependencies
The code is written directly in FORTRAN95, built using the [LAPACK](https://netlib.org/lapack/) libraries' interfaces directly (as they are also written directly in FORTRAN). It relies on a FORTRAN compiler, but was built with `gfortran` in mind (this is the compiler the included `Makefile` looks for).
