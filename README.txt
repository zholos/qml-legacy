Q Math Library

1. Introduction

The Q Math Library provides the q programming language and KDB+ database with
an interface to a number of useful mathematical functions from the FDLIBM,
Cephes, LAPACK and CONMAX libraries.

This is version 0.2.1 of the Q Math Library.


2. Licensing

The Q Math Library is free software, available under a BSD-style license.
It is provided in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranties of MERCHANTABILITY and FITNESS FOR A
PARTICULAR PURPOSE.  See the file LICENSE.txt for more details.

This library is intended to be linked, and the precompiled binaries are
linked, against several other libraries.  The copyrights and licenses for
these libraries are also listed in the file LICENSE.txt.


3. Usage

Two files make up the library: qml.dll (in a platform-specific subdirectory)
and qml.q.  Put them both somewhere where q can find them and load with

    q)\l qml.q

All functions will be in the .qml namespace.  The functions accept any
numerical arguments and convert them into floating-point.  Matrixes are in
row-major order, as usual.  Complex numbers are represented as pairs of the
real and imaginary parts.  E.g.:

    q).qml.ncdf .25 .5 .75                    / normal distribution quartiles
    -0.6744898 0 0.6744898

    q).qml.mchol (1 2 1;2 5 4;1 4 6)          / Cholesky factorization
    1 2 1
    0 1 2
    0 0 1

    q).qml.poly 2 -9 16 -15                   / solve 2x^3-9x^2+16x-15=0
    2.5
    1 1.414214
    1 -1.414214

    q).qml.conmin[{x*y+1};{1-(x*x)+y*y};0 0]  / minimize x(y+1) s.t. x^2+y^2<=1
    -0.8660254 0.5

It's recommended to run the test suite, test.q, to make sure that everything
is working correctly.

To rebuild the library, use the Makefile, which is for GNU make.  There are
configuration options near the top.  On Windows, rebuilding requires Cygwin or
another GNU environment.  On other platforms, GCC is required.  Although only
tested on the platforms for which binaries are included, it should compile on
any platform with a few tweaks.  If you rebuild the library, make sure to run
the test suite before using it.


4. Constants and functions

  pi              pi
  e               e
  eps             smallest representable step from 1.

  sin[x]          sine
  cos[x]          cosine
  tan[x]          tangent
  asin[x]         arcsine
  acos[x]         arccosine
  atan[x]         arctangent
  atan2[x;y]      atan[x%y]
  sinh[x]         hyperbolic sine
  cosh[x]         hyperbolic cosine
  tanh[x]         hyperbolic tangent
  asinh[x]        hyperbolic arcsine
  acosh[x]        hyperbolic arccosine
  atanh[x]        hyperbolic arctangent

  exp[x]          exponential
  expm1[x]        exp[x]-1
  log[x]          logarithm
  log10[x]        base-10 logarithm
  logb[x]         extract binary exponent
  log1p[x]        log[1+x]
  pow[a;x]        exponentiation
  sqrt[x]         square root
  cbrt[x]         cube root
  hypot[x;y]      sqrt[pow[x;2]+pow[y;2]]
  floor[x]        round downward
  ceil[x]         round upward
  fabs[x]         absolute value
  fmod[x;y]       remainder of x%y

  erf[x]          error function
  erfc[x]         complementary error function
  lgamma[x]       log of absolute value of gamma function
  gamma[x]        gamma function
  beta[x;y]       beta function
  pgamma[a;x]     lower incomplete gamma function (a>0)
  pgammac[a;x]    upper incomplete gamma function (a>0)
  pgammar[a;x]    regularized lower incomplete gamma function (a>0)
  pgammarc[a;x]   regularized upper incomplete gamma function (a>0)
  ipgammarc[a;p]  inverse complementary regularized incomplete gamma function
                    (a>0,p>=.5)
  pbeta[a;b;x]    incomplete beta function (a,b>0)
  pbetar[a;b;x]   regularized incomplete beta function (a,b>0)
  ipbetar[a;b;p]  inverse regularized incomplete beta function (a,b>0)
  j0[x]           order 0 Bessel function
  j1[x]           order 1 Bessel function
  y0[x]           order 0 Bessel function of the second kind
  y1[x]           order 1 Bessel function of the second kind

  ncdf[x]                 CDF of normal distribution
  nicdf[p]        inverse CDF of normal distribution
  c2cdf[k;x]              CDF of chi-squared distribution (k>=1)
  c2icdf[k;p]     inverse CDF of chi-squared distribution (k>=1)
  stcdf[k;x]              CDF of Student's t-distribution (natural k)
  sticdf[k;p]     inverse CDF of Student's t-distribution (natural k)
  fcdf[d1;d2;x]           CDF of F-distribution (d1,d2>=1,x>=0)
  ficdf[d1;d2;p]  inverse CDF of F-distribution (d1,d2>=1,x>=0)
  gcdf[k;th;x]            CDF of gamma distribution
  gicdf[k;th;p]   inverse CDF of gamma distribution
  bncdf[k;n;p]            CDF of binomial distribution
  bnicdf[k;n;x]   inverse CDF of binomial distribution for p parameter (k<n)
  pscdf[k;lambda]         CDF of Poisson distribution
  psicdf[k;p]     inverse CDF of Poisson distribution for lambda parameter
  smcdf[n;e]              CDF for one-sided Kolmogorov-Smirnov test
  smicdf[n;e]     inverse CDF for one-sided Kolmogorov-Smirnov test
  kcdf[x]                 CDF for Kolmogorov distribution
  kicdf[p]        inverse CDF for Kolmogorov distribution (p>=1e-8)

  diag[diag]      make diagonal matrix
  mdiag[matrix]   extract main diagonal
  mdet[matrix]    determinant
  mrank[matrix]   rank
  minv[matrix]    inverse
  mpinv[matrix]   pseudoinverse
  mev[matrix]     (eigenvalues; eigenvectors) sorted by decreasing modulus
  mchol[matrix]   Cholesky factorization upper matrix
  mqr[matrix]     QR factorization: (Q; R)
  mqrp[matrix]    QR factorization with column pivoting:
                    (Q; R; P), matrix@\:P=Q mmu R
  mlup[matrix]    LUP factorization with row pivoting:
                    (L; U; P), matrix[P]=L mmu U
  msvd[matrix]    singular value decomposition: (U; Sigma; V)

  poly[coef]      roots of a polynomial (highest-degree coefficient first, can
                    be complex)

  root[f;(x0;x1)]         find root on interval (f(x0)f(x1)<1)
  rootx[opt;f;(x0;x1)]    root[] with options (as dictionary or mixed list)
                           `iter:  max iterations         (default: 100)
                           `tol:   numerical tolerance    (default: ~1e-8)
                           `full:  full output            (default: only x)
                           `quiet: return null on failure (default: signal)
  solve[eqs;x0]           solve nonlinear equations (given as functions)
  solvex[opt;eqs;x0]      solve[] with options
                           `iter:  max iterations         (default: 1000)
                           `tol:   numerical tolerance    (default: ~1e-8)
                           `full:  full output            (default: only x)
                           `quiet: return null on failure (default: signal)
                           `steps: RK steps per iteration (default: 1)
                           `rk:    use RK steps only      (default: RK, SLP)
                           `slp:   use SLP steps only     (default: RK, SLP)
  line[f;base;x0]         line search for minimum from base
  linex[opt;f;base;x0]    line[] with same options as rootx[]
  min[f;x0]               find unconstrained minimum
  minx[opt;f;x0]          min[] with same options as solvex[]
  conmin[f;cons;x0]       find constrained minimum (functions cons>=0)
  conminx[opt;f;cons;x0]  min[] with same options as solvex[], plus
                           `lincon: assume linear cons    (default: nonlinear)


5. Updates and feedback

This library is hosted at http://althenia.net/qml.  It is programmed by
Andrey Zholos <aaz@althenia.net>.  Comments, bug reports and testing results
are welcome.
