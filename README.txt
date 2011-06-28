Q Math Library

1. Introduction

The Q Math Library provides the q programming language and KDB+ database with
an interface to a number of useful mathematical functions from the FDLIBM,
Cephes and LAPACK libraries.

This is the first release of the Q Math Library, version 0.1.1.  Only the two
Windows platforms are currently supported.


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

    q).qml.ncdf .25 .5 .75            / quartiles of a normal distribution
    -0.6744898 0 0.6744898

    q).qml.mchol (1 2 1;2 5 4;1 4 6)  / Cholesky decomposition
    1 2 1
    0 1 2
    0 0 1

    q).qml.poly 2 -9 16 -15           / solve 2x^3-9x^2+16x-15=0
    2.5
    1 1.414214
    1 -1.414214

It's probably a good idea to run test suite program, test.q, at some point.

Rebuilding requires Cygwin.  The top of the Makefile has some configuration
options.


4. Constants and functions

  pi              pi
  e               e
  eps             smallest representable step from 1.

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
  c2cdf[k;x]              CDF of chi-squared distribution (k>=1,x>=0)
  c2icdf[k;p]     inverse CDF of chi-squared distribution
                            (k>=1,p>=.5)
  stcdf[k;x]              CDF of Student's t-distribution (natural k)
  sticdf[k;p]     inverse CDF of Student's t-distribution (natural k)
  fcdf[d1;d2;x]           CDF of F-distribution (d1,d2>=1,x>=0)
  ficdf[d1;d2;p]  inverse CDF of F-distribution (d1,d2>=1,x>=0)
  gcdf[k;th;x]            CDF of gamma distribution (x>=0)
  gicdf[k;th;x]   inverse CDF of gamma distribution (x>=0,p>=.5)

  mdiag[diag]     make diagonal matrix
  mdet[matrix]    determinant
  minv[matrix]    inverse
  mpinv[matrix]   pseudoinverse
  mev[matrix]     (eigenvalues; eigenvectors) sorted by decreasing modulus
  mchol[matrix]   Cholesky decomposition upper matrix
  msvd[matrix]    singular value decomposition: (U; Sigma; V)

  poly[coef]      roots of a polynomial (highest-degree coefficient first, can
                    be complex)


5. Updates and feedback

This library is hosted at http://althenia.net/qml.  It is programmed by
Andrey Zholos <aaz@althenia.net>.  Comments, bug reports and testing results
are welcome and will be appreciated.
