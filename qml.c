#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <k.h>

#ifdef _WIN32
    #define DLL_EXPORT __declspec(dllexport)
#else
    #define DLL_EXPORT
#endif


//
// Utility functions
//

#if __GNUC__ >= 3
    #define unlikely(x) __builtin_expect(!!(x), 0)
#else
    #define unlikely(x) (x)
#endif

#define check(x, r, c)        \
    do                        \
        if (unlikely(!(x))) { \
            c;                \
            return (r);       \
        }                     \
    while (0)

#define check_type(x, c)    check(x, krr(ss("type")), c)
#define check_length(x, c)  check(x, krr(ss("length")), c)


#define check_lapack_return(info, c) \
    check((info) >= 0, krr(ss(info == -7 ? "wsfull" : "qml_assert")), c)


#define alloc(x, n, c)                                                 \
    do {                                                               \
        I alloc__n = (n);                                              \
        check(alloc__n >= 0 && alloc__n < wi &&                        \
              alloc__n <= SIZE_MAX / sizeof *x, krr(ss("limit")),  c); \
        check(x = malloc(alloc__n * sizeof *x), krr(ss("wsfull")), c); \
    } while (0)

static I
take_maxwork(F maxwork) {
    if (!(maxwork >= 0 && maxwork < wi))
        return wi;
    return (I)maxwork;
}

static I
add_size(I a, I m, I n) {
    check(!n || m <= (wi - a) / n, wi,);
    return a + m * n;
}


// krr() might return NULL or an error object.
// This dummy object represents no error.
static struct k0 no_error_;
#define no_error (&no_error_)

// Used internally as both a default empty list and as a flag.
static struct k0 empty_con_ = { 0 };
#define empty_con (&empty_con_)

#define bubble_error(error, c)                     \
    do {                                           \
        K bubble__error = (error);                 \
        if (unlikely(bubble__error != no_error)) { \
            c;                                     \
            return bubble__error;                  \
        }                                          \
    } while (0)


#define repeat(i, n) \
    for (I i = 0, n__##i = (n); i < n__##i; i++)


static I
min_i(I x, I y) {
    return x <= y ? x : y;
}

static I
max_i(I x, I y) {
    return x >= y ? x : y;
}

static void
swap_i(I* x, I* y) {
    I v = *x; *x = *y; *y = v;
}

static void
swap_f(F* x, F* y) {
    F v = *x; *x = *y; *y = v;
}


//
// Input value conversion
//

static int
is_i(K x) {
    return xt==-KI || xt==-KH || xt==-KB;
}

static int
is_f(K x) {
    return xt>=-KF && xt<=-KH || xt==-KB;
}

static int // if true, x->n is length
is_I(K x) {
    if (xt>0)
        return xt==KB || xt==KH || xt==KI;
    if (!xt) {
        repeat (i, xn)
            check(is_i(xK[i]), 0,);
        return 1;
    }
    return 0;
}

static int // if true, x->n is length
is_F(K x) {
    if (xt>0)
        return xt==KB || xt>=KH && xt<=KF;
    if (!xt) {
        repeat (i, xn)
            check(is_f(xK[i]), 0,);
        return 1;
    }
    return 0;
}

static I
as_i(K x) {
    switch (xt) {
    case -KI:
        return xi;
    case -KH:
        return xh==nh ? ni : xh==wh ? wi : xh==-wh ? -wi : xh;
    case -KB:
        return xg;
    }
    return ni; /* not reached */
}

static F
as_f(K x) {
    switch (xt) {
    case -KF:
        return xf;
    case -KE:
        return xe;
    case -KJ:
        return xj==nj ? nf : xj==wj ? wf : xj==-wj ? -wf : xj;
    case -KI:
        return xi==ni ? nf : xi==wi ? wf : xi==-wi ? -wf : xi;
    case -KH:
        return xh==nh ? nf : xh==wh ? wf : xh==-wh ? -wf : xh;
    case -KB:
        return xg;
    }
    return ni; /* not reached */
}

static I
item_I(K x, I i) {
    switch (xt) {
    case KI:
        return xI[i];
    case KH:
        return xH[i]==nh ? ni : xH[i]==wh ? wi : xH[i]==-wh ? -wi : xH[i];
    case KB:
        return xG[i];
    case 0:
        return as_i(xK[i]);
    }
    return ni; /* not reached */
}

static F
item_F(K x, I i) {
    switch (xt) {
    case KF:
        return xF[i];
    case KE:
        return xE[i];
    case KJ:
        return xJ[i]==nj ? nf : xJ[i]==wj ? wf : xJ[i]==-wj ? -wf : xJ[i];
    case KI:
        return xI[i]==ni ? nf : xI[i]==wi ? wf : xI[i]==-wi ? -wf : xI[i];
    case KH:
        return xH[i]==nh ? nf : xH[i]==wh ? wf : xH[i]==-wh ? -wf : xH[i];
    case KB:
        return xG[i];
    case 0:
        return as_f(xK[i]);
    }
    return ni; /* not reached */
}

static void
copy_F(K x, F* f, I step) {
    if (xt == KF)
        if (step == 1)
            memcpy(f, xF, xn * sizeof *f);
        else
            repeat (i, xn)
                f[i * step] = xF[i];
    else
        repeat (i, xn)
            f[i * step] = item_F(x, i);
}


//
// Scalar functions
//

#define wrap_last(x, call)       \
    do {                         \
        if (is_f(x)) {           \
            F v = as_f(x);       \
            return kf(call);     \
        }                        \
        check_type(is_F(x),);    \
        K we__r = ktn(KF, x->n); \
        repeat (i, x->n)  {      \
            F v = item_F(x, i);  \
            kF(we__r)[i] = call; \
        }                        \
        return we__r;            \
    } while (0)


// Wrap a function of (atom/vector F)
#define wrap_F(name, func) \
K DLL_EXPORT               \
qml_##name(K x) {          \
    wrap_last(x, func(v)); \
}

#define wrap_F_(func) \
    wrap_F(func, func)


// Wrap a function of (atom I, atom/vector F)
#define wrap_iF(name, func)     \
K DLL_EXPORT                    \
qml_##name(K x, K y) {          \
    check_type(is_i(x),);       \
    I x_i = as_i(x);            \
    wrap_last(y, func(x_i, v)); \
}


// Wrap a function of (atom F, atom/vector F)
#define wrap_fF(name, func)     \
K DLL_EXPORT                    \
qml_##name(K x, K y) {          \
    check_type(is_f(x),);       \
    F x_f = as_f(x);            \
    wrap_last(y, func(x_f, v)); \
}


// Wrap a function of (atom/vector F, atom/vector F)
#define wrap_FF(name, func)                          \
K DLL_EXPORT                                         \
qml_##name(K x, K y) {                               \
    if (is_f(y)) {                                   \
        F y_f = as_f(y);                             \
        wrap_last(x, func(v, y_f));                  \
    }                                                \
    check_type(is_F(y),);                            \
    if (is_f(x)) {                                   \
        F x_f = as_f(x);                             \
        K r = ktn(KF, y->n);                         \
        repeat (i, y->n)                             \
            kF(r)[i] = func(x_f, item_F(y, i));      \
        return r;                                    \
    }                                                \
    check_type(is_F(x),);                            \
    check_length(xn == y->n,);                       \
    K r = ktn(KF, y->n);                             \
    repeat (i, y->n)                                 \
        kF(r)[i] = func(item_F(x, i), item_F(y, i)); \
    return r;                                        \
}

#define wrap_FF_(func) \
    wrap_FF(func, func)


// Wrap a function of (atom I; atom I; atom/vector F)
#define wrap_iiF(name, func)                 \
K DLL_EXPORT                                 \
qml_##name(K x, K y, K z) {                  \
    check_type(is_i(x) && is_i(y),);         \
    I x_i = as_i(x), y_i = as_i(y);          \
    wrap_last(z, func(x_i, y_i, v));         \
}


// Wrap a function of (atom F; atom F; atom/vector F)
#define wrap_ffF(name, func)                 \
K DLL_EXPORT                                 \
qml_##name(K x, K y, K z) {                  \
    check_type(is_f(x) && is_f(y),);         \
    F x_f = as_f(x), y_f = as_f(y);          \
    wrap_last(z, func(x_f, y_f, v));         \
}


// libm functions
wrap_F_(cos)
wrap_F_(sin)
wrap_F_(tan)
wrap_F_(acos)
wrap_F_(asin)
wrap_F_(atan)
wrap_FF_(atan2)
wrap_F_(cosh)
wrap_F_(sinh)
wrap_F_(tanh)
wrap_F_(acosh)
wrap_F_(asinh)
wrap_F_(atanh)
wrap_F_(exp)
wrap_F_(log)
wrap_F_(log10)
wrap_F_(logb)
wrap_F_(expm1)
wrap_F_(log1p)
wrap_FF_(pow)
wrap_F_(floor)
wrap_F_(ceil)
wrap_F_(fabs)
wrap_FF_(fmod)
wrap_F_(erf)
wrap_F_(erfc)
wrap_F_(lgamma)
wrap_F(gamma, tgamma)
wrap_F_(j0)
wrap_F_(j1)
wrap_F_(y0)
wrap_F_(y1)
wrap_F_(sqrt)
wrap_F_(cbrt)
wrap_FF_(hypot)


static I
sgamma(F x) {
    if (x > 0)
        return 1;
    return fmod(x, 2) > -1 ? -1 : 1;
}

static F
beta(F x, F y) {
    return sgamma(x) * sgamma(y) * sgamma(x + y) *
       exp(lgamma(x) + lgamma(y) - lgamma(x + y));
}

wrap_FF_(beta)


// Cephes probability functions
double igam(double, double);
double igamc(double, double);
double igami(double, double);
double incbet(double, double, double);
double incbi(double, double, double);
double ndtr(double);
double ndtri(double);
double stdtr(int, double);
double stdtri(int, double);
double fdtr(int, int, double);
double fdtrc(int, int, double);
double fdtri(int, int, double);
double chdtrc(double, double);
double chdtr(double, double);
double chdtri(double, double);
double bdtr(int, int, double);
double bdtri(int, int, double);
double pdtr(int, double);
double pdtri(int, double);
double smirnov(int, double);
double kolmogorov(double);
double smirnovi(int, double);
double kolmogi(double);


// chdtri() returns the inverse of the complementary CDF
static F
c2icdf(F x, F y) {
    return chdtri(x, 1 - y);
}

// fdtri() returns the inverse of the complementary CDF
static F
ficdf(I x, I y, F z) {
    return fdtri(x, y, 1 - z);
}

// gdtr() interprets second parameter differently
static F
gcdf(F x, F y, F z) {
    return igam(x, z / y);
}

// no gdtri() function
static F
gicdf(F x, F y, F z) {
    return y * igami(x, 1 - z);
}

// smirnov() returns complementary CDF
static F
smcdf(I x, F y) {
    return 1 - smirnov(x, y);
}

// smirnovi() returns inverse of the complementary CDF
static F
smicdf(I x, F y) {
    return smirnovi(x, 1 - y);
}

// komogorov() returns complementary CDF
static F
kcdf(F x) {
    return 1 - kolmogorov(x);
}

// komogoi() returns inverse of the complementary CDF
static F
kicdf(F x) {
    check(x>=1e-8, nf,); // doesn't converge well for small values
    return kolmogi(1 - x);
}

wrap_fF(pgammar, igam)
wrap_fF(pgammarc, igamc)
wrap_fF(ipgammarc, igami)
wrap_ffF(pbetar,incbet)
wrap_ffF(ipbetar,incbi)
wrap_F(ncdf, ndtr)
wrap_F(nicdf, ndtri)
wrap_iF(stcdf, stdtr)
wrap_iF(sticdf, stdtri)
wrap_iiF(fcdf, fdtr)
wrap_iiF(ficdf, ficdf)
wrap_fF(c2cdf, chdtr)
wrap_fF(c2icdf, c2icdf)
wrap_ffF(gcdf, gcdf)
wrap_ffF(gicdf, gicdf)
wrap_iiF(bncdf, bdtr)
wrap_iiF(bnicdf, bdtri)
wrap_iF(pscdf, pdtr)
wrap_iF(psicdf, pdtri)
wrap_iF(smcdf, smcdf)
wrap_iF(smicdf, smicdf)
wrap_F_(kcdf)
wrap_F_(kicdf)


//
// Options processing
//

struct opt {
    const char* s;
    I t; // 0 for boolean, -KI, or -KF
};

union optv {
    I i; /* first member */
    F f;
};

static int
take_opt_get(K x, I i, int flat, const S s,
             const struct opt* opt, union optv* v)
{
    check(s, 0,);
    check(*s, 1,);

    I k;
    for (k = 0; opt[k].s; k++)
        if (!strcmp(s, opt[k].s))
            goto found;
    return 0;

found:
    if (!opt[k].t && flat) {
        v[k].i = 1;
        return 1;
    } else {
        check(x && xt>=0 && i < xn, 0,);
        if (opt[k].t==-KF) {
            if (!xt) {
                check(is_f(xK[i]), 0,);
                v[k].f = as_f(xK[i]);
            } else {
                check(is_F(x), 0,);
                v[k].f = item_F(x, i);
            }
        } else {
            if (!xt) {
                check(is_i(xK[i]), 0,);
                v[k].i = as_i(xK[i]);
            } else {
                check(is_I(x), 0,);
                v[k].i = item_I(x, i);
            }
        }
        return 2;
    }
}

// Fills in an optv array
static int // 1 on success
take_opt(K x, const struct opt* opt, union optv* v) {
    switch (xt) {
    case 0:
        repeat (i, xn) {
            S s = xK[i]->t==-KS ? xK[i]->s : NULL;
            int r = take_opt_get(x, i+1, 1, s, opt, v);
            check(r, 0,);
            i += r-1;
        }
        break;
    case -KS:
    case KS:
        repeat (i, xt>0 ? xn    : 1) {
            S s =  xt>0 ? xS[i] : xs;
            check(take_opt_get(NULL, 0, 1, s, opt, v), 0,);
        }
        break;
    case XD:
        check(xx->t==KS || !xx->t, 0,);
        repeat (i, xx->n) {
            S s = xx->t==KS         ? kS(xx)[i] :
                  kK(xx)[i]->t==-KS ? kK(xx)[i]->s : NULL;
            check(take_opt_get(xy, i, 0, s, opt, v), 0,);
        }
    }
    return 1;
}


//
// Matrix functions
//

// Fortran BLAS prototypes
int dgemv_(char* trans, int* m, int* n,
           double* alpha, double* a, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
int dgemm_(char* transa, char* transb, int* m, int* n, int* k,
           double* alpha, double* a, int* lda,
           double* b, int* ldb, double* beta, double* c, int* ldc);
int dtrsv_(char* uplo, char* trans, char* diag, int* n,
           double* a, int* lda, double* x, int* incx);
int dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

// Fortran LAPACK prototypes
int dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
int dgetri_(int* n, double* a, int* lda, int* ipiv,
            double* work, int* lwork, int* info);
int dgeqrf_(int* m, int* n, double* a, int* lda,
            double* tau, double* work, int* lwork, int* info);
int dgeqp3_(int* m, int* n, double* a, int* lda, int* jpvt,
            double* tau, double* work, int* lwork, int* info);
int dorgqr_(int* m, int* n, int* k, double* a, int* lda,
            double* tau, double* work, int* lwork, int* info);
int dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
int dgesdd_(char* jobz, int* m, int* n, double* a, int* lda,
            double* s, double* u, int* ldu, double* vt, int* ldvt,
            double* work, int* lwork, int* iwork, int* info);
int dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
           double* b, int* ldb, int* info);
int dgesvx_(char* fact, char* trans, int* n, int* nrhs,
            double* a, int* lda, double* af, int* ldaf, int* ipiv,
            char* equed, double* r, double* c,
            double* b, int* ldb, double* x, int* ldx,
            double* rcond, double* ferr, double* berr,
            double* work, int* iwork, int* info);
int dgels_(char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
           double* b, int* ldb, double* work, int* lwork, int* info);
int dgelsd_(int* m, int* n, int* nrhs,
            double* a, int* lda, double* b, int* ldb,
            double* s, double* rcond, int* rank,
            double* work, int* lwork, int *iwork, int* info);
int dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda,
           double* wr, double* wi_,
           double* vl, int* ldvl, double* vr, int* ldvr,
           double* work, int* lwork, int* info);
int zgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda,
           double* w,
           double* vl, int* ldvl, double* vr, int* ldvr,
           double* work, int* lwork, double* rwork, int* info);

// Error handler for the Fortran BLAS interface, which is presumably only
// invoked if invalid parameters are specified by us or LAPACK.
int
xerbla_() {
    abort();
}

// Error handlers for the C BLAS interface, which are presumably only invoked if
// we specify invalid parameters.
void
ATL_xerbla(int p, char* rout, char* form, ...) {
    abort();
}

void
cblas_xerbla(int p, char* rout, char* form, ...) {
    abort();
}

// This should only be called before a logic error happens, but we quiet it to
// avoid linking with printf(). The original routine also involves error code
// arithmetic; we just preemptively abort to avoid reproducing it.
int
cblas_errprn(int ierr, int info, char* form, ...) {
    abort();
}


static K // returns error object
take_square_matrix(K x, F** r_, I* n_, int* triangular) {
    I n;
    F* r;

    check_type(!xt && (n = xn),);
    repeat (j, n) {
        check_type(is_F(xK[j]),);
        check_length(xK[j]->n == n,);
    }

    alloc(r, add_size(0, n, n),);

    if (triangular) {
        int upper = 1, lower = 1;
        repeat (j, n)
            repeat (i, n) {
                F f = item_F(xK[j], i);
                if (f != 0)
                    if (i < j) // non-zero below diagonal, so not upper
                        upper = 0;
                    else if (i > j) // non-zero above diagonal, so not lower
                        lower = 0;
                r[j + i*n] = f;
            }
        *triangular = upper ? 1 : lower ? -1 : 0;
    } else
        repeat (j, n)
            copy_F(xK[j], r + j, n);

    *r_ = r;
    *n_ = n;
    return no_error;
}


static int matrix_alloc_square_; // flag
#define matrix_alloc_square (&matrix_alloc_square_)

// ldr == m is allowed, otherwise *ldr must be set on input
static K // returns error object
take_matrix(K x, F** r, I* ldr, I* m, I* n, int* column) {
    if (column && column != matrix_alloc_square) {
        if (is_F(x)) {
            *m = xn; *n = 1;
            *ldr = max_i(*ldr, *m); // if ldr == m, *ldr is set above
            alloc(*r, add_size(0, *ldr, *n),);
            copy_F(x, *r, 1);
            *column = 1;
            return no_error;
        }
        *column = 0;
    }

    check_type(!xt && (*m = xn) && is_F(xK[0]) && (*n = xK[0]->n),);
    repeat (j, *m) {
        check_type(is_F(xK[j]),);
        check_length(*n == xK[j]->n,);
    }
    *ldr = max_i(*ldr, *m); // if ldr == m, *ldr is set above

    I n1 = column == matrix_alloc_square ? max_i(*m, *n) : *n;
    alloc(*r, add_size(0, *ldr, n1),);
    repeat (j, *m)
        copy_F(xK[j], *r + j, *ldr);
    return no_error;
}


static K
make_vector(const F* r, I offset, I n, I step) {
    K x = ktn(KF, n);
    if (!r)
        repeat (i, n)
            xF[i] = nf;
    else if (step == 1)
        memcpy(xF, r + offset, n * sizeof *r);
    else
        repeat (i, n)
            xF[i] = r[offset + i * step];
    return x;
}


#define matrix_transpose (-1)

// r can be null
static K
make_matrix(const F* r, I ldr, I m, I n, int column) {
    if (column && column != matrix_transpose)
        return make_vector(r, 0, m, 1);

    K x;
    if (column == matrix_transpose) {
        x = ktn(0, n);
        repeat (j, n)
            xK[j] = make_vector(r, j * ldr, m, 1);
    } else {
        x = ktn(0, m);
        repeat (j, m)
            xK[j] = make_vector(r, j, n, ldr);
    }
    return x;
}

// r can be null
static K
make_upper_matrix(const F* r, I ldr, I m, I n) {
    K x = ktn(0, m);
    repeat (j, m) {
        F* q = kF(xK[j] = ktn(KF, n));
        repeat (i, n)
            q[i] = i < j ? 0 : r ? r[j + i * ldr] : nf;
    }
    return x;
}


static K
make_complex(F a, F b) {
    if (b == 0) // in case isnan(b), return complex pair
        return kf(a);
    else {
        K x = ktn(KF, 2);
        xF[0] = a;
        xF[1] = b;
        return x;
    }
}

// a and b can be null
static K // xt==0 when complex
make_complex_vector(const F* a, const F* b, int n, int step) {
    if (b)
        repeat (i, n)
            if (!(b[i] == 0))
                goto complex;

    return make_vector(a, 0, n, step);

complex:;
    K x = ktn(0, n);
    repeat (i, n)
        xK[i] = make_complex(a[i * step], b[i * step]);
    return x;
}


// Matrix determinant
K DLL_EXPORT
qml_mdet(K x) {
    int triangular;
    I n, info;
    F *a, r = 1;

    bubble_error(take_square_matrix(x, &a, &n, &triangular),                  );

    if (!triangular) {
        I* ipiv;
        alloc(ipiv, n,                                                 free(a));
        dgetrf_(&n, &n, a, &n, ipiv, &info);
        check_lapack_return(info,                          free(ipiv); free(a));
        check(!info, kf(0),                                free(ipiv); free(a));

        repeat (i, n)
            if (ipiv[i]-1 != i)
                r = -r;
        /*                                              */ free(ipiv);
    }

    repeat (i, n)
        r *= a[i + i*n];
    /*                                                               */ free(a);
    return kf(r);
}


// Matrix inverse
K DLL_EXPORT
qml_minv(K x) {
    I n, info, lwork, *ipiv;
    F *a, *w, maxwork;

    bubble_error(take_square_matrix(x, &a, &n, NULL),                         );
    alloc(ipiv, n,                                                     free(a));

    dgetrf_(&n, &n, a, &n, ipiv, &info);
    check_lapack_return(info,                              free(ipiv); free(a));

    if (!info) {
        lwork = -1;
        dgetri_(&n, NULL, &n, NULL, &maxwork, &lwork, &info);
        check_lapack_return(info,                          free(ipiv); free(a));

        lwork = take_maxwork(maxwork);
        alloc(w, lwork,                                    free(ipiv); free(a));

        dgetri_(&n, a, &n, ipiv, w, &lwork, &info);
        /*                                     */ free(w);
        check_lapack_return(info,                          free(ipiv); free(a));
    }

    /*                                                  */ free(ipiv);
    x = make_matrix(info ? NULL : a, n, n, n, 0);
    /*                                                              */ free(a);
    return x;
}


// Matrix multiply
K DLL_EXPORT
qml_mm(K x, K y) {
    int b_column;
    I a_m, a_n, b_m, b_n;
    F *a, *b, *r;
    int i_1 = 1;
    double f_0 = 0, f_1 = 1;

    bubble_error(take_matrix(x, &a, &a_m, &a_m, &a_n, NULL),                  );
    bubble_error(take_matrix(y, &b, &b_m, &b_m, &b_n, &b_column),      free(a));
    check_length(a_n == b_m,                                  free(b); free(a));
    alloc(r, add_size(0, a_m, b_n),                           free(b); free(a));

    if (b_n == 1)
        dgemv_("N", &a_m, &a_n, &f_1, a, &a_m, b, &i_1, &f_0, r, &i_1);
    else if (a_m == 1)
        dgemv_("T", &b_m, &b_n, &f_1, b, &b_m, a, &i_1, &f_0, r, &i_1);
    else
        dgemm_("N", "N", &a_m, &b_n, &a_n,
               &f_1, a, &a_m, b, &b_m, &f_0, r, &a_m);

    /*                                                     */ free(b); free(a);
    x = make_matrix(r, a_m, a_m, b_n, b_column);
    /*                                            */ free(r);
    return x;
}


// Matrix substitution solve
K DLL_EXPORT
qml_ms(K x, K y) {
    int a_triangular, b_column;
    I a_n, b_m, b_n, info;
    F *a, *b;
    int i_1 = 1;
    double f_1 = 1;

    bubble_error(take_square_matrix(x, &a, &a_n, &a_triangular),              );
    check(a_triangular, krr(ss("domain")),                             free(a));
    bubble_error(take_matrix(y, &b, &b_m, &b_m, &b_n, &b_column),      free(a));
    check_length(a_n == b_m,                                  free(b); free(a));

    for (info = a_n; info; info--)
        if (a[(info-1) + (info-1)*a_n] == 0)
            goto done;

    if (b_n == 1)
        dtrsv_(a_triangular > 0 ? "U" : "L", "N", "N", &a_n, a, &a_n, b, &i_1); 
    else
        dtrsm_("L", a_triangular > 0 ? "U" : "L", "N", "N",
               &b_m, &b_n, &f_1, a, &a_n, b, &b_m);

done:
    /*                                                              */ free(a);
    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column);
    /*                                                     */ free(b);
    return x;
}


// Eigenvalues and eigenvectors
K DLL_EXPORT
qml_mevu(K x) {
    I n, info, lwork;
    F *a, *w, maxwork;

    bubble_error(take_square_matrix(x, &a, &n, NULL),                         );

    lwork = -1;
    dgeev_("N", "V", &n, NULL, &n, NULL, NULL,
           NULL, &n, NULL, &n, &maxwork, &lwork, &info);
    check_lapack_return(info,                                          free(a));

    lwork = take_maxwork(maxwork);
    alloc(w, add_size(add_size(lwork, n, 2), n, n),                    free(a));
    F *lr  = w  + lwork,
      *li  = lr + n,
      *ev  = li + n;
   /* *end = ev + n*n */

    dgeev_("N", "V", &n, a, &n, lr, li,
           NULL, &n, ev, &n, w, &lwork, &info);
    /*                                                              */ free(a);
    check_lapack_return(info,                                 free(w);        );

    x = ktn(0, 2);
    xK[0] = make_complex_vector(info ? NULL : lr, info ? NULL : li, n, 1);

    if (xK[0]->t) // no complex elements (also when info)
        xK[1] = make_matrix(info ? NULL : ev, n, n, n, matrix_transpose);
    else {
        xK[1] = ktn(0, n);
        repeat (j, n) {
            if (li[j] != 0 && j+1 < n) {
                K* q1 = kK(kK(xK[1])[j]   = ktn(0, n));
                K* q2 = kK(kK(xK[1])[j+1] = ktn(0, n));
                repeat (i, n) {
                    F a = ev[i + j*n], b = ev[i + (j+1)*n];
                    q1[i] = make_complex(a,  b);
                    q2[i] = make_complex(a, -b);
                }
                j++;
            } else
                kK(xK[1])[j] = make_vector(ev, j*n, n, 1);
        }
    }

    /*                                                     */ free(w);
    return x;
}


// Cholesky decomposition
K DLL_EXPORT
qml_mchol(K x) {
    I n, info;
    F* a;

    bubble_error(take_square_matrix(x, &a, &n, NULL),                         );

    dpotrf_("U", &n, a, &n, &info);
    check_lapack_return(info,                                          free(a));

    x = make_upper_matrix(info ? NULL : a, n, n, n);
    /*                                                              */ free(a);
    return x;
}


// QR factorization
static K
mqr(K x, int pivot) {
    I n, m, min, info, lwork;
    F *a, *w, maxwork[2];

    bubble_error(take_matrix(x, &a, &m, &m, &n, matrix_alloc_square),         );

    lwork = -1;
    if (pivot)
        dgeqp3_(&m, &n, NULL, &m, NULL, NULL, &maxwork[0], &lwork, &info);
    else
        dgeqrf_(&m, &n, NULL, &m, NULL,       &maxwork[0], &lwork, &info);
    check_lapack_return(info,                                          free(a));

    min = min_i(m, n);
    dorgqr_(&m, &m, &min, NULL, &m, NULL, &maxwork[1], &lwork, &info);
    check_lapack_return(info,                                          free(a));

    lwork = max_i(take_maxwork(maxwork[0]), take_maxwork(maxwork[1]));
    alloc(w, add_size(lwork, min, 1),                                  free(a));
    F* tau = w + lwork;

    K ipiv;
    if (pivot) {
        I* q = kI(ipiv = ktn(KI, n));
        repeat (i, n)
            q[i] = 0;
        dgeqp3_(&m, &n, a, &m, q, tau, w, &lwork, &info);
        repeat (i, n)
            q[i]--;
    } else
        dgeqrf_(&m, &n, a, &m,    tau, w, &lwork, &info);
    check_lapack_return(info,            if (pivot) r0(ipiv); free(w); free(a));

    K r = make_upper_matrix(a, m, m, n);

    dorgqr_(&m, &m, &min, a, &m, tau, w, &lwork, &info);
    /*                                                     */ free(w);
    check_lapack_return(info,     r0(r); if (pivot) r0(ipiv);          free(a));

    x = ktn(0, pivot ? 3 : 2);
    xK[0] = make_matrix(a, m, m, m, 0);
    xK[1] = r;
    if (pivot)
        xK[2] = ipiv;
    /*                                                              */ free(a);
    return x;
}

K DLL_EXPORT
qml_mqr(K x) {
    return mqr(x, 0);
}

K DLL_EXPORT
qml_mqrp(K x) {
    return mqr(x, 1);
}


// LUP factorization
K DLL_EXPORT
qml_mlup(K x) {
    I m, n, min, info, *ipiv;
    F* a;

    bubble_error(take_matrix(x, &a, &m, &m, &n, NULL),                        );

    min = min_i(m, n);
    alloc(ipiv, min,                                                   free(a));

    dgetrf_(&m, &n, a, &m, ipiv, &info);
    check_lapack_return(info,                              free(ipiv); free(a));

    x = ktn(0, 3);
    xK[0] = ktn(0, m);
    repeat (j, m) {
        F* q = kF(kK(xK[0])[j] = ktn(KF, min));
        repeat (i, min)
            q[i] = i < j ? a[j + i*m] : i == j ? 1 : 0;
    }
    xK[1] = make_upper_matrix(a, m, min, n);
    I* q = kI(xK[2] = ktn(KI, m));
    repeat (i, m)
        q[i] = i;
    repeat (i, min)
        swap_i(q + i, q + ipiv[i]-1);

    /*                                                  */ free(ipiv); free(a);
    return x;
}


// Singular value decomposition
K DLL_EXPORT
qml_msvd(K x) {
    I m, n, min, lwork, info, *iw;
    F *a, *w, maxwork;

    bubble_error(take_matrix(x, &a, &m, &m, &n, NULL),                        );

    min = min_i(m, n);
    lwork = -1;
    dgesdd_("A", &m, &n, NULL, &m, NULL, NULL, &m, NULL, &n,
            &maxwork, &lwork, NULL, &info);
    check_lapack_return(info,                                          free(a));

    lwork = take_maxwork(maxwork);
    alloc(w, add_size(add_size(add_size(lwork, min, 1), m, m), n, n),  free(a));
    F *s   = w  + lwork,
      *u   = s  + min,
      *vt  = u  + m * m;
   /* *end = vt + n * n */

    alloc(iw, add_size(0, min, 8),                            free(w); free(a));

    dgesdd_("A", &m, &n, a, &m, s, u, &m, vt, &n, w, &lwork, iw, &info);
    /*                                           */ free(iw);          free(a);
    check_lapack_return(info,                                 free(w);        );

    x = ktn(0, 3);
    xK[0] = make_matrix(info ? NULL : u, m, m, m, 0);
    xK[1] = ktn(0, m);
    repeat (j, m) {
        F* q = kF(kK(xK[1])[j] = ktn(KF, n));
        repeat (i, n)
            q[i] = i == j ? info ? nf : s[j] : 0;
    }
    xK[2] = make_matrix(info ? NULL : vt, n, n, n, matrix_transpose);

    /*                                                     */ free(w);
    return x;
}


//
// Advanced matrix functions
//

// Polynomial root finding
K DLL_EXPORT
qml_poly(K x) {
    int complex;
    I n, lwork, info;
    F lca, lcb, *w, maxwork[2];

    if (is_F(x))
        complex = 0;
    else {
        check_type(!xt,);
        repeat (i, xn)
            if (!is_f(xK[i])) {
                check_type(is_F(xK[i]),);
                check_length(xK[i]->n == 2,);
            }
        complex = 1;
    }
    check_length(n = xn,);
    n--;

    if (complex && !is_f(xK[0])) {
        lca = item_F(xK[0], 0);
        lcb = item_F(xK[0], 1);
    } else {
        lca = item_F(x, 0);
        lcb = 0;
    }

    check(lca != 0 || lcb != 0, krr(ss("roots")),);
    check(n, ktn(0, 0),);

    lwork = -1;
    if (complex)
        zgeev_("N", "N", &n, NULL, &n, NULL, NULL,       &n, NULL, &n,
               maxwork, &lwork, NULL, &info);
    else
        dgeev_("N", "N", &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n,
               maxwork, &lwork,       &info);
    check_lapack_return(info,);

    lwork = take_maxwork(maxwork[0]);
    alloc(w, add_size(0, add_size(add_size(lwork, n, 2), n, n), 1+complex),);
    F *l     = w     + lwork * (1+complex),
      *rwork = l     + n * 2,
      *cm    = rwork + n * (complex*2);

    /* make companion matrix */
    if (complex) {
        repeat (i, n-1)
            repeat (j, n*2)
                cm[j + i*n*2] = j == (i + 1)*2;
        F d = lca*lca + lcb*lcb;
        if (lcb == 0) {
            d = lca;
            lca = 1;
        }
        repeat (j, n) {
            F a, b;
            if (is_f(xK[n-j])) {
                a = as_f(xK[n-j]);
                b = 0;
            } else {
                a = item_F(xK[n-j], 0);
                b = item_F(xK[n-j], 1);
            }
            cm[j*2 + (n-1)*n*2]     = (a * lca + b * lcb) / -d;
            cm[j*2 + (n-1)*n*2 + 1] = (b * lca - a * lcb) / -d;
        }
    } else {
        repeat (i, n-1)
            repeat (j, n)
                cm[j + i*n] = j == i + 1;
        repeat (j, n)
            cm[j + (n-1)*n] = item_F(x, n-j) / -lca;
    }

    if (complex)
        zgeev_("N", "N", &n, cm, &n, l,
               NULL, &n, NULL, &n, w, &lwork, rwork, &info);
    else
        dgeev_("N", "N", &n, cm, &n, l, l+n,
               NULL, &n, NULL, &n, w, &lwork,        &info);
    check_lapack_return(info,                                          free(w));
    check(!info, krr(ss("roots")),                                     free(w));

    x = make_complex_vector(l, complex ? l+1 : l+n, n, 1+complex);
    /*                                                              */ free(w);
    return x;
}


// Linear equations
static K
mls(K x, K y, int equi) {
    char equed;
    int b_column;
    I a_n, b_m, b_n, info, *ipiv;
    F *a, *b;

    bubble_error(take_square_matrix(x, &a, &a_n, NULL),                       );
    bubble_error(take_matrix(y, &b, &b_m, &b_m, &b_n, &b_column),      free(a));
    check_length(a_n == b_m,                                  free(b); free(a));

    alloc(ipiv, add_size(0, a_n, 1 + equi),                   free(b); free(a));

    if (equi) {
        F* bi = b;
        alloc(b, add_size(add_size(add_size(add_size(0,
                 a_n, b_n), a_n, a_n), a_n, 6), b_n, 2),
                                                 free(ipiv); free(bi); free(a));
        F *af   = b    + a_n*b_n, 
          *r    = af   + a_n*a_n,
          *c    = r    + a_n,
          *w    = c    + a_n,
          *ferr = w    + a_n*4,
          *berr = ferr + b_n,
       /* *end  = berr + b_n, */
          rcond;
        I* iw = ipiv + a_n;

        dgesvx_("E", "N", &a_n, &b_n, a, &a_n, af, &a_n, ipiv, &equed, r, c,
                bi, &b_m, b, &b_m, &rcond, ferr, berr, w, iw, &info);
        /*                                                */ free(bi);
    } else
        dgesv_(&a_n, &b_n, a, &a_n, ipiv, b, &b_m, &info);
    /*                                         */ free(ipiv);          free(a);
    check_lapack_return(info,                                 free(b)         );

    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column);
    /*                                                     */ free(b);
    return x;
}

static const struct opt ls_opt[] = {
    { "equi", 0 },
    { NULL }
};

K DLL_EXPORT
qml_mlsx(K opts, K x, K y) {
    union optv v[] = { { 0 } };
    check(take_opt(opts, ls_opt, v), krr(ss("opt")),);
    return mls(x, y, v[0].i);
}

K DLL_EXPORT
qml_mls(K x, K y) {
    return mls(x, y, 0);
}


// Linear least squares
static K
mlsq(K x, K y, int svd) {
    int b_column;
    I a_m, a_n, b_m, b_n, ldb, rank, info, lwork, liwork, *iw;
    F *a, *b, *w, maxwork, rcond = -1;

    bubble_error(take_matrix(x, &a, &a_m, &a_m, &a_n, NULL),                  );
    ldb = a_n;
    bubble_error(take_matrix(y, &b, &ldb, &b_m, &b_n, &b_column),      free(a));
    check_length(a_m == b_m,                                  free(b); free(a));

    lwork = -1;
    if (svd)
        dgelsd_(&a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
                NULL, &rcond, &rank, &maxwork, &lwork, &liwork, &info);
    else
        dgels_("N", &a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
               &maxwork, &lwork, &info);
    check_lapack_return(info,                                 free(b); free(a));

    lwork = take_maxwork(maxwork);
    alloc(w, add_size(lwork, min_i(a_m, a_n), svd),           free(b); free(a));
    if (svd) {
        F* s = w + lwork;
        alloc(iw, liwork,                            free(w); free(b); free(a));
        dgelsd_(&a_m, &a_n, &b_n, a, &a_m, b, &ldb,
                s, &rcond, &rank, w, &lwork, iw, &info);
        /*                              */ free(iw);
    } else
        dgels_("N", &a_m, &a_n, &b_n, a, &a_m, b, &ldb,
               w, &lwork, &info);
    /*                                            */ free(w); /*    */ free(a);
    check_lapack_return(info,                                 free(b);        );

    x = make_matrix(info ? NULL : b, ldb, a_n, b_n, b_column);
    /*                                                     */ free(b);
    return x;
}

static const struct opt lsq_opt[] = {
    { "svd", 0 },
    { NULL }
};

K DLL_EXPORT
qml_mlsqx(K opts, K x, K y) {
    union optv v[] = { { 0 } };
    check(take_opt(opts, lsq_opt, v), krr(ss("opt")),);
    return mlsq(x, y, v[0].i);
}

K DLL_EXPORT
qml_mlsq(K x, K y) {
    return mlsq(x, y, 0);
}


//
// Optimization functions
//

int conmax_(int* ioptn, int* nparm, int* numgr, int* itlim,
            double* fun, int* ifun, double* pttbl,
            int* iwork, int* liwrk, double* work, int* lwrk,
            int* iter, double* param, double* error);
int muller_(int* limmul, int* nsrch, int* ioptn, int* nparm, int* numgr,
            double* dvec, double* fun, int* ifun, double* pttbl,
            double* zwork, double* tolcon, int* iphse,
            int* iwork, int* liwrk, double* work, int* lwrk, double* parwrk,
            double* err1, double* p1, double* f1, double* procor, double* emin);
int searsl_(int* initlm, int* nadd, int* lims1, int* ioptn,
            int* numgr, int* nparm, double* prjlim, double* tol1,
            double* x, double* fun, int* ifun, double* pttbl,
            double* param, double* error,
            double* rchdwn, int* mact, int* iact, int* iphse, double* unit,
            double* tolcon, double* rchin, int* itypm1, int* itypm2,
            int* iwork, int* liwrk, double* work, int* lwrk,
            double* err1, double* parprj, double* projct,
            double* emin, double* emin1, double* parser, int* nsrch);

double
d1mach_(int* i) {
    // This should only be called by CONMAX, and only with this argument.
    if (unlikely(*i != 3))
        abort();

    // DBL_EPSILON is pow(FLT_RADIX, 1-DBL_MANT_DIG),
    // CONMAX needs   pow(FLT_RADIX,  -DBL_MANT_DIG).
    return DBL_EPSILON / FLT_RADIX;
}


static K // returns error object
count_param(K x, I* n) {
    if (!xt)
        repeat (i, xn)
            bubble_error(count_param(xK[i], n),);
    else if (xt>0) {
        check_type(is_F(x),);
        *n = add_size(*n, xn, 1);
    } else {
        check_type(is_f(x),);
        *n = add_size(*n, 1, 1);
    }
    return no_error;
}

static F*
take_param(K x, F* param) {
    if (!xt)
        repeat (i, xn)
            param = take_param(xK[i], param);
    else if (xt>0) {
        copy_F(x, param, 1);
        param += xn;
    } else
        *param++ = as_f(x);
    return param;
}

static F*
make_param(K x, F* param, K* a) {
    if (!xt) {
        *a = ktn(0, xn);
        repeat (i, xn)
            param = make_param(xK[i], param, &kK(*a)[i]);
    } else if (xt>0) {
        *a = make_vector(param, 0, xn, 1);
        param += xn;
    } else
        *a = kf(*param++);
    return param;
}


// The first member must be of the same type as pttbl. This allows casting a
// structure pointer to a first-member pointer, passing it through CONMAX, and
// casting it back to get the original pointer in a well-defined manner.
struct fnset_info {
    F base; /* first member */
    K fun, con, start;
    I contyp;
    int neg;
    K error;
};

static F
fnset_call(struct fnset_info* info, K f, int neg, F* param) {
    check(info->error == no_error, 0,);

    K a;
    if (!info->start && info->con) /* line */
        a = knk(1, kf(info->base + (info->neg ? -*param : *param)));
    else if (!info->start || info->start->t<0) /* root, or scalar param */
        a = knk(1, kf(*param));
    else
        make_param(info->start, param, &a);

    K x = dot(f, a);
    r0(a);
    if (x && is_f(x)) {
        F v = as_f(x);
        r0(x);
        if (isnan(v)) // protect optimization routines from NaNs
            return wf; // this is -wf for ">=" constraints
        return neg ? -v : v;
    }

    // function didn't return a float as we'd hoped
    if (!x || x->t==-128)
        info->error = x;
    else {
        info->error = krr(x->t>=100 ? "rank" : "type");
        r0(x);
    }
    return 0;
}

int
fnset_(int* nparm, int* numgr, double* pttbl, double* param,
       int* ipt, int* indfn, int* icntyp, double* confun)
{
    struct fnset_info* info = (struct fnset_info*)pttbl;

    icntyp += *ipt - 1;
    confun += *ipt - 1;

    K f;
    int neg;
    if (!info->con) { /* root or solve */
        f = info->fun->t ? info->fun : kK(info->fun)[*ipt-1];
        neg = info->neg;
        *icntyp = 2;
    } else /* line, min or conmin */
        if (*ipt == 1) { /* objective function */
            f = info->fun;
            neg = 0;
            *icntyp = 1;
        } else { /* constraints */
            f = info->con->t ? info->con : kK(info->con)[*ipt-2];
            neg = 1;
            *icntyp = info->contyp;
        }

    F v = *confun = fnset_call(info, f, neg, param);
    if (*indfn) {
        I m = *numgr, n = *nparm;
        if (*icntyp == -1) /* linear function */
            repeat (i, n) {
                F p = param[i];
                param[i] = p + 1;
                *(confun += m) = fnset_call(info, f, neg, param) - v;
                param[i] = p;
            }
        else { /* nonlinear function */
            F h = sqrt(DBL_EPSILON / FLT_RADIX);
            repeat (i, n) {
                F p = param[i], p1, p2, v1, v2;
                param[i] = p1 = p + h; v1 = fnset_call(info, f, neg, param);
                param[i] = p2 = p - h; v2 = fnset_call(info, f, neg, param);
                *(confun += m) = (v1 - v2) / (p1 - p2);
                param[i] = p;
            }
        }
    }

    if (info->error != no_error)
        *icntyp = 0;
    return 0;
}


static K
solvemin(K fun, K con, K start, I maxiter, F tolcon, I steps,
         int slp, int rk, int lincon, int full, int quiet)
{
    struct fnset_info info;

    check_type(fun->t>=100 || !con && !fun->t,);
    check_type(!con || con->t>=100 || !con->t,);
    check(!slp || !rk && steps<=0, krr(ss("opt")),);

    I numgr, ifun;
    if (con) { /* min or conmin */
        numgr = 1 + (!con->t ? con->n : 1);
        ifun = 1;
    } else { /* solve */
        numgr = !fun->t ? fun->n : 1;
        ifun = numgr;
    }

    I nparm = 0;
    bubble_error(count_param(start, &nparm),);

    I lwrk = add_size(add_size(add_size(add_size(
        13, nparm, 27), numgr, 11), numgr, nparm*4), nparm, nparm*2);

    I twrk = add_size(add_size(add_size(
        lwrk, ifun, 1), nparm, 1), numgr + 3, 1);

    I liwrk = add_size(add_size(3, nparm, 7), numgr, 7);

    F* work;
    I* iwork;
    alloc(work, twrk,                                                         );
    alloc(iwork, liwrk,                                             free(work));

    F* error = work + lwrk;
    F* vfun  = error + numgr + 3; // vfun[-3] will refer to error[numgr]
    F* param = vfun + ifun;

    take_param(start, param);

    I ioptn = 200;
    I itlim = max_i(0, maxiter);
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    work[1] = tolcon;
    if (steps > 0) {
        iwork[1] = steps;
        ioptn += 100;
    }
    ioptn += slp ? 1000 : rk ? 2000 : 0;
    repeat (i, ifun)
        vfun[i] = 0;

    info.fun = fun;
    info.con = con;
    info.start = start;
    info.contyp = lincon ? -1 : -2;
    info.neg = 0;
    info.error = no_error;

    I iter;
    conmax_(&ioptn, &nparm, &numgr, &itlim, vfun, &ifun, (F*)&info,
            iwork, &liwrk, work, &lwrk, &iter, param, error);
    /*                                              */ free(iwork);
    bubble_error(info.error,                                        free(work));

    S sig = NULL;
    if (iter >= maxiter)
        sig = "iter";
    else if (iter < 0 || !con && !(vfun[-3] >= -tolcon && vfun[-3] <= tolcon))
        sig = "feas";
    else
        repeat (i, nparm)
            if (isnan(param[i])) {
                sig = "nan";
                break;
            }
    check(!sig || quiet, krr(ss(sig)),                              free(work));

    K x, y;
    if (sig) {
        if (full)
            make_param(start, param, &y);
        repeat (i, nparm)
            param[i] = nf;
    }
    make_param(start, param, &x);

    if (full) {
        /*          solve       min        conmin */
        /* normal: `x`iter      `x`f`iter  `x`f`cons`iter */
        /* error:  `x`last`err  `x`last`f  `x`last`f`cons */
        /*         `iter`sig    `iter`sig  `err`iter`sig */
        I j = 2 + (sig ? 2 + (con != empty_con) : 0) +
                  (con ? 1 + (con != empty_con) : 0);
        K k = ktn(KS, j), v = ktn(0, j);
        if (sig) kS(k)[--j] = ss("sig"),  kK(v)[j] = ks(ss(sig));
        /*    */ kS(k)[--j] = ss("iter"), kK(v)[j] = ki(iter >= 0 ? iter : 0);
        if (sig && con != empty_con) {
            F err = vfun[!con ? -3 : lincon ? -2 : -1];
            /**/ kS(k)[--j] = ss("err"),  kK(v)[j] = kf(err);
        }
        if (con) {
            if (con != empty_con) {
                 K cons = ktn(KF, numgr-1);
                 repeat (i, numgr-1) {
                     F e = error[i+1];
                     kF(cons)[i] = e >= -tolcon && e <= tolcon ? 0 : -e;
                 }
                 kS(k)[--j] = ss("cons"), kK(v)[j] = cons;
            }
            /**/ kS(k)[--j] = ss("f"),    kK(v)[j] = kf(vfun[-3]);
        }
        if (sig) kS(k)[--j] = ss("last"), kK(v)[j] = y;
        /*    */ kS(k)[--j] = ss("x"),    kK(v)[j] = x;
        x = xD(k, v);
    }
    /*                                                           */ free(work);
    return x;
}


static K
root(K fun, K x, I maxiter, F tolcon, int full, int quiet) {
    I iter, nsrch, ioptn, nparm, numgr, ifun, iphse, iwork[17], liwrk, lwrk;
    F dvec, cfun, zwork, work[6], parwrk, err1[4], p1, p2, f1, f2;
    struct fnset_info info;

    check_type(fun->t>=100 && is_F(x),);
    check_length(xn == 2,);
    p1 = item_F(x, 0);
    p2 = item_F(x, 1);
    if (p1 > p2)
        swap_f(&p1, &p2);

    info.fun = fun;
    info.con = NULL; /* root/solve flag */
    info.start = NULL; /* root/line flag */
    info.contyp = -2;
    info.neg = 0;
    info.error = no_error;

    f1 = fnset_call(&info, fun, 0, &p1);
    f2 = fnset_call(&info, fun, 0, &p2);
    bubble_error(info.error,);

    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    if (f1 < -tolcon && f2 > tolcon) {
        f1 = -f1;
        f2 = -f2;
        info.neg = 1;
    }
    int sig_sign = !(f1 > tolcon && f2 < -tolcon);

    iter = maxiter >= 0 ? maxiter : 0;
    ifun = dvec = numgr = nparm = 1;
    zwork = cfun = iphse = ioptn = 0;
    liwrk = 17;
    lwrk = 6;
    iwork[15] = -2;

    muller_(&iter, &nsrch, &ioptn, &nparm, &numgr,
            &dvec, &cfun, &ifun, (F*)&info,
            &zwork, &tolcon, &iphse,
            iwork, &liwrk, work, &lwrk, &parwrk,
            err1, &p1, &f1, &p2, &f2);
    bubble_error(info.error,);

    S sig = nsrch >= maxiter                 ? "iter" :
            !(f2 >= -tolcon && f2 <= tolcon) ? "feas" :
            isnan(p2)                        ? "nan"  : NULL;
    if (sig && sig_sign)
        sig = "sign";
    check(!sig || quiet, krr(ss(sig)),);

    x = sig ? kf(nf) : kf(p2);
    if (full) {
        /* normal: `x`iter */
        /* error:  `x`last`err`iter`sig */
        I j = 2 + (sig ? 3 : 0);
        K k = ktn(KS, j), v = ktn(0, j);
        if (sig) kS(k)[--j] = ss("sig"),  kK(v)[j] = ks(ss(sig));
        /*    */ kS(k)[--j] = ss("iter"), kK(v)[j] = ki(nsrch);
        if (sig) kS(k)[--j] = ss("err"),  kK(v)[j] = kf(fabs(f2));
        if (sig) kS(k)[--j] = ss("last"), kK(v)[j] = kf(p2);
        /*    */ kS(k)[--j] = ss("x"),    kK(v)[j] = x;
        x = xD(k, v);
    }
    return x;
}


static K
line(K fun, K base, K x, I maxiter, F tolcon, int full, int quiet) {
    I initlm, nadd, lims1, ioptn, numgr, nparm, ifun, mact, iact, iphse,
      itypm1, itypm2, iwork[17], liwrk, lwrk, nsrch;
    F prjlim, tol1, cx[2], cfun, param, error[4], rchdwn, unit, rchin,
      work[42], err1[4], parprj, projct, emin, emin1, parser;
    struct fnset_info info;

    check_type(fun->t>=100 && is_f(base) && is_f(x),);

    info.base = as_f(base);
    info.fun = fun;
    info.con = fun; /* min/line flag */
    info.start = NULL; /* root/line flag */
    info.neg = 0;
    info.error = no_error;

    projct = as_f(x) - info.base;
    if (projct < 0) {
        projct = -projct;
        info.neg = 1;
    }

    /* Call sequence based on sample driver from conmax.f. */

    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    prjlim = wf;
    tol1 = 100 * DBL_EPSILON;
    if (tol1 > tolcon)
        tol1 = tolcon;
    lims1 = nadd = initlm = maxiter >= 0 ? maxiter : 0;
    unit = iact = mact = ifun = nparm = numgr = 1;
    param = cfun = itypm2 = itypm1 = iphse = ioptn = 0;
    rchin = rchdwn = 2;
    liwrk = 17;
    lwrk = 42;
    cx[0] = 1;
    iwork[6] = 1;

    searsl_(&initlm, &nadd, &lims1, &ioptn,
            &numgr, &nparm, &prjlim, &tol1,
            cx, &cfun, &ifun, (F*)&info,
            &param, error,
            &rchdwn, &mact, &iact, &iphse, &unit,
            &tolcon, &rchin, &itypm1, &itypm2,
            iwork, &liwrk, work, &lwrk, err1, &parprj, &projct,
            &emin, &emin1, &parser, &nsrch);
    bubble_error(info.error,);

    S sig = nsrch >= maxiter || nsrch >= lims1 ? "iter" :
            isnan(projct)                      ? "nan"  : NULL;
    check(!sig || quiet, krr(ss(sig)),);

    projct = info.base + (info.neg ? -projct : projct);
    x = sig ? kf(nf) : kf(projct);
    if (full) {
        /* normal: `x`f`iter */
        /* error:  `x`last`f`iter`sig */
        I j = 3 + (sig ? 2 : 0);
        K k = ktn(KS, j), v = ktn(0, j);
        if (sig) kS(k)[--j] = ss("sig"),  kK(v)[j] = ks(ss(sig));
        /*    */ kS(k)[--j] = ss("iter"), kK(v)[j] = ki(nsrch);
        /*    */ kS(k)[--j] = ss("f"),    kK(v)[j] = kf(emin);
        if (sig) kS(k)[--j] = ss("last"), kK(v)[j] = kf(projct);
        /*    */ kS(k)[--j] = ss("x"),    kK(v)[j] = x;
        x = xD(k, v);
    }
    return x;
}


Z const struct opt solve_opt[] = {
    { "iter",   -KI },
    { "tol",    -KF },
    { "steps",  -KI },
    { "slp",      0 },
    { "rk",       0 },
    { "full",     0 },
    { "quiet",    0 },
    { NULL }
};

K DLL_EXPORT
qml_solvex(K opts, K x, K y) {
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 } };
    check(take_opt(opts, solve_opt, v), krr(ss("opt")),);
    return solvemin(x, NULL, y, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K DLL_EXPORT
qml_solve(K x, K y) {
    return qml_solvex(empty_con, x, y);
}

K DLL_EXPORT
qml_minx(K opts, K x, K y) {
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 } };
    check(take_opt(opts, solve_opt, v), krr(ss("opt")),);
    return solvemin(x, empty_con, y, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K DLL_EXPORT
qml_min(K x, K y) {
    return qml_minx(empty_con, x, y);
}

static const struct opt conmin_opt[] = {
    { "iter",   -KI },
    { "tol",    -KF },
    { "steps",  -KI },
    { "slp",      0 },
    { "rk",       0 },
    { "lincon",   0 },
    { "full",     0 },
    { "quiet",    0 },
    { NULL }
};

K DLL_EXPORT
qml_conminx(K opts, K x, K y, K z) {
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    check(take_opt(opts, conmin_opt, v), krr(ss("opt")),);
    return solvemin(x, y, z, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, v[5].i, v[6].i, v[7].i);
}

K DLL_EXPORT
qml_conmin(K x, K y, K z) {
    return qml_conminx(empty_con, x, y, z);
}

static const struct opt rootline_opt[] = {
    { "iter", -KI },
    { "tol",  -KF },
    { "full",   0 },
    { "quiet",  0 },
    { NULL }
};

K DLL_EXPORT
qml_rootx(K opts, K x, K y) {
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    check(take_opt(opts, rootline_opt, v), krr(ss("opt")),);
    return root(x, y, v[0].i, v[1].f, v[2].i, v[3].i);
}

K DLL_EXPORT
qml_root(K x, K y) {
    return qml_rootx(empty_con, x, y);
}

K DLL_EXPORT
qml_linex(K opts, K x, K y, K z) {
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    check(take_opt(opts, rootline_opt, v), krr(ss("opt")),);
    return line(x, y, z, v[0].i, v[1].f, v[2].i, v[3].i);
}

K DLL_EXPORT
qml_line(K x, K y, K z) {
    return qml_linex(empty_con, x, y, z);
}


//
// Constants
//

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

K DLL_EXPORT
qml_const(K x) {
    check_type(xt==-KI,);
    switch (xi) {
    case 0:
        return ks(QUOTE(QML_VERSION));
    case 1:
        return kf(3.1415926535897932384626433832795028842);
    case 2:
        return kf(2.7182818284590452353602874713526624978);
    case 3:
        return kf(DBL_EPSILON);
    }
    return krr(ss("domain"));
}
