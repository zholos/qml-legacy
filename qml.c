#include <config.h>

#undef QML_DLLMAIN
#if defined(QML_DLLEXPORT) && !defined(__GNUC__)
    #define QML_DLLMAIN
#endif

#include <string.h>
#include <float.h>

#include <fdlibm.h>

#ifdef QML_DLLMAIN
    #include <windows.h>
#endif

#undef min
#undef max
#undef small
#undef large
#include <f2c.h>
#include <clapack.h>

#include <k.h>

#undef QML_EXPORT
#ifdef QML_DLLEXPORT
    #define QML_EXPORT __declspec(dllexport)
#elif __GNUC__>=4
    #define QML_EXPORT __attribute__((visibility("default")))
#else
    #define QML_EXPORT
#endif


// Return 1 if argument can be converted to I atom
ZI mi1p(K x) {
    R xt==-KB || xt<=-KG && xt>=-KI;
}

// Return I value for argument
ZI mi1i(K x) {
    SW(xt) {
        CS(-KI, R xi==ni ? ni : xi==wi ? wi : xi==-wi ? -wi : xi)
        CS(-KH, R xh==nh ? ni : xh==wh ? wi : xh==-wh ? -wi : xh)
        case -KB:
        CS(-KG, R xg)
        CD: R ni;
    }
}

// Make new I atom out of B, G, H or I argument
Z K1(mi1) {
    U(mi1p(x));
    R ki(mi1i(x));
}

// Make new I atom or vector out of B, G, H or I argument
Z K1(mi) {
    K n;
    P(xt<0, mi1(x));
    U(xt>=KG && xt<=KI || xt==KB || xt==0)
    n = ktn(KI, xn);
    SW(xt) {
        CS(KI, memcpy(kI(n), xI, xn*sizeof(I)))
        CS(KH, DO(xn, kI(n)[i] = xH[i]==nh ? ni :
            xH[i]==wh ? wi : xH[i]==-wh ? -wi : xH[i]))
        case KB:
        CS(KG, DO(xn, kI(n)[i] = xG[i]))
        CS(0, DO(xn, if (!mi1p(xK[i])) { r0(n); R 0; } kI(n)[i] = mi1i(xK[i])))
    }
    R n;
}

// Return 1 if argument can be converted to F atom
ZI mf1p(K x) {
    R xt==-KB || xt<=-KG && xt >= -KF;
}

// Return F value for argument
ZF mf1f(K x) {
    SW(xt) {
        CS(-KF, R xf)
        CS(-KE, R xe)
        CS(-KJ, R xj==nj ? nf : xj==wj ? wf : xj==-wj ? -wf : xj)
        CS(-KI, R xi==ni ? nf : xi==wi ? wf : xi==-wi ? -wf : xi)
        CS(-KH, R xh==nh ? nf : xh==wh ? wf : xh==-wh ? -wf : xh)
        case -KB:
        CS(-KG, R xg)
        CD: R nf;
    }
}

// Make new F atom out of B, G, H, I, J, E or F argument
Z K1(mf1) {
    U(mf1p(x));
    R kf(mf1f(x));
}

// Make new F atom or vector out of B, G, H, I, J, E or F argument
Z K1(mf) {
    K f;
    P(xt<0, mf1(x));
    U(xt>=KG && xt<=KF || xt==KB || xt==0)
    f = ktn(KF, xn);
    SW(xt) {
        CS(KF, memcpy(kF(f), xF, xn*sizeof(F)))
        CS(KE, DO(xn, kF(f)[i] = xE[i]))
        CS(KJ, DO(xn, kF(f)[i] = xJ[i]==nj ? nf :
            xJ[i]==wj ? wf : xJ[i]==-wj ? -wf : xJ[i]))
        CS(KI, DO(xn, kF(f)[i] = xI[i]==ni ? nf :
            xI[i]==wi ? wf : xI[i]==-wi ? -wf : xI[i]))
        CS(KH, DO(xn, kF(f)[i] = xH[i]==nh ? nf :
            xH[i]==wh ? wf : xH[i]==-wh ? -wf : xH[i]))
        case KB:
        CS(KG, DO(xn, kF(f)[i] = xG[i]))
        CS(0, DO(xn, if (!mf1p(xK[i])) { r0(f); R 0; } kF(f)[i] = mf1f(xK[i])))
    }
    R f;
}

// Wrap a function of F
#define WRAPf(ff) WRAPnf(ff,ff)
#define WRAPnf(ffn,ff) \
K QML_EXPORT qml_##ffn(K u) { \
    K x = mf(u); \
    P(x==0, krr(ss("type"))); \
    SW(xt) { \
        CS(-KF, xf = ff(xf)) \
        CS(KF, DO(xn, xF[i] = ff(xF[i]))) \
    } \
    R x; \
}

// Wrap a function of (F; F)
#define WRAPff(ff) \
K QML_EXPORT qml_##ff(K u, K v) { \
    K x, y; \
    P(!(u->t<0 || v->t<0 || u->n==v->n), krr(ss("length"))) \
    x = mf(u); P(x==0, krr(ss("type"))) \
    y = mf(v); if (y==0) { r0(x); R krr(ss("type")); } \
    SW(xt) { \
        case -KF: SW(y->t) { \
            CS(-KF, xf = ff(xf, y->f)) \
            case KF: DO(y->n, kF(y)[i] = ff(xf, kF(y)[i])) r0(x); R y; \
        } break; \
        case KF: SW(y->t) { \
            CS(-KF, DO(x->n, xF[i] = ff(xF[i], y->f))) \
            CS(KF, DO(x->n, xF[i] = ff(xF[i], kF(y)[i]))) \
        } \
    } \
    r0(y); R x; \
}

// Wrap a function of (F; F)
#define WRAPnff(ffn,ff) \
K QML_EXPORT qml_##ffn(K u, K v) { \
    P(!mf1p(u) || !mf1p(v), krr(ss("type"))) \
    R kf(ff(mf1f(u), mf1f(v))); \
}

// Wrap a function of (I; F)
#define WRAPnif(ffn,ff) \
K QML_EXPORT qml_##ffn(K u, K v) { \
    P(!mi1p(u) || !mf1p(v), krr(ss("type"))) \
    R kf(ff(mi1i(u), mf1f(v))); \
}

// Wrap a function of (F; F; F)
#define WRAPnfff(ffn,ff) \
K QML_EXPORT qml_##ffn(K u, K v, K w) { \
    P(!mf1p(u) || !mf1p(v) || !mf1p(w), krr(ss("type"))) \
    R kf(ff(mf1f(u), mf1f(v), mf1f(w))); \
}

// Wrap a function of (I; I; F)
#define WRAPniif(ffn,ff) \
K QML_EXPORT qml_##ffn(K u, K v, K w) { \
    P(!mi1p(u) || !mi1p(v) || !mf1p(w), krr(ss("type"))) \
    R kf(ff(mi1i(u), mi1i(v), mf1f(w))); \
}

// fdlibmm functions
WRAPf(cos)WRAPf(sin)WRAPf(tan)WRAPf(acos)WRAPf(asin)WRAPf(atan)WRAPff(atan2)
WRAPf(cosh)WRAPf(sinh)WRAPf(tanh)WRAPf(acosh)WRAPf(asinh)WRAPf(atanh)
WRAPf(exp)WRAPf(log)WRAPf(log10)WRAPf(logb)WRAPf(expm1)WRAPf(log1p)WRAPff(pow)
WRAPf(floor)WRAPf(ceil)WRAPf(fabs)WRAPff(fmod)
WRAPf(erf)WRAPf(erfc)WRAPf(lgamma)
WRAPf(j0)WRAPf(j1)WRAPf(y0)WRAPf(y1)
WRAPf(sqrt)WRAPf(cbrt)WRAPff(hypot)

// Cephes probability functions
extern double gamma(double);
extern double igam(double, double);
extern double igamc(double, double);
extern double igami(double, double);
extern double incbet(double, double, double);
extern double incbi(double, double, double);
extern double ndtr(double);
extern double ndtri(double);
extern double stdtr(int, double);
extern double stdtri(int, double);
extern double fdtr(int, int, double);
extern double fdtrc(int, int, double);
extern double fdtri(int, int, double);
extern double chdtrc(double, double);
extern double chdtr(double, double);
extern double chdtri(double, double);

// this is not available in all distributions of Cephes so we provide it here
F beta(F x, F y) {
    I sx, sy, sxy;
    F r;
    r = lgamma_r(x, &sx);
    r += lgamma_r(y, &sy);
    r -= lgamma_r(x + y, &sxy);
    R sx*sy*sxy*exp(r);
}

// chdtri returns the inverse of the complementary CDF
F c2icdf(F x, F y) {
    R chdtri(x, 1-y);
}

// fdtri returns the inverse of the complementary CDF
F ficdf(I x, I y, F z) {
    R fdtri(x, y, 1-z);
}

// gdtr interprets second parameter differently
F gcdf(F x, F y, F z) {
    R igam(x, z/y);
}

// no gdtri function
F gicdf(F x, F y, F z) {
    R y*igami(x, 1-z);
}

WRAPf(gamma)
WRAPff(beta)
WRAPnff(pgammar,igam)WRAPnff(pgammarc,igamc)WRAPnff(ipgammarc,igami)
WRAPnfff(pbetar,incbet)WRAPnfff(ipbetar,incbi)
WRAPnf(ncdf,ndtr)WRAPnf(nicdf,ndtri)
WRAPnif(stcdf,stdtr)WRAPnif(sticdf,stdtri)
WRAPniif(fcdf,fdtr)WRAPniif(ficdf,ficdf)
WRAPnff(c2cdf,chdtr)WRAPnff(c2icdf,c2icdf)
WRAPnfff(gcdf,gcdf)WRAPnfff(gicdf,gicdf)

// Utility functions
Z integer min_i(I x, I y) {
    R x<=y ? x : y;
}

Z integer max_i(I x, I y) {
    R x>=y ? x : y;
}

ZV swap_i(I* x, I* y) {
    I v = *x; *x = *y; *y = v;
}

// Make new F vector out of square matrix argument, column-major order
ZK mfms(K x, integer *n, char **err) {
    integer j;
    K f, r;
    if (xt!=0 || (*n=xn)==0)
        { *err="type"; R 0; }
    r = ktn(KF, *n**n);
    for (j=0; j<*n; ++j) {
        if ((f=mf(xK[j]))==0)
            { r0(r); *err="type"; R 0; }
        if (f->t<0 || f->n!=*n)
            { r0(f); r0(r); *err="length"; R 0; }
        DO(*n, kF(r)[*n*i+j] = kF(f)[i])
        r0(f);
    }
    R r;
}

// Make new F vector out of matrix argument, column-major order
ZK mfm(K x, integer *m, integer *n, char **err) {
    integer j;
    K f, r;
    if (xt!=0 || (*m=xn)==0 || (f=mf(xK[0]))==0)
        { *err="type"; R 0; }
    if (f->t<0 || (*n=f->n)==0)
        { r0(f); *err="length"; R 0; }
    r = ktn(KF, *m**n);
    for (j=0; j<*m; ++j) {
        if (j!=0) {
            if ((f=mf(xK[j]))==0)
                { r0(r); *err="type"; R 0; }
            if (f->t<0 || f->n!=*n)
                { r0(f); r0(r); *err="length"; R 0; }
        }
        DO(*n, kF(r)[*m*i+j] = kF(f)[i])
        r0(f);
    }
    R r;
}

// Represent a complex value
ZK mcv(F a, F b) {
    K x;
    if (b!=0)
        { x = ktn(KF, 2); xF[0] = a; xF[1] = b; R x; }
    else
        R kf(a);
}

// Matrix determinant
K QML_EXPORT qml_mdet(K x) {
    char* err;
    integer j, n;
    F r;
    K ipiv;
    x = mfms(x, &n, &err);
    P(x==0, krr(ss(err)))

    ipiv = ktn(QML_KLONG, n);
    dgetrf_(&n, &n, kF(x), &n, QML_kLONG(ipiv), &j);
    if (j!=0)
        { r0(ipiv); r0(x); R j>0 ? kf(0) : krr(ss("qml_assert")); }

    j = 0;
    DO(n, j += QML_kLONG(ipiv)[i]!=i+1);
    r0(ipiv);

    r = j%2 ? -1 : 1;
    DO(n, r *= kF(x)[i*n+i])
    r0(x);
    R kf(r);
}

// Matrix inverse
K QML_EXPORT qml_minv(K x) {
    char* err;
    integer j, n, lwork;
    F maxwork;
    K ipiv, work, a = mfms(x, &n, &err);
    P(a==0, krr(ss(err)))

    lwork = -1;
    dgetri_(&n, NULL, &n, NULL, &maxwork, &lwork, &j);
    if (j!=0)
        { r0(a); R krr(ss("qml_assert")); }
    lwork = max_i((integer)(maxwork+.5), n);

    ipiv = ktn(QML_KLONG, n);
    dgetrf_(&n, &n, kF(a), &n, QML_kLONG(ipiv), &j);
    if (j!=0)
        { r0(ipiv); r0(a); R j>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

    work = ktn(KF, lwork);
    dgetri_(&n, kF(a), &n, QML_kLONG(ipiv), kF(work), &lwork, &j);
    r0(work); r0(ipiv);
    if (j!=0)
        { r0(a); R j>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

    x = ktn(0, n);
    for (j=0; j<n; ++j) {
        xK[j] = ktn(KF, n);
        DO(n, kF(xK[j])[i] = kF(a)[i*n+j])
    }
    r0(a); R x;
}

// Eigenvalues and eigenvectors
K QML_EXPORT qml_mevu(K x) {
    char* err;
    integer j, n, lwork;
    F maxwork;
    K f, a = mfms(x, &n, &err);
    P(a==0, krr(ss(err)))

    lwork = -1;
    dgeev_("N", "V", &n, NULL, &n, NULL, NULL, NULL, &n,
        NULL, &n, &maxwork, &lwork, &j);
    if (j!=0)
        { r0(a); R krr(ss("qml_assert")); }
    lwork = max_i((integer)(maxwork+.5), 4*n);

    f = ktn(KF, n*2 + n*n + lwork);
    dgeev_("N", "V", &n, kF(a), &n, kF(f), kF(f)+n, NULL, &n,
        kF(f)+n*2, &n, kF(f)+n*2+n*n, &lwork, &j);
    r0(a);
    if (j!=0)
        { r0(f); R j>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

    x = ktn(0, 2);
    xK[1] = ktn(0, n);
    for (j=0; j<n; ++j)
        if (kF(f)[n+j]!=0)
            goto complex_result;

/*real_result:*/
    xK[0] = ktn(KF, n);
    for (j=0; j<n; ++j)
        kF(xK[0])[j] = kF(f)[j];
    for (j=0; j<n; ++j) {
        kK(xK[1])[j] = ktn(KF, n);
        DO(n, kF(kK(xK[1])[j])[i] = kF(f)[n*2+j*n+i])
    }
    goto done;

complex_result:
    xK[0] = ktn(0, n);
    for (j=0; j<n; ++j)
        kK(xK[0])[j] = mcv(kF(f)[j], kF(f)[n+j]);
    for (j=0; j<n; ++j) {
        if (kF(f)[n+j]!=0 && j+1<n) {
            kK(xK[1])[j] = ktn(0, n);
            DO(n, kK(kK(xK[1])[j])[i] =
                mcv(kF(f)[n*2+j*n+i], kF(f)[n*2+(j+1)*n+i]))
            kK(xK[1])[j+1] = ktn(0, n);
            DO(n, kK(kK(xK[1])[j+1])[i] =
                mcv(kF(f)[n*2+j*n+i], -kF(f)[n*2+(j+1)*n+i]))
            ++j;
        } else {
            kK(xK[1])[j] = ktn(KF, n);
            DO(n, kF(kK(xK[1])[j])[i] = kF(f)[n*2+j*n+i])
        }
    }

done:
    r0(f); R x;
}

// Cholesky decomposition
K QML_EXPORT qml_mchol(K x) {
    char* err;
    integer j, n;
    K a = mfms(x, &n, &err);
    P(a==0, krr(ss(err)))

    dpotrf_("U", &n, kF(a), &n, &j);
    if (j!=0)
        { r0(a); R k>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

    x = ktn(0, n);
    for (j=0; j<n; ++j) {
        xK[j] = ktn(KF, n);
        /*DO(j+1, kF(xK[j])[i] = kF(a)[i*n+j])
        DO(n-(j+1), kF(xK[j])[j+1+i] = 0)*/
        DO(j, kF(xK[j])[i] = 0)
        DO(n-j, kF(xK[j])[j+i] = kF(a)[(j+i)*n+j])
    }
    r0(a); R x;
}

// QR factorization
K QML_EXPORT qml_mqr(K x) {
    char* err;
    integer j, k, n, m, min, lwork;
    F maxwork, *q;
    K f, a = mfm(x, &m, &n, &err);
    P(a==0, krr(ss(err)))

    min = min_i(m, n);
    lwork = -1;
    dgeqrf_(&m, &n, NULL, &m, NULL, &maxwork, &lwork, &j);
    if (j!=0)
        { r0(a); R krr(ss("qml_assert")); }
    j = (integer)(maxwork+.5);
    dorgqr_(&m, &m, &min, NULL, &m, NULL, &maxwork, &lwork, &k);
    if (k!=0)
        { r0(a); R krr(ss("qml_assert")); }
    k = (integer)(maxwork+.5);
    lwork = max_i(max_i(j, k), n);

    f = ktn(KF, min + lwork + (n<m ? m*m : 0));
    dgeqrf_(&m, &n, kF(a), &m, kF(f), kF(f)+min, &lwork, &j);
    if (j!=0)
        { r0(f); r0(a); R krr(ss("qml_assert")); }

    x = ktn(0, 2);
    xK[0] = ktn(0, m);
    xK[1] = ktn(0, m);
    for (j=0; j<m; ++j) {
        kK(xK[1])[j] = ktn(KF, n);
        if (j<min) {
            DO(j, kF(kK(xK[1])[j])[i] = 0)
            DO(n-j, kF(kK(xK[1])[j])[j+i] = kF(a)[(j+i)*m+j])
        } else
            DO(n, kF(kK(xK[1])[j])[i] = 0)
    }

    if (n<m) {
        q = kF(f)+min+lwork;
        for (j=0; j<m; ++j)
            DO(n, q[i*m+j] = kF(a)[i*m+j])
    } else
        q = kF(a);
    dorgqr_(&m, &m, &min, q, &m, kF(f), kF(f)+min, &lwork, &j);
    if (j!=0)
        { r0(f); r0(a); r0(x); R krr(ss("qml_assert")); }

    for (j=0; j<m; ++j) {
        kK(xK[0])[j] = ktn(KF, m);
        DO(m, kF(kK(xK[0])[j])[i] = q[i*m+j])
    }
    r0(f); r0(a); R x;
}

// LUP factorization
K QML_EXPORT qml_mlup(K x) {
    char* err;
    integer j, m, n, min;
    K ipiv, a = mfm(x, &m, &n, &err);
    P(a==0, krr(ss(err)))

    min = min_i(m, n);
    ipiv = ktn(QML_KLONG, min);
    dgetrf_(&m, &n, kF(a), &m, QML_kLONG(ipiv), &j);
    if (j<0)
        { r0(a); R krr(ss("qml_assert")); }

    x = ktn(0, 3);
    xK[0] = ktn(0, m);
    for (j=0; j<m; ++j) {
        kK(xK[0])[j] = ktn(KF, min);
        if (j<min) {
            DO(j, kF(kK(xK[0])[j])[i] = kF(a)[i*m+j])
            kF(kK(xK[0])[j])[j] = 1;
            DO(min-(j+1), kF(kK(xK[0])[j])[j+1+i] = 0)
        } else
            DO(min, kF(kK(xK[0])[j])[i] = kF(a)[i*m+j])
    }
    xK[1] = ktn(0, min);
    for (j=0; j<min; ++j) {
        kK(xK[1])[j] = ktn(KF, n);
        DO(j, kF(kK(xK[1])[j])[i] = 0)
        DO(n-j, kF(kK(xK[1])[j])[j+i] = kF(a)[(j+i)*m+j])
    }
    xK[2] = ktn(KI, m);
    DO(m, kI(xK[2])[i] = i)
    DO(min, swap_i(kI(xK[2])+i, kI(xK[2])+(QML_kLONG(ipiv)[i]-1)))
    r0(ipiv); r0(a); R x;
}

// Singular value decomposition
K QML_EXPORT qml_msvd(K x) {
    char* err;
    integer j, m, n, min, lwork;
    F maxwork;
    K f, a = mfm(x, &m, &n, &err);
    P(a==0, krr(ss(err)))

    min = min_i(m, n);
    lwork = -1;
    dgesvd_("A", "A", &m, &n, NULL, &m, NULL, NULL, &m,
        NULL, &n, &maxwork, &lwork, &j);
    if (j!=0)
        { r0(a); R krr(ss("qml_assert")); }
    lwork = max_i((integer)(maxwork+.5), 3*min+m+n);

    f = ktn(KF, min + m*m + n*n + lwork);
    dgesvd_("A", "A", &m, &n, kF(a), &m, kF(f), kF(f)+min, &m,
        kF(f)+min+m*m, &n, kF(f)+min+m*m+n*n, &lwork, &j);
    r0(a);
    if (j!=0)
        { r0(f); R j>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

    x = ktn(0, 3);
    xK[0] = ktn(0, m);
    for (j=0; j<m; ++j) {
        kK(xK[0])[j] = ktn(KF, m);
        DO(m, kF(kK(xK[0])[j])[i] = kF(f)[min+i*m+j])
    }
    xK[1] = ktn(0, m);
    for (j=0; j<m; ++j) {
        kK(xK[1])[j] = ktn(KF, n);
        DO(n, kF(kK(xK[1])[j])[i] = 0)
        if (j<min)
            kF(kK(xK[1])[j])[j] = kF(f)[j];
    }
    xK[2] = ktn(0, n);
    for (j=0; j<n; ++j) {
        kK(xK[2])[j] = ktn(KF, n);
        DO(n, kF(kK(xK[2])[j])[i] = kF(f)[min+m*m+j*n+i])
    }
    r0(f); R x;
}


// Polynomial root finding
#ifdef QML_R_POLY
    void R_cpolyroot(double *opr, double *opi, int *degree, double *zeror,
        double *zeroi, int *fail, double* work);
#elif !defined(QML_LAPACK_POLY)
    int cpoly_(double *opr, double* opi, long* degree,
        double *zeror, double* zeroi, long* fail);
    int rpoly_(double *op, long* degree,
        double *zeror, double* zeroi, long* fail);
#endif

K QML_EXPORT qml_poly(K x) {
    integer j, n;
#ifndef QML_R_POLY
    I complex;
#endif
    K a, c, r;
#ifdef QML_LAPACK_POLY
    F d, maxwork[2];
    integer lwork;
    K work;
#elif defined(QML_R_POLY)
    K work;
#else
    integer k;
#endif

    if ((a=mf(x))==0) {
        P(xt!=0, krr(ss("type")))
        P((n=xn)==0, krr(ss("length")))
        c = ktn(KF, n*2);
        for (j=0; j<n; ++j) {
            if ((a=mf(xK[j]))==0)
                { r0(c); R krr(ss("type")); }
            if (a->t>0) {
                if (a->n!=1 && a->n!=2)
                    { r0(a); r0(c); R krr(ss("length")); }
                kF(c)[j] = kF(a)[0];
                kF(c)[n+j] = kF(a)[1];
            } else {
                kF(c)[j] = a->f;
                kF(c)[n+j] = 0;
            }
            r0(a);
        }
#ifndef QML_R_POLY
        complex = 1;
#endif
    } else {
        if (a->t<0)
            { r0(a); R krr(ss("type")); }
        if (a->n==0)
            { r0(a); R krr(ss("length")); }
        c = ktn(KF, (n=a->n)*2);
        DO(n, kF(c)[i] = kF(a)[i])
        DO(n, kF(c)[n+i] = 0)
        r0(a);
#ifndef QML_R_POLY
        complex = 0;
#endif
    }
#ifdef QML_R_POLY
    if (n > 50) /* matches constant in cpoly.c */
#else
    if (n > (complex ? 50 : 100)) /* matches constants in cpoly.c and rpoly.c */
        /* not strictly necessary for LAPACK version, but reasonable */
#endif
        { r0(c); R krr(ss("limit")); }
    --n;

    a = ktn(KF, n*2);
#ifdef QML_LAPACK_POLY
    /* compute roots as the eigenvalues of the companion matrix */
    if (complex)
        d = kF(c)[0]*kF(c)[0] + kF(c)[n+1]*kF(c)[n+1];
    else
        d = kF(c)[0];
    if (n==0)
        j = kF(c)[0]==0 && (!complex || kF(c)[n+1]==0);
    else if (d == 0)
        j = 1;
    else {
        lwork = -1;
        if (complex)
            zgeev_("N", "N", &n, NULL, &n, NULL, NULL, &n,
                NULL, &n, (doublecomplex*)maxwork, &lwork, NULL, &j);
        else
            dgeev_("N", "N", &n, NULL, &n, NULL, NULL, NULL, &n,
                NULL, &n, maxwork, &lwork, &j);
        if (j!=0)
            { r0(c); r0(a); R krr(ss("qml_assert")); }
        lwork = max_i((integer)(maxwork[0]+.5), (complex?2:3)*n);

        if (complex) {
            work = ktn(KF, 2*n*n + 2*lwork + 2*n);
            DO(2*n*n, kF(work)[i] = 0)
            for (j=0; j<n; ++j) {
                kF(work)[2*j*n] =
                    -(kF(c)[n+j+2]*kF(c)[n+1]+kF(c)[j+1]*kF(c)[0])/d;
                kF(work)[2*j*n+1] =
                    (kF(c)[j+1]*kF(c)[n+1]-kF(c)[n+j+2]*kF(c)[0])/d;
            }
            DO(n-1, kF(work)[2*(i*(n+1)+1)] = 1)
            zgeev_("N", "N", &n, (doublecomplex*)kF(work), &n,
                (doublecomplex*)kF(a), NULL, &n, NULL, &n,
                (doublecomplex*)(kF(work)+2*n*n), &lwork,
                kF(work)+2*n*n+2*lwork, &j);
        } else {
            work = ktn(KF, n*n + lwork);
            DO(n*n, kF(work)[i] = 0)
            DO(n, kF(work)[i*n] = -kF(c)[i+1]/d)
            DO(n-1, kF(work)[i*(n+1)+1] = 1)
            dgeev_("N", "N", &n, kF(work), &n, kF(a), kF(a)+n, NULL, &n,
                NULL, &n, kF(work)+n*n, &lwork, &j);
        }
        r0(work);
        if (j<0)
            { r0(c); r0(a); R krr(ss("roots")); }
    }
#elif defined(QML_R_POLY)
    work = ktn(KF, 10*(n+1)); /* matches expression in cpoly.c */
    R_cpolyroot(kF(c), kF(c)+n+1, &n, kF(a), kF(a)+n, &j, kF(work));
    r0(work);
#else
    if (!complex) {
        k = n;
        rpoly_(kF(c), &k, kF(a), kF(a)+n, &j);
        if (j)
            goto cpoly; /* fallback to cpoly, e.g. for 5 0 0 0 0 1 */
    } else
    cpoly:
        cpoly_(kF(c), kF(c)+n+1, &n, kF(a), kF(a)+n, &j);
#endif
    r0(c);
    if (j!=0)
        { r0(a); R krr(ss("roots")); }
    r = ktn(0, n);
    for (j=0; j<n; ++j)
#ifdef QML_LAPACK_POLY
        if (complex)
            kK(r)[j] = mcv(kF(a)[j*2], kF(a)[j*2+1]);
        else
#endif
            kK(r)[j] = mcv(kF(a)[j], kF(a)[n+j]);
    r0(a); R r;
}


// Initialization
static int initialized = 0;

#ifndef QML_DLLMAIN
__attribute__ ((constructor))
#endif
ZV initialize() {
    // these functions are only thread-unsafe on their first invocation
    slamch_("E");
    dlamch_("E");
    initialized = 1;
}

#ifdef QML_DLLMAIN
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved) {
    if (fdwReason == DLL_PROCESS_ATTACH)
        initialize();
    return TRUE;
}
#endif


// Constants
#define QUOTE1(x) (#x)
#define QUOTE(x) QUOTE1(x)

K QML_EXPORT qml_const(K x) {
    F r;
    P(!initialized, krr(ss("qml_assert"))) // catch possibly missed init here
    P(xt!=-KI, krr(ss("type")))
    SW(xi) {
        CS(0,R ks(QUOTE(QML_VERSION)))
        CS(1,r=3.1415926535897932384626433832795028842)
        CS(2,r=2.7182818284590452353602874713526624978)
        CS(3,r=DBL_EPSILON)
        CD: r=nf;
    }
    R kf(r);
}
