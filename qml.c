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
#include <conmax.h>

#include <k.h>

#undef QML_EXPORT
#ifdef QML_DLLEXPORT
    #define QML_EXPORT __declspec(dllexport)
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

// Return 1 if argument can be converted to I atom or vector
ZI mip(K x) {
    if (xt==0) {
        DO(xn, U(mi1p(xK[i])))
        R 1;
    } else
        R mi1p(x) || xt==KB || xt>=KG && xt<=KI;
}

// Make new I atom or vector out of B, G, H or I argument
Z K1(mi) {
    K n;
    U(mip(x))
    P(xt<0, mi1(x));
    n = ktn(KI, xn);
    SW(xt) {
        CS(KI, memcpy(kI(n), xI, xn*sizeof(I)))
        CS(KH, DO(xn, kI(n)[i] = xH[i]==nh ? ni :
            xH[i]==wh ? wi : xH[i]==-wh ? -wi : xH[i]))
        case KB:
        CS(KG, DO(xn, kI(n)[i] = xG[i]))
        CS(0, DO(xn, kI(n)[i] = mi1i(xK[i])))
    }
    R n;
}

// Return 1 if argument can be converted to F atom
ZI mf1p(K x) {
    R xt==-KB || xt>=-KF && xt<=-KG;
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

// Return 1 if argument can be converted to F atom or vector
ZI mfp(K x) {
    if (xt==0) {
        DO(xn, U(mf1p(xK[i])))
        R 1;
    } else
        R mf1p(x) || xt==KB || xt>=KG && xt<=KF;
}

// Make new F atom or vector out of B, G, H, I, J, E or F argument
Z K1(mf) {
    K f;
    U(mfp(x))
    P(xt<0, mf1(x));
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
        CS(0, DO(xn, kF(f)[i] = mf1f(xK[i])))
    }
    R f;
}

// Wrap a function of F
#define WRAPf(ff) WRAPnf(ff,ff)
#define WRAPnf(ffn,ff) \
K QML_EXPORT qml_##ffn(K u) { \
    K x = mf(u); \
    P(x==0, krr(ss("type"))) \
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
double gamma(double);
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

// this is not available in all distributions of Cephes so we provide it here
F beta(F x, F y) {
    I sx, sy, sxy;
    F r;
    r = lgamma_r(x, &sx);
    r += lgamma_r(y, &sy);
    r -= lgamma_r(x + y, &sxy);
    R sx*sy*sxy*exp(r);
}

// chdtri() returns the inverse of the complementary CDF
F c2icdf(F x, F y) {
    R chdtri(x, 1-y);
}

// fdtri() returns the inverse of the complementary CDF
F ficdf(I x, I y, F z) {
    R fdtri(x, y, 1-z);
}

// gdtr() interprets second parameter differently
F gcdf(F x, F y, F z) {
    R igam(x, z/y);
}

// no gdtri() function
F gicdf(F x, F y, F z) {
    R y*igami(x, 1-z);
}

// smirnov() returns complementary CDF
F smcdf(I x, F y) {
    R 1-smirnov(x, y);
}

// smirnovi() returns inverse of the complementary CDF
F smicdf(I x, F y) {
    R smirnovi(x, 1-y);
}

// komogorov() returns complementary CDF
F kcdf(F x) {
    R 1-kolmogorov(x);
}

// komogoi() returns inverse of the complementary CDF
F kicdf(F x) {
    P(x<1e-8,nf) // doesn't converge well for small values
    R kolmogi(1-x);
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
WRAPniif(bncdf,bdtr)WRAPniif(bnicdf,bdtri)
WRAPnif(pscdf,pdtr)WRAPnif(psicdf,pdtri)
WRAPnif(smcdf,smcdf)WRAPnif(smicdf,smicdf)
WRAPf(kcdf)WRAPf(kicdf)

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
        { r0(a); R j>0 ? ktn(0, 0) : krr(ss("qml_assert")); }

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

// QR factorization base
ZK qr(K x, const int pivot) {
    char* err;
    integer j, k, n, m, min, lwork;
    F maxwork, *q;
    K jpiv, f, a = mfm(x, &m, &n, &err);
    P(a==0, krr(ss(err)))

    min = min_i(m, n);
    lwork = -1;
    if (pivot)
        dgeqp3_(&m, &n, NULL, &m, NULL, NULL, &maxwork, &lwork, &j);
    else
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
    if (pivot) {
        jpiv = ktn(QML_KLONG, n);
        DO(n, QML_kLONG(jpiv)[i] = 0)
        dgeqp3_(&m, &n, kF(a), &m, QML_kLONG(jpiv), kF(f), kF(f)+min, &lwork,
            &j);
    } else
        dgeqrf_(&m, &n, kF(a), &m, kF(f), kF(f)+min, &lwork, &j);
    if (j!=0)
        { if (pivot) r0(jpiv); r0(f); r0(a); R krr(ss("qml_assert")); }

    x = ktn(0, pivot ? 3 : 2);
    if (pivot) {
        xK[2] = ktn(KI, n);
        DO(n, kI(xK[2])[i] = QML_kLONG(jpiv)[i]-1)
        r0(jpiv);
    }
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

// QR factorization
K QML_EXPORT qml_mqr(K x) {
    R qr(x, 0);
}

// QR factorization with column pivoting
K QML_EXPORT qml_mqrp(K x) {
    R qr(x, 1);
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
    int cpoly_(double *opr, double* opi, integer* degree,
        double *zeror, double* zeroi, logical* fail);
    int rpoly_(double *op, integer* degree,
        double *zeror, double* zeroi, logical* fail);
#endif

K QML_EXPORT qml_poly(K x) {
#ifdef QML_R_POLY
    I j, n;
#else
    integer j, n;
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
    else if (d==0)
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


// Options processing
struct opt {
    S s;
    I t;
};

union optv {
    I i; /* first member */
    F f;
};

ZI mopt_find(const S s, const struct opt* opt) {
    I k;
    for (k=0; opt[k].s!=0; k++)
        P(!strcmp(s, opt[k].s), k)
    R k;
}

// Fills in an optv array; returns 1 on success
ZI mopt(K x, const struct opt* opt, union optv* v) {
    I j, k;
    switch (xt) {
    case 0:
        for (j=0; j<xn; j++) {
            U(xK[j]->t==-KS)
            if (*xK[j]->s) {
                k = mopt_find(xK[j]->s, opt); U(opt[k].s!=0)
                SW(opt[k].t) {
                    CS(0, v[k].i=1)
                    CS(-KI, ++j; U(j<xn && mi1p(xK[j])) v[k].i=mi1i(xK[j]))
                    CS(-KF, ++j; U(j<xn && mf1p(xK[j])) v[k].f=mf1f(xK[j]))
                }
            }
        }
        break;
    case -KS:
        if (*xs) {
            k = mopt_find(xs, opt); U(opt[k].s!=0 && opt[k].t==0)
            v[k].i = 1;
        }
        break;
    case KS:
        for (j=0; j<xn; j++)
            if (*xS[j]) {
                k = mopt_find(xS[j], opt); U(opt[k].s!=0 && opt[k].t==0)
                v[k].i = 1;
            }
        break;
    case XD:
        U(xx->t==KS || xx->t==0)
        for (j=0; j<xx->n; j++) {
            S s;
            U(xx->t==KS || kK(xx)[j]->t==-KS)
            s = xx->t==KS ? kS(xx)[j] : kK(xx)[j]->s;
            if (*s) {
                k = mopt_find(s, opt);
                U(opt[k].s!=0)
                if (xy->t==0) {
                    U(j<xy->n)
                    SW(opt[k].t) {
                        case 0:
                        CS(-KI, U(mi1p(kK(xy)[j])) v[k].i=mi1i(kK(xy)[j]))
                        CS(-KF, U(mf1p(kK(xy)[j])) v[k].f=mf1f(kK(xy)[j]))
                    }
                } else {
                    K r;
                    SW(opt[k].t) {
                        case 0:
                        CS(-KI, U(mip(xy) && xy->t>=0 && j<xy->n)
                                r=mi(xy); v[k].i=kI(r)[j]; r0(r))
                        CS(-KF, U(mfp(xy) && xy->t>=0 && j<xy->n)
                                r=mf(xy); v[k].f=kF(r)[j]; r0(r))
                    }
                }
            }
        }
    }
    R 1;
}


// Nonlinear optimization
doublereal d1mach_3, d1mach_3_sqrt;

doublereal d1mach_(integer *i) {
    R d1mach_3;
}

struct fnset_info {
    K fun, con, start;
    I contyp;
    I neg, sig;
    K sigval;
    doublereal base;
};

Z doublereal* fnset_arg(K x, doublereal* param, K* a) {
    K y;
    if (xt==0) {
        y = ktn(0, xn);
        DO(xn, param=fnset_arg(xK[i], param, &kK(y)[i]))
    } else if (xt>0) {
        y = ktn(KF, xn);
        DO(xn, kF(y)[i]=param[i])
        param += xn;
    } else
        y = kf(*param++);
    *a = y;
    R param;
}

ZF fnset_call(struct fnset_info* info, K f, I neg, doublereal* param) {
    K a, x;
    F v;
    P(info->sig, 0)

    if (info->start==0 && info->con!=0) /* line */
        a = knk(1, kf(info->neg ? info->base-*param : info->base+*param));
    else if (info->start==0 || info->start->t<0) /* root or scalar param */
        a = knk(1, kf(*param));
    else
        fnset_arg(info->start, param, &a);

    x = dot(f, a); r0(a);
    if (x!=0 && mf1p(x))
        { v = mf1f(x); r0(x); R neg ? -v : v; }

    /* function didn't return a float as we'd hoped */
    info->sig = 1;
    if (x==0 || x->t==-128)
        info->sigval = x;
    else
        { S s = x->t>=100 ? "rank" : "type"; r0(x); info->sigval = krr(ss(s)); }
    R 0;
}

int fnset_(integer *nparm, integer *numgr, void *pttbl, doublereal *param,
           integer *ipt, integer *indfn, integer *icntyp, doublereal *confun)
{
    I j, neg;
    F v;
    K f;
    struct fnset_info* info = pttbl;

    icntyp += *ipt - 1;
    confun += *ipt - 1;

    if (info->con==0) { /* root or solve */
        f = info->fun->t==0 ? kK(info->fun)[*ipt-1] : info->fun;
        neg = info->neg;
        *icntyp = 2;
    } else /* line, min or conmin */
        if (*ipt==1) { /* objective function */
            f = info->fun;
            neg = 0; *icntyp = 1;
        } else { /* constraints */
            f = info->con->t==0 ? kK(info->con)[*ipt-2] : info->con;
            neg = 1; *icntyp = info->contyp;
        }

    *confun = v = fnset_call(info, f, neg, param);
    if (*indfn) {
        I m = *numgr, n = *nparm;
        if (*icntyp==-1) /* linear function */
            for (j=0; j<n; j++) {
                F p = param[j];
                param[j] = p + 1;
                *(confun+=m) = fnset_call(info, f, neg, param) - v;
                param[j] = p;
            }
        else /* nonlinear function */
            for (j=0; j<n; j++) {
                F p = param[j], h = d1mach_3_sqrt, p1, p2, v1, v2;
                param[j] = p1 = p + h; v1 = fnset_call(info, f, neg, param);
                param[j] = p2 = p - h; v2 = fnset_call(info, f, neg, param);
                *(confun+=m) = (v1-v2) / (p1-p2);
                param[j] = p;
            }
    }

    if (info->sig)
        *icntyp = 0;
    R 0;
}

ZI solvemin_count(K x) {
    if (xt==0) {
        I j, c, n = 0;
        for (j=0; j<xn; j++) {
            c = solvemin_count(xK[j]); P(c<0, c)
            n += c; P(n>1000000, -2)
        }
        R n;
    } else {
        P(!mfp(x), -1)
        R xt<0 ? 1 : xn>1000000 ? -2 : xn;
    }
}

Z doublereal* solvemin_param(K x, doublereal* param) {
    if (xt==0)
        DO(xn, param=solvemin_param(xK[i], param))
    else if (xt>0) {
        x = mf(x); /* type was checked by solvemin_count */
        DO(xn, param[i]=xF[i])
        param += xn;
        r0(x);
    } else
        *param++ = mf1f(x);
    R param;
}

Z struct k0 con_dummy; /* empty con used only by min */

ZK solvemin(K fun, K con, K start, I maxiter, F tolcon, I steps, I slp, I rk,
            I lincon, I full, I quiet)
{
    integer ioptn, nparm, numgr, itlim, ifun, liwrk, lwrk, iter;
    doublereal *param, *error;
    K iwork, work, x, y;
    S sig;
    struct fnset_info info;

    nparm = solvemin_count(start);
    P(nparm<0, krr(ss(nparm==-2 ? "limit" : "type")))
    P(!(fun->t>=100 || con==0 && fun->t==0), krr(ss("type")))
    P(!(con==0 || con->t>=100 || con->t==0), krr(ss("type")))
    P(fun->t==0 && fun->n>1000000, krr(ss("limit")))
    P(con!=0 && con->t==0 && con->n>1000000, krr(ss("limit")))
    P(slp && (rk || steps>0), krr(ss("opt")))

    if (con==0) { /* solve */
        numgr = fun->t==0 ? fun->n : 1;
        ifun = numgr;
    } else { /* min or conmin */
        numgr = 1 + (con->t==0 ? con->n : 1);
        ifun = 1;
    }

    P(2*nparm+4*numgr>2000000000/nparm, krr(ss("limit")))
    liwrk = 7*numgr + 7*nparm + 3;
    lwrk = 2*nparm*nparm + 4*numgr*nparm + 11*numgr + 27*nparm + 13;
    iwork = ktn(QML_KLONG, liwrk);
    work = ktn(KF, lwrk + ifun + (numgr+3) + nparm);

    error = kF(work) + lwrk + ifun;
    param = error + (numgr+3);
    solvemin_param(start, param);

    ioptn = 200;
    itlim = maxiter>=0 ? maxiter : 0;
    if (!(tolcon>=0))
        tolcon = d1mach_3_sqrt;
    kF(work)[1] = tolcon;
    if (steps>0) {
        QML_kLONG(iwork)[1] = steps;
        ioptn += 100;
    }
    ioptn += slp ? 1000 : rk ? 2000 : 0;
    DO(ifun, (kF(work)+lwrk)[i]=0)

    info.fun = fun;
    info.con = con;
    info.start = start;
    info.contyp = lincon ? -1 : -2;
    info.neg = 0;
    info.sig = 0;

    conmax_(&ioptn, &nparm, &numgr, &itlim, kF(work)+lwrk, &ifun, &info,
            QML_kLONG(iwork), &liwrk, kF(work), &lwrk, &iter,
            param, kF(work)+lwrk+ifun);
    r0(iwork);
    if (info.sig)
        { r0(work); R info.sigval; }

    if (iter>=maxiter)
        sig = "iter";
    else if (iter<0 || con==0 && !(param[-3]>=-tolcon && param[-3]<=tolcon))
        sig = "feas";
    else {
        I j;
        sig = NULL;
        for (j=0; j<nparm; j++)
            if (isnan(param[j]))
                { sig="nan"; break; }
    }
    if (sig && !quiet)
        { r0(work); R krr(ss(sig)); }

    if (sig) {
        if (full)
            fnset_arg(start, param, &y);
        DO(nparm, param[i]=nf)
    }
    fnset_arg(start, param, &x);
    if (full) {
        /*          solve       min        conmin */
        /* normal: `x`iter      `x`f`iter  `x`f`cons`iter */
        /* error:  `x`last`err  `x`last`f  `x`last`f`cons */
        /*         `iter`sig    `iter`sig  `err`iter`sig */
        I j = 2 + (sig ? 2 + (con!=&con_dummy) : 0) +
               (con!=0 ? 1 + (con!=&con_dummy) : 0);
        K k = ktn(KS, j), v = ktn(0, j);
        if (sig) kS(k)[--j] = ss("sig"),  kK(v)[j] = ks(ss(sig));
        /*    */ kS(k)[--j] = ss("iter"), kK(v)[j] = ki(iter>=0 ? (I)iter : 0);
        if (sig && con!=&con_dummy) {
            F err = param[con==0 ? -3 : lincon ? -2 : -1];
            /**/ kS(k)[--j] = ss("err"),  kK(v)[j] = kf(err);
        }
        if (con!=0) {
            if (con!=&con_dummy) {
                 F e; K cons = ktn(KF, numgr-1);
                 DO(numgr-1, e = error[i+1];
                             kF(cons)[i] = e>=-tolcon && e<=tolcon ? 0 : -e)
                 kS(k)[--j] = ss("cons"), kK(v)[j] = cons;
            }
            /**/ kS(k)[--j] = ss("f"),    kK(v)[j] = kf(param[-3]);
        }
        if (sig) kS(k)[--j] = ss("last"), kK(v)[j] = y;
        /*    */ kS(k)[--j] = ss("x"),    kK(v)[j] = x;
        x = xD(k, v);
    }
    r0(work);
    R x;
}

ZK root(K fun, K x, I maxiter, doublereal tolcon, I full, I quiet) {
    S sig;
    I sigsign;
    integer iter, nsrch, ioptn, nparm, numgr, ifun, iphse, iwork[17],
            liwrk, lwrk;
    doublereal dvec, cfun, zwork, work[6], parwrk, err1[4], p1, p2, f1, f2;
    struct fnset_info info;

    P(fun->t<100, krr(ss("type")))
    x = mf(x); P(x==0, krr(ss("type")))
    if (xt<0 || xn!=2)
        { r0(x); R krr(ss("length")); }
    p1 = xF[0];
    p2 = xF[1];
    r0(x);

    info.fun = fun;
    info.con = NULL; /* root/solve flag */
    info.start = NULL; /* root/line flag */
    info.contyp = -2;
    info.neg = 0;
    info.sig = 0;

    if (p1>p2)
        { f2 = p1; p1 = p2; p2 = f2; }
    f1 = fnset_call(&info, fun, 0, &p1);
    f2 = fnset_call(&info, fun, 0, &p2);
    P(info.sig, info.sigval)

    if (!(tolcon>=0))
        tolcon = d1mach_3_sqrt;
    if (f1<-tolcon && f2>tolcon)
        { f1 = -f1; f2 = -f2; info.neg = 1; }
    sigsign = !(f1>tolcon && f2<-tolcon);

    iter = maxiter>=0 ? maxiter : 0;
    ifun = dvec = numgr = nparm = 1;
    zwork = cfun = iphse = ioptn = 0;
    liwrk = 17;
    lwrk = 6;
    iwork[15] = -2;

    muller_(&iter, &nsrch, &ioptn, &nparm, &numgr, &dvec, &cfun, &ifun,
            &info, &zwork, &tolcon, &iphse, iwork, &liwrk,
            work, &lwrk, &parwrk, err1, &p1, &f1, &p2, &f2);
    if (info.sig)
        R info.sigval;

    sig = nsrch>=maxiter ? "iter" :
          !(f2>=-tolcon && f2<=tolcon) ? "feas" :
          isnan(p2) ? "nan" : NULL;
    if (sig && sigsign)
        sig = "sign";
    P(sig && !quiet, krr(ss(sig)))

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
    R x;
}

ZK line(K fun, K base, K x, I maxiter, doublereal tolcon, I full, I quiet) {
    S sig;
    integer initlm, nadd, lims1, ioptn, numgr, nparm, ifun, mact, iact, iphse,
            itypm1, itypm2, iwork[17], liwrk, lwrk, nsrch;
    doublereal prjlim, tol1, cx[2], cfun, param, error[4], rchdwn, unit, rchin,
               work[42], err1[4], parprj, projct, emin, emin1, parser;
    struct fnset_info info;

    P(fun->t<100, krr(ss("type")))
    P(!mf1p(base) || !mf1p(x), krr(ss("type")))

    info.fun = fun;
    info.con = fun; /* min/line flag */
    info.start = NULL; /* root/line flag */
    info.neg = 0;
    info.sig = 0;
    info.base = mf1f(base);

    projct = mf1f(x) - info.base;
    if (!(projct>0))
        { projct = -projct; info.neg = 1; }

    if (!(tolcon>=0))
        tolcon = d1mach_3_sqrt;
    prjlim = wf;
    tol1 = 100 * d1mach_3;
    if (tol1>tolcon)
        tol1=tolcon;
    lims1 = nadd = initlm = maxiter>=0 ? maxiter : 0;
    unit = iact = mact = ifun = nparm = numgr = 1;
    param = cfun = itypm2 = itypm1 = iphse = ioptn = 0;
    rchin = rchdwn = 2;
    liwrk = 17;
    lwrk = 42;
    cx[0] = 1;
    iwork[6] = 1;

    searsl_(&initlm, &nadd, &lims1, &ioptn, &numgr, &nparm, &prjlim, &tol1, cx,
            &cfun, &ifun, &info, &param, error, &rchdwn, &mact, &iact, &iphse,
            &unit, &tolcon, &rchin, &itypm1, &itypm2, iwork, &liwrk, work,
            &lwrk, err1, &parprj, &projct, &emin, &emin1, &parser, &nsrch);
    if (info.sig)
        R info.sigval;

    sig = nsrch>=maxiter || nsrch>=lims1 ? "iter" :
          isnan(projct) ? "nan" : NULL;
    P(sig && !quiet, krr(ss(sig)));

    projct = info.neg ? info.base-projct : info.base+projct;
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
    R x;
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

K QML_EXPORT qml_solvex(K opts, K x, K y) {
    union optv v[] =
        { { 1000 }, { .f = -1 }, { -1 }, { 0 }, { 0 }, { 0 }, { 0 } };
    P(opts!=0 && !mopt(opts, solve_opt, v), krr(ss("opt")))
    R solvemin(x, NULL, y,
               v[0].i, v[1].f, v[2].i, v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K QML_EXPORT qml_solve(K x, K y) {
    R qml_solvex(NULL, x, y);
}

K QML_EXPORT qml_minx(K opts, K x, K y) {
    union optv v[] =
        { { 1000 }, { .f = -1 }, { -1 }, { 0 }, { 0 }, { 0 }, { 0 } };
    P(opts!=0 && !mopt(opts, solve_opt, v), krr(ss("opt")))
    R solvemin(x, &con_dummy, y,
               v[0].i, v[1].f, v[2].i, v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K QML_EXPORT qml_min(K x, K y) {
    R qml_minx(NULL, x, y);
}

Z const struct opt conmin_opt[] = {
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

K QML_EXPORT qml_conminx(K opts, K x, K y, K z) {
    union optv v[] =
        { { 1000 }, { .f = -1 }, { -1 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    P(opts!=0 && !mopt(opts, conmin_opt, v), krr(ss("opt")))
    R solvemin(x, y, z,
               v[0].i, v[1].f, v[2].i, v[3].i, v[4].i, v[5].i, v[6].i, v[7].i);
}

K QML_EXPORT qml_conmin(K x, K y, K z) {
    R qml_conminx(NULL, x, y, z);
}

Z const struct opt rootline_opt[] = {
    { "iter", -KI },
    { "tol",  -KF },
    { "full",   0 },
    { "quiet",  0 },
    { NULL }
};

K QML_EXPORT qml_rootx(K opts, K x, K y) {
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    P(opts!=0 && !mopt(opts, rootline_opt, v), krr(ss("opt")))
    R root(x, y, v[0].i, v[1].f, v[2].i, v[3].i );
}

K QML_EXPORT qml_root(K x, K y) {
    R qml_rootx(NULL, x, y);
}

K QML_EXPORT qml_linex(K opts, K x, K y, K z) {
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    P(opts!=0 && !mopt(opts, rootline_opt, v), krr(ss("opt")))
    R line(x, y, z, v[0].i, v[1].f, v[2].i, v[3].i);
}

K QML_EXPORT qml_line(K x, K y, K z) {
    R qml_linex(NULL, x, y, z);
}


// Initialization
static int initialized = 0;

#ifndef QML_DLLMAIN
__attribute__ ((constructor))
#endif
ZV initialize() {
    // these functions are only thread-unsafe on their first invocation
    slamch_("E");
    d1mach_3 = dlamch_("E");
    d1mach_3_sqrt = sqrt(d1mach_3);
    con_dummy.t = 0;
    con_dummy.n = 0; /* can't initialize with designators in initializer */
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
    P(!initialized, krr(ss("qml_assert"))) // catch potentially missed init here
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
