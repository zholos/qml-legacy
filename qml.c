#include <math.h>
#include <float.h>
#include <fdlibm.h>
#include <stdio.h>
#include <k.h>

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
K __declspec(dllexport) qml_##ffn(K u) { \
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
K __declspec(dllexport) qml_##ff(K u, K v) { \
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
K __declspec(dllexport) qml_##ffn(K u, K v) { \
    P(!mf1p(u) || !mf1p(v), krr(ss("type"))) \
    R kf(ff(mf1f(u), mf1f(v))); \
}

// Wrap a function of (I; F)
#define WRAPnif(ffn,ff) \
K __declspec(dllexport) qml_##ffn(K u, K v) { \
    P(!mi1p(u) || !mf1p(v), krr(ss("type"))) \
    R kf(ff(mi1i(u), mf1f(v))); \
}

// Wrap a function of (F; F; F)
#define WRAPnfff(ffn,ff) \
K __declspec(dllexport) qml_##ffn(K u, K v, K w) { \
    P(!mf1p(u) || !mf1p(v) || !mf1p(w), krr(ss("type"))) \
    R kf(ff(mf1f(u), mf1f(v), mf1f(w))); \
}

// Wrap a function of (I; I; F)
#define WRAPniif(ffn,ff) \
K __declspec(dllexport) qml_##ffn(K u, K v, K w) { \
    P(!mi1p(u) || !mi1p(v) || !mf1p(w), krr(ss("type"))) \
    R kf(ff(mi1i(u), mi1i(v), mf1f(w))); \
}

// fdlibmm functions
/*WRAP1(cos)WRAP1(sin)*/WRAPf(tan)WRAPf(acos)WRAPf(asin)WRAPf(atan)WRAPff(atan2)
WRAPf(cosh)WRAPf(sinh)WRAPf(tanh)WRAPf(acosh)WRAPf(asinh)WRAPf(atanh)
WRAPf(exp)WRAPf(log)WRAPf(log10)WRAPf(logb)WRAPf(expm1)WRAPf(log1p)WRAPff(pow)
WRAPf(floor)WRAPf(ceil)WRAPf(fabs)WRAPff(fmod)
WRAPf(erf)WRAPf(erfc)WRAPf(lgamma)
WRAPf(j0)WRAPf(j1)WRAPf(y0)WRAPf(y1)
/*WRAP1(sqrt)*/WRAPf(cbrt)WRAPff(hypot)

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
    I s;
    F r;
    r = lgamma(x); s = signgam;
    r += lgamma(y); s *= signgam;
    r -= lgamma(x + y); s *= signgam;
    R s*exp(r);
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

// Make new F vector out of square matrix argument, column-major order
ZK mfm(K x, I *n, const char **err) {
    I j;
    K f, r;
    if (xt!=0 || (*n=xn)==0)
        { *err="type"; R 0; }
    r = ktn(KF, *n**n);
    for (j=0; j<*n; ++j) {
        if ((f=mf(xK[j]))==0)
            { r0(r); *err="type"; R 0; }
        if (f->t<0 || f->n!=*n)
            { r0(f); r0(r); *err="length"; R 0; }
        DO(*n, kF(r)[*n*i+j]=kF(f)[i])
        r0(f);
    }
    R r;
}

// Represent a complex value
ZK mcv(F a, F b) {
    K x;
    if (b!=0) {
        x = ktn(KF, 2);
        xF[0] = a;
        xF[1] = b;
        return x;
    } else
        return kf(a);
}

// Matrix determinant
int dgetrf_(int *m, int *n, double *a, int *lda, int* ipiv, int *info);

K __declspec(dllexport) qml_mdet(K x) {
    char* err;
    I j, n;
    F r;
    K ipiv;
    x = mfm(x, &n, &err);
    P(x==0, krr(ss(err)))

    ipiv = ktn(KI, n);
    dgetrf_(&n, &n, kF(x), &n, kI(ipiv), &j);
    if (j!=0)
        { r0(ipiv); r0(x); R kf(j>0 ? 0 : nf); }

    j = 0;
    DO(n, j+=kI(ipiv)[i]!=i+1);
    r0(ipiv);

    r = j%2 ? -1 : 1;
    DO(n, r*=kF(x)[i*n+i])
    r0(x);
    R kf(r);
}

// Matrix inverse
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork,
    int *info);

K __declspec(dllexport) qml_minv(K x) {
    char* err;
    I j, n, lwork;
    K ipiv, work, a = mfm(x, &n, &err);
    P(a==0, krr(ss(err)))

    ipiv = ktn(KI, n);
    dgetrf_(&n, &n, kF(a), &n, kI(ipiv), &j);
    if (j!=0)
        { r0(ipiv); r0(a); R ktn(0, 0); }

    lwork = 64*n;
    work = ktn(KF, lwork);
    dgetri_(&n, kF(a), &n, kI(ipiv), kF(work), &lwork, &j);
    r0(work); r0(ipiv);
    if (j!=0)
        { r0(a); R ktn(0, 0); }

    x = ktn(0, n);
    for (j=0; j<n; ++j) {
        xK[j] = ktn(KF, n);
        DO(n, kF(xK[j])[i]=kF(a)[i*n+j])
    }
    r0(a); R x;
}

// Eigenvalues and eigenvectors
int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
    double *wr, double* wi_, double *vl, int *ldvl, double *vr, int *ldvr,
    double *work, int *lwork, int *info);

K __declspec(dllexport) qml_mevu(K x) {
    char* err;
    I j, n, lwork;
    K f, a = mfm(x, &n, &err);
    P(a==0, krr(ss(err)))

    lwork = 50*n;
    f = ktn(KF, n*2 + n*n + lwork);
    dgeev_("N", "V", &n, kF(a), &n, kF(f), kF(f)+n, NULL, &n,
        kF(f)+n*2, &n, kF(f)+n*2+n*n, &lwork, &j);
    r0(a);
    if (j!=0)
        { r0(f); R ktn(0, 0); }

    x = ktn(0, 2);
    xK[1] = ktn(0, n);
    for (j=0; j<n; ++j)
        if (kF(f)[n+j]!=0)
            goto complex;

real:
    xK[0] = ktn(KF, n);
    for (j=0; j<n; ++j)
        kF(xK[0])[j] = kF(f)[j];
    for (j=0; j<n; ++j) {
        kK(xK[1])[j] = ktn(KF, n);
        DO(n, kF(kK(xK[1])[j])[i]=kF(f)[n*2+j*n+i])
    }
    goto done;

complex:
    xK[0] = ktn(0, n);
    for (j=0; j<n; ++j)
        kK(xK[0])[j] = mcv(kF(f)[j], kF(f)[n+j]);
    for (j=0; j<n; ++j) {
        if (kF(f)[n+j]!=0 && j+1<n) {
            kK(xK[1])[j] = ktn(0, n);
            DO(n, kK(kK(xK[1])[j])[i]=
                mcv(kF(f)[n*2+j*n+i], kF(f)[n*2+(j+1)*n+i]))
            kK(xK[1])[j+1] = ktn(0, n);
            DO(n, kK(kK(xK[1])[j+1])[i]=
                mcv(kF(f)[n*2+j*n+i], -kF(f)[n*2+(j+1)*n+i]))
            ++j;
        } else {
            kK(xK[1])[j] = ktn(KF, n);
            DO(n, kF(kK(xK[1])[j])[i]=kF(f)[n*2+j*n+i])
        }
    }

done:
    r0(f); R x;
}

// Cholesky decomposition
int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);

K __declspec(dllexport) qml_mchol(K x) {
    char* err;
    I j, n;
    K f, a = mfm(x, &n, &err);
    P(a==0, krr(ss(err)))

    dpotrf_("U", &n, kF(a), &n, &j);
    if (j!=0)
        { r0(a); R ktn(0, 0); }

    x = ktn(0, n);
    for (j=0; j<n; ++j) {
        xK[j] = ktn(KF, n);
        /*DO(j+1, kF(xK[j])[i]=kF(a)[i*n+j])
        DO(n-(j+1), kF(xK[j])[j+1+i]=0)*/
        DO(j, kF(xK[j])[i]=0)
        DO(n-j, kF(xK[j])[j+i]=kF(a)[(j+i)*n+j])
    }
    r0(a); R x;
}

// Singular value decomposition
int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda,
    double *s, double *u, int *ldu, double *vt, int *ldvt, double *work,
    int *lwork, int *info);

K __declspec(dllexport) qml_msvd(K x) {
    I j, n, m, min, lwork;
    K a, f;
    P(xt!=0 || (n=xn)==0, krr(ss("type")))
    for (j=0; j<n; ++j) {
        if (j==0) {
            if ((f=mf(xK[j]))==0 || f->t<0 || (m=f->n)==0)
                { if(f!=0) r0(f); R krr(ss("type")); }
            a = ktn(KF, n*m);
        } else {
            if ((f=mf(xK[j]))==0)
                { r0(a); R krr(ss("type")); }
            if (f->t<0 || f->n!=m)
                { r0(f); r0(a); R krr(ss("length")); }
        }
        DO(m, kF(a)[n*i+j]=kF(f)[i])
        r0(f);
    }

    min = m < n ? m : n;
    lwork = 4 * min + m + n;
    f = ktn(KF, min + n*n + m*m + lwork);
    dgesvd_("A", "A", &n, &m, kF(a), &n, kF(f), kF(f)+min, &n,
        kF(f)+min+n*n, &m, kF(f)+min+n*n+m*m, &lwork, &j);
    r0(a);
    if (j!=0)
        { r0(f); R ktn(0, 0); }

    x = ktn(0, 3);
    xK[0] = ktn(0, n);
    for (j=0; j<n; ++j) {
        kK(xK[0])[j] = ktn(KF, n);
        DO(n, kF(kK(xK[0])[j])[i]=kF(f)[min+i*n+j])
    }
    xK[1] = ktn(0, n);
    for (j=0; j<n; ++j) {
        kK(xK[1])[j] = ktn(KF, m);
        DO(m, kF(kK(xK[1])[j])[i]=0)
        if (j<min)
            kF(kK(xK[1])[j])[j] = kF(f)[j];
    }
    xK[2] = ktn(0, m);
    for (j=0; j<m; ++j) {
        kK(xK[2])[j] = ktn(KF, m);
        DO(m, kF(kK(xK[2])[j])[i]=kF(f)[min+n*n+j*m+i])
    }
    r0(f); R x;
}


// Polynomial root finding
#ifdef USE_LAPACK_POLY
    int zgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
        double *w, double *vl, int *ldvl, double *vr, int *ldvr,
        double *work, int *lwork, double* rwork, int *info);
#elif defined(USE_R_POLY)
    void R_cpolyroot(double *opr, double *opi, int *degree, double *zeror,
        double *zeroi, int *fail, double* work);
#else
    int cpoly_(double *opr, double* opi, long* degree,
        double *zeror, double* zeroi, long* fail);
    int rpoly_(double *op, long* degree,
        double *zeror, double* zeroi, long* fail);
#endif

K __declspec(dllexport) qml_poly(K x) {
    I j, n;
#ifndef USE_R_POLY
    I complex;
#endif
    K a, c, r;
#ifdef USE_LAPACK_POLY
    F d;
    I lwork;
    K work;
#elif defined(USE_R_POLY)
    K work;
#else
    I k;
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
#ifndef USE_R_POLY
        complex = 1;
#endif
    } else {
        if (a->t<0)
            { r0(a); R krr(ss("type")); }
        if (a->n==0)
            { r0(a); R krr(ss("length")); }
        c = ktn(KF, (n=a->n)*2);
        DO(n, kF(c)[i]=kF(a)[i])
        DO(n, kF(c)[n+i]=0)
        r0(a);
#ifndef USE_R_POLY
        complex = 0;
#endif
    }
#ifdef USE_R_POLY
    if (n > 50) /* matches constant in cpoly.c */
#else
    if (n > (complex ? 50 : 100)) /* matches constants in cpoly.c and rpoly.c */
        /* not strictly necessary for LAPACK version, but reasonable */
#endif
        { r0(c); R krr(ss("limit")); }
    --n;

    a = ktn(KF, n*2);
#ifdef USE_LAPACK_POLY
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
        lwork = 50*n;
        if (complex) {
            work = ktn(KF, 2*(n*n + lwork) + 2*n + lwork);
            DO(2*n*n, kF(work)[i]=0)
            for (j=0; j<n; ++j) {
                kF(work)[2*j*n] =
                    -(kF(c)[n+j+2]*kF(c)[n+1]+kF(c)[j+1]*kF(c)[0])/d;
                kF(work)[2*j*n+1] =
                    (kF(c)[j+1]*kF(c)[n+1]-kF(c)[n+j+2]*kF(c)[0])/d;
            }
            DO(n-1, kF(work)[2*(i*(n+1)+1)]=1)
            zgeev_("N", "N", &n, kF(work), &n, kF(a), NULL, &n,
                NULL, &n, kF(work)+2*n*n, &lwork, kF(work)+2*n*n+2*lwork, &j);
        } else {
            work = ktn(KF, n*n + lwork);
            DO(complex ? 2*n*n : n*n, kF(work)[i]=0)
            DO(n, kF(work)[i*n]=-kF(c)[i+1]/d)
            DO(n-1, kF(work)[i*(n+1)+1]=1)
            dgeev_("N", "N", &n, kF(work), &n, kF(a), kF(a)+n, NULL, &n,
                NULL, &n, kF(work)+n*n, &lwork, &j);
        }
        r0(work);
    }
#elif defined(USE_R_POLY)
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
    if (j)
        { r0(a); R krr(ss("roots")); }
    r = ktn(0, n);
    for (j=0; j<n; ++j)
#ifdef USE_LAPACK_POLY
        if (complex)
            kK(r)[j] = mcv(kF(a)[j*2], kF(a)[j*2+1]);
        else
#endif
            kK(r)[j] = mcv(kF(a)[j], kF(a)[n+j]);
    r0(a); R r;
}


// Constants
#define QUOTE1(x) (#x)
#define QUOTE(x) QUOTE1(x)

K __declspec(dllexport) qml_const(K x) {
    F r;
    P(xt!=-KI, krr(ss("type")))
    SW(xi) {
        CS(0,R ks(QUOTE(VERSION)))
        CS(1,r=M_PI)
        CS(2,r=M_E)
        CS(3,r=DBL_EPSILON)
        CD: r=nf;
    }
    R kf(r);
}
