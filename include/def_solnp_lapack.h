#ifndef DEF_SOLNP_LAPACK_H_GUARD
#define DEF_SOLNP_LAPACK_H_GUARD
#ifdef __cplusplus
"C"
{
#endif

/* Name Manglng for lapack functions*/
#ifndef LAPACK_GLOBAL
#if defined(LAPACK_GLOBAL_PATTERN_LC) || defined(ADD_)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#elif defined(LAPACK_GLOBAL_PATTERN_UC) || defined(UPPER)
#define LAPACK_GLOBAL(lcname,UCNAME)  UCNAME
#elif defined(LAPACK_GLOBAL_PATTERN_MC) || defined(NOCHANGE)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname
#else
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname
#endif
#endif


/* Dense direct */
#define LAPACK_RET_OK (0)
#define LAPACK_UPLOW_LOW ('L')
#define LAPACK_UPLOW_UP ('U')
#define LAPACK_NOTRANS ('N')
#define LAPACK_TRANS ('T')
#define LAPACK_SIDE_LEFT ('L')
#define LAPACK_DIAG_NONUNIT ('N')

    void dpotrf(const char *uplo, const int *n, double *a, const int *lda, int *info);
    void dpotrs(const char *uplo, const int *n, const int *nrhs, const double *a,
                const int *lda, double *b, const int *ldb, int *info);
    void dpocon(const char *uplo, const int *n, const double *a, const int *lda,
                const double *anorm, double *rcond, double *work, int *iwork,
                int *info);
    void dgetrf(const int *m, const int *n, double *a, const int *lda, int *ipiv,
                int *info);
    void dgetrs(const char *trans, const int *n, const int *nrhs, const double *a,
                const int *lda, const int *ipiv, double *b, const int *ldb,
                int *info);
    void dsytrf(const char *uplo, const int *n, double *a, const int *lda,
                int *ipiv, double *work, const int *lwork, int *info);
    void dsytrs(const char *uplo, const int *n, const int *nrhs, const double *a,
                const int *lda, const int *ipiv, double *b, const int *ldb, int *info);
    void dgelss(const int *m, const int *n, const int *nrhs, double *a, const int *lda,
                double *b, const int *ldb, double *s, const double *rcond, int *rank,
                double *work, const int *lwork, int *info);

    // for name manging
    #define dpotrf_base LAPACK_GLOBAL(dpotrf, DPOTRF)
    #define dpotrs_base LAPACK_GLOBAL(dpotrs, DPOTRS)
    #define dpocon_base LAPACK_GLOBAL(dpocon, DPOCON)
    #define dgetrf_base LAPACK_GLOBAL(dgetrf, DGETRF)
    #define dgetrs_base LAPACK_GLOBAL(dgetrs, DGETRS)
    #define dsytrf_base LAPACK_GLOBAL(dsytrf, DSYTRF)
    #define dsytrs_base LAPACK_GLOBAL(dsytrs, DSYTRS)
    #define dgelss_base LAPACK_GLOBAL(dgelss, DGELSS)

    void dpotrf_base(const char *uplo, const int *n, double *a, const int *lda, int *info);
    void dpotrs_base(const char *uplo, const int *n, const int *nrhs, const double *a,
                const int *lda, double *b, const int *ldb, int *info);
    void dpocon_base(const char *uplo, const int *n, const double *a, const int *lda,
                const double *anorm, double *rcond, double *work, int *iwork,
                int *info);
    void dgetrf_base(const int *m, const int *n, double *a, const int *lda, int *ipiv,
                int *info);
    void dgetrs_base(const char *trans, const int *n, const int *nrhs, const double *a,
                const int *lda, const int *ipiv, double *b, const int *ldb,
                int *info);
    void dsytrf_base(const char *uplo, const int *n, double *a, const int *lda,
                int *ipiv, double *work, const int *lwork, int *info);
    void dsytrs_base(const char *uplo, const int *n, const int *nrhs, const double *a,
                const int *lda, const int *ipiv, double *b, const int *ldb, int *info);
    void dgelss_base(const int *m, const int *n, const int *nrhs, double *a, const int *lda,
                double *b, const int *ldb, double *s, const double *rcond, int *rank,
                double *work, const int *lwork, int *info);

    #define dpotrf(...) dpotrf_base(__VA_ARGS__)
    #define dpotrs(...) dpotrs_base(__VA_ARGS__)
    #define dpocon(...) dpocon_base(__VA_ARGS__)
    #define dgetrf(...) dgetrf_base(__VA_ARGS__)
    #define dgetrs(...) dgetrs_base(__VA_ARGS__)
    #define dsytrf(...) dsytrf_base(__VA_ARGS__)
    #define dsytrs(...) dsytrs_base(__VA_ARGS__)
    #define dgelss(...) dgelss_base(__VA_ARGS__)

#ifdef __cplusplus
}
#endif
#endif