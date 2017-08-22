
#ifndef ATLAS_HDR__
#define ATLAS_HDR__
#include "config.h"

extern "C" {
#include "atlas/clapack.h"

// FIXME: external header (where?)
int dgelss_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int
    *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int
    *info);

void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
    double *A, const int *lda, double *b, const int *ldb, double *work, const
    int * lwork, int *info);
}

#ifndef HAVE_LIBLAPACK
#warning "this will not work, install lapack to get full functionality"
int dgelss_(int *, int *, int *, double *, int *, double *, int *, double *,
		double *, int *, double *, int *, int *){
	incomplete();
	return 1;
}
#endif

#endif
