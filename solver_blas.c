/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"
#include "cblas.h"

/* 
 * Add your BLAS implementation here
 */
/*
 * FORMULA = B x At + A^2 x B
 */
double* my_solver(int N, double *A, double *B) {

	double *fst = (double *)calloc(N * N, sizeof(double));
	double *snd = (double *)calloc(N * N, sizeof(double));
	double *trd = (double *)calloc(N * N, sizeof(double));

	int i;
	int j;
	
	i = 0;
	j = 0;

	for (i = 0; i < N * N; i++) {
		fst[i] = B[i];
		snd[i] = A[i];
		trd[i] = B[i];
	}

	// calcul B x At
	cblas_dtrmm(CblasRowMajor, CblasRight , CblasUpper, CblasTrans,
	CblasNonUnit, N, N, 1, A, N, fst, N);

	// calcul A x A
	cblas_dtrmm(CblasRowMajor, CblasLeft , CblasUpper, CblasNoTrans,
	CblasNonUnit, N, N, 1, A, N, snd, N);

	// calcul  A^2 x B
	cblas_dtrmm(CblasRowMajor, CblasLeft , CblasUpper, CblasNoTrans,
	CblasNonUnit, N, N, 1, snd, N, trd, N);
	
	/* suma celor calculate fst + trd */
	for (i = 0; i < N; i++) {
		register double *orig_fst = &fst[N * i];
		register double *orig_trd = &trd[N * i];
		for (j = 0; j < N; j++) {
			orig_fst[j] += orig_trd[j];
		}
	}

	free(snd);
	free(trd);
	return fst;
}
