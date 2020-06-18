/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 * Formula =  B * At +  A^2 * B
 */
double* my_solver(int N, double *A, double* B) {
	
	double *At = (double *)calloc(N * N, sizeof(double));
	double *fst = (double *)calloc(N * N, sizeof(double));
	double *snd = (double *)calloc(N * N, sizeof(double));
	double *trd = (double *)calloc(N * N, sizeof(double));
	int i, j, k;

	
	/* calcul Atranspus */
	for (i = 0; i < N; i++) {
		register double *orig_At = &At[N * i];
		for (j = 0; j <= i; j++){
			orig_At[j] = A[j * N + i];
		} 
	}

	/* calcul B x At */
	for (i = 0; i < N; i++) {
		register double *orig_B = &B[i * N];
		register double *orig_fst = &fst[i * N]; 
		for (k = 0; k < N; k++) {
			register double idxB = orig_B[k];
			register double *orig_At = &At[k * N];
			for (j = 0; j <= k; j++) {
				orig_fst[j] += idxB * orig_At[j];
			}
		}
	}

	/* calcul A x A */
	for (i = 0; i < N; i++) {
		register double *orig_A1 = &A[i * N];
		register double *orig_snd = &snd[i * N]; 
		for (k = i; k < N; k++) {
			register double idxA = orig_A1[k];
			register double *orig_A2 = &A[k * N];
			for (j = 0; j < N; j += 8) {
				orig_snd[j] += idxA * orig_A2[j];
				orig_snd[j + 1] += idxA * orig_A2[j + 1];
				orig_snd[j + 2] += idxA * orig_A2[j + 2];
				orig_snd[j + 3] += idxA * orig_A2[j + 3];
				orig_snd[j + 4] += idxA * orig_A2[j + 4];
				orig_snd[j + 5] += idxA * orig_A2[j + 5];
				orig_snd[j + 6] += idxA * orig_A2[j + 6];
				orig_snd[j + 7] += idxA * orig_A2[j + 7];
			
			}
		}
	}

	/* calcul A^2 x B */
	for (i = 0; i < N; i++) {
		register double *orig_A2 = &snd[i * N];
		register double *orig_trd = &trd[i * N]; 
		for (k = i; k < N; k++) {
			register double idxA2 = orig_A2[k];
			register double *orig_B = &B[k * N];
			for (j = 0; j < N; j += 8) {
				orig_trd[j] += idxA2 * orig_B[j];
				orig_trd[j + 1] += idxA2 * orig_B[j + 1];
				orig_trd[j + 2] += idxA2 * orig_B[j + 2];
				orig_trd[j + 3] += idxA2 * orig_B[j + 3];
				orig_trd[j + 4] += idxA2 * orig_B[j + 4];
				orig_trd[j + 5] += idxA2 * orig_B[j + 5];
				orig_trd[j + 6] += idxA2 * orig_B[j + 6];
				orig_trd[j + 7] += idxA2 * orig_B[j + 7];
			}
		}
	}

	/* suma intre fst si trd */
	for (i = 0; i < N; i++) {
		register double *orig_fst = &fst[N * i];
		register double *orig_trd = &trd[N * i];
		for (j = 0; j < N; j += 8) {
			orig_fst[j] += orig_trd[j];
			orig_fst[j+1] += orig_trd[j+1];
			orig_fst[j+2] += orig_trd[j+2];
			orig_fst[j+3] += orig_trd[j+3];
			orig_fst[j+4] += orig_trd[j+4];
			orig_fst[j+5] += orig_trd[j+5];
			orig_fst[j+6] += orig_trd[j+6];
			orig_fst[j+7] += orig_trd[j+7];
		}
	}

	free(At);
	free(snd);
	free(trd);
	return fst;
}
