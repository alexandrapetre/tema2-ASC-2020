/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 * Formula =  B x At + A^2 x B
 */
double* my_solver(int N, double *A, double* B) {
	
	double *At = (double *)calloc(N * N, sizeof(double));
	double *fst = (double *)calloc(N * N, sizeof(double));
	double *snd = (double *)calloc(N * N, sizeof(double));
	double *trd = (double *)calloc(N * N, sizeof(double));
	int i;
	int j;
	int k;
	
	/* calcul Atranspus */
	for (i = 0; i < N; i++) {
		register double *orig_At = &At[N * i];
		for (j = 0; j < N; j++){
			if (j <= i) {
				orig_At[j] = A[j * N + i];
			} else 
				break;
		} 
	}

	/* calcul B x At */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			register int idx = j + i * N;
			register double sum = 0.0;
			for (k = 0; k < N; k++) {
				if (k >= j) {
					register int idxB = k + N * i;
					register int idxAt = j + k * N;
					sum += B[idxB] * At[idxAt];
				} else
					continue;
			}
			fst[idx] = sum;
		}
	}

	/* calcul A x A */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			register int idx = j + i * N;
			register double sum = 0.0;
			for (k = 0; k < N; k++) {
				if (i <= j) {
					register int idxA1 = i * N + k;
					register int idxA2 = k * N + j;
					sum += A[idxA1] * A[idxA2];
				} else
					break;
			}
			snd[idx] = sum;
		}
	}

	/* calcul A^2 x B */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			register int idx = j + i * N;
			register double sum = 0.0;
			for (k = 0; k < N; k++) {
				if (i <= k) {
					register int idxA2 = i * N + k;
					register int idxB = k * N + j;
					sum += snd[idxA2] * B[idxB];
				} else 
					continue;
			}
			trd[idx] = sum;
		}
	}

	/* suma celor calculate fst + trd */
	for (i = 0; i < N; i++) {
		register double *orig_fst = &fst[N * i];
		register double *orig_trd = &trd[N * i];
		for (j = 0; j < N; j++) {
			orig_fst[j] += orig_trd[j];
		}
	}

	free(At);
	free(snd);
	free(trd);
	return fst;
}
