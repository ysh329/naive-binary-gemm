#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include "gemm_common.h"

#define ENCODE_BIT 32


void generate_matrix(unsigned int *mat, int mat_len) {
    srand( (unsigned) time(0) );
    for (int i = 0; i < mat_len; i++) {
        mat[i] = rand() % 2;
    }
}

void ones_matrix(unsigned int *mat, int mat_len) {
    for (int i = 0; i < mat_len; i++) {
        mat[i] = 1;
    }
}

void ones_matrix(int *mat, int mat_len) {
    for (int i = 0; i < mat_len; i++) {
        mat[i] = 1;
    }
}

void random_matrix(int *mat, int mat_len) {
    for (int i = 0; i < mat_len; i++) {
        mat[i] = rand()%2==1 ? 1: -1;
    }
}

void print_matrix(unsigned int *mat, int mat_len) {
    for (int i = 0; i < mat_len; i++) {
        printf("%d ", mat[i]);
    }
    printf("\n");
}

void print_matrix(unsigned int *mat, int m, int n) {
    for (int i = 0; i < m; i++) {
    	for(int j = 0; j < n; j++) {
    		printf("%d ", mat[j + i * n]);
    	}
        printf("\n");
    }
    printf("\n");
}

void print_matrix(int *mat, int m, int n) {
    for (int i = 0; i < m; i++) {
    	for(int j = 0; j < n; j++) {
    		printf("%d ", mat[j + i * n]);
    	}
        printf("\n");
    }
    printf("\n");
}

void print_matrix(float *mat, int m, int n) {
    for (int i = 0; i < m; i++) {
    	for(int j = 0; j < n; j++) {
    		printf("%.0f ", mat[j + i * n]);
    	}
        printf("\n");
    }
    printf("\n");
}

void encode_rows(int *columns, unsigned int *columns_binary, int m, int n) {
    // m: row number
    // n: col number, new col number
    int col_bin_n = n / ENCODE_BIT; // col_bin_m: new row number
    int i, j, k;

    for (i = 0; i < m; i++) {
        for (j = 0; j < col_bin_n; j++) {
            /* start index */
            int j32 = j * ENCODE_BIT;
            unsigned int rvalue = 0;
            unsigned int sign;

            for (k = 0; k < ENCODE_BIT; k++) {
                sign = columns[j32 + k+ col_bin_n * i] > 0;
                rvalue |= (sign << k);
            }
            /* store 32bit-encoded elem in i-th row j-th col*/
            columns_binary[j + col_bin_n * i] = rvalue;
        }
    }
}

void encode_cols(int *columns, unsigned int *columns_binary, int m, int n) {
    // m: row number
    // n: col number, new col number
    int col_bin_m = m / ENCODE_BIT; // col_bin_m: new row number
    int i, j, k;

    for (i = 0; i < col_bin_m; i++) {
        for (j = 0; j < n; j++) {
            /* start index */
            int i32 = i * ENCODE_BIT;
            unsigned int rvalue = 0;
            unsigned int sign;

            for (k = 0; k < ENCODE_BIT; k++) {
                sign = columns[j + n * (i32 + k)] > 0;
                rvalue |= (sign << k);
            }
            /* store 32bit-encoded elem in i-th row j-th col*/
            columns_binary[j + n * i] = rvalue;
        }
    }
}

void itos(unsigned int *i_var, float *s_var, int len) {
    for (int i = 0; i < len; i++) {
	s_var[i] = i_var[i] * 1.0;
    }
}

void itos(int *i_var, float *s_var, int len) {
    for (int i = 0; i < len; i++) {
	s_var[i] = i_var[i] * 1.0;
    }
}

void binary_gemm(int m, int n, int k, unsigned int *a, int lda, unsigned int *b, int ldb, int *c, int ldc) {
    /* column-major order */
    int i;
    unsigned int * a_ = a;
    for (i = 0; i < m; i++) {
        int j, l;
        unsigned int *b_ = b;
        for (j = 0; j < n; j++) {
            register unsigned int Cvalue = 0;
            for (l = 0; l < k; l++) {
            	int cc = popcount(b_[l] ^ a_[l * lda + i]);
                Cvalue += popcount(b_[l] ^ a_[l * lda + i]);
                // printf("%d %d %d, ",b_[l], a_[l * lda + i], cc);
            }
            // printf("%d ", Cvalue);
            b_ += ldb;
            //c[j *ldc + i] = 32 * k - 2 * Cvalue;
            c[j * ldc + i] = ENCODE_BIT * k - 2 * Cvalue;
        }
    }
}

int main(void){
    struct timeval start,finish;
    double duration;
    srand((unsigned)time(NULL));

    int m, n, k;
    m = 3;
    n = 4;
    k = 32;

    int len_a, len_b, len_c;
    len_a = m*k;
    len_b = k*n;
    len_c = m*n;

    int *a, *b, *c;
    a = (int *) malloc(len_a * sizeof(int));
    b = (int *) malloc(len_b * sizeof(int));
    c = (int *) malloc(len_c * sizeof(int));


    printf("A\n");
    // print_matrix(a, m, k);
    //generate_matrix(a, len_a);
    ones_matrix(a, len_a);
    a[0] = -1;
    a[23] = -1;
    print_matrix(a, m, k);
    printf("\n");

    printf("B\n");
    // print_matrix(b, k, n);
    //generate_matrix(b, len_b);
    ones_matrix(b, len_b);
    b[0] = -1;
    b[35] = -1;
    b[69] = -1;
    print_matrix(b, k, n);
    printf("\n");

    printf("C\n");
    // print_matrix(c, m, n);
    //generate_matrix(c, len_c);
    ones_matrix(c, len_c);
    print_matrix(c, m, n);
    printf("\n");

    /* encode a, b */
    unsigned int *encoded_a, *encoded_b, *encoded_c;
    encoded_a = (unsigned int *) malloc(len_a/ENCODE_BIT * sizeof(unsigned int));
    encoded_b = (unsigned int *) malloc(len_b/ENCODE_BIT * sizeof(unsigned int));

    encode_rows(a, encoded_a, m, k);
    encode_cols(b, encoded_b, k, n);

    print_matrix(encoded_a, m, k/ENCODE_BIT);
    print_matrix(encoded_b, k/ENCODE_BIT, n);

    /* binary gemm */
    gettimeofday(&start, NULL);
    binary_gemm(m, n, k/ENCODE_BIT,
                encoded_a, m,
                encoded_b, k/ENCODE_BIT,
                c, m);
    gettimeofday(&finish, NULL);
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + 
		(double)(finish.tv_usec-start.tv_usec)) / 1000000;
    double gflops = 2.0 * m * n * k / ENCODE_BIT / ENCODE_BIT / ENCODE_BIT;
    gflops = gflops/duration*1.0e-6;
    printf("binary gemm \n %dx%dx%d\t%lf s\t%lf MFLOPS\n", m/ENCODE_BIT, n/ENCODE_BIT, k/ENCODE_BIT, duration, gflops);

    // binary gemm result
    printf("binary_gemm result:\n");
    print_matrix(c, m, n);

    /***************
     * cblas sgemm *
     ***************/

    float *sa, *sb, *sc;
    sa = (float *) malloc(len_a * sizeof(float));
    sb = (float *) malloc(len_b * sizeof(float));
    sc = (float *) malloc(len_c * sizeof(float));

    itos(a, sa, len_a);
    itos(b, sb, len_b);

    // precision benchmark
    printf("[cblas gemm]precision benchmark for un-encoded matrix\n");
    print_matrix(sa, m, k);
    print_matrix(sb, k, n);

    cblas_sgemm(CblasColMajor,
		CblasNoTrans,
		CblasTrans,
		m, n, k,
		1, sa, m,
		   sb, n,
		0, sc, m);
    print_matrix(sc, m, n);

    // speed benchmark
    printf("[cblas gemm]speed benchmark for encoded matrix\n");
    gettimeofday(&start, NULL);
    cblas_sgemm(CblasColMajor,
		CblasNoTrans,
		CblasTrans,
		m,
		n,
		k/ENCODE_BIT,
		1, sa, m,
		   sb, k,
	        2, sc, k/ENCODE_BIT);
    gettimeofday(&finish, NULL);
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 +
                (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    gflops = gflops/duration*1.0e-6;
    printf("cblas sgemm \n %dx%dx%d\t%lf s\t%lf MFLOPS\n", m/ENCODE_BIT, n/ENCODE_BIT, k/ENCODE_BIT, duration, gflops);

    // cblas sgemm result is meanless

    /* free */
    free(a); free(sa); free(encoded_a);
    free(b); free(sb); free(encoded_b);
    free(c); free(sc);
    return 0;
}
