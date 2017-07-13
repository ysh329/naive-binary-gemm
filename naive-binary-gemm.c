#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <cblas.h>
#include "gemm_common.h"

#define ENCODE_BIT 32

void generate_matrix(unsigned int *mat, int mat_len);
void ones_matrix(unsigned int *mat, int mat_len);
void binary_gemm(int m, int n, int k, unsigned int *a, int lda, unsigned int *b, int ldb, unsigned int *c, int ldc);
void print_matrix(unsigned int *mat, int mat_len);

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

void print_matrix(unsigned int *mat, int mat_len) {
    for (int i = 0; i < mat_len; i++) {
        printf("%d ", mat[i]);
    }
    printf("\n");
}

void encode_cols(unsigned int *columns, unsigned int *columns_binary, int m, int n) {
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

            for (k = 0; k < 1; k++) {
                sign = columns[j + n * (i32 + k)];
                rvalue |= (sign << k);
            }
            /* store 32bit-encoded elem in i-th row j-th col*/
            columns_binary[j + n * i] = rvalue;
        }
    }
}

void binary_gemm(int m, int n, int k, unsigned int *a, int lda, unsigned int *b, int ldb, unsigned int *c, int ldc) {
    /* column-major order */
    int i;
    unsigned int * a_ = a;
    for (i = 0; i < m; i++) {
        int j, l;
        unsigned int *b_ = b;
        for (j = 0; j < n; j++) {
            register unsigned int Cvalue = 0;
            for (l = 0; l < k; l++) {
                Cvalue += popcount(b_[l] ^ a_[l * lda + i]);
            }
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
    n = 3;
    k = 2;

    int len_a, len_b, len_c;
    len_a = m*k;
    len_b = k*n;
    len_c = m*n/ENCODE_BIT/ENCODE_BIT;

    unsigned int *a, *b, *c, *encoded_a, *encoded_b, *encoded_c;
    a = (unsigned int *) malloc(len_a * sizeof(unsigned int));
    b = (unsigned int *) malloc(len_b * sizeof(unsigned int));
    c = (unsigned int *) malloc(len_c * sizeof(unsigned int));

    encoded_a = (unsigned int *) malloc(len_a/ENCODE_BIT * sizeof(unsigned int));
    encoded_b = (unsigned int *) malloc(len_b/ENCODE_BIT * sizeof(unsigned int));
    encoded_c = (unsigned int *) malloc(len_c/ENCODE_BIT/ENCODE_BIT * sizeof(unsigned int));

    printf("A\n");
    print_matrix(a, len_a);
    //generate_matrix(a, len_a);
    ones_matrix(a, len_a);
    print_matrix(a, len_a);
    printf("\n");

    printf("B\n");
    print_matrix(b, len_b);
    //generate_matrix(b, len_b);
    ones_matrix(b, len_b);
    print_matrix(b, len_b);
    printf("\n");

    printf("C\n");
    print_matrix(c, len_c);
    //generate_matrix(c, len_c);
    ones_matrix(c, len_c);
    print_matrix(c, len_c);
    printf("\n");

    /* encode a, b */
    encode_cols(a, encoded_a, m, k);
    encode_cols(b, encoded_b, k, n);

    print_matrix(encoded_a, len_a/ENCODE_BIT);
    print_matrix(encoded_b, len_b/ENCODE_BIT);

    /* binary gemm */
    gettimeofday(&start, NULL);
    binary_gemm(m/ENCODE_BIT, n/ENCODE_BIT, k/ENCODE_BIT,
                encoded_a, len_a/ENCODE_BIT,
                encoded_b, len_b/ENCODE_BIT,
                encoded_c, len_c/ENCODE_BIT/ENCODE_BIT);
    gettimeofday(&finish, NULL);
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + 
		(double)(finish.tv_usec-start.tv_usec)) / 1000000;
    double gflops = 2.0 * m * n * k / ENCODE_BIT / ENCODE_BIT / ENCODE_BIT;
    gflops = gflops/duration*1.0e-6;
    printf("binary gemm \n %dx%dx%d\t%lf s\t%lf MFLOPS\n", m/ENCODE_BIT, n/ENCODE_BIT, k/ENCODE_BIT, duration, gflops);

    // binary gemm result
    printf("binary_gemm result:\n");
    print_matrix(encoded_c, len_c/ENCODE_BIT/ENCODE_BIT);

    /* cblas sgemm */
    gettimeofday(&start, NULL);
    cblas_sgemm(CblasColMajor,
		CblasNoTrans,
		CblasTrans,
		m/ENCODE_BIT,
		n/ENCODE_BIT,
		k/ENCODE_BIT, alpha, encoded_a, m, encoded_b, k, beta, encoded_c, k);
    gettimeofday(&finish, NULL);
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 +
                (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    gflops = gflops/duration*1.0e-6;
    printf("cblas sgemm \n %dx%dx%d\t%lf s\t%lf MFLOPS\n", m/ENCODE_BIT, n/ENCODE_BIT, k/ENCODE_BIT, duration, gflops);

    // cblas sgemm result is meanless
    free(a); free(encoded_a);
    free(b); free(encoded_b);
    free(c); free(encoded_c);
    return 0;
}
