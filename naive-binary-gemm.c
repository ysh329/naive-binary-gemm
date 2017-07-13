#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gemm_common.h"

#define ENCONDE_BIT 32

void generate_matrix(unsigned int *mat, int mat_len);
void ones_matrix(unsigned int *mat, int mat_len);
void binary_gemm(int m, int n, int k, unsigned int *a, int lda, unsigned int *b, int ldb, int *c, int ldc);
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
    int m, n, k;
    m = 1;
    n = 2;
    k = 3;

    int len_a, len_b, len_c;
    len_a = m*n;
    len_b = n*k;
    len_c = m*k/ENCODE_BIT/ENCODE_BIT;

    unsigned int *a, *b, *c, *encoded_a, *encoded_b;
    a = (unsigned int *) malloc(len_a * sizeof(unsigned int));
    b = (unsigned int *) malloc(len_b * sizeof(unsigned int));
    c = (unsigned int *) malloc(len_c * sizeof(unsigned int));

    encoded_a = (unsigned int *) malloc(len_a/ENCODE_BIT * sizeof(unsigned int));
    encoded_b = (unsigned int *) malloc(len_b/ENCODE_BIT * sizeof(unsigned int));

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
    encode_cols(a, encoded_a, m, n);
    encode_cols(b, encoded_b, n, k);

    print_matrix(encoded_a, m*n/ENCODE_BIT);
    print_matrix(encoded_b, n*k/ENCODE_BIT);

    binary_gemm(m, n, k,
                encoded_a, len_a/ENCODE_BIT,
                encoded_b, len_b/ENCODE_BIT,
                c, len_c);

    printf("binary_gemm result:\n");
    print_matrix(c, len_c);

    return 0;
}