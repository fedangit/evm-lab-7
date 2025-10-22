#include <limits.h>
#include <immintrin.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <cblas.h>


void CopyMatrix(float *Matrix, float *MatrixC, int n) {
    cblas_scopy(n * n, Matrix, 1, MatrixC, 1);
}

void FillMatrix(float *Matrix, int n) {
    srand(time(NULL));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix[i * n + j] = rand() % 5;
        }
    }
}

void SumMatrix(float *Matrix1, float *Matrix2, float *MatrixResult, int n) {
    cblas_scopy(n * n, Matrix1, 1, MatrixResult, 1);
    cblas_saxpy(n * n, 1.0, Matrix2, 1, MatrixResult, 1);
}

void SubMatrix(float *Matrix1, float *Matrix2, float *MatrixResult, int n) {
    cblas_scopy(n * n, Matrix1, 1, MatrixResult, 1);
    cblas_saxpy(n * n, -1.0, Matrix2, 1, MatrixResult, 1);
}

void DivMatrix(const float *Matrix, float *MatrixResult, float digit, int n) {
    cblas_scopy(n * n, Matrix, 1, MatrixResult, 1);

    for (int i = 0; i < n * n; i++) {
        MatrixResult[i] /= digit;
    }
}

void InitializationMatrix(float *Matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix[i * n + j] = 0;
        }
    }
}

void FillIdentityMatrix(float *Matrix, int n) {
    InitializationMatrix(Matrix, n);
    for (int i = 0; i < n; i++) {
        Matrix[i * n + i] = 1;
    }
}

float UnitRate(float *Matrix, int n) {
    float result = FLT_MIN, sum = 0;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            sum += abs(Matrix[i * n + j]);
        }

        if (sum > result) {
            result = sum;
        }

        sum = 0;
    }

    return result;
}

float InfiniteRate(float *Matrix, int n) {
    float result = 0, sum = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += fabs(Matrix[i * n + j]);
        }

        if (sum > result) {
            result = sum;
        }

        sum = 0;
    }

    return result;
}

void TranspositionMatrix(const float *Matrix, float *MatrixTransposition, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixTransposition[j * n + i] = Matrix[i * n + j];
        }
    }
}

void MultMatrix(const float *A, const float *B, float *C, int N) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, A, N, B, N, 0.0, C, N);
}

void InversionMatrix(float *Matrix, float *MatrixResult, int n, int m) {
    float *MatrixRFix = (float *)_mm_malloc(n * n * sizeof(float), 16);
    float *MatrixB = (float *)_mm_malloc(n * n * sizeof(float), 16);
    float *MatrixIdentity = (float *)_mm_malloc(n * n * sizeof(float), 16);

    FillIdentityMatrix(MatrixIdentity, n);

    float digit = UnitRate(Matrix, n) * InfiniteRate(Matrix, n);
    float *MatrixTransposition = (float *)_mm_malloc(n * n * sizeof(float), 16);

    TranspositionMatrix(Matrix, MatrixTransposition, n);
    DivMatrix(MatrixTransposition, MatrixB, digit, n);

    float *Matrix1  = (float *)_mm_malloc(n * n * sizeof(float), 16);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, MatrixB, n, Matrix, n, 0.0, Matrix1, n);
    SubMatrix(MatrixIdentity, Matrix1, MatrixRFix, n);

    float *MatrixC  = (float *)_mm_malloc(n * n * sizeof(float), 16);
    FillIdentityMatrix(MatrixC, n);

    float *MatrixR  = (float *)_mm_malloc(n * n * sizeof(float), 16);
    CopyMatrix(MatrixRFix, MatrixR, n);

    for (int i = 1; i < m; i++) {
        SumMatrix(MatrixC, MatrixR, MatrixC, n);
        float *MatrixRCopy = (float *)_mm_malloc(n * n * sizeof(float), 16);
        CopyMatrix(MatrixR, MatrixRCopy, n);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, MatrixRCopy, n, MatrixRFix, n, 0.0, MatrixR, n);
    }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, MatrixC, n, MatrixB, n, 0.0, MatrixResult, n);

    _mm_free(MatrixRFix);
    _mm_free(MatrixB);
    _mm_free(MatrixIdentity);
    _mm_free(Matrix1);
    _mm_free(MatrixC);
    _mm_free(MatrixR);
    _mm_free(MatrixTransposition);
}

int main() {
    int n, m;
    scanf("%d%d", &n, &m);

    float *Matrix1 = (float *)_mm_malloc(n * n * sizeof(float), 16);
    float *MatrixResult = (float *)_mm_malloc(n * n * sizeof(float), 16);
    float *MatrixCheck = (float *) _mm_malloc(n * n * sizeof(float), 16);

    FillMatrix(Matrix1, n);

    time_t start = time(NULL), end;

    InversionMatrix(Matrix1, MatrixResult, n, m);

    end = time(NULL);

    printf("Time: %ld seconds\n", end - start);
    
    MultMatrix(Matrix1, MatrixResult, MatrixCheck, n);

    printf("Second standart: %f\n", UnitRate(MatrixCheck, n));

    _mm_free(MatrixResult);
    _mm_free(Matrix1);
    _mm_free(MatrixCheck);

    return 0;
}
