#include <limits.h>
#include <immintrin.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


void CopyMatrix(float *Matrix, float *MatrixC, int n) {
    __m128 a;

    for (int i = 0; i < n * n; i += 4) {
        a = _mm_load_ps(Matrix + i);
        _mm_store_ps(MatrixC + i, a);
    }
}

void FillMatrix(float *Matrix, int n) {
    srand(time(0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix[i * n + j] = rand() % 5;
        }
    }
}

float *SumMatrix(float *Matrix1, float *Matrix2, int n) {
    float *MatrixResult = (float *)_mm_malloc(n * n * sizeof(float), 16);
    __m128 a, b, c;

    for (int i = 0; i < n * n; i += 4) {
        a = _mm_load_ps(Matrix1 + i);
        b = _mm_load_ps(Matrix2 + i);
        c = _mm_add_ps(a, b);
        _mm_store_ps(MatrixResult + i, c);
    }

    return MatrixResult;
}

void SubMatrix(float *Matrix1, float *Matrix2, float *MatrixResult, int n) {
    __m128 a, b, c;
    for (int i = 0; i < n * n; i += 4)
    {
        a = _mm_load_ps(Matrix1 + i);
        b = _mm_load_ps(Matrix2 + i);
        c = _mm_sub_ps(a, b);
        _mm_store_ps(MatrixResult + i, c);
    }
}

void DivMatrix(float *Matrix1, float *MatrixResult, float digit, int n) {
    __m128 a, b = _mm_set1_ps(digit), c;

    for (int i = 0; i < n * n; i += 4) {
        a = _mm_load_ps(Matrix1 + i);
        c = _mm_div_ps(a, b);
        _mm_store_ps(MatrixResult + i, c);
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
    float result = 0, sum = 0;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            sum += fabs(Matrix[i * n + j]);
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
            sum += abs(Matrix[i * n + j]);
        }

        if (sum > result) {
            result = sum;
        }

        sum = 0;
    }

    return result;
}

void TranspositionMatrix(float *Matrix, float *MatrixResult, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixResult[j * n + i] = Matrix[i * n + j];
            int a = j * n + i;
            int b = i * n + j;
            int v = 5;
        }
    }
}

void MultMatrix(const float *A, const float *B, float *C, int N) {
    for (int i = 0; i < N; i++) {
        float *c = C + i * N;

        for (int j = 0; j < N; j += 4) {
            _mm_store_ps(c + j + 0, _mm_setzero_ps());
        }

        for (int k = 0; k < N; k++) {
            const float *b = B + k * N;

            __m128 a = _mm_set1_ps(A[i * N + k]);

            for (int j = 0; j < N; j += 4) {
                __m128 b_val = _mm_load_ps(b + j);
                __m128 c_val = _mm_load_ps(c + j);
                __m128 result = _mm_add_ps(_mm_mul_ps(a, b_val), c_val);
                _mm_store_ps(c + j, result);
            }
        }
    }
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

    float *Matrix1 = (float *)_mm_malloc(n * n * sizeof(float), 16);
    MultMatrix(MatrixB, Matrix, Matrix1, n);
    SubMatrix(MatrixIdentity, Matrix1, MatrixRFix, n);

    float *MatrixC = (float *)_mm_malloc(n * n * sizeof(float), 16);
    FillIdentityMatrix(MatrixC, n);

    float *MatrixR = (float *)_mm_malloc(n * n * sizeof(float), 16);
    CopyMatrix(MatrixRFix, MatrixR, n);

    for (int i = 1; i < m; i++) {
        MatrixC = SumMatrix(MatrixC, MatrixR, n);

        float *MatrixRCopy = (float *)_mm_malloc(n * n * sizeof(float), 16);
        CopyMatrix(MatrixR, MatrixRCopy, n);
        MultMatrix(MatrixRCopy, MatrixRFix, MatrixR, n);
    }

    MultMatrix(MatrixC, MatrixB, MatrixResult, n);

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
