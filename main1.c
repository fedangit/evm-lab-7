#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


void CopyMatrix(float* Matrix1, float* Matrix2, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix2[i * n + j] = Matrix1[i * n + j];
        }
    }
}

void FillMatrix(float *Matrix, int n) {
    srand(time(NULL));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix[i * n + j] = rand() % 5;
        }
    }
}

float* SumMatrix(float* Matrix1, float* Matrix2, int n) {
    float *MatrixResult = (float *)calloc(n * n, sizeof(float));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixResult[i * n + j] = Matrix1[i * n + j] + Matrix2[i * n + j];
        }
    }

    return MatrixResult;
}

float *SubMatrix(float *Matrix1, float *Matrix2, int n) {
    float *MatrixResult = (float *)calloc(n * n, sizeof(float));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixResult[i * n + j] = Matrix1[i * n + j] - Matrix2[i * n + j];
        }
    }

    return MatrixResult;
}

float *DivMatrix(float *Matrix1, float digit, int n) {
    float *MatrixResult = (float *)calloc(n * n, sizeof(float));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixResult[i * n + j] = Matrix1[i * n + j] / digit;
        }
    }

    return MatrixResult;
}

void MultMatrix(float *Matrix1, float *Matrix2, float* MatrixResult, int n) {
    for (int i = 0; i < n; i++) {
        float *c = MatrixResult + i * n;

        for (int j = 0; j < n; j++) {
            c[j] = 0;
        }

        for (int k = 0; k < n; k++) {
            const float *b = Matrix2 + k * n;
            float a = Matrix1[i * n + k];

            for (int j = 0; j < n; j++) {
                c[j] += a * b[j];
            }
        }
    }
}

void InitializationMatrix(float* Matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix[i * n + j] = 0;
        }
    }
}

void FillIdentityMatrix(float *Matrix, int n) {
    for (int i = 0; i < n; i++) {
       Matrix[i * n + i] = 1;
    }
}

float UnitRate(float* Matrix, int n) {
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

float* TranspositionMatrix(float* Matrix, int n){
    float* MatrixResult = (float *)calloc(n * n, sizeof(float));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MatrixResult[j * n + i] = Matrix[i * n + j];
            int a = j * n + i;
            int b = i * n + j;
        }
    }

    return MatrixResult;
}

void InversionMatrix(float *Matrix, float* MatrixResult, int n, int m) {
    float *MatrixIdentity = (float *)calloc(n * n, sizeof(float));
    FillIdentityMatrix(MatrixIdentity, n);

    float digit = UnitRate(Matrix, n) * InfiniteRate(Matrix, n);
    float *MatrixB = DivMatrix(TranspositionMatrix(Matrix, n), digit, n);

    float *Matrix1 = (float *)calloc(n * n, sizeof(float));
    MultMatrix(MatrixB, Matrix, Matrix1, n);

    float *MatrixC = (float *)calloc(n * n, sizeof(float));
    FillIdentityMatrix(MatrixC, n);

    float* MatrixR = (float *)calloc(n * n, sizeof(float));
    float *MatrixRFix = SubMatrix(MatrixIdentity, Matrix1, n);
    CopyMatrix(MatrixRFix, MatrixR, n);

    float *MatrixRCopy = (float *)calloc(n * n, sizeof(float));
    for (int i = 1; i < m; i++) {
        MatrixC = SumMatrix(MatrixC, MatrixR, n);
        CopyMatrix(MatrixR, MatrixRCopy, n);
        MultMatrix(MatrixRCopy, MatrixRFix, MatrixR, n);
    }

    MultMatrix(MatrixC, MatrixB, MatrixResult, n);

    free(MatrixRFix);
    free(MatrixB);
    free(MatrixIdentity);
    free(Matrix1);
    free(MatrixC);
    free(MatrixR);
    free(MatrixRCopy);
}

int main() {

    int n, m;
    scanf("%d%d", &n, &m);

    float *Matrix = (float *)calloc(n * n, sizeof(float));
    float* Matrix1 = (float *)calloc(n * n, sizeof(float));

    FillMatrix(Matrix1, n);

    time_t start = time(NULL), end;

    InversionMatrix(Matrix1, Matrix, n, m);

    end = time(NULL);

    printf("Time: %ld seconds\n", end - start);

    float *Matrix2 = (float *)calloc(n * n, sizeof(float));
    MultMatrix(Matrix1, Matrix, Matrix2, n);

    printf("Second standart: %f\n", UnitRate(Matrix2, n));

    free(Matrix);
    free(Matrix1);
    free(Matrix2);

    return 0;
}
