#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>

float** matrix_create(int n, int m, int initialize) {
    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix dimensions\n");
        return NULL;
    }

    float** matrix = (float**)malloc(n * sizeof(float*));
    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    float* data = (float*)calloc(n * m, sizeof(float));
    if (!data) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(matrix);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matrix[i] = data + i * m;
    }
    return matrix;
}

float** matrix_input(int n, int m) {
    float** matrix = matrix_create(n, m, 0);
    if (!matrix) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (scanf("%f", &matrix[i][j]) != 1) {
                fprintf(stderr, "Error: Invalid input\n");
                matrix_free(matrix);
                return NULL;
            }
        }
    }
    return matrix;
}

float** matrix_rcreate(int n, int m, int min, int max) {
    if (min > max) {
        fprintf(stderr, "Error: Invalid range\n");
        return NULL;
    }

    float** matrix = matrix_create(n, m, 1);
    if (!matrix) return NULL;

    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = (float)(rand() % (max - min + 1)) + min;
        }
    }
    return matrix;
}

float** matrix_copy(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** copy = matrix_create(n, m, 0);
    if (!copy) return NULL;

    for (int i = 0; i < n; i++) {
        memcpy(copy[i], matrix[i], m * sizeof(float));
    }
    return copy;
}

void matrix_free(float** matrix) {
    if (!matrix) return;
    free(matrix[0]);
    free(matrix);
}

void matrix_print(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix\n");
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%8.2f", matrix[i][j]);
        }
        printf("\n");
    }
}

int matrix_savef(float** matrix, int n, int m, const char* filename) {
    if (!matrix || n <= 0 || m <= 0 || !filename) {
        fprintf(stderr, "Error: Invalid input\n");
        return -1;
    }

    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error");
        return -1;
    }

    fprintf(file, "%d %d\n", n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fprintf(file, "%.6f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}

float** matrix_loadf(const char* filename) {
    if (!filename) {
        fprintf(stderr, "Error: Invalid filename\n");
        return NULL;
    }

    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error");
        return NULL;
    }

    int n, m;
    if (fscanf(file, "%d %d", &n, &m) != 2 || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid file format\n");
        fclose(file);
        return NULL;
    }

    float** matrix = matrix_create(n, m, 0);
    if (!matrix) {
        fclose(file);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (fscanf(file, "%f", &matrix[i][j]) != 1) {
                fprintf(stderr, "Error: Invalid file format\n");
                matrix_free(matrix);
                fclose(file);
                return NULL;
            }
        }
    }

    fclose(file);
    return matrix;
}

float** matrix_sum(float** A, float** B, int n, int m) {
    if (!A || !B || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_copy(A, n, m);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] += B[i][j];
        }
    }
    return result;
}

float** matrix_multp(float** A, float** B, int n, int p, int m) {
    if (!A || !B || n <= 0 || p <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_create(n, m, 1);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

float matrix_determin(float** matrix, int n) {
    if (!matrix || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return 0;
    }

    if (n == 1) return matrix[0][0];

    float det = 0;
    float** submatrix = matrix_create(n-1, n-1, 0);
    if (!submatrix) return 0;

    for (int x = 0; x < n; x++) {
        int subi = 0;
        for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
                if (j == x) continue;
                submatrix[subi][subj] = matrix[i][j];
                subj++;
            }
            subi++;
        }
        det += (x % 2 == 0 ? 1 : -1) * matrix[0][x] * matrix_determin(submatrix, n-1);
    }

    matrix_free(submatrix);
    return det;
}

float** matrix_inver(float** matrix, int n) {
    if (!matrix || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float det = matrix_determin(matrix, n);
    if (fabs(det) < 1e-10) {
        fprintf(stderr, "Error: Matrix is singular\n");
        return NULL;
    }

    float** inverse = matrix_create(n, n, 0);
    if (!inverse) return NULL;

    float** temp = matrix_create(n-1, n-1, 0);
    if (!temp) {
        matrix_free(inverse);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int ti = 0;
            for (int x = 0; x < n; x++) {
                if (x == i) continue;
                int tj = 0;
                for (int y = 0; y < n; y++) {
                    if (y == j) continue;
                    temp[ti][tj] = matrix[x][y];
                    tj++;
                }
                ti++;
            }
            float cofactor = ((i + j) % 2 == 0 ? 1 : -1) * matrix_determin(temp, n-1);
            inverse[j][i] = cofactor / det;
        }
    }

    matrix_free(temp);
    return inverse;
}