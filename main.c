#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

void mergeAsc(float arr[], int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    float* L = (float*)malloc(n1 * sizeof(float));
    float* R = (float*)malloc(n2 * sizeof(float));

    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];

    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k++] = L[i++];
    }

    while (j < n2) {
        arr[k++] = R[j++];
    }

    free(L);
    free(R);
}

void mergeSortAsc(float arr[], int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        mergeSortAsc(arr, left, mid);
        mergeSortAsc(arr, mid + 1, right);
        mergeAsc(arr, left, mid, right);
    }
}


void mergeDesc(float arr[], int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    float* L = (float*)malloc(n1 * sizeof(float));
    float* R = (float*)malloc(n2 * sizeof(float));

    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];

    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L[i] >= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k++] = L[i++];
    }

    while (j < n2) {
        arr[k++] = R[j++];
    }

    free(L);
    free(R);
}

void mergeSortDesc(float arr[], int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        mergeSortDesc(arr, left, mid);
        mergeSortDesc(arr, mid + 1, right);
        mergeDesc(arr, left, mid, right);
    }
}

float** matrix_input(int n, int m) {
    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix dimensions\n");
        return NULL;
    }

    float** matrix = (float**)malloc(n * sizeof(float*));
    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matrix[i] = (float*)malloc(m * sizeof(float));
        if (!matrix[i]) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }

        for (int j = 0; j < m; j++) {
            if (scanf("%f", &matrix[i][j]) != 1) {
                fprintf(stderr, "Error: Invalid input\n");
                for (int k = 0; k <= i; k++) free(matrix[k]);
                free(matrix);
                return NULL;
            }
        }
    }
    return matrix;
}

float** matrix_createv(int n, int m) {
    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix dimensions\n");
        return NULL;
    }

    float** matrix = (float**)calloc(n, sizeof(float*));
    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matrix[i] = (float*)calloc(m, sizeof(float));
        if (!matrix[i]) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

float** matrix_rcreate(int n, int m, int min, int max) {
    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix dimensions\n");
        return NULL;
    }

    if (min > max) {
        fprintf(stderr, "Error: Invalid range\n");
        return NULL;
    }

    float** matrix = (float**)malloc(n * sizeof(float*));
    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        matrix[i] = (float*)malloc(m * sizeof(float));
        if (!matrix[i]) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }

        for (int j = 0; j < m; j++) {
            matrix[i][j] = (float)(rand() % (max - min + 1)) + min;
        }
    }
    return matrix;
}

float** matrix_inputsqr(int n) {
    return matrix_input(n, n);
}

float** matrix_createsqrv(int n) {
    return matrix_createv(n, n);
}

float** matrix_rcreatesqr(int n, int min, int max) {
    return matrix_rcreate(n, n, min, max);
}

float** matrix_copy(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** copy = (float**)malloc(n * sizeof(float*));
    if (!copy) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        copy[i] = (float*)malloc(m * sizeof(float));
        if (!copy[i]) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(copy[j]);
            free(copy);
            return NULL;
        }

        for (int j = 0; j < m; j++) {
            copy[i][j] = matrix[i][j];
        }
    }
    return copy;
}

void matrix_free(float** matrix) {
    if (!matrix) return;
    for (int i = 0; matrix[i]; i++) free(matrix[i]);
    free(matrix);
}

void matrix_print(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix\n");
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%8.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int matrix_savef(float** matrix, int n, int m, char* filename) {
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
            fprintf(file, "%f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}

float** matrix_loadf(char* filename) {
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
    if (fscanf(file, "%d %d", &n, &m) != 2) {
        fprintf(stderr, "Error: Invalid file format\n");
        fclose(file);
        return NULL;
    }

    if (n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid matrix dimensions\n");
        fclose(file);
        return NULL;
    }

    float** matrix = (float**)malloc(n * sizeof(float*));
    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(file);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matrix[i] = (float*)malloc(m * sizeof(float));
        if (!matrix[i]) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            fclose(file);
            return NULL;
        }

        for (int j = 0; j < m; j++) {
            if (fscanf(file, "%f", &matrix[i][j]) != 1) {
                fprintf(stderr, "Error: Invalid file format\n");
                for (int k = 0; k <= i; k++) free(matrix[k]);
                free(matrix);
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

    float** result = matrix_createv(n, m);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

float** matrix_sub(float** A, float** B, int n, int m) {
    if (!A || !B || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_createv(n, m);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

float** matrix_multp(float** A, float** B, int n, int p, int m) {
    if (!A || !B || n <= 0 || p <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_createv(n, m);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = 0;
            for (int k = 0; k < p; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

float** matrix_multpscl(float** matrix, int n, int m, float scalar) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_copy(matrix, n, m);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] *= scalar;
        }
    }
    return result;
}

float** matrix_power(float** matrix, int n, int p) {
    if (!matrix || n <= 0 || p < 0 || p > 102) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    if (p == 0) {
        float** identity = matrix_createv(n, n);
        if (!identity) return NULL;
        for (int i = 0; i < n; i++) identity[i][i] = 1;
        return identity;
    }

    float** result = matrix_copy(matrix, n, n);
    if (!result) return NULL;

    for (int pow = 1; pow < p; pow++) {
        float** temp = matrix_multp(result, matrix, n, n, n);
        if (!temp) {
            matrix_free(result);
            return NULL;
        }
        matrix_free(result);
        result = temp;
    }
    return result;
}

float** matrix_transpose(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float** result = matrix_createv(m, n);
    if (!result) return NULL;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = matrix[j][i];
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
    float** temp = matrix_createv(n, n);
    if (!temp) return 0;

    int sign = 1;
    for (int f = 0; f < n; f++) {
        int i = 0, j = 0;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                if (row != 0 && col != f) {
                    temp[i][j++] = matrix[row][col];
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
        det += sign * matrix[0][f] * matrix_determin(temp, n - 1);
        sign = -sign;
    }

    matrix_free(temp);
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

    float** adj = matrix_createv(n, n);
    if (!adj) return NULL;

    float** temp = matrix_createv(n, n);
    if (!temp) {
        matrix_free(adj);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int ii = 0, jj = 0;
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    if (row != i && col != j) {
                        temp[ii][jj++] = matrix[row][col];
                        if (jj == n - 1) {
                            jj = 0;
                            ii++;
                        }
                    }
                }
            }
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = sign * matrix_determin(temp, n - 1);
        }
    }

    matrix_free(temp);

    float** inverse = matrix_createv(n, n);
    if (!inverse) {
        matrix_free(adj);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = adj[i][j] / det;
        }
    }

    matrix_free(adj);
    return inverse;
}

float matrix_trace(float** matrix, int n) {
    if (!matrix || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return 0;
    }

    float trace = 0;
    for (int i = 0; i < n; i++) {
        trace += matrix[i][i];
    }
    return trace;
}

float** matrix_norm(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        printf("Error: Invalid input\n");
        return NULL;
    }
    
    float min = FLT_MAX, max = -FLT_MAX;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (matrix[i][j] < min) min = matrix[i][j];
            if (matrix[i][j] > max) max = matrix[i][j];
        }
    }
    
    if (max == min) {
        printf("Error: Cannot normalize matrix with same values\n");
        return NULL;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = (matrix[i][j] - min) / (max - min);
        }
    }
    
    return matrix;
}

float* matrix_getrow(float** matrix, int n, int m, int row) {
    if (matrix == NULL || n <= 0 || m <= 0 || row < 0 || row >= n) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* row_vector = (float*)malloc(m * sizeof(float));
    if (!row_vector) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    for (int j = 0; j < m; j++) {
        row_vector[j] = matrix[row][j];
    }
    
    return row_vector;
}

float* matrix_getcol(float** matrix, int n, int m, int col) {
    if (matrix == NULL || n <= 0 || m <= 0 || col < 0 || col >= m) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* col_vector = (float*)malloc(n * sizeof(float));
    if (!col_vector) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    for (int j = 0; j < n; j++) {
        col_vector[j] = matrix[j][col];
    }
    
    return col_vector;
}

float** matrix_setrow(float** matrix, int n, int m, int row, float* mas, size_t size_mas) {
    if (matrix == NULL || mas == NULL || n <= 0 || m <= 0 || row < 0 || row >= n || size_mas <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    for (int j = 0; j < m; j++) {
        if (j < size_mas) {
            matrix[row][j] = mas[j]; 
        } else {
            matrix[row][j] = 0.0; 
        }
    }

    return matrix;
}

float** matrix_setcol(float** matrix, int n, int m, int col, float* mas, size_t size_mas) {
    if (matrix == NULL || mas == NULL || n <= 0 || m <= 0 || col < 0 || col >= m || size_mas <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    for (int j = 0; j < n; j++) {
        if (j < size_mas) {
            matrix[j][col] = mas[j];  
        } else {
            matrix[j][col] = 0.0;  
        }
    }

    return matrix;  
}

float matrix_min(float** matrix, int n, int m) {
    if (matrix == NULL || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return FLT_MAX;  
    }

    float min_value = matrix[0][0]; 

    for (int i = 0; i < n; i++) {
        if (matrix[i] == NULL) {  
            fprintf(stderr, "Error: Invalid line #%d\n", i);
            return FLT_MAX;
        }
        for (int j = 0; j < m; j++) {
            if (matrix[i][j] < min_value) {
                min_value = matrix[i][j];
            }
        }
    }

    return min_value;
}

float matrix_max(float** matrix, int n, int m) {
    if (matrix == NULL || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return FLT_MIN;  
    }

    float max_value = matrix[0][0]; 

    for (int i = 0; i < n; i++) {
        if (matrix[i] == NULL) { 
            fprintf(stderr, "Error: Invalid line #%d\n", i);
            return FLT_MIN;
        }
        for (int j = 0; j < m; j++) {
            if (matrix[i][j] > max_value) {
                max_value = matrix[i][j];
            }
        }
    }

    return max_value;
}

float* matrix_colmean(float** matrix, int n, int m) {
    if (matrix == NULL || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float* col_means = (float*)malloc(m * sizeof(float));
    if (!col_means) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    for (int j = 0; j < m; j++) {
        float sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += matrix[i][j];
        }
        col_means[j] = sum / n;  
    }

    return col_means;
}

float* matrix_rowmean(float** matrix, int n, int m) {
    if (matrix == NULL || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }

    float* row_means = (float*)malloc(n * sizeof(float));
    if (!row_means) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        float sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += matrix[i][j];
        }
        row_means[i] = sum / m;  
    }

    return row_means;
}

bool matrix_issym(float** matrix, int n) {
    if (matrix == NULL || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return false;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_isiden(float** matrix, int n) {
    if (matrix == NULL || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return false;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i == j && matrix[i][j] != 1.0) || (i != j && matrix[i][j] != 0.0)) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_isdiag(float** matrix, int n) {
    if (matrix == NULL || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return false;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && matrix[i][j] != 0.0) {
                return false;
            }
        }
    }
    return true;
}

float** matrix_sortasc(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc(n * m * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    int index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            temp[index++] = matrix[i][j];
    
    mergeSortAsc(temp, 0, n * m - 1);
    
    index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j] = temp[index++];
    
    free(temp);
    return matrix;
}

float** matrix_sortdesc(float** matrix, int n, int m) {
    if (!matrix || n <= 0 || m <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc(n * m * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    int index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            temp[index++] = matrix[i][j];
    
    mergeSortDesc(temp, 0, n * m - 1);
    
    index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j] = temp[index++];
    
    free(temp);
    return matrix;
}

float** matrix_sortcolasc(float** matrix, int n, int m, int col) {
    if (!matrix || n <= 0 || m <= 0 || col < 0 || col >= m) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc(n * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    for (int i = 0; i < n; i++)
        temp[i] = matrix[i][col];
    
    mergeSortAsc(temp, 0, n - 1);
    
    for (int i = 0; i < n; i++)
        matrix[i][col] = temp[i];
    
    free(temp);
    return matrix;
}

float** matrix_sortcoldesc(float** matrix, int n, int m, int col) {
    if (!matrix || n <= 0 || m <= 0 || col < 0 || col >= m) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc(n * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    for (int i = 0; i < n; i++)
        temp[i] = matrix[i][col];
    
    mergeSortDesc(temp, 0, n - 1);
    
    for (int i = 0; i < n; i++)
        matrix[i][col] = temp[i];
    
    free(temp);
    return matrix;
}

float** matrix_sortrowasc(float** matrix, int n, int m, int row) {
    if (!matrix || n <= 0 || m <= 0 || row < 0 || row >= n) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    mergeSortAsc(matrix[row], 0, m - 1);
    return matrix;
}

float** matrix_sortrowdesc(float** matrix, int n, int m, int row) {
    if (!matrix || n <= 0 || m <= 0 || row < 0 || row >= n) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    mergeSortDesc(matrix[row], 0, m - 1);
    return matrix;
}

float** matrix_sortrangeasc(float** matrix, int n, int m, int start, int end) {
    if (!matrix || n <= 0 || m <= 0 || start < 0 || end >= n * m || start > end) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc((end - start + 1) * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    int index = 0;
    for (int i = start; i <= end; i++)
        temp[index++] = matrix[i / m][i % m];
    
    mergeSortAsc(temp, 0, end - start);
    
    index = 0;
    for (int i = start; i <= end; i++)
        matrix[i / m][i % m] = temp[index++];
    
    free(temp);
    return matrix;
}

float** matrix_sortrangedesc(float** matrix, int n, int m, int start, int end) {
    if (!matrix || n <= 0 || m <= 0 || start < 0 || end >= n * m || start > end) {
        fprintf(stderr, "Error: Invalid input\n");
        return NULL;
    }
    
    float* temp = (float*)malloc((end - start + 1) * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }
    
    int index = 0;
    for (int i = start; i <= end; i++)
        temp[index++] = matrix[i / m][i % m];
    
    mergeSortDesc(temp, 0, end - start);
    
    index = 0;
    for (int i = start; i <= end; i++)
        matrix[i / m][i % m] = temp[index++];
    
    free(temp);
    return matrix;
}