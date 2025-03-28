#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void mergeAsc(float arr[], int left, int mid, int right);
void mergeSortAsc(float arr[], int left, int right);
void mergeDesc(float arr[], int left, int mid, int right);
void mergeSortDesc(float arr[], int left, int right);

float** matrix_norm(float** matrix, int n, int m);
float* matrix_getrow(float** matrix, int n, int m, int row);
float* matrix_getcol(float** matrix, int n, int m, int col);
float** matrix_setrow(float** matrix, int n, int m, int row, float* mas, size_t size_mas);
float** matrix_setcol(float** matrix, int n, int m, int col, float* mas, size_t size_mas);
float matrix_min(float** matrix, int n, int m);
float matrix_max(float** matrix, int n, int m);
float* matrix_colmean(float** matrix, int n, int m);
float* matrix_rowmean(float** matrix, int n, int m);
bool matrix_issym(float** matrix, int n);
bool matrix_isiden(float** matrix, int n);
bool matrix_isdiag(float** matrix, int n);
float** matrix_sortasc(float** matrix, int n, int m);
float** matrix_sortdesc(float** matrix, int n, int m);
float** matrix_sortcolasc(float** matrix, int n, int m, int col);
float** matrix_sortcoldesc(float** matrix, int n, int m, int col);
float** matrix_sortrowasc(float** matrix, int n, int m, int row);
float** matrix_sortrowdesc(float** matrix, int n, int m, int row);
float** matrix_sortrangeasc(float** matrix, int n, int m, int start, int end);
float** matrix_sortrangedesc(float** matrix, int n, int m, int start, int end);
float** matrix_input(int n, int m);
float** matrix_createv(int n, int m);
float** matrix_rcreate(int n, int m, int min, int max);
float** matrix_inputsqr(int n);
float** matrix_createsqrv(int n);
float** matrix_rcreatesqr(int n, int min, int max);
float** matrix_copy(float** matrix, int n, int m);
void matrix_free(float** matrix);
void matrix_print(float** matrix, int n, int m);
int matrix_savef(float** matrix, int n, int m, char* filename);
float** matrix_loadf(char* filename);
float** matrix_sum(float** A, float** B, int n, int m);
float** matrix_sub(float** A, float** B, int n, int m);
float** matrix_multp(float** A, float** B, int n, int p, int m);
float** matrix_multpscl(float** matrix, int n, int m, float scalar);
float** matrix_power(float** matrix, int n, int p);
float** matrix_transpose(float** matrix, int n, int m);
float matrix_determin(float** matrix, int n);
float** matrix_inver(float** matrix, int n);
float matrix_trace(float** matrix, int n);

#endif