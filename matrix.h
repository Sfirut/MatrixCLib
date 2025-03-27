#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

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