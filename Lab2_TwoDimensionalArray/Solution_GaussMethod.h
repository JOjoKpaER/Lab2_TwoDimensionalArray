#pragma once

#include <Matrix.h>

#define EPS 0.00001f

void Solution_GaussMethod(Matrix *M, Matrix *Y, Matrix *X);
float Det_GaussMethod(Matrix *M);
void Inv_GaussMethod(Matrix *M, Matrix *E);
float Delta(Matrix *M, Matrix *Y, Matrix *X);