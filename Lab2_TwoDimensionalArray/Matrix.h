#pragma once

#include <stdio.h>

enum Type{
	Identity
};

class Matrix
{
private:

	float** array;
		
public:
	int w;
	int h;

	Matrix() {
		array = new (float*[1]);
		array[0] = new float[0];
		w = 1;
		h = 1;
	}

	Matrix(int MaxH, int MaxW) {
		if (MaxH < 1 || MaxW < 1)
		{
			printf("Incorrect size values\n");
			return;
		}
		array = new (float*[MaxH]);
		for (int i = 0; i < MaxH; i++) {
			array[i] = new float[MaxW]();
		}
		h = MaxH;
		w = MaxW;
	}

	Matrix(const Matrix& MatrixPtr) {
		h = MatrixPtr.h;
		w = MatrixPtr.w;

		array = new (float*[h]);
		for (int i = 0; i < h; i++) {
			array[i] = new float[w]();
		}
		for (int i = 0; i < h; i++)
			for (int k = 0; k < w; k++)
				array[i][k] = MatrixPtr.array[i][k];
	}

	Matrix(int MaxH, int MaxW, Type t) {
		if (MaxH < 1 || MaxW < 1)
		{
			printf("Incorrect size values\n");
			return;
		}
		switch(t){
		case(Identity): {
			array = new (float*[MaxH]);
			for (int i = 0; i < MaxH; i++) {
				array[i] = new float[MaxW]();
			}
			h = MaxH;
			w = MaxW;
			for (int i = 0; i < h; i++) {
				this->set(i, i, 1);
			}
			break;
			}
		default:
		{
			printf("Cannot resolve matrix for type %i\n", t);
			return;
		}
		}
	}

	~Matrix() {
		if (w > 0)
		for (int i = 0; i < h; i++) {
			if (array[i])
			delete[] array[i];
		}
		if (h > 0)
		delete[] array;
	}

	float get(int i, int j) {
		if (i < 0 || j < 0 || i > h || j > w) {
			printf("Incorrect position values\n");
			return NULL;
		}
		else
			return array[i][j];
	}

	void set(int i, int j, float value) {
		if (i < 0 || j < 0 || i > h || j > w) {
			printf("Incorrect position values: i %i j %i\n", i, j);
			return;
		}
		else
			array[i][j] = value;
	}

	void resize(int MaxH, int MaxW) {
		if (MaxH < 1 || MaxW < 1)
		{
			printf("Incorrect size values\n");
			return;
		}
		Matrix *tempM = new Matrix(*this);
		delete this;
		array = new (float*[MaxH]);
		for (int i = 0; i < MaxH; i++) {
			array[i] = new float[MaxW]();
		}
		h = MaxH;
		w = MaxW;
		for (int i = 0; i < h; i++)
			for (int j = 0; j < w; j++)
			{
				array[i][j] = tempM->array[i][j];
			}
		delete tempM;
	}

	void swapRows(int a, int b) {
		if (a < 0 || b < 0 || a >= h || b >= h) {
			printf("Incorrect indexes while swap rows procedure\n");
			exit(0);
		}
		float* temp = array[a];
		array[a] = array[b];
		array[b] = temp;
		temp = nullptr;
	}

	void print() {
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++) {
				printf("%f ", array[i][j]);
			}
			printf("\n");
		}
	}

	Matrix *Multiply(Matrix *R) {
		if (this->w != R->h) {
			printf("Cannot multyply these matrices\n");
			return nullptr;
		}
		Matrix *ret = new Matrix(this->h, R->w);
		for (int i = 0; i < ret->h; i++) {
			for (int j = 0; j < ret->w; j++) {
				float t = 0.0f;
				for (int k = 0; k < this->w; k++) {
					t += this->get(i, k) * R->get(k, j);
				}
				ret->set(i, j, t);
			}
		}
		return ret;
	}

};

