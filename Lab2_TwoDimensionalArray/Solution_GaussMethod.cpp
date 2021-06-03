#include <cstdlib>
#include <stdio.h>

#include <Matrix.h>
#include <Solution_GaussMethod.h>

void Solution_GaussMethod(Matrix *M, Matrix *Y, Matrix *X) {
	if (M->h < 1 || M->w < 1 || M->h != M->w || Y->h != M->h || Y->w != 1 || X->h != M->h || X->w != 1) {
		printf("Incomplete matrix, cannot solve\n");
		return;
	}
	int k = 0;
	int index;
	int n = M->h;
	float max;
	while (k < n) {
		max = abs(M->get(k, k));
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(M->get(i, k) > max)) {
				max = abs(M->get(i, k));
				index = i;
			}
		}
		if (max < EPS) {
			printf("Cannot solve due to null column\n");
			for (int i = 0; i < X->h; i++) {
				X->set(i, 0, NAN);
			}
			return;
		}
		M->swapRows(k, index);
		Y->swapRows(k, index);
		for (int i = k; i < n; i++) {
			float t = M->get(i, k);
			if (abs(t) < EPS) continue;
			for (int j = 0; j < n; j++) {
				M->set(i, j, M->get(i, j) / t);
			}
			Y->set(i, 0, Y->get(i, 0) / t);
			if (i == k) continue;
			for (int j = 0; j < n; j++) {
				M->set(i, j, M->get(i, j) - M->get(k, j));
			}
			Y->set(i, 0, Y->get(i, 0) - Y->get(k, 0));
		}
		k++;
	}
	/*=======================================================*/
	for (k = n - 1; k >= 0; k--) {
		X->set(k, 0, Y->get(k, 0));
		for (int i = 0; i < k; i++) {
			Y->set(i, 0, Y->get(i, 0) - M->get(i, k) * X->get(k, 0));
		}
	}
}

float Det_GaussMethod(Matrix *M) {
	if (M->h < 1 || M->w < 1 || M->h != M->w) {
		printf("Incomplete matrix, cannot solve\n");
		return NAN;
	}
	float ret = 1;
	int n = M->h;
	for (int i = 0; i < n; i++){
		int k = i;
		for (int j = i + 1; j < n; j++) {
			if (abs(M->get(j, i)) > abs(M->get(k, j))) {
				k = j;
			}
		}
		if (abs(M->get(k, i)) < EPS) {
			return 0.0f;
		}
		M->swapRows(i, k);
		if (i != k) {
			ret *= -1;
		}
		ret *= M->get(i, i);
		for (int j = i + 1; j < n; j++) {
			M->set(i, j, M->get(i, j) / M->get(i, i));
		}
		/*=================================================*/
		for (int j = 0; j < n; j++) {
			if (j != i && abs(M->get(i, j)) > EPS) {
				for (int k = i + 1; k < n; k++) {
					M->set(i, j, M->get(j, k) - 1.0f * M->get(i, k) * M->get(j, i));
				}
			}
		}
	}
	return ret;
}

void Inv_GaussMethod(Matrix *M, Matrix *E) {
	if (M->h < 1 || M->w < 1 || M->h != M->w) {
		printf("Incomplete matrix, cannot solve\n");
		return;
	}
	int k = 0;
	int index;
	int n = M->h;
	float max;
	while (k < n) {
		max = abs(M->get(k, k));
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(M->get(i, k) > max)) {
				max = abs(M->get(i, k));
				index = i;
			}
		}
		if (max < EPS) {
			printf("Cannot solve due to null column\n");
			return;
		}
		M->swapRows(k, index);
		E->swapRows(k, index);
		for (int i = k; i < n; i++) {
			float t = M->get(i, k);
			if (abs(t) < EPS) continue;
			for (int j = 0; j < n; j++) {
				M->set(i, j, M->get(i, j) / t);
				E->set(i, j, E->get(i, j) / t);
			}
			if (i == k) continue;
			for (int j = 0; j < n; j++) {
				M->set(i, j, M->get(i, j) - M->get(k, j));
				E->set(i, j, E->get(i, j) - E->get(k, j));
			}
		}
		k++;
	}
}
float Delta(Matrix *M, Matrix *Y, Matrix *X) {
	if (M->h < 1 || M->w < 1 || M->h != M->w || Y->h != M->h || Y->w != 1 || X->h != M->h || X->w != 1) {
		printf("Incomplete matrix, cannot solve\n");
		return NAN;
	}
	Matrix *B = M->Multiply(X);
	float max = EPS;
	for (int i = 0; i < B->h; i++) {
		if (abs(B->get(i, 0) - Y->get(i, 0)) > max) {
			max = abs(B->get(i, 0) - Y->get(i, 0));
		}
	}
	if (max <= EPS) max = 0.0f;
	return max;
}