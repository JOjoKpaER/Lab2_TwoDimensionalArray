#include <iostream>
#include <Matrix.h>
#include <Solution_GaussMethod.h>

#define MAXRAND 1000

static Matrix Precise = Matrix(3, 3);
static Matrix PreciseB = Matrix(3, 1);

static Matrix Random = Matrix(3, 3);
static Matrix RandomB = Matrix(3, 1);

static Matrix Iden = Matrix(3,3);
static Matrix IdenB = Matrix(3,1);

static Matrix Collinear = Matrix(3, 3);
static Matrix CollinearB = Matrix(3, 1);

static Matrix Gilbert = Matrix(3, 3);
static Matrix GilbertB = Matrix(3, 1);

int main()
{
	Precise.set(0, 0, 2); Precise.set(0, 1, 4); Precise.set(0, 2, 1);
	Precise.set(1, 0, 5); Precise.set(1, 1, 2); Precise.set(1, 2, 1);
	Precise.set(2, 0, 2); Precise.set(2, 1, 3); Precise.set(2, 2, 4);
	PreciseB.set(0, 0, 36); PreciseB.set(1, 0, 47); PreciseB.set(2, 0, 37);

	Random.set(0, 0, rand() % MAXRAND); Random.set(0, 1, rand() % MAXRAND); Random.set(0, 2, rand() % MAXRAND);
	Random.set(1, 0, rand() % MAXRAND); Random.set(1, 1, rand() % MAXRAND); Random.set(1, 2, rand() % MAXRAND);
	Random.set(2, 0, rand() % MAXRAND); Random.set(2, 1, rand() % MAXRAND); Random.set(2, 2, rand() % MAXRAND);
	RandomB.set(0, 0, rand() % MAXRAND); RandomB.set(1, 0, rand() % MAXRAND); RandomB.set(2, 0, rand() % MAXRAND);

	Iden.set(0, 0, 1); Iden.set(0, 1, 0); Iden.set(0, 2, 0);
	Iden.set(1, 0, 0); Iden.set(1, 1, 1); Iden.set(1, 2, 0);
	Iden.set(2, 0, 0); Iden.set(2, 1, 0); Iden.set(2, 2, 1);
	IdenB.set(0, 0, 1); IdenB.set(1, 0, 2); IdenB.set(2, 0, 3);

	Collinear.set(0, 0, 2); Collinear.set(0, 1, 4); Collinear.set(0, 2, 1);
	Collinear.set(1, 0, 5); Collinear.set(1, 1, 2); Collinear.set(1, 2, 1);
	Collinear.set(2, 0, -2); Collinear.set(2, 1, -4); Collinear.set(2, 2, -1);
	CollinearB.set(0, 0, 36); CollinearB.set(1, 0, 47); CollinearB.set(2, 0, 37);

	Gilbert.set(0, 0, 1.0f); Gilbert.set(0, 1, 0.5f); Gilbert.set(0, 2, 1.0f/3.0f);
	Gilbert.set(1, 0, 0.5f); Gilbert.set(1, 1, 1.0f/3.0f); Gilbert.set(1, 2, 0.25f);
	Gilbert.set(2, 0, 1.0f/3.0f); Gilbert.set(2, 1, 0.25f); Gilbert.set(2, 2, 0.2f);
	GilbertB.set(0, 0, 36); GilbertB.set(1, 0, 47); GilbertB.set(2, 0, 37);


	printf("Precise marix:\n");
	Matrix *M = new Matrix(Precise);
	Matrix *Y = new Matrix(PreciseB);
	Matrix *X = new Matrix(3, 1);
	printf("Matrix:\n");
	M->print();
	Solution_GaussMethod(M, Y, X);
	printf("Solution:\n");
	X->print();
	printf("Delta: %f\n+", Delta(new Matrix(Precise), new Matrix(PreciseB), X));

	delete M;
	delete Y;
	delete X;
	M = new Matrix(Precise);
	printf("Matrix:\n");
	M->print();
	float det = Det_GaussMethod(M);
	printf("Determinant: %f\n", det);

	delete M;

	Matrix *Inverse = new Matrix(3, 3, Identity);
	M = new Matrix(Precise);
	printf("Matrix:\n");
	M->print();
	Inv_GaussMethod(M, Inverse);
	printf("Inverse matrix:\n");
	Inverse->print();

	delete M;
	delete Inverse;
	/*=========================================================*/
	printf("Random matrix:\n");
	M = new Matrix(Random);
	Y = new Matrix(RandomB);
	X = new Matrix(3, 1);
	printf("Matrix:\n");
	M->print();
	Solution_GaussMethod(M, Y, X);
	printf("Solution:\n");
	X->print();
	printf("Delta: %f\n+", Delta(new Matrix(Random), new Matrix(RandomB), X));

	delete M;
	delete Y;
	delete X;
	M = new Matrix(Random);
	printf("Matrix:\n");
	M->print();
	det = Det_GaussMethod(M);
	printf("Determinant: %f\n", det);

	delete M;

	Inverse = new Matrix(3, 3, Identity);
	M = new Matrix(Random);
	printf("Matrix:\n");
	M->print();
	Inv_GaussMethod(M, Inverse);
	printf("Inverse matrix:\n");
	Inverse->print();

	delete M;
	delete Inverse;
	/*=========================================================*/
	printf("Identity matrix:\n");
	M = new Matrix(Iden);
	Y = new Matrix(IdenB);
	X = new Matrix(3, 1);
	printf("Matrix:\n");
	M->print();
	Solution_GaussMethod(M, Y, X);
	printf("Solution:\n");
	X->print();
	printf("Delta: %f\n+", Delta(new Matrix(Iden), new Matrix(IdenB), X));

	delete M;
	delete Y;
	delete X;
	M = new Matrix(Iden);
	printf("Matrix:\n");
	M->print();
	det = Det_GaussMethod(M);
	printf("Determinant: %f\n", det);

	delete M;

	Inverse = new Matrix(3, 3, Identity);
	M = new Matrix(Iden);
	printf("Matrix:\n");
	M->print();
	Inv_GaussMethod(M, Inverse);
	printf("Inverse matrix:\n");
	Inverse->print();

	delete M;
	delete Inverse;
	/*=========================================================*/
	printf("Collinear matrix:\n");
	M = new Matrix(Collinear);
	Y = new Matrix(CollinearB);
	X = new Matrix(3, 1);
	printf("Matrix:\n");
	M->print();
	Solution_GaussMethod(M, Y, X);
	printf("Solution:\n");
	X->print();
	printf("Delta: %f\n+", Delta(new Matrix(Collinear), new Matrix(CollinearB), X));

	delete M;
	delete Y;
	delete X;
	M = new Matrix(Collinear);
	printf("Matrix:\n");
	M->print();
	det = Det_GaussMethod(M);
	printf("Determinant: %f\n", det);

	delete M;

	Inverse = new Matrix(3, 3, Identity);
	M = new Matrix(Collinear);
	printf("Matrix:\n");
	M->print();
	Inv_GaussMethod(M, Inverse);
	printf("Inverse matrix:\n");
	Inverse->print();

	delete M;
	delete Inverse;
	/*=========================================================*/
	printf("Gilbert matrix:\n");
	M = new Matrix(Gilbert);
	Y = new Matrix(GilbertB);
	X = new Matrix(3, 1);
	printf("Matrix:\n");
	M->print();
	Solution_GaussMethod(M, Y, X);
	printf("Solution:\n");
	X->print();
	printf("Delta: %f\n+", Delta(new Matrix(Gilbert), new Matrix(GilbertB), X));

	delete M;
	delete Y;
	delete X;
	M = new Matrix(Gilbert);
	printf("Matrix:\n");
	M->print();
	det = Det_GaussMethod(M);
	printf("Determinant: %f\n", det);

	delete M;

	Inverse = new Matrix(3, 3, Identity);
	M = new Matrix(Gilbert);
	printf("Matrix:\n");
	M->print();
	Inv_GaussMethod(M, Inverse);
	printf("Inverse matrix:\n");
	Inverse->print();

	delete M;
	delete Inverse;
	/*=========================================================*/
	return 0;
}

