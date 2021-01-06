#include "Matrix.h"

double eps = 1E-10;;
vector<double> x0;
int maxiter = 10000;


int main()
{

	Mat<double> A;
	A.Mat_init();

	//vector<vector<double>> MxGb(A.n);
	//A.MatGilbert(A.n, MxGb, A.ig, A.ggl, A.di);

	x0.resize(A.n);
	vector<double> x1(A.n);
	ifstream InitVal("x0.txt");

	for (int i = 0; i < A.n; i++)
		InitVal >> x0[i];

	//A.PowSolv(x0, maxiter, eps);
	A.InversePowSolv(x0, maxiter, eps);
	return 0;
}
