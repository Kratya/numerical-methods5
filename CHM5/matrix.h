#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

template<class mytype>
class Mat
{
public:
	int n;
	vector<int> ig;
	vector<mytype> ggl, ggu, di;

	void Mat_init();
	void Gilb_init();
	void MultMatXVec(vector<mytype>& x, vector<mytype>& res);
	mytype PowSolv(vector<mytype>& x0, int maxiter, int eps);
	mytype InversePowSolv(vector<mytype>& x0, int maxiter, int eps);
	void MatGilbert(int k, vector<vector<mytype>>& mxg, vector<int>& ia, vector<mytype>& al, vector<mytype>& di);//храним только нижний треугольнник

private:

	mytype norma(vector<mytype>& b);
	void LU();
	void LUSolver(vector<mytype>& x, vector<mytype>& b);
};

template <typename mytype>
void Mat<mytype>::Gilb_init()
{
	ifstream F_di("di.txt");
	F_di >> n;
	F_di.close();
}

template <typename mytype>
void Mat<mytype>::Mat_init()
{
	ifstream F_di("di.txt");
	ifstream F_ggl("ggl.txt");
	ifstream F_ggu("ggu.txt");
	ifstream F_ig("ig.txt");

	F_di >> n;
	di.resize(n);

	for (int i = 0; i < n; i++)	
		F_di >> di[i];

	ig.resize(n + 1);

	for (int i = 0; i < n + 1; i++)	
		F_ig >> ig[i];

	ggl.resize(ig[n]);
	ggu.resize(ig[n]);
	

	for (int i = 0; i < ig[n]; i++)
	{
		F_ggl >> ggl[i];
		F_ggu >> ggu[i];
	}
	

	if (ig[0])	for (int i = 0, nn = n + 1; i < nn; i++) ig[i]--;

	F_ggl.close();
	F_ggu.close();
	F_ig.close();
	F_di.close();
}

template <typename mytype>
mytype Mat<mytype>::norma(vector<mytype>& b)
{
	mytype nrm = 0;
	for (int i = 0; i < n; i++)	
		nrm += b[i] * b[i];

	return sqrt(nrm);
}

template <typename mytype>
void  Mat<mytype>::MultMatXVec(vector<mytype>& x, vector<mytype>& res)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = di[i] * x[i];
		int i0 = ig[i];
		int i1 = ig[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			res[i] += ggl[k] * x[j];
			res[j] += ggu[k] * x[i];
		}
	}
}

template <typename mytype>
mytype Mat<mytype>::PowSolv(vector<mytype>& x0, int maxiter, int eps)
{
	vector<mytype> xk(n);
	
	mytype normaxk = 0;
	mytype normax0 = norma(x0);
	mytype lambda0 = 0;
	mytype lambdak = 0;
	mytype nevas = 1;

	for (int k = 0; k < maxiter && nevas>eps; k++)
	{
		MultMatXVec(x0, xk);
		normaxk = norma(xk);
		lambdak = normaxk / normax0;
		nevas = abs((lambdak - lambda0) / lambdak);
		lambda0 = lambdak;
		x0 = xk;
		normax0 = normaxk;
		if (k % 1 == 0)
		{
			for (int i = 0; i < n; i++)
				x0[i] /= normax0;
			normax0 = 1;
		}
		cout << "\rN_itter = " << k + 1 << " l_max = " << scientific << lambdak << endl;
	}
	return lambdak;
}

template <typename mytype>
void Mat<mytype>::LU() // разложение lu
{
	for (int i = 0; i < n; i++)
	{
		int j = i - (ig[i + 1] - ig[i]);// первый не нулевой элемент i строки
		mytype sumdi = 0;
		for (int k = ig[i]; k < ig[i + 1]; k++, j++)
		{
			int iaj = ig[j];
			mytype sumal = 0;
			mytype sumau = 0;
			int tl = ig[i];
			int tu = ig[j];
			int a = k - ig[i] - (ig[j + 1] - ig[j]);//сдиг
			if (a < 0)	tu += abs(a);
			else tl += a;
			for (; tl < k; tl++, tu++)
			{
				sumal += ggu[tu] * ggl[tl];
				sumau += ggu[tl] * ggl[tu];
			}
			ggl[k] = (ggl[k] - sumal) / di[j];
			ggu[k] = (ggu[k] - sumau) / di[j];
			sumdi += ggl[k] * ggu[k];
		}
		di[i] = sqrt(di[i] - sumdi);
	}
}

template <typename mytype>
void Mat<mytype>::LUSolver(vector<mytype>& b, vector<mytype>& x)
{
	vector<mytype> y(n);
	for (int i = 0; i < n; i++)
	{
		mytype sum = 0;
		int j = i - (ig[i + 1] - ig[i]);
		for (int p = ig[i]; p < ig[i + 1]; p++, j++)
		{
			sum += ggl[p] * y[j];
		}
		y[i] = (b[i] - sum) / di[i];
	}
	for (int i = n - 1; i >= 0; i--)
	{
		int iai = ig[i];
		int iai1 = ig[i + 1];
		int j = i - (iai1 - iai);
		x[i] = y[i] / di[i];
		for (int p = iai; p < iai1; p++, j++)
		{
			y[j] -= ggu[p] * x[i];
		}
	}
}

template <typename mytype>
mytype Mat<mytype>::InversePowSolv(vector<mytype>& x0, int maxiter, int eps)
{
	vector<mytype> xk(n);
	mytype normaxk = 0;
	mytype normax0 = norma(x0);
	mytype lambda0 = 0;
	mytype lambdak = 0;
	mytype nevas = 1;

	LU();
	for (int k = 0; k < maxiter && nevas>eps; k++)
	{
		LUSolver(x0, xk);
		normaxk = norma(xk);

		lambdak = normaxk / normax0;
		nevas = abs((lambdak - lambda0) / lambdak);
		x0 = xk;
		lambda0 = lambdak;
		normax0 = normaxk;
		if (k % 10 == 0)
		{
			for (int i = 0; i < n; i++)
				x0[i] /= normax0;
			normax0 = 1;
		}
		cout << "\rN_Itter = " << k + 1 << ", l_min: " << scientific << 1 / lambdak << endl;
	}

	return 1 / lambdak;
}

template <typename mytype>
void Mat<mytype>::MatGilbert(int k, vector<vector<mytype>>& mxg, vector<int>& ia, vector<mytype>& al, vector<mytype>& di)//храним только нижний треугольнник
{
	int col = 0;
	int iaa = 1;

	di.resize(n);
	ggl.resize((n * n - n) / 2);
	ggu.resize((n * n - n) / 2);
	di.resize(n);
	ig.resize(n + 1);
	

	for (int i = 0; i < k; i++)
	{
		mxg[i].resize(i + 1);
		for (int j = 0; j <= i; j++)
		{
			mytype c = mytype(1) / (i + j + 1);
			mxg[i][j] = c;
		}

	}
	di[0] = mxg[0][0];
	ia[0] = 0;
	
	for (int i = 0; i < k; i++)
		for (int kn = 0; kn < i; kn++)
		{
			if (mxg[i][kn] != 0)
			{
				ia[iaa] = col;
				iaa++;
				for (int j = kn; j < i; j++)
				{
					al[col] = mxg[i][j];
					col++;
				}
				di[i] = mxg[i][i];
				kn = i;
			}
		}
	ia[iaa] = col;
	ggu = ggl;
}
