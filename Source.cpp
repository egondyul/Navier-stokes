#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <math.h>
#include<sstream>
#include <chrono>
using namespace std;
typedef std::vector<std::vector<double>> Matrix;

void read_file(Matrix& space)
{
	char name[1000] = { "test1.txt" };
	std::ifstream in;
	in.open(name);
	if (!in.is_open())
	{
		std::cout << "Sorry :(";
	}
	else
	{
		std::string str = {};

		while (std::getline(in, str))
		{
			double tmp = 0.0;
			std::stringstream is(str);
			std::vector<double> subspace;
			while (is >> tmp)
			{
				subspace.push_back(tmp);
			}
			space.push_back(subspace);
		}
	}

}

void poisson_1D(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0 && i == 0)
			{
				out[i][j] = (in[i][j + 1] - 2 * in[i][j] + p1) / (hx*hx);
			}
			else if (space[i][j] == 1 && i == 0 && j != 0 && j != (m - 1))
			{
				out[i][j] = (in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hx*hx);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				out[i][j] = (p2 - 2 * in[i][j] + in[i][j - 1]) / (hx*hx);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] = -out[i][j] + 2 / (hx*hy)*in[i][j];
		}
	}
}

void poisson_1D_withzero(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = in.size();
	int m = in[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && space[i + 1][j] == 0 && j == 0)
			{
				out[i][j]= (in[i][j + 1] - 2 * in[i][j] + p1) / (hx*hx);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && space[i + 1][j] == 0 && j != 0 && j != (m - 1))
			{
				out[i][j]= (in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hx*hx);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && space[i + 1][j] == 0 && j == (m - 1))
			{
				out[i][j] = (p2 - 2 * in[i][j] + in[i][j - 1]) / (hx*hx);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] = -out[i][j] + (2 / (hx*hy))*in[i][j];
		}
	}
}

void poisson(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = in.size();
	int m = in[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0 )
			{
				//first
				out[i][j] = (in[i + 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0)
			{
				//second
				out[i][j] = (in[i + 1][j] - in[i][j]) / hx + (p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0)
			{
				//third
				out[i][j] = (in[i - 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0)
			{
				//fourth
				out[i][j] = (in[i - 1][j] - in[i][j]) / hx + (p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0)
			{
				//dirichlet-right
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				//dirichlet-left
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//top-neiman
				out[i][j] = (in[i + 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//bottom-neiman
				out[i][j] = (in[i - 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i][j + 1] == 0 && i != 0 && i != (n - 1))
			{
				//left-neiman
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (in[i][j - 1] - in[i][j]) / hy;
			}
			else if (space[i][j] == 1 && space[i][j - 1] == 0 && i != 0 && i != (n - 1))
			{
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (in[i][j + 1] - in[i][j]) / hy;
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && j != 0 && j != (m - 1))
			{
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] = -out[i][j] + (2 / (hx*hy))*in[i][j];
		}
	}
}

void A_poisson(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = in.size();
	int m = in[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0)
			{
				//first
				out[i][j] = /*(in[i + 1][j] - in[i][j]) / (hx*hx) +*/ (in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0)
			{
				//second
				out[i][j] = /*(in[i + 1][j] - in[i][j]) / (hx*hx) +*/(p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0)
			{
				//third
				out[i][j] = /*(in[i - 1][j] - in[i][j]) / (hx*hx) +*/(in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0)
			{
				//fourth
				out[i][j] = /*(in[i - 1][j] - in[i][j]) / (hx*hx) +*/(p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0)
			{
				//dirichlet-right
				out[i][j] = /*(in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) +*/(in[i][j + 1] - 2 * in[i][j] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				//dirichlet-left
				out[i][j] = /*(in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) +*/(p2 - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//top-neiman
				out[i][j] = /*(in[i + 1][j] - in[i][j]) / (hx*hx) +*/(in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//bottom-neiman
				out[i][j] = /*(in[i - 1][j] - in[i][j]) / (hx*hx) +*/(in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i][j + 1] == 0 && i != 0 && i != (n - 1))
			{
				//left-neiman
				out[i][j] = /*(in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) +*/(in[i][j - 1] - in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i][j - 1] == 0 && i != 0 && i != (n - 1))
			{
				out[i][j] = /*(in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) +*/(in[i][j + 1] - in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && j != 0 && j != (m - 1))
			{
				out[i][j] = /*(in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) +*/(in[i][j + 1] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] = -out[i][j] + (2 / (hx*hy))*in[i][j];
		}
	}
}


void Resize(Matrix& A, int n, int m)
{
	A.resize(n);
	for (int i = 0; i < n; i++)
	{
		A[i].resize(m);
	}
}

double mlt(Matrix& a, Matrix& b)
{
	double tmp = 0.0;

	int n = a.size();
	int m = a[0].size();

	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			tmp += a[i][j] * b[i][j];
		}
	}
	return tmp;
}

double norm(Matrix& r)
{
	double tmp = 0.0;

	int n = r.size();
	int m = r[0].size();

	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			tmp += r[i][j] * r[i][j];
		}
	}
	return sqrt(tmp);
}


void min_nev(Matrix& space, Matrix& x, Matrix& rhs, int n, int m, double p11, double p22, double hx, double hy)
{
	double rescure = 1e5;
	int ii = 0;
	double error = 1e-2;
	double tau = 0.0;
	do {
		Matrix x1;
		Resize(x1, n, m);
		A_poisson(space, x, x1, p11, p22, hx, hy);
		Matrix r;
		Resize(r, n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = x1[i][j] - rhs[i][j];
			}
		}
		rescure = norm(r);
		Matrix rr1;
		Resize(rr1, n, m);
		A_poisson(space, r, rr1, p11, p22, hx, hy);
		tau = mlt(rr1, r) / (mlt(rr1, rr1));
		if (ii == 1)
		{
			std::cout << "tau = " << tau << endl;
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x1[i][j] = x[i][j] - tau * r[i][j];
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = x1[i][j];
			}
		}
		ii++;
	} while (rescure > error);
	std::cout << "min_nev" << endl;
	std::cout << "it = " << ii << endl;
	std::cout << "norm = " << rescure << endl;
	std::cout << "tau = " << tau << endl;
	std::cout << endl;
}

void gmres(Matrix& space, Matrix& x, Matrix& rhs, int n, int m, double p1, double p2, double hx, double hy)
{
	double rescure = 1e5;
	int ii = 0;
	double error = 1e-2;
	do {
		Matrix x1;
		Resize(x1, n, m);
		A_poisson(space, x, x1, p1, p2, hx, hy);
		Matrix r;
		Resize(r, n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = x1[i][j] - rhs[i][j];
			}
		}
		rescure = norm(r);
		Matrix rr1;
		Resize(rr1, n, m);
		A_poisson(space, r, rr1, p1, p2, hx, hy);
		double tau = mlt(r, r) / (mlt(rr1, r));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x1[i][j] = x[i][j] - tau * r[i][j];
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = x1[i][j];
			}
		}
		ii++;
	} while (rescure > error);
	std::cout << "gmres" << endl;
	std::cout << "it = " << ii << endl;
	std::cout << "norm = " << rescure << endl;
	std::cout << endl;
}

/* ================================================================= */
void poisson_vect(Matrix& space, std::vector<double>& in, std::vector<double>& out, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			int k = i * m + j;
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0 && i == 0)
			{
				out[k] = (in[k + 1] - 2 * in[k] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && i == 0 && j != 0 && j != (m - 1))
			{
				out[k] = (in[k + 1] - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				out[k] = (p2 - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0 && space[i - 1][j] == 0)
			{
				out[k] = (in[k + 1] - 2 * in[k] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && space[i - 1][j] == 0 && j != 0 && j != (m - 1))
			{
				out[k] = (in[k + 1] - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0 && space[i - 1][j] == 0)
			{
				out[k] = (p2 - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0)
			{
				//first
				out[k] = (in[(i + 1)*m+j] - in[k]) / hx + (in[k + 1] - 2 * in[k] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0)
			{
				//second
				out[k] = (in[(i + 1)*m + j] - in[k]) / hx + (p2 - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0)
			{
				//third
				out[k] = (in[(i - 1)*m + j] - in[k]) / hx + (in[k + 1] - 2 * in[k] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0)
			{
				//fourth
				out[k] = (in[(i - 1)*m + j] - in[k]) / hx + (p2 - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0)
			{
				//dirichlet-right
				out[k] = (in[(i + 1)*m + j] - 2 * in[k] + in[(i - 1)*m + j]) / (hx*hx) + (in[k + 1] - 2 * in[k] + p1) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				//dirichlet-left
				out[k] = (in[(i + 1)*m + j] - 2 * in[k] + in[(i - 1)*m + j]) / (hx*hx) + (p2 - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//top-neiman
				out[k] = (in[(i + 1)*m + j] - in[k]) / hx + (in[k + 1] - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j != 0 && j != (m - 1))
			{
				//bottom-neiman
				out[k] = (in[(i - 1)*m + j] - in[k]) / hx + (in[k + 1] - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i][j + 1] == 0 && i != 0 && i != (n - 1))
			{
				//left-neiman
				out[k] = (in[(i + 1)*m + j] - 2 * in[k] + in[(i - 1)*m + j]) / (hx*hx) + (in[k - 1] - in[k]) / hy;
			}
			else if (space[i][j] == 1 && space[i][j - 1] == 0 && i != 0 && i != (n - 1))
			{
				out[k] = (in[(i + 1)*m + j] - 2 * in[k] + in[(i - 1)*m + j]) / (hx*hx) + (in[k] - in[k]) / hy;
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && j != 0 && j != (m - 1))
			{
				out[k] = (in[(i + 1)*m + j] - 2 * in[k] + in[(i - 1)*m + j]) / (hx*hx) + (in[k + 1] - 2 * in[k] + in[k - 1]) / (hy*hy);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			int k = i * m + j;
			out[k] = -out[k] + 2 / (hx*hy)*in[k];
		}
	}
}

double norm_vect(std::vector<double>& a)
{
	double tmp=0.0;

	int n = a.size();
	for (int i = 0; i < n; i++)
	{
		tmp += a[i] * a[i];
	}
	return sqrt(tmp);
}

double mlt_vect(std::vector<double>& a, std::vector<double>& b)
{
	double tmp=0.0;

	int n = a.size();
	for (int i = 0; i < n; i++)
	{
		tmp += a[i] * b[i];
	}
	return tmp;
}

void min_vect(Matrix& space, std::vector<double>x, std::vector<double>rhs, double p1, double p2, double hx, double hy)
{
	int n = x.size();
	std::vector<double>x0;
	x0.resize(n);
	double error = 10 ^ 5;
	int ii = 0;
	double error_min = 1e-3;
	do
	{
		poisson_vect(space, x0, x, p1, p2, hx, hy);
		std::vector<double> r;
		r.resize(n);
		for (int i = 0; i < n; i++)
		{
			r[i] = x[i] - rhs[i];
		}
		error = norm_vect(r);
		std::vector<double> r1;
		r1.resize(n);
		poisson_vect(space, r, r1, p1, p2, hx, hy);
		double tau = mlt_vect(r1, r) / (mlt_vect(r1, r1));
		for (int i = 0; i < n; i++)
		{
			x[i] = x0[i] - tau * r[i];
		}
		for (int i = 0; i < n; i++)
		{
			x0[i] = x[i];
		}
		ii++;
	} while (error > error_min);
	std::cout << "it = " << ii << endl;
	std::cout << "norm = " << error << endl;
}
/*=======================================================================*/


void gmresnew(Matrix& space, Matrix& x0, Matrix& rhs, int n,int m, double p1, double p2, double hx, double hy)
{
	Matrix Ax; 
	Resize(Ax, n, m);

	A_poisson(space, x0, Ax, p1, p2, hx, hy);
	Matrix r;
	Resize(r, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}

	Matrix x;
	Resize(x, n, m);
	double error = norm(r);
	double eps = 1e-2;
	int it = 0;
	double tau = 0.01;
	while (error > eps) //&& it < endCycle)
	{
		it++;
		A_poisson(space, x0, Ax, p1, p2, hx, hy);

		Matrix Ar;
		Resize(Ar, n, m);
		A_poisson(space, r, Ar, p1, p2, hx, hy);
		
		if (mlt(Ar, r) != 0.0)
		{
			tau = mlt(r, r) / (mlt(Ar, r));
		}
		if (it == 1)
		{
			std::cout << "tau = " << tau << endl;
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = x0[i][j] - tau * (Ax[i][j] - rhs[i][j]);
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = rhs[i][j] - Ax[i][j];
			}
		}
		error = norm(r);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x0[i][j] = x[i][j];
			}
		}
	}
	std::cout << "gmres2" << endl;
	cout << "iterations = " << it << endl;
	std::cout << "error = " << error << std::endl;
	std::cout << "tau = " << tau << endl;
	std::cout << endl;
}

void minnev(Matrix& space, Matrix& x0, Matrix& rhs, int n, int m, double p1, double p2, double hx, double hy)
{
	Matrix Ax;
	Resize(Ax, n, m);

	poisson(space, x0, Ax, p1, p2, hx, hy);
	Matrix r;
	Resize(r, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}

	/*Matrix Ar;
	Resize(Ar, n, m);
	poisson(space, r, Ar, p1, p2, hx, hy);
	double tau = 0.01;
	if (mlt(Ar, Ar) != 0.0)
	{
		tau = mlt(Ar, r) / (mlt(Ar, Ar));
	}*/

	Matrix x;
	Resize(x, n, m);
	double error = norm(r);
	double eps = 1e-2;
	int it = 0;
	double tau = 0.01;
	while (error > eps) //&& it < endCycle)
	{
		it++;
		poisson(space, x0, Ax, p1, p2, hx, hy);
		Matrix Ar;
		Resize(Ar, n, m);
		poisson(space, r, Ar, p1, p2, hx, hy);
		
		if (mlt(Ar, Ar) != 0.0)
		{
			tau = mlt(Ar, r) / (mlt(Ar, Ar));
		}
		if (it == 1)
		{
			std::cout << "tau = " << tau << endl;
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (tau > 0.0)
				{
					x[i][j] = x0[i][j] - tau * (Ax[i][j] - rhs[i][j]);
				}
				else
				{
					x[i][j] = x0[i][j] + tau * (Ax[i][j] - rhs[i][j]);
				}
				
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = rhs[i][j] - Ax[i][j];
			}
		}
		error = norm(r);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x0[i][j] = x[i][j];
			}
		}
	}
	std::cout << "minnev2" << endl;
	std::cout << "iterations = " << it << endl;
	std::cout << "error = " << error << std::endl;
	std::cout << "tau = " << tau << endl;
	std::cout << endl;
}

/*void minres_saad(Matrix& space, Matrix& x, Matrix& rhs, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	Matrix Ax;
	Resize(Ax, n, m);
	poisson(space, x, Ax, p1, p2, hx, hy);

	Matrix r;
	Resize(r, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 1)
			{
				r[i][j] = rhs[i][j]- Ax[i][j];
			}
		}
	}
	double error;// = norm(r);

	Matrix Ar;
	Resize(Ar, n, m);
	poisson(space, r, Ar, p1, p2, hx, hy);

	double eps = 1e-2;
	double tau;
	int it = 0;
	do
	{
		it++;
		tau = mlt(Ar, r) / (mlt(Ar, Ar));
		//tau = abs(tau);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (space[i][j] == 1)
				{
					x[i][j] = x[i][j] - tau * r[i][j];
					r[i][j] = r[i][j] + tau * Ar[i][j];
				}
			}
		}
		poisson(space, r, Ar, p1, p2, hx, hy);
		error = norm(r);
	} while (error > eps);
	std::cout << "min_res_saad" << endl;
	std::cout << "it = " << it<< endl;
	std::cout << "norm of res = " << error << endl;
	std::cout << endl;
}*/

void Dpoisson(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	for (int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 1))
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j+1]) / (hx*hx);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx);
			}
		}
	}

	/*for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] += 2 * in[i][j] / (hx*hx);
		}
	}*/

}

void DDpoisson(Matrix& space, Matrix& in, Matrix& out, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			//top
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j == 0)
			{
				out[i][j]= (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j != 0 && j != (m-1))
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i - 1][j] == 0 && j == (m - 1))
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (in[i][j] - in[i + 1][j] )/ (hy*hy);
			}
			//inside
			else if (space[i][j] == 1 && j == 0 && space[i+1][j] == 1 &&space[i-1][j]==1)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j])/(hy*hy);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 1 && space[i - 1][j] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 1 && space[i - 1][j] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//bottom
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j == 0)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j != 0 && j != (m - 1))
			{
				out[i][j]= (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j == (m - 1))
			{
				out[i][j]= (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
		}
	}

	/*for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] += 2 * in[i][j] / (hx*hy);
		}
	}*/

}

void simple_iterations(Matrix& space, Matrix& x0, Matrix& rhs, double t0, double p0, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();
	Matrix Ax;
	Resize(Ax, n, m);

	Dpoisson(space, x0, Ax, p1, p2, hx, hy);
	Matrix r;
	Resize(r, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}

	Matrix x;
	Resize(x, n, m);
	double error = norm(r);
	double eps = 1e-3;
	//int endCycle = 1000;
	int it = 0;
	while (error > eps) //&& it < endCycle)
	{
		it++;
		Dpoisson(space, x0, Ax, p1, p2, hx, hy);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = x0[i][j] - t0 * (Ax[i][j] - rhs[i][j]);
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = rhs[i][j] - Ax[i][j];
			}
		}
		error = norm(r);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x0[i][j] = x[i][j];
			}
		}
	}
	std::cout << "simple_iteration" << endl;
	cout << "iterations = " << it << endl;
	std::cout << "error = " << error << std::endl;
	std::cout << endl;
}

void gmres_saad(Matrix& space, Matrix& x, Matrix& rhs,double p1, double p2, double hx, double hy )
{
	int n = space.size();
	int m = space[0].size();

	Matrix Ax;
	Resize(Ax, n, m);
	DDpoisson(space, x, Ax, p1, p2, hx, hy);

	Matrix r;
	Resize(r, n, m);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}
	Matrix p;
	Resize(p, n, m);
	p = r;

	Matrix q;
	Resize(q, n, m);
	//Dpoisson(space, r, p, p1, p2, hx, hy);
	
	double error = norm(r);
	double eps = 1e-3;
	int it = 0;
	while (error > eps)
	{
		it++;
		DDpoisson(space, p, q, p1, p2, hx, hy);
		double tau = mlt(p, r) / (mlt(p, q));
		//std::cout << "tau = " << tau << endl;
		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < n; i++)
			{
				x[i][j] += tau * p[i][j];
			}
		}
		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < n; i++)
			{
				r[i][j] -= tau * q[i][j];
			}
		}
		p = r;
		error = norm(r);
	}
	std::cout << "gmres_saad " << endl;
	std::cout << "error =  " <<error<< endl;
	std::cout << "it = " << it << endl;
}

void minres_saad(Matrix& space, Matrix& x, Matrix& rhs, double p1, double p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();

	Matrix Ax;
	Resize(Ax, n, m);
	DDpoisson(space, x, Ax, p1, p2, hx, hy);

	Matrix r;
	Resize(r, n, m);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}
	Matrix p;
	Resize(p, n, m);
	p = r;

	Matrix q;
	Resize(q, n, m);
	//Dpoisson(space, r, p, p1, p2, hx, hy);

	double error = norm(r);
	double eps = 1e-3;
	int it = 0;
	while (error > eps)
	{
		it++;
		DDpoisson(space, p, q, p1, p2, hx, hy);
		double tau = mlt(q, p) / (mlt(q, q));
		//std::cout << "tau = " << tau << endl;
		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < n; i++)
			{
				x[i][j] += tau * p[i][j];
			}
		}
		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < n; i++)
			{
				r[i][j] -= tau * q[i][j];
			}
		}
		p = r;
		error = norm(r);
	}
	std::cout << "min_nev_saad " << endl;
	std::cout << "error =  " << error << endl;
	std::cout << "it = " << it << endl;
}

int main()
{
	Matrix space;
	read_file(space);
	int n = space.size();
	int m = space[0].size();

	//----(x,y)-----//
	double x0 = 0.0;
	double xEnd = 8.0;
	double y0 = -2.0;
	double yEnd = 4.0;

	//-------step of gris-------//
	double hx = (xEnd - x0) / m;
	double hy = (yEnd - y0) / n;

	std::cout << "hx = " << hx << std::endl;
	std::cout << "hy = " << hy << std::endl;
	std::cout << endl;
	std::cout << "n = " << n << endl;
	std::cout << "m = " << m << endl;

	double p11 = 8.0;
	double p22 = 6.0;

	Matrix rhs;
	Resize(rhs, n, m);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j == 0)
			{
				rhs[i][j] = p11/(hx*hx);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				rhs[i][j] = p22/(hx*hx);
			}
		}
	}


	Matrix x;
	Resize(x, n, m);

	simple_iterations(space, x, rhs, 0.01, 0.01, p11, p22, hy, hx);
	//gmres_saad(space, x, rhs, p11, p22, hx, hy);

	/*minres_saad(space, x, rhs, p11, p22, hy, hx);
	minnev(space, x, rhs, n, m, p11, p22, hy, hx);*/

	Matrix x3;
	Resize(x3, n, m);
	//minres_saad(space, x3, rhs, p11, p22, hx, hy);

	/*Matrix x4;
	Resize(x4, n, m);
	simple_iterations(space, x4, rhs, 0.01, 0.01, p11, p22, hy, hx);
	/*std::vector<double>x2;
	x2.resize(n*m);
	std::vector<double>rhs_vect;
	rhs_vect.resize(n*m);
	min_vect(space, x2, rhs_vect, p11, p22, hx, hy);*/

	system("pause");
	return 0;
}
