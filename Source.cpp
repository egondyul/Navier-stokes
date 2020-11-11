#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <math.h>
#include<sstream>
using namespace std;

typedef std::vector<std::vector<double>> Matrix;

void read_file(Matrix& space)
{
	char name[1000] = { "test1.txt" };
	std::ifstream in;
	in.open(name);
	if (!in.is_open())
	{
		cout << "Sorry :(";
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

void poisson(Matrix& space, Matrix& in, Matrix& out, std::vector<double>& p1, std::vector<double>& p2, double hx, double hy)
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
				out[i][j] = (in[i + 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + p1[i]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0)
			{
				//second
				out[i][j] = (in[i + 1][j] - in[i][j]) / hx + (p2[i] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0)
			{
				//third
				out[i][j] = (in[i - 1][j] - in[i][j]) / hx + (in[i][j + 1] - 2 * in[i][j] + p1[i]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0)
			{
				//fourth
				out[i][j] = (in[i - 1][j] - in[i][j]) / hx + (p2[i] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == 0)
			{
				//dirichlet-right
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (in[i][j + 1] - 2 * in[i][j] + p1[i]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				//dirichlet-left
				out[i][j] = (in[i + 1][j] - 2 * in[i][j] + in[i - 1][j]) / (hx*hx) + (p2[i] - 2 * in[i][j] + in[i][j - 1]) / (hy*hy);
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
}

void Resize(Matrix& vect, int n, int m)
{
	vect.resize(n);
	for (int i = 0; i < n; i++)
	{
		vect[i].resize(m);
	}
}

double norma_remains(Matrix& a)
{
	int n = a.size();
	int m = a[0].size();
	
	double norma = 0.0;
	for (int i = 0; i < n; i++)
	{
		double s = 0.0;
		for (int j = 0; j < m; j++)
		{
			s += abs(a[i][j]);
		}
		if (s > norma)
		{
			norma = s;
		}
	}
	return norma;
}

void simple_iterations(Matrix& space, Matrix& x0, Matrix& rhs, double t0, double p0, std::vector<double>& p1, std::vector<double>& p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();
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

	Matrix x;
	Resize(x, n, m);
	double error = norma_remains(r);
	double eps = 1e-3;
	int endCycle = 1000;
	int it = 0;
	while (error > eps && it < endCycle)
	{
		/*double pi = 3.14159;
		double tn = cos((2 * it - 1)*pi / (2 * endCycle));
		double t = t0 / (1 + p0 * tn);
		cout << "t = " << t << endl;*/
		it++;
		poisson(space, x0, Ax, p1, p2, hx, hy);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = x0[i][j] - t0*(Ax[i][j] - rhs[i][j]);
			}
		}
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] = rhs[i][j] - Ax[i][j];
			}
		}
		error = norma_remains(r);
		x0 = x;
	}
	cout << "iterations = " << it << endl;
}

/*void min_remains(Matrix& space, Matrix& x0, Matrix& rhs, std::vector<double>& p1, std::vector<double>& p2, double hx, double hy)
{
	int n = space.size();
	int m = space[0].size();
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

	Matrix x;
	Resize(x, n, m);

	double error = norma_remains(r);
	double eps = 1e-3;
	int endCycle = 100;
	int it = 0;
	while (error > eps && it < endCycle)
	{
		it++;
		double t_p=
	}
}*/


int main()
{
	Matrix space; //записываем 0 и 1
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

	cout << "hx = " << hx << endl;
	cout << "hy = " << hy << endl;

	int what;
	cout << "What kind of test would you like? (1)-linear and (2)- sin: ";
	cin >> what;

	if (what == 1)
	{
		Matrix rhs;
		Resize(rhs, n, m);

		double p11 = 8.0;
		double p22 = 6.0;

		std::vector<double> p1, p2;
		p1.resize(n);
		p2.resize(n);
		for (int i = 0; i < n; i++)
		{
			p1[i] = p11;
			p2[i] = p22;
		}

		double t0 = -0.01;
		double p0 = -0.9971;

		Matrix x;
		Resize(x, n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (space[i][j] == 1)
				{
					x[i][j] = p22;
				}
			}
		}
		simple_iterations(space, x, rhs, t0, p0, p1, p2, hx, hy);
	}
	else if (what == 2)
	{
		Matrix Sin;
		Resize(Sin, n, m);
		std::vector<double>p1;
		std::vector<double>p2;

		Matrix big;
		char name[1000] = { "sin.txt" };
		std::ifstream in;
		in.open(name);
		if (!in.is_open())
		{
			cout << "Sorry :(";
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
				big.push_back(subspace);
			}
		}
		for (int i = 0; i < n; i++)
		{
			p1.push_back(big[i][0]);
		}
		for (int i = 0; i < n; i++)
		{
			p2.push_back(big[i][m + 1]);
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (space[i][j] == 1)
				{
					Sin[i][j] = -big[i][j + 1];
				}
			}
		}
		double t0 = -0.01;
		double p0 = -0.9971;
		Matrix x;
		Resize(x, n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				x[i][j] = 0.0;
			}
		}
		simple_iterations(space, x, Sin, t0, p0, p1, p2, hx, hy);
	}
	

	system("pause");
	return 0;
}