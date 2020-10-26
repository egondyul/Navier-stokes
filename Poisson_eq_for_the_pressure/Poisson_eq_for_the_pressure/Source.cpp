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

/*class Pressure {
private:
	std::vector<double> p1;//краевые условия
	std::vector<double> p2;
	Matrix p;
public:

};*/

void read_file(Matrix& space)
{
	char name[1000] = { "input.txt" };
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

/*Matrix matrix_algorithm(Matrix& space)
{
	Matrix result;
	double x_min = 0;
	double y_min = -2;
	double x_max = 8;
	double y_max = 4;

	int n = space.size();
	int m = space[0].size();

	double h = 0.2;

	Matrix rhs;
	rhs.resize(n);
	for (int i = 0; i < n; i++)
	{
		rhs[i].resize(m);
	}

	Matrix I;
	I.resize(n);
	for (int i = 0; i < n; i++)
	{
		I[i].resize(m);
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j)
			{
				I[i][j] = 1.0;
			}
		}
	}

	Matrix C;
	C.resize(n);
	for (int i = 0; i < n; i++)
	{
		C[i].resize(n);
	}
	for (int i = 0; i < n; i++)
	{
		C[i][i] = -4.0;
		C[i][i + 1] = 1.0;
		C[i][i - 1] = 1.0;
	}

	Matrix alpha;
	alpha.resize(n);
	for (int i = 0; i < n; i++)
	{
		alpha[i].resize(n);
	}


	

	result.resize(n+2);
	for (int i = 0; i < n; i++)
	{
		result[i].resize(m+2);
	}

	// граничные условия
	std::vector<double> p1,p2;
	p1.resize(n+2);
	p2.resize(n+2);
	for (int i = 0; i < n+2; i++)
	{
		p1[i] = 80000.0;
		p2[i] = 60000.0;
	}
	for (int i = 0; i < n+2; i++)
	{
		result[i][0] = p1[i];
		result[i][1] = p1[i];
		result[i][m + 1] = p2[i];
		result[i][m] = p2[i];
	}

	Matrix sol;
	sol.resize(n + 2);
	for (int i = 0; i < n; i++)
	{
		sol[i].resize(m + 2);
	}

	for (int j = 1; j < m+2; j++)
	{
		for (int i = 1; i < m; i++)
		{
			if (space[i-1][j-1] > 0.0)
			{
				if (space[i-1][j] == 0.0)
				{
					result[i][j + 1] = result[i + 1][j + 1];
				}
				if (space[i + 1][j] == 0.0)
				{
					result[i+2][j+1]=result[]
				}
				double tmp = sol[i][j];
				result[i][j + 1] = tmp*(h*h) - result[i - 1][j] - result[i][j - 1] - result[i + 1][j] + 4 * result[i][j];

			}
		}
	}


}*/

Matrix pressure(Matrix& space)
{
	Matrix result;

	double h = 0.2;

	int n_old = space.size();
	int m_old;
	for (int i = 0; i < n_old; i++)
	{
		m_old = space[i].size();
	}

	int n = n_old + 2;
	int m = m_old + 2;

	Matrix space_new;
	space_new.resize(n);
	for (int i = 0; i < n; i++)
	{
		space_new[i].resize(m);
	}
	for (int i = 1; i < n - 1; i++)
	{
		for (int j = 1; j < m - 1; j++)
		{
			space_new[i][j] = space[i - 1][j - 1];
		}
	}

	result.resize(n);
	for (int i = 0; i < n; i++)
	{
		result[i].resize(m);
	}

	std::vector<double> p1, p2;
	p1.resize(n);
	p2.resize(n);
	for (int i = 0; i < n; i++)
	{
		p1[i] = 80000.0;
		p2[i] = 60000.0;
	}
	for (int i = 0; i < n; i++)
	{
		result[i][0] = p1[i];
		result[i][1] = p1[i];
		result[i][m-1] = p2[i];
		result[i][m-2] = p2[i];
	}

	Matrix sol;
	sol.resize(n);
	for (int i = 0; i < n; i++)
	{
		sol[i].resize(m);
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			sol[i][j] = 1.0;
		}
	}

	for (int j = 1; j < m-1; j++)
	{
		for (int i = 1; i < n-1; i++)
		{
			if (space_new[i][j] > 0.0)
			{
				double tmp = sol[i][j];
				result[i][j + 1] = tmp * (h*h) - result[i - 1][j] - result[i][j - 1] - result[i + 1][j] + 4 * result[i][j];
			}
		}
	}

	return result;
}

int main()
{
	Matrix space; //записываем 0 и 1
	read_file(space);
	pressure(space);

	return 0;
}
