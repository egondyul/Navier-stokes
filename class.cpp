#include "class.h"

NS::NS()
{
	read_file();
	Resize_nm(u_tilda,n,m-1);
	Resize_nm(v_tilda,n+1,m);
	Resize(u_vect);
	Resize_nm(u,n,m-1);
	Resize_nm(v,n+1,m);
	Resize(rhs);
	Resize(p);
	dt = 0.1;
	nu = 1;
	rho = 1;

	Resize(phi);
	epsilon = 0.05;
}

void NS::read_file()
{
	std::string name;
	std::ifstream file;
	do {
		std::cout << "enter the name of txt: ";
		std::cin >> name;
		file.open(name);
		if (file.good())
		{
			break;
		}
	} while (true);

	if (file.is_open())
	{
		std::string str = {};

		while (std::getline(file, str))
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
	n = space.size();
	m = space[0].size();

	//----(x,y)-----//
	double x0 = 0.0;
	double xEnd = 8.0;
	double y0 = -2.0;
	double yEnd = 4.0;

	//-------step of gris-------//
	hx = (xEnd - x0) / m;
	hy = (yEnd - y0) / n;
	p1 = 8.0;
	p2 = 6.0;

	/*std::string name_phi;
	std::ifstream file_phi;
	do {
		std::cout << "enter the name of txt: ";
		std::cin >> name_phi;
		file_phi.open(name_phi);
		if (file_phi.good())
		{
			break;
		}
	} while (true);

	if (file_phi.is_open())
	{
		std::string str = {};

		while (std::getline(file_phi, str))
		{
			double tmp = 0.0;
			std::stringstream is(str);
			std::vector<double> subspace;
			while (is >> tmp)
			{
				subspace.push_back(tmp);
			}
			phi.push_back(subspace);
		}
	}*/

}

double NS:: u_for_second_step(double u1, double u2, double u3, double u4, double u5, double v1, double v2, double v3, double v4)
{
	double result;
	double ux, uy;
	if (u1 > 0.0)
	{
		ux = (u1 - u3) / hx;
	}
	else {
		ux = (u5 - u1) / hx;
	}

	if ((v1 + v2 + v3 + v4) > 0)
	{
		uy = (u1 - u4) / hy;
	}
	else {
		uy = (u2 - u1) / hy;
	}

	double uu_vu = u1 * ux + (1 / 4)*(v1 + v2 + v3 + v4)*uy;
	result = u1 - dt * uu_vu + 2 * dt*(nu/rho)*(u5 + u3 - 4 * u1 + u4 + u2);
	//result = u1 + 2 * dt*nu*(u5 + u3 - 4 * u1 + u4 + u2);
	return result;
}

double NS:: v_for_second_step(double v1, double v2, double v3, double v4, double v5, double u1, double u2, double u3, double u4)
{
	double result;
	double vx, vy;

	if (v1 > 0.0)
	{
		vy = (v1 - v5) / hy;
	}
	else {
		vy = (v3 - v1) / hy;
	}

	if ((u1 + u2 + u3 + u4) > 0.0)
	{
		vx = (v1 - v4) / hx;
	}
	else {
		vx = (v2 - v1) / hx;
	}

	double uv_vv = (1 / 4)*(u1 + u2 + u3 + u4)*vx + v1 * vy;
	result = v1 - dt * uv_vv + 2 * dt*(nu/rho)*(v2 + v4 - 4 * v1 + v3 + v5);
	//result = v1 + 2 * dt*nu*(v2 + v4 - 4 * v1 + v3 + v5);
	return result;
}

void NS::second_step_velocity()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			//top
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0 && space[i][j + 1] == 1 && space[i + 1][j] == 1)
			{
				//функция для верхнего левого угла
				double u_jminus = (hx / hy)*(v[i][j] - v[i + 1][j]) + u[i][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], -u[i][j], u_jminus, u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = 0.0; //v=0
			}
			//нет проверки что ячейка и верхний правый угол и (j+2) right wall
			//top-right wall 
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i - 1][j] == 0 && space[i][j + 1] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1)
			{
				u_tilda[i][j] = 0.0;
				v_tilda[i][j] = 0.0;
			}
			//top-left wall (j-1)
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i - 1][j] == 0 && space[i][j - 1] == 0 && space[i + 1][j] == 1 && space[i][j + 1] == 1)
			{
				//для верхней границы with left wall где выполняется условие u=0, v=0
				double u_jminus = 0.0;
				u_tilda[i][j] = u_for_second_step(u[i][j], -u[i][j], u_jminus, u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = 0.0;
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i - 1][j] == 0 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && space[i + 1][j] == 1)
			{
				//для верхней границы где выполняется условие u=0, v=0
				u_tilda[i][j] = u_for_second_step(u[i][j], -u[i][j], u[i][j - 1], u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = 0.0;
			}
			else if (space[i][j] == 1 && j == (m - 2) && space[i - 1][j] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1)
			{
				//для верхнего правого угла j=(m-2)
				double u_jplus = -(hx / hy)*(v[i][j + 1] - v[i + 1][j + 1]) + u[i][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], -u[i][j], u[i][j - 1], u[i + 1][j], u_jplus, v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = 0.0;
			}
			//проверить также на второй слой сверху
			//inside
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 1 && space[i - 1][j] == 1 && space[i][j + 1] == 1)
			{
				//для левой внешней границы
				double u_jminus = (hx / hy)*(v[i][j] - v[i + 1][j]) + u[i][j];
				double u_jminus_iminus = (hx / hy)*(v[i - 1][j] - v[i][j]) + u[i - 1][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u_jminus, u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], 0.0, v[i + 1][j], u[i - 1][j], u_jminus_iminus, u_jminus, u[i][j]);

			}
			//right wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i + 1][j] == 1 && space[i][j + 1] == 0 && space[i - 1][j] == 1 && space[i][j - 1] == 1)
			{
				//(j+1)==0
				u_tilda[i][j] = 0.0;
				v_tilda[i][j] = v_for_second_step(v[i][j], -v[i][j], v[i - 1][j], v[i][j - 1], v[i + 1][j], u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], 0.0);
			}
			//left wall 
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i][j - 1] == 0&&space[i+1][j]==1&&space[i-1][j]==1&&space[i][j+1]==1)
			{
				//(j-1)==0
				double u_jminus = 0.0;
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u_jminus, u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], -v[i][j], v[i + 1][j], u[i - 1][j], u[i-1][j-1], 0.0, u[i][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i+1][j+1]==0 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1)
			{
				//[i+1][j+1]==0
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], -u[i][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], v[i + 1][j], u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i - 1][j + 1] == 0 && space[i + 1][j] == 1 && space[i - 1][j] == 1 && space[i][j + 1] == 1 && space[i][j - 1] == 1)
			{
				//[i-1][j+1]==0
				u_tilda[i][j] = u_for_second_step(u[i][j], 0.0, u[i][j - 1], u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], v[i + 1][j], u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1)
			{
				//для внутренней ячейки где всё прекрасно
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], u[i + 1][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], v[i + 1][j], u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j == (m - 2) && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j + 1] == 1 && space[i][j - 1] == 1)
			{
				//для правой внешней границы j==(m-2)
				double u_jplus = (hx / hy)*(v[i + 1][j + 1] - v[i][j + 1]) + u[i][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], u[i + 1][j], u_jplus, v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], v[i + 1][j], u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);

			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i][j - 1] == 1)
			{
				//для правой внешней границы j==(m-1)
				double u_iminus_j = (hx / hy)*(v[i][j] - v[i - 1][j]) + u[i - 1][j - 1];
				double u_j_i = (hx / hy)*(v[i + 1][j] - v[i][j]) + u[i][j - 1];
				v_tilda[i][j] = v_for_second_step(v[i][j], 0.0, v[i - 1][j], v[i][j - 1], v[i + 1][j], u_iminus_j, u[i - 1][j - 1], u[i][j - 1], u_j_i);
			}
			//bottom
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i][j + 1] == 1)
			{
				//для нижнего левого угла
				double u_jminus = (hx / hy)*(v[i][j] - v[i + 1][j]) + u[i][j];
				double u_jminus_iminus = (hx / hy)*(v[i - 1][j] - v[i][j]) + u[i - 1][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u_jminus, -u[i][j], u[i][j + 1], v[i][j + 1], v[i][j], v[i + 1][j], v[i + 1][j + 1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], 0.0, v[i + 1][j], u[i - 1][j], u_jminus_iminus, u_jminus, u[i][j]);
			}
			//bottom-right wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i + 1][j] == 0 && space[i][j + 1] == 0)
			{
				//(j+1)==0
				u_tilda[i][j] = 0.0;
				v_tilda[i][j] = v_for_second_step(v[i][j], -v[i][j], v[i - 1][j], v[i][j - 1], 0.0, u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], 0.0);
			}
			//bottom-left wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i + 1][j] == 0 && space[i][j - 1] == 0)
			{
				double u_jminus = 0.0;
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u_jminus, -u[i][j], u[i][j + 1], v[i][j + 1], v[i][j], 0.0, v[i+1][j+1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], -v[i][j], 0.0, u[i - 1][j], u[i-1][j-1], 0.0, u[i][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i + 1][j] == 0 && space[i + 1][j + 1] == 1 && space[i - 1][j] == 1 && space[i][j + 1] == 1 && space[i][j - 1] == 1)
			{
				//для нижней границы [i+1][j+1]==1
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], -u[i][j], u[i][j + 1], v[i][j + 1], v[i][j], 0.0, v[i+1][j+1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], 0.0, u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 2) && j != (m - 1) && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1)
			{
				//для нижней границы
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], -u[i][j], u[i][j + 1], v[i][j + 1], v[i][j], 0.0, 0.0);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], 0.0, u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j == (m - 2) && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i][j + 1] == 1 && space[i][j - 1] == 1)
			{
				//для нижнего правого угла j+2
				double u_jplus = (hx / hy)*(v[i+1][j+1] - v[i][j + 1]) + u[i][j];
				u_tilda[i][j] = u_for_second_step(u[i][j], u[i - 1][j], u[i][j - 1], -u[i][j], u_jplus, v[i][j + 1], v[i][j], 0.0, v[i+1][j+1]);
				v_tilda[i][j] = v_for_second_step(v[i][j], v[i][j + 1], v[i - 1][j], v[i][j - 1], 0.0, u[i - 1][j], u[i - 1][j - 1], u[i][j - 1], u[i][j]);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i][j - 1] == 1)
			{
				//для нижнего правого угла j+1
				double u_ij = (hx / hy)*(v[i + 1][j] - v[i][j]) + u[i][j - 1];
				double u_iminus_j = (hx / hy)*(v[i][j] - v[i - 1][j]) + u[i - 1][j - 1];
				v_tilda[i][j] = v_for_second_step(v[i][j], 0.0, v[i - 1][j], v[i][j - 1], 0.0, u_iminus_j, u[i - 1][j - 1], u[i][j - 1], u_ij);
			}
		}
	}
}

void NS::divergence()
{
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
				double u_jminus = (hx / hy)*(v_tilda[i][j] - v_tilda[i + 1][j]) + u_tilda[i][j];
				u_vect[i][j] = (1 / hx)*(u_tilda[i][j] - u_jminus) + (1 / hy)*(v_tilda[i][j] - v_tilda[i + 1][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 1))
			{
				u_vect[i][j] = (1 / hx)*(u_tilda[i][j] - u_tilda[i][j - 1]) + (1 / hy)*(v_tilda[i][j] - v_tilda[i + 1][j]);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				double u_jplus = (hx / hy)*(v_tilda[i + 1][j] - v_tilda[i][j]) + u_tilda[i][j-1];
				u_vect[i][j] = (1 / hx)*(u_jplus - u_tilda[i][j - 1]) + (1 / hy)*(v_tilda[i][j] - v_tilda[i + 1][j]);
			}
		}
	}
}

void NS:: third_step()
{
	divergence();
	rhs = u_vect;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			rhs[i][j] = rhs[i][j] * (rho / dt);
		}
	}

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
				rhs[i][j] -= p1 / (hx*hx);
			}
			else if (space[i][j] == 1 && j == (m - 1))
			{
				rhs[i][j] -= p2 / (hx*hx);
			}
		}
	}

	conj_grad_saad();
}

void NS::curve_poisson(Matrix& in, Matrix& out)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			/*========top========*/
			//top-first
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0 && space[i][j + 1] == 1 && space[i + 1][j] == 1)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//top-right wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i - 1][j] == 0 && space[i][j + 1] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1)
			{
				out[i][j] = (-in[i][j - 1] + in[i][j]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//top-left wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i - 1][j] == 0 && space[i][j - 1] == 0 && space[i + 1][j] == 1 && space[i][j + 1] == 1)
			{
				out[i][j] = (in[i][j] - in[i][j + 1]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			// just top
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i - 1][j] == 0 && space[i][j + 1] == 1 && space[i][j - 1] == 1 && space[i + 1][j] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//top-second
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			/*=========inside=======*/
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 1 && space[i - 1][j] == 1)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 1 && space[i - 1][j] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//right wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 1 && space[i][j + 1] == 0)
			{
				out[i][j] = (-in[i][j - 1] + in[i][j]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			//left wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 1 && space[i][j - 1] == 0)
			{
				out[i][j] = (in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 1 && space[i - 1][j] == 1)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + 2 * in[i][j] - in[i + 1][j]) / (hy*hy);
			}
			/*=========bottom========*/
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j == 0)
			{
				out[i][j] = (2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && space[i + 1][j] == 0 && j == (m - 1))
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			//bottom-right wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 0 && space[i][j + 1] == 0)
			{
				out[i][j] = (-in[i][j - 1] + in[i][j]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			//bottom-left wall
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 0 && space[i][j - 1] == 0)
			{
				out[i][j] = (in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}
			else if (space[i][j] == 1 && j != 0 && j != (m - 1) && space[i + 1][j] == 0)
			{
				out[i][j] = (-in[i][j - 1] + 2 * in[i][j] - in[i][j + 1]) / (hx*hx) + (-in[i - 1][j] + in[i][j]) / (hy*hy);
			}

		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			out[i][j] = -out[i][j];
		}
	}
}

void NS::conj_grad_saad()
{
	Matrix Ax;
	Resize(Ax);
	curve_poisson(p, Ax);
	Matrix r;
	Resize(r);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			r[i][j] = rhs[i][j] - Ax[i][j];
		}
	}
	Matrix pp;
	Resize(pp);
	pp = r;
	Matrix Ap;
	Resize(Ap);

	double error = max_norm(r);
	double eps = 1e-14;
	int it = 0;
	int itmax = 100;
	while (error > eps)// && it < itmax)
	{
		it++;
		curve_poisson(pp, Ap);
		double alpha = mlt(r, r) / (mlt(Ap, pp));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				p[i][j] += alpha * pp[i][j];
			}
		}
		double for_betta = mlt(r, r);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				r[i][j] -= alpha * Ap[i][j];
			}
		}
		double betta = mlt(r, r) / for_betta;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				pp[i][j] = r[i][j] + betta * pp[i][j];
			}
		}
		error = max_norm(r);
		//std::cout << "error =  " << error << endl;
	}
	std::cout << "conj_grad_saad " << endl;
	std::cout << "error =  " << error << endl;
	std::cout << "it = " << it << endl;
}

void NS::Resize(Matrix& A)
{
	A.resize(n);
	for (int i = 0; i < n; i++)
	{
		A[i].resize(m);
	}
}

void NS::Resize_nm(Matrix& A, int nn, int mm)
{
	A.resize(nn);
	for (int i = 0; i < n; i++)
	{
		A[i].resize(mm);
	}
}

double NS::mlt(Matrix& a, Matrix& b)
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

double NS:: max_norm(Matrix& r)
{
	double tmp = 0.0;

	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			if (abs(r[i][j]) > tmp)
			{
				tmp = abs(r[i][j]);
			}
		}
	}

	return tmp;
}

void NS::fourth_step_u()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m-1; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && space[i][j + 1] == 1)
			{
				u[i][j] = u_tilda[i][j] - (dt / (rho*hx))*(p[i][j + 1] - p[i][j]);
			}
			
		}
	}
}

void NS::fourth_step_v()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && i != 0 && space[i - 1][j] == 1)
			{
				v[i][j] = v_tilda[i][j] - (dt / (rho*hx))*(p[i - 1][j] - p[i][j]);
			}
			else if (space[i][j] == 1 && i != 0 && space[i - 1][j] == 0)
			{
				v[i][j] = v_tilda[i][j];
			}
			else if (space[i][j] == 1 && i == 0)
			{
				v[i][j] = v_tilda[i][j];
			}
		}
	}
}

void NS::calculate_phi()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 0)
			{
				continue;
			}
			else if (space[i][j] == 1 && j != 0 && j != 1 && j != (m - 1) && j != (m - 2) && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i + 1][j] != 0 && space[i + 2][j] != 0 && space[i][j+1]==1&&space[i][j+2]==1&&space[i][j-1]==1&&space[i][j-2]==1)
			{
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//разобраться с полем скорости!!!!!!!!!!
			//j==0/1
			else if (space[i][j] == 1 && j == 0 && space[i - 2][j] == 1 && space[i - 1][j] == 1 && space[i + 2][j] == 1 && space[i + 1][j] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i][j + 1];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + tmp1)) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(tmp1 * tmp1 * tmp1 + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(tmp1 + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + tmp1 + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			else if (space[i][j] == 1 && j == 1 && space[i - 2][j] == 1 && space[i - 1][j] == 1 && space[i + 2][j] == 1 && space[i + 1][j] == 1)
			{
				double tmp2 = phi[i][j - 1];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//j=(m-2)/(m-1)
			else if (space[i][j] == 1 && j == (m-1) && space[i - 2][j] == 1 && space[i - 1][j] == 1 && space[i + 2][j] == 1 && space[i + 1][j] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i][j - 1];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (tmp1 + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + tmp1 * tmp1 * tmp1 -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + tmp1 - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp2 + tmp1 + phi[i + 2][j] + phi[i + 1][j]);
			}
			else if (space[i][j] == 1 && j == (m-2) && space[i - 2][j] == 1 && space[i - 1][j] == 1 && space[i + 2][j] == 1 && space[i + 1][j] == 1)
			{
				double tmp2 = phi[i][j + 1];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp2 + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=0/1
			else if (space[i][j] == 1 && j != 0 && j != 1 && space[i - 1][j] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1 && space[i][j - 2] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i + 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (tmp1 + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + tmp1 * tmp1 * tmp1) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + tmp1) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp2 + tmp1 - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			else if (space[i][j] == 1 && j != 0 && j != 1 && space[i - 1][j] == 1 && space[i - 2][j] == 0 && space[i + 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1 && space[i][j - 2] == 1)
			{
				double tmp2 = phi[i - 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp2 + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=(n-1)/(n-2)
			else if (space[i][j] == 1 && j != 0 && j != 1 && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1 && space[i][j - 2] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i - 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + tmp1)) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + tmp1 * tmp1 * tmp1 + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + tmp1 + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp2 + tmp1);
			}
			else if (space[i][j] == 1 && j != 0 && j != 1 && space[i + 2][j] == 0 && space[i + 1][j] == 1 && space[i-1][j]==1 && space[i][j - 1] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1 && space[i][j - 2] == 1)
			{
				double tmp2 = phi[i + 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp2 + phi[i + 1][j]);
			}
			//левый верхний угол
			//i=0,j=0
			else if (space[i][j] == 1 && j == 0 && space[i - 1][j] == 0 && space[i + 1][j] == 1 && space[i + 2][j] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i][j + 1];
				double tmp3 = phi[i + 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + tmp1)) -
					(dt / 2 * hy)*(v[i][j] * (tmp1 + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(tmp1 * tmp1 * tmp1 + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + tmp1 * tmp1 * tmp1) -
						(1 - phi[i][j])*(dt / hx * hx)*(tmp1 + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + tmp1) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + tmp1 + tmp3 + tmp1 - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=1,j=0
			else if (space[i][j] == 1 && j == 0 && space[i - 2][j] == 0 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i + 2][j] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i][j + 1];
				double tmp3 = phi[i - 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + tmp1)) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(tmp1 * tmp1 * tmp1 + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(tmp1 + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + tmp1 + tmp3 + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=1,j=1
			else if (space[i][j] == 1 && j == 1 && space[i - 2][j] == 0 && space[i - 1][j] == 1 && space[i + 1][j] == 1 && space[i + 2][j] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
				double tmp1 = phi[i][j - 1];
				double tmp2 = phi[i - 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp1 + phi[i][j - 1] + tmp2 + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=0,j=1
			else if (space[i][j] == 1 && j == 1 && space[i - 1][j] == 0 && space[i + 1][j] == 1 && space[i + 2][j] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
				double tmp1 = phi[i][j-1];
				double tmp2 = phi[i][j];
				double tmp3 = phi[i + 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
					(dt / 2 * hy)*(v[i][j] * (tmp2 + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + tmp2 * tmp2 * tmp2) -
						(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + tmp2) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp1 + phi[i][j - 1] + tmp3 + tmp2 - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//левый нижний угол
			//i=(n-1),j=0
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] == 0 && space[i - 1][j] == 1 && space[i - 2][j] == 1 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
				double tmp1 = phi[i][j];
				double tmp2 = phi[i][j + 1];
				double tmp3 = phi[i - 1][j];
				phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + tmp1)) -
					(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + tmp1)) +
					(1 - phi[i][j])*(dt / hx * hx)*(tmp1 * tmp1 * tmp1 + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
						4 * phi[i][j] * phi[i][j] * phi[i][j] + tmp1 * tmp1 * tmp1 + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
						(1 - phi[i][j])*(dt / hx * hx)*(tmp1 + phi[i][j + 1] - 4 * phi[i][j] + tmp1 + phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + tmp1 + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp3 + tmp1);
			}
			//i=(n-2),j=0
			else if (space[i][j] == 1 && j == 0 && space[i + 1][j] != 0 && space[i + 2][j] == 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 &&  space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
			double tmp1 = phi[i][j];
			double tmp2 = phi[i][j + 1];
			double tmp3 = phi[i + 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + tmp1)) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(tmp1 *tmp1 * tmp1 + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp2 + tmp1 + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp3 + phi[i + 1][j]);
			}
			//i=(n-2),j=1
			else if (space[i][j] == 1 && j == 1 && space[i + 1][j] != 0 && space[i + 2][j] == 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0  && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
			double tmp1 = phi[i][j - 1];
			double tmp2 = phi[i + 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp1 + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp2 + phi[i + 1][j]);
			}
			//i=(n-1),j=1
			else if (space[i][j] == 1 && j == 1 && space[i + 1][j] == 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i][j + 1] == 1 && space[i][j + 2] == 1)
			{
			double tmp1 = phi[i][j - 1];
			double tmp2 = phi[i][j];
			double tmp3 = phi[i - 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + tmp2)) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + tmp2 * tmp2 * tmp2 + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + tmp2 + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(tmp1 + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + phi[i][j + 2] + phi[i][j + 1] + tmp3 + tmp2);
			}
			//правый верхний угол
			//i=0,j=(m-1)
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] == 0 && space[i + 1][j] != 0 && space[i + 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j];
			double tmp2 = phi[i][j - 1];
			double tmp3 = phi[i + 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (tmp1 + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (tmp1 + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + tmp1 * tmp1 * tmp1 -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + tmp1 * tmp1 * tmp1) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + tmp1 - 4 * phi[i][j] + phi[i + 1][j] + tmp1) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp3 + tmp1 - 8 * phi[i][j] + tmp2 + tmp1 + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=0,j=(m-2)
			else if (space[i][j] == 1 && j == (m - 2) && space[i - 1][j] == 0 && space[i + 1][j] != 0 && space[i + 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j + 1];
			double tmp2 = phi[i][j];
			double tmp3 = phi[i + 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (tmp2 + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + tmp2 * tmp2 * tmp2) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + tmp2) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp3 + tmp2 - 8 * phi[i][j] + tmp1 + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=1,j=(m-2)
			else if (space[i][j] == 1 && j == (m-2) && space[i - 1][j] == 0 && space[i + 1][j] != 0 && space[i + 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j + 1];
			double tmp2 = phi[i-1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp2 + phi[i - 1][j] - 8 * phi[i][j] + tmp1 + phi[i][j + 1] + phi[i + 2][j] + phi[i + 1][j]);
			}
			//i=1,j=(m-1)
			else if (space[i][j] == 1 && j == (m - 1) && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i + 1][j] != 0 && space[i + 2][j] != 0 && space[i][j + 1] == 1 && space[i][j + 2] == 1 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j];
			double tmp2 = phi[i][j-1];
			double tmp3 = phi[i-1][j];

			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (tmp1 + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + tmp1 * tmp1 * tmp1 -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + tmp1 - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + tmp3 + phi[i - 1][j] - 8 * phi[i][j] + tmp2 + tmp1 + phi[i + 2][j] + phi[i + 1][j]);
			}
			//правый нижний угол
			//i=(n-1),j=(m-1)
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 1][j] == 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j];
			double tmp2 = phi[i][j - 1];
			double tmp3 = phi[i - 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (tmp1 + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + tmp1)) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + tmp1 * tmp1 *tmp1 -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + tmp1 * tmp1 * tmp1 + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + tmp1 - 4 * phi[i][j] + tmp1 + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp2 + tmp1 + tmp3 + tmp1);
			}
			//i=(n-2), j=(m-1)
			else if (space[i][j] == 1 && j == (m - 1) && space[i + 2][j] == 0 && space[i + 1][j] != 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j];
			double tmp2 = phi[i][j - 1];
			double tmp3 = phi[i + 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (tmp1 + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + tmp1 * tmp1 * tmp1 -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + tmp1 - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp2 + tmp1 + tmp3 + phi[i + 1][j]);
			}
			//i=(n-2),j=(m-2)
			else if (space[i][j] == 1 && j == (m - 2) && space[i + 2][j] == 0 && space[i + 1][j] != 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j + 1];
			double tmp2 = phi[i+1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + phi[i + 1][j])) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + phi[i + 1][j] * phi[i + 1][j] * phi[i + 1][j] + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + phi[i + 1][j] + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp1 + phi[i][j + 1] + tmp2 + phi[i + 1][j]);
			}
			//i=(n-1),j=(m-2)
			else if (space[i][j] == 1 && j == (m - 2) && space[i + 1][j] == 0 && space[i - 1][j] != 0 && space[i - 2][j] != 0 && space[i][j - 1] == 1 && space[i][j - 2] == 1)
			{
			double tmp1 = phi[i][j + 1];
			double tmp2 = phi[i][j];
			double tmp3 = phi[i - 1][j];
			phi[i][j] += -(dt / 2 * hx)*(u[i][j] * (phi[i][j + 1] + phi[i][j]) - u[i][j - 1] * (phi[i][j] + phi[i][j - 1])) -
				(dt / 2 * hy)*(v[i][j] * (phi[i - 1][j] + phi[i][j]) - v[i + 1][j] * (phi[i][j] + tmp2)) +
				(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] * phi[i][j - 1] * phi[i][j - 1] + phi[i][j + 1] * phi[i][j + 1] * phi[i][j + 1] -
					4 * phi[i][j] * phi[i][j] * phi[i][j] + tmp2 * tmp2 * tmp2 + phi[i - 1][j] * phi[i - 1][j] * phi[i - 1][j]) -
					(1 - phi[i][j])*(dt / hx * hx)*(phi[i][j - 1] + phi[i][j + 1] - 4 * phi[i][j] + tmp2 + phi[i - 1][j]) -
				(1 - phi[i][j])*(dt / hx * hx*hx*hx)*epsilon*epsilon*(phi[i][j - 2] + phi[i][j - 1] + phi[i - 2][j] + phi[i - 1][j] - 8 * phi[i][j] + tmp1 + phi[i][j + 1] + tmp3 + tmp2);
			}
		}
	}
}







//эксперименты, тесты

void NS::cycle()
{
	for (int i = 0; i < 500; i++)
	{
		//print();
		second_step_velocity();
		third_step();
		fourth_step_u();
		fourth_step_v();
		change_dt();

		std::ofstream outp;
		outp.open("p" + std::to_string(i) + ".txt");
		if (outp.is_open())
		{
			for (int it = 0; it < n; it++)
			{
				for (int jt = 0; jt < m; jt++)
				{
					outp << p[it][jt] << " ";
				}
				outp << std::endl;
			}
		}

		std::ofstream outu;
		outu.open("u" + std::to_string(i) + ".txt");
		if (outu.is_open())
		{
			for (int it = 0; it < n; it++)
			{
				for (int jt = 0; jt < m-1; jt++)
				{
					outu << u[it][jt] << " ";
				}
				outu << std::endl;
			}
		}

		std::ofstream outv;
		outv.open("v" + std::to_string(i) + ".txt");
		if (outv.is_open())
		{
			for (int it = 0; it < n; it++)
			{
				for (int jt = 0; jt < m; jt++)
				{
					outv << v[it][jt] << " ";
				}
				outv << std::endl;
			}
		}
	}
}

void NS::print()
{
	std::cout << "First step: " << std::endl;
	std::cout << "u_tilda: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m - 1; j++)
		{
			std::cout << u_tilda[i][j] << " ";
		}
		std::cout << endl;
		if (i == (n - 1))
		{
			std::cout << std::endl;
		}
	}
	std::cout << "v_tilda: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << v_tilda[i][j] << " ";
		}
		std::cout << endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Second step: " << std::endl;

	std::cout << "presuure: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << p[i][j] << " ";
		}
		std::cout << endl;
		if (i == (n - 1))
		{
			std::cout << std::endl;
			std::cout << std::endl;
		}
	}

	
	std::cout << "Third step: " << std::endl;

	std::cout << "u: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m-1; j++)
		{
			std::cout << u[i][j] << " ";
		}
		std::cout << endl;
		if (i == (n - 1))
		{
			std::cout << std::endl;
		}
	}

	std::cout << "v: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m ; j++)
		{
			std::cout << v[i][j] << " ";
		}
		std::cout << endl;
	}
	std::cout << std::endl;
	std::cout << "NEXT" << std::endl;
	
}

void NS::change_dt()
{
	double u_max = 0.0;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (space[i][j] == 1&&j!=(m-2)&&j!=(m-1))
			{
				if (sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) > u_max)
				{
					u_max = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);
				}
			}
		}
	}

	if (u_max!=0.0 && dt > (hx / (sqrt(2)*u_max)))
	{
		dt = hx / (sqrt(2)*u_max);
	}

	std::cout << "dt 1 =" << dt << endl;
	std::cout << "u_max =" << u_max;
	std::cout << endl;

	/*if (nu!=0 && dt > (hx*hx / (4 * nu)))
	{
		dt = (hx*hx / (4 * nu));
	}

	std::cout << "dt 2 =" << dt;*/
}

void NS::test_rho_nu()
{
	double epss = 10e-4;
	double k = 10.0;
	int i = 0;
	//изменение плотности и вязкости
	rho = 1;
	nu = 1;
	p1 = 8;
	p2 = 6;

	std::ofstream out_file;
	out_file.open("inform.txt");
	if (out_file.is_open())
	{
		out_file << "dp= " << (p2 - p1) << std::endl;
		out_file << "rho= " << rho << std::endl;
		out_file << "nu= " << nu << std::endl;
		out_file << std::endl;
	}

	auto dt_min = 0.0;
	auto dt_max = 0.0;

	auto dt_P_min = 0.0;
	auto dt_P_max = 0.0;

	while(k>epss)
	{
		i++;

		auto time1 = std::chrono::steady_clock::now();
		second_step_velocity();

		auto begin_time = std::chrono::steady_clock::now();
		third_step();
		auto end_time = std::chrono::steady_clock::now();
		auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time);
		//std::cout << "The time of regex = " << elapsed_ms.count() << endl;

		double u_n = 0.0;
		for (int ii = 0; ii < n; ii++)
		{
			for (int jj = 0; jj < m-1; jj++)
			{
				u_n += sqrt(u[ii][jj]* u[ii][jj]+ v[ii][jj]* v[ii][jj])*(hx*hy);
			}
		}
		fourth_step_u();
		fourth_step_v();
		change_dt();

		auto time2 = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1);
		//std::cout << "The time of regex = " << ms.count() << endl;

		double u_n_plus = 0.0;
		for (int ii = 0; ii < n; ii++)
		{
			for (int jj = 0; jj < m-1; jj++)
			{
				u_n_plus += sqrt(u[ii][jj] * u[ii][jj] + v[ii][jj] * v[ii][jj]) * (hx*hy);
			}
		}

		//критерий установления течения
		if (u_n != 0)
		{
			k = abs((u_n_plus - u_n) / u_n);
		}

		/*std::ofstream out_file;
		out_file.open("inform_1.txt");
		if (out_file.is_open())
		{
			out_file << "Number of iteration= " << i << std::endl;
			out_file << "dt of Poisson= " << elapsed_ms.count() << std::endl;
			out_file << "dt of iteration= " << ms.count() << std::endl;
			out_file << std::endl;

		}
		out_file.close();*/
		std::cout << "Number of iteration= " << i << std::endl;
		std::cout << "dt of Poisson= " << elapsed_ms.count() << std::endl;
		std::cout << "dt of iteration= " << ms.count() << std::endl;
		std::cout << std::endl;

	}

	
}

