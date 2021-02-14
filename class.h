#pragma once
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

class NS
{
private:
	Matrix space;
	int n;
	int m;
	Matrix u;
	Matrix v;
	Matrix p;

	Matrix phi;
public:
	NS();
	void cycle();
	void test_rho_nu();
private:
	Matrix u_tilda; //inter field of u
	Matrix v_tilda; //inter field of v
	Matrix u_vect; // divergence of velocity
	Matrix rhs; //rhs for third step
	double hx;
	double hy;
	double p1;
	double p2;
	double dt;
	double nu;
	double rho;

	double epsilon;

private:
	void read_file();
	void second_step_velocity();
	void third_step();
	void fourth_step_u();
	void fourth_step_v();
	void print();
private:
	//for every step
	void Resize(Matrix& A);
	void Resize_nm(Matrix& A, int nn, int mm);
	//for second step
	double u_for_second_step(double u1, double u2, double u3, double u4, double u5, double v1, double v2, double v3, double v4);
	double v_for_second_step(double v1, double v2, double v3, double v4, double v5, double u1, double u2, double u3, double u4);
	//for third step
	void divergence(); //rhs
	void conj_grad_saad(); //SLAU
	void curve_poisson(Matrix& in, Matrix& out); // poisson
	double mlt(Matrix& a, Matrix& b);
	double max_norm(Matrix& r);
	void change_dt();
	//test

	void calculate_phi();
};


