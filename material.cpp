#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

double Lambda(double sigma, double sigma_0) {
	if (abs(sigma) <= sigma_0) { return 0; }
	else return ((abs(sigma) / sigma_0) - 1);
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<double, double, double>get_Stress_Strains(double eps_p_m, double eps, double eta, double h, double t_total, double sigma_0, double E) {
	
	int s = (t_total / h);
	double Ct = E;
	double sigma_m,eps_p_next, eps_p_prev, sigma_next;
	
	// 1 step Euler Explicit

	for (int i = 0; i<s; i++) {
		sigma_m = E * (eps - eps_p_m);
		eps_p_next = eps_p_m + h * eta * signum(sigma_m) * Lambda(sigma_m, sigma_0);

		// 1000 Steps Euler Implicit

		for (int j = 0; j < 1000; j++) {
			eps_p_prev = eps_p_next;
			sigma_next = E * (eps - eps_p_next);
			eps_p_next = eps_p_m + h * eta * signum(sigma_next) * Lambda(sigma_next, sigma_0);

			if (abs(eps_p_next - eps_p_prev) < 1e-5) {
				break;
			}
		}
		
		sigma_m = E * (eps - eps_p_next);						//Updating final values of STress and Plastic Strain
		eps_p_m = eps_p_next;

	}

	if (abs(sigma_m) > sigma_0) {
		Ct = (E * sigma_0) / (sigma_0 + E * eta * t_total);		//Updating Tangent Stiffness for Plastic Case
	}
	else Ct = E;

	return { sigma_m, eps_p_next, Ct };//,sigmas, ep_p_vals
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>material_routine(int elements, double gauss_weights, MatrixXd Lengths, MatrixXd Areas, MatrixXd  E_mat, MatrixXd u, double time_step, double eta, MatrixXd eps_p, double sigma_0, double E) {

	MatrixXd F_int= MatrixXd::Zero(elements + 1, 1);
	MatrixXd stress_mat = MatrixXd::Zero(elements , 1);
	MatrixXd K_t = MatrixXd::Zero(elements+1 , elements+1);
	MatrixXd eps = MatrixXd::Zero(elements , 1);
	MatrixXd u_temp(2, 1);

	for (int ele = 0; ele < elements; ele++) {
		
		u_temp(0, 0) = u(ele,0);
		u_temp(1, 0) = u(ele + 1, 0);
		MatrixXd t = B_mat(Lengths(ele,0)) * u_temp;
		eps(ele, 0) = t(0,0);

		double smt, ept, et;
		tie(smt, ept, et) = get_Stress_Strains(eps_p(ele, 0), eps(ele, 0), eta, time_step / 100, time_step, sigma_0, E);
		stress_mat(ele, 0) = smt;
		eps_p(ele, 0) = ept;
		E_mat(ele, 0) = et;
		
		MatrixXd F_int_temp = B_mat(Lengths(ele, 0)).transpose() * Areas(ele, 0) * stress_mat(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;
		MatrixXd K_t_temp = B_mat(Lengths(ele, 0)).transpose() * B_mat(Lengths(ele, 0)) * E_mat(ele, 0) * Areas(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;

		//Converting from Local to Global at one Gauss Point

		F_int(ele, 0) = F_int(ele, 0) + F_int_temp(0, 0);
		F_int(ele + 1, 0) = F_int(ele + 1, 0) + F_int_temp(1, 0);

		MatrixXd A = assign(ele, elements);
		K_t = K_t + (A.transpose() * K_t_temp * A);

	}

	return { K_t, F_int, eps_p, eps, stress_mat };

}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/