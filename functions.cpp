#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<math.h>
#include"functions.h"
#include"fea_main.h"


using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd function_y(MatrixXd para, double area_1, double area_2, double len, double t_total, double ForceMax) {

	MatrixXd y = fea_main(para, area_1, area_2, len, t_total, ForceMax);

	return y;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd getJacobianMatrix(MatrixXd para_est, MatrixXd deflection, MatrixXd ym, double area_1, double area_2, double len, double t_total, double ForceMax) {

	MatrixXd Jacobian_Matrix(ym.rows(), para_est.rows());
	MatrixXd y = function_y(para_est, area_1, area_2, len, t_total, ForceMax);
	MatrixXd y_deflected(ym.rows(), 1);

	for (int i = 0; i < para_est.rows(); i++) {

		para_est(i, 0) = para_est(i, 0) + deflection(i, 0);		/*Changing the parameters one by one */

		y_deflected = function_y(para_est, area_1, area_2, len, t_total, ForceMax);				/*Computing the deflected function arrray */
		for (int j = 0; j < ym.rows(); j++) {

			// [f(v, p + dp) - f(v, p) ] / [dp] 

			Jacobian_Matrix(j, i) = (y_deflected(j, 0) - y(j, 0)) / deflection(i, 0);
		}
		para_est(i, 0) = para_est(i, 0) - deflection(i, 0);		/*Bringing back the parametes to original value*/
	}
	return Jacobian_Matrix;
};
/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/


MatrixXd GaussNewton(MatrixXd para_guess, MatrixXd deflection, MatrixXd ym, double area_1, double area_2, double len, double t_total, double ForceMax) {

	cout << "-> Entered Gauﬂ-Newton\n";

	//cout << "para_guess: " << para_guess << endl;
	/*cout << "deflection: " << deflection << endl;
	cout << "ym: " << ym << endl;
	cout << "input: " << input << endl;*/

	int npara = para_guess.rows(), ndata = ym.rows();

	//MatrixXd IdentityMat = MatrixXd::Identity(npara, npara);
	MatrixXd H(npara, npara);
	MatrixXd d(ndata, 1);
	MatrixXd J(ndata, npara);
	double error,error_gn;

	MatrixXd y_init = function_y(para_guess, area_1, area_2, len, t_total, ForceMax);
	//cout << "yinit: " << y_init << endl;

	//double lambda = 10;
	int updateJ = 1;
	MatrixXd para_est = para_guess;
	int maxiter = 1000, counter=0;
	

	while (counter < maxiter) {
		
		cout << "--> Iteration : " << counter << endl;
	
		J = getJacobianMatrix(para_est, deflection, ym, area_1, area_2, len, t_total, ForceMax);
		//cout << "J: \n" << J << endl;
		MatrixXd y_est = function_y(para_est, area_1, area_2, len, t_total, ForceMax);
		d = ym - y_est;
		//cout << "d: \n" << d << endl;
		H = J.transpose() * J;
		cout << "H: \n" << H << endl;

		if (counter == 0) {
			MatrixXd temp1 = d.transpose() * d;
			error = temp1(0, 0);
			//cout << "error" << error << endl;
		}
				
		//cout << "Hinverse=\n" << H_lm.completeOrthogonalDecomposition().pseudoInverse() << endl;
		MatrixXd dp = H.completeOrthogonalDecomposition().pseudoInverse() * J.transpose() * d;
		cout << "dp: \n" << dp;
		MatrixXd para_gn = para_est + dp;
		MatrixXd y_est_gn = function_y(para_gn, area_1, area_2, len, t_total, ForceMax);
		MatrixXd d_gn = ym - y_est_gn;
		MatrixXd temp2 = d_gn.transpose() * d_gn;
		error_gn = temp2(0,0);

		para_est = para_gn;
		error = error_gn;
		
		if (dp.norm() < 1e-6) {
			counter = 1000;
		}
		else counter++;
		//cout << "Para:\n" << para_est << endl;
		cout << "dp:\n" << dp << endl;
	}

	return para_est;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
