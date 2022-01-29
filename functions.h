#include<iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXd function_y(MatrixXd para, double area_1, double area_2, double len, double t_total, double ForceMax);
MatrixXd getJacobianMatrix(MatrixXd para_est, MatrixXd deflection, MatrixXd ym, double area_1, double area_2, double len, double t_total, double ForceMax);
MatrixXd GaussNewton(MatrixXd para_guess, MatrixXd deflection, MatrixXd ym, double area_1, double area_2, double len, double t_total, double ForceMax);
