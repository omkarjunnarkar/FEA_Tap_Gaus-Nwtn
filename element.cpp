#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include<math.h>
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd B_mat(double le) {
	MatrixXd b(1, 2);
	b(0, 0) = -1 / le;
	b(0, 1) = 1 / le;
	return b;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

int signum(double x) {
	if (x >= 0) return 1;
	else return -1;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*To compute the assignment matrix for each element*/
MatrixXd assign(int ele_index_cpp, int numele) {

	//ele_index is current element number, not defined here!! different for each ele !!

	MatrixXd a_element = MatrixXd::Zero(2, numele + 1);
	a_element(0, ele_index_cpp) = 1;
	a_element(1, ele_index_cpp + 1) = 1;
	return a_element;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

//tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double a_1, double a_2, double len_1, double len_2) {
//	MatrixXd areas(elements, 1);
//	MatrixXd lengths(elements, 1);
//	for (int i = 0; i < elements; i++) {
//		if (i < elements / 2) {
//			areas(i, 0) = a_1;
//			lengths(i, 0) = len_1 / (elements / 2);
//		}
//		else {
//			areas(i, 0) = a_2;
//			lengths(i, 0) = len_2 / (elements / 2);
//		}
//	}
//	
//	return { areas,lengths };
//}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double ar_1, double ar_2, double len) {
	MatrixXd areas(elements, 1);
	MatrixXd lengths(elements, 1);
	double len_per_ele = len / elements;
	double r1 = pow(ar_1 / 3.14, 0.5), r2 = pow(ar_2 / 3.14, 0.5);
	double theta = (r2 - r1) / len, y_n;
	for (int i = 0; i < elements; i++) {
		lengths(i, 0) = len_per_ele;
		y_n = r1 + theta * (i + 1) * len_per_ele;
		areas(i, 0) = 3.14 * pow(y_n, 2);
	}

	return { areas,lengths };
}


