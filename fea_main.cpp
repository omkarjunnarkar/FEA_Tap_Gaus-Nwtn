#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"
#include"fea_main.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd fea_main(MatrixXd Parameters, double A_1, double A_2, double Len, double t_tot, double F_max) {

	//cout << "Starting the Finite Element Program" << endl;
	//Creating csv files for storing stress, train and displacement data

	//ofstream mydisplacefile("displacement.csv");
	//ofstream mystressfile("stress.csv");
	//ofstream mystrainfile("strain.csv");

	//double F_max = 3900;			//Maximum Force
	//double A_1 = 6;					//csArea of left element
	//double A_2 = 12;				//csArea of right element
	//double L_1 = 20;				//Length of left element
	//double L_2 = 40;				//Length of right element
	//double E = 100000;				//Youngs Modulus
	int steps_num = 10000;			//Number of Steps to solve the problem
	double steps = steps_num;
	//double t_tot = 0.5;				//Total time for which force is applied
	//double sigma_0 = 200;			//Yield Strength
	int elements = 25;				//Number of elements in which problem is to be divided
	double force = 0;					
	int gauss_weights = 2;
	//double eta = 1;					//Viscosity
	MatrixXd Areas, Lengths;

	double E = 1000*Parameters(0, 0);
	double sigma_0 = Parameters(1, 0);
	double eta = Parameters(2, 0);

	steps = t_tot / double(steps);

	//Calling Element Routine to get Areas & Lengths of all the elements

	tie(Areas, Lengths) = get_areas_lengths(elements, A_1, A_2, Len);	
	
	MatrixXd F_ext = MatrixXd::Zero(elements + 1,1);
	MatrixXd E_mat= E*MatrixXd::Ones(elements , 1);
	MatrixXd eps_p = MatrixXd::Zero(elements , 1);
	MatrixXd eps = MatrixXd::Zero(elements , 1);
	MatrixXd u = MatrixXd::Zero(elements + 1, 1);
	MatrixXd K_t, F_int, strain, stress;
	MatrixXd final_displacement_cpp_node = MatrixXd::Zero(steps_num+2, 1);
	
	int count = 0;

	//Iterating through all force steps

	while (force < F_max) {
		if (count != 0) { force = force + (steps / t_tot) * F_max; }		//Force Increment
		
		count++;
		
		//cout << "Force = " << force << endl;

		F_ext(int(elements / 2), 0) = force;								//Assigning the node to apply force in Global External Force Matrix

		//Entering the Newton-Raphson Method

		for (int NR = 0; NR < 6; NR++) {

			//Calling Material Routine to get the Global Tangent Stiffness Matrix, Global Internal Force Matrix, Plastic Strain, Total Strain & Stress Values of All Elements

			tie(K_t, F_int, eps_p, strain, stress) = material_routine(elements, gauss_weights, Lengths, Areas, E_mat, u, steps, eta, eps_p, sigma_0, E);
			
			MatrixXd K_t_red = K_t.block(1, 1, K_t.cols() - 2, K_t.rows() - 2);		//Reducing the Tangent Stiffness
			MatrixXd K_t_red_inv = K_t_red.inverse();

			MatrixXd R = F_int - F_ext;
			//cout << "R=" << R << endl;
			MatrixXd R_red(elements - 1, 1),u_red(elements-1,1);					//Reducing the Residual 
			for (int c = 1; c < elements; c++) { R_red(c - 1, 0) = R(c, 0); }
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }

			MatrixXd del_u_red = K_t_red_inv * R_red;								//Computing Delta U
			
			u_red = u_red - del_u_red;
			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }

		}
		
		//Appending values of Center Node/First Element to the csv files

		//mydisplacefile << u((elements / 2) , 0) << endl;
		//mystressfile << stress((elements / 2)-1, 0) << endl;
		//mystrainfile << strain((elements / 2) - 1, 0) << endl;
		
		final_displacement_cpp_node(count, 0) = u((elements / 2), 0);
	}

	//mydisplacefile.close();
	//mystressfile.close();
	//mystrainfile.close();

	//cout << "Computation Done ! Exiting the Program" << endl;
	cout << "---FEA---" << endl;

	return final_displacement_cpp_node;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
