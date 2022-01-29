using namespace std;
using namespace Eigen;

double Lambda(double sigma, double sigma_0);
tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>material_routine(int elements, double gauss_weights, MatrixXd Lengths, MatrixXd Areas, MatrixXd  E_mat, MatrixXd u, double time_step, double eta, MatrixXd eps_p, double sigma_0, double E);
tuple<double, double, double>get_Stress_Strains(double eps_p_m, double eps, double eta, double h, double t_total, double sigma_0, double E);

