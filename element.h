using namespace std;
using namespace Eigen;

MatrixXd B_mat(double le);
int signum(double x);
MatrixXd assign(int ele_index_cpp, int numele);
tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double a_1, double a_2, double len);
