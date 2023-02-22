#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include<Eigen/Core> 
#include<Eigen/SVD>   

template<typename _Matrix_Type_> 
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = 
    std::numeric_limits<double>::epsilon()) 
{  
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);  
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);  
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint(); 
}   

class min_snap
{
private:
    /* data */
public:
    int dim = 2; // if 2: x, y; if 3: x, y, z
    int wp_num = 5;
    Eigen::MatrixXd path = Eigen::MatrixXd::Zero(dim, wp_num);
    
    // path.push_back(wp1);
    double v0[2] = {0, 0};
    double a0[2] = {0, 0};
    double v1[2] = {0, 0};
    double a1[2] = {0, 0};
    double T = 5;
    min_snap(/* args */);
    ~min_snap();
    void solve_bAp();
    void solve_Nseg_bAp();
    void time_arrange(std::vector<std::vector<double>> path_);
};

min_snap::min_snap(/* args */)
{
    path.col(0) = Eigen::Vector2d(0, 0);
    path.col(1) = Eigen::Vector2d(1, 2);
    path.col(2) = Eigen::Vector2d(2, -1);
    path.col(3) = Eigen::Vector2d(4, 8);
    path.col(4) = Eigen::Vector2d(5, 2);
    solve_bAp();
    // std::cout << " path: " << std::endl << path << std::endl;
    // path.push_back(wp1);
    // path.push_back(wp2);
    // path.push_back(wp3);
    // path.push_back(wp4);
    // path.push_back(wp5);
    // std::cout << " path: " << path.size() << std::endl;
    // time_arrange(path);
}

min_snap::~min_snap()
{
}

void min_snap::solve_Nseg_bAp()
{
    /* x(t) = pi * t^i + ... + p5 * t^5 + p4 * t^4 + p3 * t^3 + p2 * t^2 + p1 * t + p0
    for 5th order, i = 5, size of vector 7 is i + 1
    x = 0, 1, 2, 4, 5, so segments N = 4 */
    Eigen::VectorXd x;
    x << 0, 1, 2, 4, 5;
    int N = x.size() - 1;
    int order = 6;
    double t0 = 0;
    double t1 = 2;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(6); // Boundary condition: [x0, x1, v0, v1, a0, a1]
    b(1) = 2;
    b(3) = 1;
    Eigen::VectorXd p(order + 1); // for order = 6, p(7)
    Eigen::MatrixXd A(6, order + 1); // dimension is: (Boundary condition size: 6) * (order + 1))
    
    A.row(0) << pow(t0, 5), pow(t0, 4), pow(t0, 3), pow(t0, 2), t0, 1;
    A.row(1) << pow(t1, 5), pow(t1, 4), pow(t1, 3), pow(t1, 2), t1, 1;
    A.row(2) << 0, 0, 0, 0, 1, 0;
    A.row(3) << 5 * pow(t1, 4), 4 * pow(t1, 3), 3 * pow(t1, 2), 2 * t1, 1, 0;
    A.row(4) << 0, 0, 0, 2, 0, 0;
    A.row(5) << 20 * pow(t1, 3), 12 * pow(t1, 2), 6 * t1, 2, 0, 0;
    Eigen::FullPivLU<Eigen::MatrixXd> luA(A);
    int rank = luA.rank();
    std::cout << "Rank(A) = " << rank << std::endl;
    std::cout << "A^-1 = " << std::endl << A.inverse() << std::endl;
    p = A.inverse() * b;
    std::cout << "p = [" << p.transpose() << ']' << std::endl;
}

void min_snap::solve_bAp()
{
    /* x(t) = pi * t^i + ... + p5 * t^5 + p4 * t^4 + p3 * t^3 + p2 * t^2 + p1 * t + p0
    for 5th order, i = 5, size of vector to be determined p is i + 1 */
    int order = 6;
    double t0 = 0;
    double t1 = 2;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(6); // Boundary condition: [x0, x1, v0, v1, a0, a1]
    b(1) = 2;
    b(3) = 1;
    Eigen::VectorXd p(order + 1); // for order = 6, p(7)
    Eigen::MatrixXd A(6, order + 1); // dimension is: (Boundary condition size: 6) * (order + 1))
    
    A.row(0) << pow(t0, 6), pow(t0, 5), pow(t0, 4), pow(t0, 3), pow(t0, 2), t0, 1;
    A.row(1) << pow(t1, 6), pow(t1, 5), pow(t1, 4), pow(t1, 3), pow(t1, 2), t1, 1;
    A.row(2) << 6 * pow(t0, 5), 5 * pow(t0, 4), 4 * pow(t0, 3), 3 * pow(t0, 2), 2 * t0, 1, 0;
    A.row(3) << 6 * pow(t1, 5), 5 * pow(t1, 4), 4 * pow(t1, 3), 3 * pow(t1, 2), 2 * t1, 1, 0;
    A.row(4) << 30 * pow(t0, 4), 20 * pow(t0, 3), 12 * pow(t0, 2), 6 * t0, 2, 0, 0;
    A.row(5) << 30 * pow(t1, 4), 20 * pow(t1, 3), 12 * pow(t1, 2), 6 * t1, 2, 0, 0;
    Eigen::FullPivLU<Eigen::MatrixXd> luA(A);
    int rank = luA.rank();
    std::cout << "Rank(A) = " << rank << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    // std::cout << "A^-1 = " << std::endl << A.inverse() << std::endl;
    std::cout << "A^-1 = " << std::endl << pseudoInverse(A) << std::endl;
    // p = A.inverse() * b;
    p = pseudoInverse(A) * b;
    std::cout << "p = [" << p.transpose() << ']' << std::endl;
    /* The number of equations is less than the number of unknowns,
    so there are infinite solutions. 
    The solution here is the minimum norm solution. i.e. min||p|| */
    double x1 = A.row(1) * p;
    std::cout << "x1 = " << x1 << std::endl;
    double v1 = A.row(3) * p;
    std::cout << "v1 = " << v1 << std::endl;
}

void min_snap::time_arrange(std::vector<std::vector<double>> path_)
{
    for(int i; i < path_.size(); i++)
    {

    }
}

int main(int argc, const char** argv) {
    // Eigen::Vector3d vec(1,2,3);
    // cout << vec << endl;
    min_snap ms;
    return 0;
}