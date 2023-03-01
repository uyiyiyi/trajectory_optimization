#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include "OsqpEigen/OsqpEigen.h"

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
    static const int wp_num = 5;
    double time[wp_num];
    double max_velocity = 10;
    double total_time = 5;
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
    void solve_Nseg_bAp(int dim);
    void uniform_time_arrange(Eigen::MatrixXd path_);
    int solveQP(Eigen::MatrixXd Hessian);
};

min_snap::min_snap(/* args */)
{
    path.col(0) = Eigen::Vector2d(0, 0);
    path.col(1) = Eigen::Vector2d(1, 2);
    path.col(2) = Eigen::Vector2d(2, -1);
    path.col(3) = Eigen::Vector2d(4, 8);
    path.col(4) = Eigen::Vector2d(5, 2);
    uniform_time_arrange(path);
    solve_Nseg_bAp(0);
    solve_Nseg_bAp(1);

    // solve_bAp();
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

void min_snap::solve_Nseg_bAp(int dim)
{
    // p(t) = p0 + p1 * t + p2 * t2 + ... + pn * tn = ∑ pi * ti
    int derivative_order = 3; // pos = 0, vel = 1, acc = 2, jerk = 3, snap = 4;
    int seg_num = path.row(0).size() - 1;
    int p_order = 2 * derivative_order - 1; // Polynomial order, for jerk is 5, for snap is 7
    int p_num = p_order + 1;
    int m = 2 * wp_num - 1;
    int n = p_num * seg_num;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(p_num * seg_num, p_num * seg_num);
    Eigen::MatrixXd Qi = Eigen::MatrixXd::Zero(p_num, p_num);
    Eigen::VectorXd t = Eigen::VectorXd::Ones(p_num);
    std::cout << "0" << std::endl;

    for(int seg = 0; seg < seg_num; seg++)
    {
        for(int i = 0; i < p_num; i++)
        {
            t(i) = 1;
            for(int j = 0; j < derivative_order; j++)
            {
                int Ni = i - j;
                if(Ni < 0)
                    Ni = 0;
                t(i) *= Ni;
            }
        }
        Qi = t * t.transpose();
        Q.block(seg * p_num, seg * p_num, p_num, p_num) = Qi;
    } 
    // std::cout << t << std::endl;
    // std::cout << Q << std::endl;
    std::cout << "1" << std::endl;


    Eigen::SparseMatrix<double> hessian(n, n);      //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian = Q.sparseView();
    // std::cout << hessian << std::endl;
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(n);                    //Q: n*1向量
    Eigen::SparseMatrix<double> linearMatrix(m, n); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    Eigen::VectorXd lowerBound(m);                  //L: m*1下限向量
    Eigen::VectorXd upperBound(m);                  //U: m*1上限向量
    
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(m, n);
    Eigen::MatrixXd Ai = Eigen::MatrixXd::Zero(1, p_num);
    Eigen::MatrixXd Ai_ = Eigen::MatrixXd::Zero(1, p_num);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < p_order; j++)
        {
            int index = i / 2;
            Ai(0, j) = pow(time[index], j);
            Ai_(0, j) = pow(time[index + 1], j);
        }
        if(i % 2 == 0)
        {
            Aeq.block(i, p_num * i / 2, 1, p_num) = Ai;
            if(dim == 0)
            {
                lowerBound[i] = path.col(i / 2).x();
                upperBound[i] = path.col(i / 2).x();
            }
            else if (dim == 1)
            {
                lowerBound[i] = path.col(i / 2).y();
                upperBound[i] = path.col(i / 2).y();
            }
            else
            {
                lowerBound[i] = path.col(i / 2).z();
                upperBound[i] = path.col(i / 2).z();
            }
        }
        else
        {
            Aeq.block(i, p_num * (i - 1) / 2, 1, p_num) = Ai;
            Aeq.block(i, p_num * (i + 1) / 2, 1, p_num) = -Ai_;
            lowerBound[i] = 0;
            upperBound[i] = 0;
            // std::cout << Ai << std::endl;
            // std::cout << Ai_ << std::endl;
        }
    }
    std::cout << "2" << std::endl;
    std::cout << "m x n = " << m << " x " << n << std::endl;
    std::cout << Aeq;
    linearMatrix = Aeq.sparseView();
    std::cout << "2.1" << std::endl;
    std::cout << lowerBound << std::endl;
    
    // Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(4 * wp_num + 2, p_num * wp_num);
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, p_num);
    // // 初始条件以及结束条件的限制（包含有位置，速度，加速度等）
    // for(int i = 0; i < wp_num; i++)
    // {
    //     double t0 = time[i];
    //     for(int n = 0; n < p_order; n++)
    //     {
    //         A(0, n) = pow(time[i], n);
    //         if(n < 1)
    //             A(1, 0) = 0;
    //         else
    //             A(1, n) = n * pow(time[i], n - 1);
    //         if(n < 2)
    //             A(2, n) = 0;
    //         else
    //             A(2, n) = n * (n - 1) * pow(time[i], n - 2);
    //     }
    //     Aeq.block(3 * i, p_num * i, 3, p_num) = A;
    // }
    // std::cout << Aeq.sparseView() << std::endl;

    // 中间节点连续（包含有位置，速度，加速度）

    // instantiate the solver
    OsqpEigen::Solver solver;
    std::cout << "2.2" << std::endl;
    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver.data()->setNumberOfVariables(n);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);
    std::cout << "3" << std::endl;

    // instantiate the solver
    solver.initSolver();
    Eigen::VectorXd QPSolution;

    // solve the QP problem
    OsqpEigen::ErrorExitFlag err_flag = solver.solveProblem();
  
    // std::cout << "err_flag = " << err_flag << std::endl;
    QPSolution = solver.getSolution();
    std::cout << "QPSolution" << std::endl
              << QPSolution << std::endl; //输出为m*1的向量
    

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

void min_snap::uniform_time_arrange(Eigen::MatrixXd path_)
{
    int point_num = path_.row(0).size();
    // std::cout << path_ << std::endl;
    // std::cout << point_num << std::endl;
    int seg_num = point_num - 1;
    double dist[seg_num], total_dist = 0;
    time[0] = 0;
    for(int i = 0; i < seg_num; i++)
    {
        dist[i] = (path_.col(i + 1) - path_.col(i)).norm();
        total_dist += dist[i];
        // std::cout << "dist: " << dist[i] << std::endl;
    }
    for(int i = 0; i < point_num; i++)
    {
        time[i + 1] = time[i] + dist[i] / total_dist * total_time;
        
    }
    // for(int i = 0; i < point_num; i++)
    // {
    //     std::cout << "time: " << time[i] << std::endl;
    // }
}

int min_snap::solveQP(Eigen::MatrixXd Hessian)
{
    // allocate QP problem matrices and vectores
    
    return 0;
}

int main(int argc, const char** argv) {
    min_snap ms;
    return 0;
}