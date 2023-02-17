#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>

class min_snap
{
private:
    /* data */
public:
    std::vector<double> wp1 {0, 0};
    std::vector<double> wp2 {1, 2};
    std::vector<double> wp3 {2, -1};
    std::vector<double> wp4 {4, 8};
    std::vector<double> wp5 {5, 2};
    std::vector<std::vector<double>> path;
    
    // path.push_back(wp1);
    double v0[2] = {0, 0};
    double a0[2] = {0, 0};
    double v1[2] = {0, 0};
    double a1[2] = {0, 0};
    double T = 5;
    min_snap(/* args */);
    ~min_snap();
    void time_arrange(std::vector<std::vector<double>> path_);
};

min_snap::min_snap(/* args */)
{
    path.push_back(wp1);
    path.push_back(wp2);
    path.push_back(wp3);
    path.push_back(wp4);
    path.push_back(wp5);
    std::cout << " path: " << path.size() << std::endl;
    time_arrange(path);
}

min_snap::~min_snap()
{
}

void min_snap::time_arrange(std::vector<std::vector<double>> path_)
{
    for(int i; i < path.size(); i++)
    {
        
    }
}

int main(int argc, const char** argv) {
    // Eigen::Vector3d vec(1,2,3);
    // cout << vec << endl;
    min_snap ms;
    return 0;
}