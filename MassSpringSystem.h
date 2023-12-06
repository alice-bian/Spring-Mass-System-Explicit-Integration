#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f)
    {
        // Allocate and clear space in f
        f.clear();
        int num_nodes = m.size();
        f.resize(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            f[i] = TV::Zero();
        }
        
        // 3D convergent spring force
        int num_segs = segments.size();
        for (int i = 0; i < num_segs; i++) {
            Eigen::Matrix<int,2,1> segment = segments[i];
            TV x_1 = TV::Zero();
            x_1 = x[segment[0]];
            TV x_2 = TV::Zero();
            x_2 = x[segment[1]];
            T l = std::sqrt(std::pow(x_1[0] - x_2[0], 2) + std::pow(x_1[1] - x_2[1], 2) + std::pow(x_1[2] - x_2[2], 2));
            T l_0 = rest_length[i];
            TV n_12 = TV::Zero();
            n_12 = (x_1 - x_2) / l;
            TV f_1 = TV::Zero();
            f_1 = (-1) * youngs_modulus * (l / l_0 - 1) * n_12;
            f[segment[0]] += f_1;
            f[segment[1]] -= f_1;
        }
    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        // Allocate and clear space in f
        f.clear();
        int num_nodes = m.size();
        f.resize(num_nodes);
        for (int i = 0; i < num_nodes; i++) {
            f[i] = TV::Zero();
        }
        
        // Spring dashpot damping force
        int num_segs = segments.size();
        for (int i = 0; i < num_segs; i++) {
            Eigen::Matrix<int,2,1> segment = segments[i];
            TV x_1 = TV::Zero();
            x_1 = x[segment[0]];
            TV x_2 = TV::Zero();
            x_2 = x[segment[1]];
            T l = std::sqrt(std::pow(x_1[0] - x_2[0], 2) + std::pow(x_1[1] - x_2[1], 2) + std::pow(x_1[2] - x_2[2], 2));
            TV n_12 = TV::Zero();
            n_12 = (x_1 - x_2) / l;
            TV v_1 = TV::Zero();
            v_1 = v[segment[0]];
            TV v_2 = TV::Zero();
            v_2 = v[segment[1]];
            T v_rel = (v_1 - v_2)[0] * n_12[0] + (v_1 - v_2)[1] * n_12[1] + (v_1 - v_2)[2] * n_12[2];
            TV f_1 = TV::Zero();
            f_1 = (-1) * damping_coeff * v_rel * n_12;
            f[segment[0]] += f_1;
            f[segment[1]] -= f_1;
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
};
