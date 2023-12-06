#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unordered_set>
#include <string>
#include <math.h>

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 0;
    T damping_coeff = 0; 
    T dt = 0;

    // node data
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    // segment data
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> rest_length;

    if (argc < 2) 
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "0") == 0) // cloth case
    {
        // 1. Create node data: position, mass, velocity
        T cloth_total_mass = 2;
        int cloth_res = 64;
        int num_nodes = std::pow(cloth_res, 2);
        for (int i = 0; i < cloth_res; i++) {
            for (int j = 0; j < cloth_res; j++) {
                TV pos = TV::Zero();
                pos(0) = i / 32.f;
                pos(2) = j / 32.f;
                x.push_back(pos);
                TV vel = TV::Zero();
                v.push_back(vel);
                m.push_back(cloth_total_mass / num_nodes);
                node_is_fixed.push_back(false);
            }
        }
        
        // 2. Fill segments and rest_length, including struct springs, shearing springs and bending springs.
        int curr_idx = 0;
        for (int i = 0; i < cloth_res; i++) {
            for (int j = 0; j < cloth_res; j++) {
                // Add structural springs
                if (j != cloth_res-1) {
                    Eigen::Matrix<int,2,1> horiz_seg = {curr_idx, curr_idx+1};
                    segments.push_back(horiz_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx][0] - x[curr_idx+1][0], 2) +
                                                    std::pow(x[curr_idx][1] - x[curr_idx+1][1], 2) +
                                                    std::pow(x[curr_idx][2] - x[curr_idx+1][2], 2)));
                }
                if (i != cloth_res-1) {
                    Eigen::Matrix<int,2,1> vert_seg = {curr_idx, curr_idx+cloth_res};
                    segments.push_back(vert_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx][0] - x[curr_idx+cloth_res][0], 2) +
                                                    std::pow(x[curr_idx][1] - x[curr_idx+cloth_res][1], 2) +
                                                    std::pow(x[curr_idx][2] - x[curr_idx+cloth_res][2], 2)));
                }
                // Add shearing springs
                if (i != cloth_res-1 && j != cloth_res-1) {
                    Eigen::Matrix<int,2,1> diag1_seg = {curr_idx, curr_idx+cloth_res+1};
                    segments.push_back(diag1_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx][0] - x[curr_idx+cloth_res+1][0], 2) +
                                                    std::pow(x[curr_idx][1] - x[curr_idx+cloth_res+1][1], 2) +
                                                    std::pow(x[curr_idx][2] - x[curr_idx+cloth_res+1][2], 2)));
                    Eigen::Matrix<int,2,1> diag2_seg = {curr_idx+1, curr_idx+cloth_res};
                    segments.push_back(diag2_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx+1][0] - x[curr_idx+cloth_res][0], 2) +
                                                    std::pow(x[curr_idx+1][1] - x[curr_idx+cloth_res][1], 2) +
                                                    std::pow(x[curr_idx+1][2] - x[curr_idx+cloth_res][2], 2)));
                }
                // Add bending springs
                if (j < cloth_res-2) {
                    Eigen::Matrix<int,2,1> horiz_seg = {curr_idx, curr_idx+2};
                    segments.push_back(horiz_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx][0] - x[curr_idx+2][0], 2) +
                                                    std::pow(x[curr_idx][1] - x[curr_idx+2][1], 2) +
                                                    std::pow(x[curr_idx][2] - x[curr_idx+2][2], 2)));
                }
                if (i < cloth_res-2) {
                    Eigen::Matrix<int,2,1> vert_seg = {curr_idx, curr_idx+(2*cloth_res)};
                    segments.push_back(vert_seg);
                    rest_length.push_back(std::sqrt(std::pow(x[curr_idx][0] - x[curr_idx+(2*cloth_res)][0], 2) +
                                                    std::pow(x[curr_idx][1] - x[curr_idx+(2*cloth_res)][1], 2) +
                                                    std::pow(x[curr_idx][2] - x[curr_idx+(2*cloth_res)][2], 2)));
                }
                curr_idx++;
            }
        }
        
        // 3. Choose proper youngs_modulus, damping_coeff and dt.
        youngs_modulus = 5.3;
        damping_coeff = 0.4;
        dt = 0.0001;
        
        // 4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        node_is_fixed[cloth_res * (cloth_res-1)] = true;
        node_is_fixed[0] = true;
        
        // Move fixed nodes
        driver.helper = [&](T t, T dt) {
            if (t < 1.5) {
                TV v_displacement = TV::Zero();
                v_displacement(0) = 0;
                v_displacement(1) = 0;
                v_displacement(2) = 1;
                driver.ms.v[cloth_res * (cloth_res-1)] = v_displacement * t;
                driver.ms.v[0] = v_displacement * t;
                driver.ms.x[cloth_res * (cloth_res-1)] += driver.ms.v[cloth_res * (cloth_res-1)] * dt;
                driver.ms.x[0] += driver.ms.v[0] * dt;
            } else if (t < 2.25) {
                TV v_displacement = TV::Zero();
                v_displacement(0) = 0;
                v_displacement(1) = 0;
                v_displacement(2) = -1;
                driver.ms.v[cloth_res * (cloth_res-1)] = v_displacement * t;
                driver.ms.v[0] = v_displacement * t;
                driver.ms.x[cloth_res * (cloth_res-1)] += driver.ms.v[cloth_res * (cloth_res-1)] * dt;
                driver.ms.x[0] += driver.ms.v[0] * dt;
            } else if (t < 3.5) {
                // top right corner
                TV v_displacement_left = TV::Zero();
                v_displacement_left(0) = 0.25;
                v_displacement_left(1) = 0;
                v_displacement_left(2) = 0.25;
                driver.ms.v[cloth_res * (cloth_res-1)] = v_displacement_left * t;
                driver.ms.x[cloth_res * (cloth_res-1)] += driver.ms.v[cloth_res * (cloth_res-1)] * dt;
                // top left corner
                TV v_displacement_right = TV::Zero();
                v_displacement_right(0) = -0.25;
                v_displacement_right(1) = 0;
                v_displacement_right(2) = -0.25;
                driver.ms.v[0] = v_displacement_right * t;
                driver.ms.x[0] += driver.ms.v[0] * dt;
            } else if (t < 4.5) {
                // top right corner
                TV v_displacement_left = TV::Zero();
                v_displacement_left(0) = -0.25;
                v_displacement_left(1) = 0;
                v_displacement_left(2) = -0.25;
                driver.ms.v[cloth_res * (cloth_res-1)] = v_displacement_left * t;
                driver.ms.x[cloth_res * (cloth_res-1)] += driver.ms.v[cloth_res * (cloth_res-1)] * dt;
                // top left corner
                TV v_displacement_right = TV::Zero();
                v_displacement_right(0) = 0.25;
                v_displacement_right(1) = 0;
                v_displacement_right(2) = 0.25;
                driver.ms.v[0] = v_displacement_right * t;
                driver.ms.x[0] += driver.ms.v[0] * dt;
            }
        };
        driver.test = "cloth";
        
        // 5. Generate quad mesh for rendering.
        std::string filename = "data/" + driver.test + ".obj";
        std::ofstream cloth_quad_mesh;
        cloth_quad_mesh.open(filename);
        // Add vertices
        for (int i = 0; i < num_nodes; i++) {
            cloth_quad_mesh << "v";
            for (int j = 0; j < dim; j++) {
                cloth_quad_mesh << " " << x[i][j];
            }
            cloth_quad_mesh << "\n";
        }
        // Add faces
        curr_idx = 1;
        for (int i = 0; i < cloth_res-1; i++) {
            for (int j = 0; j < cloth_res-1; j++) {
                cloth_quad_mesh << "f";
                cloth_quad_mesh << " " << curr_idx;
                cloth_quad_mesh << " " << curr_idx+cloth_res;
                cloth_quad_mesh << " " << curr_idx+cloth_res+1;
                cloth_quad_mesh << " " << curr_idx+1 << "\n";
                curr_idx++;
            }
            curr_idx++;
        }
        cloth_quad_mesh.close();
    }

    else if (strcmp(argv[1], "1") == 0) // volumetric bunny case
    {
        // 1. Create node data from data/points: The first line indicates the number of points and dimension (which is 3).
        T bunny_total_mass = 18;
        // Set up file reader
        std::ifstream points_file;
        points_file.open("data/points");
        if (!points_file) { // Check if file opened properly
            std::cout << "Unable to open file" << std::endl;
            exit(1); // terminate with error
        }
        // Extract number of nodes and dimension
        int num_nodes = 0;
        std::string points_line;
        if (getline(points_file, points_line)) {
            std::istringstream in(points_line); // Make a stream for points_line
            int n, d;
            in >> n >> d;
            num_nodes = n;
        }
        // Extract position data one line at a time
        while (getline(points_file, points_line)) {
            std::istringstream in(points_line); // Make a stream for points_line
            T x_coord, y_coord, z_coord;
            in >> x_coord >> y_coord >> z_coord;
            TV pos = TV::Zero();
            pos(0) = x_coord;
            pos(1) = y_coord;
            pos(2) = z_coord;
            x.push_back(pos);
            TV vel = TV::Zero();
            v.push_back(vel);
            m.push_back(bunny_total_mass / num_nodes);
            node_is_fixed.push_back(false);
        }
        points_file.close();
        
        // 2. Fill segments and rest_length from data/cells: The first line indicates the number of tetrahedra and the number of vertices of each tet (which is 4). Each edge in this tetrahedral mesh will be a segment. Be careful not to create duplicate edges.
        // Set up file reader
        std::ifstream cells_file;
        cells_file.open("data/cells");
        if (!cells_file) { // Check if file opened properly
            std::cout << "Unable to open file" << std::endl;
            exit(1); // terminate with error
        }
        // Extract number of tetrahedra and dimension
        int num_tetrahedra = 0;
        std::string cells_line;
        if (getline(cells_file, cells_line)) {
            std::istringstream in(cells_line); // Make a stream for cells_line
            int n, d;
            in >> n >> d;
            num_tetrahedra = n;
        }
        // Extract tetrahedra data one line at a time
        std::unordered_set<std::string> duplicate = {};
        while (getline(cells_file, cells_line)) {
            std::istringstream in(cells_line); // Make a stream for cells_line
            T w_coord, x_coord, y_coord, z_coord;
            in >> w_coord >> x_coord >> y_coord >> z_coord;
            std::vector<T> tetrahedron_indices;
            tetrahedron_indices.push_back(w_coord);
            tetrahedron_indices.push_back(x_coord);
            tetrahedron_indices.push_back(y_coord);
            tetrahedron_indices.push_back(z_coord);
            for (int i = 0; i < 4; i++) {
                for (int j = i+1; j < 4; j++) {
                    // Generate segment label and combination of two indices, with smaller index first, to check for duplicate segments
                    std::string seg_label = std::to_string(tetrahedron_indices[i]) + std::to_string(tetrahedron_indices[j]);
                    if (tetrahedron_indices[j] < tetrahedron_indices[i]) {
                        seg_label = std::to_string(tetrahedron_indices[j]) + std::to_string(tetrahedron_indices[i]);
                    }
                    if (duplicate.find(seg_label) == duplicate.end()) { // Segment not previously found
                        duplicate.insert(seg_label);
                        Eigen::Matrix<int,2,1> seg = {tetrahedron_indices[i], tetrahedron_indices[j]};
                        segments.push_back(seg);
                        rest_length.push_back(std::sqrt(std::pow(x[tetrahedron_indices[i]][0] - x[tetrahedron_indices[j]][0], 2) +
                                                        std::pow(x[tetrahedron_indices[i]][1] - x[tetrahedron_indices[j]][1], 2) +
                                                        std::pow(x[tetrahedron_indices[i]][2] - x[tetrahedron_indices[j]][2], 2)));
                    }
                }
            }
        }
        cells_file.close();
        
        // 3. Choose proper youngs_modulus, damping_coeff, dt.
        youngs_modulus = 7.5;
        damping_coeff = 2.5;
        dt = 0.0001;
        
        // 4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        node_is_fixed[2140] = true; // ear nodes
        node_is_fixed[2346] = true;
        node_is_fixed[1036] = true; // tail node
        
        // Move fixed nodes
        driver.helper = [&](T t, T dt) {
            if (t < 1) {
                TV v_displacement = TV::Zero();
                v_displacement(0) = 1;
                v_displacement(1) = 0.1;
                v_displacement(2) = 0.1;
                driver.ms.v[1036] = v_displacement * t;
                driver.ms.x[1036] += driver.ms.v[1036] * dt;
            } else {
                driver.ms.node_is_fixed[1036] = false;
            }
        };
        driver.test = "bunny";
    }

    else {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    // simulate
    driver.dt = dt;
    driver.ms.segments = segments;
    driver.ms.m = m;
    driver.ms.v = v;
    driver.ms.x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;
    
    driver.run(120);

    return 0;
}
