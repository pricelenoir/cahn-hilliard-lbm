#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include "json.hpp"
#include "node.hpp"
#include "constants.hpp"

// Directional offsets for neighbors
/* 
    # D2Q9 Lattice structure
    #   6         3         0
    #      .      ^     .
    #         .   |   .
    #           . | . 
    #   7 <------ 4 ------> 1
    #           . |  .
    #        .    |     .
    #     .       v       .
    #   8         5         2
*/

static int e[9][2] = {
    { 1,  1}, { 1,  0}, { 1, -1},
    { 0,  1}, { 0,  0}, { 0, -1},
    {-1,  1}, {-1,  0}, {-1, -1}
};

// Weights
static double w[9] = {1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 4.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36};

class Domain {
public:
    long nX;                        // Number of nodes in X direction
    long nY;                        // Number of nodes in Y direction
    int periodicity[2];             // Periodicity in X and Y directions

    std::vector<double> resUVec;    // Residual velocity iteration series
    std::vector<double> resPhiVec;  // Residual order parameter iteration series
    std::vector<double> resPVec;    // Residual pressure iteration series

    std::vector<std::vector <Node *>> nodes;

    Domain(long nX, long nY);
    void initialize(const nlohmann::json& config, Constants &constants);
    void save(const nlohmann::json& config, int iter);
};

#endif // DOMAIN_HPP