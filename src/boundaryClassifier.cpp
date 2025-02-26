#include "boundaryClassifier.hpp"
#include <cmath>

using namespace std;

vector<Node*> identifyNeighbors(Domain &domain, Node* node) {
    long x = node->x;
    long y = node->y;
    long nx = domain.nX;
    long ny = domain.nY;

    vector<Node*> neighbors;
    neighbors.resize(9, nullptr);

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

    /* NOTE: neighbors[4] will always be a nullptr */
    if ((0 < x && x < nx - 1) && (0 < y && y < ny - 1)) {
        // Internal node
        neighbors[0] = domain.nodes[x+1][y+1]; // NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y-1]; // SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x-1][y+1]; // NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y-1]; // SW
    } else if ((x == 0) && (y == 0)) {
        // Bottom left corner of boundary
        neighbors[0] = domain.nodes[x+1][y+1]; // NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y+1]; // SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y+1];   // Down
        neighbors[6] = domain.nodes[x+1][y+1]; // NW
        neighbors[7] = domain.nodes[x+1][y];   // Left
        neighbors[8] = domain.nodes[x+1][y+1]; // SW
    } else if ((x == 0) && (y == ny - 1)) {
        // Top left corner
        neighbors[0] = domain.nodes[x+1][y-1]; // Flipped NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y-1]; // SE
        neighbors[3] = domain.nodes[x][y-1];   // Flipped up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x+1][y-1]; // Flipped NW
        neighbors[7] = domain.nodes[x+1][y];   // Flipped left
        neighbors[8] = domain.nodes[x+1][y-1]; // Flipped SW
    } else if ((x == nx - 1) && (y == 0)) {
        // Bottom right corner
        neighbors[0] = domain.nodes[x-1][y+1]; // Flipped NE
        neighbors[1] = domain.nodes[x-1][y];   // Flipped Right
        neighbors[2] = domain.nodes[x-1][y+1]; // Flipped SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y+1];   // Flipped down
        neighbors[6] = domain.nodes[x-1][y+1]; // NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y+1]; // Flipped SW
    } else if ((x == nx - 1) && (y == ny - 1)) {
        // Top right corner
        neighbors[0] = domain.nodes[x-1][y-1]; // Flipped NE
        neighbors[1] = domain.nodes[x-1][y];   // Flipped right
        neighbors[2] = domain.nodes[x-1][y-1]; // Flipped SE
        neighbors[3] = domain.nodes[x][y-1];   // Flipped up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x-1][y-1]; // Flipped NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y-1]; // SW
    } else if ((x == 0) && (0 < y && y < ny - 1)) {
        // Left edge
        neighbors[0] = domain.nodes[x+1][y+1]; // NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y-1]; // SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x+1][y+1]; // Flipped NW
        neighbors[7] = domain.nodes[x+1][y];   // Flipped left
        neighbors[8] = domain.nodes[x+1][y-1]; // Flipped SW
    } else if ((x == nx - 1) && (0 < y && y < ny - 1)) {
        // Right edge
        neighbors[0] = domain.nodes[x-1][y+1]; // Flipped NE
        neighbors[1] = domain.nodes[x-1][y];   // Flipped right
        neighbors[2] = domain.nodes[x-1][y-1]; // Flipped SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x-1][y+1]; // NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y-1]; // SW
    } else if ((0 < x && x < nx - 1) && (y == 0)) {
        // Bottom edge
        neighbors[0] = domain.nodes[x+1][y+1]; // NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y+1]; // Flipped SE
        neighbors[3] = domain.nodes[x][y+1];   // Up
        neighbors[5] = domain.nodes[x][y+1];   // Flipped down
        neighbors[6] = domain.nodes[x-1][y+1]; // NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y+1]; // Flipped SW
    } else if ((0 < x && x < nx - 1) && (y == ny - 1)) {
        // Top edge
        neighbors[0] = domain.nodes[x+1][y-1]; // Flipped NE
        neighbors[1] = domain.nodes[x+1][y];   // Right
        neighbors[2] = domain.nodes[x+1][y-1]; // SE
        neighbors[3] = domain.nodes[x][y-1];   // Flipped up
        neighbors[5] = domain.nodes[x][y-1];   // Down
        neighbors[6] = domain.nodes[x-1][y-1]; // Flipped NW
        neighbors[7] = domain.nodes[x-1][y];   // Left
        neighbors[8] = domain.nodes[x-1][y-1]; // SW
    }
    return neighbors;
}

void identifyBoundaryNodes(Domain &domain) {
    Node* node;
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            if (node->id == 0) {
                vector<Node*> neighbors = identifyNeighbors(domain, node);

                bool ne    = neighbors[0]->id > 0;
                bool right = neighbors[1]->id > 0;
                bool se    = neighbors[2]->id > 0;
                bool above = neighbors[3]->id > 0;
                bool below = neighbors[5]->id > 0;
                bool nw    = neighbors[6]->id > 0;
                bool left  = neighbors[7]->id > 0;
                bool sw    = neighbors[8]->id > 0;

                int numNeighbors = ne + right + se + above + below + nw + left + sw;


                if (numNeighbors > 0) {

                    if (numNeighbors > 3) node->isConcave = true;

                    // Electrolyte to the right
                    if (right && !(above && left && below && nw && ne && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte below
                    if (above && !(left && right && below && nw && ne && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte to the left
                    if (left && !(right && above && below && nw && ne && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte above
                    if (below && !(right && above && left && nw && ne && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte above right
                    if (nw && !(right && above && left && below && ne && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte above left
                    if (ne && !(right && above && left && below && nw && se && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte below left
                    if (se && !(right && above && left && below && nw && ne && sw)) {
                        node->isBoundary = true;
                    }
                    
                    // Electrolyte below right
                    if (sw && !(right && above && left && below && nw && ne && se)) {
                        node->isBoundary = true;
                    }
                }
            }
        }
    }
}

double euclideanDistance(Node* node, pair<int, int> neighbor) {
    double x1 = node->x;
    double x2 = neighbor.first;
    double y1 = node->y;
    double y2 = neighbor.second;

    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

vector<pair<int, int>> kClosestNeighbors(Node* node, vector<pair<int, int>> edgeNeighbors, int k) {
    if (edgeNeighbors.size() <= k) {
        return edgeNeighbors;  // If fewer than or exactly k neighbors, return all available
    }

    // Sort neighbors based on Euclidean distance
    sort(edgeNeighbors.begin(), edgeNeighbors.end(), [&](const pair<int, int>& a, const pair<int, int>& b) {
        return euclideanDistance(node, a) < euclideanDistance(node, b);
    });

    // Return the first k closest neighbors
    return vector<pair<int, int>>(edgeNeighbors.begin(), edgeNeighbors.begin() + k);
}

void classifyBC(Domain &domain) {
    Node* node;
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            if (node->isBoundary) {
                vector<Node*> neighbors = identifyNeighbors(domain, node);
                vector<pair<int, int>> edgeNeighbors;

                for (int k=0; k < neighbors.size(); k++) {
                    if (neighbors[k] != nullptr && neighbors[k]->isBoundary) {
                        edgeNeighbors.emplace_back(i+e[k][0], j+e[k][1]);
                    }
                }

                edgeNeighbors = kClosestNeighbors(node, edgeNeighbors, 2);

                // Calculate slope
                double dx = edgeNeighbors[0].first  - edgeNeighbors[1].first;
                double dy = edgeNeighbors[0].second - edgeNeighbors[1].second;
                double normal;

                if (dx == 0) {
                    normal = 0;
                } else if (dy == 0) {
                    normal = numeric_limits<double>::infinity();
                } else {
                    normal = -(dx/dy);
                }

                /* 
                Node ID Key:
                ┌────---┬───────────────────────────┬───────────┐
                │ ID    │ Description               │ Normal ID │
                ├────---┼───────────────────────────┼───────────┤
                │ 0     │ General fluid node        │           │
                │ 20    │ Solid node                │           │
                │ 21-29 │ Fluid node next to solid  │           │
                ├────---┼───────────────────────────┼───────────┤
                │ 21    │ NE                        │ -5        │
                │ 22    │ E                         │ -1        │
                │ 23    │ SE                        │ -7        │
                │ 24    │ N                         │ -3        │
                │ 25    │ Middle                    │           │
                │ 26    │ S                         │ -4        │
                │ 27    │ NW                        │ -8        │
                │ 28    │ W                         │ -2        │
                │ 29    │ SW                        │ -6        │
                └────---┴───────────────────────────┴───────────┘
                */
                if (normal == 0) {
                    if (domain.nodes[i+1][j]->id == 0 || domain.nodes[i+1][j]->isBoundary) {
                        node->normalNodeID = -1;
                        node->id = 22;
                    } else {
                        node->normalNodeID = -2;
                        node->id = 28;
                    }

                } else if (isinf(normal)) {
                    if (domain.nodes[i][j+1]->id == 0 || domain.nodes[i][j+1]->isBoundary) {
                        if (domain.nodes[i][j+1]->id == 0) {
                            node->normalNodeID = -3;
                            node->id = 24;
                        } else if (domain.nodes[i][j-1]->id != 1) {
                            node->normalNodeID = -4;
                            node->id = 26;
                        } else {
                            node->normalNodeID = -3;
                            node->id = 24;
                        }
                    } else {
                        node->normalNodeID = -4;
                        node->id = 26;
                    }

                } else if (normal == 1) {
                    if (domain.nodes[i+1][j+1]->id == 0 || domain.nodes[i+1][j+1]->isBoundary) {
                        node->normalNodeID = -5;
                        node->id = 21;
                    } else {
                        node->normalNodeID = -6;
                        node->id = 29;
                    }
                
                } else if (normal == -1) {
                    if (domain.nodes[i+1][j-1]->id == 0 || domain.nodes[i+1][j-1]->isBoundary) {
                        node->normalNodeID = -7;
                        node->id = 23;
                    } else {
                        node->normalNodeID = -8;
                        node->id = 27;
                    }
                }
            }
        }
    }
}

void identifyInletNodes(Domain &domain) {
    for (long j = 1; j < domain.nY - 1; j++) {
        domain.nodes[0][j]->isInlet = true;
        domain.nodes[0][j]->id = 22;
        domain.nodes[0][j]->normalNodeID = -1;
    }
    domain.nodes[0][1]->normalNodeID = -5;
    domain.nodes[0][domain.nY - 2]->normalNodeID = -7;
}

void identifyOutletNodes(Domain &domain) {
    for (long j = 1; j < domain.nY - 1; j++) {
        domain.nodes[domain.nX - 1][j]->isOutlet = true;
        domain.nodes[domain.nX - 1][j]->id = 28;
        domain.nodes[domain.nX - 1][j]->normalNodeID = -2;
    }
    domain.nodes[domain.nX - 1][1]->normalNodeID = -8;
    domain.nodes[domain.nX - 1][domain.nY - 2]->normalNodeID = -6;
}
