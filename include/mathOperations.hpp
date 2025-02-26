#ifndef MATH_OPERATIONS_HPP
#define MATH_OPERATIONS_HPP

#include "domain.hpp"
#include "node.hpp"

using namespace std;

template <typename T>
void derivativeX(Node* node, Domain &domain, T Node::*value, T Node::*derivative) {
    if (node->id == 20) return; // Solid wall

    long ix = node->x;
    long iy = node->y;
    long nx = domain.nX;
    long ny = domain.nY;

    // Calculate neighboring indices
    long ixp = ix + 1;
    long ixn = ix - 1;
    long iyp = iy + 1;
    long iyn = iy - 1;
    long ixpp = ixp + 1;
    long ixnn = ixn - 1;
    long iypp = iyp + 1;
    long iynn = iyn - 1;

    if (ixp > nx - 1) {
        ixp = ix - 1;
    } else if (ixp < 0) {
        ixp = ix + 1;
    }

    if (iyp > ny - 1) {
        iyp = iy - 1;
    } else if (iyp < 0) {
        iyp = iy + 1;
    }

    if (ixn > nx - 1) {
        ixn = ix - 1;
    } else if (ixn < 0) {
        ixn = ix + 1;
    }

    if (iyn > ny - 1) {
        iyn = iy - 1;
    } else if (iyn < 0) {
        iyn = iy + 1;
    }

    if (ixpp > nx- 1) {
        ixpp = (ixpp == nx) ? ix : ix - 2;
    } else if (ixpp < 0) {
        ixpp = (ixpp == -1) ? ix : ix + 2;
    }

    if (iypp > ny - 1) {
        iypp = (iypp == ny) ? iy : iy - 2;
    } else if (iypp < 0) {
        iypp = (iypp == -1) ? iy : iy + 2;
    }

    if (ixnn > nx - 1) {
        ixnn = (ixnn == nx) ? ix : ix - 2;
    } else if (ixnn < 0) {
        ixnn = (ixnn == -1) ? ix : ix + 2;
    }

    if (iynn > ny - 1) {
        iynn = (iynn == ny) ? iy : iy - 2;
    } else if (iynn < 0) {
        iynn = (iynn == -1) ? iy : iy + 2;
    }

    switch(node->id) {
        // Fluid nodes: central difference
        case 0:
            node->*derivative = (domain.nodes[ixp][iy]->*value - domain.nodes[ixn][iy]->*value)/3 + 
                                (domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value)/12 + 
                                (domain.nodes[ixp][iyn]->*value - domain.nodes[ixn][iyp]->*value)/12;
            break;

        // East boundary nodes: forward difference
        case 21:
        case 22:
        case 23:
            node->*derivative = 0.5*(-(domain.nodes[ixpp][iy]->*value) + 4*domain.nodes[ixp][iy]->*value - 3*domain.nodes[ix][iy]->*value);
            break;

        // West boundary node: backward difference
        case 27:
        case 28:
        case 29:
            node->*derivative = -0.5*(-(domain.nodes[ixnn][iy]->*value) + 4*domain.nodes[ixn][iy]->*value - 3*domain.nodes[ix][iy]->*value);
            break;

        // Other boundary nodes: regular central difference
        default:
            node->*derivative = 0.5*(domain.nodes[ixp][iy]->*value - domain.nodes[ixn][iy]->*value);
            break;
    }
}

template <typename T>
void derivativeY(Node* node, Domain &domain, T Node::*value, T Node::*derivative) {
    if (node->id == 20) return; // Solid wall

    long ix = node->x;
    long iy = node->y;
    long nx = domain.nX;
    long ny = domain.nY;

    // Calculate neighboring indices
    long ixp = ix + 1;
    long ixn = ix - 1;
    long iyp = iy + 1;
    long iyn = iy - 1;
    long ixpp = ixp + 1;
    long ixnn = ixn - 1;
    long iypp = iyp + 1;
    long iynn = iyn - 1;

    if (ixp > nx - 1) {
        ixp = ix - 1;
    } else if (ixp < 0) {
        ixp = ix + 1;
    }

    if (iyp > ny - 1) {
        iyp = iy - 1;
    } else if (iyp < 0) {
        iyp = iy + 1;
    }

    if (ixn > nx - 1) {
        ixn = ix - 1;
    } else if (ixn < 0) {
        ixn = ix + 1;
    }

    if (iyn > ny - 1) {
        iyn = iy - 1;
    } else if (iyn < 0) {
        iyn = iy + 1;
    }

    if (ixpp > nx- 1) {
        ixpp = (ixpp == nx) ? ix : ix - 2;
    } else if (ixpp < 0) {
        ixpp = (ixpp == -1) ? ix : ix + 2;
    }

    if (iypp > ny - 1) {
        iypp = (iypp == ny) ? iy : iy - 2;
    } else if (iypp < 0) {
        iypp = (iypp == -1) ? iy : iy + 2;
    }

    if (ixnn > nx - 1) {
        ixnn = (ixnn == nx) ? ix : ix - 2;
    } else if (ixnn < 0) {
        ixnn = (ixnn == -1) ? ix : ix + 2;
    }

    if (iynn > ny - 1) {
        iynn = (iynn == ny) ? iy : iy - 2;
    } else if (iynn < 0) {
        iynn = (iynn == -1) ? iy : iy + 2;
    }

    switch(node->id) {
        // Fluid nodes: central difference
        case 0:
            node->*derivative = (domain.nodes[ix][iyp]->*value - domain.nodes[ix][iyn]->*value)/3 + 
                                (domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value)/12 + 
                                (domain.nodes[ixn][iyp]->*value - domain.nodes[ixp][iyn]->*value)/12;
            break;

        // East boundary nodes: forward difference
        case 21:
        case 24:
        case 27:
            node->*derivative = 0.5*(-(domain.nodes[ix][iypp]->*value) + 4*domain.nodes[ix][iyp]->*value - 3*domain.nodes[ix][iy]->*value);
            break;

        // West boundary node: backward difference
        case 23:
        case 26:
        case 29:
            node->*derivative = -0.5*(-(domain.nodes[ix][iynn]->*value) + 4*domain.nodes[ix][iyn]->*value - 3*domain.nodes[ix][iy]->*value);
            break;

        // Other boundary nodes: regular central difference
        default:
            node->*derivative = 0.5*(domain.nodes[ix][iyp]->*value - domain.nodes[ix][iyn]->*value);
            break;
    }
}

template <typename T>
void laplace(Node* node, Domain &domain, T Node::*value, T Node::*solution) {
    if (node->id == 20) return; // Solid wall

    long ix = node->x;
    long iy = node->y;
    long nx = domain.nX;
    long ny = domain.nY;

    // Define weights for calculation
    array<double, 9> w1 = {1, 4, 1, 4, -20, 4, 1, 4, 1};

    switch(node->id) {
        case 0: {
            double sum = 0.0;
            for (int i = 0; i < 9; ++i) {
                long ixp = ix + e[i][0];
                long iyp = iy + e[i][1];

                // Apply boundary checks
                if (ixp >= nx) {
                    ixp = ix - e[i][0];
                } else if (ixp < 0) {
                    ixp = ix - e[i][0];
                }

                if (iyp >= ny) {
                    iyp = iy - e[i][1];
                } else if (iyp < 0) {
                    iyp = iy - e[i][1];
                }

                sum += w1[i] * domain.nodes[ixp][iyp]->*value / 6.0;
            }
            node->*solution = sum;
            break;
        }

        case 21:
            node->*solution = (2 * node->*value - 5 * domain.nodes[ix + 1][iy]->*value + 4 * domain.nodes[ix + 2][iy]->*value - domain.nodes[ix + 3][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy + 1]->*value + 4 * domain.nodes[ix][iy + 2]->*value - domain.nodes[ix][iy + 3]->*value);
            break;

        case 22: {
            long iyp = iy + 1 >= ny ? iy - 1 : iy + 1;
            long iyn = iy - 1 < 0 ? iy + 1 : iy - 1;

            node->*solution = (2 * node->*value - 5 * domain.nodes[ix + 1][iy]->*value + 4 * domain.nodes[ix + 2][iy]->*value - domain.nodes[ix + 3][iy]->*value) +
                                (domain.nodes[ix][iyp]->*value - 2 * node->*value + domain.nodes[ix][iyn]->*value);
            break;
        }

        case 23:
            node->*solution = (2 * node->*value - 5 * domain.nodes[ix + 1][iy]->*value + 4 * domain.nodes[ix + 2][iy]->*value - domain.nodes[ix + 3][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy - 1]->*value + 4 * domain.nodes[ix][iy - 2]->*value - domain.nodes[ix][iy - 3]->*value);
            break;

        case 24: {
            long ixp = ix + 1 >= nx ? ix - 1 : ix + 1;
            long ixn = ix - 1 < 0 ? ix + 1 : ix - 1;

            node->*solution = (domain.nodes[ixp][iy]->*value - 2 * node->*value + domain.nodes[ixn][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy + 1]->*value + 4 * domain.nodes[ix][iy + 2]->*value - domain.nodes[ix][iy + 3]->*value);
            break;
        }

        case 26: {
            long ixp = ix + 1 >= nx ? ix - 1 : ix + 1;
            long ixn = ix - 1 < 0 ? ix + 1 : ix - 1;

            node->*solution = (domain.nodes[ixp][iy]->*value - 2 * node->*value + domain.nodes[ixn][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy - 1]->*value + 4 * domain.nodes[ix][iy - 2]->*value - domain.nodes[ix][iy - 3]->*value);
            break;
        }

        case 27:
            node->*solution = (2 * node->*value - 5 * domain.nodes[ix - 1][iy]->*value + 4 * domain.nodes[ix - 2][iy]->*value - domain.nodes[ix - 3][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy + 1]->*value + 4 * domain.nodes[ix][iy + 2]->*value - domain.nodes[ix][iy + 3]->*value);
            break;

        case 28: {
            long iyp = iy + 1 >= ny ? iy - 1 : iy + 1;
            long iyn = iy - 1 < 0 ? iy + 1 : iy - 1;

            node->*solution = (2 * node->*value - 5 * domain.nodes[ix - 1][iy]->*value + 4 * domain.nodes[ix - 2][iy]->*value - domain.nodes[ix - 3][iy]->*value) +
                                (domain.nodes[ix][iyp]->*value - 2 * node->*value + domain.nodes[ix][iyn]->*value);
            break;
        }

        case 29:
            node->*solution = (2 * node->*value - 5 * domain.nodes[ix - 1][iy]->*value + 4 * domain.nodes[ix - 2][iy]->*value - domain.nodes[ix - 3][iy]->*value) +
                                (2 * node->*value - 5 * domain.nodes[ix][iy - 1]->*value + 4 * domain.nodes[ix][iy - 2]->*value - domain.nodes[ix][iy - 3]->*value);
            break;
    }
}

template <typename T1, typename T2>
void eGrad(Node* node, Domain &domain, T1 Node::*value, T2 Node::*gradient) {
    if (node->id == 20) return; // Solid wall

    long ix = node->x;
    long iy = node->y;
    long nx = domain.nX;
    long ny = domain.nY;

    // Calculate neighboring indices
    long ixp = ix + 1;
    long ixn = ix - 1;
    long iyp = iy + 1;
    long iyn = iy - 1;
    long ixpp = ixp + 1;
    long ixnn = ixn - 1;
    long iypp = iyp + 1;
    long iynn = iyn - 1;

    for (int i = 0; i < 9; i++) {
        if (i == 4) {
            (*node.*gradient)[i] = 0.0;
        } else {
            ixp = ix + e[i][0];
            ixn = ix - e[i][0];
            iyp = iy + e[i][1];
            iyn = iy - e[i][1];
            ixpp = ix + 2*e[i][0];
            ixnn = ix - 2*e[i][0];
            iypp = iy + 2*e[i][1];
            iynn = iy - 2*e[i][1];

            if (ixp > nx - 1) {
                ixp = ix - e[i][0];
            } else if (ixp < 0) {
                ixp = ix - e[i][0];
            }
            if (ixn > nx - 1) {
                ixn = ix + e[i][0];
            } else if (ixn < 0) {
                ixn = ix + e[i][0];
            }
            if (iyp > ny - 1) {
                iyp = iy - e[i][1];
            } else if (iyp < 0) {
                iyp = iy - e[i][1];
            }
            if (iyn > ny - 1) {
                iyn = iy + e[i][1];
            } else if (iyn < 0) {
                iyn = iy + e[i][1];
            }
            if (ixpp > nx - 1) {
                ixpp = (ixpp == nx) ? ix : ix - 2 * e[i][0];
            } else if (ixpp < 0) {
                ixpp = (ixpp == -1) ? ix : ix - 2 * e[i][0];
            }
            if (ixnn > nx - 1) {
                ixnn = (ixnn == nx) ? ix : ix + 2 * e[i][0];
            } else if (ixnn < 0) {
                ixnn = (ixnn == -1) ? ix : ix + 2 * e[i][0];
            }
            if (iypp > ny - 1) {
                iypp = (iypp == ny) ? iy : iy - 2 * e[i][1];
            } else if (iypp < 0) {
                iypp = (iypp == -1) ? iy : iy - 2 * e[i][1];
            }
            if (iynn > ny - 1) {
                iynn = (iynn == ny) ? iy : iy + 2 * e[i][1];
            } else if (iynn < 0) {
                iynn = (iynn == -1) ? iy : iy + 2 * e[i][1];
            }

            if (node->id == 0) {
                (*node.*gradient)[i] = 0.5*(domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value);
            } else if (node->id > 20) {
                if (1) { // if domain_info[ix,iy][i+3] == 1:
                    (*node.*gradient)[i] = 0.0;
                } else {
                    (*node.*gradient)[i] = 0.5 * (domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value);
                }
            }
        }
    }
}

/* The biased derivative functions are not used in the most recent solver.
template <typename T>
void biasedDerivativeX(Domain &domain, T Node::*value, T Node::*derivative, T Node::*biased_derivative) {
    long ixp, ixn, iyp, iyn, ixpp, ixnn, iypp, iynn;
    Node* node;

    long nx = domain.nX;
    long ny = domain.nY;

    for (long ix = 0; ix < nx; ix++) {
        for (long iy = 0; iy < ny; iy++) {
            node = domain.nodes[ix][iy];

            if (node->id == 20) continue; // Solid wall

            // Calculate neighboring indices
            ixp = ix + 1;
            ixn = ix - 1;
            iyp = iy + 1;
            iyn = iy - 1;
            ixpp = ixp + 1;
            ixnn = ixn - 1;
            iypp = iyp + 1;
            iynn = iyn - 1;

            if (ixp > nx - 1) {
                ixp = ix - 1;
            } else if (ixp < 0) {
                ixp = ix + 1;
            }

            if (iyp > ny - 1) {
                iyp = iy - 1;
            } else if (iyp < 0) {
                iyp = iy + 1;
            }

            if (ixn > nx - 1) {
                ixn = ix - 1;
            } else if (ixn < 0) {
                ixn = ix + 1;
            }

            if (iyn > ny - 1) {
                iyn = iy - 1;
            } else if (iyn < 0) {
                iyn = iy + 1;
            }

            if (ixpp > nx - 1) {
                ixpp = (ixpp == nx) ? ix : ix - 2;
            } else if (ixpp < 0) {
                ixpp = (ixpp == -1) ? ix : ix + 2;
            }

            if (iypp > ny - 1) {
                iypp = (iypp == ny) ? iy : iy - 2;
            } else if (iypp < 0) {
                iypp = (iypp == -1) ? iy : iy + 2;
            }

            if (ixnn > nx - 1) {
                ixnn = (ixnn == nx) ? ix : ix - 2;
            } else if (ixnn < 0) {
                ixnn = (ixnn == -1) ? ix : ix + 2;
            }

            if (iynn > ny - 1) {
                iynn = (iynn == ny) ? iy : iy - 2;
            } else if (iynn < 0) {
                iynn = (iynn == -1) ? iy : iy + 2;
            }

            double solution = 0.0;

            switch (node->id) {
                // Fluid nodes: central difference
                case 0:
                    solution = (-(domain.nodes[ixpp][iy]->*value) + 4 * domain.nodes[ixp][iy]->*value - 4 * domain.nodes[ixn][iy]->*value + domain.nodes[ixnn][iy]->*value) / 6 +
                                      (-(domain.nodes[ixpp][iypp]->*value) + 4 * domain.nodes[ixp][iyp]->*value - 4 * domain.nodes[ixn][iyn]->*value + domain.nodes[ixnn][iynn]->*value
                                      - domain.nodes[ixpp][iynn]->*value + 4 * domain.nodes[ixp][iyn]->*value - 4 * domain.nodes[ixn][iyp]->*value + domain.nodes[ixnn][iypp]->*value) / 24;
                    break;

                // East boundary nodes: forward difference
                case 21: 
                case 22: 
                case 23: 
                    solution = 0.5 * (-(domain.nodes[ixpp][iy]->*value) + 4 * domain.nodes[ixp][iy]->*value - 3 * domain.nodes[ix][iy]->*value);
                    break;

                case 27: 
                case 28: 
                case 29:
                    solution = -0.5 * (-(domain.nodes[ixnn][iy]->*value) + 4 * domain.nodes[ixn][iy]->*value - 3 * domain.nodes[ix][iy]->*value);
                    break;

                default: // Central difference for other boundary nodes
                    solution = 0.5 * (domain.nodes[ixp][iy]->*value - domain.nodes[ixn][iy]->*value);
                    break;
            }

            node->*biased_derivative = 0.5 * solution + 0.5 * node->*derivative;
        }
    }
}

template <typename T>
void biasedDerivativeY(Domain &domain, T Node::*value, T Node::*derivative, T Node::*biased_derivative) {
    long ixp, ixn, iyp, iyn, ixpp, ixnn, iypp, iynn;
    Node* node;

    long nx = domain.nX;
    long ny = domain.nY;

    for (long ix = 0; ix < nx; ix++) {
        for (long iy = 0; iy < ny; iy++) {
            node = domain.nodes[ix][iy];

            if (node->id == 20) continue; // Solid wall

            // Calculate neighboring indices
            ixp = ix + 1;
            ixn = ix - 1;
            iyp = iy + 1;
            iyn = iy - 1;
            ixpp = ixp + 1;
            ixnn = ixn - 1;
            iypp = iyp + 1;
            iynn = iyn - 1;

            if (ixp > nx - 1) {
                ixp = ix - 1;
            } else if (ixp < 0) {
                ixp = ix + 1;
            }

            if (iyp > ny - 1) {
                iyp = iy - 1;
            } else if (iyp < 0) {
                iyp = iy + 1;
            }

            if (ixn > nx - 1) {
                ixn = ix - 1;
            } else if (ixn < 0) {
                ixn = ix + 1;
            }

            if (iyn > ny - 1) {
                iyn = iy - 1;
            } else if (iyn < 0) {
                iyn = iy + 1;
            }

            if (ixpp > nx - 1) {
                ixpp = (ixpp == nx) ? ix : ix - 2;
            } else if (ixpp < 0) {
                ixpp = (ixpp == -1) ? ix : ix + 2;
            }

            if (iypp > ny - 1) {
                iypp = (iypp == ny) ? iy : iy - 2;
            } else if (iypp < 0) {
                iypp = (iypp == -1) ? iy : iy + 2;
            }

            if (ixnn > nx - 1) {
                ixnn = (ixnn == nx) ? ix : ix - 2;
            } else if (ixnn < 0) {
                ixnn = (ixnn == -1) ? ix : ix + 2;
            }

            if (iynn > ny - 1) {
                iynn = (iynn == ny) ? iy : iy - 2;
            } else if (iynn < 0) {
                iynn = (iynn == -1) ? iy : iy + 2;
            }

            double solution = 0.0;

            switch (node->id) {
                // Fluid nodes: central difference
                case 0:
                    solution = (-(domain.nodes[ix][iypp]->*value) + 4 * domain.nodes[ix][iyp]->*value - 4 * domain.nodes[ix][iyn]->*value + domain.nodes[ix][iynn]->*value)/6 +
                                      (-(domain.nodes[ixpp][iypp]->*value) + 4 * domain.nodes[ixp][iyp]->*value - 4 * domain.nodes[ixn][iyn]->*value + domain.nodes[ixnn][iynn]->*value
                                      - domain.nodes[ixpp][iynn]->*value + 4 * domain.nodes[ixp][iyn]->*value - 4 * domain.nodes[ixn][iyp]->*value + domain.nodes[ixnn][iypp]->*value)/24;
                    break;

                // East boundary nodes: forward difference
                case 21: 
                case 24:
                case 27: 
                    solution = 0.5 * (-(domain.nodes[ixpp][iy]->*value) + 4 * domain.nodes[ixp][iy]->*value - 3 * domain.nodes[ix][iy]->*value);
                    break;

                case 23:
                case 26: 
                case 29:
                    solution = -0.5 * (-(domain.nodes[ixnn][iy]->*value) + 4 * domain.nodes[ixn][iy]->*value - 3 * domain.nodes[ix][iy]->*value);
                    break;

                default: // Central difference for other boundary nodes
                    solution = 0.5 * (domain.nodes[ixp][iy]->*value - domain.nodes[ixn][iy]->*value);
                    break;
            }

            node->*biased_derivative = 0.5 * solution + 0.5 * node->*derivative;
        }
    }
}

template <typename T>
void biasedEGrad(Domain &domain, T Node::*value, T Node::*gradient) {
  long ixp, ixn, iyp, iyn, ixpp, ixnn, iypp, iynn;
    Node* node;

    long nx = domain.nX;
    long ny = domain.nY;

    for (long ix = 0; ix < nx; ix++) {
        for (long iy = 0; iy < ny; iy++) {
            node = domain.nodes[ix][iy];

            if (node->id == 20) continue; // Solid wall

            for (int i = 0; i < 9; i++) {
                if (i == 4) {
                    (*node.*gradient)[i] = 0.0;
                } else {
                    ixp = ix + e[i][0];
                    ixn = ix - e[i][0];
                    iyp = iy + e[i][1];
                    iyn = iy - e[i][1];
                    ixpp = ix + 2*e[i][0];
                    ixnn = ix - 2*e[i][0];
                    iypp = iy + 2*e[i][1];
                    iynn = iy - 2*e[i][1];

                    if (ixp > nx - 1) {
                        ixp = ix - e[i][0];
                    } else if (ixp < 0) {
                        ixp = ix - e[i][0];
                    }
                    if (ixn > nx - 1) {
                        ixn = ix + e[i][0];
                    } else if (ixn < 0) {
                        ixn = ix + e[i][0];
                    }
                    if (iyp > ny - 1) {
                        iyp = iy - e[i][1];
                    } else if (iyp < 0) {
                        iyp = iy - e[i][1];
                    }
                    if (iyn > ny - 1) {
                        iyn = iy + e[i][1];
                    } else if (iyn < 0) {
                        iyn = iy + e[i][1];
                    }
                    if (ixpp > nx - 1) {
                        ixpp = (ixpp == nx) ? ix : ix - 2 * e[i][0];
                    } else if (ixpp < 0) {
                        ixpp = (ixpp == -1) ? ix : ix - 2 * e[i][0];
                    }
                    if (ixnn > nx - 1) {
                        ixnn = (ixnn == nx) ? ix : ix + 2 * e[i][0];
                    } else if (ixnn < 0) {
                        ixnn = (ixnn == -1) ? ix : ix + 2 * e[i][0];
                    }
                    if (iypp > ny - 1) {
                        iypp = (iypp == ny) ? iy : iy - 2 * e[i][1];
                    } else if (iypp < 0) {
                        iypp = (iypp == -1) ? iy : iy - 2 * e[i][1];
                    }
                    if (iynn > ny - 1) {
                        iynn = (iynn == ny) ? iy : iy + 2 * e[i][1];
                    } else if (iynn < 0) {
                        iynn = (iynn == -1) ? iy : iy + 2 * e[i][1];
                    }

                    // THIS SECTION BELOW IS NOT FINISHED
                    if (node->id == 0) {
                        if (domain.nodes[ixpp][iypp]->id == 20) {
                            (*node.*gradient)[i] = (0.5*(domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value))/2;
                        } else {
                            (*node.*gradient)[i] = (0.5*-(domain.nodes[ixpp][iypp]->*value) + 4*domain.nodes[ixp][iyp]->*value - 3*domain.nodes[ix][iy]->*value + 0.5*(domain.nodes[ixp][iyp]->*value - domain.nodes[ixn][iyn]->*value))/2;
                        }
                    } else if (node->id > 20) {
                        if (1) {
                            if(1) {

                            } else {

                            }
                        } else {

                        }
                    }
                }
            }
        }
    }
    return;
}
*/

#endif // MATH_OPERATIONS_HPP
