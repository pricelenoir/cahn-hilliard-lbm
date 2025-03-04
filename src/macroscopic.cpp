#include "macroscopic.hpp"
#include "mathOperations.hpp"
#include <cmath>

using namespace std;

void macroscopic(Domain &domain, Constants &constants) {
    // Definition of constants
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double beta = constants.beta;
    double nu1 = constants.nu1;
    double nu2 = constants.nu2;
    double gX = constants.gX;
    double gY = constants.gY;
    double k = constants.k;

    Node* node;
    for (long i = 0; i < domain.nX; i++) {
        for (long j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            // Calculate gradient of phi
            derivativeX(node, domain, &Node::phi, &Node::dphidx);
            derivativeY(node, domain, &Node::phi, &Node::dphidy);

            // Gradient calculations
            derivativeX(node, domain, &Node::uX, &Node::dudx);
            derivativeY(node, domain, &Node::uX, &Node::dudy);
            derivativeX(node, domain, &Node::uY, &Node::dvdx);
            derivativeY(node, domain, &Node::uY, &Node::dvdy);

            node->pStar = 0;
            node->phi = 0;
            node->uX = 0;
            node->uY = 0;

            for (long k = 0; k < 9; k++) {
                // Calculates composition from CH equilibrium functions
                node->phi += node->hIn[k];

                // Calculates velocity from NS equilibrium functions
                node->uX += e[k][0] * node->gIn[k];
                node->uY += e[k][1] * node->gIn[k];

                // Calculates pressure from NS equilibrium functions
                node->pStar += node->gIn[k];
            }

            // Set limits for composition to remain between 0 and 1
            if (node->phi <= min(phi1, phi2)) {
                node->phi = min(phi1, phi2);
            } else if (node->phi >= max(phi1, phi2)) {
                node->phi = max(phi1, phi2);
            }

            // Calculate density from composition
            node->rho = (rho1 * (node->phi - phi2) + rho2 * (phi1 - node->phi)) / (phi1 - phi2);
        }
    }

    // Need to iterate through entire domain before using laplace() (depends on neighboring values)
    for (long i = 0; i < domain.nX; i++) {
        for (long j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            // Calculate mu
            node->mu0 = 4 * beta * (node->phi - phi2) * (node->phi - phi1) * (node->phi - 0.5 * (phi2 + phi1));
            laplace(node, domain, &Node::phi, &Node::d2phidx2);
            node->mu = node->mu0 - k * node->d2phidx2;
            node->mu = 2.2 * node->mu;

            node->nu = (nu1 * (node->rho - rho2) + nu2 * (rho1 - node->rho)) / (rho1 - rho2);
            node->p = node->pStar * node->rho / 3;

            // Calculate intermediate values for force calculation below
            node->tmp = node->rho * node->nu * (node->dudx + node->dvdy + node->dvdx + node->dudy);
            derivativeX(node, domain, &Node::tmp, &Node::forceX);
            derivativeY(node, domain, &Node::tmp, &Node::forceY);

            /// Calculate force (forceX and forceY currently hold intermediate gradient values)
            node->forceX = gX * node->rho + node->forceX - (1/3) * node->pStar * (rho1 - rho2) * node->dphidx + node->mu * node->dphidx;
            node->forceY = gY * node->rho + node->forceY - (1/3) * node->pStar * (rho1 - rho2) * node->dphidy + node->mu * node->dphidy;
            
            // Adjusting pressure and velocities based on Lee paper suggestions (macroscopic equations)
            node->uX = node->uX + (0.5 * node->forceX) / node->rho;
            node->uY = node->uY + (0.5 * node->forceY) / node->rho;

            // Thermodynamic pressure calculation
            node->pThermo = node->p + (node->phi * node->mu0 - beta * pow((pow(node->phi, 2) - node->phi), 2)) - k * node->phi * node->d2phidx2 + 0.5 * k * (pow(node->dphidx, 2) + pow(node->dphidy, 2));
        }
    }
}