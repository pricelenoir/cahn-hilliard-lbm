#include "equilibrium.hpp"
#include "mathOperations.hpp"
#include <cmath>

void equilibriumG(Domain &domain) {
    double fEq, gammaU;
    long duE = 0;
    long eU = 0;

    Node* node;
    for (long i; i < domain.nX; i++) {
        for (long j; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node->pStar = node->p / (node->rho / 3);

            if (node->id != 20) {
                for (int k; k < 9; k++) {
                    /* Gamma equilibrium calculation */
                    eU = e[k][0] * node->uX + e[k][1] * node->uY;
                    duE = node->eDvdx[k] * e[k][1] + node->eDudy[k] * e[k][0];
                    gammaU = w[k] * (1 + 3*eU + 4.5*pow(eU, 2) - 1.5*node->uSqr  + 3*(node->tau * duE));

                    /* Equilibrium f calculation */
                    fEq = node->rho * gammaU;
                    
                    /* Equilibrium g calculation */
                    node->gEq[k] = w[k] * node->pStar + (fEq - w[k]);
                }
            }
        }
    }
}

void equilibriumH(Domain &domain) {
    double fEq, gammaU;
    long duE = 0;
    long eU = 0;

    Node* node;
    for (long i = 0; i < domain.nX; i++) {
        for (long j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

             if (node->id != 20) {
                for (int k; k < 9; k++) {
                    /* Gamma equilibrium calculation */
                    eU = e[k][0] * node->uX + e[k][1] * node->uY;
                    duE = node->eDvdx[k] * e[k][1] + node->eDudy[k] * e[k][0];
                    gammaU = w[k] * (1 + 3*eU + 4.5*pow(eU, 2) - 1.5*node->uSqr  + 3*(node->tau * duE));

                    /* Equilibrium f calculation */
                    fEq = node->rho * gammaU;

                    /* Equilibrium h calculation */
                    node->hEq[k] = fEq * node->phi;
                }
            }
        }
    }
}

void sourceG(Domain &domain, Constants &constants) {
    // Definition of constants
    double nu1 = constants.nu1;
    double nu2 = constants.nu2;
    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double gX = constants.gX;
    double gY = constants.gY;

    Node* node;
    for (long i = 0; i < domain.nX; i++) {
        for (long j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node-> nu = (nu1 * (node->rho - rho2) + nu2 * (rho1 - node->rho)) / (rho1 - rho2);

            // Gradient calculations of phi
            derivativeX(node, domain, &Node::phi, &Node::dphidx);
            derivativeY(node, domain, &Node::phi, &Node::dphidy);

            // Gradient calculations of pStar
            derivativeX(node, domain, &Node::pStar, &Node::dpStardx);
            derivativeY(node, domain, &Node::pStar, &Node::dpStardy);

            // Calculate intermediate value for force calculation below
            node->tmp = node->rho * node->nu * (node->dudx + node->dvdy + node->dvdx + node->dudy);
            derivativeX(node, domain, &Node::tmp, &Node::forceX);
            derivativeY(node, domain, &Node::tmp, &Node::forceY);

            // Calculate force (forceX and forceY currently hold intermediate gradient values)
            node->forceX = gX * node->rho + node->forceX - (1/3) * node->pStar * (rho1 - rho2) * node->dphidx + node->mu * node->dphidx;
            node->forceY = gY * node->rho + node->forceY - (1/3) * node->pStar * (rho1 - rho2) * node->dphidy + node->mu * node->dphidy;

            long eU = 0;
            if (node->id != 20) {
                for (int k; k < 9; k++) {
                    eU = e[k][0] * node->uX + e[k][1] * node->uY;
                    node->sourceG[k] = 0.5 * w[k] * (3 * (e[k][0] - node->uX) * node->forceX + 3 * (e[k][1] - node->uY) * node->forceY + 9 * eU * (e[k][0] * node->forceX + e[k][1] * node->forceY));
                }
            }
        }
    }
}

void sourceH(Domain &domain, Constants &constants) {
    // Definition of constants
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double beta = constants.beta;
    double k = constants.k;

    long eU = 0;
    double L2NormEGradPhi, interfaceNormalX, interfaceNormalY;
    double L2NormThrshld = 0.00005;
    Node* node;
    for (long i = 0; i < domain.nX; i++) {
        for (long j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            // Gradient calculations of phi
            derivativeX(node, domain, &Node::phi, &Node::dphidx);
            derivativeY(node, domain, &Node::phi, &Node::dphidy);

            double L2NormEGradPhi = sqrt(pow(node->dphidx, 2) + pow(node->dphidy, 2));
            
            if (node->id != 20) {
                if (L2NormEGradPhi >= L2NormThrshld) {
                    interfaceNormalX = node->dphidx / L2NormEGradPhi;
                    interfaceNormalY = node->dphidy / L2NormEGradPhi;
                }

                node->forceX = -4 * (node->phi - phi1) * (node->phi - phi2) * (0.25 * sqrt(2 * beta / k)) * interfaceNormalX / (phi1 - phi2);
                node->forceY = -4 * (node->phi - phi1) * (node->phi - phi2) * (0.25 * sqrt(2 * beta / k)) * interfaceNormalY / (phi1 - phi2);
                                
                for (int k; k < 9; k++) {
                    eU = e[k][0] * node->uX + e[k][1] * node->uY;
                    node->sourceH[k] = 0.5 * w[k] * node->mbl * (3 * (e[k][0] - node->uX) * node->forceX + 3 * (e[k][1] - node->uY) * node->forceY + 9 * eU * (e[k][0] * node->forceX + e[k][1] * node->forceY));
                }
            }
        }
    }
}