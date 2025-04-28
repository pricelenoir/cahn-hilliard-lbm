#include "equilibrium.hpp"
#include "mathOperations.hpp"
#include <cmath>

using namespace std;

void equilibriumG(Node* node, Domain &domain, array<double, 9> Node::*g) {
    double fEq, gammaU;
    int duE = 0;
    int eU = 0;

    node->pStar = node->p / (node->rho / 3);

    if (node->id != 20) {
        for (int k = 0; k < 9; k++) {
            /* Gamma equilibrium calculation */
            eU = e[k][0] * node->uX + e[k][1] * node->uY;
            duE = node->eDvdx[k] * e[k][1] + node->eDudy[k] * e[k][0];
            gammaU = w[k] * (1.0 + 3.0*eU + 4.5*pow(eU, 2) - 1.5*node->uSqr  + 3.0*(node->tau * duE));

            /* Equilibrium f calculation */
            fEq = node->rho * gammaU;
            fEq = 1.0 * gammaU;
            
            /* Equilibrium g calculation */
            (node->*g)[k] = w[k] * node->pStar + (fEq - w[k]);
        }
    }
}

void equilibriumH(Node* node, Domain &domain, array<double, 9> Node::*h) {
    if (node->id == 20) return; // Solid wall

    double fEq, gammaU;
    int duE = 0;
    int eU = 0;

    for (int k = 0; k < 9; k++) {
        /* Gamma equilibrium calculation */
        eU = e[k][0] * node->uX + e[k][1] * node->uY;
        duE = node->eDvdx[k] * e[k][1] + node->eDudy[k] * e[k][0];
        gammaU = w[k] * (1.0 + 3.0*eU + 4.5*pow(eU, 2) - 1.5*node->uSqr  + 3.0*(node->tau * duE));

        /* Equilibrium f calculation */
        fEq = node->rho * gammaU;
        fEq = 1.0 * gammaU;

        /* Equilibrium h calculation */
        (node->*h)[k] = fEq * node->phi;
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
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node->nu = (nu1 * (node->rho - rho2) + nu2 * (rho1 - node->rho)) / (rho1 - rho2);

            // Gradient calculations of phi
            derivativeX(node, domain, &Node::phi, &Node::dphidx);
            derivativeY(node, domain, &Node::phi, &Node::dphidy);
        
            // Gradient calculations of pStar
            derivativeX(node, domain, &Node::pStar, &Node::dpStardx);
            derivativeY(node, domain, &Node::pStar, &Node::dpStardy);
        
            // Calculate intermediate value for force calculation below
            node->tmp = node->rho * node->nu * (node->dudx + node->dvdy + node->dvdx + node->dudy);
        }
    }

    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            derivativeX(node, domain, &Node::tmp, &Node::forceX);
            derivativeY(node, domain, &Node::tmp, &Node::forceY);
        
            // Calculate force (forceX and forceY currently hold intermediate gradient values)
            node->forceX = gX * node->rho + node->forceX - (1.0/3.0) * node->pStar * (rho1 - rho2) * node->dphidx + node->mu * node->dphidx;
            node->forceY = gY * node->rho + node->forceY - (1.0/3.0) * node->pStar * (rho1 - rho2) * node->dphidy + node->mu * node->dphidy;

            int eU = 0;
            if (node->id != 20) {
                for (int k = 0; k < 9; k++) {
                    eU = e[k][0] * node->uX + e[k][1] * node->uY;
                    node->sourceG[k] = 0.5 * w[k] * (3.0 * (e[k][0] - node->uX) * node->forceX + 3.0 * (e[k][1] - node->uY) * node->forceY + 9.0 * eU * (e[k][0] * node->forceX + e[k][1] * node->forceY));
                }
            }
        }
    }
}

void sourceH(Node* node, Domain &domain, Constants &constants) {
    // Definition of constants
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double beta = constants.beta;
    double k = constants.k;
    double mbl = constants.mbl;

    int eU = 0;
    double L2NormEGradPhi, interfaceNormalX, interfaceNormalY;
    double L2NormThrshld = 0.00005;

    // Gradient calculations of phi
    derivativeX(node, domain, &Node::phi, &Node::dphidx);
    derivativeY(node, domain, &Node::phi, &Node::dphidy);

    L2NormEGradPhi = sqrt(pow(node->dphidx, 2) + pow(node->dphidy, 2));
    
    if (node->id != 20) {
        if (L2NormEGradPhi >= L2NormThrshld) {
            interfaceNormalX = node->dphidx / L2NormEGradPhi;
            interfaceNormalY = node->dphidy / L2NormEGradPhi;
        }

        node->forceX = -4.0 * (node->phi - phi1) * (node->phi - phi2) * (0.25 * sqrt(2.0 * beta / k)) * interfaceNormalX / (phi1 - phi2);
        node->forceY = -4.0 * (node->phi - phi1) * (node->phi - phi2) * (0.25 * sqrt(2.0 * beta / k)) * interfaceNormalY / (phi1 - phi2);
                        
        for (int k = 0; k < 9; k++) {
            eU = e[k][0] * node->uX + e[k][1] * node->uY;
            node->sourceH[k] = 0.5 * w[k] * mbl * (3.0 * (e[k][0] - node->uX) * node->forceX + 3.0 * (e[k][1] - node->uY) * node->forceY + 9.0 * eU * (e[k][0] * node->forceX + e[k][1] * node->forceY));
        }
    }
}