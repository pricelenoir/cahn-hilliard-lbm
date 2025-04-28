#include "boundaryConditions.hpp"
#include "domain.hpp"
#include "mathOperations.hpp"

using namespace std;
#include <iostream>

void macroscopicInflowBC(Domain &domain, Constants &constants) {
    // Definition of constants
    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double deltaM = constants.deltaM;
    double deltaX = constants.deltaX;
    double deltaT = constants.deltaT;

    int xp, xpp, yp, ypp;
    Node* node;
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node->oldUX = node->uX;
            node->oldUY = node->uY;
        }
    }

    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            if (node->isInlet) {
                if (node->inletNormalNodeID == -1) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j;
                    ypp = j;
                } else if (node->inletNormalNodeID == -2) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j;
                    ypp = j;
                } else if (node->inletNormalNodeID == -3) {
                    xp = i;
                    xpp = i;
                    yp = j + 1;
                    ypp = j + 2;
                } else if (node->inletNormalNodeID == -4) {
                    xp = i;
                    xpp = i;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->inletNormalNodeID == -5) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j + 1;
                    ypp = j + 2;
                } else if (node->inletNormalNodeID == -6) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->inletNormalNodeID == -7) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->inletNormalNodeID == -8) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j + 1;
                    ypp = j + 2;
                }

                // Neumann boundary condition
                node->uX = (-2.0 * 0.0 - domain.nodes[xpp][ypp]->oldUX + 4.0 * domain.nodes[xp][yp]->oldUX) / 3.0;
                node->uY = (-2.0 * 0.0 - domain.nodes[xpp][ypp]->oldUY + 4.0 * domain.nodes[xp][yp]->oldUY) / 3.0;

                node->phi = phi1;
                if (node->phi <= min(phi1, phi2)) {
                    node->phi = min(phi1, phi2);
                } else if (node->phi >= max(phi1, phi2)) {
                    node->phi = max(phi1, phi2);
                }

                node->mu = 0.0;
                node->p = 1256 / (deltaM / (deltaX * deltaT * deltaT)); // Boundary pressure of 1200Pa + 10% than the equilibrium value
            }

            node->rho = (rho1 * (node->phi - phi2) + rho2 * (phi1 - node->phi)) / (phi1 - phi2);
        }
    }
}

void macroscopicOutflowBC(Domain &domain, Constants &constants) {
    // Definition of constants
    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double deltaM = constants.deltaM;
    double deltaX = constants.deltaX;
    double deltaT = constants.deltaT;


    int xp, xpp, yp, ypp;
    Node* node;
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node->oldUX = node->uX;
            node->oldUY = node->uY;
        }
    }

    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            if (node->isOutlet) {
                if (node->outletNormalNodeID == -1) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j;
                    ypp = j;
                } else if (node->outletNormalNodeID == -2) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j;
                    ypp = j;
                } else if (node->outletNormalNodeID == -3) {
                    xp = i;
                    xpp = i;
                    yp = j + 1;
                    ypp = j + 2;
                } else if (node->outletNormalNodeID == -4) {
                    xp = i;
                    xpp = i;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->outletNormalNodeID == -5) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j + 1;
                    ypp = j + 2;
                } else if (node->outletNormalNodeID == -6) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->outletNormalNodeID == -7) {
                    xp = i + 1;
                    xpp = i + 2;
                    yp = j - 1;
                    ypp = j - 2;
                } else if (node->outletNormalNodeID == -8) {
                    xp = i - 1;
                    xpp = i - 2;
                    yp = j + 1;
                    ypp = j + 2;
                }

                // Neumann boundary condition
                node->uX = (-2.0 * 0.0 - domain.nodes[xpp][ypp]->oldUX + 4.0 * domain.nodes[xp][yp]->oldUX) / 3.0;
                node->uY = (-2.0 * 0.0 - domain.nodes[xpp][ypp]->oldUY + 4.0 * domain.nodes[xp][yp]->oldUY) / 3.0;

                node->phi = phi2;
                if (node->phi <= min(phi1, phi2)) {
                    node->phi = min(phi1, phi2);
                } else if (node->phi >= max(phi1, phi2)) {
                    node->phi = max(phi1, phi2);
                }

                node->mu = 0.0;
                node->p = 0.0;
            }

            node->rho = (rho1 * (node->phi - phi2) + rho2 * (phi1 - node->phi)) / (phi1 - phi2);
        }
    }
}

void macroscopicWallBC(Domain &domain, Constants &constants) {
    // Definition of constants
    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double nu1 = constants.nu1;
    double nu2 = constants.nu2;
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;
    double beta = constants.beta;
    double k = constants.k;
    double contactAngle = constants.contactAngle;
    double gX = constants.gX;
    double gY = constants.gY;
    double phiContactAngleThreshold = 0.5 * constants.W * constants.deltaX / constants.ly;

    double cubicPhiFunction, gradPhiNew, gradPStarNew;
    int xp, xpp, xppp, yp, ypp, yppp;
    int dx, dy;

    Node* node;
    for (int i = 0; i < domain.nX; i++) {
        for (int j = 0; j < domain.nY; j++) {
            node = domain.nodes[i][j];

            node->oldPhi = node->phi;
            node->oldMu = node->mu;
            node->oldUX = node->uX;
            node->oldUY = node->uY;
            node->oldP = node->p;
            node->pStar = 3 * node->p / node->rho;

            derivativeY(node, domain, &Node::uX, &Node::dudy);
            derivativeX(node, domain, &Node::uY, &Node::dvdx);
            derivativeY(node, domain, &Node::uX, &Node::dudx);
            derivativeX(node, domain, &Node::uY, &Node::dvdy);

            node->nu = (nu1 * (node->rho - rho2) + nu2 * (rho1 - node->rho)) / (rho1 - rho2);

            // Gradient calculation of phi
            derivativeX(node, domain, &Node::phi, &Node::dphidx);
            derivativeY(node, domain, &Node::phi, &Node::dphidy);

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

            if (node->isBoundary) {
                if (node->normalNodeID == -1) {
                    xp = i + 1;
                    xpp = i + 2;
                    xppp = i + 3;
                    yp = j;
                    ypp = j;
                    yppp = j;
                    dx = 1;
                    dy = 0;
                } else if (node->normalNodeID == -2) {
                    xp = i - 1;
                    xpp = i - 2;
                    xppp = i - 3;
                    yp = j;
                    ypp = j;
                    yppp = j;
                    dx = -1;
                    dy = 0;
                } else if (node->normalNodeID == -3) {
                    xp = i;
                    xpp = i;
                    xppp = i;
                    yp = j + 1;
                    ypp = j + 2;
                    yppp = j + 3;
                    dx = 0;
                    dy = 1;
                } else if (node->normalNodeID == -4) {
                    xp = i;
                    xpp = i;
                    xppp = i;
                    yp = j - 1;
                    ypp = j - 2;
                    yppp = j - 3;
                    dx = 0;
                    dy = -1;
                } else if (node->normalNodeID == -5) {
                    xp = i + 1;
                    xpp = i + 2;
                    xppp = i + 3;
                    yp = j + 1;
                    ypp = j + 2;
                    yppp = j + 3;
                    dx = 1;
                    dy = 1;
                } else if (node->normalNodeID == -6) {
                    xp = i - 1;
                    xpp = i - 2;
                    xppp = i - 3;
                    yp = j - 1;
                    ypp = j - 2;
                    yppp = j - 3;
                    dx = -1;
                    dy = -1;
                } else if (node->normalNodeID == -7) {
                    xp = i + 1;
                    xpp = i + 2;
                    xppp = i + 3;
                    yp = j - 1;
                    ypp = j - 2;
                    yppp = j - 3;
                    dx = 1;
                    dy = -1;
                } else if (node->normalNodeID == -8) {
                    xp = i - 1;
                    xpp = i - 2;
                    xppp = i - 3;
                    yp = j + 1;
                    ypp = j + 2;
                    yppp = j + 3;
                    dx = -1;
                    dy = 1;
                }

                cubicPhiFunction = -(phiContactAngleThreshold - domain.nodes[xp][yp]->oldPhi) * ((1.0 - phiContactAngleThreshold) - domain.nodes[xp][yp]->oldPhi) / ((1.0 - 2.0 * phiContactAngleThreshold) * (1.0 - 2.0 * phiContactAngleThreshold));
                if (cubicPhiFunction < 0.0) cubicPhiFunction = 0.0;

                // Neumann boundary condition
                gradPhiNew = (-1.0) * sqrt(2.0 * beta / k) * cos(contactAngle) * cubicPhiFunction;
                domain.nodes[xp][yp]->phi = (-2.0 * gradPhiNew - domain.nodes[xppp][yppp]->oldPhi + 4.0 * domain.nodes[xpp][ypp]->oldPhi) / 3.0;

                cubicPhiFunction = -(phiContactAngleThreshold - node->oldPhi) * ((1.0 - phiContactAngleThreshold) - node->oldPhi) / ((1.0 - 2.0 * phiContactAngleThreshold) * (1.0 - 2.0 * phiContactAngleThreshold));
                if (cubicPhiFunction <= 0) cubicPhiFunction = 0.0;

                // Neumann boundary condition
                gradPhiNew = (-1.0) * sqrt(2.0 * beta / k) * cos(contactAngle) * cubicPhiFunction;
                node->phi = (-2.0 * gradPhiNew - domain.nodes[xpp][ypp]->oldPhi + 4.0 * domain.nodes[xp][yp]->oldPhi) / 3.0;

                if (node->phi <= min(phi1, phi2)) {
                    node->phi = min(phi1, phi2);
                } else if (node->phi >= max(phi1, phi2)) {
                    node->phi = max(phi1, phi2);
                }

                if (domain.nodes[xp][yp]->phi <= min(phi1, phi2)) {
                    domain.nodes[xp][yp]->phi = min(phi1, phi2);
                } else if (domain.nodes[xp][yp]->phi >= max(phi1, phi2)) {
                    domain.nodes[xp][yp]->phi = max(phi1, phi2);
                }

                node->rho = (rho1 * (node->phi - phi2) + rho2 * (phi1 - node->phi)) / (phi1 - phi2);
                domain.nodes[xp][yp]->rho = (rho1 * (domain.nodes[xp][yp]->phi - phi2) + rho2 * (phi1 - domain.nodes[xp][yp]->phi)) / (phi1 - phi2);

                node->mu = 0.0;
                domain.nodes[xp][yp]->mu = 0.0;

                node->uX = 0.0;
                node->uY = 0.0;
                domain.nodes[xp][yp]->uX = 0.0;
                domain.nodes[xp][yp]->uY = 0.0;

                // Neumann boundary condition
                gradPStarNew = sqrt(pow(domain.nodes[xp][yp]->forceX * dx, 2) + pow(domain.nodes[xp][yp]->forceY * dy, 2));
                domain.nodes[xp][yp]->p = ((-2 * gradPStarNew - domain.nodes[xppp][yppp]->pStar + 4.0 * domain.nodes[xpp][ypp]->pStar) / 3.0) * (domain.nodes[xp][yp]->rho / 3.0);

                // Neumann boundary condition
                gradPStarNew = sqrt(pow(node->forceX * dx, 2) + pow(node->forceY * dy, 2));
                node->p = ((-2 * gradPStarNew - domain.nodes[xpp][ypp]->pStar + 4.0 * domain.nodes[xp][yp]->pStar) / 3.0) * (node->rho / 3.0);
            }
        }
    }
}

void zouHeBC(Node* node, Domain &domain) {
    if (node->isBoundary || node->isInlet || node->isOutlet) {
        array<int, 3> bounceBackIndices;
        array<int, 2> concaveIndices;
        int normalNodeID;

        if (node->isBoundary) normalNodeID = node->normalNodeID;
        if (node->isInlet)    normalNodeID = node->inletNormalNodeID;
        if (node->isOutlet)   normalNodeID = node->outletNormalNodeID;

        if (normalNodeID == -1) {
            bounceBackIndices = {6, 7, 8};
        } else if (normalNodeID == -2) {
            bounceBackIndices = {0, 1, 2};
        } else if (normalNodeID == -3) {
            bounceBackIndices = {2, 5, 8};
        } else if (normalNodeID == -4) {
            bounceBackIndices = {0, 3, 6};
        } else if (normalNodeID == -5) {
            bounceBackIndices = {5, 7, 8};
        } else if (normalNodeID == -6) {
            bounceBackIndices = {0, 1, 3};
        } else if (normalNodeID == -7) {
            bounceBackIndices = {3, 6, 7};
        } else if (normalNodeID == -8) {
            bounceBackIndices = {1, 2, 5};
        }

        // Apply Zou-He bounce-back condition
        for (int idx : bounceBackIndices) {
            node->gIn[idx] = node->gEq[idx] + node->gIn[8 - idx] - node->gEq[8 - idx];
            node->hIn[idx] = node->hEq[idx] + node->hIn[8 - idx] - node->hEq[8 - idx];
        }

        if (node->isConcave) {
            if (node->normalNodeID == -5) {
                concaveIndices = {6, 2};
            } else if (node->normalNodeID == -6) {
                concaveIndices = {6, 2};
            } else if (node->normalNodeID == -7) {
                concaveIndices = {0, 8};
            } else if (node->normalNodeID == -8) {
                concaveIndices = {0, 8};
            }

            // Apply concave bounce-back condition
            for (int idx : concaveIndices) {
                node->gIn[idx] = node->gEq[idx] + node->gIn[8 - idx] - node->gEq[8 - idx];
                node->hIn[idx] = node->hEq[idx] + node->hIn[8 - idx] - node->hEq[8 - idx];
            }
        }
    } else return;
}
