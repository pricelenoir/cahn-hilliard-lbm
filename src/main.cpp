#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "json.hpp"
#include "boundaryClassifier.hpp"
#include "boundaryConditions.hpp"
#include "collideAndStream.hpp"
#include "constants.hpp"
#include "domain.hpp"
#include "equilibrium.hpp"
#include "macroscopic.hpp"
#include "mathOperations.hpp"

using json = nlohmann::json;
using namespace std;

int main(int argc, char** argv) {
    // Create a JSON object and read from config file
    json config;
    ifstream inFile("config.json");

    if (!inFile.is_open()) {
        cerr << "Failed to open config.json for reading." << endl;
        return 1;
    }
    inFile >> config;
    inFile.close();

    // Initialize simulation constants from config
    Constants constants(config);

    long nX = constants.nx;
    long nY = constants.ny;
    double deltaX = constants.deltaX;
    double deltaT = constants.deltaT;
    double deltaM = constants.deltaM;

    // Create domain and initialize node variables (u, p, phi, rho, mu0, tau)
    Domain domain(nX, nY);
    domain.initialize(config, constants);
    boundaryClassification(domain);

    Node* node;
    // Initial calculation of mu
    for (long i = 0; i < nX; i++) {
        for (long j = 0; j < nY; j++) {
            node = domain.nodes[i][j];
            laplace(node, domain, &Node::phi, &Node::oldMu);
            node->oldMu = node->mu0 - (constants.k * node->oldMu);
        }
    }

    for (long i = 0; i < nX; i++) {
        for (long j = 0; j < nY; j++) {
            node = domain.nodes[i][j];

            // Gradient calculations
            laplace(node, domain, &Node::oldMu, &Node::d2mu);
            derivativeY(node, domain, &Node::uX, &Node::dudy);
            derivativeX(node, domain, &Node::uY, &Node::dvdx);
            derivativeX(node, domain, &Node::uX, &Node::dudx);
            derivativeY(node, domain, &Node::uY, &Node::dvdy);
            eGrad(node, domain, &Node::uX, &Node::eDudy);
            eGrad(node, domain, &Node::uY, &Node::eDvdx);

            node->uSqr = pow(node->uX, 2) + pow(node->uY, 2);

            // Initial equilibrium values
            equilibriumG(node, domain);
            equilibriumH(node, domain);
        }
    }
    // Add some print stuff here about the domain and sim
    cout << "Setup complete. Starting simulation..." << endl;

    int maxIter = config["simulation"]["max_iter"];
    int saveDomainIter = config["simulation"]["save_domain_iter"];
    double resU, resPhi, resP;
    try {
        for (int iter = 0; iter < maxIter; iter++) {
            // Update macroscopic values
            macroscopic(domain, constants);

            // Update boundary conditions
            macroscopicInflowBC(domain, constants);
            macroscopicOutflowBC(domain, constants);
            macroscopicWallBC(domain, constants);

            for (long i = 0; i < nX; i++) {
                for (long j = 0; j < nY; j++) {
                    node = domain.nodes[i][j];
                    node->oldMu = domain.nodes[i][j]->mu;
                }
            }

            for (long i = 0; i < nX; i++) {
                for (long j = 0; j < nY; j++) {
                    node = domain.nodes[i][j];

                    // Gradient calculations
                    laplace(node, domain, &Node::oldMu, &Node::d2mu);
                    derivativeY(node, domain, &Node::uX, &Node::dudy);
                    derivativeX(node, domain, &Node::uY, &Node::dvdx);
                    derivativeX(node, domain, &Node::uX, &Node::dudx);
                    derivativeY(node, domain, &Node::uY, &Node::dvdy);
                    eGrad(node, domain, &Node::uX, &Node::eDudy);
                    eGrad(node, domain, &Node::uY, &Node::eDvdx);

                    node->uSqr = pow(node->uX, 2) + pow(node->uY, 2);

                    // Update equilibrium values and sources
                    equilibriumG(node, domain);
                    sourceG(node, domain, constants);
                    equilibriumH(node, domain);
                    sourceH(node, domain, constants);

                    // Zou He boundary condition
                    zouHeBC(node, domain);

                    // Collide and stream steps
                    collide(node, domain);
                    stream(node, domain);

                    // Calculate residuals
                    if (node->id == 0) {
                        resU   = sqrt(pow(abs(node->oldUX - node->uX), 2) + pow(abs(node->oldUY - node->uY), 2)) / (deltaT / deltaX);
                        resPhi = abs(node->oldPhi - node->phi);
                        resP   = abs(node->oldP - node->p) / ((deltaX * pow(deltaT, 2)) / deltaM);
    
                        node->oldUX  = node->uX;
                        node->oldUY  = node->uY;
                        node->oldPhi = node->phi;
                        node->oldP   = node->p;
    
                        if (domain.resU < resU) domain.resU = resU;
                        if (domain.resPhi < resPhi) domain.resPhi = resPhi;
                        if (domain.resP < resP) domain.resP = resP;
                    }
                }
            }

            if (iter % saveDomainIter == 0 && iter != 0) {
                cout << "Iteration: " << iter << " / " << maxIter << ". Saving domain..." << endl;
                domain.save(config, iter);
            }

            // Check for convergence
            if (domain.resU <= 1e-6 && domain.resPhi <= 1e-6) {
                cout << "Convergence reached. Stopping simulation at iteration " << iter << "." << endl;
                domain.save(config, iter);
                return 0;
            } else if (iter > 1000 && isinf(domain.resU)) {
                cout << "Simulation failed. Solution diverged." << endl;
                return 1;
            }
        }
    } catch (runtime_error &e) {
        cout << "Simulation stopped. Error." << endl;
        return 1;
    }
    cout << "Simulation complete." << endl;
    return 0;
}