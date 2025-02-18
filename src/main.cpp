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

    // Create domain and initialize node variables (u, p, phi, rho, mu0, tau)
    Domain domain(nX, nY);
    domain.initialize(config, constants);

    // Identify and classify nodes with boundary conditions
    identifyBoundaryNodes(domain);
    classifyBC(domain);
    identifyInletNodes(domain);
    identifyOutletNodes(domain);

    // Initial calculation of mu
    laplace(domain, &Node::phi, &Node::oldMu);
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            Node* node = domain.nodes[i][j];
            node->oldMu = node->mu0 - (constants.k * node->oldMu);
        }
    }

    // Gradient calculations
    laplace(domain, &Node::oldMu, &Node::d2mu);
    derivativeY(domain, &Node::uX, &Node::dudy);
    derivativeX(domain, &Node::uY, &Node::dvdx);
    derivativeX(domain, &Node::uX, &Node::dudx);
    derivativeY(domain, &Node::uY, &Node::dvdy);
    eGrad(domain, &Node::uX, &Node::eDudy);
    eGrad(domain, &Node::uY, &Node::eDvdx);
    uSqr(domain);

    // Initial equilibrium distribution based on the macroscopic values
    equilibriumG(domain);
    equilibriumH(domain);

    cout << "Setup complete. Starting simulation..." << endl;

    int maxIter = config["simulation"]["max_iter"];
    int saveDomainIter = config["simulation"]["save_domain_iter"];
    try {
        for (int iter = 0; iter < maxIter; iter++) {
            // Update macroscopic values
            macroscopic(domain, constants);

            // Update boundary conditions
            macroscopicInflowBC(domain);
            macroscopicOutflowBC(domain);
            macroscopicWallBC(domain);

            // Update old mu values
            for (long i = 0; i < nX; i++) {
                for (long j = 0; j < nY; j++) {
                    domain.nodes[i][j]->oldMu = domain.nodes[i][j]->mu;
                }
            }

            // Gradient calculations
            laplace(domain, &Node::mu, &Node::d2mu);
            derivativeY(domain, &Node::uX, &Node::dudy);
            derivativeX(domain, &Node::uY, &Node::dvdx);
            derivativeX(domain, &Node::uX, &Node::dudx);
            derivativeY(domain, &Node::uY, &Node::dvdy);
            eGrad(domain, &Node::uX, &Node::eDudy);
            eGrad(domain, &Node::uY, &Node::eDvdx);
            uSqr(domain);

            // Update equilibrium values and sources
            equilibriumG(domain);
            sourceG(domain, constants);

            equilibriumH(domain);
            sourceH(domain, constants);
            
            // Zou He boundary condition
            zouHeBC(domain);

            // Collide and stream steps
            collide(domain);
            stream(domain);

            domain.calcResiduals();

            if (iter % saveDomainIter == 0 && iter != 0) {
                cout << "Iteration: " << iter << " / " << maxIter << ". Saving domain..." << endl;
                domain.save(config, iter);
            }
        }
    } catch (runtime_error &e) {
        cout << "Simulation stopped. Error." << endl;
        return 1;
    }
    cout << "Simulation complete." << endl;
    return 0;
}