#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
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
using namespace std::chrono;

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

    int nX = constants.nx;
    int nY = constants.ny;
    double deltaX = constants.deltaX;
    double deltaT = constants.deltaT;
    double deltaM = constants.deltaM;

    // Create domain and initialize node variables (u, p, phi, rho, mu0, tau)
    Domain domain(nX, nY);
    domain.initialize(config, constants);
    boundaryClassification(domain);

    Node* node;
    // Initial calculation of mu
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            node = domain.nodes[i][j];
            laplace(node, domain, &Node::phi, &Node::oldMu);
            node->oldMu = node->mu0 - (constants.k * node->oldMu);
        }
    }

    // Initial gradient and equilibrium calculations
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            node = domain.nodes[i][j];

            derivativeY(node, domain, &Node::uX, &Node::dudy);
            derivativeX(node, domain, &Node::uY, &Node::dvdx);
            derivativeX(node, domain, &Node::uX, &Node::dudx);
            derivativeY(node, domain, &Node::uY, &Node::dvdy);
            eGrad(node, domain, &Node::uX, &Node::eDudy);
            eGrad(node, domain, &Node::uY, &Node::eDvdx);

            node->uSqr = pow(node->uX, 2) + pow(node->uY, 2);

            equilibriumG(node, domain, &Node::gIn);
            equilibriumH(node, domain, &Node::hIn);
        }
    }
    
    string domainDir = config["domain"]["domain_dir"];
    string outputDir = config["simulation"]["output_dir"];
    int maxIter = config["simulation"]["max_iter"];
    int saveDomainIter = config["simulation"]["save_domain_iter"];

    cout << "\033[1;34m";
    cout << "========================================\n";
    cout << " Cahn-Hilliard Lattice Boltzmann Solver\n";
    cout << "========================================\033[0m\n\n";
    cout << "\033[1;36m" << left;
    cout << setw(22) << "Domain:"            << "\033[0m" << domainDir << '\n';
    cout << "\033[1;36m" << setw(22) << "Output:"            << "\033[0m" << outputDir << '\n';
    cout << "\033[1;36m" << setw(22) << "Iterations:"        << "\033[0m" << maxIter << '\n';
    cout << "\033[1;36m" << setw(22) << "Save interval:"     << "\033[0m" << saveDomainIter << " iterations\n";
    cout << "\033[1;36m" << setw(22) << "Grid size:"         << "\033[0m" << nX << " x " << nY << '\n';
    cout << "\033[1;36m" << setw(22) << "Periodicity:"       << "\033[0m"
        << (config["domain"]["periodicity"]["x"] ? "[1" : "[0") << ","
        << (config["domain"]["periodicity"]["y"] ? "1]" : "0]") << '\n';
    cout << "\033[1;36m" << setw(22) << "Density ratio:"     << "\033[0m" << constants.densityRatio << '\n';
    cout << "\033[1;36m" << setw(22) << "Viscosity ratio:"   << "\033[0m" << constants.viscosityRatio << '\n';
    cout << "\033[1;36m" << setw(22) << "Gravity:"           << "\033[0m[" << constants.gX << "," << constants.gY << "]" <<'\n';
    cout << "\n\033[1;32m✔ Setup complete:\033[0m Starting simulation...\n";

    double resU, resPhi, resP;
    auto start = high_resolution_clock::now();
    try {
        for (int iter = 1; iter <= maxIter; iter++) {
            // Reset residuals
            resU = 0;
            resPhi = 0;
            resP = 0;

            // Update macroscopic values
            macroscopic(domain, constants);

            // Update boundary conditions
            macroscopicInflowBC(domain, constants);
            macroscopicOutflowBC(domain, constants);
            macroscopicWallBC(domain, constants);

            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nY; j++) {
                    node = domain.nodes[i][j];
                    node->oldMu = domain.nodes[i][j]->mu;
                }
            }

            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nY; j++) {
                    node = domain.nodes[i][j];

                    // Gradient calculations
                    derivativeY(node, domain, &Node::uX, &Node::dudy);
                    derivativeX(node, domain, &Node::uY, &Node::dvdx);
                    derivativeX(node, domain, &Node::uX, &Node::dudx);
                    derivativeY(node, domain, &Node::uY, &Node::dvdy);
                    eGrad(node, domain, &Node::uX, &Node::eDudy);
                    eGrad(node, domain, &Node::uY, &Node::eDvdx);

                    node->uSqr = pow(node->uX, 2) + pow(node->uY, 2);

                    // Update equilibrium values
                    equilibriumG(node, domain, &Node::gEq);
                    equilibriumH(node, domain, &Node::hEq);
                    
                    node->pStar = node->p / (node->rho / 3.0);
                }
            }
            
            // Update source terms
            sourceG(domain, constants);

            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nY; j++) {
                    node = domain.nodes[i][j];

                    sourceH(node, domain, constants);

                    // Zou He boundary condition
                    zouHeBC(node, domain);

                    // Collide and stream steps
                    collide(node, domain);
                    stream(node, domain);

                    // Calculate residuals and update old values
                    if (node->id == 0) {
                        resU   = max(sqrt(pow(abs(node->uX0 - node->uX), 2) + pow(abs(node->uY0 - node->uY), 2)) / (deltaT / deltaX), resU);
                        resPhi = max(abs(node->phi0 - node->phi), resPhi);
                        resP   = max(abs(node->p0 - node->p) / ((deltaX * pow(deltaT, 2)) / deltaM), resP);
    
                        node->uX0  = node->uX;
                        node->uY0  = node->uY;
                        node->phi0 = node->phi;
                        node->p0   = node->p;
                    }
                }
            }

            // Store residual at each timestep
            if (saveDomainIter != 0) {
                domain.resUVec.push_back(resU);
                domain.resPhiVec.push_back(resPhi);
                domain.resPVec.push_back(resP);

                if (iter % saveDomainIter == 0 && iter != 1) {
                    cout << "\033[1;36m[Iteration " << right << setw(7) << iter << " / " << maxIter << "]\033[0m Saving domain...\n";
                    domain.save(config, iter);
                }
            }

            // Check for convergence
            if (resU <= 1e-6 && resPhi <= 1e-6) {
                cout << "\033[1;32m✔ Convergence reached.\033[0m Stopping simulation at iteration " << iter << ".\n";
                domain.save(config, iter);
                break;
            } else if (iter > 1000 && isinf(resU)) {
                cout << "\033[1;31m✘ Simulation failed:\033[0m Solution diverged at iteration " << iter << ".\n";
                return 1;
            }
        }
    } catch (runtime_error &e) {
        cout << "\033[1;31m✘ Simulation stopped:\033[0m Error occurred during execution.\n";
        cout << e.what() << endl;
        return 1;
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<chrono::duration<double>>(end - start);
    auto elapsed = duration_cast<seconds>(end - start).count();

    int hours   = elapsed / 3600;
    int minutes = (elapsed % 3600) / 60;
    int seconds = elapsed % 60;
    
    cout << "\033[1;33m⚑ Simulation complete.\033[0m\n";
    cout << "\033[1;35m⏱ Total time: \033[0m" << right
              << setfill('0') << setw(2) << hours << ":"
              << setfill('0') << setw(2) << minutes << ":"
              << setfill('0') << setw(2) << seconds << '\n';
    return 0;
}