#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <functional>
#include "json.hpp"
#include "domain.hpp"
#include "node.hpp"

using namespace std;
namespace fs = std::filesystem;

Domain::Domain(int nx, int ny) : nX(nx), nY(ny) {
    nodesChunk.resize(nX * nY); // Allocate all Nodes in one block
    nodes.resize(nX);

    for (int i = 0; i < nX; ++i) {
        nodes[i].resize(nY, nullptr);
        for (int j = 0; j < nY; ++j) {
            Node* nodePtr = &nodesChunk[i * nY + j]; // Get address inside contiguous block
            nodePtr->x = i;
            nodePtr->y = j;
            nodes[i][j] = nodePtr; // Point the 2D vector to the right Node
        }
    }
}

// Lambda function to read in CSV input files
template <typename T>
void readInputFile(int nX, int nY, const string& filename, function<void(int, int, T)> populateFunc) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    int i = 0;
    while (getline(fin, line)) {
        stringstream ss(line);
        T value;
        int j = 0;

        while (ss >> value) {
            if (i < nX && j < nY) {
                populateFunc(i, j, value);
                j++;
                if (ss.peek() == ',') ss.ignore();
            }
        }
        i++;
        if (i >= nX) break;
    }
    fin.close();
}

// Lambda function to write CSV output files
template <typename T>
void writeOutputFile(int nX, int nY, const string& filename, function<T(int, int)> retrieveFunc) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    fout << fixed; // Ensure fixed-point notation for floating-point numbers

    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            if constexpr (is_same_v<T, double>) {
                fout << setprecision(8) << retrieveFunc(i, j); // Set precision for doubles
            } else {
                fout << retrieveFunc(i, j); // Default behavior for other types
            }

            if (j < nY - 1) fout << ",";
        }
        fout << endl;
    }
    fout.close();
}

void Domain::initialize(const nlohmann::json& config, Constants &constants) {
    // Read periodicity from config
    periodicity[0] = config["domain"]["periodicity"]["x"] ? 1 : 0;
    periodicity[1] = config["domain"]["periodicity"]["y"] ? 1 : 0;

    string domainDir = config["domain"]["domain_dir"];

    // Read in domain file
    string inputFile = domainDir + "/domain.txt";
    readInputFile<int>(nX, nY, inputFile, [this](int i, int j, int value) {
        nodes[i][j]->id = value;
    });

    // Read in pressure file
    inputFile = domainDir + "/p.txt";
    readInputFile<double>(nX, nY, inputFile, [this](int i, int j, double value) {
        nodes[i][j]->p = value;
    });

    // Read in order parameter file
    inputFile = domainDir + "/phi.txt";
    readInputFile<double>(nX, nY, inputFile, [this](int i, int j, double value) {
        nodes[i][j]->phi = value;
    });

    // Read in velocity (x-direction) file
    inputFile = domainDir + "/uX.txt";
    readInputFile<double>(nX, nY, inputFile, [this](int i, int j, double value) {
        nodes[i][j]->uX = value;
    });

    // Read in velocity (y-direction) file
    inputFile = domainDir + "/uY.txt";
    readInputFile<double>(nX, nY, inputFile, [this](int i, int j, double value) {
        nodes[i][j]->uY = value;
    });

    double rho1 = constants.rho1;
    double rho2 = constants.rho2;
    double beta = constants.beta;
    double phi1 = constants.phi1;
    double phi2 = constants.phi2;

    // Set properties for each node
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            Node* node = nodes[i][j];

            // Calculate density from composition and mu0
            node->rho = (rho1 * (node->phi - phi2) + rho2 * (phi1 - node->phi)) / (phi1 - phi2);
            node->mu0 = 2 * beta * (node->phi - phi2) * (node->phi - phi1) * (2 * node->phi - (phi2 + phi1));
            node->tau = config["simulation"]["tau"];
        }
    }
}

void Domain::save(const nlohmann::json& config, int iter) {
    // Create output directory
    string outputDir = config["simulation"]["output_dir"];
    if (!fs::exists(outputDir)) {
        fs::create_directories(outputDir);
    }

    string residualsFile = outputDir + "/residuals.txt";

    // Create iteration output directory
    string iterDir = outputDir + "/iter" + to_string(iter) + "/";
    if (!fs::exists(iterDir)) {
        fs::create_directory(iterDir);
    }

    string outputFile = iterDir + "domain.txt";
    writeOutputFile<int>(nX, nY, outputFile, [this](int i, int j) {
        return nodes[i][j]->id;
    });

    // Write out pressure file
    outputFile = iterDir + "p.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](int i, int j) {
        return nodes[i][j]->p;
    });

    // Write out order parameter file
    outputFile = iterDir + "phi.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](int i, int j) {
        return nodes[i][j]->phi;
    });

    // Write out velocity (x-direction) file
    outputFile = iterDir + "uX.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](int i, int j) {
        return nodes[i][j]->uX;
    });

    // Write out velocity (y-direction) file
    outputFile = iterDir + "uY.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](int i, int j) {
        return nodes[i][j]->uY;
    });

    // Check if this is the first save of the run
    bool isFirstSave = (iter == config["simulation"]["save_domain_iter"]);

    ofstream resFile;
    if (isFirstSave) {
        resFile.open(residualsFile);
    } else {
        resFile.open(residualsFile, ios::app);
    }

    if (resFile.is_open()) {
        size_t firstIteration = iter - resUVec.size() + 1;

        for (size_t i = 0; i < resUVec.size(); i++) {
            size_t iterationNumber = firstIteration + i;
            resFile << "Iteration " << iterationNumber << ":\n";
            resFile << "------------------\n";
            resFile << fixed << setprecision(10);
            resFile << "Residual U: " << resUVec[i] << "\n";
            resFile << "Residual Phi: " << resPhiVec[i] << "\n";
            resFile << "Residual P: " << resPVec[i] << "\n";
            resFile << "\n";
        }
        // Clear the vectors after writing to the file
        resUVec.clear();
        resPhiVec.clear();
        resPVec.clear();
    } else {
        cerr << "Error opening file: " << residualsFile << endl;
    }
}