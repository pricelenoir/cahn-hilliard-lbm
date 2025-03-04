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

Domain::Domain(long nx, long ny) : nX(nx), nY(ny) {
    nodes.resize(nX);

    for (long i = 0; i < nX; ++i) {
        nodes[i].resize(nY, nullptr);
        for (long j = 0; j < nY; ++j) {
            nodes[i][j] = new Node(i, j);
        }
    }
    resU = 1;
    resP = 1;
    resPhi = 1;
}

// Lambda function to read in CSV input files
template <typename T>
void readInputFile(long nX, long nY, const string& filename, function<void(long, long, T)> populateFunc) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    long i = 0;
    while (getline(fin, line)) {
        stringstream ss(line);
        T value;
        long j = 0;

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
void writeOutputFile(long nX, long nY, const string& filename, function<T(long, long)> retrieveFunc) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    fout << fixed; // Ensure fixed-point notation for floating-point numbers

    for (long i = 0; i < nX; i++) {
        for (long j = 0; j < nY; j++) {
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

    if (!fs::exists("domain/")) {
        cout << "Domain directory not found." << endl;
    }

    string domainDir = config["domain"]["domain_dir"];
    domainDir = "domain/" + domainDir;

    // Read in domain file
    string inputFile = domainDir + "/domain.txt";
    readInputFile<int>(nX, nY, inputFile, [this](long i, long j, int value) {
        nodes[i][j]->id = value;
    });

    // Read in pressure file
    inputFile = domainDir + "/p.txt";
    readInputFile<double>(nX, nY, inputFile, [this](long i, long j, double value) {
        nodes[i][j]->p = value;
    });

    // Read in order parameter file
    inputFile = domainDir + "/phi.txt";
    readInputFile<double>(nX, nY, inputFile, [this](long i, long j, double value) {
        nodes[i][j]->phi = value;
    });

    // Read in velocity (x-direction) file
    inputFile = domainDir + "/uX.txt";
    readInputFile<double>(nX, nY, inputFile, [this](long i, long j, double value) {
        nodes[i][j]->uX = value;
    });

    // Read in velocity (y-direction) file
    inputFile = domainDir + "/uY.txt";
    readInputFile<double>(nX, nY, inputFile, [this](long i, long j, double value) {
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
    if (!fs::exists("output/")) {
        fs::create_directory("output/");
    }

    string outputDir = config["simulation"]["output_dir"];
    outputDir = "output/" + outputDir;
    if (!fs::exists(outputDir)) {
        fs::create_directory(outputDir);
    }

    // Create iteration output directory
    outputDir = outputDir + "/iter" + to_string(iter) + "/";
    if (!fs::exists(outputDir)) {
        fs::create_directory(outputDir);
    }

    string outputFile = outputDir + "domain.txt";
    writeOutputFile<int>(nX, nY, outputFile, [this](long i, long j) {
        return nodes[i][j]->id;
    });

    // Write out pressure file
    outputFile = outputDir + "p.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](long i, long j) {
        return nodes[i][j]->p;
    });

    // Write out order parameter file
    outputFile = outputDir + "phi.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](long i, long j) {
        return nodes[i][j]->phi;
    });

    // Write out velocity (x-direction) file
    outputFile = outputDir + "uX.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](long i, long j) {
        return nodes[i][j]->uX;
    });

    // Write out velocity (y-direction) file
    outputFile = outputDir + "uY.txt";
    writeOutputFile<double>(nX, nY, outputFile, [this](long i, long j) {
        return nodes[i][j]->uY;
    });
}
