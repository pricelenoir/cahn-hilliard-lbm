#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

#include "constants.hpp"

Constants::Constants(const nlohmann::json& config) {
    // Define domain properties
    lx = config["domain"]["lx"];
    ly = config["domain"]["ly"];
    nLB = config["domain"]["nLB"];

    ref_len = ly;               
    deltaX = ly / (nLB - 1);
    deltaY = deltaX; // Uniform grid

    nx = lx / deltaX;
    ny = ly / deltaY;

    gX = config["gravity"]["gx"];
    gY = config["gravity"]["gy"];

    tau = config["simulation"]["tau"];
    axisSymmetry = config["simulation"]["axis_symmetry"];

    // Define fluid properties
    nuP = config["fluid"]["nuP"];
    rhoP = config["fluid"]["rhoP"];
    rho0 = config["fluid"]["rho0"];
    Re = config["fluid"]["Re"];
    
    phi1 = config["order"]["phi1"];
    phi2 = config["order"]["phi2"];

    sigmaP = config["fluid"]["sigmaP"];
    sigma = config["fluid"]["sigma"]; // Lattice units (must be less than 1e-3)

    deltaM = (rhoP / rho0) * pow(deltaX, 3);
    deltaT = sqrt((sigma / sigmaP) * deltaM);

    uP = (nuP * Re) / ref_len;
    nulb = nuP / (pow(deltaX, 2) / deltaT);

    uLB = (uP != 0) ? uP / (deltaX / deltaT) : deltaT / deltaX; // Velocity in lattice units (should be much smaller than the speed of sound, 1/3)

    densityRatio = config["ratios"]["density"];
    viscosityRatio = config["ratios"]["viscosity"];

    // Species 1 (water), Species 2 (air)
    rho1 = rho0;
    rho2 = densityRatio * rho0;
    nu1 = nulb;
    nu2 = viscosityRatio * nulb;

    mu1 = rho1 * nu1; // Dynamic viscosity
    mu2 = rho2 * nu2;

    // Cahn-Hilliard parameters define surface tension
    W = 5;

    if (axisSymmetry) {
        k = (3.0 / 2.0) * sigma * W / pow(phi1 - phi2, 2);
        beta = (12.0 / pow(phi1 - phi2, 4)) * (sigma / W);
    } else {
        k = (3.0 / 2.0) * (0.5 * sigma) * W / pow(phi1 - phi2, 2);
        beta = (12.0 / pow(phi1 - phi2, 4)) * ((0.5 * sigma) / W);
    }

    D_2 = 2.42e-7;
    D_1 = 2.42e-9;

    We = pow(uLB, 2) * (ref_len / deltaX) * (rho1 - rho2) / sigma;
    Pe = (uP != 0) ? (1.1 / beta) : (1.1 / beta);

    eta = 2;
    mbl = (eta * tau) / 3;
    contactAngle = M_PI - (M_PI / 3.0 + M_PI / 2.0);
}

#endif // CONSTANTS_CPP
