#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "json.hpp"

struct Constants {
    // Lattice definition
    int e[9][2];
    double weight[9];
    double h[9];

    // Define domain properties
    double lx;             // Physical domain size [meters] in x and y directions
    double ly;
    long nLB;              // Number of lattice units (characteristic length scale)

    double ref_len;        // Reference length in physical units
    double deltaX;         // Lattice spacing in the x and y directions
    double deltaY;         

    long nx;               // Number of lattice points in the x and y directions
    long ny;

    // Define fluid properties
    double nuP;            // Viscosity of water at 80°C [m^2/s]
    double rhoP;           // Density of water at 80°C [kg/m^3]
    double rho0;           // Density of water in lattice units
    double Re;             // Reynolds number
    
    double phi1;          // Limits of the order parameter
    double phi2;
    
    double sigmaP;        // Surface tension of the air-water interface in physical units [kg/s^2]
    double sigma;         // Lattice units (must be less than 1e-3)
    double deltaM;        // Mass per lattice volume
    double deltaT;        // Lattice time step
    double uP;            // Velocity in physical units
    double nulb;          // Kinematic viscosity in lattice units
    double uLB;           // Velocity in lattice units (should be much smaller than the speed of sound, 1/3)
    double densityRatio;
    double viscosityRatio;

    // Species 1 (water), Species 2 (air)
    double rho1, rho2;
    double nu1, nu2; 
    double mu1, mu2;     // Dynamic viscosity

    // Cahn-Hilliard parameters define surface tension
    double diff;         // Diffusive interface thickness
    double k;            // Surface tension coefficient
    double beta;         // Energy parameter

    double D_2;          // Diffusion coefficient of O2 in water [m^2/s]
    double D_1;          // Diffusion coefficient of water in water [m^2/s]

    double We;           // Weber number
    double Pe;           // Pecelet number definition
    double mlb;          // Mobility
    double contactAngle;

    Constants(const nlohmann::json& config); // Constructor
};

#endif // CONSTANTS_HPP
