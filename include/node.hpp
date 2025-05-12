#ifndef NODE_HPP
#define NODE_HPP

#include <array>

class Node {
public:
    int x;
    int y;

    /* 
    Node ID Key:
    ┌────---┬───────────────────────────┬───────────┐
    │ ID    │ Description               │ Normal ID │
    ├────---┼───────────────────────────┼───────────┤
    │ 0     │ General fluid node        │           │
    │ 20    │ Solid node                │           │
    │ 21-29 │ Fluid node next to solid  │           │
    ├────---┼───────────────────────────┼───────────┤
    │ 21    │ NE                        │ -5        │
    │ 22    │ E                         │ -1        │
    │ 23    │ SE                        │ -7        │
    │ 24    │ N                         │ -3        │
    │ 25    │ Middle                    │           │
    │ 26    │ S                         │ -4        │
    │ 27    │ NW                        │ -8        │
    │ 28    │ W                         │ -2        │
    │ 29    │ SW                        │ -6        │
    └────---┴───────────────────────────┴───────────┘
    */

    // Node classification
    int id;
    int normalNodeID;
    int inletNormalNodeID;
    int outletNormalNodeID;
    bool isBoundary; // Boundary node (wall or internal solid)
    bool isConcave;  // Concave solid node
    bool isInlet;    // Inlet boundary condition
    bool isOutlet;   // Outlet boundary condition

    // Physical properties
    double uX;       // Velocity in x direction
    double uY;       // Velocity in y direction
    double phi;      // Order parameter
    double rho;      // Density
    double tau;      // Relaxation time
    double p;        // Pressure
    double pStar;    // Corrected pressure value from paper
    double pThermo;  // Thermodynamic pressure
    double nu;       // Kinematic viscosity
    double mu;       // Dynamic viscosity
    double forceX;   // Summation of forces term: gravity, surface, pressure, and viscous forces
    double forceY;   // Used in NS eqn
    double tmp;      // Temporary variable to store intermediate values

    // Residual values
    double phi0;
    double p0;
    double uX0;
    double uY0;
    
    double oldUX;
    double oldUY;
    double oldPhi;
    double oldP;
    double oldMu;
    double mu0;

    // Gradient values
    double dudx;
    double dudy;
    double dvdx;
    double dvdy;
    double dpdx;
    double dpdy;
    double dpStardx;
    double dpStardy;
    double dphidx;
    double dphidy;
    double d2phidx2;
    double uSqr;
    std::array<double, 9> eDudy;
    std::array<double, 9> eDvdx;

    // Particle distribution functions
    std::array<double, 9> gIn;
    std::array<double, 9> gOut;
    std::array<double, 9> gEq;
    std::array<double, 9> sourceG;

    std::array<double, 9> hIn;
    std::array<double, 9> hOut;
    std::array<double, 9> hEq;
    std::array<double, 9> sourceH;
    
    std::array<bool, 9> neighborLookUp; // Neighbor look up table for boundary nodes

    Node();
};

#endif // NODE_HPP