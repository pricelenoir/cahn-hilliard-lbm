#ifndef NODE_HPP
#define NODE_HPP

#include <vector>

class Node {
public:
    long x;
    long y;

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
    double mbl;      // Mobility
    double nu;       // Kinematic viscosity
    double mu;       // Dynamic viscosity
    double forceX;   // Summation of forces term: gravity, surface, pressure, and viscous forces
    double forceY;   // Used in NS eqn
    double tmp;      // Temporary variable to store intermediate values

    // Residual values
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
    std::vector<double> eDvdx;
    std::vector<double> eDudy;

    // Particle distribution functions
    std::vector<double> gIn;
    std::vector<double> gOut;
    std::vector<double> gEq;
    std::vector<double> sourceG;

    std::vector<double> hIn;
    std::vector<double> hOut;
    std::vector<double> hEq;
    std::vector<double> sourceH;

    std::vector<bool> neighborLookUp; // Neighbor look up table for boundary nodes

    Node(int posX, int posY);
};

#endif // NODE_HPP