#include "node.hpp"

Node::Node() {
	x = 0;
	y = 0;
	id = 0;

	isBoundary = false;
	isConcave = false;
	isInlet = false;
	isOutlet = false;
	normalNodeID = 0;
	inletNormalNodeID = 0;
	outletNormalNodeID = 0;
	
	uX = 0;
	uY = 0;
	
	phi = 0;
	rho = 0;
	tau = 0;
	p = 0;
	pStar = 0;
	pThermo = 0;
	nu = 0;
	mu = 0;
	forceX = 0;
	forceY = 0;
	tmp = 0;

	phi0 = 1;
	p0 = 1;
	uX0 = 1;
	uY0 = 1;

	oldUX = 0;
	oldUY = 0;
	oldPhi = 0;
	oldP = 0;
	oldMu = 0;
	mu0 = 0;

	dudx = 0;
	dudy = 0;
	dvdx = 0;
	dvdy = 0;
	dpdx = 0;
	dpdy = 0;
	dpStardx = 0;
	dpStardy = 0;
	dphidx = 0;
	dphidy = 0;
	d2phidx2 = 0;
	uSqr = 0;

    eDudy.fill(0.0);
    eDvdx.fill(0.0);

    gIn.fill(0.0);
    gOut.fill(0.0);
    gEq.fill(0.0);
    sourceG.fill(0.0);
	
    hIn.fill(0.0);
    hOut.fill(0.0);
    hEq.fill(0.0);
    sourceH.fill(0.0);
}
