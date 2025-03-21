#include "node.hpp"

Node::Node(int x_coord, int y_coord) : x(x_coord), y(y_coord) {
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
	mbl = 0;
	nu = 0;
	mu = 0;
	forceX = 0;
	forceY = 0;
	tmp = 0;

	oldUX = 0;
	oldUY = 0;
	oldPhi = 1;
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

	eDudy.resize(9, 0);
	eDvdx.resize(9, 0);

	gIn.resize(9, 0);
	gOut.resize(9, 0);
	gEq.resize(9, 0);
	sourceG.resize(9, 0);

	hIn.resize(9, 0);
	hOut.resize(9, 0);
	hEq.resize(9, 0);
	sourceH.resize(9, 0);
}
