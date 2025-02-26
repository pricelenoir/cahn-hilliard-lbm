#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "domain.hpp"

void macroscopicInflowBC(Domain &domain, Constants &constants);
void macroscopicOutflowBC(Domain &domain, Constants &constants);
void macroscopicWallBC(Domain &domain, Constants &constants);
void zouHeBC(Node* node, Domain &domain);

#endif // BOUNDARY_CONDITIONS_HPP