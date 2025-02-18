#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "domain.hpp"

void neumannBC(Domain &domain);
void macroscopicInflowBC(Domain &domain);
void macroscopicOutflowBC(Domain &domain);
void macroscopicWallBC(Domain &domain);
void zouHeBC(Domain &domain);

#endif // BOUNDARY_CONDITIONS_HPP