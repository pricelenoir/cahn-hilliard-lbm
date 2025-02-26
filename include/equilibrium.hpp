#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include "domain.hpp"

void equilibriumG(Node* node, Domain &domain);
void equilibriumH(Node* node, Domain &domain);
void sourceG(Node* node, Domain &domain, Constants &constants);
void sourceH(Node* node, Domain &domain, Constants &constants);

#endif // EQUILIBRIUM_HPP