#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include "domain.hpp"
#include <vector>

void equilibriumG(Node* node, Domain &domain, std::array<double, 9> Node::*g);
void equilibriumH(Node* node, Domain &domain, std::array<double, 9> Node::*h);
void sourceG(Domain &domain, Constants &constants);
void sourceH(Node* node, Domain &domain, Constants &constants);

#endif // EQUILIBRIUM_HPP