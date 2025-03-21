#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include "domain.hpp"
#include <vector>

void equilibriumG(Node* node, Domain &domain, std::vector<double> Node::*g);
void equilibriumH(Node* node, Domain &domain, std::vector<double> Node::*h);
void sourceG(Node* node, Domain &domain, Constants &constants);
void sourceH(Node* node, Domain &domain, Constants &constants);

#endif // EQUILIBRIUM_HPP