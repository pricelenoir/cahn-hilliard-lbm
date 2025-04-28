#ifndef MATH_OPERATIONS_HPP
#define MATH_OPERATIONS_HPP

#include "domain.hpp"

void derivativeX(Node* node, Domain &domain, double Node::*value, double Node::*derivative);
void derivativeY(Node* node, Domain &domain, double Node::*value, double Node::*derivative);
void laplace(Node* node, Domain &domain, double Node::*value, double Node::*solution);
void eGrad(Node* node, Domain &domain, double Node::*value, std::array<double, 9> Node::*gradient);

#endif // MATH_OPERATIONS_HPP
