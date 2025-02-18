#ifndef BOUNDARY_CLASSIFIER_HPP
#define BOUNDARY_CLASSIFIER_HPP

#include "domain.hpp"

std::vector<Node*> identifyNeighbors(Domain &domain, Node* node);
void identifyBoundaryNodes(Domain &domain);
double euclideanDistance(Node* node, std::pair<int, int> neighbor);
std::vector<std::pair<int, int>> kClosestNeighbors(Node* node, std::vector<std::pair<int, int>> edgeNeighbors, int k);
void classifyBC(Domain &domain);
void identifyInletNodes(Domain &domain);
void identifyOutletNodes(Domain &domain);

#endif // BOUNDARY_CLASSIFIER_HPP