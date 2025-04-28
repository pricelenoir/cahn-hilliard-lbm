#include "collideAndStream.hpp"

using namespace std;

void collide(Node* node, Domain &domain) {
    if (node->id == 20) return;

    for (int i = 0; i < 9; i++) {
        node->gOut[i] = node->gIn[i] - (node->gIn[i] - node->gEq[i])/(node->tau + 0.5) + node->sourceG[i];
        node->hOut[i] = node->hIn[i] - (node->hIn[i] - node->hEq[i])/(node->tau + 0.5) + node->sourceH[i];
    }
}

void stream(Node* node, Domain &domain) {
    if (node->id == 20) return; // Solid wall

    int nx = domain.nX;
    int ny = domain.nY;
    int x = node->x;
    int y = node->y;
    int nextX, nextY;

    for (int i = 0; i < 9; i++) {
        nextX = x + e[i][0];
        nextY = y + e[i][1];

        // Handle x-boundary conditions
        if (domain.periodicity[0] == 1) {
            // Periodic boundary
            if (nextX < 0) {
                nextX = nx - 1;
            } else if (nextX >= nx) {
                nextX = 0;
            }
        } else {
            // Stream to infinity boundary
            if (nextX < 0 || nextX >= nx) return;
        }

        // Handle y-boundary conditions
        if (domain.periodicity[1] == 1) {
            // Periodic boundary
            if (nextY < 0) {
                nextY = ny - 1;
            } else if (nextY >= ny) {
                nextY = 0;
            }
        } else {
            // Stream to infinity boundary
            if (nextY < 0 || nextY >= ny) return;
        }

        // Skip if the target node is an obstacle
        if (domain.nodes[nextX][nextY]->id == 20) return;

        domain.nodes[nextX][nextY]->gIn[i] = node->gOut[i];
        domain.nodes[nextX][nextY]->hIn[i] = node->hOut[i];
    }
}
