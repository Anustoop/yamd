#include "verlet.h"

// First step of the Verlet integration: update positions
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    // Update positions based on current velocities and forces
    x += vx * timestep + 0.5 * fx * timestep * timestep;
    y += vy * timestep + 0.5 * fy * timestep * timestep;
    z += vz * timestep + 0.5 * fz * timestep * timestep;

    // Update velocities based on current forces (half-step)
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
}

// Second step of the Verlet integration: update velocities
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    // Complete the velocity update with the new forces
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
}
