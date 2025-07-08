#include "global.h"

int updateParticleBackground(){
    // Update position and velocity from background acceleration
    Particle *ptcl = nullptr;
    double dt = global_variable->time_step;

    for (int i = 0; i <= global_variable->LastParticleIndex; i++)
    {
        ptcl = &particles[i];
        if (!ptcl->isActive)
            continue;
        for (int dim = 0; dim < Dim; dim++)
        {
            ptcl->Position[dim] += (ptcl->BackgroundAcceleration[dim]*dt/2 + ptcl->Velocity[dim])*dt;
            ptcl->Velocity[dim] += ptcl->BackgroundAcceleration[dim]*dt;
        }
    }

    return 1;
}