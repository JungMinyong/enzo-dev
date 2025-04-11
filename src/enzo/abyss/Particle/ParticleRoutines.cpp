#include "../global.h"
#include <cmath>

void Particle::set(int *PID, double *Mass, double *CreationTime, double *DynamicalTime, double *Metallicity,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim], int &i) {

	__initialize__();
	this->isActive 					 = true;
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass*mass_unit; // [Msol unit] by EW 2025.4.3
	this->CreationTime               = CreationTime[i]*EnzoTime*1e4; // [Myr unit] by EW 2025.4.3
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->InitialMetallicity		 = Metallicity[i]; // [Absolute unit] No unit conversion here! by EW 2025.4.3
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]*EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]*EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]*EnzoAcceleration;
	this->RadiusOfNeighbor           = InitialNeighborRadius2; // fixed by EW 2025.3.12
	fprintf(stdout, "PID: %d. T_ini: %e Myr, M_ini: %e Msol, Z_ini: %e\n", this->PID, this->CreationTime, this->InitialMass, this->InitialMetallicity);
#ifndef SEVN
	this->ParticleType = NoFeedbackStar;
	this->radius = 2.25461e-8/position_unit*pow(this->Mass*mass_unit, 1./3); // stellar radius in code unit
	/*
	if (this->Mass*1e9 > 8) {
		this->ParticleType = Blackhole+SingleStar;
		this->radius = 2*this->Mass*1e9/mass_unit/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzshild radius in code unit
		// initialBHspin(this);
	}
	else {
		this->ParticleType = NormalStar+SingleStar;
		this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 1./3); // stellar radius in code unit
	}
	*/
#endif
}

void Particle::update(double *Mass, double *BackgroundAcceleration[Dim], int &i)
{
#ifndef SEVN
	this->Mass                       = Mass[i]*EnzoMass;
#endif
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]*EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]*EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]*EnzoAcceleration;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
#ifdef SEVN
	this->dm						 = 0.0;
	this->SNEjectedMass				 = 0.0;
	if (this->StellarEvolution != nullptr) {
		if (this->StellarEvolution->vkick[3] > 0.0) {
			fprintf(SEVNout, "PID: %d. Kicked velocity: (%e, %e, %e) [km/s]\n", this->PID, this->StellarEvolution->vkick[0], this->StellarEvolution->vkick[1], this->StellarEvolution->vkick[2]);
			for(int i=0; i<Dim; i++)
				this->Velocity[i] += this->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
			if (this->CMPtclIndex != -1)
				this->setBinaryInterruptState(BinaryInterruptState::kicked);
		}
	}
#endif
}