#include "../global.h"
#include <cmath>

void Particle::set(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim], int &i) {

	__initialize__();
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	//this->RadiusOfNeighbor           = InitialRadiusOfAC;
}

void Particle::set(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
                   double *Position[Dim], double *Velocity[Dim],
                   double *BackgroundAcceleration[Dim], int ParticleType, int &i)
{
    __initialize__();
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType               = ParticleType;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	//this->RadiusOfNeighbor           = InitialRadiusOfAC;
}

void Particle::update(double *Mass, double *BackgroundAcceleration[Dim], int &i)
{
	this->Mass                       = Mass[i]*EnzoMass;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}