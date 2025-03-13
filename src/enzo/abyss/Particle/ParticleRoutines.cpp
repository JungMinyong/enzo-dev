#include "../global.h"
#include <cmath>

void Particle::set(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim], int &i) {

	__initialize__();
	this->isActive 					 = true;
	this->PID                        = PID[i];
	this->Mass                       = Mass[i];
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i];
	this->DynamicalTime              = DynamicalTime[i];
	this->Position[0]                = Position[0][i];
	this->Position[1]                = Position[1][i];
	this->Position[2]                = Position[2][i];
	this->Velocity[0]                = Velocity[0][i];
	this->Velocity[1]                = Velocity[1][i];
	this->Velocity[2]                = Velocity[2][i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i];
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i];
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i];
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->RadiusOfNeighbor           = InitialNeighborRadius2;
}

void Particle::set(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
                   double *Position[Dim], double *Velocity[Dim],
                   double *BackgroundAcceleration[Dim], int ParticleType, int &i)
{
    __initialize__();
	this->isActive 					 = true;
	this->PID                        = PID[i];
	this->Mass                       = Mass[i];
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i];
	this->DynamicalTime              = DynamicalTime[i];
	this->Position[0]                = Position[0][i];
	this->Position[1]                = Position[1][i];
	this->Position[2]                = Position[2][i];
	this->Velocity[0]                = Velocity[0][i];
	this->Velocity[1]                = Velocity[1][i];
	this->Velocity[2]                = Velocity[2][i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i];
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i];
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i];
	this->ParticleType               = ParticleType;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->RadiusOfNeighbor           = InitialNeighborRadius2;
}

void Particle::update(double *Mass, double *BackgroundAcceleration[Dim], int &i)
{
	this->Mass                       = Mass[i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i];
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i];
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i];
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}