#include "global.h"
#include "def.h"
#include <cmath>
#include <algorithm>
#include <iostream>




double getNewTimeStepReg(double v[3], double df[3][4]) {

	double v2, F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     = df[0][0]*df[0][0] + df[1][0]*df[1][0] + df[2][0]*df[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];
	v2     = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];


	double dta;
  dta  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	dta  = std::sqrt(eta*dta);
	return dta; //dta;
}


double getNewTimeStepIrr(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


	/*
	else if  (F2 == 0 && F2dot2 != 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
	else if (F2 != 0 && F2dot2 == 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
													 */
	//fprintf(stdout, "NBODY+: f dots: %e, %e, %e, %e\n", F2, Fdot2, F2dot2, F3dot2);
	/*
	if (F2 != 0 && F2dot2 != 0 && dt != 0) {
		DivergentPrevent = std::sqrt(F2/Fdot2)/dt;
		fprintf(stdout, "F2 = %e, F2dot2 = %e, dt = %e\n", F2, F2dot2, dt);
		fprintf(stdout, "DivergentPrevent = %e\n", DivergentPrevent);
	}
	else {
		DivergentPrevent = 1;
	}
	*/

	TimeStep = (std::sqrt(F2 * F2dot2) + Fdot2) / (std::sqrt(Fdot2 * F3dot2) + F2dot2);
	// TimeStep  = F2/Fdot2;
	// TimeStep  = std::sqrt(DivergentPrevent*eta*TimeStep);
	TimeStep  = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

double getNewTimeStep(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


	/*
	else if  (F2 == 0 && F2dot2 != 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
	else if (F2 != 0 && F2dot2 == 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
													 */
	//fprintf(stdout, "NBODY+: f dots: %e, %e, %e, %e\n", F2, Fdot2, F2dot2, F3dot2);
	/*
	if (F2 != 0 && F2dot2 != 0 && dt != 0) {
		DivergentPrevent = std::sqrt(F2/Fdot2)/dt;
		fprintf(stdout, "F2 = %e, F2dot2 = %e, dt = %e\n", F2, F2dot2, dt);
		fprintf(stdout, "DivergentPrevent = %e\n", DivergentPrevent);
	}
	else {
		DivergentPrevent = 1;
	}
	*/

  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
  //TimeStep  = F2/Fdot2;
	//TimeStep  = std::sqrt(DivergentPrevent*eta*TimeStep);
	TimeStep  = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

void getBlockTimeStep(double dt, int &TimeLevel, ULL &TimeBlock, double &TimeStep) {
	TimeLevel = static_cast<int>(floor(log(dt/global_variable->EnzoTimeStep)/log(2.0)));
	//TimeLevel = static_cast<int>(ceil(log(dt/EnzoTimeStep)/log(2.0)));
	//std::cout << "NBODY+: TimeLevel = " << TimeLevel << std::endl;
	//std::cout << "NBODY+: TimeStep = " << TimeStep << std::endl;
	
	if (TimeLevel < global_variable->time_block) {
		//std::cerr << "TimeLevel is less than time block!!" << std::endl;
		TimeLevel = global_variable->time_block;
	}

	TimeStep = static_cast<double>(pow(2, TimeLevel));
	TimeBlock = static_cast<ULL>(pow(2, TimeLevel-global_variable->time_block));
}


void Particle::setNewTimeStepWithNewEnzoTimeStep(double &OldEnzoTimeStep, double &NewEnzoTimeStep) {
	double dt = this->TimeStepReg*OldEnzoTimeStep; //*1e10/1e6;
	this->TimeLevelReg = static_cast<int>(floor(log(dt / NewEnzoTimeStep) / log(2.0)));
	this->TimeLevelReg = std::min(0,this->TimeLevelReg);
	this->TimeStepReg = static_cast<double>(pow(2, this->TimeLevelReg));
	this->TimeBlockReg = static_cast<ULL>(pow(2, this->TimeLevelReg-global_variable->time_block));

	dt = this->TimeStepIrr*OldEnzoTimeStep; //*1e10/1e6;
	this->TimeLevelIrr = static_cast<int>(floor(log(dt / NewEnzoTimeStep) / log(2.0)));
	this->TimeLevelIrr = std::min(0,this->TimeLevelIrr);
	this->TimeStepIrr = static_cast<double>(pow(2, this->TimeLevelIrr));
	this->TimeBlockIrr = static_cast<ULL>(pow(2, this->TimeLevelIrr-global_variable->time_block));
}

