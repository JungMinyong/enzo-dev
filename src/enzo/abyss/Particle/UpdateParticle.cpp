#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include "../global.h"
#include "../def.h"

#define SMALL_TIMESTEP_TEST




#ifdef CUDA_FLOAT
void Particle::predictParticleSecondOrder(double dt, CUDA_REAL pos[], CUDA_REAL vel[]) {
	// Doubling check
	// temporary variables for calculation

	// only predict the positions if necessary
	// how about using polynomial correction here?
	
#ifdef COMOVE
	// (Query) Should I use CurrentTimeReg instead of CurrentTimeIrr if NumNeighbor == 0?
	double vel_mid;
	double a = getScaleFactor(dt, global_variable->a_i, global_variable->a_f); //global_variable->a_i + (global_variable->a_f-global_variable->a_i)*(CurrentTimeIrr+dt*0.5); //dt is range from 0 to 1
	double dadt = getScaleFactorDot(dt, global_variable->a_i, global_variable->a_f); //global_variable->dadt_i + (global_variable->dadt_f-global_variable->dadt_i)*(CurrentTimeIrr+dt*0.5);
	double H   = dadt / a;
	double a2  = a*a;
#endif

	dt = dt*global_variable->EnzoTimeStep;
	fprintf(nbpout, "dt = %e, a, dadt, H = %e, %e, %e\n", dt, a, dadt, H);

	if (dt == 0 || std::isnan(a_tot[0][0])) {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = (CUDA_REAL)Position[dim];
			vel[dim] = (CUDA_REAL)Velocity[dim];
		}
	}
	else {
		for (int dim=0; dim<Dim; dim++) {
#ifdef COMOVE
			// a_tot[dim] is obtained from comoving separation. So we have to devided it by a^2
			// vel_mid = (a_tot[dim][1]*dt/2 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt*0.5 + Velocity[dim];
			double acc = a_tot[dim][0] + BackgroundAcceleration[dim];// g/a^2
			double jerk= a_tot[dim][1]; // - H * a_tot[dim][0]/a2;       // j/a^2
			double dudt   = acc - H * Velocity[dim];   // total peculiar acceleration

			fprintf(nbpout, "acc = %e, jerk = %e, dudt = %e\n", acc, jerk, dudt);

			// 2nd‑order velocity
			vel[dim] = (CUDA_REAL) Velocity[dim] + (dudt + 0.5*jerk*dt)*dt;
			double xdot  = Velocity[dim]/a;
			double xddot = (dudt - H*Velocity[dim]) / a;   // from d(u/a)/dt
			pos[dim] = (CUDA_REAL) Position[dim] + xdot*dt + 0.5*xddot*dt*dt + jerk*dt*dt*dt/(6.0*a);
			// fprintf(stderr, "PID = %d, a_tot[dim][0] = %e, a_tot[dim][1] = %e, BackgroundAcceleration[dim] = %e, vel_mid = %e, dt = %e\n",  PID, a_tot[dim][0], a_tot[dim][1], BackgroundAcceleration[dim], vel_mid, dt);
			// fprintf(stderr, "PID = %d, a = %e, dadt = %e, dt = %e\n", PID, a, dadt, dt);
#else
			pos[dim] = (CUDA_REAL) ((a_tot[dim][1]*dt/3 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt/2 + Velocity[dim])*dt + Position[dim];
			vel[dim] = (CUDA_REAL)  (a_tot[dim][1]*dt/2 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt   + Velocity[dim];
#endif
		}
	}

	/*
	if (StarParticleFeedback != 0) {
		PredMass = Mass + evolveStarMass(CurrentTimeIrr, CurrentTimeIrr+TimeStepIrr*1.01);
		Mdot     = (PredMass - Mass)/TimeStepIrr*1e-2;
	}
	else {
		PredMass = Mass;
		Mdot     = 0;
	}
	*/
	return;
}
#endif

void Particle::predictParticleSecondOrder(double dt, double pos[], double vel[]) {
	// Doubling check
	// temporary variables for calculation

	// only predict the positions if necessary
	// how about using polynomial correction here?
#ifdef COMOVE
	// (Query) Should I use CurrentTimeReg instead of CurrentTimeIrr if NumNeighbor == 0?
	double vel_mid;
	double a = getScaleFactor(dt, global_variable->a_i, global_variable->a_f); //global_variable->a_i + (global_variable->a_f-global_variable->a_i)*(CurrentTimeIrr+dt*0.5); //dt is range from 0 to 1
	double dadt = getScaleFactorDot(dt, global_variable->a_i, global_variable->a_f); //global_variable->dadt_i + (global_variable->dadt_f-global_variable->dadt_i)*(CurrentTimeIrr+dt*0.5);
	double H   = dadt / a;
	double a2  = a*a;
#endif

	dt = dt*global_variable->EnzoTimeStep;

	if (dt == 0) {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = Position[dim];
			vel[dim] = Velocity[dim];
		}
	}
	else {
		for (int dim=0; dim<Dim; dim++) {
#ifdef COMOVE
			// a_tot[dim] is obtained from comoving separation. So we have to devided it by a^2
			// vel_mid = (a_tot[dim][1]*dt/2 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt*0.5 + Velocity[dim];
			double acc = a_tot[dim][0] + BackgroundAcceleration[dim];// g/a^2
			double jerk= a_tot[dim][1]; // - H * a_tot[dim][0]/a2;       // j/a^2
			double dudt   = acc - H * Velocity[dim];   // total peculiar acceleration

			// 2nd‑order velocity
			vel[dim] = Velocity[dim] + (dudt + 0.5*jerk*dt)*dt;
			double xdot  = Velocity[dim]/a;
			double xddot = (dudt - H*Velocity[dim]) / a;   // from d(u/a)/dt
			pos[dim] = Position[dim] + xdot*dt + 0.5*xddot*dt*dt + jerk*dt*dt*dt/(6.0*a);
			// fprintf(stderr, "PID = %d, a_tot[dim][0] = %e, a_tot[dim][1] = %e, BackgroundAcceleration[dim] = %e, vel_mid = %e, dt = %e\n",  PID, a_tot[dim][0], a_tot[dim][1], BackgroundAcceleration[dim], vel_mid, dt);
			// fprintf(stderr, "PID = %d, a = %e, dadt = %e, dt = %e\n", PID, a, dadt, dt);
#else
			pos[dim] = ((a_tot[dim][1] * dt / 3 + a_tot[dim][0] + BackgroundAcceleration[dim]) * dt / 2 + Velocity[dim]) * dt + Position[dim];
			vel[dim] = (a_tot[dim][1] * dt / 2 + a_tot[dim][0] + BackgroundAcceleration[dim]) * dt + Velocity[dim];
#endif
		}
	}
	return;
}


/*
 *  Purporse: Correct particle positions and velocities up to fourth order
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */
void Particle::correctParticleFourthOrder(double dt, double pos[], double vel[], double a[3][4]) {
	double dt3,dt4,dt5;

	dt = dt*global_variable->EnzoTimeStep;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<Dim; dim++) {
		NewPosition[dim] = pos[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		NewVelocity[dim] = vel[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}


/*
void Particle::polynomialPrediction(double current_time) {

}
*/


void Particle::updateParticle() {
	
	for (int dim=0; dim<Dim; dim++) {
		this->Position[dim] = this->NewPosition[dim];
		this->Velocity[dim] = this->NewVelocity[dim];
	}
	
	//updateTimeStep();
}



void Particle::updateRadius() {

	/* exponential (aggressive) */
	/*
		const double c = 0.5;
		const double b = std::log(2) / (NumNeighborMax);  // ln(2) / 40
		double exp = a * (std::exp(b * NumberOfAC) - 1);
		*/

	/* n=2 polynomial (mild) as n increases it grows mild */

	if (this->NumberOfNeighbor > FixNumNeighbor) {
		const int n = 2;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.9 / std::pow(c,n);  // ln(2) / 40
		double x = this->NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->RadiusOfNeighbor *= (1.-a);
	}
	else if (this->NumberOfNeighbor < FixNumNeighbor) {
		const int n = 3;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.5 / std::pow(c,n);  // ln(2) / 40
		double x = this->NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->RadiusOfNeighbor *= (1.-a);
	}
}







double getNewTimeStepReg(double v[3], double f[3][4]);
double getNewTimeStepIrr(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ULL &TimeBlock, double &TimeStep);

void Particle::calculateTimeStepIrr() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ULL TimeBlockTmp;


#ifdef SMALL_TIMESTEP_TEST
	fprintf(stderr, "In calculateTimeStepyIrr, before SMALL_TIMESTEP_TEST, TimeLevelIrr=%d, TimeStepIrr=%e\n", TimeLevelIrr, TimeStepIrr);
	TimeStepIrr  = 0.25; // 0.5 Myr
	TimeLevelIrr  = static_cast<int>(log(TimeStepIrr)/log(2.)); // 0.5 Myr
	TimeBlockIrr = static_cast<ULL>(pow(2,TimeLevelIrr-global_variable->time_block));
	return;
#endif


	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
		TimeStepIrr = static_cast<double>(pow(2, TimeLevelReg));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelReg-global_variable->time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (NextRegTimeBlock < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, NextRegTimeBlock=%llu, TimeBlockReg=%llu\n", PID, NextRegTimeBlock, TimeBlockReg);
	}
	*/

	if (TimeLevelTmp > TimeLevelIrr) {
		if (fmod(NewCurrentBlockIrr, 2*TimeBlockIrr)==0) {
			TimeLevelTmp = TimeLevelIrr+1;
			TimeBlockTmp = 2*TimeBlockIrr;
		}
		else {
			TimeLevelTmp = TimeLevelIrr;
			TimeBlockTmp = TimeBlockIrr;
		}
	}
	else if (TimeLevelTmp < TimeLevelIrr) {
		if (TimeLevelTmp < TimeLevelIrr-1) {
			TimeLevelTmp = TimeLevelIrr - 2;
			TimeBlockTmp = TimeBlockIrr/4;
		}
		else {
			TimeLevelTmp = TimeLevelIrr - 1;
			TimeBlockTmp = TimeBlockIrr/2;
		}
	} else {
		TimeLevelTmp = TimeLevelIrr;
		TimeBlockTmp = TimeBlockIrr;
	}


	while (((NewCurrentBlockIrr < CurrentBlockReg+TimeBlockReg) && (NewCurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg)) || (TimeLevelTmp >= TimeLevelReg)) {
		/*
		fprintf(stderr,"CurrentBlockIrr = %llu\n",
				CurrentBlockIrr);
		fprintf(stderr,"CurrentBlockReg = %llu\n",
				CurrentBlockReg);
		fprintf(stderr,"TimeBlockReg    = %llu\n",
				TimeBlockReg);
		fprintf(stderr,"TimeBlockIrr    = %llu\n",
				TimeBlockIrr);
				*/
		//if (TimeLevelTmp > TimeLevelReg)
			//fprintf(stderr, "PID=%d, Irr=%d, Reg=%d\n", PID, TimeLevelTmp0, TimeLevelReg);

		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < global_variable->time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		TimeLevelIrr = std::max(global_variable->time_block, TimeLevelIrr);
	}


	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));

	if (TimeStepIrr*global_variable->EnzoTimeStep*1e4<1e-10) {
		fprintf(stderr, "Too small TimeStepIrr! PID: %d (NN: %d, rad: %e pc), TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, NumberOfNeighbor, sqrt(RadiusOfNeighbor)*position_unit,
				TimeStepIrr*global_variable->EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*global_variable->EnzoTimeStep*1e4);
		while (TimeStepIrr*global_variable->EnzoTimeStep*1e4<1e-10) {
			TimeLevelIrr++;
			TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
			TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));
		}
	}

	if (TimeStepIrr > 1) {
		fprintf(stderr, "Why TimeStepIrr is too large? PID: %d, pos: (%e, %e, %e), vel: (%e, %e, %e)\n",
				PID, Position[0], Position[1], Position[2], Velocity[0], Velocity[1], Velocity[2]);
		fprintf(stderr, "airr0: (%e, %e, %e), airr1: (%e, %e, %e), airr2: (%e, %e, %e), airr3: (%e, %e, %e)\n",
				a_irr[0][0], a_irr[1][0], a_irr[2][0], a_irr[0][1], a_irr[1][1], a_irr[2][1],
				a_irr[0][2], a_irr[1][2], a_irr[2][2], a_irr[0][3], a_irr[1][3], a_irr[2][3]);
		fprintf(stderr, "areg0: (%e, %e, %e), areg1: (%e, %e, %e), areg2: (%e, %e, %e), areg3: (%e, %e, %e)\n",
				a_reg[0][0], a_reg[1][0], a_reg[2][0], a_reg[0][1], a_reg[1][1], a_reg[2][1],
				a_reg[0][2], a_reg[1][2], a_reg[2][2], a_reg[0][3], a_reg[1][3], a_reg[2][3]);
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*EnzoTimeStep*1e4 < KSTime && (this->isCMptcl == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}

void Particle::calculateTimeStepIrr2() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ULL TimeBlockTmp;

#ifdef SMALL_TIMESTEP_TEST
fprintf(stderr, "In calculateTimeStepIrr2, before SMALL_TIMESTEP_TEST, TimeLevelIrr=%d, TimeStepIrr=%e\n", TimeLevelIrr, TimeStepIrr);
	TimeStepIrr  = 0.25; // 0.5 Myr
	TimeLevelIrr  = static_cast<int>(log(TimeStepIrr)/log(2.)); // 0.5 Myr
	TimeBlockIrr = static_cast<ULL>(pow(2,TimeLevelIrr-global_variable->time_block));
	return;
#endif

	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
		TimeStepIrr = static_cast<double>(pow(2, TimeLevelReg));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelReg-global_variable->time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (NextRegTimeBlock < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, NextRegTimeBlock=%llu, TimeBlockReg=%llu\n", PID, NextRegTimeBlock, TimeBlockReg);
	}
	*/


	while (((NewCurrentBlockIrr < CurrentBlockReg+TimeBlockReg) && (NewCurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg)) || (TimeLevelTmp >= TimeLevelReg)) {
		/*
		fprintf(stderr,"CurrentBlockIrr = %llu\n",
				CurrentBlockIrr);
		fprintf(stderr,"CurrentBlockReg = %llu\n",
				CurrentBlockReg);
		fprintf(stderr,"TimeBlockReg    = %llu\n",
				TimeBlockReg);
		fprintf(stderr,"TimeBlockIrr    = %llu\n",
				TimeBlockIrr);
				*/
		//if (TimeLevelTmp > TimeLevelReg)
			//fprintf(stderr, "PID=%d, Irr=%d, Reg=%d\n", PID, TimeLevelTmp0, TimeLevelReg);

		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < global_variable->time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		TimeLevelIrr = std::max(global_variable->time_block, TimeLevelIrr);
	}


	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));

	if (TimeStepIrr*global_variable->EnzoTimeStep*1e4<1e-10) {
		fprintf(stderr, "Too small TimeStepIrr! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, TimeStepIrr*global_variable->EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*global_variable->EnzoTimeStep*1e4);
		exit(1);
		while (TimeStepIrr*global_variable->EnzoTimeStep*1e4<1e-10) {
			TimeLevelIrr++;
			TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
			TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));
		}
	}

	if (TimeStepIrr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*EnzoTimeStep*1e4 < KSTime && (this->isCMptcl == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}



// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

#ifdef SMALL_TIMESTEP_TEST
	TimeStepReg  = 0.5; 
	TimeLevelReg  = static_cast<int>(log(TimeStepReg)/log(2.)); // 0.5 Myr
	TimeBlockReg = static_cast<ULL>(pow(2,TimeLevelReg-global_variable->time_block));
	fprintf(stderr, "In calculateTimeStepReg, before SMALL_TIMESTEP_TEST, TimeLevelReg=%d, TimeStepReg=%e\n", TimeLevelReg, TimeStepReg);
	return;
#endif

	getBlockTimeStep(getNewTimeStepReg(Velocity, a_reg), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);


	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepTmp << std::endl;

	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp >= TimeLevelReg+1) {
		if (fmod(CurrentBlockReg, 2*TimeBlockReg)==0 \
				&& CurrentTimeReg != 0) {
			TimeLevelTmp   = TimeLevelReg+1;
			TimeBlockTmp   = TimeBlockReg*2;

			while ((TimeLevelTmp0 > TimeLevelTmp) \
					&& (fmod(CurrentBlockReg, 2*TimeBlockTmp) == 0)) {
				TimeBlockTmp = 2*TimeBlockTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeLevelTmp < TimeLevelReg) {
		TimeLevelTmp = TimeLevelReg - 1;
		if (TimeLevelTmp0 < TimeLevelTmp)
			TimeLevelTmp--;
	}
	else {
		TimeLevelTmp = TimeLevelReg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
	/*
	if (TimeStepRegTmp > 0.1 && TimeStepRegTmp > TimeStepReg) {
		double v2 = 0., a2=0., dt;
		for (int dim=0; dim<Dim; dim++) {
			v2 += (PredVelocity[dim]-NewVelocity[dim])*(PredVelocity[dim]-NewVelocity[dim]);
			a2 += a_reg[dim][0]*a_reg[dim][0];
		}
		dt = TimeStepReg*std::pow((1e-4*TimeStepReg*TimeStepReg*a2/v2),0.1);
		if (dt < TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg;
		}	
	}
	*/

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);

	TimeLevelReg = std::max(global_variable->time_block,TimeLevelTmp);
	//TimeLevelReg = std::max(time_block, TimeLevelReg);

	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
	}

	TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
	TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-global_variable->time_block));

	if (TimeStepReg*global_variable->EnzoTimeStep*1e4 < 1e-7) {
		fprintf(stderr, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
				PID, TimeStepReg * global_variable->EnzoTimeStep * 1e4, static_cast<double>(pow(2, TimeLevelTmp0)) * global_variable->EnzoTimeStep * 1e4);
		fflush(stderr);
	}

	if (CurrentTimeReg+TimeStepReg > 1 && CurrentTimeReg != 1.0) {
		TimeStepReg = 1 - CurrentTimeReg;
		TimeBlockReg = global_variable->block_max-CurrentBlockReg;
	}
	/*
	if (TimeStepReg*EnzoTimeStep*1e4<1e-9) {
		fprintf(stderr, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
				PID, TimeStepReg*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		throw std::runtime_error("TimeStepReg is too small.");
	}
	*/
	if (TimeStepReg*global_variable->EnzoTimeStep*1e4<1e-7) {
		fprintf(stderr, "Too small TimeStepReg! PID: %d (NN: %d, rad: %e pc), TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, NumberOfNeighbor, sqrt(RadiusOfNeighbor)*position_unit, 
				TimeStepReg*global_variable->EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*global_variable->EnzoTimeStep*1e4);
		while (TimeStepReg*global_variable->EnzoTimeStep*1e4<1e-7) {
			TimeLevelReg++;
			TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
			TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-global_variable->time_block));
		}
	}
	/* // original code
	if (TimeStepReg > 1) {
		fprintf(stderr, "TimeStepReg=%e, TimeLevelReg=%d, TimeLevelTmp0=%d\n",TimeStepReg, TimeLevelReg, TimeLevelTmp0);
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}
	*/
	// /* // modified version by EW 2025.5.12
	while (TimeStepReg > 1) {
		fprintf(stderr, "Too large TimeStepReg! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, TimeStepReg*global_variable->EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*global_variable->EnzoTimeStep*1e4);
		TimeLevelReg--;
		TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
		TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-global_variable->time_block));
	}
	// */


	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}


// Use this function when OnlyIrregularRoutines is true (RadiusOfNeighbor == 1e20)
void Particle::calculateTimeStepOnlyIrr() {

#ifdef SMALL_TIMESTEP_TEST
	TimeStepIrr  = 0.25; // 0.5 Myr
	TimeLevelIrr  = static_cast<int>(log(TimeStepIrr)/log(2.)); // 0.5 Myr
	TimeBlockIrr = static_cast<ULL>(pow(2,TimeBlockIrr-global_variable->time_block));
	fprintf(stderr, "In calculateTimeStepOnlyIrr, before SMALL_TIMESTEP_TEST, TimeLevelIrr=%d, TimeStepIrr=%e\n", TimeLevelIrr, TimeStepIrr);
	return;
#endif

	assert(this->NumberOfNeighbor > 0);
	assert(this->RadiusOfNeighbor == 1e20);

	double TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStepReg(Velocity, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp > TimeLevelIrr+1) {
		if (fmod(NewCurrentBlockIrr, 2*TimeBlockIrr)==0) {
			TimeLevelTmp = TimeLevelIrr+1;
			TimeBlockTmp = 2*TimeBlockIrr;
		}
		else {
			TimeLevelTmp = TimeLevelIrr;
			TimeBlockTmp = TimeBlockIrr;
		}
	}
	else if (TimeLevelTmp < TimeLevelIrr) {
		if (TimeLevelTmp < TimeLevelIrr-1) {
			TimeLevelTmp = TimeLevelIrr - 2;
			TimeBlockTmp = TimeBlockIrr/4;
		}
		else {
			TimeLevelTmp = TimeLevelIrr - 1;
			TimeBlockTmp = TimeBlockIrr/2;
		}
	} else {
		TimeLevelTmp = TimeLevelIrr;
		TimeBlockTmp = TimeBlockIrr;
	}

	TimeLevelIrr = std::max(global_variable->time_block,TimeLevelTmp);
	TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));

	if (CurrentTimeIrr+TimeStepIrr > 1 && CurrentTimeIrr != 1.0) {
		TimeStepIrr = 1 - CurrentTimeIrr;
		TimeBlockIrr = global_variable->block_max-CurrentBlockIrr;
	}

	while (TimeStepIrr > 1) {
		TimeLevelIrr--;
		TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));
	}

}



// Use this function in FBInitialization when OnlyIrregularRoutines is true (RadiusOfNeighbor == 1e20)
void Particle::calculateTimeStepOnlyIrr2() {

#ifdef SMALL_TIMESTEP_TEST
	TimeStepIrr  = 0.25; // 0.5 Myr
	TimeLevelIrr  = static_cast<int>(log(TimeStepIrr)/log(2.)); // 0.5 Myr
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeBlockIrr-global_variable->time_block));
fprintf(stderr, "In calculateTimeStepOnlyIrr2, before SMALL_TIMESTEP_TEST, TimeLevelIrr=%d, TimeStepIrr=%e\n", TimeLevelIrr, TimeStepIrr);
return;
#endif

	assert(this->NumberOfNeighbor > 0);
	assert(this->RadiusOfNeighbor == 1e20);

	double TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStepReg(Velocity, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	TimeLevelIrr = std::max(global_variable->time_block,TimeLevelTmp);
	TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));

	if (CurrentTimeIrr+TimeStepIrr > 1 && CurrentTimeIrr != 1.0) {
		TimeStepIrr = 1 - CurrentTimeIrr;
		TimeBlockIrr = global_variable->block_max-CurrentBlockIrr;
	}

	while (TimeStepIrr > 1) {
		TimeLevelIrr--;
		TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-global_variable->time_block));
	}


}


