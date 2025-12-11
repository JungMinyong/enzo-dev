#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

struct GlobalVariable {
	int NumberOfSingleParticle;
	int LastParticleIndex;
	ULL NextRegTimeBlock;
	
	int time_block;
	double time_step;
	ULL block_max;
	double EnzoTimeStep;
	double OldEnzoTimeStep;
	double EnzoCurrentTime; // Current Enzo Time in Myr
	// int NumberOfParticle;
// #ifdef COMOVE
	double a_i = 1.0;
	double a_f = 1.0;
	double dadt_i = 0.0;
	double dadt_f = 0.0;
// #endif
};

#endif
