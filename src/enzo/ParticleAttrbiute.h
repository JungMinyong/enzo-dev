#include <iostream>
#include <string>



struct ParticleAttribute {

	double creation_time; // ParticleAttributes[0] in the old version
	double dynamical_time; // ParticleAttributes[1] in the old version
	double life_time;
	double mass_p;	 // ParticleAttributes[2] in the old version
	double metallicity;	
	double birth_mass;	
	double dummy[5];



#ifdef INDIVIDUALSTARS
	double se_table_position[2];
	double rad_table_position[3];
	double yield_table_position[2];
	double wind_mass_ejected;
	double sn_mass_ejected;
	double effective_temp;
#endif
#ifdef NBODY
	double background_acceleration[3];
	int AbyssParticleType;
#endif
#ifdef SEVN

#endif

	// Constructor
	//MyStruct(T val, const std::string& n) : value(val), name(n) {}

	// Member function
	//void print() const {
		//std::cout << "Name: " << name << ", Value: " << value << std::endl;
	//}
};
