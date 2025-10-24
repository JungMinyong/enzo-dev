#ifdef SEVN

#include "global.h"
#include <random>
#include <map>
#include "supernova.h"

void UpdateEvolution(Particle* ptcl);
void setBHspin(Particle* ptcl);
void adjustEnzoSEVNParticle(int index, bool delete_sevn);

// Start SEVN stellar evolution
// Set stellar radii, BH spin, etc
void initializeStellarEvolution(double *MHECurrent, double *MCOCurrent, double *SEVNFrac, double *SEVNType, double *WorldTime) {

    std::ostringstream oss;
	oss << sevnio->svpar;
	fprintf(SEVNout, "%s\n\n\n\n\n", oss.str().c_str());
    // fprintf(SEVNout, "PID\tPhase\tMass (Msun)\tRadius (pc)\tTime (Myr)\tWorldtime (Myr)\n")
	fflush(SEVNout);

    assert(SEVNList.empty());

    for (int i=0; i<NumberOfSingleParticle; i++) {

        Particle* ptcl = &particles[i];

        ptcl->WorldTime = global_variable->EnzoCurrentTime;

        if (ptcl->CreationTime < 0.0) {
            ptcl->radius = 2.25461e-8/position_unit*pow(ptcl->Mass*mass_unit, 1/3);
            ptcl->ParticleType = NO_FEEDBACK_STAR;
            continue;
        }

        if (ptcl->InitialMass < SEVNLowerMassLimit) {
			ptcl->radius = 2.25461e-8/position_unit*pow(ptcl->Mass*mass_unit, 1/3);
            ptcl->ParticleType = NO_FEEDBACK_STAR;
            // stellar radius in code unit
            // extrapolated from the solar radius (it is assumed that low-mass stars have the same stellar density to the sun)

			// fprintf(SEVNout, "PID: %d. Mass: %e Msol, Radius: %e pc\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->radius*position_unit);
			continue;
		}

        ptcl->ParticleType = static_cast<int>(SEVNType[i]);
        ptcl->WorldTime = WorldTime[i];

        // Mass, metallicity, spin, sn model, tini, tf, dtout, random seed(optional)
        double Metallicity = ptcl->InitialMetallicity > SEVNMetallicityUpperLimit ? SEVNMetallicityUpperLimit : ptcl->InitialMetallicity;
        Metallicity = Metallicity < SEVNMetallicityLowerLimit ? SEVNMetallicityLowerLimit : Metallicity;

        std::stringstream Mass;
        std::stringstream tini;
        if (ptcl->ParticleType > 6) {

            std::string remnant_suffix;

            if (ptcl->ParticleType == 7 || ptcl->ParticleType == 8) {
                fprintf(stderr, "Error: Particle %d is already a remnant (type %d) at birth! Setting to NO_FEEDBACK_STAR.\n", ptcl->PID, ptcl->ParticleType);
                ptcl->ParticleType = NO_FEEDBACK_STAR;
                continue;
            }

            if (ptcl->ParticleType == 9)
                remnant_suffix = "HEWD";
            else if (ptcl->ParticleType == 10)
                remnant_suffix = "COWD";
            else if (ptcl->ParticleType == 11)
                remnant_suffix = "ONEWD";
            else if (ptcl->ParticleType == 12)
                remnant_suffix = "NSEC";
            else if (ptcl->ParticleType == 13)
                remnant_suffix = "NS";
            else if (ptcl->ParticleType == 14 || ptcl->ParticleType == 15)
                remnant_suffix = "BH";

            Mass << std::setprecision(17) << ptcl->Mass*mass_unit << remnant_suffix;

            tini << "zams";
        } else {
            Mass << std::setprecision(17) << "(" << ptcl->InitialMass << ","
                << ptcl->Mass*mass_unit << ","
                << MHECurrent[i] << "," 
                << MCOCurrent[i] << ")";

            tini << std::setprecision(17) << "%" << SEVNFrac[i] * 100
                << ":" << ptcl->ParticleType;
        }

        std::vector<std::string> init_params{Mass.str(), std::to_string(Metallicity), "0.0", "delayed", tini.str(), "end", "events"};

        size_t id = ptcl->PID;
        ptcl->StellarEvolution = new StarSEVN(sevnio, init_params, id, false);
        fprintf(stderr, "New star in SEVN! PID: %d, Mass: %e Msun, Z: %e\n", ptcl->PID, ptcl->Mass * mass_unit, ptcl->InitialMetallicity);
        // fprintf(stdout, "New star in SEVN! PID: %d, Initial mass: %e Msun, Z: %e\n", ptcl->PID, ptcl->InitialMass, ptcl->InitialMetallicity);

        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // stellar radius in code unit
        ptcl->T_eff = ptcl->StellarEvolution->getp(Temperature::ID);

        if (ptcl->StellarEvolution->amiremnant()) {
            UpdateEvolution(ptcl);
            ptcl->dm = 0.0;
            ptcl->SNEjectedMass = 0.0;
        }
        else
            SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
    }
}

void initializeStellarEvolution(int ParticleIndex) {

    Particle* ptcl = &particles[ParticleIndex];

    ptcl->WorldTime = global_variable->EnzoCurrentTime;
    ptcl->ParticleType = NO_FEEDBACK_STAR;

    if (PIDtoIndexMap_SEVN.find(ptcl->PID) != PIDtoIndexMap_SEVN.end()) {
        int index_SEVN = PIDtoIndexMap_SEVN[ptcl->PID];

        ptcl->StellarEvolution = SEVNList_Enzo[index_SEVN];
        // ptcl->CreationTime = creation_time_Enzo[index_SEVN]; // These two should be the same by EW 2025.5.1
        ptcl->WorldTime = world_time_Enzo[index_SEVN];
        ptcl->InitialMass = ptcl->StellarEvolution->get_zams();

        adjustEnzoSEVNParticle(index_SEVN, false);

        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // stellar radius in code unit
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*mass_unit/ptcl->StellarEvolution->get_max_zams(), 1./3);

        if (ptcl->StellarEvolution->amiremnant()) {
            ptcl->ParticleType = REMNANT + (int)ptcl->StellarEvolution->getp(RemnantType::ID);
            if (ptcl->StellarEvolution->amiBH())
                setBHspin(ptcl);
            if (ptcl->Mass*mass_unit > MassiveBlackHoleCutoff) {
                ptcl->ParticleType = MASSIVE_BLACK_HOLE;
            }
        }
        else
            ptcl->ParticleType = (int)ptcl->StellarEvolution->getp(Phase::ID);
        return;
    }

    double Metallicity = ptcl->InitialMetallicity > SEVNMetallicityUpperLimit ? SEVNMetallicityUpperLimit : ptcl->InitialMetallicity;
    Metallicity = Metallicity < SEVNMetallicityLowerLimit ? SEVNMetallicityLowerLimit : Metallicity;
    std::vector<std::string> init_params{std::to_string(double(ptcl->InitialMass)), std::to_string(Metallicity), "0.0", "delayed", "zams", "end", "events"};

    size_t id = ptcl->PID;
    ptcl->StellarEvolution = new StarSEVN(sevnio, init_params, id, false);
    fprintf(stderr, "New star in SEVN! PID: %d, Initial mass: %.2e Msun, Z: %e\n", ptcl->PID, ptcl->InitialMass, ptcl->InitialMetallicity);
    fprintf(stdout, "New star in SEVN! PID: %d, Initial mass: %.2e Msun, Z: %e\n", ptcl->PID, ptcl->InitialMass, ptcl->InitialMetallicity);

    SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});

    ptcl->ParticleType = static_cast<int>(ptcl->StellarEvolution->getp(Phase::ID));
    ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // stellar radius in code unit
    ptcl->T_eff = ptcl->StellarEvolution->getp(Temperature::ID);
}

void setBHspin(Particle* ptcl) {
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    double phi = 2 * M_PI * distr(mt);
    double theta = M_PI * distr(mt);

    ptcl->a_spin[0] = ptcl->StellarEvolution->getp(Xspin::ID) * sin(theta)*cos(phi);
    ptcl->a_spin[1] = ptcl->StellarEvolution->getp(Xspin::ID) * sin(theta)*sin(phi);
    ptcl->a_spin[2] = ptcl->StellarEvolution->getp(Xspin::ID) * cos(theta);
}

void StellarEvolution() {

    Particle* ptcl;
    while (!SEVNList.empty()) {
         
        auto it = SEVNList.begin();
        ptcl = &particles[it->second];

        ptcl->WorldTime += ptcl->StellarEvolution->getp(Timestep::ID);
        ptcl->StellarEvolution->evolve();
        if (!ptcl->StellarEvolution->amiremnant())
            ptcl->T_eff = ptcl->StellarEvolution->getp(Temperature::ID);

        while (ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID) <= global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime) {
            ptcl->WorldTime += ptcl->StellarEvolution->getp(Timestep::ID);
            ptcl->StellarEvolution->evolve();
            if (!ptcl->StellarEvolution->amiremnant())
                ptcl->T_eff = ptcl->StellarEvolution->getp(Temperature::ID);
        }

        it = SEVNList.erase(it);
        UpdateEvolution(ptcl);

        if (SEVNList.empty() || SEVNList.begin()->first > global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime)
            break;
    }
    fflush(SEVNout);
}

void StellarEvolution_Enzo(double *Mass, double *Wind, double *SN, double *Temperature,
                            double *newMass, double *newWind, double *newSN, double *newTemperature) {


    if (NumberOfEnzoSEVNParticle == 0 && newNumberOfEnzoSEVNParticle == 0)
        return;

    int index;
    StarSEVN* star_sevn;
    bool evolved;
    double CurrentMass;
    double WorldTime;
    if (NumberOfEnzoSEVNParticle - newNumberOfEnzoSEVNParticle > 0) {
        for (int i=0; i<NumberOfEnzoSEVNParticle - newNumberOfEnzoSEVNParticle; i++) {
            index = PIDtoIndexMap_SEVN[EnzoPIDs_SEVN[i]];
            star_sevn = SEVNList_Enzo[index];

            if (star_sevn->amiremnant()) {
                if (star_sevn->amiempty()) {
                    Mass[i] = -1.0;
                    adjustEnzoSEVNParticle(index, true);
                }
                else
                    Mass[i] = star_sevn->getp(Mass::ID)/mass_unit/EnzoMass;
                Wind[i] = 0.0;
                SN[i] = 0.0;
                Temperature[i] = 0.0;
            }
            else {
                evolved = false;
                CurrentMass = star_sevn->getp(Mass::ID);
                WorldTime = world_time_Enzo[index];
                Mass[i] = CurrentMass/mass_unit/EnzoMass;
                Wind[i] = 0.0;
                SN[i] = 0.0;
                Temperature[i] = star_sevn->getp(Temperature::ID);
                while (WorldTime + star_sevn->getp(Timestep::ID) <= global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime) {
                    evolved = true;
                    WorldTime += star_sevn->getp(Timestep::ID);
                    star_sevn->evolve();

                    if (!star_sevn->amiremnant())
                        Temperature[i] = star_sevn->getp(Temperature::ID);
                    else
                        break;
                }
                world_time_Enzo[index] = WorldTime;
                if (evolved) {

                    fprintf(SEVNout, "(Enzo) PID: %d, Phase: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, T_eff: %e K, Worldtime: %e Myr, CurrentTime: %e Myr\n", 
                        EnzoPIDs_SEVN[i], int(star_sevn->getp(Phase::ID)), star_sevn->getp(Mass::ID), star_sevn->get_zams(), star_sevn->get_Z(),
                        star_sevn->getp(Radius::ID)/(utilities::parsec_to_Rsun), star_sevn->getp(Temperature::ID), star_sevn->getp(Worldtime::ID), 
                        global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime);

                    Mass[i] = star_sevn->getp(Mass::ID)/mass_unit/EnzoMass;
                    if (star_sevn->amiempty()) {
                        Mass[i] = -1.0;
                        SN[i] = CurrentMass;
                        adjustEnzoSEVNParticle(index, true);
                    }
                    else if (star_sevn->amiremnant()) {
                        if (star_sevn->get_supernova()->get_fallback_frac() == 1.0) // direct collapse BH
                            SN[i] = 0.0;
                        else
                            SN[i] = star_sevn->get_supernova()->get_Mejected();
                        Wind[i] = CurrentMass - star_sevn->getp(Mass::ID) - SN[i];
                    }
                    else {
                        Wind[i] = CurrentMass - star_sevn->getp(Mass::ID);
                        SN[i] = 0.0;
                    }
                    fprintf(stderr, "(Enzo) Feedback info send to Enzo...\n");
                    fprintf(stderr, "\tPID: %d. WindEjectedMass: %e Msun, SNEjectedMass: %e Msun, T_eff: %e K\n", 
                            EnzoPIDs_SEVN[i], Wind[i], SN[i], Temperature[i]);
                }
            }
        }
    }
    if (newNumberOfEnzoSEVNParticle > 0) {
        int offset = NumberOfEnzoSEVNParticle - newNumberOfEnzoSEVNParticle;
        for (int i=0; i<newNumberOfEnzoSEVNParticle; i++) {
            index = PIDtoIndexMap_SEVN[EnzoPIDs_SEVN[i+offset]];
            star_sevn = SEVNList_Enzo[index];

            evolved = false;
            CurrentMass = star_sevn->getp(Mass::ID);
            WorldTime = world_time_Enzo[index];
            newMass[i] = CurrentMass/mass_unit/EnzoMass;
            newWind[i] = 0.0;
            newSN[i] = 0.0;
            newTemperature[i] = star_sevn->getp(Temperature::ID);
            if (star_sevn->amiremnant()) {
                fprintf(stderr, "Why new EnzoSEVN particle (PID: %d) is remnant?\n", EnzoPIDs_SEVN[i+offset]);
                assert(!star_sevn->amiremnant());
            }
            while (WorldTime + star_sevn->getp(Timestep::ID) <= global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime) {
                evolved = true;
                WorldTime += star_sevn->getp(Timestep::ID);
                star_sevn->evolve();
                if (!star_sevn->amiremnant())
                    newTemperature[i] = star_sevn->getp(Temperature::ID);
                else
                    break;
            }
            world_time_Enzo[index] = WorldTime;
            if (evolved) {

                fprintf(SEVNout, "(Enzo) PID: %d, Phase: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, T_eff: %e K, Worldtime: %e Myr, CurrentTime: %e Myr\n", 
                        EnzoPIDs_SEVN[i+offset], int(star_sevn->getp(Phase::ID)), star_sevn->getp(Mass::ID), star_sevn->get_zams(), star_sevn->get_Z(),
                        star_sevn->getp(Radius::ID)/(utilities::parsec_to_Rsun), star_sevn->getp(Temperature::ID), star_sevn->getp(Worldtime::ID), 
                        global_time * global_variable->EnzoTimeStep * 1e4 + global_variable->EnzoCurrentTime);

                newMass[i] = star_sevn->getp(Mass::ID)/mass_unit/EnzoMass;
                if (star_sevn->amiempty()) {
                    newMass[i] = -1.0;
                    newSN[i] = CurrentMass;
                    adjustEnzoSEVNParticle(index, true);
                }
                else if (star_sevn->amiremnant()) {
                    if (star_sevn->get_supernova()->get_fallback_frac() == 1.0) // direct collapse BH
                        newSN[i] = 0.0;
                    else
                        newSN[i] = star_sevn->get_supernova()->get_Mejected();
                    newWind[i] = CurrentMass - star_sevn->getp(Mass::ID) - newSN[i];
                }
                else {
                    newWind[i] = CurrentMass - star_sevn->getp(Mass::ID);
                    newSN[i] = 0.0;
                }
                fprintf(stderr, "(Enzo) Feedback info send to Enzo...\n");
                fprintf(stderr, "\tPID: %d. WindEjectedMass: %e Msun, SNEjectedMass: %e Msun, T_eff: %e K\n", 
                        EnzoPIDs_SEVN[i+offset], newWind[i], newSN[i], newTemperature[i]);
            }
        }
    }
}

void UpdateEvolution(Particle* ptcl) {

    if (ptcl->StellarEvolution->amiremnant())
        ptcl->ParticleType = REMNANT + (int)ptcl->StellarEvolution->getp(RemnantType::ID);
    else
        ptcl->ParticleType = (int)ptcl->StellarEvolution->getp(Phase::ID);

    if (!ptcl->StellarEvolution->amiremnant()) {
        SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*mass_unit/ptcl->StellarEvolution->get_max_zams(), 1./3);
        fprintf(SEVNout, "PID: %d, Phase: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, T_eff: %e K, Time: %e Myr, Worldtime: %e Myr\n", 
            ptcl->PID, int(ptcl->StellarEvolution->getp(Phase::ID)), ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->StellarEvolution->get_Z(),
            ptcl->radius*position_unit, ptcl->T_eff, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
    }
    else if (ptcl->StellarEvolution->amiWD()) {

        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "WD. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", 
            ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->StellarEvolution->get_Z(),
            ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        /* // WD is not kicked in SEVN by EW 2025.4.1
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        */
    }
    else if (ptcl->StellarEvolution->amiNS()) {

        ptcl->SNEjectedMass = ptcl->StellarEvolution->get_supernova()->get_Mejected()/mass_unit;
        ptcl->dm += ptcl->Mass - (ptcl->StellarEvolution->getp(Mass::ID)/mass_unit + ptcl->SNEjectedMass);
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "NS. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", 
            ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->StellarEvolution->get_Z(),
            ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        /* // (SEVN Query) This particle will be kicked after SN feedback in Enzo by EW 2025.4.1
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        */
    }
    else if (ptcl->StellarEvolution->amiBH()) {

        setBHspin(ptcl);

        // Direct collapse case! newly implemented by EW 2025.5.21
        if (ptcl->StellarEvolution->get_supernova()->get_fallback_frac() == 1.0)
            ptcl->SNEjectedMass = 0.0;
        else // CCSN case!
            ptcl->SNEjectedMass = ptcl->StellarEvolution->get_supernova()->get_Mejected()/mass_unit;

        ptcl->dm += ptcl->Mass - (ptcl->StellarEvolution->getp(Mass::ID)/mass_unit + ptcl->SNEjectedMass);
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        if (ptcl->Mass*mass_unit > MassiveBlackHoleCutoff)
            ptcl->ParticleType = MASSIVE_BLACK_HOLE; // This particle does feedback & accretion in Enzo by EW 2025.6.25
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "BH. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Z: %e, Radius: %e pc, Fallback_frac: %e, Time: %e Myr, Worldtime: %e Myr\n", 
            ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->StellarEvolution->get_Z(),
            ptcl->radius*position_unit, ptcl->StellarEvolution->get_supernova()->get_fallback_frac(),
            ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        fprintf(SEVNout, "\tDimless spin. mag: %e, (%e, %e, %e)\n", ptcl->StellarEvolution->getp(Xspin::ID), ptcl->a_spin[0], ptcl->a_spin[1], ptcl->a_spin[2]);
        /* // (SEVN Query) This particle will be kicked after SN feedback in Enzo by EW 2025.4.1
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        */
    }
    else if (ptcl->StellarEvolution->amiempty()) {

        /* // (SEVN Query) dm set to be 0 as SN is processed in Enzo by EW 2025.4.1
        ptcl->dm += ptcl->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        */
        ptcl->SNEjectedMass = ptcl->Mass;
        ptcl->Mass = -1.0;
        fprintf(SEVNout, "Empty. PID: %d, ZAMS Mass: %e Msol, Z: %e, Time: %e Myr, Worldtime: %e Myr\n", 
            ptcl->PID, ptcl->StellarEvolution->get_zams(), ptcl->StellarEvolution->get_Z(), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        if (ptcl->CMPtclIndex != -1) {
            ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        else {
            if (!ptcl->isActive) {
                fprintf(stderr, "Particle %d is not active\n", ptcl->PID);
                fflush(stderr);
                assert(ptcl->isActive);
            }
            // assert(ptcl->isActive); // for debugging by EW 2025.1.20
            ptcl->isActive = false;
            NumberOfParticle--;
        }
        // /* // (SEVN Query) Let's delete later... by EW 2025.4.1
        delete ptcl->StellarEvolution;
        ptcl->StellarEvolution = nullptr;
        // */
    }
}

// Reference: int Mix::special_evolve(Binstar *binstar) in Processes.cpp of SEVN
void Mix(StarSEVN* star1, StarSEVN* star2) {

    // utilities::wait("Hey I have to mix",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);

    StarSEVN *donor = star1; //Donor is the star that will set to empty
    StarSEVN *accretor = star2; //Accretor is the star that will remain as results of the mix

    ///Choose the star that remains and the one that is set to empty
    //First handle the general case: the accretor is the more evolved star,
    // or the more compact (massive) remnant if both are remnant (handled in get_id_more_evolved_star).

    //If one of the star is empty return the other one
    int get_id_more_evolved_star;

    //In the stars have the same phase (not remnant handled before), return the more evolved is the one with the larger plife
    if (donor->getp(Phase::ID)==accretor->getp(Phase::ID))
        get_id_more_evolved_star = donor->plife()>=accretor->plife() ? 0 : 1;
    else
        get_id_more_evolved_star = donor->getp(Phase::ID)>=accretor->getp(Phase::ID) ? 0 : 1;
    if(get_id_more_evolved_star==0)
        //Now the star 0 becomes the accretor and the star 1 the donor
        utilities::swap_stars(donor,accretor);

    //Now handle a special case: If the accretor is a naked helium and the donor not and they have the same innermost core
    //we swap so that the accretor will be the normal star.  This is done because jumping to a normal
    //track from a pureHe is more dangerous since if there are problems on finding a good match the star
    //does not jump and the code raises an error.  In case of a normal track, if we don't find a match
    //we can just continue following the original track.
    //GI 18/02/21: bug fix, to really swap the donor have to have at least a Helium core (the CO core is handled inside).
    //We check >1E-3 instead of >0 to avoid to select the accretor star as the one that is just starting to grow its HE core.
    //GI 02/09/22: We further extend the bug fix to disable the donor swab if the H-star is growing its core from 0, i.e. it is the TMS phase
    //
    //if (accretor->aminakedhelium() and (!donor->aminakedhelium() and donor->getp(MHE::ID)>1E-3)){
    if (accretor->aminakedhelium() and (!donor->aminakedhelium() and
        donor->getp(Phase::ID)>Lookup::TerminalMainSequence and
            donor->getp(Phase::ID)!=Lookup::TerminalCoreHeBurning)){
        unsigned int id_inner_core_accretor = accretor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
        unsigned int id_inner_core_donor =  donor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
        if (id_inner_core_accretor==id_inner_core_donor)
            utilities::swap_stars(donor,accretor);
    }



    ///Let's mix
    //Case 1 - Mix between not remnant stars
    if (!accretor->amiremnant() and !donor->amiremnant()){

        //Set new maximumCO and minmum HE
        accretor->set_MCO_max_aftermerge(donor);
        accretor->set_MHE_min_aftermerge(donor);

        //Mass from the donor
        double DM_new     = donor->getp(Mass::ID);
        double DMHE_new   = donor->getp(MHE::ID);
        double DMCO_new   = donor->getp(MCO::ID);
        double MCO_old    = accretor->getp(MCO::ID);
        double MHE_old    = accretor->getp(MHE::ID);

        //Update masses of the accretor
        accretor->update_from_binary(Mass::ID, DM_new);
        accretor->update_from_binary(dMcumul_binary::ID, DM_new);
        accretor->update_from_binary(MHE::ID, DMHE_new);
        accretor->update_from_binary(MCO::ID, DMCO_new);

        ///Jump to a new track
        //Sanity check on MCO
        if (MCO_old==0 and DMCO_new>0){
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
        }
        // Sanity check on MHE
        else if (MHE_old==0 and DMHE_new>0){
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
        }
        // The accretor is a nakedHelium and the donor not, so we have to trigger a jump to a normal track
        // There is no possibility that the donor has a CO core and the naked helium not because we already check this
        // at the beginning
        else if (accretor->aminakedhelium() and !donor->aminakedhelium()){
            accretor->jump_to_normal_tracks();
        }
        else{
            accretor->find_new_track_after_merger();
        }
    }
    //Case 2 - Mix between WDs
    //TODO Notice that here mixing with a WD is treated as mixing with a NS/BH, in Hurley we have a more complicate outcomes
    else if (accretor->amiWD() and donor->amiWD()) {

        //Check if we have two HeWD
        bool check_double_HeWD = donor->getp(RemnantType::ID)==Lookup::Remnants::HeWD && accretor->getp(RemnantType::ID)==Lookup::Remnants::HeWD;

        //Case 2a - SNI explosion
        if (check_double_HeWD){
            accretor->explode_as_SNI(); // Eunwoo: this should be treated carefully!
            // if (binstar->onesurvived) binstar->set_onesurvived(); //Why this condition? I forgot, we have to check
        }
        else{
            //Update the Mass
            accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
            // binstar->set_onesurvived();
        }

    }
    //Case 3 -  NS, BH or WD mixing with a NS/BH
    else if (accretor->amiCompact() and donor->amiremnant()){
        //Update the Mass
        accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
    }
    //Case 4 - A not remnant donor over a NS, BH, WD
    //-> Do not accreate nothing, just set the donor to empty, see below
    //From commonenvelope::main_stellar_collision in common_envelope.cpp in SEVN1


    ///LOGPRINT AND EVENT SET
    // std::string w =Mix::log_message(binstar,accretor,donor);
    // binstar->print_to_log(w);
    //utilities::hardwait(" Event before ",get_event());
    // set_event((double)Lookup::EventsList::Merger);
    //utilities::hardwait(" Event After ",get_event());

    ///LAST STUFF
    //In any case set the binary to broken and put the donor to empty
    //donor->set_empty_in_bse(binstar);
    donor->set_empty();

}

// Use this function when merger happened
void SetRadius(Particle* ptcl) {

    if (ptcl->ParticleType < REMNANT) {
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*mass_unit/ptcl->StellarEvolution->get_max_zams(), 1./3);
    }
    else if (ptcl->ParticleType <= WHITE_DWARF_ONE) {
        double RNS = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
        double Mch = 1.41/mass_unit;
        double RWD = 0.0115*std::sqrt(pow(Mch/ptcl->Mass,0.6666666667) -  pow(ptcl->Mass/Mch,0.6666666667));
        
        ptcl->radius = std::max(RNS,RWD);
    }
    else if (ptcl->ParticleType <= NEUTRON_STAR_CCSN)
        ptcl->radius = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
    else
        ptcl->radius = 2*ptcl->Mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius in code unit
}

void adjustEnzoSEVNParticle(int index, bool delete_sevn) {

    StarSEVN* star_sevn = SEVNList_Enzo[index];

    PIDtoIndexMap_SEVN.erase(int(star_sevn->get_ID()));

    if (delete_sevn) {
        assert(star_sevn->amiempty());
        delete star_sevn;
    }
    star_sevn = nullptr;
        
    if (index != SEVNList_Enzo.size() - 1) {
        std::swap(SEVNList_Enzo[index], SEVNList_Enzo.back());
        std::swap(creation_time_Enzo[index], creation_time_Enzo.back());
        std::swap(world_time_Enzo[index], world_time_Enzo.back());
        PIDtoIndexMap_SEVN[int(SEVNList_Enzo[index]->get_ID())] = index;
    }
    SEVNList_Enzo.pop_back();
    creation_time_Enzo.pop_back();
    world_time_Enzo.pop_back();
}

#endif