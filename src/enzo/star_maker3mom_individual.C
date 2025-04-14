#ifdef NBODY

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm> 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "phys_constants.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"


int GetUnits(float *DensityUnits, float *LengthUnits,
    float *TemperatureUnits, float *TimeUnits,
    float *VelocityUnits, FLOAT Time);

void grid::individual_star_feedback3mom(const float &dx, const float &kinf_in, float *mu, const float &yield) {

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
		Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}


    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum, MetalIaNum, MetalIINum;
    if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
        ENZO_FAIL("Error in grid->IdentifyColourFields.\n");


    // dx = CellWidthTemp == float(CellWidth[0][0])

    /* Get Units */
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits");
    }
    MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
    EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;


    int nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
    int ibuff = NumberOfGhostZones;
    float xstart = this->CellLeftEdge[0][0];
    float ystart = this->CellLeftEdge[1][0];
    float zstart = this->CellLeftEdge[2][0];

    float xfc, yfc, zfc, xfcshift, yfcshift, zfcshift;
    float fbuff;
    float face_shift;
    int iface, jface, kface;
    float dxf, dyf, dzf;
    float xface, yface, zface, xpos, ypos, zpos;
    int ic, jc, kc;
    float dxc, dyc, dzc;

    // Assuming 3x3x3 cube, calculate cell distribution
    int distrad = 3; // feedback distribution radius in cells
    float dist_mass_cells = distrad*distrad*distrad;

    float mass_per_cell;
    float energy, energy_per_cell;
    
    float kinf;
    float energy51;
    float Zsol, num_d, d_ave;
    float mu_cell;
    float t_PDS, R_PDS;
    float R_resolve;
    float thermal_energy_per_cell;
    float mom_per_cell;
    float ke_after, delta_ke;
    float asum, bsum, csum;
    int idir;
    float m_eject;
    float energy_before, energy_after, mass_before, mass_after, kin_energy_before, kin_energy_after;
    float ke_injected;

    float mass_ejected_Msun;
    float T_eff;
    float v_wind_low_cgs = 20.0 * km_cm;
    float v_wind_high_cgs = 100.0 * km_cm;

    // Allocate a 3D array using pointers to pointers
    float*** u1 = new float**[4];
    float*** v1 = new float**[4];
    float*** w1 = new float**[4];
    float*** d1 = new float**[4];
    float*** ge1 = new float**[4];
    float*** te1 = new float**[4];
    float*** metal1 = new float**[4];
    float*** ke_before = new float**[4];

    // Allocate the second dimension (NY)
    for (int i = 0; i < 4; ++i) {
        u1[i] = new float*[4];
        v1[i] = new float*[4];
        w1[i] = new float*[4];
        d1[i] = new float*[4];
        ge1[i] = new float*[4];
        te1[i] = new float*[4];
        metal1[i] = new float*[4];
        ke_before[i] = new float*[4];
        
        // Allocate the third dimension (NZ)
        for (int j = 0; j < 4; ++j) {
            u1[i][j] = new float[4];
            v1[i][j] = new float[4];
            w1[i][j] = new float[4];
            d1[i][j] = new float[4];
            ge1[i][j] = new float[4];
            te1[i][j] = new float[4];
            metal1[i][j] = new float[4];
            ke_before[i][j] = new float[4];
        }
    }

    for (int n=0; n < this->NumberOfParticles; n++) {
        // Feedback condition: creation time >= 0 && (WindEjectedMass > 0.0 || SNejectedMass > 0.0)
        // (SEVN Query) This conditions might need to be fixed later!
        if (this->ParticleAttribute[0][n] >= 0 && (
                    ParticleAttribute[NumberOfParticleAttributes-8+1][n] > 0.0 || ParticleAttribute[NumberOfParticleAttributes-8+2][n] > 0.0)) {

            if (this->ParticlePosition[0][n] < xstart || this->ParticlePosition[0][n] > xstart + dx * nx ||
                this->ParticlePosition[1][n] < ystart || this->ParticlePosition[1][n] > ystart + dx * ny ||
                this->ParticlePosition[2][n] < zstart || this->ParticlePosition[2][n] > zstart + dx * nz) {

                ENZO_FAIL("star particle out of grid");
            }

            // Set center of feedback zone
            xfc = this->ParticlePosition[0][n];
            yfc = this->ParticlePosition[1][n];
            zfc = this->ParticlePosition[2][n];
            fbuff = ibuff + 2.0;

            // Check bounds - if star particle is near grid edge then shift center of feedback region
            if (xfc < xstart + fbuff*dx || xfc > xstart + dx * (nx - fbuff) ||
                yfc < ystart + fbuff*dx || yfc > ystart + dx * (ny - fbuff) ||
                zfc < zstart + fbuff*dx || zfc > zstart + dx * (nz - fbuff)) {

                xfcshift = xfc;
                yfcshift = yfc;
                zfcshift = zfc;

                xfc = max(xfc, xstart+fbuff*dx);
                yfc = max(yfc, ystart+fbuff*dx);
                zfc = max(zfc, zstart+fbuff*dx);

                xfc = min(xfc, xstart + dx * (nx - fbuff - 1));
                yfc = min(yfc, ystart + dx * (ny - fbuff - 1));
                zfc = min(zfc, zstart + dx * (nz - fbuff - 1));

                xfcshift = xfcshift - xfc;
                yfcshift = yfcshift - yfc;
                zfcshift = zfcshift - zfc;
            }

            // If using zeus, then velocities are face-centered so shift
            face_shift = 0.0;
            if (HydroMethod == 2) face_shift = 0.5;

            // Compute index of the first cell to add momentum. accounting for possible face-centering
            xface = (xfc - xstart)/dx - face_shift;
            yface = (yfc - ystart)/dx - face_shift;
            zface = (zfc - zstart)/dx - face_shift;

            iface = int(xface);
            jface = int(yface);
            kface = int(zface);

            dxf = iface + 1.0 - xface;
            dyf = jface + 1.0 - yface;
            dzf = kface + 1.0 - zface;

            // Compute index of the first cell to add mass, assuming cell-centering
            xpos = (xfc - xstart)/dx;
            ypos = (yfc - ystart)/dx;
            zpos = (zfc - zstart)/dx;

            ic = int(xpos);
            jc = int(ypos);
            kc = int(zpos);

            dxc = ic + 1.0 - xpos;
            dyc = jc + 1.0 - ypos;
            dzc = kc + 1.0 - zpos;

            // Stellar wind feedback
            if (this->ParticleAttribute[NumberOfParticleAttributes-8+1][n] > 0.0) {

                mass_ejected_Msun = this->ParticleAttribute[NumberOfParticleAttributes-8+1][n];
                T_eff = this->ParticleAttribute[NumberOfParticleAttributes-8+3][n];
                assert(T_eff > 0.0);

                // Ejected mass fraction
                m_eject = (mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx)) / (this->ParticleMass[n] + mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx));
                assert(m_eject > 0.0 && m_eject <= 1.0);

                // Stellar mass for AGB stars (M <= 8 Msun, v_wind = 20 km/s)
                if (this->ParticleAttribute[NumberOfParticleAttributes-8+0][n] <= 8.0) {

                    fprintf(stderr, "PID: %d. Weak wind!!! WindEjectedMass: %e Msun\n", this->ParticleNumber[n], mass_ejected_Msun);

                    // Calculate mass per cell (total 27) ejected
                    mass_per_cell = (mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx)) / dist_mass_cells; // (SEVN Query) This should be code unit!!!

                    // 1.5 * k_B * N * T +  kinetic energy of wind ==> Let's assume fully thermalized wind as AEOS did!
                    energy = 1.5 * T_eff * (mass_ejected_Msun * SolarMass / (mh)) * kboltz; // current T of wind
                    energy += 0.5 * (mass_ejected_Msun * SolarMass) * v_wind_low_cgs * v_wind_low_cgs; // assume 100% KE thermalization
                    energy = energy / EnergyUnits / (dx*dx*dx); // (SEVN Query) This should be code unit!!!
                    energy_per_cell = energy / dist_mass_cells;
                    
                    // Fully thermalized wind energy
                    kinf = 0.0;

                    // Finished computing kinf; now compute the amount of thermal energy that needs to be added to each cell
                    thermal_energy_per_cell = (1.0 - kinf) * energy_per_cell;

                    // Inject energy, mass & metals
                    // Zero local dummy field and kinetic energy field
                    for (int k = 0; k < 4; k++) {
                        for (int j = 0; j < 4; j++) {
                            for (int i = 0; i < 4; i++) {
                                u1[i][j][k] = 0.0f;
                                v1[i][j][k] = 0.0f;
                                w1[i][j][k] = 0.0f;
                                d1[i][j][k] = 0.0f;
                                ge1[i][j][k] = 0.0f;
                                te1[i][j][k] = 0.0f;
                                metal1[i][j][k] = 0.0f;
                            }
                        }
                    }

                    // Compute the kinetic energy in the affected region before momentum is added (except for ZEUS).
                    // This is needed at the end of the calculation to update the total energy (te) field.
                    if (HydroMethod != 2) {
                        for (int k = -1; k <= 2; k++) {
                            for (int j = -1; j <= 2; j++) {
                                for (int i = -1; i <= 2; i++) {
                                    int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                    ke_before[i+1][j+1][k+1] = 0.5 * this->BaryonField[DensNum][index] * (
                                        this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                        this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                        this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                    );
                                }
                            }
                        }
                    }

                    // First convert velocities to momenta and transform into frame comoving with particle
                    idir = +1;
                    this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                    // Sum mass and energy before
                    this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_before, kin_energy_before);

                    this->add_feedback1(u1, v1, w1, d1, ge1, te1, metal1, 
                        dxf, dyf, dzf, dxc, dyc, dzc, m_eject, yield, this->ParticleAttribute[2][n],
                        mass_per_cell, 1.0, 0.0);

                    // Fully thermalized wind energy
                    mom_per_cell = 0.0;

                    // Now add mass and momentum, using three-point CIC
                    // Note that if the te field is being used, it is only updated with the thermal energy
                    this->add_feedback2(nx, ny, nz,
                        ic, jc, kc, iface, jface, kface,
                        dxf, dyf, dzf, dxc, dyc, dzc,
                        m_eject, yield, this->ParticleAttribute[2][n],
                        mass_per_cell, mom_per_cell, thermal_energy_per_cell);

                    // Sum mass and energy after
                    this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_after, kin_energy_after);

                    if (kinf != 0.0 && abs(kin_energy_after - kin_energy_before - kinf*energy)/(kinf*energy) > 0.01) {
                        fprintf(stderr, "Kinetic energy added to mesh does not match!!!\n");
                        fprintf(stderr, "kinf: %e\n", kinf);
                        fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                        fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                        fprintf(stderr, "diff: %e %\n", abs(kin_energy_after - kin_energy_before)/kin_energy_before);
                    }
                    if (abs(mass_after - mass_before - mass_per_cell*dist_mass_cells)/(mass_per_cell*dist_mass_cells) > 0.01) {
                        fprintf(stderr, "Mass added to mesh does not match!!!\n");
                        fprintf(stderr, "kinf: %e\n", kinf);
                        fprintf(stderr, "mass_after: %e, mass_before: %e\n", mass_after, mass_before);
                        fprintf(stderr, "mass_per_cell*dist_mass_cells: %e\n", mass_per_cell*dist_mass_cells);
                        fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                        fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                        fprintf(stderr, "diff: %e %\n", abs(mass_after - mass_before)/mass_before);
                    }

                    // Convert momenta back to velocities and transform back to lab frame
                    idir = -1;
                    this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                    // Add the increase in the kinetic energy to the total energy field (unless we're using Zeus).
                    // If using dual energy formalism, we might want to enforce consistency.
                    if (HydroMethod != 2) {
                        for (int k = -1; k <= 2; k++) {
                            for (int j = -1; j <= 2; j++) {
                                for (int i = -1; i <= 2; i++) {
                                    int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                    
                                    ke_after = 0.5 * this->BaryonField[DensNum][index] * (
                                        this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                        this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                        this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                    );
                                    delta_ke = ke_after - ke_before[i+1][j+1][k+1];
                                    this->BaryonField[TENum][index] += delta_ke/this->BaryonField[DensNum][index];
                                    ke_injected += delta_ke;
                                }
                            }
                        }
                    }
                } // Low-mass star wind feedback
                // Stellar mass for massive stars (M > 8 Msun, v_wind = 100 km/s)
                else {

                    fprintf(stderr, "PID: %d. Strong wind!!! WindEjectedMass: %e Msun\n", this->ParticleNumber[n], mass_ejected_Msun);

                    // Calculate mass per cell (total 27) ejected
                    mass_per_cell = (mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx)) / dist_mass_cells; // (SEVN Query) This should be code unit!!!

                    // 1.5 * k_B * N * T +  kinetic energy of wind ==> Let's assume fully thermalized wind as AEOS did!
                    energy = 1.5 * T_eff * (mass_ejected_Msun * SolarMass / (mh)) * kboltz; // current T of wind
                    energy += 0.5 * (mass_ejected_Msun * SolarMass) * v_wind_high_cgs * v_wind_high_cgs; // assume 100% KE thermalization
                    energy = energy / EnergyUnits / (dx*dx*dx); // (SEVN Query) This should be code unit!!!
                    energy_per_cell = energy / dist_mass_cells;
                    
                    // Fully thermalized wind energy
                    kinf = 0.0;

                    // Finished computing kinf; now compute the amount of thermal energy that needs to be added to each cell
                    thermal_energy_per_cell = (1.0 - kinf) * energy_per_cell;

                    // Inject energy, mass & metals
                    // Zero local dummy field and kinetic energy field
                    for (int k = 0; k < 4; k++) {
                        for (int j = 0; j < 4; j++) {
                            for (int i = 0; i < 4; i++) {
                                u1[i][j][k] = 0.0f;
                                v1[i][j][k] = 0.0f;
                                w1[i][j][k] = 0.0f;
                                d1[i][j][k] = 0.0f;
                                ge1[i][j][k] = 0.0f;
                                te1[i][j][k] = 0.0f;
                                metal1[i][j][k] = 0.0f;
                            }
                        }
                    }

                    // Compute the kinetic energy in the affected region before momentum is added (except for ZEUS).
                    // This is needed at the end of the calculation to update the total energy (te) field.
                    if (HydroMethod != 2) {
                        for (int k = -1; k <= 2; k++) {
                            for (int j = -1; j <= 2; j++) {
                                for (int i = -1; i <= 2; i++) {
                                    int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                    ke_before[i+1][j+1][k+1] = 0.5 * this->BaryonField[DensNum][index] * (
                                        this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                        this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                        this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                    );
                                }
                            }
                        }
                    }

                    // First convert velocities to momenta and transform into frame comoving with particle
                    idir = +1;
                    this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                    // Sum mass and energy before
                    this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_before, kin_energy_before);

                    this->add_feedback1(u1, v1, w1, d1, ge1, te1, metal1, 
                        dxf, dyf, dzf, dxc, dyc, dzc, m_eject, yield, this->ParticleAttribute[2][n],
                        mass_per_cell, 1.0, 0.0);

                    // Fully thermalized wind energy
                    mom_per_cell = 0.0;

                    // Now add mass and momentum, using three-point CIC
                    // Note that if the te field is being used, it is only updated with the thermal energy
                    this->add_feedback2(nx, ny, nz,
                        ic, jc, kc, iface, jface, kface,
                        dxf, dyf, dzf, dxc, dyc, dzc,
                        m_eject, yield, this->ParticleAttribute[2][n],
                        mass_per_cell, mom_per_cell, thermal_energy_per_cell);

                    // Sum mass and energy after
                    this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_after, kin_energy_after);

                    if (kinf != 0.0 && abs(kin_energy_after - kin_energy_before - kinf*energy)/(kinf*energy) > 0.01) {
                        fprintf(stderr, "Kinetic energy added to mesh does not match!!!\n");
                        fprintf(stderr, "kinf: %e\n", kinf);
                        fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                        fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                        fprintf(stderr, "diff: %e %\n", abs(kin_energy_after - kin_energy_before)/kin_energy_before);
                    }
                    if (abs(mass_after - mass_before - mass_per_cell*dist_mass_cells)/(mass_per_cell*dist_mass_cells) > 0.01) {
                        fprintf(stderr, "Mass added to mesh does not match!!!\n");
                        fprintf(stderr, "kinf: %e\n", kinf);
                        fprintf(stderr, "mass_after: %e, mass_before: %e\n", mass_after, mass_before);
                        fprintf(stderr, "mass_per_cell*dist_mass_cells: %e\n", mass_per_cell*dist_mass_cells);
                        fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                        fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                        fprintf(stderr, "diff: %e %\n", abs(mass_after - mass_before)/mass_before);
                    }

                    // Convert momenta back to velocities and transform back to lab frame
                    idir = -1;
                    this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                    // Add the increase in the kinetic energy to the total energy field (unless we're using Zeus).
                    // If using dual energy formalism, we might want to enforce consistency.
                    if (HydroMethod != 2) {
                        for (int k = -1; k <= 2; k++) {
                            for (int j = -1; j <= 2; j++) {
                                for (int i = -1; i <= 2; i++) {
                                    int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                    
                                    ke_after = 0.5 * this->BaryonField[DensNum][index] * (
                                        this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                        this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                        this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                    );
                                    delta_ke = ke_after - ke_before[i+1][j+1][k+1];
                                    this->BaryonField[TENum][index] += delta_ke/this->BaryonField[DensNum][index];
                                    ke_injected += delta_ke;
                                }
                            }
                        }
                    }
                } // High-mass star wind feedback
                this->ParticleAttribute[NumberOfParticleAttributes-8+1][n] = 0.0; // Set wind mass to 0.0
            }

            // Supernova feedback
            // /*
            if (this->ParticleAttribute[NumberOfParticleAttributes-8+2][n] > 0.0) {

                mass_ejected_Msun = this->ParticleAttribute[NumberOfParticleAttributes-8+2][n];

                // Ejected mass fraction
                m_eject = (mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx)) / (this->ParticleMass[n] + mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx));
                assert(m_eject > 0.0 && m_eject <= 1.0);

                fprintf(stderr, "PID: %d. Supernova!!! SNEjectedMass: %e Msun\n", this->ParticleNumber[n], mass_ejected_Msun);

                // Calculate mass per cell (total 27) ejected
                mass_per_cell = (mass_ejected_Msun * SolarMass / MassUnits / (dx*dx*dx)) / dist_mass_cells; // (SEVN Query) This should be code unit!!!

                // Calculate how much of the star formation in this timestep would have gone into supernova energy
                // (SEVN Query) Let's use 10^51 erg here!
                // energy = sn_param * mform * (clight/vunits) * (clight/vunits);
                // energy_per_cell = energy / dist_mass_cells;

                energy = 1e51 / EnergyUnits / (dx*dx*dx); // (SEVN Query) This should be code unit!!!
                energy_per_cell = energy / dist_mass_cells;

                // Use fixed kinf value, unless kinf < 0 - then compute variable kinf
                kinf = kinf_in;
                if (kinf < 0) {

                    energy51 = 1.0; // (SEVN Query) We are setting supernova energy into fixed value of 10**51 erg

                    Zsol = 0.0;
                    num_d = 0.0;
                    d_ave = 0.0;

                    // Compute averave density and metallicity around cell containing particle
                    for (int k = -1; k <= 1; k++) {
                        for (int j = -1; j <= 1; j++) {
                            for (int i = -1; i <= 1; i++) {
                                // Compute the 1D flattened index
                                int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;

                                // Access mu_field just like mu(ic+i, jc+j, kc+k) in Fortran
                                mu_cell = mu[index];
                                Zsol += this->BaryonField[MetalNum][index] / 0.02;
                                num_d += this->BaryonField[DensNum][index] * DensityUnits / (mu_cell * mh);
                                d_ave += this->BaryonField[DensNum][index] * DensityUnits;
                            }
                        }
                    }

                    Zsol = Zsol / 27.0;
                    num_d = num_d / 27.0;
                    d_ave = d_ave / 27.0;

                    // Compute time and radius of transition to PDS phase for gas with the computed properties
                    // Note: t_PDS has units of 1e3 yrs and R_PDS has units of pc

                    // For metal poor gas
                    if (Zsol < 0.01) {
                        t_PDS = 3.06e2 * pow(energy51, 1.0/8.0) * pow(num_d, -3.0/4.0);
                        R_PDS = 49.3 * pow(energy51, 0.25) * pow(num_d, -0.5);
                    } 
                    // For metal rich gas
                    else {
                        t_PDS = 26.5 * pow(energy51, 3.0/14.0) * pow(Zsol, -5.0/14.0) * pow(num_d, -4.0/7.0);
                        R_PDS = 18.5 * pow(energy51, 2.0/7.0) * pow(num_d, -3.0/7.0) * pow(Zsol, -1.0/7.0);
                    }
                    R_resolve = dx*LengthUnits/pc_cm; // now in pc scale as R_PDS
                    if (R_PDS > 4.5*R_resolve) {
                        kinf = 0.0;
                    } else {
                        kinf = 3.97133e-6 * (d_ave/mh) * pow(R_resolve, -2.0) * pow(R_PDS, 7.0) * pow(t_PDS, -2.0) * pow(energy51, -1.0);
                    }
                }
                
                // Check kinf. If it is too small, set it to 0
                if (kinf <= 1e-10) kinf = 0.0;
                // kinf = 0.0;

                // Finished computing kinf; now compute the amount of thermal energy that needs to be added to each cell
                thermal_energy_per_cell = (1.0 - kinf) * energy_per_cell;

                // Inject energy, mass & metals
                // Zero local dummy field and kinetic energy field
                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        for (int i = 0; i < 4; i++) {
                            u1[i][j][k] = 0.0f;
                            v1[i][j][k] = 0.0f;
                            w1[i][j][k] = 0.0f;
                            d1[i][j][k] = 0.0f;
                            ge1[i][j][k] = 0.0f;
                            te1[i][j][k] = 0.0f;
                            metal1[i][j][k] = 0.0f;
                        }
                    }
                }

                // Compute the kinetic energy in the affected region before momentum is added (except for ZEUS).
                // This is needed at the end of the calculation to update the total energy (te) field.
                if (HydroMethod != 2) {
                    for (int k = -1; k <= 2; k++) {
                        for (int j = -1; j <= 2; j++) {
                            for (int i = -1; i <= 2; i++) {
                                int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                ke_before[i+1][j+1][k+1] = 0.5 * this->BaryonField[DensNum][index] * (
                                    this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                    this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                    this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                );
                            }
                        }
                    }
                }

                // First convert velocities to momenta and transform into frame comoving with particle
                idir = +1;
                this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                // Sum mass and energy before
                this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_before, kin_energy_before);

                this->add_feedback1(u1, v1, w1, d1, ge1, te1, metal1, 
                    dxf, dyf, dzf, dxc, dyc, dzc, m_eject, yield, this->ParticleAttribute[2][n],
                    mass_per_cell, 1.0, 0.0);

                // The kinetic energy after the feedback event is a quadratic equation with delta_p,
                // where delta_p is the total amount of momentum added to the feedback region:
                // E_k, a = asum + bsum*delta_p + csum*delta_p*delta_p;

                // Sum a, b, and c terms to get momentum normalization
                if (kinf > 0.0) {
                    this->sum_abc(u1, v1, d1, iface, jface, kface,
                        ic, jc, kc, asum, bsum, csum);

                    asum = asum - (kin_energy_before + kinf*energy);
                    
                    // Calculate momentum contribution
                    mom_per_cell = (-bsum + sqrt(bsum*bsum - 4.0*asum*csum)) / (2.0*csum);
                } else {
                    // If the kinetic fraction is set to 0, then 0 momentum is added to each cell.
                    // Note that this is different from adding 0 kinetic energy to
                    mom_per_cell = 0.0;
                }

                // Now add mass and momentum, using three-point CIC
                // Note that if the te field is being used, it is only updated with the thermal energy
                this->add_feedback2(nx, ny, nz,
                    ic, jc, kc, iface, jface, kface,
                    dxf, dyf, dzf, dxc, dyc, dzc,
                    m_eject, yield, this->ParticleAttribute[2][n],
                    mass_per_cell, mom_per_cell, thermal_energy_per_cell);

                // Sum mass and energy after
                this->sum_mass_kinetic_energy(n, iface, jface, kface, ic, jc, kc, mass_after, kin_energy_after);

                if (kinf != 0.0 && abs(kin_energy_after - kin_energy_before - kinf*energy)/(kinf*energy) > 0.01) {
                    fprintf(stderr, "Kinetic energy added to mesh does not match!!!\n");
                    fprintf(stderr, "kinf: %e\n", kinf);
                    fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                    fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                    fprintf(stderr, "diff: %e %\n", abs(kin_energy_after - kin_energy_before)/kin_energy_before);
                }
                if (abs(mass_after - mass_before - mass_per_cell*dist_mass_cells)/(mass_per_cell*dist_mass_cells) > 0.01) {
                    fprintf(stderr, "Mass added to mesh does not match!!!\n");
                    fprintf(stderr, "kinf: %e\n", kinf);
                    fprintf(stderr, "mass_after: %e, mass_before: %e\n", mass_after, mass_before);
                    fprintf(stderr, "mass_per_cell*dist_mass_cells: %e\n", mass_per_cell*dist_mass_cells);
                    fprintf(stderr, "kin_energy_after: %e, kin_energy_before: %e\n", kin_energy_after, kin_energy_before);
                    fprintf(stderr, "energy*kinf: %e\n", kinf*energy);
                    fprintf(stderr, "diff: %e %\n", abs(mass_after - mass_before)/mass_before);
                }

                // Convert momenta back to velocities and transform back to lab frame
                idir = -1;
                this->momentum(n, ic, jc, kc, iface, jface, kface, idir);

                // Add the increase in the kinetic energy to the total energy field (unless we're using Zeus).
                // If using dual energy formalism, we might want to enforce consistency.
                if (HydroMethod != 2) {
                    for (int k = -1; k <= 2; k++) {
                        for (int j = -1; j <= 2; j++) {
                            for (int i = -1; i <= 2; i++) {
                                int index = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                                
                                ke_after = 0.5 * this->BaryonField[DensNum][index] * (
                                    this->BaryonField[Vel1Num][index] * this->BaryonField[Vel1Num][index] +
                                    this->BaryonField[Vel2Num][index] * this->BaryonField[Vel2Num][index] +
                                    this->BaryonField[Vel3Num][index] * this->BaryonField[Vel3Num][index]
                                );
                                delta_ke = ke_after - ke_before[i+1][j+1][k+1];
                                this->BaryonField[TENum][index] += delta_ke/this->BaryonField[DensNum][index];
                                ke_injected += delta_ke;
                            }
                        }
                    }
                }
                this->ParticleAttribute[NumberOfParticleAttributes-8+2][n] = 0.0; // Set SN mass to 0.0 
            } // SN feedback
            // */
            // this->ParticleAttribute[NumberOfParticleAttributes-8+2][n] = 0.0; // Set SN mass to 0.0 
        } // for proper particles only
    } // every particle loop

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            delete[] u1[i][j];
            delete[] v1[i][j];
            delete[] w1[i][j];
            delete[] d1[i][j];
            delete[] ge1[i][j];
            delete[] te1[i][j];
            delete[] metal1[i][j];
            delete[] ke_before[i][j];
        }
        delete[] u1[i];
        delete[] v1[i];
        delete[] w1[i];
        delete[] d1[i];
        delete[] ge1[i];
        delete[] te1[i];
        delete[] metal1[i];
        delete[] ke_before[i];
    }
    
    // Deallocate the first dimension
    delete[] u1;
    delete[] v1;
    delete[] w1;
    delete[] d1;
    delete[] ge1;
    delete[] te1;
    delete[] metal1;
    delete[] ke_before;
}

void grid::momentum(const int &ParticleIndex,
    const int &ic, const int &jc, const int &kc,
    const int &iface, const int &jface, const int &kface,
    const int &idir) {

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
		Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}

    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum, MetalIaNum, MetalIINum;
    if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
        ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

    int MetallicityField = (MetalNum != -1 || SNColourNum != -1);

    int nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];

    int index1, index2;

    if (idir != -1 && idir != 1) {
        ENZO_FAIL("Invalid direction for momentum calculation");
    }

    for (int k = -1; k <= 2; k++) {
        for (int j = -1; j <= 2; j++) {
            for (int i = -1; i <= 2; i++) {

                // idir = +1: convert vel -> mom
                if (idir == 1) {
                    if (HydroMethod == 2) {

                        index1 = (iface + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        index2 = (iface + i + 1) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[Vel1Num][index1] = (this->BaryonField[Vel1Num][index1] - this->ParticleVelocity[0][ParticleIndex]) * 0.5 * 
                                                                        (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2]);
                        
                        index1 = (ic + i) + (jface + j) * nx + (kc + k) * nx * ny;
                        index2 = (ic + i) + (jface + j + 1) * nx + (kc + k) * nx * ny;                                  
                        this->BaryonField[Vel2Num][index1] = (this->BaryonField[Vel2Num][index1] - this->ParticleVelocity[1][ParticleIndex]) * 0.5 * 
                                                                        (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2]);
                        
                        index1 = (ic + i) + (jc + j) * nx + (kface + k) * nx * ny;
                        index2 = (ic + i) + (jc + j) * nx + (kface + k + 1) * nx * ny;
                        this->BaryonField[Vel3Num][index1] = (this->BaryonField[Vel3Num][index1] - this->ParticleVelocity[2][ParticleIndex]) * 0.5 * 
                                                                        (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2]);
                    } else {
                        index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[Vel1Num][index1] = (this->BaryonField[Vel1Num][index1] - this->ParticleVelocity[0][ParticleIndex]) * 
                                                                        this->BaryonField[DensNum][index1];
                        this->BaryonField[Vel2Num][index1] = (this->BaryonField[Vel2Num][index1] - this->ParticleVelocity[1][ParticleIndex]) * 
                                                                        this->BaryonField[DensNum][index1];
                        this->BaryonField[Vel3Num][index1] = (this->BaryonField[Vel3Num][index1] - this->ParticleVelocity[2][ParticleIndex]) *
                                                                        this->BaryonField[DensNum][index1];
                    }
                    if (MetallicityField == 1) {
                        index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[MetalNum][index1] = this->BaryonField[MetalNum][index1] * 
                                                                        this->BaryonField[DensNum][index1];
                    }
                }
                // idir = -1: convert mom - >vel
                else {
                    if (HydroMethod == 2) {
                        index1 = (iface + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        index2 = (iface + i + 1) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[Vel1Num][index1] = this->BaryonField[Vel1Num][index1] / 
                                                                        (0.5 * (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2])) + 
                                                                        this->ParticleVelocity[0][ParticleIndex];
                        
                        index1 = (ic + i) + (jface + j) * nx + (kc + k) * nx * ny;
                        index2 = (ic + i) + (jface + j + 1) * nx + (kc + k) * nx * ny;                                  
                        this->BaryonField[Vel2Num][index1] = this->BaryonField[Vel2Num][index1] / 
                                                                        (0.5 * (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2])) + 
                                                                        this->ParticleVelocity[1][ParticleIndex];

                        index1 = (ic + i) + (jc + j) * nx + (kface + k) * nx * ny;
                        index2 = (ic + i) + (jc + j) * nx + (kface + k + 1) * nx * ny;
                        this->BaryonField[Vel3Num][index1] = this->BaryonField[Vel3Num][index1] / 
                                                                        (0.5 * (this->BaryonField[DensNum][index1] + this->BaryonField[DensNum][index2])) + 
                                                                        this->ParticleVelocity[2][ParticleIndex];
                    } else {
                        index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[Vel1Num][index1] = this->BaryonField[Vel1Num][index1] / 
                                                                        this->BaryonField[DensNum][index1] + this->ParticleVelocity[0][ParticleIndex];
                        this->BaryonField[Vel2Num][index1] = this->BaryonField[Vel2Num][index1] / 
                                                                        this->BaryonField[DensNum][index1] + this->ParticleVelocity[1][ParticleIndex];
                        this->BaryonField[Vel3Num][index1] = this->BaryonField[Vel3Num][index1] /
                                                                        this->BaryonField[DensNum][index1] + this->ParticleVelocity[2][ParticleIndex];
                    }
                    if (MetallicityField == 1) {
                        index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                        this->BaryonField[MetalNum][index1] = this->BaryonField[MetalNum][index1] / 
                                                                        this->BaryonField[DensNum][index1];

                    }
                }
            }
        }
    }
}

void grid::sum_mass_kinetic_energy(const int &ParticleIndex,
    const int &iface, const int &jface, const int &kface,
    const int &ic, const int &jc, const int &kc,
    float &mass_sum, float &kin_energy_sum) {

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
		Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}

    int nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];

    float mass_term, mom_term, kin_energy;
    int index1, index2, index3, index4;

    mass_sum = 0.0;
    kin_energy_sum = 0.0;

    for (int k = -1; k <= 2; k++) {
        for (int j = -1; j <= 2; j++) {
            for (int i = -1; i <= 2; i++) {

                index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                index2 = (iface + i) + (jc + j) * nx + (kc + k) * nx * ny;
                index3 = (ic + i) + (jface + j) * nx + (kc + k) * nx * ny;
                index4 = (ic + i) + (jc + j) * nx + (kface + k) * nx * ny;

                mass_term = this->BaryonField[DensNum][index1];
                mom_term = this->BaryonField[Vel1Num][index2] * this->BaryonField[Vel1Num][index2] +
                           this->BaryonField[Vel2Num][index3] * this->BaryonField[Vel2Num][index3] +
                           this->BaryonField[Vel3Num][index4] * this->BaryonField[Vel3Num][index4];

                kin_energy = mom_term / (2.0 * mass_term);
                mass_sum += mass_term;
                kin_energy_sum += kin_energy;
            }
        }
    }
}

void grid::sum_abc(float ***u1, float ***v1, float *** d1,
    const int &iface, const int &jface, const int &kface,
    const int &ic, const int &jc, const int &kc,
    float &asum, float &bsum, float &csum) {

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
		Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}

    int nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
    
    asum = 0.0, bsum = 0.0, csum = 0.0;
    float mass_sum = 0.0, energy_sum = 0.0;
    int istart = 1, jstart = 1, kstart = 1;

    float mass_term, mom_term, aterm, bterm, cterm;
    int index1, index2, index3, index4;

    for (int k = -1; k <= 2; k++) {
        for (int j = -1; j <= 2; j++) {
            for (int i = -1; i <= 2; i++) {
                
                index1 = (ic + i) + (jc + j) * nx + (kc + k) * nx * ny;
                index2 = (iface + i) + (jc + j) * nx + (kc + k) * nx * ny;
                index3 = (ic + i) + (jface + j) * nx + (kc + k) * nx * ny;
                index4 = (ic + i) + (jc + j) * nx + (kface + k) * nx * ny;

                mass_term = this->BaryonField[DensNum][index1];
                mom_term = this->BaryonField[Vel1Num][index2] * this->BaryonField[Vel1Num][index2] +
                           this->BaryonField[Vel2Num][index3] * this->BaryonField[Vel2Num][index3] +
                           this->BaryonField[Vel3Num][index4] * this->BaryonField[Vel3Num][index4];

                mass_term += d1[istart+i][jstart+j][kstart+k];

                asum += mom_term/(2.0 * mass_term);

                bterm = this->ParticleVelocity[0][index1]*u1[istart+i][jstart+j][kstart+k] + 
                        this->ParticleVelocity[1][index2]*v1[istart+i][jstart+j][kstart+k] +
                        this->ParticleVelocity[2][index3]*d1[istart+i][jstart+j][kstart+k];
                bsum += bterm/mass_term;
                
                cterm = u1[istart+i][jstart+j][kstart+k] * u1[istart+i][jstart+j][kstart+k] +
                        v1[istart+i][jstart+j][kstart+k] * v1[istart+i][jstart+j][kstart+k] +
                        d1[istart+i][jstart+j][kstart+k] * d1[istart+i][jstart+j][kstart+k];
                csum += cterm/(2.0 * mass_term);
            }
        }
    }
}

void grid::add_feedback1(float ***u1, float ***v1, float ***w1, float ***d1, float ***ge1, float ***te1, float ***metal1,
    const float &dxf, const float &dyf, const float &dzf, 
    const float &dxc, const float &dyc, const float &dzc,
    const float &m_eject, const float &yield, const float &metalf,
    const float &mass_per_cell, const float mom_per_cell, const float therm_per_cell) {

    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum, MetalIaNum, MetalIINum;
    if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
        ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

    int MetallicityField = (MetalNum != -1 || SNColourNum != -1);

    const int ic = 1, jc = 1, kc = 1, iface = 1, jface = 1, kface = 1;

    float delta_mass, delta_pu, delta_pv, delta_pw, delta_therm, dratio;
    float dxf1, dyf1, dzf1, dxc1, dyc1, dzc1;

    if (MultiMetals == 1) {
        ENZO_FAIL("momentum: not supported");
    }

    float total_mass = 0.0;

    for (int k = -1; k <= 1; k++) {
        for (int j = -1; j <= 1; j++) {
            for (int i = -1; i <= 1; i++) {
                for (int i1 = i; i1 <= i + 1; i1++) {
                    dxf1 = dxf;
                    dxc1 = dxc;
                    if (i1 == i + 1) {
                        dxf1 = 1.0 - dxf;
                        dxc1 = 1.0 - dxc;
                    }
                    for (int j1 = j; j1 <= j + 1; j1++) {
                        dyf1 = dyf;
                        dyc1 = dyc;
                        if (j1 == j + 1) {
                            dyf1 = 1.0 - dyf;
                            dyc1 - 1.0 - dyc;
                        }
                        for (int k1 = k; k1 <= k + 1; k1++) {
                            dzf1 = dzf;
                            dzc1 = dzc;
                            if (k1 == k + 1) {
                                dzf1 = 1.0 - dzf;
                                dzc1 = 1.0 - dzc;
                            }

                            delta_mass = mass_per_cell*dxc1*dyc1*dzc1;
                            delta_pu = i*mom_per_cell*dxf1*dyc1*dzc1;
                            delta_pv = j*mom_per_cell*dxc1*dyf1*dzc1;
                            delta_pw = k*mom_per_cell*dxc1*dyc1*dzf1;
                            delta_therm = therm_per_cell*dxc1*dyc1*dzc1;

                            dratio = d1[ic+i1][jc+j1][kc+k1] / (d1[ic+i1][jc+j1][kc+k1] + delta_mass);
                            d1[ic+i1][jc+j1][kc+k1] += delta_mass;
                            u1[iface+i1][jc+j1][kc+k1] += delta_pu;
                            v1[ic+i1][jface+j1][kc+k1] += delta_pv;
                            w1[ic+i1][jc+j1][kface+k1] += delta_pw;
                            total_mass += delta_mass;

                            te1[ic+i1][jc+j1][kc+k1] = te1[ic+i1][jc+j1][kc+k1] * dratio + delta_therm / d1[ic+i1][jc+j1][kc+k1];

                            if (DualEnergyFormalism == 1) {
                                ge1[ic+i1][jc+j1][kc+k1] = ge1[ic+i1][jc+j1][kc+k1] * dratio + delta_therm / d1[ic+i1][jc+j1][kc+k1];
                            }

                            if (MetallicityField == 1) {
                                metal1[ic+i1][jc+j1][kc+k1] += (delta_mass/m_eject) * (yield * (1.0 - metalf) + m_eject * metalf);
                            }
                        }
                    }
                }
            }
        }
    }
}

void grid::add_feedback2(const int &nx, const int &ny, const int &nz,
    const int &ic, const int &jc, const int &kc, const int &iface, const int &jface, const int &kface,
    const float &dxf, const float &dyf, const float &dzf, 
    const float &dxc, const float &dyc, const float &dzc,
    const float &m_eject, const float &yield, const float &metalf,
    const float &mass_per_cell, const float mom_per_cell, const float therm_per_cell) {

    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum, MetalIaNum, MetalIINum;
    if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
        ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

    int MetallicityField = (MetalNum != -1 || SNColourNum != -1);

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
		Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}
		
    float delta_mass, delta_pu, delta_pv, delta_pw, delta_therm, dratio;
    float dxf1, dyf1, dzf1, dxc1, dyc1, dzc1;
    int index1, index2, index3, index4;

    if (MultiMetals == 1) {
        ENZO_FAIL("momentum: not supported");
    }

    float total_mass = 0.0;

    for (int k = -1; k <= 1; k++) {
        for (int j = -1; j <= 1; j++) {
            for (int i = -1; i <= 1; i++) {
                for (int i1 = i; i1 <= i + 1; i1++) {
                    dxf1 = dxf;
                    dxc1 = dxc;
                    if (i1 == i + 1) {
                        dxf1 = 1.0 - dxf;
                        dxc1 = 1.0 - dxc;
                    }
                    for (int j1 = j; j1 <= j + 1; j1++) {
                        dyf1 = dyf;
                        dyc1 = dyc;
                        if (j1 == j + 1) {
                            dyf1 = 1.0 - dyf;
                            dyc1 - 1.0 - dyc;
                        }
                        for (int k1 = k; k1 <= k + 1; k1++) {
                            dzf1 = dzf;
                            dzc1 = dzc;
                            if (k1 == k + 1) {
                                dzf1 = 1.0 - dzf;
                                dzc1 = 1.0 - dzc;
                            }

                            delta_mass = mass_per_cell*dxc1*dyc1*dzc1;
                            delta_pu = i*mom_per_cell*dxf1*dyc1*dzc1;
                            delta_pv = j*mom_per_cell*dxc1*dyf1*dzc1;
                            delta_pw = k*mom_per_cell*dxc1*dyc1*dzf1;
                            delta_therm = therm_per_cell*dxc1*dyc1*dzc1;

                            index1 = (ic + i1) + (jc + j1) * nx + (kc + k1) * nx * ny;
                            index2 = (iface + i1) + (jc + j1) * nx + (kc + k1) * nx * ny;
                            index3 = (ic + i1) + (jface + j1) * nx + (kc + k1) * nx * ny;
                            index4 = (ic + i1) + (jc + j1) * nx + (kface + k1) * nx * ny;

                            dratio = this->BaryonField[DensNum][index1] / (this->BaryonField[DensNum][index1] + delta_mass);
                            this->BaryonField[DensNum][index1] += delta_mass;
                            this->ParticleVelocity[0][index2] += delta_pu;
                            this->ParticleVelocity[1][index3] += delta_pv;
                            this->ParticleVelocity[2][index4] += delta_pw;
                            total_mass += delta_mass;

                            this->BaryonField[TENum][index1] = this->BaryonField[TENum][index1] * dratio + delta_therm / this->BaryonField[DensNum][index1];

                            if (DualEnergyFormalism == 1) {
                                this->BaryonField[GENum][index1] = this->BaryonField[GENum][index1] * dratio + delta_therm / this->BaryonField[DensNum][index1];
                            }

                            if (MetallicityField == 1) {
                                this->BaryonField[MetalNum][index1] += (delta_mass/m_eject) * (yield * (1.0 - metalf) + m_eject * metalf);
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif