/***********************************************************************
/
	/  Individual star galaxy for Enzo-Abyss-AEOS
/
	/  written by: Yongseok Jo & Eunwoo Chung
	/  date:       May, 2025
/
/  PURPOSE:
/
************************************************************************/
#undef INDIVIDUALSTAR
#ifdef INDIVIDUALSTAR
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "ProblemType.h"
#include "EventHooks.h"
#include "phys_constants.h"

#include "IndividualStarProperties.h"
#include "StellarYieldsRoutines.h"

#define VCIRC_TABLE_LENGTH 10000

void mt_init(unsigned_int seed);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
					 LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
			 float *TemperatureUnits, float *TimeUnits,
			 float *VelocityUnits, double *MassUnits, FLOAT Time);
inline int nlines(const char *fname);
char *ChemicalSpeciesBaryonFieldLabel(const int &atomic_number, int element_set = 1);

int nlines(const char *fname)
{

	FILE *fptr = fopen(fname, "r");
	int ch, n = 0;

	do
	{
		ch = fgetc(fptr);
		if (ch == '\n')
			n++;
	} while (ch != EOF);

	fclose(fptr);
	if (debug)
		fprintf(stderr, "Read %" ISYM " lines \n", n);
	return n;
}

class ProblemType_IsolatedAEOS;

class IsolatedAEOSGrid : public grid
{
	friend class ProblemType_IsolatedAEOS;
};

class ProblemType_IsolatedAEOS : public EnzoProblemType
{
private:

public:
	FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
	FLOAT CenterPosition[MAX_DIMENSION];
	float Bfield[MAX_DIMENSION];
	FLOAT ScaleLength;
	FLOAT ScaleHeight;
	FLOAT ScaleLengthPlummer; // by Jo YS
	float DiskMass;
	float GasFraction;
	float DiskTemperature;
	float DiskMassPlummer; // by Jo YS
	float DiskMetallicity;
	float HaloMass;
	float HaloTemperature;
	float HaloMetallicity;
	FLOAT VCircRadius[VCIRC_TABLE_LENGTH], VCircRadius1[VCIRC_TABLE_LENGTH];	 // by Ahram Lee
	float VCircVelocity[VCIRC_TABLE_LENGTH], VCircVelocity1[VCIRC_TABLE_LENGTH]; // by Ahram Lee
	int RefineAtStart;
	int UseGasParticles;
	int UseGasParticlesEqualizePressure;
	FLOAT *GasParticlePosition[MAX_DIMENSION];
	FLOAT *GasParticleVelocity[MAX_DIMENSION];
	float *GasParticleMass;
	int NumberOfGasParticles;
	float GasHaloDensity;
	float GasHaloRadius;
	int isPlummer = 0;
	ProblemType_IsolatedAEOS() : EnzoProblemType()
	{
		if (MyProcessorNumber == 0)
			std::cout << "Creating problem type Agora Restart" << std::endl;
	}

	~ProblemType_IsolatedAEOS() {}

	virtual int InitializeFromRestart(
		HierarchyEntry &TopGrid, TopGridData &MetaData)
	{
		return SUCCESS;
	}

	virtual int InitializeSimulation(
		FILE *fptr, FILE *Outfptr,
		HierarchyEntry &TopGrid, TopGridData &MetaData)
	{
		if (debug)
		{
			printf("Entering IsolatedAEOSInitialize\n");
			fflush(stdout);
		}

		char *DensName = "Density";
		char *TEName = "TotalEnergy";
		char *GEName = "GasEnergy";
		char *Vel1Name = "x-velocity";
		char *Vel2Name = "y-velocity";
		char *Vel3Name = "z-velocity";
  	char *CRName      = "CREnergyDensity";
  	char *GravPotentialName = "GravPotential";

		char *ElectronName = "Electron_Density";
		char *HIName = "HI_Density";
		char *HIIName = "HII_Density";
		char *HeIName = "HeI_Density";
		char *HeIIName = "HeII_Density";
		char *HeIIIName = "HeIII_Density";

		char *HMName = "HM_Density";
		char *H2IName = "H2I_Density";
		char *H2IIName = "H2II_Density";
		char *DIName = "DI_Density";
		char *DIIName = "DII_Density";
		char *HDIName = "HDI_Density";
		char *MetalName = "Metal_Density";
		char *MetalIaName = "MetalSNIa_Density";

		char *MetalSNIaName = "MetalSNIa_Density";
		char *MetalSNIIName = "MetalSNII_Density";
		char *BxName = "Bx";
		char *ByName = "By";
		char *BzName = "Bz";
		char *PhiName = "Phi";

		char *AGBMetalName    = "AGB_Metal_Density";
		char *PopIIIMetalName = "PopIII_Metal_Density";
		char *PopIIIPISNeMetalName = "PopIII_PISNe_Metal_Density";
		char *SNIIMetalName   = "SNII_Metal_Density";
		char *SNIaMetalName   = "SNIa_Metal_Density";
		char *RProcMetalName  = "RProcess_Metal_Density";

		char *ExtraMetalName0    = "SNIa_sCh_Metal_Density";
		char *ExtraMetalName1    = "SNIa_SDS_Metal_Density";
		char *ExtraMetalName2    = "SNIa_HeRS_Metal_Density";

		/* local declarations */

		char line[MAX_LINE_LENGTH];
		int i, ret, level;

		/* make sure it is 3D */

		if (MetaData.TopGridRank != 3)
		{
			printf("Cannot do AcoraRestart in %" ISYM " dimension(s)\n",
				   MetaData.TopGridRank);
			ENZO_FAIL("Agora Restart simulations must be 3D!");
		}

		for (i = 0; i < MAX_DIMENSION; i++)
		{
			this->CenterPosition[i] = 0.5;
			this->Bfield[i] = 0.;
		}

		// These come from Oscar's sample output.  The units are:
		// Velocity: km/s
		// Mass: 10^9 Msun
		// Length: kpc
		// Temperature: K
		this->UseGasParticles = 0; // by default, do not use gas particles to init gas
		this->UseGasParticlesEqualizePressure = 1;
		this->ScaleLength = .0343218;
		this->ScaleHeight = .00343218;
		this->DiskMass = 42.9661;
		this->GasFraction = 0.2;
		this->DiskTemperature = 1e4;
		this->DiskMetallicity = 0.0;
		this->HaloMass = 0.10000;
		this->HaloTemperature = this->DiskTemperature;
		this->HaloMetallicity = 0.0;
		this->RefineAtStart = TRUE;

		this->NumberOfGasParticles = 0;
		this->GasParticleMass = NULL;
		for (int i = 0; i < MAX_DIMENSION; i++)
		{
			this->GasParticlePosition[i] = NULL;
			this->GasParticleVelocity[i] = NULL;
		}
		this->GasHaloDensity = 0.0;
		this->GasHaloRadius = 0.0; // ignore if zero
								   // set this from global data (kind of a hack)
		TestProblemData.MultiSpecies = MultiSpecies;

		/* read input from file */
		while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
		{
			ret = 0;
			ret += sscanf(line, "IsolatedAEOSUseGasParticles = %" ISYM,
						  &UseGasParticles);
			ret += sscanf(line, "IsolatedAEOSUseGasParticlesEqualizePressure = %" ISYM,
						  &UseGasParticlesEqualizePressure);
			ret += sscanf(line, "IsolatedAEOSCenterPosition = %" PSYM " %" PSYM " %" PSYM,
						  CenterPosition, CenterPosition + 1, CenterPosition + 2);
			ret += sscanf(line, "IsolatedAEOSScaleLength = %" PSYM, &ScaleLength);
			ret += sscanf(line, "IsolatedAEOSScaleLengthPlummer = %" PSYM, &ScaleLengthPlummer); // by YS Jo
			ret += sscanf(line, "IsolatedAEOSScaleHeight = %" PSYM, &ScaleHeight);
			ret += sscanf(line, "IsolatedAEOSDiskMass = %" FSYM, &DiskMass);
			ret += sscanf(line, "IsolatedAEOSDiskMassPlummer = %" FSYM, &DiskMassPlummer); // by YS Jo
			ret += sscanf(line, "isPlummer = %" ISYM, &isPlummer);						   // by YS Jo
			ret += sscanf(line, "IsolatedAEOSGasFraction = %" FSYM, &GasFraction);
			ret += sscanf(line, "IsolatedAEOSDiskTemperature = %" FSYM,
						  &DiskTemperature);
			ret += sscanf(line, "IsolatedAEOSDiskMetallicity = %" FSYM,
						  &DiskMetallicity);
			ret += sscanf(line, "IsolatedAEOSHaloMass = %" FSYM, &HaloMass);
			ret += sscanf(line, "IsolatedAEOSHaloTemperature = %" FSYM,
						  &HaloTemperature);
			ret += sscanf(line, "IsolatedAEOSHaloMetallicity = %" FSYM,
						  &HaloMetallicity);
			ret += sscanf(line, "IsolatedAEOSGasHaloDensity = %" FSYM,
						  &GasHaloDensity);
			ret += sscanf(line, "IsolatedAEOSGasHaloRadius = %" FSYM,
						  &GasHaloRadius);
			ret += sscanf(line, "IsolatedAEOSMagneticField = %" FSYM " %" FSYM " %" FSYM,
						  Bfield, Bfield + 1, Bfield + 2);

			ret += sscanf(line, "IsolatedAEOSRefineAtStart = %" ISYM,
						  &RefineAtStart);
			ret += sscanf(line, "IsolatedAEOSMultiMetals = %" ISYM,
						  &TestProblemData.MultiMetals);
			ret += sscanf(line, "IsolatedAEOSHydrogenFractionByMass = %" FSYM,
						  &TestProblemData.HydrogenFractionByMass);
			ret += sscanf(line, "IsolatedAEOSHeliumFractionByMass = %" FSYM,
						  &TestProblemData.HeliumFractionByMass);
			ret += sscanf(line, "IsolatedAEOSMetalFractionByMass = %" FSYM,
						  &TestProblemData.MetalFractionByMass);
			ret += sscanf(line, "IsolatedAEOSDeuteriumToHydrogenRatio = %" FSYM,
						  &TestProblemData.DeuteriumToHydrogenRatio);
			ret += sscanf(line, "IsolatedAEOSInitialHIFraction  = %" FSYM,
						  &TestProblemData.HI_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHIIFraction  = %" FSYM,
						  &TestProblemData.HII_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHeIFraction  = %" FSYM,
						  &TestProblemData.HeI_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHeIIFraction  = %" FSYM,
						  &TestProblemData.HeII_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHeIIIFraction  = %" FSYM,
						  &TestProblemData.HeIII_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHMFraction  = %" FSYM,
						  &TestProblemData.HM_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialH2IFraction  = %" FSYM,
						  &TestProblemData.H2I_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialH2IIFraction  = %" FSYM,
						  &TestProblemData.H2II_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialDIFraction  = %" FSYM,
						  &TestProblemData.DI_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialDIIFraction  = %" FSYM,
						  &TestProblemData.DII_Fraction);
			ret += sscanf(line, "IsolatedAEOSInitialHDIFraction  = %" FSYM,
						  &TestProblemData.HDI_Fraction);
			ret += sscanf(line, "IsolatedAEOSUseMetallicityField  = %" ISYM,
						  &TestProblemData.UseMetallicityField);


    /* Read in chemical tracer IC's */
    ret += sscanf(line, "IsolatedAEOSInitialCIFraction = %"FSYM,
                        &TestProblemData.CI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialNIFraction = %"FSYM,
                        &TestProblemData.NI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialOIFraction = %"FSYM,
                        &TestProblemData.OI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialMgIFraction = %"FSYM,
                        &TestProblemData.MgI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialSiIFraction = %"FSYM,
                        &TestProblemData.SiI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialFeIFraction = %"FSYM,
                        &TestProblemData.FeI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialYIFraction = %"FSYM,
                        &TestProblemData.YI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialBaIFraction = %"FSYM,
                        &TestProblemData.BaI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialLaIFraction = %"FSYM,
                        &TestProblemData.LaI_Fraction);
    ret += sscanf(line, "IsolatedAEOSInitialEuIFraction = %"FSYM,
                        &TestProblemData.EuI_Fraction);


    /* Initial Chemical tracer values for halo */
    ret += sscanf(line, "IsolatedAEOSInitialCIFractionHalo = %"FSYM,
                        &TestProblemData.CI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialNIFractionHalo = %"FSYM,
                        &TestProblemData.NI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialOIFractionHalo = %"FSYM,
                        &TestProblemData.OI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialMgIFractionHalo = %"FSYM,
                        &TestProblemData.MgI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialSiIFractionHalo = %"FSYM,
                        &TestProblemData.SiI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialFeIFractionHalo = %"FSYM,
                        &TestProblemData.FeI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialYIFractionHalo = %"FSYM,
                        &TestProblemData.YI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialBaIFractionHalo = %"FSYM,
                        &TestProblemData.BaI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialLaIFractionHalo = %"FSYM,
                        &TestProblemData.LaI_Fraction_2);
    ret += sscanf(line, "IsolatedAEOSInitialEuIFractionHalo = %"FSYM,
                        &TestProblemData.EuI_Fraction_2);

    ret += sscanf(line, "IsolatedAEOSInitialSpeciesFractionsDisk = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 0,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 1,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 2,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 3,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 4,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 5,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 6,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 7,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 8,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 9,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 10,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 11,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 12,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 13,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 14,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 15,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 16 );

    ret += sscanf(line, "IsolatedAEOSInitialSpeciesFractionsHalo = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 0,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 1,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 2,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 3,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 4,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 5,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 6,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 7,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 8,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 9,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 10,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 11,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 12,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 13,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 14,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 15,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 16 );


			if (ret == 0 && strstr(line, "=") &&
				(strstr(line, "IsolatedAEOS") || strstr(line, "TestProblem")) &&
				line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
				fprintf(stderr,
						"*** warning: the following parameter line from IsolatedAEOS was not interpreted:\n%s\n",
						line);

		} // end input from parameter file




		// Read in circular velocity table

		this->ReadInVcircData();

		if (UseGasParticles)
		{
			this->ReadInGasParticleData();
		}

		/* set up top grid */

		float dummy_density = 1.0;
		float dummy_gas_energy = 1.0; // Only used if DualEnergyFormalism is True
		float dummy_total_energy = 1.0;
		float dummy_velocity[3] = {0.0, 0.0, 0.0};
		float dummy_b_field[3] = {1e-20, 1e-20, 1e-20}; // Only set if HydroMethod = mhd_rk

		if (this->InitializeUniformGrid(
				TopGrid.GridData, TopGrid, MetaData, dummy_density, dummy_total_energy,
				dummy_gas_energy, dummy_velocity, dummy_b_field, this) == FAIL)
		{
			ENZO_FAIL("Error in InitializeUniformGrid");
		}

		/*
		if (UseGasParticles)
		{
			this->InitializeGridWithParticles(TopGrid.GridData, TopGrid, MetaData);
		}
		else
		{
			this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData);
		}*/

		this->InitializeParticles(TopGrid.GridData, TopGrid, MetaData);

		/* Convert minimum initial overdensity for refinement to mass
			 (unless MinimumMass itself was actually set). */

		if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED)
		{
			MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
			for (int dim = 0; dim < MetaData.TopGridRank; dim++)
				MinimumMassForRefinement[0] *= (DomainRightEdge[dim] - DomainLeftEdge[dim]) /
											   float(MetaData.TopGridDims[dim]);
		}

		/* If requested, refine the grid to the desired level. */

		if (RefineAtStart)
		{
			/* Declare, initialize, and fill out the first level of the LevelArray. */
			LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
			for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
				LevelArray[level] = NULL;
			AddLevel(LevelArray, &TopGrid, 0);

			/* Add levels to the maximum depth or until no new levels are created,
				 and re-initialize the level after it is created. */
			for (level = 0; level < MaximumRefinementLevel; level++)
			{
				if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL)
				{
					fprintf(stderr, "Error in RebuildHierarchy.\n");
					return FAIL;
				}
				if (LevelArray[level + 1] == NULL)
					break;
				LevelHierarchyEntry *Temp = LevelArray[level + 1];
				while (Temp != NULL)
				{
					if (UseGasParticles)
					{
						if (this->InitializeGridWithParticles(Temp->GridData, TopGrid, MetaData) == FAIL)
						{
							ENZO_FAIL("Error in AgoraReseart->InitializeGridWithParticles");
						}
					}
					else
					{
						if (this->InitializeGrid(Temp->GridData, TopGrid, MetaData) == FAIL)
						{
							ENZO_FAIL("Error in IsolatedAEOS->InitializeGrid");
						}
					}
					Temp = Temp->NextGridThisLevel;
				} // end: loop over grids on this level
			} // end: loop over levels
			if (UseGasParticles)
			{
				for (level = MaximumRefinementLevel; level > 0; level--)
				{
					LevelHierarchyEntry *Temp = LevelArray[level];
					while (Temp != NULL)
					{
						if (Temp->GridData->ProjectSolutionToParentGrid(
								*LevelArray[level - 1]->GridData) == FAIL)
						{
							fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
							return FAIL;
						}
						Temp = Temp->NextGridThisLevel;
					}
				}
			} // end use particles project

			if (debug1)
				fprintf(stderr, "Why1?\n"); // by YS
		}

		/* set up field names and units */
		int count = 0;
		DataLabel[count++] = DensName;
		DataLabel[count++] = Vel1Name;
		if (MetaData.TopGridRank > 1)
			DataLabel[count++] = Vel2Name;
		if (MetaData.TopGridRank > 2)
			DataLabel[count++] = Vel3Name;
		DataLabel[count++] = TEName;
		if (DualEnergyFormalism)
			DataLabel[count++] = GEName;

		if (HydroMethod == MHD_RK)
		{
			DataLabel[count++] = (char *)BxName;
			DataLabel[count++] = (char *)ByName;
			DataLabel[count++] = (char *)BzName;
			DataLabel[count++] = (char *)PhiName;
		}

		if (MultiSpecies)
		{
			DataLabel[count++] = ElectronName;
			DataLabel[count++] = HIName;
			DataLabel[count++] = HIIName;
			DataLabel[count++] = HeIName;
			DataLabel[count++] = HeIIName;
			DataLabel[count++] = HeIIIName;
			if (MultiSpecies > 1)
			{
				DataLabel[count++] = HMName;
				DataLabel[count++] = H2IName;
				DataLabel[count++] = H2IIName;
			}
			if (MultiSpecies > 2)
			{
				DataLabel[count++] = DIName;
				DataLabel[count++] = DIIName;
				DataLabel[count++] = HDIName;
			}
		}
		if (TestProblemData.UseMetallicityField)
			DataLabel[count++] = MetalName;
		if (StarMakerTypeIaSNe)
			DataLabel[count++] = MetalSNIaName;
		if (StarMakerTypeIISNeMetalField)
			DataLabel[count++] = MetalSNIIName;

		/* Chemical tracer set ups */
		if (TestProblemData.MultiMetals)
		{
			MultiMetals = TestProblemData.MultiMetals;
		}
		else if (MultiMetals)
		{
			TestProblemData.MultiMetals = MultiMetals;
		}

		if (TestProblemData.MultiMetals == 2)
		{

			for (int i = 0; i < StellarYieldsNumberOfSpecies; i++)
			{
				if (StellarYieldsAtomicNumbers[i] > 2)
				{
					DataLabel[count++] = ChemicalSpeciesBaryonFieldLabel(StellarYieldsAtomicNumbers[i]);
				}
			} // yields loop
		}
		for (i = 0; i < count; i++)
			DataUnits[i] = NULL;

		if (MyProcessorNumber == ROOT_PROCESSOR)
		{
			fprintf(Outfptr, "IsolatedAEOSUseGasParticles         = %" ISYM "\n",
					UseGasParticles);
			fprintf(Outfptr, "IsolatedAEOSUseGasParticlesEqualizePressure = %" ISYM "\n",
					UseGasParticlesEqualizePressure);
			fprintf(Outfptr, "IsolatedAEOSCenterPosition          = %" PSYM " %" PSYM " %" PSYM "\n",
					CenterPosition[0], CenterPosition[1], CenterPosition[2]);
			fprintf(Outfptr, "IsolatedAEOSMagneticField           = %" FSYM " %" FSYM " %" FSYM,
					Bfield[0], Bfield[1], Bfield[2]);
			fprintf(Outfptr, "IsolatedAEOSScaleLength             = %" PSYM "\n",
					ScaleLength);
			fprintf(Outfptr, "IsolatedAEOSScaleHeight             = %" PSYM "\n",
					ScaleHeight);
			fprintf(Outfptr, "IsolatedAEOSDiskMass                = %" FSYM "\n",
					DiskMass);
			/** by YS Jo **/
			fprintf(Outfptr, "isPlummer                           = %" ISYM "\n",
					isPlummer);
			fprintf(Outfptr, "IsolatedAEOSScaleLengthPlummer      = %" PSYM "\n",
					ScaleLengthPlummer);
			fprintf(Outfptr, "IsolatedAEOSDiskMassPlummer	        = %" FSYM "\n",
					DiskMassPlummer);
			fprintf(stdout, "isPlummer                           = %" ISYM "\n",
					isPlummer);
			fprintf(stdout, "IsolatedAEOSScaleLengthPlummer      = %" FSYM "\n",
					ScaleLengthPlummer);
			fprintf(stdout, "IsolatedAEOSDiskMassPlummer	        = %" FSYM "\n",
					DiskMassPlummer);
			/**          **/
			fprintf(Outfptr, "IsolatedAEOSGasFraction             = %" FSYM "\n",
					GasFraction);
			fprintf(Outfptr, "IsolatedAEOSDiskTemperature         = %" FSYM "\n",
					DiskTemperature);
			fprintf(Outfptr, "IsolatedAEOSHaloMass                = %" FSYM "\n",
					HaloMass);
			fprintf(Outfptr, "IsolatedAEOSHaloTemperature         = %" FSYM "\n",
					HaloTemperature);
			fprintf(Outfptr, "IsolatedAEOSGasHaloDensity          = %" FSYM "\n",
					GasHaloDensity);
			fprintf(Outfptr, "IsolatedAEOSGasHaloRadius           = %" FSYM "\n",
					GasHaloRadius);
			fprintf(Outfptr, "IsolatedAEOSRefineAtStart           = %" ISYM "\n",
					RefineAtStart);
			fprintf(Outfptr, "IsolatedAEOSHydrogenFractionByMass = %" FSYM "\n",
					TestProblemData.HydrogenFractionByMass);
			fprintf(Outfptr, "IsolatedAEOSHeliumFractionByMass = %" FSYM "\n",
					TestProblemData.HeliumFractionByMass);
			fprintf(Outfptr, "IsolatedAEOSMetalFractionByMass = %" FSYM "\n",
					TestProblemData.MetalFractionByMass);
			fprintf(Outfptr, "IsolatedAEOSInitialHIFraction  = %" FSYM "\n",
					TestProblemData.HI_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHIIFraction  = %" FSYM "\n",
					TestProblemData.HII_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHeIFraction  = %" FSYM "\n",
					TestProblemData.HeI_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHeIIFraction  = %" FSYM "\n",
					TestProblemData.HeII_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHeIIIIFraction  = %" FSYM "\n",
					TestProblemData.HeIII_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHMFraction  = %" FSYM "\n",
					TestProblemData.HM_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialH2IFraction  = %" FSYM "\n",
					TestProblemData.H2I_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialH2IIFraction  = %" FSYM "\n",
					TestProblemData.H2II_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialDIFraction  = %" FSYM "\n",
					TestProblemData.DI_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialDIIFraction  = %" FSYM "\n",
					TestProblemData.DII_Fraction);
			fprintf(Outfptr, "IsolatedAEOSInitialHDIFraction  = %" FSYM "\n",
					TestProblemData.HDI_Fraction);
			fprintf(Outfptr, "IsolatedAEOSUseMetallicityField  = %" ISYM "\n",
					TestProblemData.UseMetallicityField);
			fprintf(Outfptr, "IsolatedAEOSMultiMetals = %" ISYM "\n",
					TestProblemData.MultiMetals);
   // galaxy chemical tracers
   fprintf(Outfptr, "IsolatedAEOSInitialCIFraction = %"GOUTSYM"\n",
           TestProblemData.CI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialNIFraction = %"GOUTSYM"\n",
           TestProblemData.NI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialOIFraction = %"GOUTSYM"\n",
           TestProblemData.OI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialMgIFraction = %"GOUTSYM"\n",
           TestProblemData.MgI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialSiIFraction = %"GOUTSYM"\n",
           TestProblemData.SiI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialFeIFraction = %"GOUTSYM"\n",
           TestProblemData.FeI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialYIFraction = %"GOUTSYM"\n",
           TestProblemData.YI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialBaIFraction = %"GOUTSYM"\n",
           TestProblemData.BaI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialLaIFraction = %"GOUTSYM"\n",
           TestProblemData.LaI_Fraction);
   fprintf(Outfptr, "IsolatedAEOSInitialEuIFraction = %"GOUTSYM"\n",
           TestProblemData.EuI_Fraction);

   // halo chemical tracers
   fprintf(Outfptr, "IsolatedAEOSInitialCIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.CI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialNIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.NI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialOIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.OI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialMgIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.MgI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialSiIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.SiI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialFeIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.FeI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialYIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.YI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialBaIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.BaI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialLaIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.LaI_Fraction_2);
   fprintf(Outfptr, "IsolatedAEOSInitialEuIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.EuI_Fraction_2);
		}

		return SUCCESS;

	} // InitializeSimulation

	int InitializeGridWithParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid,
									TopGridData &MetaData)
	{
		// optional initialization scheme to read in gas directly from MakeDisk
		// and and them to cells with NO smoothing
		//
		// Cells with no deposited particles will be set to a low, constant density
		// with temperature such that they have the average pressure of adjacent cells.
		// (and average velocity).... ???? Make sure this is sensible.

		// need a flagging field to mark cells that don't host particles

		if (debug)
			printf("Entering IsolatedAEOS InitializeGridWithParticles\n");

		IsolatedAEOSGrid *thisgrid =
			static_cast<IsolatedAEOSGrid *>(thisgrid_orig);

		if (thisgrid->ProcessorNumber != MyProcessorNumber)
			return SUCCESS;

		/*
			 int nGas = 0, count = 0;
			 nGas = nlines("gas.dat");
			 if (debug) fprintf(stderr, "InitializeGridWithParticles: Number of Gas Particles %"ISYM"\n", nGas);

		// Initialize particle arrays and read gas particles
		PINT *Number = new PINT[nGas];
		int *Type = new int[nGas];
		FLOAT *Position[MAX_DIMENSION];
		float *Velocity[MAX_DIMENSION];
		for (int i = 0; i < thisgrid->GridRank; i++)
		{
		Position[i] = new FLOAT[nGas];
		Velocity[i] = new FLOAT[nGas];
		}
		float *Mass = new float[nGas];

		FLOAT dx = thisgrid->CellWidth[0][0];

		this->ReadParticlesFromFile(
		Number, Type, Position, Velocity, Mass,
		"gas.dat", -1, count, dx); // particle type is irrelevant here
		*/
		FLOAT dx = thisgrid->CellWidth[0][0];

		/* Get Units */
		float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
			  TemperatureUnits = 1;
		double MassUnits = 1;

		if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
					 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		/* Identify physical quantities */
		int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, MetalNum;

		int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
			DINum, DIINum, HDINum;

		if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
												 Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL)
		{
			fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
			ENZO_FAIL("");
		}

		if (TestProblemData.MultiSpecies)
			if (thisgrid->IdentifySpeciesFields(
					DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
					HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
				ENZO_FAIL("Error in grid->IdentifySpeciesFields.");

		int MetallicityField = FALSE;
		if ((MetalNum = FindField(
				 Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)) != -1)
			MetallicityField = TRUE;
		else
			MetalNum = 0;

		//
		int dim, i, j, k, n, size, index = 0, nx, ny, nz;
		float DiskGasEnergy, HaloGasEnergy, HaloDensity, BoxVolume, vcirc, mu;
		FLOAT xstart, ystart, zstart, xend, yend, zend, xp, yp, zp;

		nx = *(thisgrid->GridDimension);
		ny = *(thisgrid->GridDimension + 1);
		nz = *(thisgrid->GridDimension + 2);
		size = nx * ny * nz;

		// flagging field for empty cells
		int *iflag = new int[size];
		for (i = 0; i < size; i++)
			iflag[i] = 0;

		BoxVolume = 1.;
		for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
			BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

		/* Find the mean molecular weight */

		if (TestProblemData.MultiSpecies == FALSE)
			mu = Mu;
		else
		{
			// Atomic hydrogen
			mu = TestProblemData.HydrogenFractionByMass *
				 (TestProblemData.HI_Fraction + 2.0 * TestProblemData.HII_Fraction);

			// Helium
			mu += TestProblemData.HeliumFractionByMass / 4.0 *
				  (TestProblemData.HeI_Fraction + 2.0 * TestProblemData.HeII_Fraction +
				   3.0 * TestProblemData.HeIII_Fraction);

			// Molecular hydrogen, ignore Deuterium
			if (TestProblemData.MultiSpecies > 1)
				mu += TestProblemData.HydrogenFractionByMass / 2.0 *
					  (TestProblemData.H2I_Fraction + 2.0 * TestProblemData.H2II_Fraction);

			// Metals
			if (TestProblemData.UseMetallicityField)
				mu += TestProblemData.MetalFractionByMass / 16.0;

			mu = POW(mu, -1);
		}

		HaloGasEnergy = this->HaloTemperature / mu / (Gamma - 1) /
						TemperatureUnits;

		HaloDensity = this->HaloMass / BoxVolume;

		DiskGasEnergy = this->DiskTemperature / mu / (Gamma - 1) /
						TemperatureUnits;

		// loop over and deposit all particles
		int ibuff = NumberOfGhostZones;

		xstart = thisgrid->CellLeftEdge[0][0];
		ystart = thisgrid->CellLeftEdge[1][0];
		zstart = thisgrid->CellLeftEdge[2][0];
		xend = xstart + dx * nx;
		yend = ystart + dx * ny;
		zend = zstart + dx * nz;

		// initialize ALL cells to the DM halo background
		FLOAT x, y, z, radius;
		index = 0;
		for (k = 0; k < thisgrid->GridDimension[2]; k++)
		{
			for (j = 0; j < thisgrid->GridDimension[1]; j++)
			{
				for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
				{
					/* Compute position */

					x = (thisgrid->CellLeftEdge[0][i] + 0.5 * thisgrid->CellWidth[0][i]);
					y = (thisgrid->CellLeftEdge[1][j] + 0.5 * thisgrid->CellWidth[1][j]);
					z = (thisgrid->CellLeftEdge[2][k] + 0.5 * thisgrid->CellWidth[2][k]);

					x -= this->CenterPosition[0];
					y -= this->CenterPosition[1];
					z -= this->CenterPosition[2];

					radius = sqrt(POW(x, 2) +
								  POW(y, 2) +
								  POW(z, 2));

					//    for (index = 0; index < size; index++)
					//    {
					if (radius < GasHaloRadius)
					{
						thisgrid->BaryonField[DensNum][index] = GasHaloDensity * mu * mh / MassUnits * POW(LengthUnits, 3);
						thisgrid->BaryonField[TENum][index] = HaloGasEnergy;

						if (DualEnergyFormalism)
							thisgrid->BaryonField[GENum][index] = HaloGasEnergy;
					}
					else
					{
						thisgrid->BaryonField[DensNum][index] = GasHaloDensity * mu * mh / MassUnits * POW(LengthUnits, 3) / 1000.0;

						thisgrid->BaryonField[TENum][index] = HaloGasEnergy * 10;
						if (DualEnergyFormalism)
							thisgrid->BaryonField[GENum][index] = HaloGasEnergy * 10;
					}

					thisgrid->BaryonField[Vel1Num][index] = 0;
					thisgrid->BaryonField[Vel2Num][index] = 0;
					thisgrid->BaryonField[Vel3Num][index] = 0;

					if (TestProblemData.UseMetallicityField)
						thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
																 TestProblemData.MetalFractionByMass * HaloMetallicity;
				} // i
			} // j
		} // k
		  //    } // end initialize

		// Deposit gas particles and sum momentum in cells
		int number_on_grid = 0;
		for (n = 0; n < this->NumberOfGasParticles; n++)
		{

			if (this->GasParticlePosition[0][n] < xstart || // + ibuff*dx  ||
				this->GasParticlePosition[0][n] > xend ||	// - ibuff*dx  ||
				this->GasParticlePosition[1][n] < ystart || // + ibuff*dx  ||
				this->GasParticlePosition[1][n] > yend ||	// - ibuff*dx  ||
				this->GasParticlePosition[2][n] < zstart || // + ibuff*dx  ||
				this->GasParticlePosition[2][n] > zend)
			{			  // + ibuff*dx){
				continue; // particle is off of this grid
			}
			xp = (this->GasParticlePosition[0][n] - xstart) / dx;
			yp = (this->GasParticlePosition[1][n] - ystart) / dx;
			zp = (this->GasParticlePosition[2][n] - zstart) / dx;

			i = ((int)floor(xp));
			j = ((int)floor(yp));
			k = ((int)floor(zp));

			index = i + (j + k * ny) * nx;

			thisgrid->BaryonField[DensNum][index] = thisgrid->BaryonField[DensNum][index] * iflag[index] + this->GasParticleMass[n] / POW(dx, 3);
			// add momentum first, then divide by total this->GasParticleMass in cell later
			thisgrid->BaryonField[Vel1Num][index] += this->GasParticleMass[n] * this->GasParticleVelocity[0][n] / POW(dx, 3);
			thisgrid->BaryonField[Vel2Num][index] += this->GasParticleMass[n] * this->GasParticleVelocity[1][n] / POW(dx, 3);
			thisgrid->BaryonField[Vel3Num][index] += this->GasParticleMass[n] * this->GasParticleVelocity[2][n] / POW(dx, 3);

			iflag[index] = 1;
			number_on_grid++;
		} // end particle deposition
		fprintf(stderr, "Number of Particles on this Grid : %" ISYM "\n", number_on_grid);

		// correct units in disk density and velocity
		//   properly set energy and metal fraction for gas disk
		for (index = 0; index < size; index++)
		{
			if (iflag[index] == 0)
				continue;

			thisgrid->BaryonField[Vel1Num][index] /= (thisgrid->BaryonField[DensNum][index]);
			thisgrid->BaryonField[Vel2Num][index] /= (thisgrid->BaryonField[DensNum][index]);
			thisgrid->BaryonField[Vel3Num][index] /= (thisgrid->BaryonField[DensNum][index]);
			//      thisgrid->BaryonField[DensNum][index] /= (MassUnits*POW(dx,3));

			thisgrid->BaryonField[TENum][index] = DiskGasEnergy;
			if (HydroMethod != Zeus_Hydro)
			{
				thisgrid->BaryonField[TENum][index] += 0.5 *
													   (POW(thisgrid->BaryonField[Vel1Num][index], 2) +
														POW(thisgrid->BaryonField[Vel2Num][index], 2) +
														POW(thisgrid->BaryonField[Vel3Num][index], 2));
			}

			if (DualEnergyFormalism)
			{
				thisgrid->BaryonField[GENum][index] = DiskGasEnergy;
			}

			if (TestProblemData.UseMetallicityField)
			{
				thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
														 TestProblemData.MetalFractionByMass * DiskMetallicity;
			}
		}

		/// fill in cells near disk particles

		/* Now compute the temperature of everything */

		///
		/// do generic field initialization
		int xo, yo, zo;
		xo = 1;
		yo = nx;
		zo = nx * ny;
		index = 0;
		for (k = 0; k < thisgrid->GridDimension[2]; k++)
		{
			for (j = 0; j < thisgrid->GridDimension[1]; j++)
			{
				for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
				{

					if (iflag[index] == 0)
					{
						int xlow, xhigh, ylow, yhigh, zlow, zhigh;
						xlow = enzo_max(index - xo, 0);
						xhigh = enzo_min(index + xo, size - 1);
						ylow = enzo_max(index - yo, 0);
						yhigh = enzo_min(index + yo, size - 1);
						zlow = enzo_max(index - zo, 0);
						zhigh = enzo_min(index + zo, size - 1);

						float total_mass = 0.0; // density
						total_mass = iflag[xhigh] * thisgrid->BaryonField[DensNum][xhigh] +
									 iflag[xlow] * thisgrid->BaryonField[DensNum][xlow] +
									 iflag[yhigh] * thisgrid->BaryonField[DensNum][yhigh] +
									 iflag[ylow] * thisgrid->BaryonField[DensNum][ylow] +
									 iflag[zhigh] * thisgrid->BaryonField[DensNum][zhigh] +
									 iflag[zlow] * thisgrid->BaryonField[DensNum][zlow];

						if (total_mass > 0)
						{
							// average velocity with surrounding cells (mass weighted):
							thisgrid->BaryonField[Vel1Num][index] = (iflag[xhigh] * thisgrid->BaryonField[Vel1Num][xhigh] * thisgrid->BaryonField[DensNum][xhigh] +
																	 iflag[xlow] * thisgrid->BaryonField[Vel1Num][xlow] * thisgrid->BaryonField[DensNum][xlow] +
																	 iflag[yhigh] * thisgrid->BaryonField[Vel1Num][yhigh] * thisgrid->BaryonField[DensNum][yhigh] +
																	 iflag[ylow] * thisgrid->BaryonField[Vel1Num][ylow] * thisgrid->BaryonField[DensNum][ylow] +
																	 iflag[zhigh] * thisgrid->BaryonField[Vel1Num][zhigh] * thisgrid->BaryonField[DensNum][zhigh] +
																	 iflag[zlow] * thisgrid->BaryonField[Vel1Num][zlow] * thisgrid->BaryonField[DensNum][zlow]) /
																	total_mass;
							thisgrid->BaryonField[Vel2Num][index] = (iflag[xhigh] * thisgrid->BaryonField[Vel2Num][xhigh] * thisgrid->BaryonField[DensNum][xhigh] +
																	 iflag[xlow] * thisgrid->BaryonField[Vel2Num][xlow] * thisgrid->BaryonField[DensNum][xlow] +
																	 iflag[yhigh] * thisgrid->BaryonField[Vel2Num][yhigh] * thisgrid->BaryonField[DensNum][yhigh] +
																	 iflag[ylow] * thisgrid->BaryonField[Vel2Num][ylow] * thisgrid->BaryonField[DensNum][ylow] +
																	 iflag[zhigh] * thisgrid->BaryonField[Vel2Num][zhigh] * thisgrid->BaryonField[DensNum][zhigh] +
																	 iflag[zlow] * thisgrid->BaryonField[Vel2Num][zlow] * thisgrid->BaryonField[DensNum][zlow]) /
																	total_mass;
							thisgrid->BaryonField[Vel3Num][index] = (iflag[xhigh] * thisgrid->BaryonField[Vel3Num][xhigh] * thisgrid->BaryonField[DensNum][xhigh] +
																	 iflag[xlow] * thisgrid->BaryonField[Vel3Num][xlow] * thisgrid->BaryonField[DensNum][xlow] +
																	 iflag[yhigh] * thisgrid->BaryonField[Vel3Num][yhigh] * thisgrid->BaryonField[DensNum][yhigh] +
																	 iflag[ylow] * thisgrid->BaryonField[Vel3Num][ylow] * thisgrid->BaryonField[DensNum][ylow] +
																	 iflag[zhigh] * thisgrid->BaryonField[Vel3Num][zhigh] * thisgrid->BaryonField[DensNum][zhigh] +
																	 iflag[zlow] * thisgrid->BaryonField[Vel3Num][zlow] * thisgrid->BaryonField[DensNum][zlow]) /
																	total_mass;

							// compute pressure of surrounding cells (mass weighted)
							// and set temperature accordingly
							//     P   = rho * k * T / mu
							//
							//     T_new = P * mu / (rho * k)
							//
							//     removing mu and kboltz from computation (they get divided out anyway)
							if (UseGasParticlesEqualizePressure)
							{
								float pressure = this->DiskTemperature * (iflag[xhigh] * thisgrid->BaryonField[DensNum][xhigh] * thisgrid->BaryonField[DensNum][xhigh] + iflag[xlow] * thisgrid->BaryonField[DensNum][xlow] * thisgrid->BaryonField[DensNum][xlow] + iflag[yhigh] * thisgrid->BaryonField[DensNum][yhigh] * thisgrid->BaryonField[DensNum][yhigh] + iflag[ylow] * thisgrid->BaryonField[DensNum][ylow] * thisgrid->BaryonField[DensNum][ylow] + iflag[zhigh] * thisgrid->BaryonField[DensNum][zhigh] * thisgrid->BaryonField[DensNum][zhigh] + iflag[zlow] * thisgrid->BaryonField[DensNum][zlow] * thisgrid->BaryonField[DensNum][zlow]) / total_mass;

								thisgrid->BaryonField[TENum][index] = (pressure / thisgrid->BaryonField[DensNum][index]) /
																	  mu / (Gamma - 1) / TemperatureUnits;

								if (DualEnergyFormalism)
									thisgrid->BaryonField[GENum][index] = (pressure / thisgrid->BaryonField[DensNum][index]) /
																		  mu / (Gamma - 1) / TemperatureUnits;
							}

						} // endif cell is adjacent to any particle deposition

					} // equalize the pressure with the surroundings

					if (StarMakerTypeIaSNe)
					{
						int SNIaNum = FindField(MetalSNIaDensity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
						if (SNIaNum != -1)
						{
							thisgrid->BaryonField[SNIaNum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
						}
					}
					if (StarMakerTypeIISNeMetalField)
					{
						int SNIINum = FindField(MetalSNIIDensity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
						if (SNIINum != -1)
						{
							thisgrid->BaryonField[SNIINum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
						}
						else
						{
							ENZO_FAIL("Thought we would find a SNII field but did not.");
						}
					}

					if (TestProblemData.MultiSpecies)
					{
						thisgrid->BaryonField[HINum][index] = TestProblemData.HI_Fraction *
															  TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HIINum][index] = TestProblemData.HII_Fraction *
															   TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeINum][index] = TestProblemData.HeI_Fraction *
															   TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeIINum][index] = TestProblemData.HeII_Fraction *
																TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeIIINum][index] = TestProblemData.HeIII_Fraction *
																 TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						if (TestProblemData.MultiSpecies > 1)
						{
							thisgrid->BaryonField[HMNum][index] = TestProblemData.HM_Fraction *
																  TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

							thisgrid->BaryonField[H2INum][index] = 2 * TestProblemData.H2I_Fraction *
																   TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

							thisgrid->BaryonField[H2IINum][index] = 2 * TestProblemData.H2II_Fraction *
																	TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];
						}

						if (TestProblemData.MultiSpecies > 1)
							thisgrid->BaryonField[HIINum][index] -=
								(thisgrid->BaryonField[HMNum][index] + thisgrid->BaryonField[H2IINum][index] + thisgrid->BaryonField[H2INum][index]);

						// Electron "density" (remember, this is a factor of m_p/m_e scaled
						// from the 'normal' density for convenience) is calculated by
						// summing up all of the ionized species.  The factors of 0.25 and
						// 0.5 in front of HeII and HeIII are to fix the fact that we're
						// calculating mass density, not number density (because the
						// thisgrid->BaryonField values are 4x as heavy for helium for a single
						// electron)
						thisgrid->BaryonField[DeNum][index] = thisgrid->BaryonField[HIINum][index] +
															  0.25 * thisgrid->BaryonField[HeIINum][index] +
															  0.5 * thisgrid->BaryonField[HeIIINum][index];
						if (TestProblemData.MultiSpecies > 1)
							thisgrid->BaryonField[DeNum][index] += 0.5 * thisgrid->BaryonField[H2IINum][index] -
																   thisgrid->BaryonField[HMNum][index];

						// Set deuterium species (assumed to be a negligible fraction of the
						// total, so not counted in the conservation)
						if (TestProblemData.MultiSpecies > 2)
						{
							thisgrid->BaryonField[DINum][index] =
								CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HINum][index];
							thisgrid->BaryonField[DIINum][index] =
								CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HIINum][index];
							thisgrid->BaryonField[HDINum][index] = 0.75 *
																   CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[H2INum][index];
						}
					} // if(TestProblemData.MultiSpecies)

					/* set chemical tracers to small density */
					/* For now, init halo chemical tracers density to zero */
					if (TestProblemData.MultiMetals == 2)
					{
						for (int yield_i = 0; yield_i < StellarYieldsNumberOfSpecies; yield_i++)
						{
							if (StellarYieldsAtomicNumbers[yield_i] > 2)
							{
								float fraction = 0.0;
								int field_num = 0;

								thisgrid->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[yield_i]);
								fraction = tiny_number;

								//               for (i = 0; i  < size; i ++){
								thisgrid->BaryonField[field_num][index] = fraction * thisgrid->BaryonField[DensNum][index];
								//               }
							}
						} // end for loop
					} // end MM == 2 check

					if (HydroMethod == MHD_RK)
					{
						thisgrid->BaryonField[B1Num][index] = Bfield[0];
						thisgrid->BaryonField[B2Num][index] = Bfield[1];
						thisgrid->BaryonField[B3Num][index] = Bfield[2];

						thisgrid->BaryonField[TENum][index] +=
							0.5 * (POW(thisgrid->BaryonField[B1Num][index], 2) + POW(thisgrid->BaryonField[B2Num][index], 2) + POW(thisgrid->BaryonField[B3Num][index], 2)) / thisgrid->BaryonField[DensNum][index];
					}
				} // k
			} // j
		} // i

		for (index = 0; index < size; index++)
			thisgrid->BaryonField[TENum][index] = DiskGasEnergy; // debugging
																 //
																 // make sure to delete particle arrays from gas particles
																 // as well as flagging field
																 //
																 //    delete [] Mass;
																 //   delete [] Type;
		delete[] iflag;
		//    delete [] Number;
		//    for (i = 0; i < thisgrid->GridRank; i ++)
		//    {
		//      delete [] Position[i];
		//      Position[i] = NULL;
		//      delete [] Velocity[i];
		//      Velocity[i] = NULL;
		//    }
		//    Mass = NULL;
		//    Type = NULL;
		iflag = NULL;
		//    Number = NULL;
		//   Position = NULL;
		//   Velocity = NULL;

		return SUCCESS;

	} // InitializeGridWithParticles

	int InitializeGrid(grid *thisgrid_orig, HierarchyEntry &TopGrid,
					   TopGridData &MetaData)
	{

		if (debug)
			printf("Entering IsolatedAEOS InitializeGrid\n");

		IsolatedAEOSGrid *thisgrid =
			static_cast<IsolatedAEOSGrid *>(thisgrid_orig);

		if (thisgrid->ProcessorNumber != MyProcessorNumber)
			return SUCCESS;

		/* Get units */
		float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
			  TemperatureUnits = 1;
		double MassUnits = 1;

		if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
					 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		/* Identify physical quantities */
		int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, MetalNum;

		int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
			DINum, DIINum, HDINum;

    int CINum, NINum, OINum, MgINum, SiINum, FeINum, YINum, BaINum, LaINum, EuINum;

  int ExtraField[9] = {0};

	  float H_Fraction, HII_Fraction, HeII_Fraction, HeIII_Fraction, HM_Fraction,
        H2I_Fraction, H2II_Fraction, D_to_H_ratio;

	  float CI_Fraction, NI_Fraction, OI_Fraction, MgI_Fraction, SiI_Fraction, FeI_Fraction,
        YI_Fraction, LaI_Fraction, BaI_Fraction, EuI_Fraction, metal_fraction;

		if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
												 Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL)
		{
			fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
			ENZO_FAIL("");
		}

		if (TestProblemData.MultiSpecies)
			if (thisgrid->IdentifySpeciesFields(
					DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
					HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
				ENZO_FAIL("Error in grid->IdentifySpeciesFields.");

		int MetallicityField = FALSE;
		if ((MetalNum = FindField(
				 Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)) != -1)
			MetallicityField = TRUE;
		else
			MetalNum = 0;

		int dim, i, j, k, size, index = 0;
		float RhoZero, DiskGasEnergy, DiskDensity, HaloGasEnergy, HaloDensity,
			BoxVolume, vcirc, mu,
			RhoZeroPlummer; // by YS Jo
		FLOAT x, y, z, radius, xy_radius, cellwidth;

		/* Compute size of this grid */
		size = 1;
		for (dim = 0; dim < thisgrid->GridRank; dim++)
			size *= thisgrid->GridDimension[dim];
		cellwidth = thisgrid->CellWidth[0][0];

		/* Compute the size of the box */
		BoxVolume = 1.;
		for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
			BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

		/* Find the mean molecular weight */

		if (TestProblemData.MultiSpecies == FALSE)
			mu = Mu;
		else
		{
			// Atomic hydrogen
			mu = TestProblemData.HydrogenFractionByMass *
				 (TestProblemData.HI_Fraction + 2.0 * TestProblemData.HII_Fraction);

			// Helium
			mu += TestProblemData.HeliumFractionByMass / 4.0 *
				  (TestProblemData.HeI_Fraction + 2.0 * TestProblemData.HeII_Fraction +
				   3.0 * TestProblemData.HeIII_Fraction);

			// Molecular hydrogen, ignore Deuterium
			if (TestProblemData.MultiSpecies > 1)
				mu += TestProblemData.HydrogenFractionByMass / 2.0 *
					  (TestProblemData.H2I_Fraction + 2.0 * TestProblemData.H2II_Fraction);

			// Metals
			if (TestProblemData.UseMetallicityField)
				mu += TestProblemData.MetalFractionByMass / 16.0;

			mu = POW(mu, -1);
		}

		/* Find global physical properties */
		RhoZero = this->DiskMass * this->GasFraction / (4. * pi) /
				  (POW((this->ScaleLength), 2) * (this->ScaleHeight));

		if (isPlummer)
			RhoZeroPlummer = this->DiskMassPlummer / (4. * pi) / // by Jo YS
							 POW((this->ScaleLengthPlummer), 3) * 3;

		HaloGasEnergy = this->HaloTemperature / mu / (Gamma - 1) /
						TemperatureUnits;

		HaloDensity = this->HaloMass / BoxVolume;

		DiskGasEnergy = this->DiskTemperature / mu / (Gamma - 1) /
						TemperatureUnits;

		/* Loop over the mesh. */

		for (k = 0; k < thisgrid->GridDimension[2]; k++)
		{
			for (j = 0; j < thisgrid->GridDimension[1]; j++)
			{
				for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
				{
					/* Compute position */

					x = (thisgrid->CellLeftEdge[0][i] + 0.5 * thisgrid->CellWidth[0][i]) *
						LengthUnits;
					y = (thisgrid->CellLeftEdge[1][j] + 0.5 * thisgrid->CellWidth[1][j]) *
						LengthUnits;
					z = (thisgrid->CellLeftEdge[2][k] + 0.5 * thisgrid->CellWidth[2][k]) *
						LengthUnits;

					x -= this->CenterPosition[0] * LengthUnits;
					y -= this->CenterPosition[1] * LengthUnits;
					z -= this->CenterPosition[2] * LengthUnits;

					radius = sqrt(POW(x, 2) +
								  POW(y, 2) +
								  POW(z, 2));

					xy_radius = sqrt(POW(x, 2) +
									 POW(y, 2));

					/* Find disk density, halo density and internal energy */

					// by Jo YS
					if (isPlummer)
					{
						DiskDensity = plummer_in_gauss_mass(RhoZero, RhoZeroPlummer, x / LengthUnits, y / LengthUnits,
															z / LengthUnits, cellwidth) /
									  POW(cellwidth, 3);
					}
					else
					{

						DiskDensity = gauss_mass(RhoZero, x / LengthUnits, y / LengthUnits,
												 z / LengthUnits, cellwidth) /
									  POW(cellwidth, 3);
					}

					if (HaloDensity * HaloTemperature > DiskDensity * DiskTemperature)
					{
						thisgrid->BaryonField[DensNum][index] = HaloDensity;
						thisgrid->BaryonField[TENum][index] = HaloGasEnergy;
						if (DualEnergyFormalism)
							thisgrid->BaryonField[GENum][index] = HaloGasEnergy;

						thisgrid->BaryonField[Vel1Num][index] = 0;
						thisgrid->BaryonField[Vel2Num][index] = 0;
						thisgrid->BaryonField[Vel3Num][index] = 0;

						if (TestProblemData.UseMetallicityField)
							thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
																	 TestProblemData.MetalFractionByMass * HaloMetallicity;
					}
					else // Ok, we're in the disk
					{
						thisgrid->BaryonField[DensNum][index] = DiskDensity;

						vcirc = this->InterpolateVcircTable(xy_radius);

						thisgrid->BaryonField[Vel1Num][index] =
							-vcirc * y / xy_radius / VelocityUnits;
						thisgrid->BaryonField[Vel2Num][index] =
							vcirc * x / xy_radius / VelocityUnits;
						thisgrid->BaryonField[Vel3Num][index] = 0;

						thisgrid->BaryonField[TENum][index] = DiskGasEnergy;
						if (HydroMethod != Zeus_Hydro)
						{
							thisgrid->BaryonField[TENum][index] += 0.5 *
															   (POW(thisgrid->BaryonField[Vel1Num][index], 2) +
																	POW(thisgrid->BaryonField[Vel2Num][index], 2) +
																	POW(thisgrid->BaryonField[Vel3Num][index], 2));
						}

						if (DualEnergyFormalism)
						{
							thisgrid->BaryonField[GENum][index] = DiskGasEnergy;
						}

						if (TestProblemData.UseMetallicityField)
						{
							thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
																	 TestProblemData.MetalFractionByMass * DiskMetallicity;
						}
					}
					if (StarMakerTypeIaSNe)
					{
						int SNIaNum = FindField(MetalSNIaDensity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
						if (SNIaNum != -1)
						{
							thisgrid->BaryonField[SNIaNum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
						}
					}
					if (StarMakerTypeIISNeMetalField)
					{
						int SNIINum = FindField(MetalSNIIDensity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
						if (SNIINum != -1)
						{
							thisgrid->BaryonField[SNIINum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
						}
						else
						{
							ENZO_FAIL("Thought we would find a SNII field but did not.");
						}
					}

					if (TestProblemData.MultiSpecies)
					{
						thisgrid->BaryonField[HINum][index] = TestProblemData.HI_Fraction *
															  TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HIINum][index] = TestProblemData.HII_Fraction *
															   TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeINum][index] = TestProblemData.HeI_Fraction *
															   TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeIINum][index] = TestProblemData.HeII_Fraction *
																TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						thisgrid->BaryonField[HeIIINum][index] = TestProblemData.HeIII_Fraction *
																 TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

						if (TestProblemData.MultiSpecies > 1)
						{
							thisgrid->BaryonField[HMNum][index] = TestProblemData.HM_Fraction *
																  TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

							thisgrid->BaryonField[H2INum][index] = 2 * TestProblemData.H2I_Fraction *
																   TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

							thisgrid->BaryonField[H2IINum][index] = 2 * TestProblemData.H2II_Fraction *
																	TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];
						}

						if (TestProblemData.MultiSpecies > 1)
							thisgrid->BaryonField[HIINum][index] -=
								(thisgrid->BaryonField[HMNum][index] + thisgrid->BaryonField[H2IINum][index] + thisgrid->BaryonField[H2INum][index]);

						// Electron "density" (remember, this is a factor of m_p/m_e scaled
						// from the 'normal' density for convenience) is calculated by
						// summing up all of the ionized species.  The factors of 0.25 and
						// 0.5 in front of HeII and HeIII are to fix the fact that we're
						// calculating mass density, not number density (because the
						// thisgrid->BaryonField values are 4x as heavy for helium for a single
						// electron)
						thisgrid->BaryonField[DeNum][index] = thisgrid->BaryonField[HIINum][index] +
															  0.25 * thisgrid->BaryonField[HeIINum][index] +
															  0.5 * thisgrid->BaryonField[HeIIINum][index];

						if (TestProblemData.MultiSpecies > 1)
							thisgrid->BaryonField[DeNum][index] += 0.5 * thisgrid->BaryonField[H2IINum][index] -
																   thisgrid->BaryonField[HMNum][index];

						// Set deuterium species (assumed to be a negligible fraction of the
						// total, so not counted in the conservation)
						if (TestProblemData.MultiSpecies > 2)
						{
							thisgrid->BaryonField[DINum][index] =
								CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HINum][index];
							thisgrid->BaryonField[DIINum][index] =
								CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HIINum][index];
							thisgrid->BaryonField[HDINum][index] = 0.75 *
																   CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[H2INum][index];
						}
					} // if(TestProblemData.MultiSpecies)

					if (HydroMethod == MHD_RK)
					{
						thisgrid->BaryonField[B1Num][index] = Bfield[0];
						thisgrid->BaryonField[B2Num][index] = Bfield[1];
						thisgrid->BaryonField[B3Num][index] = Bfield[2];

						thisgrid->BaryonField[TENum][index] +=
							0.5 * (POW(thisgrid->BaryonField[B1Num][index], 2) + POW(thisgrid->BaryonField[B2Num][index], 2) + POW(thisgrid->BaryonField[B3Num][index], 2)) / thisgrid->BaryonField[DensNum][index];
					}

				} // i
			} // j
		} // k

		return SUCCESS;

	} // InitializeGrid

	void InitializeParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid,
							 TopGridData &MetaData)
	{
		IsolatedAEOSGrid *thisgrid =
			static_cast<IsolatedAEOSGrid *>(thisgrid_orig);

		mt_init(thisgrid->ID);

		if (debug)
			printf("Entering IsolatedAEOS InitializeParticles\n");

		// Determine the number of particles of each type
		int nBulge, nDisk, nHalo, nNbody, nParticles;
		nBulge = nlines("bulge.dat");
		if (debug)
			fprintf(stderr, "InitializeParticles: Number of Bulge Particles %" ISYM "\n", nBulge);
		nDisk = nlines("disk.dat");
		if (debug)
			fprintf(stderr, "InitializeParticles: Number of Disk Particles %" ISYM "\n", nDisk);
		nHalo = nlines("halo.dat");
		if (debug)
			fprintf(stderr, "InitializeParticles: Number of Halo Particles %" ISYM "\n", nHalo);
#ifdef NBODY
			nNbody = nlines(NbodyDir);
			if(debug) fprintf(stderr, "InitializeParticles: Number of Nbody Particles %"ISYM"\n", nNbody);
			nParticles = nBulge + nDisk + nHalo + nNbody;
#else
		nParticles = nBulge + nDisk + nHalo;
#endif
		if (debug)
			fprintf(stderr, "InitializeParticles: Total Number of Particles %" ISYM "\n", nParticles);

		// Initialize particle arrays
		PINT *Number = new PINT[nParticles];
		int *Type = new int[nParticles];
		FLOAT *Position[MAX_DIMENSION];
		float *Velocity[MAX_DIMENSION];
		float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
			  TemperatureUnits = 1;
		double MassUnits = 1;

		if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
					 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		for (int i = 0; i < thisgrid->GridRank; i++)
		{
			Position[i] = new FLOAT[nParticles];
			Velocity[i] = new float[nParticles];
		}
		float *Mass = new float[nParticles];
		float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
		for (int i = 0; i < NumberOfParticleAttributes; i++)
		{
			Attribute[i] = new float[nParticles];
			for (int j = 0; j < nParticles; j++)
			{
				Attribute[i][j] = FLOAT_UNDEFINED;
			}
		}
#ifdef NBODY
			// just for tests
			for (int j = 0; j < nParticles; j++) {
				// Attribute[0][j] = 0.0; // (SEVN Query) modified by EW 2025.4.3 // It should be changed if we consider restart case...
				Attribute[1][j] = StarMakerMinimumDynamicalTime*3.15e7/TimeUnits; // Dynamical time
				Attribute[2][j] = TestProblemData.MetalFractionByMass; // test by EW 2025.4.3
			}
#endif

			FLOAT dx = thisgrid->CellWidth[0][0];

			// Read them in and assign them as we go
			int count = 0;
			if (nBulge > 0)
				this->ReadParticlesFromFile(
					Number, Type, Position, Velocity, Mass,
					"bulge.dat", PARTICLE_TYPE_STAR, count, dx);
			if (nDisk > 0)
				this->ReadParticlesFromFile(
					Number, Type, Position, Velocity, Mass,
					"disk.dat", PARTICLE_TYPE_STAR, count, dx);
			if (nHalo > 0)
				this->ReadParticlesFromFile(
					Number, Type, Position, Velocity, Mass,
					"halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx);
			if (nNbody > 0)
				this->ReadParticlesFromFile(
					Number, Type, Position, Velocity, Mass,
					NbodyDir, PARTICLE_TYPE_NBODY, count, dx);

#ifdef NBODY
#ifdef SEVN
			// just for tests
			for (int j = 0; j < nParticles; j++)
				Attribute[NumberOfParticleAttributes-8+0][j] = Mass[j] / (SolarMass / MassUnits / dx / dx / dx); // InitialMass [Msol]
#endif
#endif


#ifdef NBODY_old

			this->ReadNbodyParticles(
					Number, Type, Position, Velocity, Mass, nNbody,
					NbodyDir, PARTICLE_TYPE_NBODY, count, dx);
			fprintf(stdout, "NbodyCluster = (%.3e,%.3e,%.3e), r2 = %.3e \n", 
					NbodyClusterPosition[0][0], NbodyClusterPosition[1][0], NbodyClusterPosition[2][0], NbodyClusterPosition[3][0]);
#endif

		thisgrid->SetNumberOfParticles(count);
		thisgrid->SetParticlePointers(Mass, Number, Type, Position,
									  Velocity, Attribute);
		MetaData.NumberOfParticles = count;
		if (debug)
			fprintf(stderr, "InitializeParticles: Set Number of Particles %" ISYM "\n", count);
	}

	/*** by YS Jo ***/
	float plummer_in_gauss_mass(
		float RhoZero, float RhoZeroPlummer, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT cellwidth)
	{
		// As of now, I hard-coded the Plummer model with 1e6 Msun by YS Jo.
		// Computes the total mass in a given cell by integrating the density
		// profile using 5-point Gaussian quadrature.
		// http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
		FLOAT EvaluationPoints[5] = {-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985};
		FLOAT Weights[5] = {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689};
		FLOAT xResult[5];
		FLOAT yResult[5];
		FLOAT r, z;
		float Mass_g = 0, Mass_p = 0;
		int i, j, k;

		// Cell mass for Exponential Disk by Jo YS
		for (i = 0; i < 5; i++)
		{
			xResult[i] = 0.0;
			for (j = 0; j < 5; j++)
			{
				yResult[j] = 0.0;
				for (k = 0; k < 5; k++)
				{
					r = sqrt((POW(xpos + EvaluationPoints[i] * cellwidth / 2.0, 2.0) +
							  POW(ypos + EvaluationPoints[j] * cellwidth / 2.0, 2.0)));
					z = fabs(zpos + EvaluationPoints[k] * cellwidth / 2.0);
					yResult[j] +=
						cellwidth / 2.0 * Weights[k] * RhoZero *
						PEXP(-r / this->ScaleLength) *
						PEXP(-fabs(z) / this->ScaleHeight);
				}
				xResult[i] += cellwidth / 2.0 * Weights[j] * yResult[j];
			}
			Mass_g += cellwidth / 2.0 * Weights[i] * xResult[i];
		}

		// Cell mass for Plummer Profile by Jo YS
		for (i = 0; i < 5; i++)
		{
			xResult[i] = 0.0;
			for (j = 0; j < 5; j++)
			{
				yResult[j] = 0.0;
				for (k = 0; k < 5; k++)
				{
					r = sqrt((POW(xpos + EvaluationPoints[i] * cellwidth / 2.0, 2.0) +
							  POW(ypos + EvaluationPoints[j] * cellwidth / 2.0, 2.0)));
					z = fabs(zpos + EvaluationPoints[k] * cellwidth / 2.0);
					yResult[j] +=
						cellwidth / 2.0 * Weights[k] * RhoZeroPlummer *
						POW(1 + (r * r + z * z) / (this->ScaleLengthPlummer * this->ScaleLengthPlummer), -5 / 2);
				}
				xResult[i] += cellwidth / 2.0 * Weights[j] * yResult[j];
			}
			Mass_p += cellwidth / 2.0 * Weights[i] * xResult[i];
		}

		// return the larger value by Jo YS
		if (Mass_g > Mass_p)
		{
			return Mass_g;
		}
		else
		{
			return Mass_p;
		}
		// return Mass;
	}

	float gauss_mass(
		float RhoZero, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT cellwidth)
	{
		// Computes the total mass in a given cell by integrating the density
		// profile using 5-point Gaussian quadrature.
		// http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
		FLOAT EvaluationPoints[5] = {-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985};
		FLOAT Weights[5] = {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689};
		FLOAT xResult[5];
		FLOAT yResult[5];
		FLOAT r, z;
		float Mass = 0;
		int i, j, k;

		for (i = 0; i < 5; i++)
		{
			xResult[i] = 0.0;
			for (j = 0; j < 5; j++)
			{
				yResult[j] = 0.0;
				for (k = 0; k < 5; k++)
				{
					r = sqrt((POW(xpos + EvaluationPoints[i] * cellwidth / 2.0, 2.0) +
							  POW(ypos + EvaluationPoints[j] * cellwidth / 2.0, 2.0)));
					z = fabs(zpos + EvaluationPoints[k] * cellwidth / 2.0);
					yResult[j] +=
						cellwidth / 2.0 * Weights[k] * RhoZero *
						PEXP(-r / this->ScaleLength) *
						PEXP(-fabs(z) / this->ScaleHeight);
				}
				xResult[i] += cellwidth / 2.0 * Weights[j] * yResult[j];
			}
			Mass += cellwidth / 2.0 * Weights[i] * xResult[i];
		}
		return Mass;
	}

	void ReadInGasParticleData(void)
	{

		int count = 0;
		this->NumberOfGasParticles = nlines("gas.dat");
		if (debug)
			fprintf(stderr, "ReadInGasParticleData: Number of Gas Particles %" ISYM "\n", this->NumberOfGasParticles);

		// Initialize particle arrays and read gas particles
		PINT *Number = new PINT[this->NumberOfGasParticles];
		int *Type = new int[this->NumberOfGasParticles];
		for (int i = 0; i < MAX_DIMENSION; i++)
		{
			this->GasParticlePosition[i] = new FLOAT[this->NumberOfGasParticles];
			this->GasParticleVelocity[i] = new FLOAT[this->NumberOfGasParticles];
		}
		this->GasParticleMass = new float[this->NumberOfGasParticles];

		this->ReadParticlesFromFile(
			Number, Type, this->GasParticlePosition, this->GasParticleVelocity, this->GasParticleMass,
			"gas.dat", -1, count, 1.0); // particle type is irrelevant here

		delete[] Number;
		delete[] Type;
		Number = NULL;
		Type = NULL;

		return;
	}

	void ReadInVcircData(void)
	{
		FILE *fptr;
		char line[MAX_LINE_LENGTH];
		int i = 0, ret;
		float vcirc;
		FLOAT rad;

		fptr = fopen("vcirc.dat", "r");

		while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
		{
			ret += sscanf(line, "%" PSYM " %" FSYM, &rad, &vcirc);
			this->VCircRadius[i] = rad * kpc_cm;  // 3.08567758e21 = kpc/cm
			this->VCircVelocity[i] = vcirc * 1e5; // 1e5 = (km/s)/(cm/s)
			i += 1;
		}

		// by Ahram Lee
		if (this->isPlummer)
		{
			fptr = fopen("vcirc_nsc.dat", "r");
			i = 0;
			while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
			{
				ret += sscanf(line, "%" PSYM " %" FSYM, &rad, &vcirc);
				this->VCircRadius1[i] = rad * kpc_cm;  // 3.08567758e21 = kpc/cm
				this->VCircVelocity1[i] = vcirc * 1e5; // 1e5 = (km/s)/(cm/s)
				i += 1;
			}
		}

		fclose(fptr);
	} // ReadInVcircData

	float InterpolateVcircTable(FLOAT radius)
	{
		int i;

		for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
		{
			//      fprintf(stderr, "xx R = %ESYM ,  VCmax = %ESYM\n",radius, this->VCircRadius[i] );

			if (radius < this->VCircRadius[i])
				break;
			if (i == 0)
			{
				/*** Plummer ***/
				if (this->isPlummer)
				{
					int j;
					for (j = 0; j < VCIRC_TABLE_LENGTH; j++)
					{
						if (radius < this->VCircRadius1[j])
							break;
					}
					if (j == 0)
						return (VCircVelocity1[j]) * (radius - VCircRadius1[0]) / VCircRadius1[0];
					else if (j == VCIRC_TABLE_LENGTH)
						return (VCircVelocity[i]) * (radius - VCircRadius[0]) / VCircRadius[0];

					return VCircVelocity1[j - 1] +
						   (VCircVelocity1[j] - VCircVelocity1[j - 1]) *
							   (radius - VCircRadius1[j - 1]) /
							   (VCircRadius1[j] - VCircRadius1[j - 1]);
				}

				else
				{
					return (VCircVelocity[i]) * (radius - VCircRadius[0]) / VCircRadius[0];
				}
			}
			else if (i == VCIRC_TABLE_LENGTH)
				ENZO_FAIL("Fell off the circular velocity interpolation table");

			// we know the radius is between i and i-1
			return VCircVelocity[i - 1] +
				   (VCircVelocity[i] - VCircVelocity[i - 1]) *
					   (radius - VCircRadius[i - 1]) /
					   (VCircRadius[i] - VCircRadius[i - 1]);
		}
	}

	int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
							  float *Velocity[], float *Mass, const char *fname,
							  Eint32 particle_type, int &c, FLOAT dx)
	{
		FILE *fptr;
		char line[MAX_LINE_LENGTH];
		int ret;
		FLOAT x, y, z;
		float vx, vy, vz;
		double mass;

		float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
			  TemperatureUnits = 1;
		double MassUnits = 1;

		if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
					 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		fptr = fopen(fname, "r");

			while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
			{
				ret +=
					sscanf(line,
							"%"PSYM" %"PSYM" %"PSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
							&x, &y, &z, &vx, &vy, &vz, &mass);

				Position[0][c] = x * kpc_cm / LengthUnits + this->CenterPosition[0];
				Position[1][c] = y * kpc_cm / LengthUnits + this->CenterPosition[1];
				Position[2][c] = z * kpc_cm / LengthUnits + this->CenterPosition[2];

				Velocity[0][c] = vx * km_cm / VelocityUnits;
				Velocity[1][c] = vy * km_cm / VelocityUnits;
				Velocity[2][c] = vz * km_cm / VelocityUnits;

				// Particle masses are actually densities.
				Mass[c] = mass * 1e9 * SolarMass / MassUnits / dx / dx / dx;
				Type[c] = particle_type;
				Number[c] = c;
				c++;
			}

		fclose(fptr);

		return c;
	} // ReadParticlesFromFile

	int ReadNbodyParticles(PINT *Number, int *Type, FLOAT *Position[],
						   float *Velocity[], float *Mass, int count, const char *fname,
						   Eint32 particle_type, int &c, FLOAT dx)
	{

		std::ifstream file(fname);

		if (!file.is_open())
		{
			std::cerr << "Error opening " << fname << "!" << std::endl;
			return 1;
		}

		std::string line;
		int line_number = 0;
		FLOAT x, y, z;
		float vx, vy, vz;
		double mass;

		float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
			  TemperatureUnits = 1;
		double MassUnits = 1;

		if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
					 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL)
		{
			ENZO_FAIL("Error in GetUnits.");
		}

		while (std::getline(file, line) && line_number < count)
		{
			std::istringstream iss(line); // Create a string stream from the line

			// Read the type ('position' or 'velocity') and the x, y, z components
			iss >> x >> y >> z >> vx >> vy >> vz >> mass;

			Position[0][c] = x * kpc_cm / LengthUnits + this->CenterPosition[0];
			Position[1][c] = y * kpc_cm / LengthUnits + this->CenterPosition[1];
			Position[2][c] = z * kpc_cm / LengthUnits + this->CenterPosition[2];

			Velocity[0][c] = vx * km_cm / VelocityUnits;
			Velocity[1][c] = vy * km_cm / VelocityUnits;
			Velocity[2][c] = vz * km_cm / VelocityUnits;

			// Particle masses are actually densities.
			Mass[c] = mass * 1e9 * SolarMass / MassUnits / dx / dx / dx;
			Type[c] = particle_type;
			Number[c] = c++;
		}

		file.clear();
		file.seekg(0, std::ios::beg);

		while (std::getline(file, line))
		{
			std::istringstream iss(line); // Create a string stream from the line

			float r;
			// Read the type ('position' or 'velocity') and the x, y, z components
			iss >> x >> y >> z >> r; // >> vx >> vy >> vz >> mass;

			NbodyClusterPosition[0] = x * kpc_cm / LengthUnits + this->CenterPosition[0];
			NbodyClusterPosition[1] = y * kpc_cm / LengthUnits + this->CenterPosition[1];
			NbodyClusterPosition[2] = z * kpc_cm / LengthUnits + this->CenterPosition[2];
			if (r < 0)
			{
				NbodyClusterPosition[3] = r;
			}
			else
			{
				NbodyClusterPosition[3] = (r * kpc_cm / LengthUnits) * (r * kpc_cm / LengthUnits);
			}
		}

		if (NbodyClusterPosition[3] < 0)
			isNbodyParticleIdentification = false;
		else
			isNbodyParticleIdentification = true;

		file.close();

		return 0;
	}
}; // class declaration

//.. register:
namespace
{
	EnzoProblemType_creator_concrete<ProblemType_IsolatedAEOS>
		reg_isoAOES("IsolatedAEOS");
}




#ifdef NEW_PROBLEM_TYPES
#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#include "ProblemType.h"

class ProblemType_IsolatedAEOS;

/*
class IsolatedAEOSGrid : public grid
{
	friend class ProblemType_IsolatedAEOS;
};*/
//int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number);
//int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number);
int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number, int element_set = 1);
int GetUnits(float *DensityUnits, float *LengthUnits,
			 float *TemperatureUnits, float *TimeUnits,
			 float *VelocityUnits, double *MassUnits, FLOAT Time);

EnzoProblemMap& get_problem_types()
{
    static EnzoProblemMap problem_type_map;
    return problem_type_map;
}

/* This takes a string, grabs the (static) plugin map defined above, and
   returns the plugin creator for that. */
EnzoProblemType *select_problem_type( std::string problem_type_name)
{
    EnzoProblemType_creator *ept_creator = get_problem_types()
            [ problem_type_name ];

    /* Simply throw an error if no such plugin exists... */

    if( !ept_creator )
    {   
        EnzoProblemMap mymap = get_problem_types();

        for (EnzoProblemMap::const_iterator it 
                = mymap.begin(); it != mymap.end(); ++it) {
          std::cout << "Available: " << it->first << std::endl;
        }
        ENZO_FAIL("Unknown output plug-in.");
    }

    EnzoProblemType *ptype = ept_creator->create();

    return ptype;
}

EnzoProblemType::EnzoProblemType() { this->DataLabelCount = 0; }

int EnzoProblemType::AddDataLabel(const char *FieldName) {
    /* We allocate a new copy of FieldName */
    /* Include NUL-terminator */
    int slen = strlen(FieldName) + 1;
    char *fcopy = new char[slen];
    strncpy(fcopy, FieldName, slen);
    DataLabel[this->DataLabelCount] = fcopy;
  std::cout << "Adding " << DataLabel[this->DataLabelCount] << " in place "
            << this->DataLabelCount << std::endl;
    return DataLabelCount++;
}

grid *EnzoProblemType::CreateNewUniformGrid(
    grid *ParentGrid, HierarchyEntry &TopGrid, TopGridData &MetaData, int Rank,
    int Dimensions[], FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles,
    float UniformDensity, float UniformTotalEnergy, float UniformInternalEnergy,
    float UniformVelocity[], float UniformBField[]) {
    grid *grid_data = new grid;
    grid_data->InheritProperties(ParentGrid);
  grid_data->PrepareGrid(Rank, Dimensions, LeftEdge, RightEdge, NumParticles);
  this->InitializeUniformGrid(grid_data, TopGrid, MetaData, UniformDensity,
                              UniformTotalEnergy, UniformInternalEnergy,
                              UniformVelocity, UniformBField);
    return grid_data;
}

int EnzoProblemType::InitializeUniformGrid(
    grid *tg, HierarchyEntry &TopGrid, TopGridData &MetaData,
    float UniformDensity, float UniformTotalEnergy, float UniformInternalEnergy,
    float UniformVelocity[], float UniformBField[],
    ProblemType_IsolatedAEOS *ag)

{
  /* declarations */
 
  int dim, i, size, field, GCM;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, B1Num, B2Num, B3Num, PhiNum;

  int CINum, CIINum, OINum, OIINum, SiINum, SiIINum, SiIIINum, CHINum, CH2INum, 
    CH3IINum, C2INum, COINum, HCOIINum, OHINum, H2OINum, O2INum;

  int ExtraField[9] = {0};

  int NINum, MgINum, FeINum, YINum, BaINum, LaINum, EuINum, MetalIaNum;

  float H_Fraction, HII_Fraction, HeII_Fraction, HeIII_Fraction, HM_Fraction,
      H2I_Fraction, H2II_Fraction, D_to_H_ratio;

  float CI_Fraction, NI_Fraction, OI_Fraction, MgI_Fraction, SiI_Fraction,
      FeI_Fraction, YI_Fraction, LaI_Fraction, BaI_Fraction, EuI_Fraction,
      metal_fraction;

  /* create fields */
 
  tg->NumberOfBaryonFields = 0;
  tg->FieldType[tg->NumberOfBaryonFields++] = Density;
  int vel = tg->NumberOfBaryonFields;
  tg->FieldType[tg->NumberOfBaryonFields++] = Velocity1;
  if (tg->GridRank > 1 || HydroMethod > 2)
    tg->FieldType[tg->NumberOfBaryonFields++] = Velocity2;
  if (tg->GridRank > 2 || HydroMethod > 2)
    tg->FieldType[tg->NumberOfBaryonFields++] = Velocity3;
  tg->FieldType[tg->NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    tg->FieldType[tg->NumberOfBaryonFields++] = InternalEnergy;
  if (UseMHD) {
    tg->FieldType[B1Num = tg->NumberOfBaryonFields++] = Bfield1;
    tg->FieldType[B2Num = tg->NumberOfBaryonFields++] = Bfield2;
    tg->FieldType[B3Num = tg->NumberOfBaryonFields++] = Bfield3;
  }
  if (HydroMethod == MHD_RK) {
    tg->FieldType[PhiNum = tg->NumberOfBaryonFields++] = PhiField;
    if (UsePoissonDivergenceCleaning) {
      tg->FieldType[tg->NumberOfBaryonFields++] = Phi_pField;
    }
  }

  int colorfields = tg->NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (TestProblemData.MultiSpecies || MultiSpecies) {
    tg->FieldType[DeNum = tg->NumberOfBaryonFields++] = ElectronDensity;
    tg->FieldType[HINum = tg->NumberOfBaryonFields++] = HIDensity;
    tg->FieldType[HIINum = tg->NumberOfBaryonFields++] = HIIDensity;
    tg->FieldType[HeINum = tg->NumberOfBaryonFields++] = HeIDensity;
    tg->FieldType[HeIINum = tg->NumberOfBaryonFields++] = HeIIDensity;
    tg->FieldType[HeIIINum = tg->NumberOfBaryonFields++] = HeIIIDensity;
    if (TestProblemData.MultiSpecies > 1 || MultiSpecies > 1) {
      tg->FieldType[HMNum = tg->NumberOfBaryonFields++] = HMDensity;
      tg->FieldType[H2INum = tg->NumberOfBaryonFields++] = H2IDensity;
      tg->FieldType[H2IINum = tg->NumberOfBaryonFields++] = H2IIDensity;
    }
    if (TestProblemData.MultiSpecies > 2 || MultiSpecies > 2) {
      tg->FieldType[DINum = tg->NumberOfBaryonFields++] = DIDensity;
      tg->FieldType[DIINum = tg->NumberOfBaryonFields++] = DIIDensity;
      tg->FieldType[HDINum = tg->NumberOfBaryonFields++] = HDIDensity;
    }
  }

  if (TestProblemData.UseMetallicityField)
    tg->FieldType[MetalNum = tg->NumberOfBaryonFields++] = Metallicity;

  if (StarMakerTypeIaSNe)
    tg->FieldType[MetalIaNum = tg->NumberOfBaryonFields++] = MetalSNIaDensity;

  // Initialize chemical species tracer fields
  if (MultiMetals == 2) {
    for (int yield_i = 0; yield_i < StellarYieldsNumberOfSpecies; yield_i++) {
      if (StellarYieldsAtomicNumbers[yield_i] > 2) {
        tg->FieldType[tg->NumberOfBaryonFields++] =
            ChemicalSpeciesBaryonFieldNumber(
                StellarYieldsAtomicNumbers[yield_i]);
        }
    } // loop over yeilds

    if (IndividualStarTrackAGBMetalDensity) {
      tg->FieldType[ExtraField[0] = tg->NumberOfBaryonFields++] = ExtraType0;
    }

    if (IndividualStarPopIIIFormation) {
      tg->FieldType[ExtraField[1] = tg->NumberOfBaryonFields++] = ExtraType1;
      tg->FieldType[ExtraField[2] = tg->NumberOfBaryonFields++] =
          MetalPISNeDensity;
    }

    if (IndividualStarTrackSNMetalDensity) {
      tg->FieldType[ExtraField[3] = tg->NumberOfBaryonFields++] =
          MetalSNIaDensity;
      if (IndividualStarSNIaModel == 2) {
        tg->FieldType[ExtraField[4] = tg->NumberOfBaryonFields++] =
            ExtraMetalField0;
        tg->FieldType[ExtraField[5] = tg->NumberOfBaryonFields++] =
            ExtraMetalField1;
        tg->FieldType[ExtraField[6] = tg->NumberOfBaryonFields++] =
            ExtraMetalField2;
      }
      tg->FieldType[ExtraField[7] = tg->NumberOfBaryonFields++] =
          MetalSNIIDensity;
    }

    if (IndividualStarRProcessModel) {
      tg->FieldType[ExtraField[8] = tg->NumberOfBaryonFields++] =
          MetalRProcessDensity;
    }
    }

  /* Return if this doesn't concern us. */
  if (tg->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* compute size of fields */
 
  size = 1;
  for (dim = 0; dim < tg->GridRank; dim++)
    size *= tg->GridDimension[dim];
 
  /* allocate fields */
 
  // for (field = 0; field < tg->NumberOfBaryonFields; field++)
  //   if (tg->BaryonField[field] == NULL)
  //     tg->BaryonField[field] = new float[size];
  tg->AllocateGrids();  

  /* set density, total energy */
 
  for (i = 0; i < size; i++) {
    tg->BaryonField[0][i] = UniformDensity;
    tg->BaryonField[1][i] = UniformTotalEnergy;
  }
 
  /* set velocities */
  for (dim = 0; dim < tg->GridRank; dim++)
    for (i = 0; i < size; i++)
      tg->BaryonField[vel + dim][i] = UniformVelocity[dim];
 
  /* Set internal energy if necessary. */
 
  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      tg->BaryonField[2][i] = UniformInternalEnergy;

  if (UseMHD) {
    for (dim = 0; dim < 3; dim++) 
      for (i = 0; i < size; i++)
        tg->BaryonField[B1Num + dim][i] = UniformBField[dim];
  }

  /* set density of color fields to user-specified values (if user doesn't
    specify, the defaults are set in SetDefaultGlobalValues.  Do some minimal
    amount of error checking to try to ensure charge conservation when
    appropriate */
  for (i = 0; i < size; i++) {

    // Set multispecies fields!
    // this attempts to set them such that species conservation is maintained,
    // using the method in CosmologySimulationInitializeGrid.C
    if (TestProblemData.UseMetallicityField) {
      for (i = 0; i < size; i++)
        tg->BaryonField[MetalNum][i] = tiny_number * UniformDensity;
    }

    if (TestProblemData.MultiSpecies) {

      tg->BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
                                   TestProblemData.HydrogenFractionByMass *
                                   UniformDensity;
 
      tg->BaryonField[HeIINum][i] =
          TestProblemData.HeII_Fraction * UniformDensity *
          (1.0 - TestProblemData.HydrogenFractionByMass);

      tg->BaryonField[HeIIINum][i] =
          TestProblemData.HeIII_Fraction * UniformDensity *
          (1.0 - TestProblemData.HydrogenFractionByMass);

      tg->BaryonField[HeINum][i] =
          (1.0 - TestProblemData.HydrogenFractionByMass) * UniformDensity -
	tg->BaryonField[HeIINum][i] - tg->BaryonField[HeIIINum][i];

      if (TestProblemData.MultiSpecies > 1) {
        tg->BaryonField[HMNum][i] =
            TestProblemData.HM_Fraction * tg->BaryonField[HIINum][i];

	tg->BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
                                     tg->BaryonField[0][i] *
                                     TestProblemData.HydrogenFractionByMass;

        tg->BaryonField[H2IINum][i] =
            TestProblemData.H2II_Fraction * 2.0 * tg->BaryonField[HIINum][i];
      }

      // HI density is calculated by subtracting off the various ionized
      // fractions from the total
      tg->BaryonField[HINum][i] = 
          TestProblemData.HydrogenFractionByMass * tg->BaryonField[0][i] -
          tg->BaryonField[HIINum][i];
      if (MultiSpecies > 1)
        tg->BaryonField[HINum][i] -=
            (tg->BaryonField[HMNum][i] + tg->BaryonField[H2IINum][i] +
             tg->BaryonField[H2INum][i]);

      // Electron "density" (remember, this is a factor of m_p/m_e scaled from
      // the 'normal' density for convenience) is calculated by summing up all
      // of the ionized species. The factors of 0.25 and 0.5 in front of HeII
      // and HeIII are to fix the fact that we're calculating mass density,
      // not number density (because the BaryonField values are 4x as heavy
      // for helium for a single electron)
      tg->BaryonField[DeNum][i] = tg->BaryonField[HIINum][i] +
                                  0.25 * tg->BaryonField[HeIINum][i] +
                                  0.5 * tg->BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
        tg->BaryonField[DeNum][i] +=
            0.5 * tg->BaryonField[H2IINum][i] - tg->BaryonField[HMNum][i];

      // Set deuterium species (assumed to be a negligible fraction of the
      // total, so not counted in the conservation)
      if (TestProblemData.MultiSpecies > 2) {
        tg->BaryonField[DINum][i] = TestProblemData.DeuteriumToHydrogenRatio *
                                    tg->BaryonField[HINum][i];
        tg->BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio *
                                     tg->BaryonField[HIINum][i];
        tg->BaryonField[HDINum][i] = 0.75 *
                                     TestProblemData.DeuteriumToHydrogenRatio *
                                     tg->BaryonField[H2INum][i];
      }

    } // if(TestProblemData.MultiSpecies)

    // metallicity fields (including 'extra' metal fields)
    if (TestProblemData.UseMetallicityField || MultiMetals) {
      tg->BaryonField[MetalNum][i] =
          TestProblemData.MetallicityField_Fraction * UniformDensity;

      if (TestProblemData.MultiMetals == 1 || MultiMetals == 1) {
        tg->BaryonField[ExtraField[0]][i] =
            TestProblemData.MultiMetalsField1_Fraction * UniformDensity;
        tg->BaryonField[ExtraField[1]][i] =
            TestProblemData.MultiMetalsField2_Fraction * UniformDensity;
      }

      if (TestProblemData.MultiMetals == 2 || MultiMetals == 2) {
        float fraction = 0.0;
        for (int yield_i = 0; yield_i < StellarYieldsNumberOfSpecies;
             yield_i++) {
          if (StellarYieldsAtomicNumbers[yield_i] > 2) {
            int field_num = 0;

            tg->IdentifyChemicalTracerSpeciesFieldsByNumber(
                field_num, StellarYieldsAtomicNumbers[yield_i]);
            fraction = tiny_number;

            for (i = 0; i < size; i++) {
            tg->BaryonField[field_num][i] = fraction * UniformDensity;
          }
        }
        }

        if (IndividualStarTrackAGBMetalDensity) {
          for (i = 0; i < size; i++)
            tg->BaryonField[ExtraField[0]][i] = tiny_number * UniformDensity;
        }

        if (IndividualStarPopIIIFormation) {
          for (i = 0; i < size; i++) {
            tg->BaryonField[ExtraField[1]][i] = tiny_number * UniformDensity;
            tg->BaryonField[ExtraField[2]][i] = tiny_number * UniformDensity;
          }
        }

        if (IndividualStarTrackSNMetalDensity) {
          for (i = 0; i < size; i++) {
            tg->BaryonField[ExtraField[3]][i] = tiny_number * UniformDensity;
            tg->BaryonField[ExtraField[7]][i] = tiny_number * UniformDensity;
          }

          if (IndividualStarSNIaModel == 2) {
            for (i = 0; i < size; i++) {
              tg->BaryonField[ExtraField[4]][i] = tiny_number * UniformDensity;
              tg->BaryonField[ExtraField[5]][i] = tiny_number * UniformDensity;
              tg->BaryonField[ExtraField[6]][i] = tiny_number * UniformDensity;
            }
          }
        }

        if (IndividualStarRProcessModel) {
          for (i = 0; i < size; i++) {
            tg->BaryonField[ExtraField[8]][i] = tiny_number * UniformDensity;
          }
        }

    } // if(TestProblemData.UseMetallicityField)
    } // for (i = 0; i < size; i++)

    /***********************
     ***  InitializeGrid ***
     ***********************/
    if (debug)
      printf("Entering IsolatedAEOS InitializeGrid\n");

    IsolatedAEOSGrid *thisgrid = static_cast<IsolatedAEOSGrid *>(tg);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    /* Get units */
    float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
          TemperatureUnits = 1;
    double MassUnits = 1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
                 &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,
        PhiNum, MetalNum;

    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                             Vel3Num, TENum, B1Num, B2Num,
                                             B3Num, PhiNum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    /*
                if (TestProblemData.MultiSpecies)
                        if (thisgrid->IdentifySpeciesFields(
                                        DeNum, HINum, HIINum, HeINum, HeIINum,
       HeIIINum, HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
                                ENZO_FAIL("Error in
       grid->IdentifySpeciesFields.");

                int MetallicityField = FALSE;
                if ((MetalNum = FindField(
                                 Metallicity, thisgrid->FieldType,
       thisgrid->NumberOfBaryonFields)) != -1) MetallicityField = TRUE; else
                        MetalNum = 0;
    */ // by YS Jo

    int dim, i, j, k, size, index = 0;
    float RhoZero, DiskGasEnergy, DiskDensity, HaloGasEnergy, HaloDensity,
        BoxVolume, vcirc, mu,
        RhoZeroPlummer; // by YS Jo
    FLOAT x, y, z, radius, xy_radius, cellwidth;

    float chemical_species_fraction[MAX_STELLAR_YIELDS] = {tiny_number};

    /* Compute size of this grid */
    size = 1;
    for (dim = 0; dim < thisgrid->GridRank; dim++)
      size *= thisgrid->GridDimension[dim];
    cellwidth = thisgrid->CellWidth[0][0];

    /* Compute the size of the box */
    BoxVolume = 1.;
    for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
      BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

    /* Find the mean molecular weight */

    if (TestProblemData.MultiSpecies == FALSE)
      mu = Mu;
    else {
      // Atomic hydrogen
      mu = TestProblemData.HydrogenFractionByMass *
           (TestProblemData.HI_Fraction + 2.0 * TestProblemData.HII_Fraction);

      // Helium
      mu +=
          TestProblemData.HeliumFractionByMass / 4.0 *
          (TestProblemData.HeI_Fraction + 2.0 * TestProblemData.HeII_Fraction +
           3.0 * TestProblemData.HeIII_Fraction);

      // Molecular hydrogen, ignore Deuterium
      if (TestProblemData.MultiSpecies > 1)
        mu += TestProblemData.HydrogenFractionByMass / 2.0 *
              (TestProblemData.H2I_Fraction +
               2.0 * TestProblemData.H2II_Fraction);

      // Metals
      if (TestProblemData.UseMetallicityField)
        mu += TestProblemData.MetalFractionByMass / 16.0;

      mu = POW(mu, -1);
    }

    /* Find global physical properties */
    RhoZero = ag->DiskMass * ag->GasFraction / (4. * pi) /
              (POW((ag->ScaleLength), 2) * (ag->ScaleHeight));

    if (ag->isPlummer)
      RhoZeroPlummer = ag->DiskMassPlummer / (4. * pi) / // by Jo YS
                       POW((ag->ScaleLengthPlummer), 3) * 3;

    HaloGasEnergy = ag->HaloTemperature / mu / (Gamma - 1) / TemperatureUnits;

    HaloDensity = ag->HaloMass / BoxVolume;

    DiskGasEnergy = ag->DiskTemperature / mu / (Gamma - 1) / TemperatureUnits;

    /* Loop over the mesh. */

    for (k = 0; k < thisgrid->GridDimension[2]; k++) {
      for (j = 0; j < thisgrid->GridDimension[1]; j++) {
        for (i = 0; i < thisgrid->GridDimension[0]; i++, index++) {
          /* Compute position */

          x = (thisgrid->CellLeftEdge[0][i] + 0.5 * thisgrid->CellWidth[0][i]) *
              LengthUnits;
          y = (thisgrid->CellLeftEdge[1][j] + 0.5 * thisgrid->CellWidth[1][j]) *
              LengthUnits;
          z = (thisgrid->CellLeftEdge[2][k] + 0.5 * thisgrid->CellWidth[2][k]) *
              LengthUnits;

          x -= ag->CenterPosition[0] * LengthUnits;
          y -= ag->CenterPosition[1] * LengthUnits;
          z -= ag->CenterPosition[2] * LengthUnits;

          radius = sqrt(POW(x, 2) + POW(y, 2) + POW(z, 2));

          xy_radius = sqrt(POW(x, 2) + POW(y, 2));

          /* Set default tracer fraction */
          metal_fraction = ag->HaloMetallicity;

          if (StellarYieldsScaledSolarInitialAbundances) {
            for (int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++) {
              if (StellarYieldsAtomicNumbers[ii] > 2) {
                chemical_species_fraction[ii] =
                    StellarYields_ScaledSolarMassFractionByNumber(
                        metal_fraction, StellarYieldsAtomicNumbers[ii]);
              }
            }
          } else {
            for (int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++) {
              if (StellarYieldsAtomicNumbers[ii] > 2) {
                chemical_species_fraction[ii] =
                    TestProblemData.ChemicalTracerSpecies_Fractions_2[ii];
              }
            }
          }

          /* Find disk density, halo density and internal energy */
          // by Jo YS
          if (ag->isPlummer) {
            DiskDensity = ag->plummer_in_gauss_mass(
                              RhoZero, RhoZeroPlummer, x / LengthUnits,
                              y / LengthUnits, z / LengthUnits, cellwidth) /
                          POW(cellwidth, 3);
          } else {

            DiskDensity =
                ag->gauss_mass(RhoZero, x / LengthUnits, y / LengthUnits,
                               z / LengthUnits, cellwidth) /
                POW(cellwidth, 3);
          }

          if (HaloDensity * ag->HaloTemperature >
              DiskDensity * ag->DiskTemperature) {
            thisgrid->BaryonField[DensNum][index] = HaloDensity;
            thisgrid->BaryonField[TENum][index] = HaloGasEnergy;
            if (DualEnergyFormalism)
              thisgrid->BaryonField[GENum][index] = HaloGasEnergy;

            thisgrid->BaryonField[Vel1Num][index] = 0;
            thisgrid->BaryonField[Vel2Num][index] = 0;
            thisgrid->BaryonField[Vel3Num][index] = 0;

            if (TestProblemData.UseMetallicityField)
              thisgrid->BaryonField[MetalNum][index] =
                  thisgrid->BaryonField[DensNum][index] *
                  TestProblemData.MetalFractionByMass * ag->HaloMetallicity;
          } else // Ok, we're in the disk
          {
            thisgrid->BaryonField[DensNum][index] = DiskDensity;

            vcirc = ag->InterpolateVcircTable(xy_radius);

            thisgrid->BaryonField[Vel1Num][index] =
                -vcirc * y / xy_radius / VelocityUnits;
            thisgrid->BaryonField[Vel2Num][index] =
                vcirc * x / xy_radius / VelocityUnits;
            thisgrid->BaryonField[Vel3Num][index] = 0;

            thisgrid->BaryonField[TENum][index] = DiskGasEnergy;
            if (HydroMethod != Zeus_Hydro) {
              thisgrid->BaryonField[TENum][index] +=
                  0.5 * (POW(thisgrid->BaryonField[Vel1Num][index], 2) +
                         POW(thisgrid->BaryonField[Vel2Num][index], 2) +
                         POW(thisgrid->BaryonField[Vel3Num][index], 2));
      }

            if (DualEnergyFormalism) {
              thisgrid->BaryonField[GENum][index] = DiskGasEnergy;
            }

            if (TestProblemData.UseMetallicityField) {
              thisgrid->BaryonField[MetalNum][index] =
                  thisgrid->BaryonField[DensNum][index] *
                  TestProblemData.MetalFractionByMass * ag->DiskMetallicity;
            }
          }
          if (StarMakerTypeIaSNe) {
            int SNIaNum = FindField(MetalSNIaDensity, thisgrid->FieldType,
                                    thisgrid->NumberOfBaryonFields);
            if (SNIaNum != -1) {
              thisgrid->BaryonField[SNIaNum][index] =
                  thisgrid->BaryonField[DensNum][index] * 1.0e-8;
            }
          }
          if (StarMakerTypeIISNeMetalField) {
            int SNIINum = FindField(MetalSNIIDensity, thisgrid->FieldType,
                                    thisgrid->NumberOfBaryonFields);
            if (SNIINum != -1) {
              thisgrid->BaryonField[SNIINum][index] =
                  thisgrid->BaryonField[DensNum][index] * 1.0e-8;
            } else {
              ENZO_FAIL("Thought we would find a SNII field but did not.");
            }
          }

          if (TestProblemData.MultiSpecies) {
            thisgrid->BaryonField[HINum][index] =
                TestProblemData.HI_Fraction *
                TestProblemData.HydrogenFractionByMass *
                thisgrid->BaryonField[DensNum][index];

            thisgrid->BaryonField[HIINum][index] =
                TestProblemData.HII_Fraction *
                TestProblemData.HydrogenFractionByMass *
                thisgrid->BaryonField[DensNum][index];

            thisgrid->BaryonField[HeINum][index] =
                TestProblemData.HeI_Fraction *
                TestProblemData.HeliumFractionByMass *
                thisgrid->BaryonField[DensNum][index];

            thisgrid->BaryonField[HeIINum][index] =
                TestProblemData.HeII_Fraction *
                TestProblemData.HeliumFractionByMass *
                thisgrid->BaryonField[DensNum][index];

            thisgrid->BaryonField[HeIIINum][index] =
                TestProblemData.HeIII_Fraction *
                TestProblemData.HeliumFractionByMass *
                thisgrid->BaryonField[DensNum][index];

            if (TestProblemData.MultiSpecies > 1) {
              thisgrid->BaryonField[HMNum][index] =
                  TestProblemData.HM_Fraction *
                  TestProblemData.HydrogenFractionByMass *
                  thisgrid->BaryonField[DensNum][index];

              thisgrid->BaryonField[H2INum][index] =
                  2 * TestProblemData.H2I_Fraction *
                  TestProblemData.HydrogenFractionByMass *
                  thisgrid->BaryonField[DensNum][index];

              thisgrid->BaryonField[H2IINum][index] =
                  2 * TestProblemData.H2II_Fraction *
                  TestProblemData.HydrogenFractionByMass *
                  thisgrid->BaryonField[DensNum][index];
            }

            if (TestProblemData.MultiSpecies > 1)
              thisgrid->BaryonField[HIINum][index] -=
                  (thisgrid->BaryonField[HMNum][index] +
                   thisgrid->BaryonField[H2IINum][index] +
                   thisgrid->BaryonField[H2INum][index]);

            // Electron "density" (remember, this is a factor of m_p/m_e
            // scaled from the 'normal' density for convenience) is calculated
            // by summing up all of the ionized species.  The factors of 0.25
            // and 0.5 in front of HeII and HeIII are to fix the fact that
            // we're calculating mass density, not number density (because the
            // thisgrid->BaryonField values are 4x as heavy for helium for a
            // single electron)
            thisgrid->BaryonField[DeNum][index] =
                thisgrid->BaryonField[HIINum][index] +
                0.25 * thisgrid->BaryonField[HeIINum][index] +
                0.5 * thisgrid->BaryonField[HeIIINum][index];

            if (TestProblemData.MultiSpecies > 1)
              thisgrid->BaryonField[DeNum][index] +=
                  0.5 * thisgrid->BaryonField[H2IINum][index] -
                  thisgrid->BaryonField[HMNum][index];

            // Set deuterium species (assumed to be a negligible fraction of
            // the total, so not counted in the conservation)
            if (TestProblemData.MultiSpecies > 2) {
              thisgrid->BaryonField[DINum][index] =
                  CoolData.DeuteriumToHydrogenRatio *
                  thisgrid->BaryonField[HINum][index];
              thisgrid->BaryonField[DIINum][index] =
                  CoolData.DeuteriumToHydrogenRatio *
                  thisgrid->BaryonField[HIINum][index];
              thisgrid->BaryonField[HDINum][index] =
                  0.75 * CoolData.DeuteriumToHydrogenRatio *
                  thisgrid->BaryonField[H2INum][index];
            }

          } // if(TestProblemData.MultiSpecies)

          if (MultiMetals == 2) {
            for (int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++) {
              if (StellarYieldsAtomicNumbers[ii] > 2) {
                int field_num;
                tg->IdentifyChemicalTracerSpeciesFieldsByNumber(
                    field_num, StellarYieldsAtomicNumbers[ii]);

                thisgrid->BaryonField[field_num][index] =
                    HaloDensity * chemical_species_fraction[ii];
              }
            }
            if (IndividualStarTrackAGBMetalDensity) {
              // for (i = 0; i < size; i ++)
              thisgrid->BaryonField[ExtraField[0]][index] =
                  tiny_number * HaloDensity;
            }

            if (IndividualStarPopIIIFormation) {
              // for (i = 0; i < size; i ++){
              thisgrid->BaryonField[ExtraField[1]][index] =
                  tiny_number * HaloDensity;
              thisgrid->BaryonField[ExtraField[2]][index] =
                  tiny_number * HaloDensity;
              //}
            }

            if (IndividualStarTrackSNMetalDensity) {
              // for (i = 0; i < size; i++){
              thisgrid->BaryonField[ExtraField[3]][index] =
                  tiny_number * HaloDensity;
              thisgrid->BaryonField[ExtraField[7]][index] =
                  tiny_number * HaloDensity;
              //}
              if (IndividualStarSNIaModel == 2) {
                // for (i = 0; i < size; i++){
                thisgrid->BaryonField[ExtraField[4]][index] =
                    tiny_number * HaloDensity;
                thisgrid->BaryonField[ExtraField[5]][index] =
                    tiny_number * HaloDensity;
                thisgrid->BaryonField[ExtraField[6]][index] =
                    tiny_number * HaloDensity;
                //}
              }
            }

            if (IndividualStarRProcessModel) {
              for (i = 0; i < size; i++) {
                thisgrid->BaryonField[ExtraField[8]][i] =
                    tiny_number * HaloDensity;
              }
            }
          } // end chemical tracer value set

          if (HydroMethod == MHD_RK) {
            thisgrid->BaryonField[B1Num][index] = ag->Bfield[0];
            thisgrid->BaryonField[B2Num][index] = ag->Bfield[1];
            thisgrid->BaryonField[B3Num][index] = ag->Bfield[2];

            thisgrid->BaryonField[TENum][index] +=
                0.5 *
                (POW(thisgrid->BaryonField[B1Num][index], 2) +
                 POW(thisgrid->BaryonField[B2Num][index], 2) +
                 POW(thisgrid->BaryonField[B3Num][index], 2)) /
                thisgrid->BaryonField[DensNum][index];
          }

        } // i
      } // j
    } // k
  }
  return SUCCESS;
}

void EnzoProblemType::FinalizeGrids(HierarchyEntry **RefLevels,
                                    HierarchyEntry &TopGrid,
                                    TopGridData &MetaData) {

  /* set up subgrids from level 1 to max refinement level -1 */
  
  int lev;
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (RefLevels[lev]->GridData->ProjectSolutionToParentGrid(
            *(RefLevels[lev - 1]->GridData)) == FAIL) {
      ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  /* set up the root grid */

  if (MaximumRefinementLevel > 0) {
    if (RefLevels[0]->GridData->ProjectSolutionToParentGrid(
            *(TopGrid.GridData)) == FAIL) {
      ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  } else if (this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData) ==
             FAIL) {
      ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
    }
}
#endif
#endif