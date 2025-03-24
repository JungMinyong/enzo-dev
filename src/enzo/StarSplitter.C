/***********************************************************************
/
/  GRID CLASS (SPLIT PARTICLES INTO CHILDREN PARTICLES)
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:  

/  PURPOSE: This routine splits particles into 13 (=12+1) children particles 
/           when requested.  See Kitsionas & Whitworth (2002) for the
/           technical details of particle splitting,  which was already 
/           implemented and used in SPH/Gadget.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#define NO_DEBUG_PS 

/* function prototypes */
int generate_plummer(int N, double **star, double Rh, double rtide,
                double xp, double yp, double zp, double up, double vp, double wp);

int generate_random(int N, double **star, double xp, double yp, double zp, double up, double vp, double wp);

int grid::StarSplitter(int np, int* nnp)
{
    // (1) Split particles into a fixed number of particles.
    const int Nsplit = 30;
    const double Rh = 0.3;    /* pc */
    const double rtide = 10.0;/* pc */
    int numinput = *nnp - np;
    int numoutput = 0;

    /* Allocate a temporary array (VLA) to store the number of splits per parent.
       (Requires C99 or later.) */
    double nsplits[numinput];

    /* Allocate the star array on the stack.
       We'll use a 2D array for the actual data and an array of pointers to its rows.
       Note: Nsplit is a compile-time constant (8). */
    double star_data[Nsplit][7];
    double *star[Nsplit];
    for (int i = 0; i < Nsplit; i++){
        star[i] = star_data[i];
    }

    /* Seed the random generator once before the loop */
    srand48((unsigned) time(NULL));

    /* Loop over each parent particle */
    for (int ii = 0; ii < numinput; ii++){
        // fprintf(stderr, "Inside StarSplitter... np: %d, nnp: %d\n", np, *nnp);
        // fprintf(stderr, "ii: %d\n", ii);
        // Initialize the mass for each child particle based on the parent's mass.
        for (int i = 0; i < Nsplit; i++){
            star[i][0] = this->ParticleMass[np+ii] / Nsplit;
        }
        // Call generate_plummer to fill star[][1..6]
        // generate_plummer(Nsplit, star, Rh, rtide, 
        //             this->ParticlePosition[0][np+ii],
        //             this->ParticlePosition[1][np+ii],
        //             this->ParticlePosition[2][np+ii],
        //             this->ParticleVelocity[0][np+ii],  // Correct velocity component
        //             this->ParticleVelocity[1][np+ii],
        //             this->ParticleVelocity[2][np+ii]);
        generate_random(Nsplit, star, 
                    this->ParticlePosition[0][np+ii],
                    this->ParticlePosition[1][np+ii],
                    this->ParticlePosition[2][np+ii],
                    this->ParticleVelocity[0][np+ii],  // Correct velocity component
                    this->ParticleVelocity[1][np+ii],
                    this->ParticleVelocity[2][np+ii]);

        // Write the generated child data into the grid's arrays.
        for (int i = 0; i < Nsplit; i++){
            int destIndex = np + numoutput + i;
            this->ParticlePosition[0][destIndex] = star[i][1];
            this->ParticlePosition[1][destIndex] = star[i][2];
            this->ParticlePosition[2][destIndex] = star[i][3];
            this->ParticleVelocity[0][destIndex] = star[i][4];
            this->ParticleVelocity[1][destIndex] = star[i][5];
            this->ParticleVelocity[2][destIndex] = star[i][6];
            this->ParticleMass[destIndex]       = star[i][0];
        }
        nsplits[ii] = Nsplit;  // Save the number of splits for this parent
        numoutput += Nsplit;
        // fprintf(stderr, "\tnumoutput: %d\n", numoutput);
    } // end of parent loop

    *nnp = np + numoutput;

    // Broadcast parent's ParticleAttribute values into the child particles.
    // Process in reverse order to avoid overwriting parent's attributes.
    for (int ii = numinput - 1; ii >= 0; ii--){
        fprintf(stderr, "ii: %d\n", ii);
        int Nsplitcurr = nsplits[ii];  // (Currently always equal to Nsplit.)
        numoutput -= Nsplitcurr;
        for (int i = 0; i < Nsplitcurr; i++){
            for (int j = 0; j < NumberOfParticleAttributes; j++){
                this->ParticleAttribute[j][np + numoutput + i] =
                    this->ParticleAttribute[j][np + ii];
            }
        }
        // fprintf(stderr, "\tnumoutput: %d\n", numoutput);
    } // end of reverse attribute broadcast loop


    // fflush(stderr);
    return 0;
}

int generate_random(int N, double **star, double xp, double yp, double zp, double up, double vp, double wp)
{
    for (int i = 0; i < N; i++){
        star[i][1] = 2.0 * 1e-6 * drand48() - 1e-6 + xp;  // x
        star[i][2] = 2.0 * 1e-6 * drand48() - 1e-6 + yp;  // y
        star[i][3] = 2.0 * 1e-6 * drand48() - 1e-6 + zp;  // z
        star[i][4] = 2.0 * 1e-6 * drand48() - 1e-6 + up;  // vx
        star[i][5] = 2.0 * 1e-6 * drand48() - 1e-6 + vp;  // vy
        star[i][6] = 2.0 * 1e-6 * drand48() - 1e-6 + wp;  // vz
    }
    return 0;
}

// #define TWOPI (2.0 * pi_val)
#define TWOPI   8.0*atan(1.0)   /* 2PI */
/***************************************************************************
 *   Copyright (C) 2009 by Andreas Kuepper                                 *
 *   akuepper@astro.uni-bonn.de                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

int generate_plummer(int N, double **star, double Rh, double rtide,
            double xp, double yp, double zp, double up, double vp, double wp)
 {
    double pc = 3.08567758149137e18;
    double yr = 3.1536e7;
    Rh = Rh * pc;  // Convert to cgs units
    rtide = rtide * pc;  // Convert to cgs units
    
    int i;
    double a = Rh / 1.305;  // Plummer scale length in pc
    rtide = rtide;  // Tidal radius in cgs units
    double r, theta, phi;
    double v, v_esc;
    double trial, f, f_max;
    
    for(i = 0; i < N; i++){
        /* --- Generate Position --- */
        do {
            double X;
            /* Generate a random number X uniformly in (0,1), avoiding X==0 */
            do { X = drand48(); } while(X < 1e-10);
            /* Inversion of the Plummer cumulative mass distribution:
             *
             *   r = a / sqrt( X^(-2/3) - 1 )
             */
            r = a / sqrt(pow(X, -2.0/3.0) - 1.0);
            /* Generate isotropic angles */
            theta = acos(1.0 - 2.0 * drand48());
            phi   = TWOPI * drand48();
            /* Compute Cartesian coordinates */
            star[i][1] = r * sin(theta) * cos(phi) + xp;
            star[i][2] = r * sin(theta) * sin(phi) + yp;
            star[i][3] = r * cos(theta) + zp;
        } while(r > rtide);  // Reject stars beyond the tidal cutoff
        
        /* --- Generate Velocity --- */
        /* The escape velocity at radius r in a Plummer model is:
         *
         *   v_esc = sqrt(2) * [1 + (r/a)^2]^(-1/4)
         */
        v_esc = sqrt(2.0) * pow(1.0 + (r * r) / (a * a), -0.25);
        /* Use rejection sampling to draw a speed from:
         *   f(v) ~ v^2 * (1 - (v/v_esc)^2)^(7/2)
         *
         * We need an envelope function that is always larger than f(v).
         * (Here an empirical constant is chosen; you may need to tweak it for efficiency.)
         */
        f_max = 0.1 * v_esc * v_esc;  // Empirical envelope constant
        for(;;) {
            double v_candidate = drand48() * v_esc;  // Candidate speed in [0, v_esc]
            double y = drand48() * f_max;
            f = v_candidate * v_candidate * pow(1.0 - (v_candidate * v_candidate) / (v_esc * v_esc), 7.0/2.0);
            if(y < f) { 
                v = v_candidate;
                break;
            }
        }
        /* Assign a random velocity direction (isotropic) */
        theta = acos(1.0 - 2.0 * drand48());
        phi   = TWOPI * drand48();
        star[i][4] = v * sin(theta) * cos(phi) + up;  // vx
        star[i][5] = v * sin(theta) * sin(phi) + vp;  // vy
        star[i][6] = v * cos(theta) + wp;             // vz
    }
    
    return 0;
}


