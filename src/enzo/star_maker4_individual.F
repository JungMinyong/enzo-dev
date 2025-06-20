#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker4_individual(nx, ny, nz,
     &                      d, dm, temp, u, v, w,
     &                      dt, r, metal, dx, t, z, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa,
     &                      initial_mass)

c
c  CREATES STAR PARTICLES
c
c  written by:  Brian O'Shea
c  date:  29 March 2005
c    this version of star_maker creates star particles using the method
c    of Andrey Kravtsov as described in "On the origin of the global
c    Schmidt law of star formation", 2003 ApJ 590:L1-L4.  As in this paper,
c    I hold the star formation timescale constant and have a fixed 
c    overdensity threshold.
c  
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold (some number * avg. density)
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    nmax     - particle array size specified by calling routine
c    metalfSNIa - metallicity fraction of particle (from SN Ia) ! MKRJ
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC imetalSNIa
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      R_PREC    metalSNIa(nx,ny,nz), metalfSNIa(nmax)
      R_PREC    odthresh, masseff, smthresh, mintdyn
      R_PREC    initial_mass(nmax)
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      R_PREC   div, tdyn, dtot
      R_PREC   sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   densthresh, timeconstant, gasfrac
      R_PREC   x_rand, y_rand, z_rand
      R_PREC   possible_mass, temp_starmass
      parameter (sndspdC=1.3095e8_RKIND)
c
      ii = np
c
c  calculate density threshold.  odthresh is in proper particles per cc
c  and (d1/mass_h) gives the mean density of the universe in 
c  particles/cm^3 (assuming hydrogen is dominant)
c
c
      densthresh = odthresh / (d1 / mass_h)

c
c  calculate time constant for star formation.  This assumes that the
c  user input is in units of years
c

c      timeconstant = mintdyn * yr_s / t1
c
c  for each zone, : "star" particle is created if the density exceeds
c  some threshold and this is the highest level of refinement.  That's
c  it.
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) is this finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) is density greater than threshold?
c
               if (d(i,j,k) .lt. densthresh) goto 10

c            make sure that we never give put than 90% of the cell's mass
c            into a star particle
c           ##### CHANGED TO BELOW FOR 080112_AGORA-DISK-ENZO (t_ff INSTEAD OF mintdyn)
c               gasfrac = min( 0.9_RKIND, dt / timeconstant )
c
c           ##### The line below is inserted to be within the loop (cut from a few lines above)
               timeconstant = mintdyn * 3.156e7_RKIND / t1

               dtot = ( d(i,j,k) )*d1
               tdyn = sqrt(3._RKIND*pi_val/32._RKIND/GravConst/dtot)/t1
               timeconstant = max(tdyn, timeconstant)
               gasfrac = min( 0.9_RKIND, masseff * dt / timeconstant )
c            calculate star mass in solar masses.  If this is less than the
c            user-defined threshold mass, do NOT make a star in this cell.
c            This is not exactly in keeping with the spirit of the Kravtsov
c            algorithm, and is somewhat degenerate with the density threshold,
c            but we really don't want millions and millions of star particles.

               starmass = gasfrac*d(i,j,k)*dble(d1)*dble(x1*dx)**3
     &                   / SolarMass

#define STOCHASTIC_STAR_FORMATION  
c
#ifdef STOCHASTIC_STAR_FORMATION
c
c                 Keep global count of "unfullfilled" star formation
c                 and when total is larger than threshold, then create
c                 a star particle with the threshold mass or 1/2 the
c                 gas in the cell, whichever is smaller.
c
               if (starmass .lt. smthresh) then
                  sformsum = sformsum + starmass
                  if (sformsum .lt. smthresh) goto 10
                  bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / SolarMass
                  gasfrac = min(smthresh/bmass, 0.5_RKIND)

c                 below is inserted to remove stars less than threshold         
                  if (gasfrac .lt. smthresh/bmass) goto 10

                  sformsum = sformsum - gasfrac*bmass
c                 in Stochastic star formation, starmass becomes
c                 gasfrac*bmass
                  starmass = gasfrac*bmass
               endif
#else
               if(starmass .lt. smthresh) goto 10
#endif

c
c              If both of these criteria are met, create a star particle
c
c               possible_mass = starfraction * bmass
               possible_mass = starmass
               call sampling_Kroupa_IMF2(temp_starmass)
               do while (possible_mass .ge. temp_starmass)
c
c                 Do not generate more star particles than available
                  if (ii .eq. nmax) goto 20
c
                  possible_mass = possible_mass - temp_starmass
                  ii = ii + 1
c  original             mp(ii)  = starfraction * d(i,j,k)
                  mp(ii) = temp_starmass
                  initial_mass(ii) = temp_starmass
                  mp(ii) = mp(ii) * SolarMass
     &                  / (dble(d1)*dble(x1*dx)**3)
                  tcp(ii) = t
                  tdp(ii) = tdyn
                  call random_number(x_rand)
                  xp(ii) = xstart + 
     &                        (REAL(i,RKIND)-1.0_RKIND+x_rand)*dx
                  call random_number(y_rand)
                  yp(ii) = ystart + 
     &                        (REAL(j,RKIND)-1.0_RKIND+y_rand)*dx
                  call random_number(z_rand)
                  zp(ii) = zstart + 
     &                        (REAL(k,RKIND)-1.0_RKIND+z_rand)*dx

                  if (imethod .eq. 2) then
                     up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                     vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                     wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
                  else
                     up(ii) = u(i,j,k)
                     vp(ii) = v(i,j,k)
                     wp(ii) = w(i,j,k)
                  endif
                  call gaussian_random_vel(up(ii), up(ii), v1)
                  call gaussian_random_vel(vp(ii), vp(ii), v1)
                  call gaussian_random_vel(wp(ii), wp(ii), v1)
c
c                 Set the particle metal fraction
c
                  if (imetal .eq. 1) then
                     metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
                  else
                     metalf(ii) = 0._RKIND
                  endif
                  if (imetalSNIa .eq. 1) then
                     metalfSNIa(ii) = metalSNIa(i,j,k)
                  endif
c
c                 Remove mass from grid
c
c  original             d(i,j,k) = (1.0_RKIND - starfraction)*d(i,j,k)
                  d(i,j,k) = d(i,j,k) - mp(ii)
                  call sampling_Kroupa_IMF2(temp_starmass)
               enddo

10          continue

            enddo
         enddo
      enddo
 20   continue
c
      if (ii .ge. nmax) then
         write(6,*) 'star_maker4: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii

cc
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback4_individual(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, m_eject, yield,
     &                      tdelay, distrad, diststep, distcells,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Brian O'Shea
c  date:       29 March 2005
c
c    This is a simplified version of stellar feedback, using the method
c    of Andrey Kravtsov as described in "On the origin of the global
c    Schmidt law of star formation", 2003 ApJ 590:L1-L4.  As in this paper,
c    the metals and energy are deposited instantaneously into the gas.
c

c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    r     - refinement field (0 if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, imethod,
     &            distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep
      R_PREC mform, tfactor, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio,
     &       tdelay, temp_starmass, distmass
      parameter (msolar_e51 = 1800.0_RKIND)
c
c-----------------------------------------------------------------------
c  SN delay time = 0 if tdelay < 0
      if (tdelay .lt. 0) tdelay = 0
c
c     Loop over particles
c
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then

c
c         Compute index of the cell that the star particle
c           resides in.
c 
            i = int((xp(n) - xstart)/dx,IKIND) + 1
            j = int((yp(n) - ystart)/dx,IKIND) + 1
            k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               goto 100
            endif

c
c         skip if we're past the first timestep (only feed back
c         on the first timestep)
c

c           if( (t-tcp(n)) .gt. dt ) goto 10   ##### CHANGED TO BELOW FOR 080112_AGORA-DISK-ENZO (5 Myr TIME DELAY INSTEAD OF IMMEDIATE INJECTION)
            if( (t-tcp(n)) .lt. tdelay*(1e6*3.15e7_RKIND/t1) 
     &         .or. (t-tcp(n)) .gt. dt + tdelay*(1e6*3.15e7_RKIND/t1) )
     &         goto 10

c
c          calculate fraction of mass from star particle put back
c          into gas in this cell
c
            temp_starmass = mp(n)
            temp_starmass = temp_starmass / SolarMass
     &           * (dble(d1)*dble(x1*dx)**3)
            if ((temp_starmass .lt. 8.0_RKIND)
     &         .or. (temp_starmass .gt. 100.0_RKIND)) goto 10
c           
c           remnant mass is 2 Msun
            mform = mp(n) - 2.0_RKIND * SolarMass
     &            / (dble(d1)*dble(x1*dx)**3)

c
c          calculate mass added to each cell
c
            distmass = mform / distcells
c
c          if using distributed feedback, check if particle is
c          too close to the boundary
c
            if (distrad .gt. 0) then
               i = max((1 + ibuff + distrad), 
     &              min((nx - ibuff - distrad), i))
               j = max((1 + ibuff + distrad), 
     &              min((ny - ibuff - distrad), j))
               k = max((1 + ibuff + distrad), 
     &              min((nz - ibuff - distrad), k))
            endif
c
c
c           Calculate how much of the star's energy should have
c           gone into supernovae energy:  Do this instantaneously.
c           use the star mass before ejected mass is removed
c           energy = sn_param * mp(n) * (c_light/v1)**2 / 
c           &                 (d(i,j,k)+mform)
c           0.00055957 Msun is 1e51 erg rest mass energy
c           No /(d(i,j,k)+mform) term as in star_maker3, 
c           we will devide those when te update
            energy = 0.00055957_RKIND * SolarMass /
     &            (dble(d1)*dble(x1*dx)**3)
            energy = energy * (c_light/v1)**2/distcells
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae, etc.)
c
            mp(n) = mp(n) - mform
c
c           Add energy to energy field
c
            do kc = k-distrad,k+distrad
               stepk = abs(kc-k)
               do jc = j-distrad,j+distrad
                  stepj = stepk + abs(jc-j)
                  do ic = i-distrad,i+distrad
                     cellstep = stepj + abs(ic-i)
                     if (cellstep .le. diststep) then
                        dratio = 1._RKIND/(d(ic,jc,kc) + distmass)
                        te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                         energy) * dratio
                        if (idual .eq. 1)
     &                      ge(ic,jc,kc) = 
     &                         ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                         energy) * dratio
                        endif
c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
                        if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c
                        yield = 0.3975_RKIND 
     &                        * exp(-17.94_RKIND/temp_starmass) 
     &                        + 57.8094_RKIND 
     &                        * exp(-51.81_RKIND/temp_starmass)
                        yield = yield / temp_starmass  
                        metal(ic,jc,kc) = 
     &                        (metal(ic,jc,kc)*d(ic,jc,kc)
     &                        + distmass * metalf(n) 
     &                        + yield * (mp(n)+mform)/distcells)
     &                        / (d(i,j,k)+distmass)
c
                        endif
c
c           Mass and momentum feedback
c
                        u(ic,jc,kc) = u(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * up(n)
                        v(ic,jc,kc) = v(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * vp(n)
                        w(ic,jc,kc) = w(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * wp(n)
                        d(ic,jc,kc) = d(ic,jc,kc) + distmass
                        u(ic,jc,kc) = u(ic,jc,kc)/d(ic,jc,kc)
                        v(ic,jc,kc) = v(ic,jc,kc)/d(ic,jc,kc)
                        w(ic,jc,kc) = w(ic,jc,kc)/d(ic,jc,kc)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
                        if (imethod .ne. 2 .and. idual .eq. 1) then
                           te(ic,jc,kc) = 0.5_RKIND*(u(ic,jc,kc)**2 + 
     &                          v(ic,jc,kc)**2 + w(ic,jc,kc)**2) +
     &                          ge(ic,jc,kc)
                        endif
                  enddo
               enddo
            enddo
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback3: end'
      return
      end
c=======================================================================
c/////////////////////// SUBROUTINE SAMPLING_KROUPA_IMF \\\\\\\\\\\\\\\\\\
c
      subroutine sampling_Kroupa_IMF2(mass)
c
c     SAMPLES A STELLAR MASS FROM THE KROUPA IMF
c     reference from mcluster by Andreas Kuepper
c
c     written by: Eunwoo Chung
c     date:       25 March 2025
c
c OUTPUTS:
c     mass - sampled stellar mass [in solar masses]
c
c------------------------------------------------------------------------
      implicit none
#include "fortran_types.def"
c------------------------------------------------------------------------
c   
      R_PREC mlow, mup, mass
      R_PREC alpha1, alpha2, c1, c2, k1, k2, xx, mth
c
      parameter (mlow = 0.1_RKIND)
      parameter (mup = 40.0_RKIND)
c
      parameter (alpha1 = 1.3_RKIND)
      parameter (alpha2 = 2.3_RKIND)
c
      c1 = 1.0_RKIND - alpha1
      c2 = 1.0_RKIND - alpha2
c
      k1 = 2.0_RKIND/c1*(0.5_RKIND**c1 - mlow**c1) 
      if (mlow .gt. 0.5_RKIND) then
         k1 = 0.0_RKIND
         k2 = 1.0_RKIND/c2 * (mup**c2 - mlow**c2) 
      else 
         k2 = k1 + 1.0_RKIND/c2 * (mup**c2 - 0.5_RKIND**c2)
      endif
      if (mup .lt. 0.5_RKIND) then
         k1 = 2.0_RKIND/c1 * (mup**c1 - mlow**c1)
         k2 = k1
      endif
c
      call random_number(xx)
      if (xx .lt. k1/k2) then   
         mass = (0.5_RKIND*c1*xx*k2 + mlow**c1) ** (1.0_RKIND/c1)
      else 
         mass=(c2*(xx*k2-k1)+max(0.5_RKIND, mlow)**c2)**(1.0_RKIND/c2)
      endif
c
      return
      end


c=======================================================================
c////////////////// SUBROUTINE GAUSSIAN_RANDOM_POS \\\\\\\\\\\\\\\\\\\\\
c
      subroutine gaussian_random_pos(boxwidth, length_unit, result)
c
c     Generates a random number following a Gaussian distribution
c     using the Box-Muller transform.
c
c INPUTS:
c     boxwidth - Width of the box in which the random number is generated [in code unit]
c     length_unit - code unit to cgs
c
c OUTPUTS:
c     result - Random number following the Gaussian distribution [in code unit]
c
c------------------------------------------------------------------------
      implicit none
#include "fortran_types.def"
c------------------------------------------------------------------------
c
      R_PREC boxwidth, length_unit, result
      R_PREC boxwidth_pc
      R_PREC u1, u2, z0
      R_PREC stddev, pi
c
c     starndard deviation - [in parsec]
      parameter (stddev = 0.1_RKIND)
      parameter (pi = 3.141593_RKIND)
c
      boxwidth_pc = boxwidth * length_unit / pc_cm
c
      do 
c        Generate two uniform random numbers in the range (0, 1)
         call random_number(u1)
         call random_number(u2)
c
c        Apply the Box-Muller transform
         z0 = sqrt(-2.0_RKIND*log(u1))*cos(2.0_RKIND*pi*u2)
c
c        Scale and shift to match the desired mean and standard deviation
         result = stddev * z0
         if (abs(result) .le. boxwidth_pc/2.0_RKIND) exit
      enddo
      result = result * pc_cm / length_unit
c
      return
      end

c=======================================================================
c////////////////// SUBROUTINE GAUSSIAN_RANDOM_POS \\\\\\\\\\\\\\\\\\\\\
c
      subroutine gaussian_random_vel(result, mean, vel_unit)
c
c     Generates a random number following a Gaussian distribution
c     using the Box-Muller transform.
c
c INPUTS:
c     vel_unit - code unit to cgs
c     mean - mean velocity [in km/s]
c
c OUTPUTS:
c     result - Random number following the Gaussian distribution [in code unit]
c
c------------------------------------------------------------------------
      implicit none
#include "fortran_types.def"
c------------------------------------------------------------------------
c
      R_PREC result, mean, vel_unit
      R_PREC u1, u2, z0
      R_PREC stddev_kms, pi, km_cm
c
c     starndard deviation - [in km/s]
      parameter (stddev_kms = 1.0_RKIND/sqrt(3.0_RKIND))
      parameter (pi = 3.141593_RKIND)
      parameter (km_cm = 1.0e5_RKIND)
c
c     Generate two uniform random numbers in the range (0, 1)
      call random_number(u1)
      call random_number(u2)
c
c     Apply the Box-Muller transform
      z0 = sqrt(-2.0_RKIND*log(u1))*cos(2.0_RKIND*pi*u2)
c
c     Scale and shift to match the desired mean and standard deviation
      result = mean + (stddev_kms * z0 * km_cm / vel_unit)
c
      return
      end
