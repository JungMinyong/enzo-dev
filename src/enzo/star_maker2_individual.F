#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker2_individual(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, tindsf, mintdyn,
     &                      odthresh, masseff, smthresh, usejeans,
     &                      level, np, 
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      imetalSNIa, metalSNIa, metalfSNIa)

c
c  CREATES GALAXY PARTICLES
c
c  written by: Chris Loken
c  date:       3 March 1997
c  modified1: 2 March 1999 by Brian O''Shea
c    2 inputs were added: odthresh and masseff, and the star creation
c    code was implemented
c  modified2: 18 May 1999 by Brian O''Shea
c    1 input added: smthresh, and star particle selection code added
c  modified3: 26 July 2002 by BWO
c    this version of star_maker2.src is hereby certified to be free of
c    the various bugs that Greg Bryan and I discovered in April 2002
c  modified4: 13 November 2002 by BWO
c    turned OFF stochastic star formation (changed ifdef) and cleaned
c    up code, added comments.
c  modified5: 20 Feb 2008 by M. Ryan Joung
c    added metalSNIa & metalfSNIa; included feedback from SN Ia/PN
c
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    tindsf   - if true, remove dt/t_dyn term in starfraction
c    odthresh - overdensity threshold (some number * avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    usejeans - use jeans mass check 
c    level - current level of refinement
c    procnum - processor number (for output)
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
c    metalfSNIa - metallicity fraction of particle (from SN Ia)
c    nmax     - particle array size specified by calling routine
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC imetalSNIa, procnum, tindsf, usejeans
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax),
     &     metalfSNIa(nmax)
      R_PREC    odthresh, masseff, smthresh, mintdyn
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
      R_PREC   x_rand, y_rand, z_rand
      R_PREC   possible_mass, temp_starmass
      parameter (sndspdC=1.3095e8_RKIND)
c
      ii = np
c
c  for each zone, : "star" particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c    is the gas mass greater than the Jeans mass?
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
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) is divergence negative?
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
               if (imethod .eq. 2) then
                  div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
               else
                  div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
               endif
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi_val/32._RKIND/
     &                      GravConst/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              5) is M > M_Jeans? (this definition involves only baryons under
c                 the assumption that the dark matter is stable, which
c                 implies that the dark matter velocity dispersion is >> 
c                 the sound speed.  This will be true for small perturbations
c                 within large halos).
c
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / SolarMass
               isosndsp2 = sndspdC * temp(i,j,k)
               if (usejeans .eq. 1) then
                  jeanmass = pi_val/(6._RKIND*
     &                 sqrt(d(i,j,k)*dble(d1)))*
     &                 dble(pi_val * isosndsp2 /
     &                 GravConst)**1.5_RKIND / SolarMass

                  if (bmass .lt. jeanmass) goto 10
               endif
c
c              6) Check to see if star is above threshold (given
c                 in units of M_solar)
c

               if (tindsf .eq. 1) then
                  starfraction = masseff
                  tdyn = mintdyn*3.15d7/t1
               else
                  starfraction = min(masseff*dt/tdyn, 0.9_RKIND)
                  tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)
               endif

c
c  13 Nov. 2002:  stochastic star formation has been turned OFF by
c    bwo, since the algorithm is a bit suspect.  This encourages
c    rational settings of the minimum star mass!
c
#define NO_STOCHASTIC_STAR_FORMATION
c
#ifdef STOCHASTIC_STAR_FORMATION
c
c                 Keep global count of "unfullfilled" star formation
c                 and when total is larger than threshold, then create
c                 a star particle with the threshold mass or 1/2 the
c                 gas in the cell, whichever is smaller.
c
               if (starfraction*bmass .lt. smthresh) then
                  sformsum = sformsum + starfraction*bmass
                  if (sformsum .lt. smthresh) goto 10
                  starfraction = min(smthresh/bmass, 0.5_RKIND)
                  sformsum = sformsum - starfraction*bmass
               endif
#else
c
c              is star mass greater than threshold, then make it.
c              if it's less than threshold, go to the next cell.
c
               if (starfraction*bmass .lt. smthresh) goto 10
#endif
c
c              Create a star particle
c
               possible_mass = starfraction * bmass
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
c
c
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)
c
 1000          format(i5,1x,6(1pe10.3,1x))
c
10          continue

            enddo
         enddo
      enddo
 20   continue
c
      if (ii .ge. nmax) then
         write(6,*) 'star_maker2: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
cc
      return
      end
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback2_individual(nx, ny, nz,
     &                      d, dm, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod, dt, r, dx, t, z,
     &                      d1, x1, v1, t1, sn_param, m_eject, yield,
     &                      distrad, diststep, distcells,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, type, justburn, 
     &                      crmodel, crfeedback, cr)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Chris Loken & Greg Bryan
c  date:       3 March 1997
c  modified1:  BWO
c              13 Nov 2002
c              Many changes have been made between the date of
c              initial creation and today - all of them unlogged,
c              unfortunately.  Bugs were fixed in calculation of 
c              total energy and conservation of metal and gas density.
c              The code has been cleaned up in general to enhance
c              readability.  This is the stable version - star_maker1
c              and star_maker3 should be used for experimentation.
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
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c    crmodel   - boolean: 1 = include cosmic ray physics, 0 = no cosmic ray physics                                                
c    crfeedback - fraction of thermal supernova energy to return as cosmic rays                                                    
c    cr         - Cosmic ray BaryonField  
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
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
     &      distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
      INTG_PREC crmodel
      R_PREC cr(nx,ny,nz)
      R_PREC crfeedback
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep
      R_PREC mform, tfactor, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio,
     &     distmass
      parameter (msolar_e51 = 1800._RKIND)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback2: start'
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c


c        The star particle creation algorithm partnered with this 
c          feedback algorithm creates a star particle instantaneously.
c          However, we do feedback as if the star particles are created 
c          over a long period of time (in code units), so the particle
c          actually loses mass over time in an exponentially decaying
c          way.

c        Determine how much of a given star particle would have been 
c          turned into stars during this timestep.  Then calculate the mass
c          which should have formed during this timestel dt using the integral
c          form of the Cen & Ostriker formula.

            xv1 = (t      - tcp(n))/tdp(n)
            if (xv1 .gt. 12._RKIND) goto 10     ! t-tcp >> tdp so ignore
            xv2 = (t + dt - tcp(n))/tdp(n)

c        First calculate the initial mass of the star particle 
c          in question.
            minitial = mp(n) / 
     &           (1._RKIND - m_eject*(1._RKIND - (1._RKIND + xv1)
     &           *exp(-xv1)))
c
c       Then, calculate the amount of mass that would have formed in
c         this timestep.
c
            mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &                          (1._RKIND + xv2)*exp(-xv2))
            mform = max(min(mform, mp(n)), 0._RKIND)
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
c          skip if very little mass is formed.
c
            if (mform/d(i,j,k) .lt. 1.e-10_RKIND) goto 10
c
c          calculate mass added to each cell
c
            distmass = (mform * m_eject) / distcells
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
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c
c           Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3
c
c           Calculate how much of the star formation in this
c           timestep would have gone into supernova energy.
c
            energy = sn_param * mform * (c_light/v1)**2/distcells
c
c             (aug 7/00 addition: share ejected energy between gas
c              being ejected and that remaining in star particle)
c
#define NO_SHARE_ENERGY
#ifdef SHARE_ENERGY
            energy = energy*mform*m_eject/
     &           (mform*m_eject + minitial*exp(-xv2)*(1._RKIND + xv2))
#endif /* SHARE_ENERGY */
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
                        if( crmodel .gt. 0 ) then
                          te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                          energy*(1._RKIND-crfeedback)) * dratio
                          if (idual .eq. 1)
     &                      ge(ic,jc,kc) =
     &                         ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                         energy*(1._RKIND-crfeedback)) * dratio
                          cr(ic,jc,kc) = cr(ic,jc,kc)+energy*crfeedback
                        else
                           te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                          energy) * dratio
                           if (idual .eq. 1)
     &                          ge(ic,jc,kc) =
     &                          ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                          energy) * dratio
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
                           metal(ic,jc,kc) = 
     &                          (metal(ic,jc,kc)*d(ic,jc,kc) + 
     &                          (mform / distcells) * 
     &                          (yield * (1._RKIND-metalf(n)) + 
     &                          m_eject * metalf(n))) * dratio
            !metal is a fraction
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
c      write(6,*) 'star_feedback2: end'
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
      parameter (mlow = 1.0_RKIND)
      parameter (mup = 150.0_RKIND)
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