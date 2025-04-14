#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback3mom_individual(nx, ny, nz,
     &               d, mu, dm, te, ge, u, v, w,
     &               metal, zfield1, zfield2,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, r, dx, t, z,
     &               dunits, x1, vunits, t1, yield,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tcp, metalf, type,
     &               kinf_in, minitial_arr, mwind_arr, 
     &               msn_arr, teff_arr)

c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Greg Bryan, Christine Simpson
c  date:      July 2015
c  	this is a copy of star_maker3.F that has been modifed to 
c       include the deposition of kinetic as well as thermal energy
c       and metals into a patch surrounding the star particle via a 
c       CIC approach.
c
c  INPUTS:
c
c    nx,ny,nz - dimensions of field arrays
c    d     - density field
c    dm    - dark matter field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    dt    - current timestep
c    r     - refinement field (0 if zone is further refined)
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c
c    npart    - particle array size specified by calling routine
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c    kinf_fixed - fixed kinetic energy fraction for use when var is
c                 false
c    minitial_arr - initial mass of star particle
c    mwind_arr    - ejected mass by wind
c    msn_arr      - ejected mass by supernova
c    teff_arr     - effective temperature of the star
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
#define NO_FORTRAN_DEBUG
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod,
     &      distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    metal(nx,ny,nz)
      R_PREC    mu(nx,ny,nz)
      R_PREC    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      R_PREC    dt, dx, t, z
      R_PREC    dunits, x1, vunits, t1, justburn
      P_PREC    xstart, ystart, zstart
      P_PREC    xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
      R_PREC 	  kinf,kinf_in
      R_PREC    minitial_arr(npart), mwind_arr(npart)
      R_PREC    msn_arr(npart), teff_arr(npart)
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     iface, jface, kface, i, j, k
      R_PREC mform, energy,
     &     yield, minitial,
     &     dist_mass_cells, dist_mom_cells, mass_per_cell, te_per_cell,
     &     thermal_energy_per_cell, kinetic_energy_per_cell,
     &	   energy_per_cell, mom_per_cell, face_shift,
     &     xface, yface, zface, xpos, ypos, zpos,
     &     energy_before, energy_after, mass_before, mass_after,
     &     kin_energy_before, kin_energy_after,
     &     asum, bsum, csum, delta_ke, ke_injected, ke_after
      R_PREC dxf, dyf, dzf, dxc, dyc, dzc, xfc, yfc, zfc,
     &     xfcshift,yfcshift,zfcshift
      R_PREC u1(4,4,4), v1(4,4,4), w1(4,4,4), d1(4,4,4),
     &     te1(4,4,4), ge1(4,4,4), metal1(4,4,4),
     &     ke_before(4,4,4)
      R_PREC Zsol, d_ave, num_d, energy_51, eunits, R_PDS, t_PDS, 
     &	   R_resolve, KE_rem, realmass, mu_cell, fbuff
      R_PREC m1, e1, v_wind_low, v_wind_high
      parameter(v_wind_low = 20.0_RKIND * 1.0e5_RKIND,
     &           v_wind_high = 100.0_RKIND * 1.0e5_RKIND)
      R_PREC mwind, msn, teff, m_eject
c
c-----------------------------------------------------------------------
c
c     Set mass unit and energy unit
      m1 = dunits * x1**3_IKIND
      e1 = m1 * vunits**2_IKIND
c
c     Assuming 3x3x3 cube, calculate cell distribution
c
      distrad = 3_IKIND
      dist_mass_cells = distrad**3_IKIND
c
c     Loop over particles
c
c      write(6,*) 'star_feedback3: start',npart
      do n=1_IKIND, npart
         if (tcp(n) .ge. 0._RKIND .and. (mwind_arr(n) .gt. 0._RKIND .or. 
     &        msn_arr(n) .gt. 0_IKIND)) then
c
            minitial = minitial_arr(n) ! in solar mass
            mwind = mwind_arr(n) ! in solar mass
            msn = msn_arr(n) ! in solar mass
            teff = teff_arr(n) ! in K
c
c     Do instantaneous injeciton
c
            ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND
c
c     
            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
               write(6,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
               goto 100
            endif
c
c	Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
c
c         check bounds - if star particle is near grid edge
c         then shift center of feedback region
c
            if (xfc .lt. xstart+fbuff*dx .or. 
     &          xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &          yfc .lt. ystart+fbuff*dx .or. 
     &          yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &          zfc .lt. zstart+fbuff*dx .or. 
     &          zfc .gt. zstart+dx*nz-fbuff*dx) then
#ifdef FORTRAN_DEBUG
               write(6,*) 'warning1: star feedback zone shifted',
     &              xfc,yfc,zfc, xstart, ystart, zstart,fbuff
#endif
c               
	       xfcshift = xfc
	       yfcshift = yfc 
	       zfcshift	= zfc
c
	       xfc = max(xfc,xstart+fbuff*dx) 
	       yfc = max(yfc,ystart+fbuff*dx) 
	       zfc = max(zfc,zstart+fbuff*dx) 
c              
	       xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
	       yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
	       zfc = min(zfc,zstart+dx*nz-fbuff*dx)
c
	       xfcshift = xfcshift - xfc
	       yfcshift = yfcshift - yfc 
	       zfcshift = zfcshift - zfc
c
#ifdef FORTRAN_DEBUG
	       write(6,*) 'warning2: star feedback zone shifted',
     &              xfc,yfc,zfc,dx,xfcshift,yfcshift,zfcshift, 
     &		    nx, ny,nz
#endif
            endif
c
c           If using zeus, then velocities are face-centered so shift
c
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
c
c         Compute index of the first cell to add momentum,
c             accounting for possible face-centering
c
            xface = (xfc - xstart)/dx - 0.5_RKIND - face_shift
            yface = (yfc - ystart)/dx - 0.5_RKIND - face_shift
            zface = (zfc - zstart)/dx - 0.5_RKIND - face_shift
c
            iface = int(xface + 0.5_RKIND)
            jface = int(yface + 0.5_RKIND)
            kface = int(zface + 0.5_RKIND)
c
            dxf = real(iface) + 0.5_RKIND - xface
            dyf = real(jface) + 0.5_RKIND - yface
            dzf = real(kface) + 0.5_RKIND - zface
c
c         Compute index of the first cell to add mass, assuming cell-centering
c 
            xpos = (xfc - xstart)/dx - 0.5_RKIND
            ypos = (yfc - ystart)/dx - 0.5_RKIND
            zpos = (zfc - zstart)/dx - 0.5_RKIND
c
            ic = int(xpos + 0.5_RKIND)
            jc = int(ypos + 0.5_RKIND)
            kc = int(zpos + 0.5_RKIND)
c
            dxc = real(ic) + 0.5_RKIND - xpos
            dyc = real(jc) + 0.5_RKIND - ypos
            dzc = real(kc) + 0.5_RKIND - zpos
c
c           Stellar wind feedback
            if (mwind .gt. 0._RKIND) then

               m_eject = (mwind * SolarMass / m1 / dx**3_IKIND)
               m_eject = m_eject / (m_eject + mp(n))
               kinf = 0.0_RKIND
c
               if (minitial .lt. 8.0_RKIND) then !AGB stellar wind
c
                  mass_per_cell = mwind * SolarMass / m1 / dx**3_IKIND / 
     &               dist_mass_cells;
                  energy = 1.5_RKIND * 
     &               (mwind * SolarMass / mass_h) *
     &               kboltz * 
     &               teff
                  energy = energy + 0.5_RKIND * (mwind * SolarMass) * 
     &               v_wind_low**2_IKIND
                  energy = energy / e1 / dx**3_IKIND
                  energy_per_cell = energy / dist_mass_cells
c
                  thermal_energy_per_cell = (1._RKIND - kinf) * 
     &                                   energy_per_cell
c
                  do k = 1_IKIND, 4_IKIND
                     do j = 1_IKIND, 4_IKIND
                        do i = 1_IKIND, 4_IKIND
                           u1(i,j,k) = 0._RKIND
                           v1(i,j,k) = 0._RKIND
                           w1(i,j,k) = 0._RKIND
                           d1(i,j,k) = 0._RKIND
                           ge1(i,j,k) = 0._RKIND
                           te1(i,j,k) = 0._RKIND
                           metal1(i,j,k) = 0._RKIND
                        enddo
                     enddo
                  enddo
c
c                 Compute the kinetic energy in the affected region
c                 before momentum is added (except for ZEUS).  This
c                 is needed at the end of the calculation to update the
c                 total energy (te) field.
c
                  if (imethod .ne. 2_IKIND) then
                     do k = -1_IKIND, +2_IKIND
                        do j = -1_IKIND, +2_IKIND
                           do i = -1_IKIND, +2_IKIND
                              ke_before(i+2, j+2, k+2) = 
     &                             0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                  (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                   v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                 w(ic+i ,jc+j ,kc+k)**2_IKIND)                     
                           enddo
                        enddo
                     enddo
                  endif
c
                  call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                          nx, ny, nz, ic, jc, kc,
     &                          iface, jface, kface, imethod, imetal,
     &                          +1_IKIND)
c
                  call sum_mass_kinetic_energy2(u, v, w, d, ge, 
     &                 te, nx, ny, nz,
     &                 iface, jface, kface, ic, jc, kc,
     &                 mass_before, kin_energy_before, idual, imethod)
c
                  call add_feedback2(u1, v1, w1, d1, ge1, te1, metal1, 
     &                              4_IKIND, 4_IKIND, 4_IKIND,
     &                              2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND,
     &                                 2_IKIND, 2_IKIND,
     &                              dxf, dyf, dzf, dxc, dyc, dzc,
     &                              imethod, imetal, imulti_metals, 
     &                              idual, m_eject, yield, metalf(n),
     &                              mass_per_cell, 1._RKIND, 0._RKIND)
c
                  mom_per_cell = 0._RKIND
c
                  call add_feedback2(u, v, w, d, ge, te, 
     &                              metal, nx, ny, nz,
     &                              ic, jc, kc, iface, jface, kface,
     &                              dxf, dyf, dzf, dxc, dyc, dzc,
     &                              imethod, imetal, 
     &                              imulti_metals, idual,
     &                              m_eject, yield, metalf(n),
     &                              mass_per_cell, mom_per_cell,
     &                              thermal_energy_per_cell)
c
                  call sum_mass_kinetic_energy2(u, v, w, d, ge, te, 
     &                 nx, ny, nz,
     &                 iface, jface, kface, ic, jc, kc,
     &                 mass_after, kin_energy_after, idual, imethod)
c
                  call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                          nx, ny, nz, ic, jc, kc,
     &                          iface, jface, kface, imethod, imetal,
     &                          -1_IKIND)
c
                  if (imethod .ne. 2_IKIND) then
                     ke_injected = 0._RKIND
                     do k = -1_IKIND, +2_IKIND
                        do j = -1_IKIND, +2_IKIND
                           do i = -1_IKIND, +2_IKIND
                              ke_after = 0.5_RKIND*
     &                                   d(ic+i, jc+j, kc+k)*
     &                                  (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                   v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                   w(ic+i ,jc+j ,kc+k)**2_IKIND)
                              delta_ke = ke_after-ke_before(i+2,j+2,k+2)
                              te(ic+i,jc+j,kc+k) =
     &                                  te(ic+i,jc+j,kc+k) +
     &                                  delta_ke/d(ic+i,jc+j,kc+k)
                              ke_injected = ke_injected + delta_ke
                           enddo
                        enddo
                     enddo
                  endif !end AGB case
c
               else !Massive stellar wind
c
                  mass_per_cell = mwind * SolarMass / m1 / 
     &               dx**3_IKIND / 
     &               dist_mass_cells;
                  energy = 1.5_RKIND * 
     &               (mwind * SolarMass / mass_h) *
     &               kboltz * 
     &               teff
                  energy = energy + 0.5_RKIND * 
     &               (mwind * SolarMass) * 
     &               v_wind_high**2_IKIND
                  energy = energy / e1 / dx**3_IKIND
                  energy_per_cell = energy / 
     &               dist_mass_cells
c
                  thermal_energy_per_cell = 
     &                (1._RKIND - kinf) *
     &                energy_per_cell
c
                  do k = 1_IKIND, 4_IKIND
                     do j = 1_IKIND, 4_IKIND
                        do i = 1_IKIND, 4_IKIND
                           u1(i,j,k) = 0._RKIND
                           v1(i,j,k) = 0._RKIND
                           w1(i,j,k) = 0._RKIND
                           d1(i,j,k) = 0._RKIND
                           ge1(i,j,k) = 0._RKIND
                           te1(i,j,k) = 0._RKIND
                           metal1(i,j,k) = 0._RKIND
                        enddo
                     enddo
                  enddo
c
c                 Compute the kinetic energy in the affected region
c                 before momentum is added (except for ZEUS).  This
c                 is needed at the end of the calculation to update the
c                 total energy (te) field.
c
                  if (imethod .ne. 2_IKIND) then
                     do k = -1_IKIND, +2_IKIND
                        do j = -1_IKIND, +2_IKIND
                           do i = -1_IKIND, +2_IKIND
                              ke_before(i+2, j+2, k+2) =
     &                             0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                  (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                   v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                 w(ic+i ,jc+j ,kc+k)**2_IKIND)
                           enddo
                        enddo
                     enddo
                  endif
c
                  call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                          nx, ny, nz, ic, jc, kc,
     &                          iface, jface, kface, imethod, imetal,
     &                          +1_IKIND)
c
                  call sum_mass_kinetic_energy2(u, v, w, d, ge, te, 
     &                 nx, ny, nz,
     &                 iface, jface, kface, ic, jc, kc,
     &                 mass_before, kin_energy_before, idual, imethod)
c
                  call add_feedback2(u1, v1, w1, d1, ge1, te1, metal1, 
     &                              4_IKIND, 4_IKIND, 4_IKIND,
     &                              2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND,
     &                              2_IKIND, 2_IKIND,
     &                              dxf, dyf, dzf, dxc, dyc, dzc,
     &                              imethod, imetal, imulti_metals, 
     &                              idual, m_eject, yield, metalf(n),
     &                              mass_per_cell, 1._RKIND, 0._RKIND)
c
                  mom_per_cell = 0._RKIND
c
                  call add_feedback2(u, v, w, d, ge, te, metal, 
     &                              nx, ny, nz,
     &                              ic, jc, kc, iface, jface, kface,
     &                              dxf, dyf, dzf, dxc, dyc, dzc,
     &                              imethod, imetal, imulti_metals, 
     &                              idual, m_eject, yield, metalf(n),
     &                              mass_per_cell, mom_per_cell,
     &                              thermal_energy_per_cell)
c
                  call sum_mass_kinetic_energy2(u, v, w, d, ge, te, 
     &                 nx, ny, nz,
     &                 iface, jface, kface, ic, jc, kc,
     &                 mass_after, kin_energy_after, idual, imethod)
c
                  call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                          nx, ny, nz, ic, jc, kc,
     &                          iface, jface, kface, imethod, imetal,
     &                          -1_IKIND)
c
                  if (imethod .ne. 2_IKIND) then
                     ke_injected = 0._RKIND
                     do k = -1_IKIND, +2_IKIND
                        do j = -1_IKIND, +2_IKIND
                           do i = -1_IKIND, +2_IKIND
                              ke_after = 0.5_RKIND*
     &                                   d(ic+i, jc+j, kc+k)*
     &                                  (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                   v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                   w(ic+i ,jc+j ,kc+k)**2_IKIND)
                              delta_ke = ke_after - 
     &                                  ke_before(i+2,j+2,k+2)
                              te(ic+i,jc+j,kc+k) =
     &                                  te(ic+i,jc+j,kc+k) +
     &                                  delta_ke/d(ic+i,jc+j,kc+k)
                              ke_injected = ke_injected + delta_ke
                           enddo
                        enddo
                     enddo
                  endif
               endif !end massive case
c
               if (abs(kin_energy_after - 
     &                 kin_energy_before - kinf*energy)
     &              /(kinf*energy) > 0.01_RKIND .and. 
     &              kinf .ne. 0._RKIND) then
c
                  write(6,*) 'Kinetic energy added to mesh',
     &                 ' does not match',
     &                 ' kinetic energy feedback parameters'
                  write(6,*) 'kinf:',kinf
                  write(6,*) 'kin_energy_after - kin_energy_before',
     &                 kin_energy_after - kin_energy_before
                  write(6,*) 'energy*fkin',kinf*energy
                  write(6,*) '%diff',abs(kin_energy_after - 
     &                 kin_energy_before)
     &                 /(kin_energy_before)
                  ERROR_MESSAGE
               endif

               if (abs(mass_after - mass_before 
     &                 - mass_per_cell*dist_mass_cells)
     &              /(mass_per_cell*dist_mass_cells) > 0.01_RKIND) then
c
                  write(6,*) 'Mass added to mesh does not match', 
     &                 ' kinetic energy feedback parameters'
                  write(6,*) 'kinf:',kinf
                  write(6,*) 'mass_after, mass_before:',
     &                        mass_after, mass_before
                  write(6,*) 'kin_energy_after - kin_energy_before',
     &                        kin_energy_after - kin_energy_before
                  write(6,*) 'energy*fkin',kinf*energy
                  write(6,*) '%diff',abs(mass_after - 
     &                 mass_before)
     &                 /(mass_before)
                  ERROR_MESSAGE
               endif
c
            endif !end wind feedback
c
c           Supernova feedback
            if (msn .gt. 0._RKIND) then
c
               m_eject = (msn * SolarMass / m1 / dx**3_IKIND)
               m_eject = m_eject / (m_eject + mp(n))
c
               mass_per_cell = msn * SolarMass / m1 / dx**3_IKIND / 
     &            dist_mass_cells
c            
               energy = 1e51_RKIND / e1 / dx**3_IKIND
               energy_per_cell = energy / dist_mass_cells
c
c              Use fixed kinf value, unless kinf < 0 - then compute
c              variable kinf
               kinf = kinf_in
               if (kinf .lt. 0._RKIND) then
c              Calculate variable kinf based on R_PDS from
c              Cioffi et al. 1988
c
                  energy_51 = 1.0_RKIND
c
                  Zsol = 0._RKIND
                  num_d = 0._RKIND
                  d_ave = 0._RKIND
c     
c                 Compute average density and metallcity around cell containing
c                 particle
c     
c                write(6,*) 'averaging over patch'
                  do k = -1_IKIND,1_IKIND
                     do j = -1_IKIND,1_IKIND
                        do i = -1_IKIND,1_IKIND
                           mu_cell = mu(ic+i,jc+j,kc+k)
                           Zsol = Zsol + 
     &                        metal(ic+i,jc+j,kc+k)/0.02_RKIND
                           num_d = num_d + 
     &                        d(ic+i,jc+j,kc+k)*dunits/mu_cell/mass_h
                           d_ave = d_ave + d(ic+i,jc+j,kc+k)*dunits
                        enddo
                     enddo
                  enddo
c               
                  Zsol = Zsol/27._RKIND
                  num_d = num_d/27._RKIND
                  d_ave = d_ave/27._RKIND
c   
c                 compute time and radius of transition to PDS phase for gas
c                 with the computed properties
c                 Note: t_PDS has units of 1e3yrs and R_PDS has units of pc
c          
c                 For metal poor gas
c     
                  if (Zsol .lt. 0.01_RKIND) then
                     t_PDS = 3.06e2_RKIND * 
     &                  (energy_51**(1._RKIND/8._RKIND)) *
     &                  (num_d**(-3._RKIND/4._RKIND))
                     R_PDS = 49.3_RKIND * (energy_51**0.25_RKIND) * 
     &                  (num_d**(-0.5_RKIND))
                  endif
c
c                 For metal rich gas
c
                  if (Zsol .ge. 0.01_RKIND) then
                     t_PDS = 26.5_RKIND * 
     &                  (energy_51**(3._RKIND/14._RKIND)) * 
     &                  (Zsol**(-5._RKIND/14._RKIND)) * 
     &                  (num_d**(-4._RKIND/7._RKIND))
                     R_PDS = 18.5_RKIND * energy_51**(2._RKIND/7._RKIND) 
     &                     * num_d**(-3._RKIND/7._RKIND)
     &                     * Zsol**(-1._RKIND/7._RKIND)
                  endif
                  R_resolve = dx*x1/pc_cm
                  if (R_PDS .gt. 4.5_RKIND*R_resolve) then
                     kinf = 0._RKIND
#ifdef FORTRAN_DEBUG
                     write(6,*) 'Sedov phase resolved:',Zsol,num_d,
     &                  energy_51,R_PDS,R_resolve, realmass
#endif
                  else                  
                     kinf = 3.97133e-6_RKIND * (d_ave/mass_h) 
     &                    * (R_resolve**(-2._RKIND)) 
     &                    * (R_PDS**7._RKIND) * (t_PDS**(-2._RKIND)) 
     &                    * (energy_51**(-1._RKIND))

#ifdef FORTRAN_DEBUG
                     write(6,*) 'Sedov phase not resolved:',kinf,Zsol,
     &                    num_d, d_ave, 
     &                    energy_51,t_PDS,R_PDS,R_resolve
#endif
                  
                  endif
               endif
c
c     Check kinf.  If it is too small, set it to 0
c
               if (kinf .le. 1e-10_RKIND) then 
#ifdef FORTRAN_DEBUG
                  write(6,*) 'kinf = ',kinf
                  write(6,*) 'kinf too small!  Setting to 0'
#endif
                  kinf = 0._RKIND
               endif
c
c              Finished computing kinf; 
c              now compute the amount of thermal energy
c              that needs to be added to each cell
               thermal_energy_per_cell = (1._RKIND - kinf) * 
     &                                   energy_per_cell

c              write(6,*) 'thermal_energy_per_cell:',
c     &           thermal_energy_per_cell,kinf,energy_per_cell
c
c              Inject energy, mass & metals
c
c              Zero local dummy field and kinetic energy field
c
               do k = 1_IKIND, 4_IKIND
                  do j = 1_IKIND, 4_IKIND
                     do i = 1_IKIND, 4_IKIND
                        u1(i,j,k) = 0._RKIND
                        v1(i,j,k) = 0._RKIND
                        w1(i,j,k) = 0._RKIND
                        d1(i,j,k) = 0._RKIND
                        ge1(i,j,k) = 0._RKIND
                        te1(i,j,k) = 0._RKIND
                        metal1(i,j,k) = 0._RKIND
                     enddo
                  enddo
               enddo
c
c              Compute the kinetic energy in the affected region
c              before momentum is added (except for ZEUS).  This
c              is needed at the end of the calculation to update the
c              total energy (te) field.
c
               if (imethod .ne. 2_IKIND) then
                  do k = -1_IKIND, +2_IKIND
                     do j = -1_IKIND, +2_IKIND
                        do i = -1_IKIND, +2_IKIND
                           ke_before(i+2, j+2, k+2) = 
     &                          0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                               (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                w(ic+i ,jc+j ,kc+k)**2_IKIND)                     
                        enddo
                     enddo
                  enddo
               endif
c
c              First convert velocities to momenta and transform
c              into frame comoving with particle
c
               call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                       nx, ny, nz, ic, jc, kc, 
     &                       iface, jface, kface, imethod, imetal, 
     &                       +1_IKIND)
c
c              Sum mass and energy before
c
               call sum_mass_kinetic_energy2(u, v, w, d, ge,
     &              te, nx, ny, nz,
     &              iface, jface, kface, ic, jc, kc,
     &              mass_before, kin_energy_before, idual, imethod)
c               write(6,*) "m,e before:", mass_before, kinetic_energy_before
c
c              Now add mass and momentum terms (normalization 1.0) to
c                 local dummy fields
c
               call add_feedback2(u1, v1, w1, d1, ge1, te1, metal1, 
     &                           4_IKIND, 4_IKIND, 4_IKIND, 
     &                           2_IKIND, 2_IKIND, 2_IKIND, 2_IKIND, 
     &                           2_IKIND, 2_IKIND,
     &                           dxf, dyf, dzf, dxc, dyc, dzc,
     &                           imethod, imetal, imulti_metals, idual,
     &                           m_eject, yield, metalf(n),
     &                           mass_per_cell, 1._RKIND, 0._RKIND)
c
c              The kinetic energy after the feedback event is a quadratic equation
c              with delta_p, where delta_p is the total amount of momentum added
c              to the feedback region:
c              E_k,a = asum + bsum*delta_p + csum*delta_p*delta_p
c
c              Sum a, b and c terms to get momentum normalization
c
               if (kinf .gt. 0._RKIND) then
                  call sum_abc2(u, v, w, d, ge, u1, v1, w1, d1, 
     &                 nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &                 asum, bsum, csum)
#ifdef FORTRAN_DEBUG
                  write(6,*) 'abc1',asum, bsum, csum, kinf
#endif
                  asum = asum - (kin_energy_before + kinf*energy)
#ifdef FORTRAN_DEBUG
c                  write(6,*) 'kin_energy_before,kinf*energy,energy:',
c     &                 kin_energy_before,kinf*energy,energy
                  write(6,*) 'abc2',asum, bsum, csum
#endif
c
c              Calculate momentum contribution
c
                  mom_per_cell = (-bsum + 
     &                 sqrt(bsum**2_IKIND - 4._RKIND*asum*csum))/
     &                 (2._RKIND*csum)
               else
c
c              If the kinetic fraction is set to 0, then 0 momentum is added
c              to each cell.  Note that this is different from adding 0 kinetic
c              energy to 
c
                  mom_per_cell = 0._RKIND
               endif
c           
#ifdef FORTRAN_DEBUG
               write(6,*) 'mass_per_cell:',mass_per_cell
               write(6,*) 'mom_per_cell:',mom_per_cell
               write(6,*) 'thermal_energy_per_cell:',
     &              thermal_energy_per_cell
#endif
c
c              Now add mass and momentum, using three-point CIC
c              Note that if the te field is being used, it is 
c              only updated with the thermal energy
c
               call add_feedback2(u, v, w, d, ge, te, metal, nx, ny, nz, 
     &                           ic, jc, kc, iface, jface, kface, 
     &                           dxf, dyf, dzf, dxc, dyc, dzc,
     &                           imethod, imetal, imulti_metals, idual,
     &                           m_eject, yield, metalf(n),
     &                           mass_per_cell, mom_per_cell, 
     &			         thermal_energy_per_cell)
c
c              Sum mass and energy before
c
               call sum_mass_kinetic_energy2(u, v, w, d, ge, te, 
     &              nx, ny, nz,
     &              iface, jface, kface, ic, jc, kc,
     &              mass_after, kin_energy_after, idual, imethod)
c
#ifdef FORTRAN_DEBUG
               write(6,*) 'energy:',energy
               write(6,*) 'mass (b,a,preda):', mass_before, mass_after, 
     &              mass_before + mform*m_eject,mform*m_eject
               write(6,*) 'kinetic energy (b,a,preda):',kin_energy_before,
     &              kin_energy_after, kin_energy_before + kinf*energy, 
     &              kinf
#endif
c
c              Error checks to see if the kinetic energy added in the particle frame
c              matches the values set by the input parameters
c
               if (abs(kin_energy_after - 
     &                 kin_energy_before - kinf*energy)
     &              /(kinf*energy) > 0.01_RKIND .and. 
     &              kinf .ne. 0._RKIND) then
c
                  write(6,*) 'Kinetic energy added to mesh',
     &                 ' does not match',
     &                 ' kinetic energy feedback parameters'
                  write(6,*) 'kinf:',kinf
                  write(6,*) 'kin_energy_after - kin_energy_before',
     &                 kin_energy_after - kin_energy_before
                  write(6,*) 'energy*fkin',kinf*energy
                  write(6,*) '%diff',abs(kin_energy_after - 
     &                 kin_energy_before)
     &                 /(kin_energy_before)
                  ERROR_MESSAGE
               endif

               if (abs(mass_after - mass_before 
     &                 - mass_per_cell*dist_mass_cells)
     &              /(mass_per_cell*dist_mass_cells) > 0.01_RKIND) then
c
                  write(6,*) 'Mass added to mesh does not match', 
     &                 ' kinetic energy feedback parameters'
                  write(6,*) 'kinf:',kinf
                  write(6,*) 'mass_after, mass_before:',
     &                        mass_after, mass_before
                  write(6,*) 'kin_energy_after - kin_energy_before',
     &                        kin_energy_after - kin_energy_before
                  write(6,*) 'energy*fkin',kinf*energy
                  write(6,*) '%diff',abs(mass_after - 
     &                 mass_before)
     &                 /(mass_before)
                  ERROR_MESSAGE
               endif
c
c              Convert momenta back to velocities and transform back to lab frame
c
               call momentum2(u, v, w, d, metal, up(n), vp(n), wp(n),
     &                       nx, ny, nz, ic, jc, kc, 
     &                       iface, jface, kface, imethod, imetal,
     &                       -1_IKIND)
c
c              Add the increase in the kinetic energy to the total energy field
c                (unless we're using Zeus).  If using dual energy formalism, we
c                might want to enforce consistency.
c     
               if (imethod .ne. 2_IKIND) then
                  ke_injected = 0._RKIND
                  do k = -1_IKIND, +2_IKIND
                     do j = -1_IKIND, +2_IKIND
                        do i = -1_IKIND, +2_IKIND
                           ke_after = 0.5_RKIND*
     &                           d(ic+i, jc+j, kc+k) *
     &                          (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                           v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                           w(ic+i ,jc+j ,kc+k)**2_IKIND) 
                           delta_ke = ke_after - ke_before(i+2,j+2,k+2)
                           te(ic+i ,jc+j ,kc+k) = 
     &                          te(ic+i ,jc+j ,kc+k) + 
     &                          delta_ke/d(ic+i, jc+j, kc+k)
                           ke_injected = ke_injected + delta_ke
                        enddo
                     enddo
                  enddo
               endif
            endif !end supervanova feedback
c
 10         continue
         endif !for proper particles only
c
 100     continue
c
      enddo !loop over particles
c
c      write(6,*) 'star_feedback3: end'
      return
      end
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine momentum2(u, v, w, d, metal, up, vp, wp, 
     &                    nx, ny, nz, ic, jc, kc, 
     &                    iface, jface, kface, imethod, imetal, idir)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC iface, jface, kface, imethod, imetal, idir
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp
c
c     Locals
c
      INTG_PREC i, j, k
c
c     Error check
c
      if (idir .ne. -1_IKIND .and. idir .ne. 1_IKIND) then
         write(6,*) 'incorrect idir value in momentum call'
         ERROR_MESSAGE
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
c
c              idir = +1: convert vel -> mom
c               
               if (idir. eq. +1_IKIND) then
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = 
     &                    (u(iface+i ,jc+j ,kc+k) - up)
     &                    * 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                    d(iface+i+1_IKIND, jc+j, kc+k))
                     v(ic+i ,jface+j ,kc+k) = 
     &                    (v(ic+i ,jface+j ,kc+k) - vp)
     &                    * 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                    d(ic+i, jface+j+1_IKIND, kc+k))
                     w(ic+i ,jc+j ,kface+k) = 
     &                    (w(ic+i ,jc+j ,kface+k) - wp)
     &                    * 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                    d(ic+i, jc+j, kface+k+1_IKIND))
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                                      d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                                      d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                                      d(ic+i, jc+j, kc+k)
                  endif
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) *
     &                    d(ic+i,jc+j,kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = u(iface+i ,jc+j ,kc+k)
     &                    /( 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                             d(iface+i+1_IKIND, jc+j, kc+k))) + up
                     v(ic+i ,jface+j ,kc+k) = v(ic+i ,jface+j ,kc+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                             d(ic+i, jface+j+1_IKIND, kc+k))) + vp
                     w(ic+i ,jc+j ,kface+k) = w(ic+i ,jc+j ,kface+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                             d(ic+i, jc+j, kface+k+1_IKIND))) + wp
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + wp
                  endif
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = metal(ic+i,jc+j,kc+k) /
     &                    d(ic+i,jc+j,kc+k)
                  endif
               endif
c
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_mass_kinetic_energy2(pu, pv, pw, d, ge, te, 
     &                           nx, ny, nz,
     &                           iface, jface, kface, ic, jc, kc,
     &                           mass_sum, kin_energy_sum,
     &                           idual, imethod)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, iface, jface, kface, ic, jc, kc, idual
      INTG_PREC imethod
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz), ge(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_sum, kin_energy_sum
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    mass_term, mom_term, kin_energy
c
      mass_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc   +k)
               mom_term = pu(iface+i ,jc   +j, kc   +k)**2_IKIND + 
     &                    pv(ic   +i ,jface+j, kc   +k)**2_IKIND + 
     &                    pw(ic   +i ,jc   +j, kface+k)**2_IKIND
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)
               mass_sum = mass_sum + mass_term
               kin_energy_sum = kin_energy_sum + kin_energy
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine sum_abc2(pu, pv, pw, d, ge, pu1, pv1, pw1, d1, 
     &                   nx, ny, nz, iface, jface, kface, ic, jc, kc,
     &                   asum, bsum, csum)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, iface, jface, kface, ic, jc, kc
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu1(4,4,4), pv1(4,4,4), pw1(4,4,4), d1(4,4,4)
      R_PREC    asum, bsum, csum
c
c     Locals
c
      INTG_PREC i, j, k, istart, jstart, kstart
      R_PREC    mass_term, mom_term, aterm, bterm, cterm
      R_PREC    mass_sum, energy_sum
c
      asum = 0._RKIND
      bsum = 0._RKIND
      csum = 0._RKIND
      mass_sum = 0._RKIND
      energy_sum = 0._RKIND
c
      istart = 2_IKIND
      jstart = 2_IKIND
      kstart = 2_IKIND
c
c     Loop over all affected cells
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
               mass_term = d(ic   +i ,jc   +j ,kc+k)
               mom_term = pu(iface+i ,jc   +j, kc+k   )**2_IKIND + 
     &                    pv(ic   +i ,jface+j, kc+k   )**2_IKIND +
     &                    pw(ic   +i ,jc   +j, kface+k)**2_IKIND

               mass_term = mass_term + d1(istart+i,jstart+j,kstart+k)
               asum = asum + mom_term/(2._RKIND * mass_term)
               bterm = 
     &           pu(iface+i,jc+j,kc+k)*pu1(istart+i,jstart+j,kstart+k) +
     &           pv(ic+i,jface+j,kc+k)*pv1(istart+i,jstart+j,kstart+k) +
     &           pw(ic+i,jc+j,kface+k)*pw1(istart+i,jstart+j,kstart+k)
               bsum = bsum + bterm/mass_term ! factors of 2 cancel
               cterm = 
     &           pu1(istart+i,jstart+j,kstart+k)**2_IKIND +
     &           pv1(istart+i,jstart+j,kstart+k)**2_IKIND +
     &           pw1(istart+i,jstart+j,kstart+k)**2_IKIND
c               
               csum = csum + cterm/(2._RKIND*mass_term)
            enddo
         enddo
      enddo
c
c      write(6,*) "m,e in sum:", mass_sum, energy_sum
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback2(pu, pv, pw, d, ge, te, metal, nx, ny, nz, 
     &                        ic, jc, kc, iface, jface, kface,
     &                        dxf, dyf, dzf, dxc, dyc, dzc,
     &                        imethod, imetal, imulti_metals, idual,
     &                        m_eject, yield, metalf,
     &                        mass_per_cell, mom_per_cell, 
     &                        therm_per_cell)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual
      INTG_PREC iface, jface, kface, imethod, imetal, imulti_metals
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell, therm_per_cell
      R_PREC    m_eject, yield, metalf
      R_PREC    dxf, dyf, dzf, dxc, dyc, dzc,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw, delta_therm,
     &     dratio, tot_mass
c
c     Error check
c
      if (imulti_metals .eq. 1_IKIND) then
         write(6,*) "momentum: not supported"
         ERROR_MESSAGE
      endif
c
c     Loop "cells" in particle-frame
c
      tot_mass = 0._RKIND
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
c     For each particle "cell", do CIC-like deposit
c     compute zone and face centered weight factors
c
               do i1 = i, i+1_IKIND
                  dxf1 = dxf
                  dxc1 = dxc
                  if (i1 .eq. i+1_IKIND) dxf1 = 1._RKIND - dxf
                  if (i1 .eq. i+1_IKIND) dxc1 = 1._RKIND - dxc
                  do j1 = j, j+1_IKIND
                     dyf1 = dyf
                     dyc1 = dyc
                     if (j1 .eq. j+1_IKIND) dyf1 = 1._RKIND - dyf
                     if (j1 .eq. j+1_IKIND) dyc1 = 1._RKIND - dyc
                     do k1 = k, k+1_IKIND
                        dzf1 = dzf
                        dzc1 = dzc
                        if (k1 .eq. k+1_IKIND) dzf1 = 1._RKIND - dzf
                        if (k1 .eq. k+1_IKIND) dzc1 = 1._RKIND - dzc

                        delta_mass =   mass_per_cell*dxc1*dyc1*dzc1
                        delta_pu   =  i*mom_per_cell*dxf1*dyc1*dzc1
                        delta_pv   =  j*mom_per_cell*dxc1*dyf1*dzc1
                        delta_pw   =  k*mom_per_cell*dxc1*dyc1*dzf1
                        delta_therm = therm_per_cell*dxc1*dyc1*dzc1
c
c     Add mass, momentum
c     (add thermal energy here)
c
                        dratio = d(ic+i1, jc+j1, kc+k1)/
     &                       (d(ic+i1, jc+j1, kc+k1) + delta_mass)
                        d(ic+i1 ,jc+j1 ,kc+k1) = d(ic+i1 ,jc+j1 ,kc+k1)
     &                       + delta_mass
                        pu(iface+i1 ,jc+j1 ,kc+k1) = 
     &                       pu(iface+i1 ,jc+j1 ,kc+k1)+ delta_pu
                        pv(ic+i1 ,jface+j1 ,kc+k1) = 
     &                       pv(ic+i1 ,jface+j1 ,kc+k1) + delta_pv
                        pw(ic+i1 ,jc+j1 ,kface+k1) = 
     &                       pw(ic+i1 ,jc+j1 ,kface+k1)+ delta_pw
c                        
                        tot_mass = tot_mass + delta_mass
c
c     Add thermal energy
c     
                        te(ic+i1 ,jc+j1 ,kc+k1) = 
     &                       te(ic+i1 ,jc+j1 ,kc+k1)*dratio + 
     &                       delta_therm / d(ic+i1, jc+j1, kc+k1)
c                        
                        if (idual .eq. 1_IKIND)
     &                       ge(ic+i1 ,jc+j1 ,kc+k1) = 
     &                       ge(ic+i1, jc+j1, kc+k1)*dratio +
     &                       delta_therm / d(ic+i1, jc+j1, kc+k1)
                        
c     
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas) rather than a density.
c     The conversion has been done in the handling routine)
c     
                        if (imetal .eq. 1_IKIND) then
c     
c     "Cen method".  This takes into account gas recycling.
c     (metal is not a fraction here)
c     
                           metal(ic+i1, jc+j1, kc+k1) = 
     &                          metal(ic+i1, jc+j1, kc+k1) + 
     &                          (delta_mass/m_eject)
     &                          * (yield * (1._RKIND-metalf) + 
     &                          m_eject * metalf)
                           
                        endif
c     
c     End loop over CIC-deposit
c     
                     enddo
                  enddo
               enddo
c     
c     End loop over "cells" in particle-frame
c     
            enddo
         enddo
      enddo
c
c     write(6,*) 'tot_mass:',tot_mass
c     
      return
      end
