!*******************************************************************
!>
!  This program checks out many parts of NOVAS Fortran, producing
!  an output file checkout.out that can be compared with validation
!  files in the NOVAS distribution.

program checkout

use novas_module

implicit none

integer nstars, ntimes, i, j, k, isun, iearth, idss

double precision tjd, starda, deltat, glat, glon, ht, &
 ujd, rai, deci, pmra, pmdec, parlax, radvel, &
 rag, decg, disg, rat, dect, dist, racio, ujd1, ujd2, gast 

character starid*20, method*8

dimension starid(4), starda(6,4), tjd(4), deltat(4)

data nstars, ntimes / 4, 4 /

data starid / 'Polaris    HIP 11767', &
              'Delta Ori  HIP 25930', &
              'Theta Car  HIP 52419', &
              'Grb 1830   FK6  1307'/

data starda / &
 37.94614689d0, 89.26413805d0, 44.22d0, -11.74d0,  7.56d0,-17.4d0, &    !
 83.00166562d0, -0.29909340d0,  1.67d0,   0.56d0,  3.56d0, 16.0d0, &    !
160.73927802d0,-64.39447937d0,-18.87d0,  12.06d0,  7.43d0, 24.0d0, &    !
11.88299133d0,37.71867646d0,4003.27d0,-5815.07d0,109.21d0,-98.8d0/      !

data tjd / 2453737.5d0, 2455013.7d0, 2456232.2d0, 2457208.9d0 /

data deltat / 64.8452d0, 66.2d0,  68.d0,  69.d0 / 

data glat /45.0d0/, glon /-75.0d0/, ht /0.0d0/

!-----Get body numbers for Sun and Earth and open output file.

isun = idss ( "SUN" )
iearth = idss ( "EARTH" )

open ( unit=8, file="checkout.out", status="unknown" )

!-----Compute geocentric places of stars.

do i = 1, nstars

   if ( index (starid(i),"HIP") /= 0 ) then
       call gethip ( starda(1,i), starda(2,i), starda(3,i), &
                     starda(4,i), starda(5,i), starda(6,i), &
                     rai, deci, pmra, pmdec, parlax, radvel )
   else
       rai = starda(1,i)
       deci = starda(2,i)
       pmra = starda(3,i)
       pmdec = starda(4,i)
       parlax = starda(5,i)
       radvel = starda(6,i)
   end if

   write (8, 101) starid(i),rai,deci,pmra,pmdec,parlax,radvel

!-----Loop through the times.

   do j = 1, ntimes
      
      call apstar (tjd(j),iearth, &
        rai,deci,pmra,pmdec,parlax,radvel, rag,decg)
   
      write (8, 201) tjd(j), rag, decg
      
   end do

end do

!------Compute topocentric Sun.

   write (8, 102) "Topocentric Sun:"

   do j = 1, ntimes
         
      ujd = tjd(j) - (deltat(j) / 86400.0d0)
      
      call applan (tjd(j),isun,iearth, rag,decg,disg) 
      call tpplan (ujd,glon,glat,ht, rat,dect,dist)
   
      write (8, 201) tjd(j), rat, dect, dist
    
   end do

!------Compute Grenwich apparent sidereal time.

   call ciotio
   call cioloc ( tjd(1),  racio, k )
   if ( k == 1 ) then
     method = "external"
   else
     method = "internal" 
   end if 

   write (8, 102) "Sidereal Time (using "//method//" CIO):" 

   do j = 1, ntimes
   
      call setdt ( deltat(j) )
         
      ujd1 = tjd(j)
      ujd2 = 0.d0

      call sidtim (ujd1,ujd2,1,gast)
      
      write (8, 202) ujd1, gast
    
   end do

end file ( unit=8 )
close ( unit=8 )
   
stop

101 format (//,1x, "Star: ", a20,": ", / &
  f12.9, 1x, f13.9, 1x, f8.2, 1x, f8.2, 1x, f6.2, 1x, f6.1)
102 format (//,1x, a, / )
201 format (1x, "JD = ", f9.1, 2x, "RA = ", f12.9, 2x, "Dec = ", &
 f12.8, :, 2x, "Dis = ", f10.8)
202 format (1x, "JD = ", f9.1, 2x, "ST = ", f12.9 )   

!*******************************************************************
  end program checkout
!*******************************************************************