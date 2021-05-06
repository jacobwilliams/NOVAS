      PROGRAM checkout
      
      IMPLICIT NONE
      
C-----This program checks out many parts of NOVAS Fortran, producing
C     an output file checkout.out that can be compared with validation
C     files in the NOVAS distribution.

      INTEGER nstars, ntimes, i, j, k, isun, iearth, idss
      
      DOUBLE PRECISION tjd, starda, deltat, glat, glon, ht,
     . ujd, rai, deci, pmra, pmdec, parlax, radvel,
     . rag, decg, disg, rat, dect, dist, racio, ujd1, ujd2, gast 
     
      CHARACTER starid*20, method*8

      DIMENSION starid(4), starda(6,4), tjd(4), deltat(4)

      data nstars, ntimes / 4, 4 /
  
      DATA starid / 'Polaris    HIP 11767',
     .              'Delta Ori  HIP 25930',
     .              'Theta Car  HIP 52419',
     .              'Grb 1830   FK6  1307'/
      
      DATA starda /
     . 37.94614689d0, 89.26413805d0, 44.22d0, -11.74d0,  7.56d0,-17.4d0,
     . 83.00166562d0, -0.29909340d0,  1.67d0,   0.56d0,  3.56d0, 16.0d0,
     .160.73927802d0,-64.39447937d0,-18.87d0,  12.06d0,  7.43d0, 24.0d0,
     .11.88299133d0,37.71867646d0,4003.27d0,-5815.07d0,109.21d0,-98.8d0/
      
      DATA tjd / 2453737.5d0, 2455013.7d0, 2456232.2d0, 2457208.9d0 /
     
      DATA deltat / 64.8452d0, 66.2d0,  68.d0,  69.d0 / 
     
      DATA glat /45.0d0/, glon /-75.0d0/, ht /0.0d0/

C-----Get body numbers for Sun and Earth and open output file.

      isun = idss ( "SUN" )
      iearth = idss ( "EARTH" )
      
      OPEN ( unit=8, file="checkout.out", status="unknown" )
      
     
C-----Compute geocentric places of stars.

      DO 10 i = 1, nstars

         IF ( index (starid(i),"HIP") .NE. 0 ) THEN
             CALL gethip ( starda(1,i), starda(2,i), starda(3,i),
     .                     starda(4,i), starda(5,i), starda(6,i),
     .                     rai, deci, pmra, pmdec, parlax, radvel )
         ELSE
             rai = starda(1,i)
             deci = starda(2,i)
             pmra = starda(3,i)
             pmdec = starda(4,i)
             parlax = starda(5,i)
             radvel = starda(6,i)
         END IF
 
         WRITE (8, 101) starid(i),rai,deci,pmra,pmdec,parlax,radvel

C-----Loop through the times.

         DO 20 j = 1, ntimes
            
            CALL apstar (tjd(j),iearth,
     .        rai,deci,pmra,pmdec,parlax,radvel, rag,decg)
         
            WRITE (8, 201) tjd(j), rag, decg
            
   20    CONTINUE
   
   10 CONTINUE

C------Compute topocentric Sun.

         WRITE (8, 102) "Topocentric Sun:"

         DO 30 j = 1, ntimes
               
            ujd = tjd(j) - (deltat(j) / 86400.0d0)
            
            CALL applan (tjd(j),isun,iearth, rag,decg,disg) 
            CALL tpplan (ujd,glon,glat,ht, rat,dect,dist)
         
            WRITE (8, 201) tjd(j), rat, dect, dist
          
   30 CONTINUE
 
C------Compute Grenwich apparent sidereal time.

         CALL ciotio
         CALL cioloc ( tjd(1),  racio, k )
         IF ( k .eq. 1 ) THEN
           method = "external"
         ELSE
           method = "internal" 
         END IF 

         WRITE (8, 102) "Sidereal Time (using "//method//" CIO):" 

         DO 40 j = 1, ntimes
         
            call setdt ( deltat(j) )
               
            ujd1 = tjd(j)
            ujd2 = 0.d0

            CALL sidtim (ujd1,ujd2,1,gast)
            
            WRITE (8, 202) ujd1, gast
          
   40 CONTINUE
   
      END FILE ( unit=8 )
      CLOSE ( unit=8 )
         
      STOP
   
  101 FORMAT (//,1X, "Star: ", A20,": ", / 
     .  F12.9, 1X, F13.9, 1X, F8.2, 1X, F8.2, 1X, F6.2, 1X, F6.1)
  102 FORMAT (//,1X, A, / )
  201 FORMAT (1X, "JD = ", F9.1, 2X, "RA = ", F12.9, 2X, "Dec = ", 
     . F12.8, :, 2X, "Dis = ", F10.8)
  202 FORMAT (1X, "JD = ", F9.1, 2X, "ST = ", F12.9 )   
   
      END
