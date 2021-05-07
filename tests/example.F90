!
!     NOVAS SAMPLE CALCULATIONS
!
!     SEE CHAPTER 3 OF USER'S GUIDE FOR EXPLANATION
!     REQUIRES SOLSYS VERSION 2 
!     FOR USE WITH SOLSYS VERSION 1, SEE COMENTS BELOW
!
program example

use novas_module

implicit double precision ( a-h, o-z )
implicit integer ( i-n )

dimension star(6), observ(6), skypos(7), &
     pos(3), vel(3), pose(3), vter(3), vcel(3), &
     sss(3), values(500) 

character names(500)*6, txteph*158, filnam*80

integer denum

parameter ( pi     = 3.14159265358979324d0 )
parameter ( degrad = pi / 180.d0           ) 

data star, observ / 12 * 0.d0 /

data iyear, month, iday, hour, leaps, ut1utc, xp, yp / &
     2008, 4, 24, 10.605d0, 33, -0.387845d0, -0.002d0, +0.529d0 /

data glon, glat, ht / -70.d0, 42.d0, 0.d0 /

data lu, filnam / 20, 'SS_EPHEM.TXT' /


!     FORMAT STATEMENTS FOR OUTPUT

9010 format ( 1x, 'JPL ephemeris DE', i3, ' open. Start JD = ', f10.2, & !
         '  End JD = ', f10.2 ) 
9020 format ( a158 )
9030 format ( 1x, 3 ( f15.10, 8x ) )
9040 format ( 1x, 2 ( f15.6, 8x ), 1f16.11 )
9050 format ( 1x, 2 ( f15.10, 8x ), 1f15.12 )
9060 format ( 1x, 2 ( f16.11, 8x ), 1f15.10 )

#ifdef usejpl
!     CHECK WHICH JPL EPHEMERIS IS AVAILABLE
!     REMOVE THIS BLOCK FOR USE WITH SOLSYS 1
call const ( names, values, sss, n )
denum = values ( 1 )
ephbeg = sss ( 1 )
ephend = sss ( 2 )
write ( *, 9010 )  denum, ephbeg, ephend  
#else
!     VERIFY THAT TEXT EPHEMERIS FILE IS AVAILABLE
!     UNCOMMENT THIS BLOCK FOR USE WITH SOLSYS 1
open ( unit=lu, file=filnam, status='UNKNOWN' )
read ( lu, 9020 ) txteph
write ( *, 9020 ) txteph
close ( lu )
#endif

write ( *, * )
write ( *, * ) 'NOVAS Sample Calculations'
write ( *, * ) '-------------------------'

!     SETUP CALLS
!     HIGH ACCURACY AND EQUINOX MODE ARE NOVAS DEFAULTS.
call hiacc
call eqinox
iearth = idss ( 'EARTH' )

!     WRITE OUT ASSUMED LONGITUDE, LATITUDE, HEIGHT (ITRS = WGS-84)
write ( *, * )
write ( *, * ) 'Geodetic location:'
write ( *, 9030 ) glon, glat, ht

!     ESTABLISH TIME ARGUMENTS
call juldat  ( iyear, month, iday, hour,   utcjd )
ttjd = utcjd + ( leaps + 32.184d0 ) / 86400.d0
ut1jd = utcjd + ut1utc / 86400.d0
deltat = 32.184d0 + leaps - ut1utc
call setdt ( deltat )
write ( *, * )
write ( *, * ) 'TT and UT1 Julian Dates and Delta-T:'
write ( *, 9040 ) ttjd, ut1jd, deltat

!     APPARENT AND TOPOCENTRIC PLACE OF STAR FK6 1307 = GRB 1830
ra2000 = 11.88299133d0
dc2000 = 37.71867646d0
pmra   =  4003.27d0
pmdec  = -5815.07d0
parx   =   109.21d0
rv     =   -98.8d0
call apstar ( ttjd, iearth, ra2000, dc2000, pmra, pmdec, parx, rv, &    !
     ra, dec )
call tpstar ( ut1jd, glon, glat, ht,   rat, dect )
write ( *, * )
write ( *, * ) 'FK6 1307 geocentric and topocentric positions:'
write ( *, 9030 ) ra, dec
write ( *, 9030 ) rat, dect

!     APPARENT AND TOPOCENTRIC PLACE OF THE MOON
moon = idss ( 'MOON' )
call applan ( ttjd, moon, iearth,   ra, dec, dis )
call tpplan ( ut1jd, glon, glat, ht,   rat, dect, dist )
write ( *, * )
write ( *, * ) 'Moon geocentric and topocentric positions:'
write ( *, 9050 ) ra, dec, dis
write ( *, 9050 ) rat, dect, dist

!     TOPOCENTRIC PLACE OF THE MOON -- SHOULD BE SAME AS ABOVE
observ(1) = glon
observ(2) = glat
observ(3) = ht
!      OBSERV(4) = ( TTJD - UT1JD ) * 86400.D0 
call place ( ttjd, 'MOON', 1, 1, star, observ, skypos )
rat  = skypos(4)
dect = skypos(5)
dist = skypos(6)
write ( *, 9050 ) rat, dect, dist

!     POSITION OF THE MOON IN LOCAL HORIZON COORDINATES
!     (POLAR MOTION IGNORED HERE)
call zdaz ( ut1jd, 0.d0, 0.d0, glon, glat, ht, rat, dect, 1, &
              zd, az, rar, decr )
write ( *, * )
write ( *, * ) 'Moon zenith distance and azimuth:'
write ( *, 9030 ) zd, az

!     GREENWICH APPARENT SIDEREAL TIME AND EARTH ROTATION ANGLE
call sidtim ( ut1jd, 0.d0, 1, gast )
st = gast + glon / 15.d0
if ( st .ge. 24.d0 ) st = st - 24.d0
if ( st .lt.  0.d0 ) st = st + 24.d0
call erot ( ut1jd, 0.d0,   era )
write ( *, * )
write ( *, * ) 'Greenwich and local sidereal time and ', &
               'Earth Rotation Angle:'
write ( *, 9060 ) gast, st, era

!     HELIOCENTRIC POSITION OF MARS IN BCRS      
tdbjd = ttjd
!     ABOVE IS AN APPROXIMATION GOOD TO 0.0017 SECONDS 
!     APPROXIMATION COULD LEAD TO ERROR OF ~50 M IN POSITION OF MARS
mars = idss ( 'MARS' )
call solsys ( tdbjd, mars, 1,   pos, vel, ierr )
if ( ierr .eq. 0 ) then
    call eqec ( 0.d0, 0, pos,   pose )
    call angles ( pose,   elon, elat )
    write ( *, * )
    write ( *, * ) 'Mars heliocentric ecliptic longitude and ', &
                   'latitude and radius vector:'
    r = dsqrt (pose(1)**2 + pose(2)**2 +pose(3)**2)
    write ( *, 9050 ) elon * 15.d0, elat, r 
else
    write ( *, * ) 'For Mars, SOLSYS returns error ', ierr
end if

!     TERRESTRIAL TO CELESTIAL TRANSFORMATION
sinlon = dsin ( glon * degrad ) 
coslon = dcos ( glon * degrad ) 
sinlat = dsin ( glat * degrad ) 
coslat = dcos ( glat * degrad ) 
!     FORM VECTOR TOWARD LOCAL ZENITH (ORTHOGONAL TO ELLIPSOID) IN ITRS
vter(1) = coslat * coslon
vter(2) = coslat * sinlon
vter(3) = sinlat
!     TRANSFORM VECTOR TO GCRS
call tercel ( ut1jd, 0.d0, xp, yp, vter,   vcel )
call angles ( vcel,   ra, dec )
write ( *, * ) 
write ( *, * ) 'Direction of zenith vector (RA & Dec) in GCRS:'
write ( *, 9030 ) ra, dec

write ( *, * )
stop
end program example


