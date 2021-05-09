!***********************************************************************
!>
!  NOVAS : NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE
!
!### Author
!  * G. H. KAPLAN, U.S. NAVAL OBSERVATORY
!
!### Version
!  * NOVAS FORTRAN VERS F3.1 of 2011 MARCH 21

    module novas_module

    !implicit none

    contains
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE APPARENT DIRECTION OF A STAR OR SOLAR
!  SYSTEM BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE
!  SYSTEM.  BASED ON KAPLAN, ET AL. (1989), ASTRONOMICAL JOURNAL 97,
!  1197-1210, WITH SOME ENHANCEMENTS FROM KLIONER (2003),
!  ASTRONOMICAL JOURNAL 125, 1580-1597.
!
!       TJD    = TT JULIAN DATE FOR PLACE (IN)
!       OBJECT = CHARACTER STRING IDENTIFYING OBJECT OF INTEREST (IN)
!                FOR SOLAR SYSTEM
!                BODY,             SPECIFY THE NAME USING ALL UPPER-
!                                  CASE LETTERS ('SUN', 'MOON',
!                                  'JUPITER', ETC.),
!                                  - OR -
!                                  SPECIFY THE BODY ID NUMBER
!                                  IN A 4-CHARACTER STRING OF THE
!                                  FORM '=NNN', WHERE NNN IS THE
!                                  BODY ID NUMBER
!                FOR STAR,         PROVIDE A BLANK STRING, THE WORD
!                                  'STAR', OR ANY STRING BEGINNING
!                                  WITH '*'
!       LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
!                SET LOCATN=0 FOR OBSERVER AT GEOCENTER
!                SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
!                SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
!       ICOORD = INTEGER CODE SPECIFYING COORDINATE SYSTEM OF OUTPUT
!                POSITION (IN)
!                SET ICOORD=0 FOR GCRS (OR 'LOCAL GCRS')
!                SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
!                SET ICOORD=3 FOR ASTROMETRIC COORDINATES, I.E.,
!                             WITHOUT LIGHT DEFLECTION OR ABERRATION
!       STAR   = ARRAY OF CATALOG DATA FOR STAR (IN)
!                (NOT USED IF SOLAR SYSTEM BODY REQUESTED)
!                STAR(1) = ICRS RIGHT ASCENSION IN HOURS
!                STAR(2) = ICRS DECLINATION IN DEGREES
!                STAR(3) = ICRS PROPER MOTION IN RA IN
!                          MILLIARCSECONDS/YEAR
!                STAR(4) = ICRS PROPER MOTION IN DEC IN
!                          MILLIARCSECONDS/YEAR
!                STAR(5) = PARALLAX IN MILLIARCSECONDS
!                STAR(6) = RADIAL VELOCITY IN KILOMETERS/SECOND
!                FURTHER STAR ARRAY ELEMENTS ARE NOT USED HERE
!                BUT ARE RESERVED FOR FUTURE USE
!       OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
!                (NOT USED IF LOCATN=0)
!                FOR LOCATN=1,
!                OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
!                            (EAST +) IN DEGREES
!                OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
!                            (NORTH +) IN DEGREES
!                OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
!                            IN METERS
!                OBSERV(4) = VALUE OF DELTA-T IN SECONDS
!                            (DELTA-T=TT-UT1)
!                OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
!                OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
!                FOR LOCATN=2,
!                OBSERV(1) = GEOCENTRIC X IN KILOMETERS
!                OBSERV(2) = GEOCENTRIC Y IN KILOMETERS
!                OBSERV(3) = GEOCENTRIC Z IN KILOMETERS
!                OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND
!                OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND
!                OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
!                WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
!       SKYPOS = ARRAY OF OUTPUT DATA SPECIFYING OBJECT'S PLACE
!                ON THE SKY AT TIME TJD, WITH RESPECT TO THE
!                SPECIFIED OUTPUT COORDINATE SYSTEM (OUT)
!                SKYPOS(1) = X, DIMENSIONLESS      UNIT VECTOR
!                SKYPOS(2) = Y, DIMENSIONLESS      TOWARD OBJECT
!                SKYPOS(3) = Z, DIMENSIONLESS
!                SKYPOS(4) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
!                            RIGHT ASCENSION IN HOURS
!                SKYPOS(5) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
!                            DECLINATION IN DEGREES
!                SKYPOS(6) = TRUE (GEOMETRIC, EUCLIDIAN) DISTANCE
!                            TO SOLAR SYSTEM BODY IN AU AT TIME TJD,
!                            OR 0.D0 FOR STAR
!                SKYPOS(7) = RADIAL VELOCITY IN KILOMETERS/SECOND
!                FURTHER SKYPOS ARRAY ELEMENTS ARE NOT USED HERE
!                BUT ARE RESERVED FOR FUTURE USE
!
!  NOTE 1: VALUES OF LOCATN AND ICOORD FOR VARIOUS STANDARD KINDS
!  OF PLACE:
!  LOCATN=0 AND ICOORD=1 APPARENT PLACE
!  LOCATN=1 AND ICOORD=1 TOPOCENTRIC PLACE
!  LOCATN=0 AND ICOORD=0 VIRTUAL PLACE
!  LOCATN=1 AND ICOORD=0 LOCAL PLACE
!  LOCATN=0 AND ICOORD=3 ASTROMETRIC PLACE
!  LOCATN=1 AND ICOORD=3 TOPOCENTRIC ASTROMETRIC PLACE
!
!  NOTE 2: ARRAYS STAR AND SKYPOS MAY BE EXPANDED IN THE FUTURE, AND
!  THIS CAN BE ALLOWED FOR IN THE CALLING CODE BY DIMENSIONING
!  THESE ARRAYS WITH 20 AND 10 ELEMENTS, RESPECTIVELY, EVEN THOUGH
!  ELEMENTS BEYOND STAR(6) AND SKYPOS(7) ARE NOT NOW REFERRED TO IN
!  THIS SUBROUTINE.
!
!  NOTE 3: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
!  BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
!  DEFINED BY THE USER VIA CALL TO SETDT.
!
!  NOTE 4: SKYPOS(7), THE RADIAL VELOCITY, IS THE PREDICTED
!  RADIAL VELOCITY MEASURE (Z) TIMES THE SPEED OF LIGHT, AN
!  INHERENTLY SPECTROSCOPIC MEASURE.  FOR A STAR, IT
!  INCLUDES ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT,
!  CONTAINED IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE,
!  WHICH IS ASSUMED GIVEN IN STAR(6).  FOR A SOLAR SYSTEM
!  BODY, IT APPLIES TO A FICTITIOUS EMITTER AT THE CENTER OF THE
!  OBSERVED OBJECT, ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT),
!  AND DOES NOT IN GENERAL APPLY TO REFLECTED LIGHT.

    subroutine place ( tjd, object, locatn, icoord, star, observ, &
        skypos )

! --- INITIAL DECLARATIONS---------------------------------------------

implicit none
integer locatn,icoord,ntimes,iearth,isun,idbody,ierr,loc,j,kcio, &
     idss
double precision tjd,star,observ,skypos, &
     t0,tlast1,tlast2,ttjd,tdbjd,c,x,secdif,tlight,dis,dt, &
     frlimb,rcio,peb,veb,psb,vsb,pog,vog,pob,vob, &
     pos1,vel1,pos2,pos3,pos4,pos5,pos6,pos7,pos8, &
     px,py,pz,ra,dec,rvs,rvd,rv,dabs,dsqrt
character*(*) object

dimension star(*), observ(6), skypos(*), &
     peb(3), veb(3), psb(3), vsb(3), &
     pog(3), vog(3), pob(3), vob(3), &
     pos1(3), vel1(3), pos2(3), pos3(3), pos4(3), pos5(3), &
     pos6(3), pos7(3), pos8(3), &
     px(3), py(3), pz(3), rvs(3), rvd(3)

save

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast1, tlast2 / 0.d0, 0.d0 /,   ntimes / 0 /

3 format ( ' PLACE: CANNOT OBTAIN COORDINATES OF ', a, ' AT JD ', &
     f10.1 )
4 format ( ' PLACE: WILL NOT PROCESS EARTH AS OBSERVED OBJECT ', &
     'EXCEPT WHEN LOCATN=2' )

! --- GET CONSTANTS, FIRST TIME ONLY ----------------------------------

ntimes = ntimes + 1

if ( ntimes == 1 ) then
    iearth = idss ( 'EARTH' )
    isun = idss ( 'SUN' )
!         GET C, THE SPEED OF LIGHT IN AU/DAY
    call astcon ( 'C(AU/DAY)', 1.d0,   c )
end if

! --- CHECK ON EARTH AS AN OBSERVED OBJECT ----------------------------

if ( object == 'EARTH' .and. locatn /= 2 ) then
    write ( *, 4 )
    return
end if

! --- GET POSITION AND VELOCITY OF EARTH (GEOCENTER) AND SUN ----------

if ( dabs ( tjd - tlast1 ) > 1.d-8 ) then

!         COMPUTE TDBJD, THE TDB JULIAN DATE CORRESPONDING TO TTJD
    ttjd = tjd
    tdbjd = tjd
    call times ( tdbjd, x,   secdif )
    tdbjd = ttjd + secdif / 86400.d0

!         GET POSITION AND VELOCITY OF THE EARTH WRT BARYCENTER OF
!         SOLAR SYSTEM, IN ICRS
    call solsys ( tdbjd, iearth, 0,   peb, veb, ierr )
    if ( ierr /= 0 ) then
        write ( *, 3 ) 'EARTH', tjd
        return
    end if

!         GET POSITION AND VELOCITY OF THE SUN WRT BARYCENTER OF
!         SOLAR SYSTEM, IN ICRS
    call solsys ( tdbjd, isun, 0,   psb, vsb, ierr )
    if ( ierr /= 0 ) then
        write ( *, 3 ) 'SUN', tjd
        return
    end if

    tlast1 = tjd

end if

! --- GET POSITION AND VELOCITY OF OBSERVER ---------------------------

if ( locatn == 1 .or. locatn == 2 ) then

!         FOR TOPOCENTRIC PLACE, GET GEOCENTRIC POSITION AND VELOCITY
!         VECTORS OF OBSERVER
    call geopos ( ttjd, locatn, observ,   pog, vog )
    loc = 1

else

!         FOR GEOCENTRIC PLACE, THERE IS NOTHING TO DO
    do j = 1, 3
       pog(j) = 0.d0
       vog(j) = 0.d0
    end do
    loc = 0

end if

!     COMPUTE POSITION AND VELOCITY OF OBSERVER WRT BARYCENTER OF
!     SOLAR SYSTEM (GALILEAN TRANSFORMATION FROM GCRS to BCRS)
do j=1,3
    pob(j) = peb(j) + pog(j)
    vob(j) = veb(j) + vog(j)
end do

! --- FIND GEOMETRIC POSITION OF OBSERVED OBJECT ----------------------

if ( object == 'STAR' .or. object == ' ' .or. &
     object(1:1) == '*' ) then

!         OBSERVED OBJECT IS STAR
    idbody = -9999

!         GET POSITION OF STAR UPDATED FOR ITS SPACE MOTION
    call vectrs ( star(1), star(2), star(3), star(4), &
                  star(5), star(6),   pos1, vel1 )
    call dlight ( pos1, pob,   dt )
    call propmo ( t0, pos1, vel1, tdbjd + dt,   pos2 )
!         GET POSITION OF STAR WRT OBSERVER (CORRECTED FOR PARALLAX)
    call geocen ( pos2, pob,   pos3, tlight )
    dis = 0.d0

else

!         OBSERVED OBJECT IS SOLAR SYSTEM BODY

!         GET ID NUMBER OF BODY
    if ( object(1:1) == '=' ) then
        read ( object, '(1X,I3)' ) idbody
    else
        idbody = idss ( object )
        if ( idbody == -9999 ) then
            write ( *, 3 ) object, tjd
            return
        end if
    end if

!         GET POSITION OF BODY WRT BARYCENTER OF SOLAR SYSTEM
    call solsys ( tdbjd, idbody, 0,   pos1, vel1, ierr )
    if ( ierr /= 0 ) then
        write ( *, 3 ) object, tjd
        return
    end if

!         GET POSITION OF BODY WRT OBSERVER, AND TRUE (EUCLIDIAN)
!         DISTANCE
    call geocen ( pos1, pob,   pos2, tlight )
    dis = tlight * c

!         GET POSITION OF BODY WRT OBSERVER, ANTEDATED FOR LIGHT-TIME
    call littim ( tdbjd, idbody, pob, tlight,   pos3, tlight )

end if

! --- APPLY GRAVITATIONAL DEFLECTION OF LIGHT AND ABERRATION ----------

if ( icoord == 3 ) then

!         THESE CALCULATIONS ARE SKIPPED FOR ASTROMETRIC PLACE
    do j = 1, 3
        pos5(j) = pos3(j)
    end do

else

!         VARIABLE LOC DETERMINES WHETHER EARTH DEFLECTION IS INCLUDED
    if ( loc == 1 ) then
        call limang ( pos3, pog,   x, frlimb )
        if ( frlimb < 0.8d0 ) loc = 0
    end if

!         COMPUTE GRAVITATIONAL DEFLECTION AND ABERRATION
    call grvdef ( tdbjd, loc, pos3, pob,   pos4 )
    call aberat ( pos4, vob, tlight,   pos5 )
!         POSITION VECTOR IS NOW IN GCRS

end if

! --- TRANSFORM, IF NECESSARY, TO OUTPUT COORDINATE SYSTEM ------------

if ( icoord == 1 ) then

!         TRANSFORM TO EQUATOR AND EQUINOX OF DATE
    call frame  ( pos5, 1,   pos6 )
    call preces ( t0, pos6, tdbjd,   pos7 )
    call nutate ( tdbjd, pos7,   pos8 )

else if ( icoord == 2 ) then

!         TRANSFORM TO EQUATOR AND CIO OF DATE
    if ( dabs ( tdbjd - tlast2 ) > 1.d-8 ) then
!             OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
!             INTERMEDIATE SYSTEM
        call cioloc ( tdbjd,   rcio, kcio )
        call ciobas ( tdbjd, rcio, kcio,   px, py, pz )
        tlast2 = tdbjd
    end if
!         TRANSFORM POSITION VECTOR TO CELESTIAL INTERMEDIATE SYSTEM
    pos8(1) = px(1) * pos5(1) + px(2) * pos5(2) + px(3) * pos5(3)
    pos8(2) = py(1) * pos5(1) + py(2) * pos5(2) + py(3) * pos5(3)
    pos8(3) = pz(1) * pos5(1) + pz(2) * pos5(2) + pz(3) * pos5(3)

else

!         NO TRANSFORMATION -- KEEP COORDINATES IN GCRS
!         (OR ICRS FOR ASTROMETRIC PLACE)
    do j = 1, 3
        pos8(j) = pos5(j)
    end do

end if

! --- GET RADIAL VELOCITY ---------------------------------------------

!     SET UP STAR DATA, IF APPLICABLE
if ( idbody == -9999 ) then
    rvs(1) = star(1)
    rvs(2) = star(2)
    rvs(3) = star(6)
    if ( star(5) <= 0.d0 ) then
        vel1(1) = 0.d0
        vel1(2) = 0.d0
        vel1(3) = 0.d0
    end if
else
    rvs(1) = 0.d0
    rvs(2) = 0.d0
    rvs(3) = 0.d0
end if

!     COMPUTE DISTANCES: OBSERVER-GEOCENTER, OBSERVER-SUN, OBJECT-SUN
rvd(1) = dsqrt ( ( pob(1)  - peb(1) ) ** 2 &
               + ( pob(2)  - peb(2) ) ** 2 &
               + ( pob(3)  - peb(3) ) ** 2 )
rvd(2) = dsqrt ( ( pob(1)  - psb(1) ) ** 2 &
               + ( pob(2)  - psb(2) ) ** 2 &
               + ( pob(3)  - psb(3) ) ** 2 )
rvd(3) = dsqrt ( ( pos1(1) - psb(1) ) ** 2 &
               + ( pos1(2) - psb(2) ) ** 2 &
               + ( pos1(3) - psb(3) ) ** 2 )

call radvl ( pos3, vel1, vob, rvs, rvd,   rv )

! --- FINISH UP -------------------------------------------------------

call angles ( pos8,   ra, dec )

x = dsqrt ( pos8(1)**2 + pos8(2)**2 + pos8(3)**2 )

do j = 1, 3
    skypos(j) = pos8(j) / x
end do

skypos(4) = ra
skypos(5) = dec
skypos(6) = dis
skypos(7) = rv

call setvec ( pos8 )

end

!***********************************************************************
!>
!  THE ENTRIES TO THIS SUBROUTINE PROVIDE 'FRONT ENDS' TO
!  SUBROUTINE PLACE, TAILORED TO SPECIFIC PLACE TYPES.  THEY
!  PROVIDE COMPATIBILITY WITH PREVIOUSLY SUPPORTED CALLING
!  SEQUENCES.

subroutine places()

implicit none
integer l,n,locatn,icoord
double precision tjd,rai,deci,pmra,pmdec,parlax,radvel, &
     ujd,glon,glat,ht,ra,dec,dis, &
     ttjd,gast,deltat,star,observ,skypos,dmod
character*4 object

dimension star(20), observ(6), skypos(10)

save

data star, observ, skypos / 36 * 0.d0 /


! --- GEOCENTRIC PLACES OF STARS --------------------------------------
!
!     ARGUMENTS COMMON TO THESE ENTRIES:
!
!          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
!          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
!                   (NO LONGER USED)
!          RAI    = ICRS RIGHT ASCENSION IN HOURS (IN)
!          DECI   = ICRS DECLINATION IN DEGREES (IN)
!          PMRA   = ICRS PROPER MOTION IN RA IN MILLIARCSECONDS/YEAR
!                   (IN)
!          PMDEC  = ICRS PROPER MOTION IN DEC IN MILLIARCSECONDS/YEAR
!                   (IN)
!          PARLAX = PARALLAX IN MILLIARCSECONDS (IN)
!          RADVEL = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
!          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
!          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
!
!     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
!     VPSTAR, ICRS FOR ASSTAR, AND EQUATOR AND EQUINOX OF DATE FOR
!     APSTAR.
!
!
!     THIS ENTRY PROVIDES THE APPARENT PLACE OF A STAR.
entry apstar (tjd,n,rai,deci,pmra,pmdec,parlax,radvel, ra,dec)
locatn = 0
icoord = 1
go to 10

!     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A STAR.
entry vpstar (tjd,n,rai,deci,pmra,pmdec,parlax,radvel, ra,dec)
locatn = 0
icoord = 0
go to 10

!     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A STAR.
entry asstar (tjd,n,rai,deci,pmra,pmdec,parlax,radvel, ra,dec)
locatn = 0
icoord = 3

10 ttjd = tjd
object = 'STAR'
star(1) = rai
star(2) = deci
star(3) = pmra
star(4) = pmdec
star(5) = parlax
star(6) = radvel

call place (ttjd,object,locatn,icoord,star,observ, skypos)

ra  = skypos(4)
dec = skypos(5)

return


! --- TOPOCENTRIC PLACES OF STARS -------------------------------------
!
!     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
!     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
!     ARGUMENTS COMMON TO THESE ENTRIES:
!
!          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
!          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
!                   IN DEGREES (IN)
!          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
!                   IN DEGREES (IN)
!          HT     = HEIGHT OF OBSERVER IN METERS (IN)
!          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
!          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
!
!     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
!     FOR LPSTAR, AND EQUATOR AND EQUINOX OF DATE FOR TPSTAR.
!
!     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
!     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
!     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
!
!
!     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A STAR.
entry tpstar (ujd,glon,glat,ht, ra,dec)
locatn = 1
icoord = 1
go to 20

!     THIS ENTRY PROVIDES THE LOCAL PLACE OF A STAR.
entry lpstar (ujd,glon,glat,ht, ra,dec)
locatn = 1
icoord = 0

20 if ( ujd > 100.d0 ) then
    deltat = ( ttjd - ujd ) * 86400.d0
else
    gast = dmod ( ujd, 24.d0 )
    if ( gast < 0.d0 ) gast = gast + 24.d0
    call placst ( gast )
    deltat = 0.d0
end if
observ(1) = glon
observ(2) = glat
observ(3) = ht
observ(4) = deltat

call place (ttjd,object,locatn,icoord,star,observ, skypos)

ra  = skypos(4)
dec = skypos(5)

return


! --- GEOCENTRIC PLACES OF PLANETS ------------------------------------
!
!     ARGUMENTS COMMON TO THESE ENTRIES:
!
!          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
!          L      = BODY IDENTIFICATION NUMBER FOR DESIRED PLANET (IN)
!          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
!                   (NO LONGER USED)
!          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
!          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
!          DIS    = TRUE DISTANCE FROM EARTH TO PLANET IN AU (OUT)
!
!     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
!     VPPLAN, ICRS FOR ASPLAN, AND EQUATOR AND EQUINOX OF DATE FOR
!     APPLAN.
!
!     NOTE: 'PLANET' IS USED GENERICALLY FOR ANY SOLAR SYSTEM BODY.
!
!     THIS ENTRY PROVIDES THE APPARENT PLACE OF A PLANET.
entry applan (tjd,l,n, ra,dec,dis)
locatn = 0
icoord = 1
go to 30

!     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A PLANET.
entry vpplan (tjd,l,n, ra,dec,dis)
locatn = 0
icoord = 0
go to 30

!     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A PLANET.
entry asplan (tjd,l,n, ra,dec,dis)
locatn = 0
icoord = 3

30 ttjd = tjd
if (l == -9999) then
    object = '    '  ! JW : bug in original?
else
    write ( object, '(''='',I3)') l
end if

call place (ttjd,object,locatn,icoord,star,observ, skypos)

ra  = skypos(4)
dec = skypos(5)
dis = skypos(6)

return


! --- TOPOCENTRIC PLACES OF PLANETS -----------------------------------
!
!     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
!     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
!     ARGUMENTS COMMON TO THESE ENTRIES:
!
!          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
!          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
!                   IN DEGREES (IN)
!          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
!                   IN DEGREES (IN)
!          HT     = HEIGHT OF OBSERVER IN METERS (IN)
!          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
!          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
!          DIS    = TRUE DISTANCE FROM OBSERVER TO PLANET IN AU (OUT)
!
!     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
!     FOR LPPLAN, AND EQUATOR AND EQUINOX OF DATE FOR TPPLAN.
!
!     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
!     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
!     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
!
!
!     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A PLANET.
entry tpplan (ujd,glon,glat,ht, ra,dec,dis)
locatn = 1
icoord = 1
go to 40

!     THIS ENTRY PROVIDES THE LOCAL PLACE OF A PLANET.
entry lpplan (ujd,glon,glat,ht, ra,dec,dis)
locatn = 1
icoord = 0

40 if ( ujd > 100.d0 ) then
    deltat = ( ttjd - ujd ) * 86400.d0
else
    gast = dmod ( ujd, 24.d0 )
    if ( gast < 0.d0 ) gast = gast + 24.d0
    call placst ( gast )
    deltat = 0.d0
end if
observ(1) = glon
observ(2) = glat
observ(3) = ht
observ(4) = deltat

call place (ttjd,object,locatn,icoord,star,observ, skypos)

ra  = skypos(4)
dec = skypos(5)
dis = skypos(6)

return

end

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE ICRS POSITION OF A STAR,
!  GIVEN ITS APPARENT PLACE AT DATE TJD.  PROPER MOTION, PARALLAX,
!  AND RADIAL VELOCITY ARE ASSUMED TO BE ZERO.
!
!      TJD    = TT JULIAN DATE OF APPARENT PLACE (IN)
!      N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
!               (NO LONGER USED)
!      RA     = APPARENT RIGHT ASCENSION IN HOURS, REFERRED TO
!               TRUE EQUATOR AND EQUINOX OF DATE (IN)
!      DEC    = APPARENT DECLINATION IN DEGREES, REFERRED TO
!               TRUE EQUATOR OF DATE (IN)
!      RAI    = ICRS RIGHT ASCENSION IN HOURS (OUT)
!      DECI   = ICRS DECLINATION IN DEGREES (OUT)

subroutine mpstar (tjd,n,ra,dec, rai,deci)
double precision tjd,ra,dec,rai,deci,t0,t1,rainew,dcinew, &
     raiold,dciold,star,observ,skypos,r,d,p,v,delra,deldec, &
     dabs
dimension star(20), observ(6), skypos(10), p(3), v(3)
save

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data star, observ, skypos / 36 * 0.d0 /

3 format ( ' MPSTAR: NO CONVERGENCE AT COORDINATES ', &
     f9.5, 1x, sp,f9.4, ' AT JD ', ss, f10.1 )

t1 = tjd

!     GET INITIAL APPROXIMATION
iter = 0
call vectrs (ra,dec,0.d0,0.d0,0.d0,0.d0,p,v)
call preces (t1,p,t0,v)
call angles (v,rainew,dcinew)

!     ITERATIVELY FIND ICRS COORDINATES THAT PRODUCE INPUT
!     APPARENT PLACE OF STAR AT DATE TJD
do
    iter = iter + 1
    raiold = rainew
    dciold = dcinew
    star(1) = raiold
    star(2) = dciold
    star(3) = 0.d0
    star(4) = 0.d0
    star(5) = 0.d0
    star(6) = 0.d0
    call place (t1,'STAR',0,1,star,observ, skypos)
    r = skypos(4)
    d = skypos(5)
    delra = r - ra
    deldec = d - dec
    if (delra<-12.d0) delra = delra + 24.d0
    if (delra>+12.d0) delra = delra - 24.d0
    rainew = raiold - delra
    dcinew = dciold - deldec
    if (iter>30) then
        write ( *, 3 ) ra, dec, tjd
        exit
    end if
    if (dabs(rainew-raiold)>1.d-12) cycle
    if (dabs(dcinew-dciold)>1.d-11) cycle
    exit
end do

rai = rainew
deci = dcinew
if (rai< 0.d0) rai = rai + 24.d0
if (rai>=24.d0) rai = rai - 24.d0

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call vectrs (rai,deci,0.d0,0.d0,0.d0,0.d0,p,v)
call setvec (p)

end

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE GREENWICH SIDEREAL TIME
!  (EITHER MEAN OR APPARENT) AT JULIAN DATE TJDH + TJDL.
!
!       TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
!       TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
!                THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
!                FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
!                PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
!                FRACTIONAL PART
!       K      = TIME SELECTION CODE (IN)
!                SET K=0 FOR GREENWICH MEAN SIDEREAL TIME
!                SET K=1 FOR GREENWICH APPARENT SIDEREAL TIME
!       GST    = GREENWICH (MEAN OR APPARENT) SIDEREAL TIME
!                IN HOURS (OUT)
!
!  NOTE:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
!  (DELTA-T = TT - UT1) TO BE USED HERE.

subroutine sidtim ( tjdh, tjdl, k,   gst )

double precision tjdh,tjdl,gst,pi,degcon,deltat, &
     t0,utjd,ttjd,tdbjd,secdif,a,theta,rcio, &
     unitx,w1,w2,x,y,z,eq,haeq,ee,dmod,datan2
dimension unitx(3), w1(3), w2(3), x(3), y(3), z(3), eq(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( degcon = 180.d0 / pi           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data unitx / 1.d0, 0.d0, 0.d0 /

3 format ( ' SIDTIM ERROR: CANNOT RETURN SIDEREAL TIME FOR ', &
     'JD ', f10.1 )

call getdt ( deltat )

!     TIME ARGUMENT FOR PRECESSION AND NUTATION COMPONENTS OF SIDEREAL
!     TIME IS TDB
utjd = tjdh + tjdl
ttjd = utjd + deltat
tdbjd = ttjd
call times ( tdbjd, a,   secdif )
tdbjd = ttjd + secdif / 86400.d0

!     GET METHOD/ACCURACY MODE
call getmod ( mode )

if ( mode >= 2 ) then
!          EQUINOX-BASED MODE
!          SEE USNO CIRCULAR 179, SECTION 2.6.2

!         GET -1 TIMES THE MEAN OR TRUE RIGHT ASCENSION OF THE CIO
    call eqxra ( tdbjd, k,   rcio )
!         GET EARTH ROTATION ANGLE
    call erot ( tjdh, tjdl,   theta )
!         COMBINE TO OBTAIN SIDEREAL TIME
    gst = dmod ( theta / 15.d0 - rcio, 24.d0 )
    if ( gst < 0.d0 ) gst = gst + 24.d0

else
!         CIO-BASED MODE
!         SEE USNO CIRCULAR 179, SECTION 6.5.4

!         GET EARTH ROTATION ANGLE
    call erot ( tjdh, tjdl,   theta )
!         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
!         INTERMEDIATE SYSTEM
    call cioloc ( tdbjd,   rcio, kcio )
    if ( rcio == 99.d0 ) then
        write ( *, 3 ) tdbjd
        gst = 99.d0
        return
    end if
    call ciobas ( tdbjd, rcio, kcio,   x, y, z )
!         COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
    call nutate ( -tdbjd, unitx,   w1 )
    call preces ( tdbjd, w1, t0,   w2 )
    call frame ( w2, -1,    eq )
!         COMPUTE THE HOUR ANGLE OF THE EQUINOX WRT THE TIO MERIDIAN
!         (NEAR GREENWICH, BUT PASSES THROUGH THE CIP AND TIO)
    haeq = theta - datan2 ( eq(1)*y(1) + eq(2)*y(2) + eq(3)*y(3), &
                            eq(1)*x(1) + eq(2)*x(2) + eq(3)*x(3) ) &    !
                   * degcon

!         FOR MEAN SIDEREAL TIME, OBTAIN THE EQUATION OF THE EQUINOXES
!         AND SUBTRACT IT
    if ( k == 0 ) then
        call etilt ( tdbjd,   a, a, ee, a, a )
        haeq = haeq - ee / 240.d0
    end if

    haeq = dmod ( haeq, 360.d0 ) / 15.d0
    if ( haeq < 0.d0 ) haeq = haeq + 24.d0
    gst = haeq

end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE TRUE RIGHT ASCENSION OF THE CELESTIAL
!  INTERMEDIATE ORIGIN (CIO) AT A GIVEN TT JULIAN DATE.  THIS IS
!  -(EQUATION OF THE ORIGINS).
!
!       TJD    = TT JULIAN DATE (IN)
!       RACIO  = RIGHT ASCENSION OF THE CIO, WITH RESPECT TO THE
!                TRUE EQUINOX OF DATE, IN HOURS (+ OR -) (OUT)

subroutine ciora ( tjd,   racio )
double precision tjd,racio,pi,degcon,t0,t,secdif,tdbjd,rcio, &
     unitx,w1,w2,x,y,z,eq,az,datan2
dimension unitx(3), w1(3), w2(3), x(3), y(3), z(3), eq(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( degcon = 180.d0 / pi           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data unitx / 1.d0, 0.d0, 0.d0 /

3 format ( ' CIORA ERROR: CANNOT RETURN CIO RA VALUE FOR JD ', &
     f10.1 )

!     TDBJD IS THE TDB JULIAN DATE
tdbjd = tjd
call times ( tdbjd, t,   secdif )
tdbjd = tjd + secdif / 86400.d0

!     OBTAIN THE BASIS VECTORS, IN THE GCRS, FOR THE CELESTIAL
!     INTERMEDIATE SYSTEM DEFINED BY THE CIP (IN THE Z DIRECTION) AND
!     THE CIO (IN THE X DIRECTION)
call cioloc ( tdbjd,   rcio, kcio )
if ( rcio == 99.d0 ) then
    write ( *, 3 ) tdbjd
    racio = 99.d0
    return
end if
call ciobas ( tdbjd, rcio, kcio,   x, y, z )

!     COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
call nutate ( -tdbjd, unitx,   w1 )
call preces ( tdbjd, w1, t0,   w2 )
call frame ( w2, -1,    eq )

!     COMPUTE THE INTERMEDIATE RA OF THE TRUE EQUINOX (EQUATION OF
!     THE ORIGINS)
az = datan2 ( eq(1)*y(1) + eq(2)*y(2) + eq(3)*y(3), &
              eq(1)*x(1) + eq(2)*x(2) + eq(3)*x(3) ) * degcon

!     THE RA OF THE CIO IS MINUS THIS COORDINATE
racio = -az / 15.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE ROTATES A VECTOR FROM THE TERRESTRIAL TO THE
!  CELESTIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
!  ITRS (A ROTATING EARTH-FIXED SYSTEM) TO THE GCRS (A LOCAL
!  SPACE-FIXED SYSTEM) BY APPLYING ROTATIONS FOR POLAR MOTION,
!  EARTH ROTATION, NUTATION, PRECESSION, AND THE DYNAMICAL-TO-GCRS
!  FRAME TIE.
!
!       TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
!       TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
!                THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
!                FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
!                PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
!                FRACTIONAL PART
!       XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
!                IN ARCSECONDS (IN)
!       YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
!                IN ARCSECONDS (IN)
!       VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
!                SYSTEM) (IN)
!       VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
!                SYSTEM) (OUT)
!
!  NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
!
!  NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
!  (DELTA-T = TT - UT1) TO BE USED HERE.
!
!  NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
!  (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
!  SPECIAL OPTION WHERE THE OUTPUT VECTOR IS PRODUCED WITH RESPECT
!  TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH THE
!  TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
!  WORKS ONLY IN 'EQUINOX' MODE.
!
!  NOTE 4: INPUT PARAMETERS XP, YP WERE XPOLE, YPOLE IN NOVAS F3.0.
!  THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
!  IERS CONVENTIONS.

subroutine tercel ( tjdh, tjdl, xp, yp, vec1,   vec2 )

double precision tjdh,tjdl,xp,yp,vec1,vec2, &
     t0,deltat,utjdh,utjdl,utjd,ttjd,tdbjd,t,secdif, &
     gast,rcio,theta,v1,v2,v3,v4,x,y,z
dimension vec1(3), vec2(3), v1(3), v2(3), v3(3), v4(3), &
     x(3), y(3), z(3)

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

call getdt ( deltat )

if ( tjdh >= 0.d0 ) then
    utjdh = tjdh
    utjdl = tjdl
else
    utjdh = tjdl
    utjdl = 0.d0
end if
utjd = utjdh + utjdl

!     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
ttjd = utjd + deltat
tdbjd = ttjd
call times ( tdbjd, t,   secdif )
tdbjd = ttjd + secdif / 86400.d0

!     GET METHOD/ACCURACY MODE
call getmod ( mode )

if ( mode >= 2 ) then
!         'EQUINOX' MODE

!         APPLY POLAR MOTION
    if ( xp == 0.d0 .and. yp == 0.d0 ) then
        v1(1) = vec1(1)
        v1(2) = vec1(2)
        v1(3) = vec1(3)
    else
        call wobble ( tdbjd, xp, yp, vec1,   v1 )
    end if

!         APPLY EARTH ROTATION
    call sidtim ( utjdh, utjdl, 1,   gast )
    call spin ( -gast * 15.d0, v1,   v2 )

!         SPECIAL OPTION SKIPS REMAINING TRANSFORMATIONS
    if ( tjdh < 0.d0 ) then
        vec2(1) = v2(1)
        vec2(2) = v2(2)
        vec2(3) = v2(3)
    else

!         APPLY NUTATION AND PRECESSION
    call nutate ( -tdbjd, v2,   v3 )
    call preces ( tdbjd, v3, t0,   v4 )

!         APPLY FRAME-TIE MATRIX
    call frame ( v4, -1, vec2 )

    end if

else
!         'CIO-TIO-THETA' MODE
!         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
!         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
!         EQ. (3) AND (4).

!         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE TERRESTRIAL
!         INTERMEDIATE SYSTEM
    if ( xp == 0.d0 .and. yp == 0.d0 ) then
        v1(1) = vec1(1)
        v1(2) = vec1(2)
        v1(3) = vec1(3)
    else
        call wobble ( tdbjd, xp, yp, vec1,   v1 )
    end if

!         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
!         INTERMEDIATE SYSTEM
    call cioloc ( tdbjd,   rcio, kcio )
    call ciobas ( tdbjd, rcio, kcio,   x, y, z )

!         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
!         THE VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
    call erot ( utjdh, utjdl,   theta )
    call spin ( -theta, v1,   v2 )

!         TRANSFORM THE VECTOR FROM THE CELESTIAL INTERMEDIATE SYSTEM
!         TO THE GCRS
    vec2(1) = x(1) * v2(1) + y(1) * v2(2) + z(1) * v2(3)
    vec2(2) = x(2) * v2(1) + y(2) * v2(2) + z(2) * v2(3)
    vec2(3) = x(3) * v2(1) + y(3) * v2(2) + z(3) * v2(3)

end if

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( vec2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE ROTATES A VECTOR FROM THE CELESTIAL TO THE
!  TERRESTRIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
!  GCRS (A LOCAL SPACE-FIXED SYSTEM) TO THE ITRS (A ROTATING
!  EARTH-FIXED SYSTEM) BY APPLYING ROTATIONS FOR THE GCRS-TO-
!  DYNAMICAL FRAME TIE, PRECESSION, NUTATION, EARTH ROTATION,
!  AND POLAR MOTION.
!
!       TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
!       TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
!                THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
!                FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
!                PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
!                FRACTIONAL PART
!       XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
!                IN ARCSECONDS (IN)
!       YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
!                IN ARCSECONDS (IN)
!       VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
!                SYSTEM) (IN)
!       VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
!                SYSTEM) (OUT)
!
!  NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
!
!  NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
!  (DELTA-T = TT - UT1) TO BE USED HERE.
!
!  NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
!  (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
!  SPECIAL OPTION WHERE THE INPUT VECTOR IS ASSUMED TO BE WITH
!  RESPECT TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH
!  THE TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
!  WORKS ONLY IN 'EQUINOX' MODE.

subroutine celter ( tjdh, tjdl, xp, yp, vec1,   vec2 )

double precision tjdh,tjdl,xp,yp,vec1,vec2, &
     t0,deltat,utjdh,utjdl,utjd,ttjd,tdbjd,t,secdif, &
     gast,rcio,theta,v1,v2,v3,v4,x,y,z
dimension vec1(3), vec2(3), v1(3), v2(3), v3(3), v4(3), &
     x(3), y(3), z(3)

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

call getdt ( deltat )

if ( tjdh >= 0.d0 ) then
    utjdh = tjdh
    utjdl = tjdl
else
    utjdh = tjdl
    utjdl = 0.d0
end if
utjd = utjdh + utjdl

!     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
ttjd = utjd + deltat
tdbjd = ttjd
call times ( tdbjd, t,   secdif )
tdbjd = ttjd + secdif / 86400.d0

!     GET METHOD/ACCURACY MODE
call getmod ( mode )

if ( mode >= 2 ) then
!         'EQUINOX' MODE

!         SPECIAL OPTION SKIPS INITIAL TRANSFORMATIONS
    if ( tjdh < 0.d0 ) then
        v3(1) = vec1(1)
        v3(2) = vec1(2)
        v3(3) = vec1(3)
    else

!         APPLY FRAME-TIE MATRIX
    call frame ( vec1, 1, v1 )

!         APPLY PRECESSION AND NUTATION
    call preces ( t0, v1, tdbjd,   v2 )
    call nutate ( tdbjd, v2,   v3 )

    end if

!         APPLY EARTH ROTATION
    call sidtim ( utjdh, utjdl, 1,   gast )
    call spin ( gast * 15.d0, v3,   v4 )

!         APPLY POLAR MOTION
    if ( xp == 0.d0 .and. yp == 0.d0 ) then
        vec2(1) = v4(1)
        vec2(2) = v4(2)
        vec2(3) = v4(3)
    else
        call wobble ( -tdbjd, xp, yp, v4,   vec2 )
    end if

else
!         'CIO-TIO-THETA' MODE
!         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
!         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
!         EQ. (3) AND (4).

!         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
!         INTERMEDIATE SYSTEM
    call cioloc ( tdbjd,   rcio, kcio )
    call ciobas ( tdbjd, rcio, kcio,   x, y, z )

!         TRANSFORM THE VECTOR FROM THE GCRS TO THE
!         CELESTIAL INTERMEDIATE SYSTEM
    v1(1) = x(1) * vec1(1) + x(2) * vec1(2) + x(3) * vec1(3)
    v1(2) = y(1) * vec1(1) + y(2) * vec1(2) + y(3) * vec1(3)
    v1(3) = z(1) * vec1(1) + z(2) * vec1(2) + z(3) * vec1(3)

!         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
!         THE VECTOR TO THE TERRESTRIAL INTERMEDIATE SYSTEM
    call erot ( utjdh, utjdl,   theta )
    call spin ( theta, v1,   v2 )

!         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE ITRS
    if ( xp == 0.d0 .and. yp == 0.d0 ) then
        vec2(1) = v2(1)
        vec2(2) = v2(2)
        vec2(3) = v2(3)
    else
        call wobble ( -tdbjd, xp, yp, v2,   vec2 )
    end if

end if

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( vec2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS HIPPARCOS DATA AT EPOCH J1991.25
!  TO EPOCH J2000.0.  TO BE USED ONLY FOR HIPPARCOS OR TYCHO STARS
!  WITH LINEAR SPACE MOTION.  BOTH INPUT AND OUTPUT DATA IS IN THE
!  ICRS.
!
!       RAH    = HIPPARCOS RIGHT ASCENSION IN DEGREES (IN)
!       DECH   = HIPPARCOS DECLINATION IN DEGREES (IN)
!       PMRAH  = HIPPARCOS PROPER MOTION IN RA
!                IN MILLIARCSECONDS/YEAR (IN)
!       PMDECH = HIPPARCOS PROPER MOTION IN DEC
!                IN MILLIARCSECONDS/YEAR (IN)
!       PARXH  = HIPPARCOS PARALLAX IN MILLIARCSECONDS (IN)
!       RVH    = RADIAL VELOCITY AT HIPPARCOS EPOCH
!                IN KILOMETERS/SECOND (IN)
!       RA2    = RIGHT ASCENSION AT J2000.0 IN HOURS (OUT)
!       DEC2   = DECLINATION AT J2000.0 IN DEGREES (OUT)
!       PMRA2  = PROPER MOTION IN RA AT J2000.0
!                IN MILLIARCSECONDS/YEAR (OUT)
!       PMDEC2 = PROPER MOTION IN DEC AT J2000.0
!                IN MILLIARCSECONDS/YEAR (OUT)
!       PARX2  = PARALLAX AT J2000.0 IN MILLIARCSECONDS (OUT)
!       RV2    = RADIAL VELOCITY AT J2000.0 IN KILOMETERS/SECOND
!                (OUT)
!
!  NOTE:  INPUT RA IS IN DEGREES, AS PER HIPPARCOS, BUT OUTPUT RA
!  IS IN HOURS.

subroutine gethip ( rah, dech, pmrah, pmdech, parxh, rvh, &
    ra2, dec2, pmra2, pmdec2, parx2, rv2 )

double precision rah,dech,pmrah,pmdech,parxh,rvh, &
     ra2,dec2,pmra2,pmdec2,parx2,rv2,epoch1,epoch2

data epoch1, epoch2 / 2448349.0625d0, 2451545.0000d0 /

call catran ( 1,epoch1,rah/15.d0,dech,pmrah,pmdech,parxh,rvh, &
                epoch2,ra2,      dec2,pmra2,pmdec2,parx2,rv2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE TRANSFORMS A STAR'S CATALOG QUANTITIES FOR
!  A CHANGE OF EPOCH AND/OR EQUATOR AND EQUINOX.  IT CAN ALSO BE
!  USED TO ROTATE CATALOG QUANTITIES ON THE DYNAMICAL EQUATOR AND
!  EQUINOX OF J2000.0 TO THE ICRS OR VICE VERSA.
!
!       IT     = TRANSFORMATION OPTION (IN)
!                SET IT=1 TO CHANGE EPOCH (SAME EQUATOR AND EQUINOX)
!                SET IT=2 TO CHANGE EQUATOR AND EQUINOX (SAME EPOCH)
!                SET IT=3 TO CHANGE EQUATOR AND EQUINOX AND EPOCH
!                SET IT=4 TO CHANGE EQUATOR AND EQUINOX OF J2000.0
!                         TO ICRS
!                SET IT=5 TO CHANGE ICRS TO EQUATOR AND EQUINOX OF
!                         J2000.0
!       DATE1  = TT JULIAN DATE, OR YEAR, OF ORIGINAL CATALOG
!                DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
!       RA1    = ORIGINAL MEAN RIGHT ASCENSION IN HOURS (IN)
!       DEC1   = ORIGINAL MEAN DECLINATION IN DEGREES (IN)
!       PMRA1  = ORIGINAL PROPER MOTION IN RA
!                IN MILLIARCSECONDS/YEAR (IN)
!       PMDEC1 = ORIGINAL PROPER MOTION IN DEC
!                IN MILLIARCSECONDS/YEAR (IN)
!       PARX1  = ORIGINAL PARALLAX IN MILLIARCSECONDS (IN)
!       RV1    = ORIGINAL RADIAL VELOCITY IN KILOMETERS/SECOND
!                (IN)
!       DATE2  = TT JULIAN DATE, OR YEAR, FOR TRANSFORMED
!                OUTPUT DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
!       RA2    = TRANSFORMED MEAN RIGHT ASCENSION IN HOURS (OUT)
!       DEC2   = TRANSFORMED MEAN DECLINATION IN DEGREES (OUT)
!       PMRA2  = TRANSFORMED PROPER MOTION IN RA
!                IN MILLIARCSECONDS/YEAR (OUT)
!       PMDEC2 = TRANSFORMED PROPER MOTION IN DEC
!                IN MILLIARCSECONDS/YEAR (OUT)
!       PARX2  = TRANSFORMED PARALLAX IN MILLIARCSECONDS (OUT)
!       RV2    = TRANSFORMED RADIAL VELOCITY IN KILOMETERS/SECOND
!                (OUT)
!
!  NOTE 1:  DATE1 AND DATE2 MAY BE SPECIFIED EITHER AS A JULIAN
!  DATE (E.G., 2433282.5D0) OR A JULIAN YEAR AND FRACTION
!  (E.G., 1950.0D0).  VALUES LESS THAN 10000 ARE ASSUMED TO
!  BE YEARS.  FOR IT=2 OR IT=3, EITHER DATE1 OR DATE2 MUST BE
!  2451545.0 OR 2000.0 (J2000.0).  FOR IT=4 AND IT=5, DATE1 AND
!  DATE2 ARE IGNORED.
!
!  NOTE 2:  IT=1 UPDATES THE STAR'S DATA TO ACCOUNT FOR
!  THE STAR'S SPACE MOTION BETWEEN THE FIRST AND SECOND DATES,
!  WITHIN A FIXED REFERENCE SYSTEM.  IT=2 APPLIES A ROTATION
!  OF THE REFERENCE SYSTEM CORRESPONDING TO PRECESSION BETWEEN
!  THE FIRST AND SECOND DATES, BUT LEAVES THE STAR FIXED IN SPACE.
!  IT=3 PROVIDES BOTH TRANSFORMATIONS.  IT=4 AND IT=5 PROVIDE A
!  A FIXED ROTATION ABOUT VERY SMALL ANGLES (<0.1 ARCSECOND) TO
!  TAKE DATA FROM THE DYNAMICAL SYSTEM OF J2000.0 TO THE ICRS (IT=4)
!  OR VICE VERSA (IT=5).
!
!  NOTE 3:  FOR IT=1, INPUT DATA CAN BE IN ANY FIXED REFERENCE
!  SYSTEM. FOR IT=2 OR IT=3, THIS SUBROUTINE ASSUMES THE INPUT DATA
!  IS IN THE DYNAMICAL SYSTEM AND PRODUCES OUTPUT IN THE DYNAMICAL
!  SYSTEM.  FOR IT=4, THE INPUT DATA MUST BE ON THE DYNAMICAL EQUATOR
!  AND EQUINOX OF J2000.0.  FOR IT=5, THE INPUT DATA MUST BE IN THE
!  ICRS.
!
!  NOTE 4:  THIS SUBROUTINE CANNOT BE PROPERLY USED TO BRING DATA
!  FROM OLD STAR CATALOGS INTO THE MODERN SYSTEM, BECAUSE
!  OLD CATALOGS WERE COMPILED USING A SET OF CONSTANTS THAT ARE
!  INCOMPATIBLE WITH MODERN VALUES.  IN PARTICULAR, IT SHOULD NOT
!  BE USED FOR CATALOGS WHOSE POSITIONS AND PROPER MOTIONS WERE
!  DERIVED BY ASSUMING A PRECESSION CONSTANT SIGNIFICANTLY DIFFERENT
!  FROM THE VALUE IMPLICIT IN SUBROUTINE PRECES.

subroutine catran ( it, &
                    date1, ra1, dec1, pmra1, pmdec1, parx1, rv1, &
                    date2, ra2, dec2, pmra2, pmdec2, parx2, rv2 )

double precision date1,ra1,dec1,pmra1,pmdec1,parx1,rv1, &
     date2,ra2,dec2,pmra2,pmdec2,parx2,rv2, &
     pi,seccon,aukm,c,tjd1,pos1,vel1,tjd2,pos2,vel2, &
     paralx,dist,r,d,cra,sra,cdc,sdc,k,pmr,pmd,rvl, &
     xyproj,dcos,dsin,datan2
integer it,j

dimension pos1(3), vel1(3), pos2(3), vel2(3)

save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

data ntimes / 0 /

ntimes = ntimes + 1
if ( ntimes == 1 ) then
!         GET LENGTH OF AU IN KILOMETERS
    call astcon ( 'AU', 1.d-3, aukm )
!         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
    call astcon ( 'C', 1.d-3,   c )
end if

! --- IF NECESSARY, COMPUTE JULIAN DATES ------------------------------

!     SUBROUTINE USES TDB JULIAN DATES INTERNALLY, BUT NO
!     DISTINCTION BETWEEN TDB AND TT IS NECESSARY

if ( date1 < 10000.d0 ) then
     tjd1 = 2451545.0d0 + ( date1 - 2000.d0 ) * 365.25d0
else
     tjd1 = date1
end if
if ( date2 < 10000.d0 ) then
     tjd2 = 2451545.0d0 + ( date2 - 2000.d0 ) * 365.25d0
else
     tjd2 = date2
end if

! --- CONVERT INPUT ANGULAR COMPONENTS TO VECTORS ---------------------

!     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
!     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
paralx = parx1
if ( paralx <= 0.d0 ) paralx = 1.d-6

!     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
!     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
dist = 1.d0 / dsin ( paralx * 1.d-3 / seccon )
r = ra1 * 54000.d0 / seccon
d = dec1 * 3600.d0 / seccon
cra = dcos ( r )
sra = dsin ( r )
cdc = dcos ( d )
sdc = dsin ( d )
pos1(1) = dist * cdc * cra
pos1(2) = dist * cdc * sra
pos1(3) = dist * sdc

!     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
!     LIGHT TRAVEL TIME TO STAR
k = 1.d0 / ( 1.d0 - rv1 / c )

!     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL
!     COMPONENTS OF MOTION, IN SPHERICAL POLAR SYSTEM AT STAR'S
!     ORIGINAL POSITION, WITH UNITS OF AU/DAY
pmr = pmra1  / ( paralx * 365.25d0 ) * k
pmd = pmdec1 / ( paralx * 365.25d0 ) * k
rvl = rv1 * 86400.d0 / aukm          * k

!     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
vel1(1) = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra
vel1(2) =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra
vel1(3) =               pmd * cdc       + rvl * sdc

! --- UPDATE STAR'S POSITION VECTOR FOR SPACE MOTION ------------------
!     (ONLY IF IT=1 OR IT=3)

if ( it == 1 .or. it == 3 ) then
    do j=1,3
        pos2(j) = pos1(j) + vel1(j) * ( tjd2 - tjd1 )
        vel2(j) = vel1(j)
    end do
else
    do j=1,3
        pos2(j) = pos1(j)
        vel2(j) = vel1(j)
    end do
end if

! --- PRECESS POSITION AND VELOCITY VECTORS ---------------------------
!     (ONLY IF IT=2 OR IT=3)

if ( it == 2 .or. it == 3 ) then
    do j=1,3
        pos1(j) = pos2(j)
        vel1(j) = vel2(j)
    end do
    call preces ( tjd1, pos1, tjd2,    pos2 )
    call preces ( tjd1, vel1, tjd2,    vel2 )
end if

! --- ROTATE DYNAMICAL J2000.0 POSITION AND VELOCITY VECTORS TO ICRS --
!     (ONLY IF IT=4)

if ( it == 4 ) then
    call frame ( pos1, -1,    pos2 )
    call frame ( vel1, -1,    vel2 )
end if

! --- ROTATE ICRS POSITION AND VELOCITY VECTORS TO DYNAMICAL J2000.0 --
!     (ONLY IF IT=5)

if ( it == 5 ) then
    call frame ( pos1, 1,    pos2 )
    call frame ( vel1, 1,    vel2 )
end if

! --- CONVERT VECTORS BACK TO ANGULAR COMPONENTS FOR OUTPUT -----------

!     FROM UPDATED POSITION VECTOR, OBTAIN STAR'S NEW POSITION
!     EXPRESSED AS ANGULAR QUANTITIES
xyproj = dsqrt ( pos2(1)**2 + pos2(2)**2 )
r = 0.d0
if ( xyproj > 0.d0 ) r = datan2 ( pos2(2), pos2(1) )
ra2 = r * seccon / 54000.d0
if ( ra2 <  0.d0 ) ra2 = ra2 + 24.d0
if ( ra2 >= 24.d0 ) ra2 = ra2 - 24.d0
d = datan2 ( pos2(3), xyproj  )
dec2 = d * seccon / 3600.d0
dist = dsqrt ( pos2(1)**2 + pos2(2)**2 + pos2(3)**2 )
paralx = dasin ( 1.d0 / dist ) * seccon * 1.d3
parx2 = paralx

!     TRANSFORM MOTION VECTOR BACK TO SPHERICAL POLAR SYSTEM AT STAR'S
!     NEW POSITION
cra = dcos ( r )
sra = dsin ( r )
cdc = dcos ( d )
sdc = dsin ( d )
pmr = - vel2(1) * sra       + vel2(2) * cra
pmd = - vel2(1) * cra * sdc - vel2(2) * sra * sdc + vel2(3) * cdc
rvl =   vel2(1) * cra * cdc + vel2(2) * sra * cdc + vel2(3) * sdc

!     CONVERT COMPONENTS OF MOTION FROM AU/DAY TO NORMAL
!     CATALOG UNITS
pmra2  = pmr * paralx * 365.25d0   / k
pmdec2 = pmd * paralx * 365.25d0   / k
rv2    = rvl * ( aukm / 86400.d0 ) / k

!     TAKE CARE OF ZERO-PARALLAX CASE
if ( parx2 <= 1.01d-6 ) then
    parx2 = 0.d0
    rv2 = rv1
end if

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( pos2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE TRANSFORMS TOPOCENTRIC RIGHT ASCENSION AND
!  DECLINATION TO ZENITH DISTANCE AND AZIMUTH.  THIS ROUTINE USES
!  A METHOD THAT PROPERLY ACCOUNTS FOR POLAR MOTION, WHICH IS
!  SIGNIFICANT AT THE SUB-ARCSECOND LEVEL.  THIS SUBROUTINE
!  CAN ALSO ADJUST COORDINATES FOR ATMOSPHERIC REFRACTION.
!
!       UJD    = UT1 JULIAN DATE (IN)
!       XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
!                ARCSECONDS (IN)
!       YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
!                ARCSECONDS (IN)
!       GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
!                IN DEGREES (IN)
!       GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
!                IN DEGREES (IN)
!       HT     = HEIGHT OF OBSERVER IN METERS (IN)
!       RA     = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
!                IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
!                OF DATE (IN)
!       DEC    = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
!                IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE (IN)
!       IREFR  = ATMOSPHERIC REFRACTION OPTION (IN)
!                SET IREFR=0 FOR NO REFRACTION
!                SET IREFR=1 TO INCLUDE REFRACTION
!       ZD     = TOPOCENTRIC ZENITH DISTANCE IN DEGREES,
!                AFFECTED BY REFRACTION IF IREFR=1 (OUT)
!       AZ     = TOPOCENTRIC AZIMUTH (MEASURED EAST FROM NORTH)
!                IN DEGREES (OUT)
!       RAR    = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
!                IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
!                OF DATE, AFFECTED BY REFRACTION IF IREFR=1 (OUT)
!       DECR   = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
!                IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE,
!                AFFECTED BY REFRACTION IF IREFR=1 (OUT)
!
!  NOTE 1:  XP AND YP CAN BE SET TO ZERO IF SUB-ARCSECOND ACCURACY IS
!  NOT NEEDED.  HT IS USED ONLY FOR REFRACTION, IF IREFR=1.  RA AND
!  DEC CAN BE OBTAINED FROM TPSTAR, TPPLAN, OR PLACE.
!
!  NOTE 2:  THE DIRECTONS ZD=0 (ZENITH) AND AZ=0 (NORTH) ARE
!  HERE CONSIDERED FIXED IN THE TERRESTRIAL SYSTEM.  SPECIFICALLY,
!  THE ZENITH IS ALONG THE GEODETIC NORMAL, AND NORTH IS TOWARD
!  THE ITRS POLE.
!
!  NOTE 3:  IF IREFR=0, THEN RAR=RA AND DECR=DEC.
!
!  NOTE 4: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
!  THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
!  IERS CONVENTIONS.

subroutine zdaz ( ujd, xp, yp, glon, glat, ht, ra, dec, irefr, &
    zd, az, rar, decr )

double precision ujd,xp,yp,glon,glat,ht,ra,dec,zd,az,rar,decr, &
     pi,degrad,raddeg, &
     sinlat,coslat,sinlon,coslon,sindc,cosdc,sinra,cosra, &
     uze,une,uwe,uz,un,uw,p,pr,pz,pn,pw,proj, &
     zd0,zd1,refr,sinzd,coszd,sinzd0,coszd0, &
     dsin,dcos,dsqrt,datan2
dimension uze(3), une(3), uwe(3), uz(3), un(3), uw(3), &
     p(3), pr(3)

parameter ( pi     = 3.14159265358979324d0 )
parameter ( degrad = pi / 180.d0           )
parameter ( raddeg = 180.d0 / pi           )

rar    = ra
decr   = dec
sinlat = dsin ( glat * degrad )
coslat = dcos ( glat * degrad )
sinlon = dsin ( glon * degrad )
coslon = dcos ( glon * degrad )
sindc  = dsin ( dec * degrad )
cosdc  = dcos ( dec * degrad )
sinra  = dsin ( ra * 15.d0 * degrad )
cosra  = dcos ( ra * 15.d0 * degrad )

! --- SET UP ORTHONORMAL BASIS VECTORS IN LOCAL EARTH-FIXED SYSTEM ----

!     DEFINE VECTOR TOWARD LOCAL ZENITH IN EARTH-FIXED SYSTEM (Z AXIS)
uze(1) =  coslat * coslon
uze(2) =  coslat * sinlon
uze(3) =  sinlat

!     DEFINE VECTOR TOWARD LOCAL NORTH IN EARTH-FIXED SYSTEM (X AXIS)
une(1) = -sinlat * coslon
une(2) = -sinlat * sinlon
une(3) =  coslat

!     DEFINE VECTOR TOWARD LOCAL WEST IN EARTH-FIXED SYSTEM (Y AXIS)
uwe(1) =  sinlon
uwe(2) = -coslon
uwe(3) =  0.d0

! --- OBTAIN VECTORS IN CELESTIAL SYSTEM ------------------------------

!     ROTATE EARTH-FIXED ORTHONORMAL BASIS VECTORS TO CELESTIAL SYSTEM
!     (WRT EQUATOR AND EQUINOX OF DATE)
call eqinox
call tercel ( -1.d0, ujd, xp, yp, uze,   uz )
call tercel ( -1.d0, ujd, xp, yp, une,   un )
call tercel ( -1.d0, ujd, xp, yp, uwe,   uw )
call resume

!     DEFINE UNIT VECTOR P TOWARD OBJECT IN CELESTIAL SYSTEM
!     (WRT EQUATOR AND EQUINOX OF DATE)
p(1) = cosdc * cosra
p(2) = cosdc * sinra
p(3) = sindc

! --- COMPUTE COORDINATES OF OBJECT WRT ORTHONORMAL BASIS -------------

!     COMPUTE COMPONENTS OF P -- PROJECTIONS OF P ONTO ROTATED
!     EARTH-FIXED BASIS VECTORS
pz = p(1) * uz(1) + p(2) * uz(2) + p(3) * uz(3)
pn = p(1) * un(1) + p(2) * un(2) + p(3) * un(3)
pw = p(1) * uw(1) + p(2) * uw(2) + p(3) * uw(3)

!     COMPUTE AZIMUTH AND ZENITH DISTANCE
proj = dsqrt ( pn**2 + pw**2 )
az = 0.d0
if ( proj > 0.d0 ) az = -datan2 ( pw, pn ) * raddeg
if ( az <   0.d0 ) az = az + 360.d0
if ( az >= 360.d0 ) az = az - 360.d0
zd = datan2 ( proj, pz ) * raddeg

! --- APPLY ATMOSPHERIC REFRACTION IF REQUESTED -----------------------

if ( irefr == 1 ) then

!         GET REFRACTION IN ZENITH DISTANCE
!         ITERATIVE PROCESS REQUIRED BECAUSE REFRACTION ALGORITHMS ARE
!         ALWAYS A FUNCTION OF OBSERVED (NOT COMPUTED) ZENITH DISTANCE
    zd0 = zd
    do
        zd1 = zd
        call refrac ( ht, zd,   refr )
        zd = zd0 - refr
    !         REQUIRE CONVERGENCE TO 0.1 ARCSEC (ACTUAL ACCURACY LESS)
        if ( dabs ( zd - zd1 ) <= 3.d-5 ) exit
    end do

!         APPLY REFRACTION TO CELESTIAL COORDINATES OF OBJECT
    if ( refr > 0.d0 .and. zd > 3.d-4 ) then

!             SHIFT POSITION VECTOR OF OBJECT IN CELESTIAL SYSTEM
!             TO ACCOUNT FOR FOR REFRACTION (SEE USNO/AA TECHNICAL
!             NOTE 1998-09)
        sinzd  = dsin ( zd * degrad )
        coszd  = dcos ( zd * degrad )
        sinzd0 = dsin ( zd0 * degrad )
        coszd0 = dcos ( zd0 * degrad )
!             COMPUTE REFRACTED POSITION VECTOR
        do j = 1, 3
           pr(j) = ( ( p(j) - coszd0 * uz(j) ) / sinzd0 ) * sinzd &
                    +                  uz(j)              * coszd
        end do

!             COMPUTE REFRACTED RIGHT ASCENSION AND DECLINATION
        proj = dsqrt ( pr(1)**2 + pr(2)**2 )
        rar = 0.d0
        if ( proj > 0.d0 ) rar = datan2 ( pr(2), pr(1) ) &
                                    * raddeg / 15.d0
        if ( rar <  0.d0 ) rar = rar + 24.d0
        if ( rar >= 24.d0 ) rar = rar - 24.d0
        decr = datan2 ( pr(3), proj ) * raddeg

    end if

end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS GCRS RIGHT ASCENSION AND DECLINATION
!  TO COORDINATES WITH RESPECT TO THE EQUATOR OF DATE (MEAN OR TRUE).
!  FOR COORDINATES WITH RESPECT TO THE TRUE EQUATOR OF DATE, THE
!  ORIGIN OF RIGHT ASCENSION CAN BE EITHER THE TRUE EQUINOX OR THE
!  CELESTIAL INTERMEDIATE ORIGIN (CIO).
!
!       TJD    = TT JULIAN DATE OF EQUATOR TO BE USED FOR
!                OUTPUT COORDINATES (IN)
!       ICOORD = COORDINATE SYSTEM SELECTION FOR OUTPUT
!                COORDINATES (IN)
!                SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
!       RAG    = GCRS RIGHT ASCENSION IN HOURS (IN)
!       DECG   = GCRS DECLINATION IN DEGREES (IN)
!       RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
!                EQUATOR AND RIGHT ASCENSION ORIGIN OF DATE (OUT)
!       DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
!                   EQUATOR OF DATE (OUT)

subroutine gcrseq ( tjd, icoord, rag, decg,   ra, dec )

double precision tjd,rag,decg,ra,dec,pi,radcon,t0,t1,t,secdif,r,d, &    !
     pos1,pos2,pos3,pos4,rcio,x,y,z,dsin,dcos
dimension pos1(3), pos2(3), pos3(3), pos4(3), x(3), y(3), z(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

!     T1 IS THE TDB JULIAN DATE
call times ( tjd, t,   secdif )
t1 = tjd + secdif / 86400.d0

!     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
r = rag * 15.d0 * radcon
d = decg * radcon
pos1(1) = dcos ( d ) * dcos ( r )
pos1(2) = dcos ( d ) * dsin ( r )
pos1(3) = dsin ( d )

if ( icoord <= 1 ) then

!         TRANSFORM THE POSITION VECTOR FROM GCRS TO MEAN EQUATOR AND
!         EQUINOX OF DATE
    call frame  ( pos1, 1,   pos2 )
    call preces ( t0, pos2, t1,   pos3 )
!         IF REQUESTED, TRANSFORM FURTHER TO TRUE EQUATOR AND EQUINOX
!         OF DATE
    if ( icoord == 1 ) then
        call nutate ( t1, pos3,   pos4)
    else
        pos4(1) = pos3(1)
        pos4(2) = pos3(2)
        pos4(3) = pos3(3)
    end if

else

!         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
!         INTERMEDIATE SYSTEM
    call cioloc ( t1,   rcio, kcio )
    call ciobas ( t1, rcio, kcio,   x, y, z )

!         TRANSFORM POSITION VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
!         (WHICH HAS THE CIO AS ITS ORIGIN OF RIGHT ASCENSION)
    pos4(1) = x(1) * pos1(1) + x(2) * pos1(2) + x(3) * pos1(3)
    pos4(2) = y(1) * pos1(1) + y(2) * pos1(2) + y(3) * pos1(3)
    pos4(3) = z(1) * pos1(1) + z(2) * pos1(2) + z(3) * pos1(3)

end if

call angles ( pos4,   ra, dec )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS RIGHT ASCENSION AND DECLINATION
!  TO ECLIPTIC LONGITUDE AND LATITUDE.
!
!       TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
!                USED FOR COORDINATES (IN)
!       ICOORD = COORDINATE SYSTEM SELECTION (IN)
!                SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
!                (ECLIPTIC IS ALWAYS THE MEAN PLANE)
!       RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
!                EQUATOR AND EQUINOX OF DATE (IN)
!       DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
!                EQUATOR AND EQUINOX OF DATE (IN)
!       ELON   = ECLIPTIC LONGITUDE IN DEGREES, REFERRED TO SPECIFIED
!                ECLIPTIC AND EQUINOX OF DATE (OUT)
!       ELAT   = ECLIPTIC LATITUDE IN DEGREES, REFERRED TO SPECIFIED
!                ECLIPTIC AND EQUINOX OF DATE (OUT)
!
!  NOTE:  TO CONVERT ICRS RA AND DEC TO ECLIPTIC COORDINATES (MEAN
!  ECLIPTIC AND EQUINOX OF J2000.0), SET TJD = 0.D0 AND ICOORD = 0.
!  EXCEPT FOR THE INPUT TO THIS CASE, ALL COORDINATES ARE DYNAMICAL.

subroutine eqecl ( tjd, icoord, ra, dec,   elon, elat )

double precision tjd,ra,dec,elon,elat,pi,radcon,r,d,xyproj,e, &
     pos1,pos2,dsin,dcos,dsqrt,datan2
dimension pos1(3), pos2(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
r = ra * 15.d0 * radcon
d = dec * radcon
pos1(1) = dcos ( d ) * dcos ( r )
pos1(2) = dcos ( d ) * dsin ( r )
pos1(3) = dsin ( d )

!     CONVERT THE VECTOR FROM EQUATORIAL TO ECLIPTIC SYSTEM
call eqec ( tjd, icoord, pos1,   pos2 )

!     DECOMPOSE ECLIPTIC VECTOR INTO ECLIPTIC LONGITUDE AND LATITUDE
xyproj = dsqrt ( pos2(1)**2 + pos2(2)**2 )
e = 0.d0
if ( xyproj > 0.d0 ) e = datan2 ( pos2(2), pos2(1) )
elon = e / radcon
if ( elon < 0.d0 ) elon = elon + 360.d0
e = datan2 ( pos2(3), xyproj )
elat = e / radcon

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS AN EQUATORIAL POSITION VECTOR TO
!  AN ECLIPTIC POSITION VECTOR.
!
!       TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
!                USED FOR COORDINATES (IN)
!       ICOORD = COORDINATE SYSTEM SELECTION (IN)
!                SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
!                (ECLIPTIC IS ALWAYS THE MEAN PLANE)
!       POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
!                EQUATOR AND EQUINOX OF DATE (IN)
!       POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
!                ECLIPTIC AND EQUINOX OF DATE (OUT)
!
!  NOTE:  TO CONVERT ICRS VECTORS TO ECLIPTIC VECTORS (MEAN ECLIPTIC
!  AND EQUINOX OF J2000.0 ONLY), SET TJD = 0.D0 AND ICOORD = 0.
!  EXCEPT FOR THE INPUT TO THIS CASE, ALL VECTORS ARE ASSUMED TO
!  BE WITH RESPECT TO A DYNAMICAL SYSTEM.

subroutine eqec ( tjd, icoord, pos1,   pos2 )

double precision tjd,pos1,pos2,pos0,pi,radcon,t0,t1,t,secdif, &
     tlast,ob2000,oblm,oblt,obl,x,dsin,dcos
dimension pos1(3), pos2(3), pos0(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /,   ob2000 / 0.d0 /

!     T1 IS THE TDB JULIAN DATE
call times ( tjd, t,   secdif )
t1 = tjd + secdif / 86400.d0

if ( tjd == 0.d0 ) then
!         CASE WHERE INPUT VECTOR IS IN ICRS SYSTEM
    call frame ( pos1, 1,   pos0 )
!         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY
    if ( ob2000 == 0.d0 ) call etilt ( t0,  ob2000, x, x, x, x )      !
    obl = ob2000 * radcon
else
!         CASE WHERE INPUT VECTOR IS IN EQUATOR OF DATE SYSTEM
    pos0(1) = pos1(1)
    pos0(2) = pos1(2)
    pos0(3) = pos1(3)
!         GET MEAN AND TRUE OBLIQUITY
    if ( dabs ( tjd - tlast ) > 1.d-8 ) then
        call etilt ( t1,   oblm, oblt, x, x, x )
        tlast = tjd
    end if
!         SELECT MEAN OR TRUE OBLIQUITY
    obl = oblm * radcon
    if ( icoord == 1 ) obl = oblt * radcon
end if

!     ROTATE EQUATORIAL POSITION VECTOR TO ECLIPTIC SYSTEM
pos2(1) =  pos0(1)
pos2(2) =  pos0(2) * dcos ( obl ) + pos0(3) * dsin ( obl )
pos2(3) = -pos0(2) * dsin ( obl ) + pos0(3) * dcos ( obl )

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( pos2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS AN ECLIPTIC POSITION VECTOR TO
!  AN EQUATORIAL POSITION VECTOR.
!
!       TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
!                USED FOR COORDINATES (IN)
!       ICOORD = COORDINATE SYSTEM SELECTION (IN)
!                SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
!                SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
!                (ECLIPTIC IS ALWAYS THE MEAN PLANE)
!       POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
!                ECLIPTIC AND EQUINOX OF DATE (IN)
!       POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
!                EQUATOR AND EQUINOX OF DATE (OUT)
!
!  NOTE:  TO CONVERT ECLIPTIC VECTORS (MEAN ECLIPTIC AND EQUINOX OF
!  OF J2000.0 ONLY) TO ICRS VECTORS, SET TJD = 0.D0 AND ICOORD = 0.
!  EXCEPT FOR THE OUTPUT FROM THIS CASE, ALL VECTORS ARE ASSUMED TO
!  BE WITH RESPECT TO A DYNAMICAL SYSTEM.

subroutine eceq ( tjd, icoord, pos1,   pos2 )

double precision tjd,pos1,pos2,pos0,pi,radcon,t0,t1,t,secdif, &
     tlast,ob2000,oblm,oblt,obl,x,dsin,dcos
dimension pos1(3), pos2(3), pos0(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /,   ob2000 / 0.d0 /

!     T1 IS THE TDB JULIAN DATE
call times ( tjd, t,   secdif )
t1 = tjd + secdif / 86400.d0

if ( tjd == 0.d0 ) then
!         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
!         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY
    if ( ob2000 == 0.d0 ) call etilt ( t0,  ob2000, x, x, x, x )      !
    obl = ob2000 * radcon
else
!         CASE WHERE OUTPUT VECTOR IS TO BE IN EQUATOR OF DATE SYSTEM
!         GET MEAN AND TRUE OBLIQUITY
    if ( dabs ( tjd - tlast ) > 1.d-8 ) then
        call etilt ( t1,   oblm, oblt, x, x, x )
        tlast = tjd
    end if
!         SELECT MEAN OR TRUE OBLIQUITY
    obl = oblm * radcon
    if ( icoord == 1 ) obl = oblt * radcon
end if

!     ROTATE ECLIPTIC POSITION VECTOR TO EQUATORIAL SYSTEM
pos2(1) =  pos1(1)
pos2(2) =  pos1(2) * dcos ( obl ) - pos1(3) * dsin ( obl )
pos2(3) =  pos1(2) * dsin ( obl ) + pos1(3) * dcos ( obl )

if ( tjd == 0.d0 ) then
!         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
    pos0(1) = pos2(1)
    pos0(2) = pos2(2)
    pos0(3) = pos2(3)
    call frame ( pos0, -1,   pos2 )
end if

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( pos2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS ICRS RIGHT ASCENSION AND DECLINATION
!  TO GALACTIC LONGITUDE AND LATITUDE.  IT USES THE TRANSFORMATION
!  GIVEN IN THE HIPPARCOS AND TYCHO CATALOGUES, VOL. 1,
!  SECTION 1.5.3.
!
!       RA     = ICRS RIGHT ASCENSION IN HOURS (IN)
!       DEC    = ICRS DECLINATION IN DEGREES (IN)
!       GLON   = GALACTIC LONGITUDE IN DEGREES (OUT)
!       GLAT   = GALACTIC LATITUDE IN DEGREES (OUT)

subroutine eqgal ( ra, dec,   glon, glat )

double precision ra,dec,glon,glat,pi,radcon,ag,r,d,xyproj,g, &
     pos1,pos2,dsin,dcos,dsqrt,datan2
dimension pos1(3), pos2(3), ag(3,3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     ROTATION MATRIX A_G FROM HIPPARCOS DOCUMENTATION EQ. 1.5.11
data ag / &
     -0.0548755604d0, +0.4941094279d0, -0.8676661490d0, &
     -0.8734370902d0, -0.4448296300d0, -0.1980763734d0, &
     -0.4838350155d0, +0.7469822445d0, +0.4559837762d0 /

!     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
r = ra * 15.d0 * radcon
d = dec * radcon
pos1(1) = dcos ( d ) * dcos ( r )
pos1(2) = dcos ( d ) * dsin ( r )
pos1(3) = dsin ( d )

!     ROTATE POSITION VECTOR TO GALACTIC SYSTEM, USING HIPPARCOS
!     DOCUMENTATION EQ. 1.5.13
pos2(1) = ag(1,1)*pos1(1) + ag(1,2)*pos1(2) + ag(1,3)*pos1(3)
pos2(2) = ag(2,1)*pos1(1) + ag(2,2)*pos1(2) + ag(2,3)*pos1(3)
pos2(3) = ag(3,1)*pos1(1) + ag(3,2)*pos1(2) + ag(3,3)*pos1(3)

!     DECOMPOSE GALACTIC VECTOR INTO LONGITUDE AND LATITUDE
xyproj = dsqrt ( pos2(1)**2 + pos2(2)**2 )
g = 0.d0
if ( xyproj > 0.d0 ) g = datan2 ( pos2(2), pos2(1) )
glon = g / radcon
if ( glon < 0.d0 ) glon = glon + 360.d0
g = datan2 ( pos2(3), xyproj )
glat = g / radcon

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec ( pos2 )

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS ANGULAR QUANTITIES RELATED TO A STAR'S
!  POSITION AND MOTION TO VECTORS.
!
!       RA     = RIGHT ASCENSION IN HOURS (IN)
!       DEC    = DECLINATION IN DEGREES (IN)
!       PMRA   = PROPER MOTION IN RA IN MILLIARCSECONDS PER YEAR
!                (IN)
!       PMDEC  = PROPER MOTION IN DEC IN MILLIARCSECONDS PER YEAR
!                (IN)
!       PARLLX = PARALLAX IN MILLIARCSECONDS (IN)
!       RV     = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
!       POS    = POSITION VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
!                WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
!                IN AU (OUT)
!       VEL    = VELOCITY VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
!                WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
!                IN AU/DAY (OUT)

subroutine vectrs (ra,dec,pmra,pmdec,parllx,rv,pos,vel)

double precision ra,dec,pmra,pmdec,parllx,rv,pos,vel, &
     pi,seccon,aukm,c,paralx,dist,r,d,cra,sra,cdc,sdc,k, &
     pmr,pmd,rvl,dcos,dsin
dimension pos(3), vel(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

data ntimes / 0 /

ntimes = ntimes + 1
if ( ntimes == 1 ) then
!         GET LENGTH OF AU IN KILOMETERS
    call astcon ( 'AU', 1.d-3,   aukm )
!         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
    call astcon ( 'C', 1.d-3,   c )
end if

!     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
!     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
paralx = parllx
if ( paralx <= 0.d0 ) paralx = 1.d-6

!     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
!     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
dist = 1.d0 / dsin ( paralx * 1.d-3 / seccon )
r = ra * 54000.d0 / seccon
d = dec * 3600.d0 / seccon
cra = dcos ( r )
sra = dsin ( r )
cdc = dcos ( d )
sdc = dsin ( d )
pos(1) = dist * cdc * cra
pos(2) = dist * cdc * sra
pos(3) = dist * sdc

!     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
!     LIGHT TRAVEL TIME TO STAR
k = 1.d0 / ( 1.d0 - rv / c )

!     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL COMPONENTS
!     OF MOTION WITH UNITS OF AU/DAY
pmr = pmra  / ( paralx * 365.25d0 ) * k
pmd = pmdec / ( paralx * 365.25d0 ) * k
rvl = rv * 86400.d0 / aukm          * k

!     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
vel(1) = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra
vel(2) =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra
vel(3) =               pmd * cdc       + rvl * sdc

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CONVERTS A VECTOR TO ANGULAR QUANTITIES.
!
!       POS = POSITION VECTOR, EQUATORIAL RECTANGULAR
!             COORDINATES (IN)
!       RA  = RIGHT ASCENSION IN HOURS (OUT)
!       DEC = DECLINATION IN DEGREES (OUT)

subroutine angles (pos,ra,dec)

double precision pos,ra,dec,pi,seccon,xyproj,r,d,dsqrt,datan2
dimension pos(3)

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

xyproj = dsqrt(pos(1)**2 + pos(2)**2)
r = 0.d0
if (xyproj>0.d0) r = datan2(pos(2),pos(1))
ra = r * seccon / 54000.d0
if (ra< 0.d0) ra = ra + 24.d0
if (ra>=24.d0) ra = ra - 24.d0
d = datan2(pos(3),xyproj)
dec = d * seccon / 3600.d0

!     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
call setvec (pos)

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE APPLIES PROPER MOTION, INCLUDING FORESHORTENING
!  EFFECTS, TO A STAR'S POSITION.
!
!       TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
!       POS1 = POSITION VECTOR OF STAR AT FIRST EPOCH (IN)
!       VEL1 = VELOCITY VECTOR OF STAR AT FIRST EPOCH (IN)
!       TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
!       POS2 = POSITION VECTOR OF STAR AT SECOND EPOCH (OUT)

subroutine propmo (tjd1,pos1,vel1,tjd2,pos2)

double precision tjd1,pos1,vel1,tjd2,pos2
dimension pos1(3), vel1(3), pos2(3)

do j=1,3
    pos2(j) = pos1(j) + vel1(j) * (tjd2 - tjd1)
end do

end subroutine propmo
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE MOVES THE ORIGIN OF COORDINATES FROM THE
!  BARYCENTER OF THE SOLAR SYSTEM TO THE OBSERVER (OR THE
!  GEOCENTER).  I.E., THIS SUBROUTINE ACCOUNTS FOR PARALLAX
!  (ANNUAL+GEOCENTRIC OR JUST ANNUAL).
!
!       POS1   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
!                ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
!                IN AU (IN)
!       PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                COMPONENTS IN AU (IN)
!       POS2   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
!                IN AU (OUT)
!       TLIGHT = LIGHT-TIME FROM STAR OR PLANET TO OBSERVER (OR THE
!                GEOCENTER) IN DAYS (OUT)
!
!  NOTE: STAR AND PLANET ARE USED GENERICALLY FOR BODIES OUTSIDE AND
!        INSIDE THE SOLAR SYSTEM, RESPECTIVELY.

subroutine geocen (pos1,pe,pos2,tlight)

double precision pos1,pe,pos2,tlight,c,dsqrt
dimension pos1(3), pe(3), pos2(3)
save

data ntimes / 0 /

ntimes = ntimes + 1
if (ntimes==1) then
    ! GET C, THE SPEED OF LIGHT IN AU/DAY
    call astcon ('C(AU/DAY)',1.d0,c)
end if

do j=1,3
    pos2(j) = pos1(j) - pe(j)
end do
tlight = dsqrt(pos2(1)**2 + pos2(2)**2 + pos2(3)**2) / c

end subroutine geocen
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE GEOCENTRIC POSITION AND VELOCITY
!  OF AN OBSERVER ON THE SURFACE OF THE EARTH OR ON A NEAR-EARTH
!  SPACECRAFT.  THE FINAL VECTORS ARE EXPRESSED IN THE GCRS.
!
!       TJD    = TT JULIAN DATE (IN)
!       LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
!                SET LOCATN=0 FOR OBSERVER AT GEOCENTER
!                SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
!                SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
!       OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
!                FOR LOCATN=0, THIS ARRAY NOT USED
!                FOR LOCATN=1,
!                OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
!                            (EAST +) IN DEGREES (IN)
!                OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
!                            (NORTH +) IN DEGREES (IN)
!                OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
!                            IN METERS (IN)
!                OBSERV(4) = VALUE OF DELTA-T IN SECONDS (IN)
!                            (DELTA-T=TT-UT1)
!                OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
!                OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
!                FOR LOCATN=2,
!                OBSERV(1) = GEOCENTRIC X IN KILOMETERS
!                OBSERV(2) = GEOCENTRIC Y IN KILOMETERS
!                OBSERV(3) = GEOCENTRIC Z IN KILOMETERS
!                OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND
!                OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND
!                OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
!                WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
!       POS    = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
!                AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
!                IN AU (OUT)
!       VEL    = VELOCITY VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
!                AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
!                IN AU/DAY (OUT)
!
!  NOTE 1: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
!  BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
!  DEFINED BY USER VIA CALL TO SETDT.
!
!  NOTE 2: THIS SUBROUTINE CALLS SUBROUTINE TERRA FOR AN OBSERVER
!  ON THE SURFACE OF THE EARTH.  TERRA NEGLECTS POLAR MOTION, AN
!  APPROXIMATION WHICH MAY YIELD UP TO 15 METERS ERROR IN POSITION
!  AND SEVERAL MILLIMETERS/SEC ERROR IN VELOCITY.

subroutine geopos (tjd,locatn,observ,pos,vel)

double precision tjd,observ,pos,vel,t0,tlast, &
     au,deltat,ttjd,tdbjd,ut1jd,st,gst,gmst,gast,eqeq,x, &
     pos1,vel1,pos2,vel2,pos3,vel3

dimension observ(6), pos(3), vel(3), &
     pos1(3), vel1(3), pos2(3), vel2(3), pos3(3), vel3(3)

save

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /,   gst / -99.d0 /,   ntimes / 0 /

ntimes = ntimes + 1
if ( ntimes == 1 ) then
!         GET AU, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
    call astcon ( 'AU', 1.d-3,   au )
end if

if ( locatn == 0 ) then
    pos(1) = 0.d0
    pos(2) = 0.d0
    pos(3) = 0.d0
    vel(1) = 0.d0
    vel(2) = 0.d0
    vel(3) = 0.d0
    return
end if

ttjd  = tjd
!     TDB IS APPROXIMATED BY TT
tdbjd = tjd

!     GET GEOCENTRIC POSITION AND VELOCITY VECTORS OF OBSERVER WRT
!     EQUATOR AND EQUINOX OF DATE

if ( locatn == 1 ) then

!         OBSERVER ON SURFACE OF EARTH

    if ( gst /= -99.d0 ) then
        ! TEMPORARY CODE TO USE SIDEREAL TIME PREVIOUSLY PROVIDED
        gast = gst
        gst = -99.d0
        ! END OF TEMPROARY CODE
    else

    !         GET DELTA-T VALUE
        if ( observ(4) /= 0.d0 ) then
            deltat = observ(4) / 86400.d0
        else
            call getdt ( deltat )
        end if

    !         USING DELTA-T VALUE, COMPUTE UT1 AND SIDEREAL TIME
        if ( ttjd == 0.d0 ) then
            ut1jd = tdbjd - deltat
        else
            ut1jd = ttjd - deltat
        end if
        if ( dabs ( ut1jd - tlast ) > 1.d-8 ) then
            call eqinox
            call sidtim ( ut1jd, 0.d0, 0,   gmst )
            call etilt ( tdbjd,   x, x, eqeq, x, x )
            call resume
            tlast = ut1jd
        end if
        gast = gmst + eqeq / 3600.d0

    end if

!         SUBROUTINE TERRA DOES THE HARD WORK, GIVEN SIDEREAL TIME
    call terra ( observ(1), observ(2), observ(3), gast, &
                 pos1, vel1 )

else if ( locatn == 2 ) then

!         OBSERVER ON NEAR-EARTH SPACECRAFT

!         CONVERT UNITS TO AU AND AU/DAY
    do j = 1, 3
        pos1(j) = observ(j)   / au
        vel1(j) = observ(j+3) / au * 86400.d0
    end do

end if

!     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
call nutate ( -tdbjd, pos1,   pos2 )
call preces ( tdbjd, pos2, t0,   pos3 )
call frame ( pos3, -1,   pos )

!     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
call nutate ( -tdbjd, vel1,   vel2 )
call preces ( tdbjd, vel2, t0,   vel3 )
call frame ( vel3, -1,   vel )

return


!     TEMPORARY CODE FOR COMPATIBILITY WITH OLD ROUTINES
entry placst ( st )
gst = st
return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE POSITION OF A SOLAR SYSTEM BODY,
!  AS ANTEDATED FOR LIGHT-TIME.
!
!       TJD    = TDB JULIAN DATE OF OBSERVATION (IN)
!       IDBODY = ID NUMBER OF BODY, USED IN CALLS TO SOLSYS (IN)
!       POSE   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
!       TLITE  = FIRST APPROXIMATION TO LIGHT-TIME, IN DAYS (IN)
!                (CAN BE SET TO 0.D0 IF UNKNOWN)
!       POS    = POSITION VECTOR OF BODY, WITH RESPECT TO ORIGIN AT
!                OBSERVER (OR THE GEOCENTER), REFERRED TO ICRS AXES,
!                COMPONENTS IN AU (OUT)
!       TLIGHT = FINAL LIGHT-TIME, IN DAYS (OUT)

subroutine littim (tjd,idbody,pose,tlite,pos,tlight)

double precision tjd,pose,tlite,pos,tlight,t0,t1,t2,t3,tol, &
     pos1,vel1,dint,dabs
logical split

dimension pose(3), pos(3), pos1(3), vel1(3)

save ntimes, split

data ntimes / 0 /,   split / .false. /

3 format ( ' LITTIM: PROBLEM WITH BODY NUMBER ', i3, ' AT JD ', &
     f10.1 )

ntimes = ntimes + 1

!     ON FIRST CALL, CHECK WHETHER SOLSYS SUPPORTS SPLIT JULIAN DATES
if ( ntimes == 1 ) split = idss ('JD') == 2

!     SET LIGHT-TIME CONVERGENCE TOLERANCE
tol = 1.d-9
if ( split .and. tlite < 0.01d0 ) tol = 1.d-12

!     IF SOLSYS SUPPORTS SPLIT JULIAN DATES, SPLIT THE JULIAN DATE
!     INTO WHOLE DAYS + FRACTION OF DAY
t0 = 0.d0
if ( split ) t0 = dint ( tjd )
t1 = tjd - t0
t2 = t1 - tlite
if ( split ) call solsys ( t0, idbody, 0,   pos1, vel1, ierr)
iter = 0

!     ITERATE TO OBTAIN CORRECT LIGHT-TIME (USUALLY CONVERGES RAPIDLY)
do
    call solsys ( t2, idbody, 0,   pos1, vel1, ierr )
    call geocen ( pos1, pose,   pos, tlight )
    if ( ierr /= 0 ) then
        write ( *, 3 ) idbody, t0 + t2
        return
    end if
    t3 = t1 - tlight
    if ( dabs ( t3 - t2 ) > tol ) then
        iter = iter + 1
        if ( iter > 10 ) then
            write ( *, 3 ) idbody, t0 + t3
            return
        end if
        t2 = t3
        cycle
    else
        exit
    end if
end do

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE RETURNS THE DIFFERENCE IN LIGHT-TIME, FOR A STAR,
!  BETWEEN THE BARYCENTER OF THE SOLAR SYSTEM AND THE OBSERVER (OR
!  THE GEOCENTER).
!
!       POS1   = POSITION VECTOR OF STAR, WITH RESPECT TO ORIGIN AT
!                SOLAR SYSTEM BARYCENTER (IN)
!       PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                COMPONENTS IN AU (IN)
!       DIFLT  = DIFFERENCE IN LIGHT TIME, IN THE SENSE STAR TO
!                BARYCENTER MINUS STAR TO EARTH, IN DAYS (OUT)
!
!  -OR-
!
!  THIS SUBROUTINE RETURNS THE LIGHT-TIME FROM THE OBSERVER (OR THE
!  GEOCENTER) TO A POINT ON A LIGHT RAY THAT IS CLOSEST TO A
!  SPECIFIC SOLAR SYSTEM BODY.
!
!       POS1   = POSITION VECTOR TOWARD OBSERVED OBJECT, WITH RESPECT
!                TO ORIGIN AT OBSERVER (OR THE GEOCENTER) (IN)
!       PE     = POSITION VECTOR OF SOLAR SYSTEM BODY, WITH RESPECT
!                TO ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
!                IN AU (IN)
!       DIFLT  = LIGHT TIME TO POINT ON LINE DEFINED BY POS1 THAT IS
!                CLOSEST TO SOLAR SYSTEM BODY (POSITIVE IF LIGHT
!                PASSES BODY BEFORE HITTING OBSERVER, I.E., IF
!                POS1 IS WITHIN 90 DEGREES OF PE)(OUT)

subroutine dlight (pos1,pe,diflt)

double precision pos1,pe,diflt,c,dis,u1,dsqrt
dimension pos1(3), pe(3), u1(3)
save

data ntimes / 0 /

ntimes = ntimes + 1
if (ntimes==1) then
!         GET C, THE SPEED OF LIGHT IN AU/DAY
    call astcon ('C(AU/DAY)',1.d0,c)
end if

!     FROM POS1, FORM UNIT VECTOR U1 IN DIRECTION OF STAR OR
!     LIGHT SOURCE
dis = dsqrt ( pos1(1)**2 + pos1(2)**2 + pos1(3)**2 )
do j=1,3
    u1(j) = pos1(j) / dis
end do

!     LIGHT-TIME RETURNED IS THE PROJECTION OF VECTOR PE ONTO THE UNIT
!     VECTOR U1 (FORMED FROM POS1), DIVIDED BY THE SPEED OF LIGHT
diflt = ( pe(1)*u1(1) + pe(2)*u1(2) + pe(3)*u1(3) ) / c

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE TOTAL GRAVITATIONAL DEFLECTION OF
!  LIGHT FOR THE OBSERVED OBJECT DUE TO THE MAJOR GRAVITATING BODIES
!  IN THE SOLAR SYSTEM.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY
!  WITHIN THE SOLAR SYSTEM AS WELL AS FOR A STAR.  SEE KLIONER
!  (2003), ASTRONOMICAL JOURNAL 125, 1580-1597, SECTION 6.
!
!       TJD    = TDB JULIAN DATE OF OBSERVATION
!       LOC    = CODE FOR LOCATION OF OBSERVER, DETERMINING
!                WHETHER THE GRAVITATIONAL DEFLECTION DUE TO THE
!                EARTH ITSELF IS APPLIED (IN)
!                SET LOC=0 FOR NO EARTH DEFLECTION (NORMALLY MEANS
!                          OBSERVER IS AT GEOCENTER)
!                SET LOC=1 TO ADD IN EARTH DEFLECTION (NORMALLY
!                          MEANS OBSERVER IS ON OR ABOVE SURFACE
!                          OF EARTH, INCLUDING EARTH ORBIT)
!       POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
!                TO ICRS AXES, COMPONENTS IN AU (IN)
!       POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
!       POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
!                TO ICRS AXES, CORRECTED FOR GRAVITATIONAL
!                DEFLECTION, COMPONENTS IN AU (OUT)

subroutine grvdef (tjd,loc,pos1,pobs,pos2)

double precision tjd,pos1,pobs,pos2,c,rmass,rmasse,pbody,vbody, &
     pbodyo,x,tlt,dlt,tclose,dsqrt
character*3 name
dimension pos1(3), pobs(3), pos2(3), name(10), id(10), rmass(10), &
     pbody(3), vbody(3), pbodyo(3)
save

!     THE FOLLOWING LIST OF NAMES IDENTIFIES WHICH GRAVITATING BODIES
!     (ASIDE FROM THE EARTH) ARE POTENTIALLY USED -- LIST IS TAKEN FROM
!     KLIONER'S TABLE 1, THE ORDER BASED ON AREA OF SKY AFFECTED (COL 2)
data name / 'SUN', 'JUP', 'SAT', 'MOO', 'VEN', 'URA', 'NEP', &
     3*'   ' /
!     CHANGE VALUE OF NBODY TO INCLUDE OR EXCLUDE GRAVITATING BODIES
!     (NBODY=0 MEANS NO DEFLECTION CALCULATED, NBODY=1 MEANS SUN ONLY,
!     NBODY=2 MEANS SUN + JUPITER, ETC.)
data nbody / 3 /

data ntimes / 0 /

ntimes = ntimes + 1
if ( ntimes == 1 ) then
!         GET C, THE SPEED OF LIGHT IN AU/DAY
    call astcon ( 'C(AU/DAY)', 1.d0, c )
!         GET ID NUMBERS AND RECIPROCAL MASSES OF GRAVITATING BODIES
    do i = 1, nbody
        id(i) = idss ( name(i) )
        call astcon ( 'MASS_'//name(i), 1.d0,   rmass(i) )
    end do
    ide = idss ( 'EARTH' )
    call astcon ( 'MASS_EARTH', 1.d0,   rmasse )
end if

!     INITIALIZE OUTPUT VECTOR OF OBSERVED OBJECT TO EQUAL INPUT VECTOR
do j = 1, 3
    pos2(j) = pos1(j)
end do
!     OPTION FOR NO DEFLECTION
if ( nbody <= 0 ) return

!     COMPUTE LIGHT-TIME TO OBSERVED OBJECT
tlt = dsqrt ( pos1(1)**2 + pos1(2)**2 + pos1(3)**2 ) / c

!     CYCLE THROUGH GRAVITATING BODIES
do i = 1, nbody

    if ( id(i) == -9999 ) cycle

!         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME TJD
    call solsys ( tjd, id(i), 0,   pbody, vbody, ierr )

!         GET POSITION OF GRAVITATING BODY WRT OBSERVER AT TIME TJD
    call geocen ( pbody, pobs,   pbodyo, x )

!         COMPUTE LIGHT-TIME FROM POINT ON INCOMING LIGHT RAY THAT
!         IS CLOSEST TO GRAVITATING BODY
    call dlight ( pos2, pbodyo,   dlt )

!         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME
!         WHEN INCOMING PHOTONS WERE CLOSEST TO IT
    tclose = tjd
    if ( dlt > 0.d0 ) tclose = tjd - dlt
    if ( tlt < dlt  ) tclose = tjd - tlt
    call solsys ( tclose, id(i), 0,   pbody, vbody, ierr )

!         COMPUTE DEFLECTION DUE TO GRAVITATING BODY
    call grvd ( pos2, pobs, pbody, rmass(i),   pos2 )

end do

!     IF OBSERVER IS NOT AT GEOCENTER, ADD IN DEFLECTION DUE TO EARTH
if ( loc /= 0 ) then

!         GET POSITION OF EARTH WRT SS BARYCENTER AT TIME TJD
    call solsys ( tjd, ide, 0,   pbody, vbody, ierr )

!         COMPUTE DEFLECTION DUE TO EARTH
    call grvd ( pos2, pobs, pbody, rmasse,   pos2 )

end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CORRECTS POSITION VECTOR FOR THE DEFLECTION
!  OF LIGHT IN THE GRAVITATIONAL FIELD OF AN ARBITRARY BODY.  ADAPTED
!  FROM MURRAY (1981) MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
!  SEE ALSO FORMULAE IN SECTION B OF THE ASTRONOMICAL ALMANAC, OR
!  KAPLAN ET AL. (1989) ASTRONOMICAL JOURNAL 97, 1197-1210, SECTION
!  III F.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY WITHIN THE
!  SOLAR SYSTEM AS WELL AS FOR A STAR.
!
!       POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
!                IN AU (IN)
!       POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                COMPONENTS IN AU (IN)
!       PBODY  = POSITION VECTOR OF GRAVITATING BODY, WITH RESPECT TO
!                ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
!                IN AU (IN)
!       RMASS  = RECIPROCAL MASS OF GRAVITATING BODY IN SOLAR MASS
!                UNITS, THAT IS, SUN MASS / BODY MASS (IN)
!       POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED FOR
!                GRAVITATIONAL DEFLECTION, COMPONENTS IN AU (OUT)

subroutine grvd (pos1,pobs,pbody,rmass,pos2)

double precision pos1,pobs,pbody,rmass,pos2,c,mau,gs,pq,pe, &
     pmag,emag,qmag,phat,ehat,qhat,pdotq,edotp,qdote, &
     fac1,fac2,p2j,dabs,dsqrt
dimension pos1(3), pobs(3), pbody(3), pos2(3), pq(3), pe(3), &
     phat(3), ehat(3), qhat(3)
save

data ntimes / 0 /

ntimes = ntimes + 1
if (ntimes==1) then
!         GET C, THE SPEED OF LIGHT IN METERS/SECOND
    call astcon ( 'C', 1.d0, c )
!         GET MAU, THE LENGTH OF THE AU IN METERS
    call astcon ( 'AU', 1.d0, mau )
!         GET GS, THE HELIOCENTRIC GRAVITATIONAL CONSTANT
    call astcon ( 'GS', 1.d0, gs )
end if

!     CONSTRUCT VECTOR PQ FROM GRAVITATING BODY TO OBSERVED OBJECT AND
!     CONSTRUCT VECTOR PE FROM GRAVITATING BODY TO OBSERVER
do j=1,3
    pq(j) = pobs(j) + pos1(j) - pbody(j)
    pe(j) = pobs(j) - pbody(j)
end do

!     COMPUTE VECTOR MAGNITUDES AND UNIT VECTORS
pmag = dsqrt (pos1(1)**2 + pos1(2)**2 + pos1(3)**2)
emag = dsqrt (  pe(1)**2 +   pe(2)**2 +   pe(3)**2)
qmag = dsqrt (  pq(1)**2 +   pq(2)**2 +   pq(3)**2)
do j = 1, 3
    phat(j) = pos1(j) / pmag
    ehat(j) =   pe(j) / emag
    qhat(j) =   pq(j) / qmag
end do

!     COMPUTE DOT PRODUCTS OF VECTORS
pdotq = phat(1)*qhat(1) + phat(2)*qhat(2) + phat(3)*qhat(3)
edotp = ehat(1)*phat(1) + ehat(2)*phat(2) + ehat(3)*phat(3)
qdote = qhat(1)*ehat(1) + qhat(2)*ehat(2) + qhat(3)*ehat(3)

!     IF GRAVITATING BODY IS OBSERVED OBJECT, OR IS ON A STRAIGHT LINE
!     TOWARD OR AWAY FROM OBSERVED OBJECT TO WITHIN 1 ARCSEC,
!     DEFLECTION IS SET TO ZERO
if ( dabs ( edotp ) > 0.99999999999d0 ) then
    do j=1,3
        pos2(j) = pos1(j)
    end do
    return
end if

!     COMPUTE SCALAR FACTORS
fac1 = 2.0d0 * gs / (c * c * emag * mau * rmass)
fac2 = 1.0d0 + qdote

!     CONSTRUCT CORRECTED POSITION VECTOR POS2
do j = 1, 3
    p2j = phat(j) + fac1 * (pdotq*ehat(j) - edotp*qhat(j)) / fac2
    pos2(j) = p2j * pmag
end do

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CORRECTS POSITION VECTOR FOR ABERRATION OF LIGHT.
!  ALGORITHM INCLUDES RELATIVISTIC TERMS.  ADAPTED FROM MURRAY (1981)
!  MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
!
!       POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH REESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
!                IN AU (IN)
!       VE     = VELOCITY VECTOR OF OBSERVER (OR THE GEOCENTER),
!                WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
!                COMPONENTS IN AU/DAY (IN)
!       TLIGHT = LIGHT TIME FROM BODY TO OBSERVER (OR THE GEOCENTER)
!                IN DAYS (IN)
!                IF TLIGHT = 0.D0, THIS SUBROUTINE WILL COMPUTE
!       POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED
!                FOR ABERRATION, COMPONENTS IN AU (OUT)

subroutine aberat (pos1,ve,tlight,pos2)

double precision pos1,ve,tlight,pos2,c,tl,p1mag,vemag, &
     beta,dot,cosd,gammai,p,q,r,dsqrt
dimension pos1(3), ve(3), pos2(3)
save

data ntimes / 0 /

ntimes = ntimes + 1
if (ntimes==1) then
!         GET C, THE SPEED OF LIGHT IN AU/DAY
    call astcon ('C(AU/DAY)',1.d0,c)
end if

tl = tlight
p1mag = tl * c
if (tl==0.d0) then
    p1mag = dsqrt(pos1(1)**2 + pos1(2)**2 + pos1(3)**2)
    tl = p1mag / c
end if
vemag = dsqrt(ve(1)**2 + ve(2)**2 + ve(3)**2)
beta = vemag / c
dot = pos1(1)*ve(1) + pos1(2)*ve(2) + pos1(3)*ve(3)
cosd = dot / (p1mag * vemag)
gammai = dsqrt(1.d0 - beta**2)
p = beta * cosd
q = (1.d0 + p / (1.d0 + gammai)) * tl
r = 1.d0 + p

do j=1,3
    pos2(j) = (gammai * pos1(j) + q * ve(j)) / r
end do

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE PREDICTS THE RADIAL VELOCITY OF THE OBSERVED
!  OBJECT AS IT WOULD BE MEASURED BY SPECTROSCOPIC MEANS.  RADIAL
!  VELOCITY IS HERE DEFINED AS THE RADIAL VELOCITY MEASURE (Z)
!  TIMES THE SPEED OF LIGHT.  FOR A SOLAR SYSTEM BODY, IT APPLIES
!  TO A FICTITIOUS EMITTER AT THE CENTER OF THE OBSERVED OBJECT,
!  ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT), AND DOES NOT
!  IN GENERAL APPLY TO REFLECTED LIGHT.  FOR STARS, IT INCLUDES
!  ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT, CONTAINED
!  IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE, A SCALAR
!  DERIVED FROM SPECTROSCOPY.  NEARBY STARS WITH A KNOWN KINEMATIC
!  VELOCITY VECTOR (OBTAINED INDEPENDENTLY OF SPECTROSCOPY) CAN BE
!  TREATED LIKE SOLAR SYSTEM OBJECTS.  SEE LINDEGREN & DRAVINS
!  (2003), ASTRONOMY & ASTROPHYSICS 401, 1185-1201.
!
!       POS    = GEOMETRIC POSITION VECTOR OF OBJECT WITH RESPECT TO
!                OBSERVER, CORRECTED FOR LIGHT-TIME, IN AU (IN)
!       VEL    = VELOCITY VECTOR OF OBJECT WITH RESPECT TO SOLAR
!                SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
!       VELOBS = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO SOLAR
!                SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
!       STAR   = 3-ELEMENT ARRAY OF CATALOG DATA FOR A STAR, TO BE
!                NON-ZERO IF OBSERVED OBJECT IS A STAR FOR WHICH THE
!                CATALOG RADIAL VELOCITY IS CONSISTENT WITH
!                THE IAU DEFINITION OF BARYCENTRIC RADIAL VELOCITY
!                MEASURE (OTHERWISE ALL ELEMENTS SHOULD BE SET TO
!                0.D0 EXACTLY) (IN)
!                STAR(1) = CATALOG RA IN HOURS
!                STAR(2) = CATALOG DEC IN DEGREES
!                STAR(3) = Z*C, THE CATALOG BARYCENTRIC RADIAL
!                          VELOCITY MEASURE TIMES THE SPEED OF LIGHT,
!                          IN KILOMETERS/SECOND
!                ALL THREE DATA ELEMENTS MUST APPLY TO THE SAME
!                EPOCH (USUALLY J2000.0 = JD 2451545.0 TT)
!       DIST   = 3-ELEMENT ARRAY OF DISTANCES IN AU (IN)
!                DIST(1) = DISTANCE OF OBSERVER FROM THE GEOCENTER
!                DIST(2) = DISTANCE OF OBSERVER FROM THE SUN
!                DIST(3) = DISTANCE OF OBJECT FROM THE SUN
!       RV     = THE OBSERVED RADIAL VELOCITY MEASURE TIMES
!                THE SPEED OF LIGHT, IN KILOMETERS/SECOND (OUT)
!
!  NOTE 1:  ALL THE INPUT ARGUMENTS ARE BCRS QUANTITIES, EXPRESSED
!  WITH RESPECT TO THE ICRS AXES.  VEL AND VELOBS ARE KINEMATIC
!  VELOCITIES -- DERIVED FROM GEOMETRY OR DYNAMICS, NOT SPECTROSCOPY.
!
!  NOTE 2:  IF ANY ELEMENT OF ARRAY STAR IS NON-ZERO, THE ALGORITHM
!  USED WILL BE CONSISTENT WITH THE IAU DEFINITION OF STELLAR
!  RADIAL VELOCITY, SPECIFICALLY, THE BARYCENTRIC RADIAL VELOCITY
!  MEASURE, WHICH IS DERIVED FROM SPECTROSCOPY.  IN THAT CASE,
!  THE VECTOR VEL CAN BE VERY APPROXIMATE -- OR, FOR DISTANT STARS
!  OR GALAXIES, ZERO -- AS IT WILL BE USED ONLY FOR A SMALL GEOMETRIC
!  CORRECTION THAT IS PROPORTIONAL TO PROPER MOTION.
!
!  NOTE 3:  ANY OF THE DISTANCES IN ARRAY DIST CAN BE SET TO ZERO
!  (0.D0) IF THE CORRESPONDING GENERAL RELATIVISTIC GRAVITATIONAL
!  POTENTIAL TERM IS NOT TO BE EVALUATED.  THESE TERMS
!  GENERALLY ARE IMPORTANT ONLY AT THE METER/SECOND LEVEL.  IF
!  THE FIRST TWO DISTANCES ARE BOTH ZERO, AN AVERAGE VALUE
!  WILL BE USED FOR THE RELATIVISTIC TERM FOR THE OBSERVER,
!  APPROPRIATE FOR AN OBSERVER ON THE SURFACE OF THE EARTH.  THE
!  THIRD DISTANCE, IF GIVEN, IS USED ONLY FOR SOLAR SYSTEM OBJECTS.

subroutine radvl ( pos, vel, velobs, star, dist,   rv )

double precision pos,vel,velobs,star,dist,rv,pi,radcon, &
     au,c,gs,ge,c2,toms,toms2, &
     posmag,uk,v2,vo2,r,phigeo,phisun,rel,zc,ra,dc,du, &
     zb1,kvobs,kv,zobs1,dsqrt,dcos,dsin

logical dostar

dimension pos(3), vel(3), velobs(3), star(3), dist(3), uk(3), &
     du(3)

save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

data ntimes / 0 /

ntimes = ntimes + 1

if ( ntimes == 1 ) then
!         GET AU, LENGTH OF ASTRONOMICAL UNIT IN METERS
    call astcon ( 'AU', 1.d0,   au )
!         GET C, THE SPEED OF LIGHT IN METERS/SECOND
    call astcon ( 'C', 1.d0,   c )
!         GET GS, HELIOCENTRIC GRAVITATIONAL CONSTANT
    call astcon ( 'GS', 1.d0,   gs )
!         GET GE, GEOCENTRIC GRAVITATIONAL CONSTANT
    call astcon ( 'GE', 1.d0,   ge )
!         (GS AND GE ARE IN METERS**3/SECOND**2)
    c2 = c**2
    toms = au / 86400.d0
    toms2 = toms**2
end if

rv = 0.d0

!     COMPUTE LENGTH OF POSITION VECTOR = DISTANCE TO OBJECT IN AU
posmag = dsqrt ( pos(1)**2 + pos(2)**2 + pos(3)**2 )
if ( posmag < 1.d-8 ) return

!     DETERMINE HOW OBJECT IS TO BE PROCESSED
dostar = star(1) /= 0.d0 .or. &
         star(2) /= 0.d0 .or. &
         star(3) /= 0.d0

!     COMPUTE UNIT VECTOR TOWARD OBJECT
do j = 1, 3
    uk(j) = pos(j) / posmag
end do

!     COMPUTE VELOCITY-SQUARED FACTORS
v2  = ( vel(1)   **2 + vel(2)   **2 + vel(3)   **2 ) * toms2
vo2 = ( velobs(1)**2 + velobs(2)**2 + velobs(3)**2 ) * toms2

!     COMPUTE GEOPOTENTIAL AT OBSERVER, UNLESS OBSERVER IS GEOCENTRIC
r = dist(1) * au
phigeo = 0.d0
if ( r > 1.d6 ) phigeo = ge / r

!     COMPUTE SOLAR POTENTIAL AT OBSERVER
r = dist(2) * au
phisun = 0.d0
if ( r > 1.d8 ) phisun = gs / r

!     COMPUTE RELATIVISTIC POTENTIAL AND VELOCITY FACTOR FOR OBSERVER
if ( dist(1) /= 0.d0 .or. dist(2) /= 0.d0 ) then
!         LINDEGREN & DRAVINS EQ. (41), SECOND FACTOR IN PARENTHESES
    rel = 1.d0 - ( phigeo + phisun ) / c2 - 0.5d0 * vo2 / c2
else
!         LINDEGREN & DRAVINS EQ. (42), INVERSE
    rel = 1.d0 - 1.550d-8
end if

if ( dostar ) then

!         FOR STARS, UPDATE BARYCENTRIC RADIAL VELOCITY MEASURE FOR
!         CHANGE IN VIEW ANGLE
    ra = star(1) * 15.d0 * radcon
    dc = star(2) * radcon
    du(1) = uk(1) - ( dcos ( dc ) * dcos ( ra ) )
    du(2) = uk(2) - ( dcos ( dc ) * dsin ( ra ) )
    du(3) = uk(3) - ( dsin ( dc )               )
    zc = star(3) * 1.d3 + &
       ( vel(1) * du(1) + vel(2) * du(2) + vel(3) * du(3) ) * toms      !

!         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A STAR (INVERSE OF
!         LINDEGREN & DRAVINS EQ. (41))
    zb1 = 1.d0 + zc / c
    kvobs = ( uk(1) * velobs(1) &
            + uk(2) * velobs(2) &
            + uk(3) * velobs(3) ) * toms
    zobs1 = zb1 * rel / ( 1.d0 + kvobs / c )

else

!         COMPUTE SOLAR POTENTIAL AT OBJECT, IF WITHIN SOLAR SYSTEM
    r = dist(3) * au
    phisun = 0.d0
    if ( r > 1.d8 .and. r < 1.d16 ) phisun = gs / r

!         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A PLANET OR OTHER
!         OBJECT -- INCLUDING A NEARBY STAR -- WHERE KINEMATIC
!         BARYCENTRIC VELOCITY VECTOR IS KNOWN AND GRAVITATIONAL
!         RED SHIFT IS NEGLIGIBLE (LINDEGREN & DRAVINS EQ. (40),
!         APPLIED AS PER S. KLIONER PRIVATE COMMUNICATION (2006))
    kv = ( uk(1) * vel(1) &
         + uk(2) * vel(2) &
         + uk(3) * vel(3) ) * toms
    zb1 = ( 1.d0 + kv / c ) / &
          ( 1.d0 - phisun / c2 - 0.5d0 * v2 / c2 )
    kvobs = ( uk(1) * velobs(1) &
            + uk(2) * velobs(2) &
            + uk(3) * velobs(3) ) * toms
    zobs1 = zb1 * rel / ( 1.d0 + kvobs / c )

end if

!     CONVERT OBSERVED RADIAL VELOCITY MEASURE TO KILOMETERS/SECOND
rv = ( zobs1 - 1.d0 ) * c / 1000.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE PRECESSES EQUATORIAL RECTANGULAR COORDINATES FROM
!  ONE EPOCH TO ANOTHER.  THE COORDINATES ARE REFERRED TO THE MEAN
!  DYNAMICAL EQUATOR AND EQUINOX OF THE TWO RESPECTIVE EPOCHS.  SEE
!  EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, PP. 103-104,
!  AND CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412,
!  567-586.
!
!       TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
!       POS1 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
!              EQUINOX OF FIRST EPOCH (IN)
!       TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
!       POS2 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
!              EQUINOX OF SECOND EPOCH (OUT)
!
!  NOTE:  EITHER TJD1 OR TJD2 MUST BE 2451545.0 (J2000.0) TDB.

subroutine preces (tjd1,pos1,tjd2,pos2)

double precision tjd1,tjd2,pos1,pos2,pi,seccon,t0,tlast,t, &
     eps0,psia,omegaa,chia,sa,ca,sb,cb,sc,cc,sd,cd, &
     xx,yx,zx,xy,yy,zy,xz,yz,zz,dabs,dcos,dsin
dimension pos1(3), pos2(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /

!     INITIALIZE PRECESSION ROTATION MATRIX AS IDENTITY MATRIX
data xx, xy, xz / 1.d0, 0.d0, 0.d0 /
data yx, yy, yz / 0.d0, 1.d0, 0.d0 /
data zx, zy, zz / 0.d0, 0.d0, 1.d0 /

3 format ( ' PRECES ERROR: PRECESSION FROM JD ', f10.1, ' TO ', &
     f10.1, ' NOT TO/FROM J2000' )

if ( tjd1 /= t0 .and. tjd2 /= t0 ) then
    write ( *, 3 ) tjd1, tjd2
    return
end if

!     T IS TIME IN TDB CENTURIES BETWEEN THE TWO EPOCHS
t = ( tjd2 - tjd1 ) / 36525.d0
if ( tjd2 == t0 ) t = -t
if ( dabs ( t - tlast ) >= 1.d-15 ) then

    !     NUMERICAL COEFFICIENTS OF PSI_A, OMEGA_A, AND CHI_A, ALONG WITH
    !     EPSILON_0, THE OBLIQUITY AT J2000.0, ARE 4-ANGLE FORMULATION
    !     FROM CAPITAINE ET AL. (2003), EQS. (4), (37), & (39)
    eps0   = 84381.406d0
    psia   = ( ( ( ( -    0.0000000951d0   * t &
                    +    0.000132851d0  ) * t &
                    -    0.00114045d0   ) * t &
                    -    1.0790069d0    ) * t &
                    + 5038.481507d0     ) * t
    omegaa = ( ( ( ( +    0.0000003337d0   * t &
                    -    0.000000467d0  ) * t &
                    -    0.00772503d0   ) * t &
                    +    0.0512623d0    ) * t &
                    -    0.025754d0     ) * t + eps0
    chia   = ( ( ( ( -    0.0000000560d0   * t &
                    +    0.000170663d0  ) * t &
                    -    0.00121197d0   ) * t &
                    -    2.3814292d0    ) * t &
                    +   10.556403d0     ) * t
    eps0 = eps0 / seccon
    psia = psia / seccon
    omegaa = omegaa / seccon
    chia = chia / seccon
    sa = dsin ( eps0 )
    ca = dcos ( eps0 )
    sb = dsin ( -psia )
    cb = dcos ( -psia )
    sc = dsin ( -omegaa )
    cc = dcos ( -omegaa )
    sd = dsin ( chia )
    cd = dcos ( chia )

    !     COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
    !     EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
    xx =  cd * cb - sb * sd * cc
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc
    xy = -sd * cb - sb * cd * cc
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc
    xz =  sb * sc
    yz = -sc * cb * ca - sa * cc
    zz = -sc * cb * sa + cc * ca

    tlast = t

end if

if ( tjd2 == t0 ) then
    !     PERFORM ROTATION FROM EPOCH TO J2000.0
    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3)
    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3)
    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3)
else
    !     PERFORM ROTATION FROM J2000.0 TO EPOCH
    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3)
    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3)
    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3)
end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE NUTATES EQUATORIAL RECTANGULAR COORDINATES FROM
!  THE MEAN DYNAMICAL EQUATOR AND EQUINOX OF EPOCH TO THE TRUE
!  EQUATOR AND EQUINOX OF EPOCH.  SEE EXPLANATORY SUPPLEMENT TO THE
!  ASTRONOMICAL ALMANAC, PP. 114-115.
!
!       TJD    = TDB JULIAN DATE OF EPOCH (IN)
!       POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
!                EQUINOX OF EPOCH (IN)
!       POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO TRUE EQUATOR AND EQUINOX
!                OF EPOCH (OUT)
!
!  NOTE:  IF TJD IS NEGATIVE, INVERSE NUTATION (TRUE TO MEAN)
!  IS APPLIED.

subroutine nutate (tjd,pos1,pos2)
double precision tjd,pos1,pos2,tjd1,pi,seccon,oblm,oblt,eqeq, &
     dpsi,deps,cobm,sobm,cobt,sobt,cpsi,spsi, &
     xx,yx,zx,xy,yy,zy,xz,yz,zz,dabs,dcos,dsin
dimension pos1(3), pos2(3)

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

tjd1 = dabs(tjd)

call etilt ( tjd1,   oblm, oblt, eqeq, dpsi, deps )
oblm = oblm * 3600.d0 / seccon
oblt = oblt * 3600.d0 / seccon
dpsi = dpsi / seccon
deps = deps / seccon
cobm = dcos ( oblm )
sobm = dsin ( oblm )
cobt = dcos ( oblt )
sobt = dsin ( oblt )
cpsi = dcos ( dpsi )
spsi = dsin ( dpsi )

!     COMPUTE ELEMENTS OF NUTATION ROTATION MATRIX
xx =  cpsi
yx = -spsi * cobm
zx = -spsi * sobm
xy =  spsi * cobt
yy =  cpsi * cobm * cobt + sobm * sobt
zy =  cpsi * sobm * cobt - cobm * sobt
xz =  spsi * sobt
yz =  cpsi * cobm * sobt - sobm * cobt
zz =  cpsi * sobm * sobt + cobm * cobt

if ( tjd < 0.d0 ) then
    !     PERFORM ROTATION FROM TRUE TO MEAN
    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3)
    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3)
    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3)
else
    !     PERFORM ROTATION FROM MEAN TO TRUE
    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3)
    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3)
    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3)
end if


end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE TRANSFORMS A VECTOR FROM ONE COORDINATE SYSTEM
!  TO ANOTHER WITH SAME ORIGIN AND AXES ROTATED ABOUT THE
!  Z AXIS.
!
!       ANGL   = ANGLE OF COORDINATE SYSTEM ROTATION, POSITIVE
!                COUNTERCLOCKWISE WHEN VIEWED FROM +Z,
!                IN DEGREES (IN)
!       POS1   = POSITION VECTOR (IN)
!       POS2   = POSITION VECTOR EXPRESSED IN NEW COORDINATE
!                SYSTEM ROTATED ABOUT Z BY ANGLE ANG (OUT)

subroutine spin (angl,pos1,pos2)

double precision angl,pos1,pos2,pi,alast,ang,cosang,sinang, &
     xx,yx,zx,xy,yy,zy,xz,yz,zz,dabs,dcos,dsin
dimension pos1(3), pos2(3)
save

parameter ( pi     = 3.14159265358979324d0 )

data alast / -999.d0 /

if ( dabs ( angl - alast ) > 1.d-12 ) then

    ang = angl / 180.d0 * pi
    cosang = dcos ( ang )
    sinang = dsin ( ang )

!         ROTATION MATRIX FOLLOWS
    xx =  cosang
    yx =  sinang
    zx =  0.d0
    xy = -sinang
    yy =  cosang
    zy =  0.d0
    xz =  0.d0
    yz =  0.d0
    zz =  1.d0

    alast = angl

end if

!     PERFORM ROTATION
pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3)
pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3)
pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3)

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE CORRECTS A VECTOR IN THE ITRS (A ROTATING EARTH-
!  FIXED SYSTEM) FOR POLAR MOTION, AND ALSO CORRECTS THE LONGITUDE
!  ORIGIN (BY A TINY AMOUNT) TO THE TERRESTRIAL INTERMEDIATE ORIGIN
!  (TIO).  THE ITRS VECTOR IS THEREBY TRANSFORMED TO THE TERRESTRIAL
!  INTERMEDIATE SYSTEM, BASED ON THE TRUE (ROTATIONAL) EQUATOR AND
!  THE TERRESTRIAL INTERMEDIATE ORIGIN (TIO).  SINCE THE TRUE EQUATOR
!  IS THE PLANE ORTHOGONAL TO THE DIRECTION OF THE CELESTIAL
!  INTERMEDIATE POLE (CIP), THE COMPONENTS OF THE OUTPUT VECTOR ARE
!  REFERRED TO Z AND X AXES TOWARD THE CIP AND TIO, RESPECTIVELY.
!
!       TJD    = TT OR UT1 JULIAN DATE (IN)
!       XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
!                ARCSECONDS (IN)
!       YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
!                INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
!                ARCSECONDS (IN)
!       POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO ITRS AXES (IN)
!       POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO TRUE EQUATOR AND TIO (OUT)
!
!  NOTE 1:  IF TJD IS NEGATIVE, THE INVERSE TRANSFORMATION (TERRESTRIAL
!  INTERMEDIATE SYSTEM TO ITRS) IS APPLIED.
!
!  NOTE 2: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
!  THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
!  IERS CONVENTIONS.

subroutine wobble (tjd,xp,yp,pos1,pos2)

double precision tjd,xp,yp,pos1,pos2,pi,seccon,t0,t,xpole,ypole, &
     sprime,tiolon,sinx,cosx,siny,cosy,sinl,cosl, &
     xx,yx,zx,xy,yy,zy,xz,yz,zz,dabs,dsin,dcos
dimension pos1(3), pos2(3)

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

!     T0 = TT JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.0d0 /

xpole = xp / seccon
ypole = yp / seccon

t = ( dabs(tjd) - t0 ) / 36525.d0

!     COMPUTE APPROXIMATE LONGITUDE OF TIO, USING EQ. (10) OF
!     LAMBERT & BIZOUARD (2002), ASTRONOMY AND ASTROPHYSICS 394,
!     317-321
sprime = -47.0d-6 * t
tiolon = -sprime / seccon
!     NOTE THAT TIOLON, THE LONGITUDE CORRECTION, IS NEGLIGIBLE FOR
!     MOST ASTRONOMICAL PURPOSES

!     COMPUTE ELEMENTS OF ROTATION MATRIX
!     EQUIVALENT TO R3(-S')R2(X)R1(Y) AS PER IERS CONVENTIONS (2003)
sinx = dsin ( xpole )
cosx = dcos ( xpole )
siny = dsin ( ypole )
cosy = dcos ( ypole )
sinl = dsin ( tiolon )
cosl = dcos ( tiolon )
xx =  cosx * cosl
yx =  sinx * siny * cosl + cosy * sinl
zx = -sinx * cosy * cosl + siny * sinl
xy = -cosx * sinl
yy =  sinx * siny * sinl + cosy * cosl
zy =  sinx * cosy * sinl + siny * cosl
xz =  sinx
yz = -cosx * siny
zz =  cosx * cosy

if ( tjd < 0.d0 ) then
    !     PERFORM ROTATION FROM TERRESTRIAL INTERMEDIATE SYSTEM TO ITRS
    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3)
    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3)
    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3)
else
    !     PERFORM ROTATION FROM ITRS TO TERRESTRIAL INTERMEDIATE SYSTEM
    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3)
    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3)
    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3)
end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE TRANSFORMS A VECTOR FROM THE DYNAMICAL REFERENCE
!  SYSTEM TO THE INTERNATIONAL CELESTIAL REFERENCE SYSTEM (ICRS),
!  OR VICE VERSA.  THE DYNAMICAL REFERENCE SYSTEM IS BASED ON THE
!  DYNAMICAL MEAN EQUATOR AND EQUINOX OF J2000.0.  THE ICRS IS
!  BASED ON THE SPACE-FIXED ICRS AXES DEFINED BY THE RADIO CATALOG
!  POSITIONS OF SEVERAL HUNDRED EXTRAGALACTIC OBJECTS.  THE ROTATION
!  MATRIX USED HERE IS EQUIVALENT TO THAT GIVEN BY HILTON AND
!  HOHENKERK (2004), ASTRONOMY AND ASTROPHYSICS 413, 765-770,
!  EQ. (6) AND (8).
!
!       POS1   = POSITION VECTOR, EQUATORIAL RECTANGULAR
!                COORDINATES (IN)
!       K      = DIRECTION OF ROTATION (IN)
!                SET K < 0 FOR DYNAMICAL TO ICRS
!                SET K > 0 FOR ICRS TO DYNAMICAL
!       POS2   = POSITION VECTOR, EQUATORIAL RECTANGULAR
!                COORDINATES (OUT)
!
!  NOTE:  FOR GEOCENTRIC COORDINATES, THE SAME TRANSFORMATION IS
!  USED BETWEEN THE DYNAMICAL REFERENCE SYSTEM AND THE GCRS.

subroutine frame (pos1,k,pos2)

double precision pos1,pos2,pi,seccon,xi0,eta0,da0, &
     xx,yx,zx,xy,yy,zy,xz,yz,zz
dimension pos1(3), pos2(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

!     XI0, ETA0, AND DA0 ARE ICRS FRAME BIASES IN ARCSECONDS TAKEN
!     FROM IERS CONVENTIONS (2003), CHAPTER 5
data xi0, eta0, da0 / -0.0166170d0, -0.0068192d0, -0.01460d0 /

data ntimes / 0 /

ntimes = ntimes + 1

!     COMPUTE ELEMENTS OF ROTATION MATRIX (TO FIRST ORDER)
if ( ntimes <= 1 ) then
    xx =  1.d0
    yx = -da0  / seccon
    zx =  xi0  / seccon
    xy =  da0  / seccon
    yy =  1.d0
    zy =  eta0 / seccon
    xz = -xi0  / seccon
    yz = -eta0 / seccon
    zz =  1.d0
    !     INCLUDE SECOND-ORDER CORRECTIONS TO DIAGONAL ELEMENTS
    xx = 1.d0 - 0.5d0 * ( yx**2 + zx**2 )
    yy = 1.d0 - 0.5d0 * ( yx**2 + zy**2 )
    zz = 1.d0 - 0.5d0 * ( zy**2 + zx**2 )
end if

if ( k >= 0 ) then
    !     PERFORM ROTATION FROM ICRS TO DYNAMICAL SYSTEM
    pos2(1) = xx * pos1(1) + xy * pos1(2) + xz * pos1(3)
    pos2(2) = yx * pos1(1) + yy * pos1(2) + yz * pos1(3)
    pos2(3) = zx * pos1(1) + zy * pos1(2) + zz * pos1(3)
else
    !     PERFORM ROTATION FROM DYNAMICAL SYSTEM TO ICRS
    pos2(1) = xx * pos1(1) + yx * pos1(2) + zx * pos1(3)
    pos2(2) = xy * pos1(1) + yy * pos1(2) + zy * pos1(3)
    pos2(3) = xz * pos1(1) + yz * pos1(2) + zz * pos1(3)
end if

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE POSITION AND VELOCITY VECTORS OF
!  A TERRESTRIAL OBSERVER WITH RESPECT TO THE GEOCENTER.
!
!       GLON   = LONGITUDE OF OBSERVER WITH RESPECT TO REFERENCE
!                MERIDIAN (EAST +) IN DEGREES (IN)
!       GLAT   = GEODETIC LATITUDE (NORTH +) OF OBSERVER
!                IN DEGREES (IN)
!       HT     = HEIGHT OF OBSERVER IN METERS (IN)
!       ST     = LOCAL APPARENT SIDEREAL TIME AT REFERENCE MERIDIAN
!                IN HOURS (IN)
!       POS    = POSITION VECTOR OF OBSERVER WITH RESPECT TO
!                GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
!                REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
!                COMPONENTS IN AU (OUT)
!       VEL    = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO
!                GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
!                REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
!                COMPONENTS IN AU/DAY (OUT)
!
!  NOTE 1:  IF REFERENCE MERIDIAN IS GREENWICH AND ST=0.D0, POS
!  IS EFFECTIVELY REFERRED TO EQUATOR AND GREENWICH.
!
!  NOTE 2:  THIS SUBROUTINE IGNORES POLAR MOTION, UNLESS THE
!  OBSERVER'S LONGITUDE AND LATITUDE HAVE BEEN CORRECTED FOR IT,
!  AND VARIATION IN THE LENGTH OF DAY (ANGULAR VELOCITY OF EARTH).
!  NEGLECT OF POLAR MOTION MAY YIELD 15 METERS ERROR IN POSITION
!  AND OF ORDER 1 MILLIMETER/SEC ERROR IN VELOCITY.  NEGLECT OF
!  VARIATIONS IN LENGTH OF DAY RESULTS IN EVEN SMALLER VELOCITY
!  ERRORS.
!
!  NOTE 3:  THE TRUE EQUATOR AND EQUINOX OF DATE DO NOT FORM AN
!  INERTIAL SYSTEM.  THEREFORE, WITH RESPECT TO AN INERTIAL SYSTEM,
!  THE SMALL VELOCITY COMPONENT, OF ORDER 0.1 MILLIMETER/SEC,
!  DUE TO THE PRECESSION AND NUTATION OF THE EARTH'S AXIS, IS NOT
!  ACCOUNTED FOR HERE.

subroutine terra (glon,glat,ht,st,pos,vel)

double precision glon,glat,ht,st,pos,vel,pi,seccon,erad,f,omega, &
     aukm,df2,phi,sinphi,cosphi,c,s,ach,ash,stlocl,sinst,cosst, &
     dsqrt,dcos,dsin
dimension pos(3), vel(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

data ntimes / 0 /

ntimes = ntimes + 1
if (ntimes==1) then
!         GET ERAD, THE EQUATORIAL RADIUS OF EARTH IN KILOMETERS
    call astcon ('ERAD',1.d-3,erad)
!         GET F, THE FLATTENING FACTOR OF EARTH ELLIPSOID
    call astcon ('F',1.d0,f)
!         GET OMEGA, THE NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF
!         EARTH IN RADIANS/SECOND
    call astcon ('ANGVEL',1.d0,omega)
!         GET AUKM, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
    call astcon ('AU',1.d-3,aukm)
end if

!     COMPUTE PARAMETERS RELATING TO GEODETIC TO GEOCENTRIC CONVERSION
df2 = (1.d0 - f)**2
phi = glat * 3600.d0 / seccon
sinphi = dsin(phi)
cosphi = dcos(phi)
c = 1.d0 / dsqrt ( cosphi**2 + df2 * sinphi**2 )
s = df2 * c
ach = erad * c + ht/1000.d0
ash = erad * s + ht/1000.d0

!     COMPUTE LOCAL SIDEREAL TIME FACTORS
stlocl = (st * 54000.d0 + glon * 3600.d0) / seccon
sinst = dsin(stlocl)
cosst = dcos(stlocl)

!     COMPUTE POSITION VECTOR COMPONENTS IN KM
pos(1) = ach * cosphi * cosst
pos(2) = ach * cosphi * sinst
pos(3) = ash * sinphi

!     COMPUTE VELOCITY VECTOR COMPONENTS IN KM/SEC
vel(1) = -omega * ach * cosphi * sinst
vel(2) =  omega * ach * cosphi * cosst
vel(3) =  0.d0

!     CONVERT POSITION AND VELOCITY COMPONENTS TO AU AND AU/DAY
do j=1,3
    pos(j) = pos(j) / aukm
    vel(j) = vel(j) / aukm * 86400.d0
end do

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE TERRESTRIAL TIME (TT) JULIAN DATE
!  CORRESPONDING TO A BARYCENTRIC DYNAMICAL TIME (TDB) JULIAN DATE.
!  THE EXPRESSION USED IN THIS VERSION IS A TRUNCATED FORM OF A
!  LONGER AND MORE PRECISE SERIES GIVEN BY FAIRHEAD & BRETAGNON
!  (1990) A&A 229, 240.  THE RESULT IS GOOD TO ABOUT 10 MICROSECONDS.
!
!       TDBJD  = TDB JULIAN DATE (IN)
!       TTJD   = TT JULIAN DATE (OUT)
!       SECDIF = DIFFERENCE TDBJD-TTJD, IN SECONDS (OUT)

subroutine times (tdbjd,ttjd,secdif)

double precision tdbjd,ttjd,secdif,t,t0,dsin

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

t = ( tdbjd - t0 ) / 36525.d0

!     EXPRESSION GIVEN IN USNO CIRCULAR 179, EQ. 2.6
secdif = 0.001657d0 * dsin (  628.3076d0 * t + 6.2401d0) &              !
       + 0.000022d0 * dsin (  575.3385d0 * t + 4.2970d0) &
       + 0.000014d0 * dsin ( 1256.6152d0 * t + 6.1969d0) &              !
       + 0.000005d0 * dsin (  606.9777d0 * t + 4.0212d0) &              !
       + 0.000005d0 * dsin (   52.9691d0 * t + 0.4444d0) &
       + 0.000002d0 * dsin (   21.3299d0 * t + 5.5431d0) &
       + 0.000010d0 * t * dsin ( 628.3076d0 * t + 4.2490d0)

ttjd = tdbjd - secdif / 86400.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES QUANTITIES RELATED TO THE ORIENTATION
!  OF THE EARTH'S ROTATION AXIS AT JULIAN DATE TJD.
!
!       TJD    = TDB JULIAN DATE FOR ORIENTATION PARAMETERS (IN)
!       OBLM   = MEAN OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
!                DATE TJD (OUT)
!       OBLT   = TRUE OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
!                DATE TJD (OUT)
!       EQEQ   = EQUATION OF THE EQUINOXES IN TIME SECONDS AT
!                DATE TJD (OUT)
!       DPSI   = NUTATION IN LONGITUDE IN ARCSECONDS AT
!                DATE TJD (OUT)
!       DEPS   = NUTATION IN OBLIQUITY IN ARCSECONDS AT
!                DATE TJD (OUT)
!
!  NOTE:  THE EQUATION OF THE EQUINOXES INCLUDES THE COMPLEMENTARY
!  TERMS.

subroutine etilt (tjd,oblm,oblt,eqeq,dpsi,deps)

double precision tjd,oblm,oblt,eqeq,dpsi,deps,pi,seccon, &
     t0,tlast,t,psi,eps,psicor,epscor,cterms,delpsi,deleps, &
     el,elp,f,d,omega,obm,obt,ee, &
     dpole1,dpole2,dx,dy,dz,sine,x,dp1,dp2,dp3, &
     obliq,dabs,dsin,dcos !,eect2000
integer accdif
dimension dp1(3), dp2(3), dp3(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /,   mlast / 0 /
data delpsi, deleps, cterms, psicor, epscor / 5 * 0.d0 /

!     FUNCTION TO COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
!     CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412, 567-586,
!     EXPRESSION FROM EQ. (39) WITH OBLIQUITY AT J2000.0 TAKEN FROM
!     EQ. (37) OR TABLE 8
obliq(t) = ( ( ( ( -  0.0000000434d0   * t &
                   -  0.000000576d0  ) * t &
                   +  0.00200340d0   ) * t &
                   -  0.0001831d0    ) * t &
                   - 46.836769d0     ) * t + 84381.406d0

!     GET METHOD/ACCURACY MODE
call getmod ( mode )

!     CHECK FOR DIFFERENCE IN ACCURACY MODE FROM LAST CALL
accdif = mod ( mode, 2 ) - mod ( mlast, 2 )

t = ( tjd - t0 ) / 36525.d0

if ( dabs ( tjd - tlast ) > 1.d-8 .or. accdif /= 0 ) then

!         OBTAIN NUTATION PARAMETERS IN ARCSECONDS
    call nod ( t,   psi, eps )

!         OBTAIN COMPLEMENTARY TERMS FOR EQUATION OF THE EQUINOXES
!         IN ARCSECONDS
    if ( mod ( mode, 2 ) == 0 ) then
!             HIGH-ACCURACY MODE
        cterms = eect2000 ( tjd, 0.d0 ) * seccon
    else
!             LOW-ACCURACY MODE
        call funarg ( t,   el, elp, f, d, omega )
!             SERIES FROM IERS CONVENTIONS (2003), CHAPTER 5,
!             TABLE 5.2C, WITH SOME ADJUSTMENTS TO COEFFICIENT VALUES
!             COPIED FROM IERS FUNCTION EECT2000, WHICH HAS A MORE
!             COMPLETE SERIES
        cterms = &
          2640.96d-6 * dsin ( omega ) &
        +   63.52d-6 * dsin ( 2.d0 * omega ) &
        +   11.75d-6 * dsin ( 2.d0 * f - 2.d0 * d + 3.d0 * omega ) &    !
        +   11.21d-6 * dsin ( 2.d0 * f - 2.d0 * d +        omega ) &    !
        -    4.55d-6 * dsin ( 2.d0 * f - 2.d0 * d + 2.d0 * omega ) &    !
        +    2.02d-6 * dsin ( 2.d0 * f            + 3.d0 * omega ) &    !
        +    1.98d-6 * dsin ( 2.d0 * f            +        omega ) &    !
        -    1.72d-6 * dsin ( 3.d0 * omega ) &
        -    0.87d-6 * t * dsin ( omega )
!             (TERMS SMALLER THAN 2 MICROARCSECONDS OMITTED)
    end if
    tlast = tjd
    mlast = mode

end if

delpsi = psi + psicor
deleps = eps + epscor

!     COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
obm = obliq(t)

!     COMPUTE TRUE OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
obt = obm + deleps

!     COMPUTE EQUATION OF THE EQUINOXES IN ARCSECONDS
ee = delpsi * dcos ( obm / seccon ) + cterms

!     CONVERT TO OUTPUT UNITS
oblm = obm / 3600.d0
oblt = obt / 3600.d0
eqeq = ee  / 15.d0
dpsi = delpsi
deps = deleps

return


entry celpol (tjd,itype,dpole1,dpole2)
!
!     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF CELESTIAL POLE
!     OFFSETS FOR HIGH-PRECISION APPLICATIONS.  EACH SET OF OFFSETS IS
!     A CORRECTION TO THE MODELED POSITION OF THE POLE FOR A SPECIFIC
!     DATE, DERIVED FROM OBSERVATIONS AND PUBLISHED BY THE IERS.
!     THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE ANY OTHER ROUTINES
!     FOR A GIVEN DATE.  VALUES OF THE POLE OFFSETS SPECIFIED VIA A CALL
!     TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
!
!          TJD    = TDB OR TT JULIAN DATE FOR POLE OFFSETS (IN)
!          ITYPE  = TYPE OF POLE OFFSET (IN)
!                   SET ITYPE=1 FOR CORRECTIONS TO ANGULAR COORDINATES
!                               OF MODELED POLE REFERRED TO MEAN
!                               ECLIPTIC OF DATE, THAT IS,
!                               DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON
!                   SET ITYPE=2 FOR CORRECTIONS TO COMPONENTS OF
!                               MODELED POLE UNIT VECTOR WITH REFERRED
!                               TO GCRS AXES, THAT IS, DX AND DY
!          DPOLE1 = VALUE OF CELESTIAL POLE OFFSET IN FIRST COORDINATE,
!                   (DELTA-DELTA-PSI OR DX) IN MILLIARCSECONDS (IN)
!          DPOLE2 = VALUE OF CELESTIAL POLE OFFSET IN SECOND COORDINATE,
!                   (DELTA-DELTA-EPSILON OR DY) IN MILLIARCSECONDS (IN)
!
!     NOTE 1:  TJD IS USED ONLY FOR ITYPE=2, TO TRANSFORM DX AND DY TO
!     THE EQUIVALENT DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON VALUES.
!
!     NOTE 2:  FOR ITYPE=2, DX AND DY ARE UNIT VECTOR COMPONENT
!     CORRECTIONS, BUT ARE EXPRESSED IN MILLIARCSECONDS SIMPLY BY
!     MULTIPLYING BY 206264806, THE NUMBER OF MILLIARCSECONDS IN ONE
!     RADIAN.
!
!
if ( itype == 1 ) then

    psicor = dpole1 * 1.d-3
    epscor = dpole2 * 1.d-3

else

    dx = dpole1
    dy = dpole2

    t = ( tjd - t0 ) / 36525.d0
!         COMPUTE SINE OF MEAN OBLIQUITY OF DATE
    sine = dsin ( obliq(t) / seccon )

!         THE FOLLOWING ALGORITHM, TO TRANSFORM DX AND DY TO
!         DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON, IS FROM G. KAPLAN
!         (2003), USNO/AA TECHNICAL NOTE 2003-03, EQS. (7)-(9).

!         TRIVIAL MODEL OF POLE TRAJECTORY IN GCRS ALLOWS COMPUTATION
!         OF DZ
    x = ( 2004.19d0 * t ) / seccon
    dz = - ( x + 0.5d0 * x**3 ) * dx

!         FORM POLE OFFSET VECTOR (OBSERVED - MODELED) IN GCRS
    dp1(1) = dx * 1.d-3 / seccon
    dp1(2) = dy * 1.d-3 / seccon
    dp1(3) = dz * 1.d-3 / seccon

!         PRECESS POLE OFFSET VECTOR TO MEAN EQUATOR AND EQUINOX OF DATE
    call frame ( dp1, 1, dp2 )
    call preces ( t0, dp2, tjd, dp3 )

!         COMPUTE DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON IN ARCSECONDS
    psicor = ( dp3(1) / sine ) * seccon
    epscor = ( dp3(2)        ) * seccon

end if

return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES FUNDAMENTAL ARGUMENTS (MEAN ELEMENTS)
!  OF THE SUN AND MOON.  SEE SIMON ET AL. (1994) ASTRONOMY AND
!  ASTROPHYSICS 282, 663-683, ESPECIALLY SECTIONS 3.4-3.5.
!
!       T      = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
!       EL     = MEAN ANOMALY OF THE MOON IN RADIANS
!                AT DATE TJD (OUT)
!       ELP    = MEAN ANOMALY OF THE SUN IN RADIANS
!                AT DATE TJD (OUT)
!       F      = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
!                OF THE MOON'S ASCENDING NODE IN RADIANS
!                AT DATE TJD (OUT)
!       D      = MEAN ELONGATION OF THE MOON FROM THE SUN IN
!                RADIANS AT DATE TJD (OUT)
!       OMEGA  = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
!                   IN RADIANS AT DATE TJD (OUT)

subroutine funarg (t,el,elp,f,d,omega)

double precision t,el,elp,f,d,omega,pi,seccon,rev,dmod

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )
parameter ( rev    = 360.d0 * 3600.d0      )

!     FUNDAMENTAL (DELAUNAY) ARGUMENTS FROM SIMON ET AL. (1994)

!     MEAN ANOMALY OF THE MOON
el    = dmod (         485868.249036d0 + &
               t*( 1717915923.2178d0 + &
               t*(         31.8792d0 + &
               t*(          0.051635d0 + &
               t*(        - 0.00024470d0 )))), rev ) / seccon

!     MEAN ANOMALY OF THE SUN
elp   = dmod (        1287104.79305d0 + &
               t*(  129596581.0481d0 + &
               t*(        - 0.5532d0 + &
               t*(          0.000136d0 + &
               t*(        - 0.00001149d0 )))), rev ) / seccon

!     MEAN ARGUMENT OF THE LATITUDE OF THE MOON
f     = dmod (         335779.526232d0 + &
               t*( 1739527262.8478d0 + &
               t*(       - 12.7512d0 + &
               t*(       -  0.001037d0 + &
               t*(          0.00000417d0 )))), rev ) / seccon

!     MEAN ELONGATION OF THE MOON FROM THE SUN
d     = dmod (        1072260.70369d0 + &
               t*( 1602961601.2090d0 + &
               t*(        - 6.3706d0 + &
               t*(          0.006593d0 + &
               t*(        - 0.00003169d0 )))), rev ) / seccon

!     MEAN LONGITUDE OF THE ASCENDING NODE OF THE MOON (FROM SIMON
!     SECTION 3.4(b.3), PRECESSION=5028.8200 ARCSEC/CY)
omega = dmod (         450160.398036d0 + &
               t*(  - 6962890.5431d0 + &
               t*(          7.4722d0 + &
               t*(          0.007702d0 + &
               t*(        - 0.00005939d0 )))), rev ) / seccon

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES ATMOSPHERIC REFRACTION IN ZENITH
!  DISTANCE.  THIS VERSION COMPUTES APPROXIMATE REFRACTION FOR
!  OPTICAL WAVELENGTHS.  IT CAN BE USED FOR PLANNING OBSERVATIONS
!  OR TELESCOPE POINTING, BUT SHOULD NOT BE USED FOR THE REDUCTION
!  OF PRECISE OBSERVATIONS.  BASIC ALGORITHM IS DESCRIBED IN THE
!  EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, P. 144,
!  AND IS AN ADAPTATION OF A FORMULA IN BENNETT (1982), JOURNAL
!  OF NAVIGATION (ROYAL INSTITUTE) 35, 255-259.
!
!       HEIGHT = HEIGHT OF OBSERVER IN METERS (IN)
!       ZDOBS  = OBSERVED ZENITH DISTANCE IN DEGREES (IN)
!       REFR   = ATMOSPHERIC REFRACTION IN DEGREES (OUT)
!
!  NOTE:  HEIGHT IS NOT USED IF ENTRY REFDAT HAS BEEN CALLED
!  TO SPECIFY ATMOSPHERIC PRESSURE.

subroutine refrac (height,zdobs,refr)

double precision height,zdobs,refr,pi,degrad,s, &
     pobs,tobs,dobs,wlobs,obsp,obst,obsd,obswl,p,t,d,wl,h,r, &
     dexp,dtan
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( degrad = pi / 180.d0           )

!     S IS APPROXIMATE SCALE HEIGHT OF ATMOSPHERE IN METERS
data s / 9.1d3 /
data pobs, tobs, dobs, wlobs / 4 * -999.d0 /

!     COMPUTE REFRACTION ONLY FOR ZENITH DISTANCES
!     BETWEEN 0.1 AND 91 DEGREES
if ( zdobs < 0.1d0 .or. zdobs > 91.d0 ) then
    refr = 0.d0
    return
end if

!     IF OBSERVED WEATHER DATA ARE AVAILABLE, USE THEM
!     OTHERWISE, USE CRUDE ESTIMATES OF AVERAGE CONDITIONS
if ( pobs >= 1.d0 .and. tobs > -100.d0 ) then
    p  = pobs
    t  = tobs
    d  = dobs
    wl = wlobs
else
    p  = 1010.d0 * dexp ( -height / s )
    t  = 10.d0
    d  =  0.d0
    wl =  0.5d0
end if
!     D AND WL NOT USED IN THIS VERSION

h = 90.d0 - zdobs
r = 0.016667d0 / dtan ( ( h +  7.31d0 / ( h + 4.4d0 ) ) * degrad )      !
refr = r * ( 0.28d0 * p / ( t + 273.d0 ) )

return

entry refdat (obsp,obst,obsd,obswl)
!
!     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF WEATHER OBSERVATIONS
!     AND OTHER DATA TO BE USED IN THE ATMOSPHERIC REFRACTION
!     CALCULATION.  THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE
!     SUBROUTINE REFRAC OR ZDAZ FOR A GIVEN DATE/TIME.  DATA SPECIFIED
!     VIA A CALL TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
!
!          OBSP   = OBSERVED ATMOSPHERIC PRESSURE IN MILLIBARS (IN)
!          OBST   = OBSERVED TEMPERATURE IN DEGREES CELSIUS (IN)
!          OBSD   = OBSERVED DEW POINT IN DEGREES CELSIUS (IN)
!          OBSWL  = OBSERVING WAVELENGTH IN MICRONS (IN)
!
!     NOTE:  OBSD AND OBSWL ARE NOT USED IN THIS VERSION'S REFRACTION
!     ALGORITHM, AND CAN BE SET TO ANY VALUE.
!
!
pobs  = obsp
tobs  = obst
dobs  = obsd
wlobs = obswl
return

end
!***********************************************************************

!***********************************************************************
!>
!
!  THIS FUNCTION DETERMINES THE ANGLE OF AN OBJECT ABOVE OR BELOW
!  THE EARTH'S LIMB (HORIZON).  THE GEOMETRIC LIMB IS COMPUTED,
!  ASSUMING THE EARTH TO BE AN AIRLESS SPHERE (NO REFRACTION OR
!  OBLATENESS IS INCLUDED).  THE OBSERVER CAN BE ON OR ABOVE THE
!  EARTH.  FOR AN OBSERVER ON THE SURFACE OF THE EARTH, THIS
!  SUBROUTINE RETURNS THE APPROXIMATE UNREFRACTED ALTITUDE.
!
!       POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT GEOCENTER, COMPONENTS IN AU (IN)
!       POSO   = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
!                AT GEOCENTER, COMPONENTS IN AU (IN)
!       ALIMB  = ANGLE OF OBSERVED OBJECT ABOVE (+) OR BELOW (-) LIMB
!                IN DEGREES (OUT)
!       AFRAC  = NADIR ANGLE OF OBSERVED OBJECT AS A FRACTION OF
!                APPARENT RADIUS OF LIMB (OUT)
!                AFRAC<1.D0 MEANS BELOW THE LIMB
!                AFRAC=1.D0 MEANS ON THE LIMB
!                AFRAC>1.D0 MEANS ABOVE THE LIMB

subroutine limang (pos1,poso,alimb,afrac)

double precision pos1,poso,alimb,afrac,pi,halfpi,degcon, &
     erad,au,rade,disobj,disobs,aprad,zdlim,coszd,zdobj, &
     dsqrt,dacos
dimension pos1(3), poso(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( halfpi = 0.5d0 * pi            )
parameter ( degcon = 180.d0 / pi           )

data ntimes / 0 /

ntimes = ntimes + 1
if ( ntimes == 1 ) then
     call astcon ( 'ERAD', 1.d0,   erad )
     call astcon ( 'AU', 1.d0,   au )
     rade = erad / au
end if

disobj = dsqrt ( pos1(1)**2 + pos1(2)**2 + pos1(3)**2 )
disobs = dsqrt ( poso(1)**2 + poso(2)**2 + poso(3)**2 )

!     COMPUTE APPARENT ANGULAR RADIUS OF EARTH'S LIMB
if ( disobs > rade ) then
    aprad = dasin ( rade / disobs )
else
    aprad = halfpi
end if

!     COMPUTE ZENITH DISTANCE OF EARTH'S LIMB
zdlim = pi - aprad

!     COMPUTE ZENITH DISTANCE OF OBSERVED OBJECT
coszd = ( pos1(1)*poso(1) + pos1(2)*poso(2) + pos1(3)*poso(3) ) &
      / ( disobj * disobs )
if ( coszd <= -1.d0 ) then
   zdobj = pi
else if ( coszd >= 1.d0 ) then
   zdobj = 0.d0
else
   zdobj = dacos ( coszd )
end if

!     ANGLE OF OBJECT WRT LIMB IS DIFFERENCE IN ZENITH DISTANCES
alimb = ( zdlim - zdobj ) * degcon

!     NADIR ANGLE OF OBJECT AS A FRACTION OF ANGULAR RADIUS OF LIMB
afrac = ( pi - zdobj ) / aprad

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE RETURNS THE LOCATION OF THE CELESTIAL
!  INTERMEDIATE ORIGIN (CIO) FOR A GIVEN JULIAN DATE, AS A
!  RIGHT ASCENSION WITH RESPECT TO EITHER THE GCRS (GEOCENTRIC ICRS)
!  ORIGIN OR THE TRUE EQUINOX OF DATE.  THE CIO IS ALWAYS LOCATED ON
!  THE TRUE EQUATOR (=INTERMEDIATE EQUATOR) OF DATE.
!
!       TJD    = TDB JULIAN DATE (IN)
!       RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (OUT)
!       K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
!                GIVEN (OUT)
!                K=1 MEANS GCRS
!                K=2 MEANS TRUE EQUATOR AND EQUINOX OF DATE
!
!  NOTE:  IF AN EXTERNAL FILE OF CIO RIGHT ASCENSIONS IS AVAILABLE,
!  IT WILL BE USED AND K WILL BE SET TO 1.  OTHERWISE AN INTERNAL
!  COMPUTATION WILL BE USED AND K WILL BE SET TO 2.

subroutine cioloc ( tjd,   racio, k )

double precision tjd,racio,a,tlast,rlast,jd,ra,p,eqor,dabs
logical usefil
dimension jd(8), ra(8), a(1)
save

!     NUMBER OF VALUES IN ARRAYS FOR LAGRANGIAN INTERPOLATION
data m / 6 /

data tlast, rlast, klast / 0.d0, 0.d0, 0 /

3 format ( ' CIOLOC ERROR: CANNOT RETURN CIO RA VALUE FOR JD ', &
     f10.1 )

!     CHECK IF EXTERNAL FILE EXISTS
call ciord ( 0.d0, 1,   a, a, ierr )
usefil = ierr == 0

!     CHECK IF PREVIOUSLY COMPUTED RA VALUE CAN BE USED
if ( dabs ( tjd - tlast ) <= 1.d-8 ) then
    racio = rlast
    k = klast
    return
end if

! --- IF EXTERNAL FILE EXISTS, INTERPOLATE RA VALUES FROM FILE --------

if ( usefil ) then

    k = 1

!         GET ARRAYS OF VALUES TO INTERPOLATE
    call ciord ( tjd, m,   jd, ra, ierr )
    if ( ierr /= 0 ) then
        write ( *, 3 ) tjd
        racio = 99.d0
        return
    end if

!         PERFORM LAGRANGIAN INTERPOLATION FOR RA AT TJD
    racio = 0.d0
    do j = 1, m
        p = 1.d0
        do i = 1, m
            if ( i == j ) cycle
            p = p * ( tjd - jd(i) ) / ( jd(j) - jd(i) )
        end do
        racio = racio + p * ra(j)
    end do

    racio = racio / 54000.d0

! --- OTHERWISE, USE INTERNAL COMPUTATION ----------------------------

else

    k = 2

!         GET EQUATION OF THE ORIGINS
    call eqxra ( tjd, 1,   eqor )

    racio = -eqor

end if

! ---------------------------------------------------------------------

tlast = tjd
rlast = racio
klast = k

end
!***********************************************************************

!***********************************************************************
!>
!  GIVEN AN INPUT TDB JULIAN DATE AND THE NUMBER OF DATA POINTS
!  DESIRED, THIS SUBROUTINE RETURNS A SET OF JULIAN DATES AND
!  CORRESPONDING VALUES OF THE GCRS RIGHT ASCENSION OF THE CELESTIAL
!  INTERMEDIATE ORIGIN (CIO).  THE RANGE OF DATES IS CENTERED (AT LEAST
!  APPROXIMATELY) ON THE REQUESTED DATE.  THE SUBROUTINE OBTAINS THE
!  DATA FROM AN EXTERNAL DATA FILE.
!
!      TJD    = TDB JULIAN DATE (IN)
!      NVALS  = NUMBER OF JULIAN DATES AND RIGHT ASCENSION VALUES
!               REQUESTED (NOT LESS THAN 2 OR MORE THAN 20) (IN)
!      TLIST  = ARRAY OF TDB JULIAN DATES (OUT)
!      RALIST = ARRAY OF GCRS RIGHT ASCENSIONS OF THE CIO, FOR THE
!               JULIAN DATES IN TLIST, IN ARCSECONDS (OUT)
!      IERR   = ERROR INDICATOR (OUT)
!               IERR=0 MEANS EVERYTHING OK
!               IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
!               IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
!               IERR=3 MEANS BAD VALUE OF NVALS
!               IERR=4 MEANS EXTERNAL FILE CANNOT BE FOUND
!
!  NOTE:  TJD=0.D0 WITH NVALS=1 INDICATES A SPECIAL CALL JUST TO
!  DETERMINE IF EXTERNAL FILE EXISTS.

subroutine ciord ( tjd, nvals,   tlist, ralist, ierr )

double precision tjd,tlist,ralist,t,t1,r,tbeg,tend,tint,dif
character filnam*40, fileid*(*)
logical fileok
dimension tlist(nvals), ralist(nvals), t(20), r(20)
save

!     LOGICAL UNIT NUMBER AND FILE ID OF TIME SERIES OF CIO RA VALUES
data lu, ityp / 24, 1 /
data filnam / 'CIO_RA.TXT                              ' /

data ntimes, tbeg, tend, fileok / 0, 0.d0, 1.d10, .false. /

1 format ( a )
2 format ( ' CIORD ERROR: CANNOT FIND FILE ', a )
3 format ( ' CIORD ERROR: REQUESTED JD ', f10.1, 1x, a, &
      ' ALLOWED JD ', f10.1 )
4 format ( f16.6, f24.14 )

!     SPECIAL CALL JUST TO DETERMINE IF FILE EXITS
!     (NO PRINTED ERROR MESSAGE IF NOT)
if ( tjd == 0.d0 .and. nvals == 1 ) then
    ierr = 4
    if ( ntimes == 0 ) inquire ( file=filnam, exist=fileok )
    if ( fileok ) ierr = 0
    go to 77
end if

!     IF EXTERNAL FILE IS ALREADY KNOWN NOT TO EXIST, IMMEDIATELY
!     EXIT WITH IERR=4
if ( ntimes > 0 .and. .not. fileok ) then
    write ( *, 2 ) filnam
    ierr = 4
    go to 77
end if

!     CHECK FOR VALID VALUE OF NVALS
if ( nvals < 2 .or. nvals > 20 ) then
    ierr = 3
    go to 77
end if

middl = nvals / 2

!     CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF FILE
10 if ( tjd < tbeg ) then
    write ( *, 3 )  tjd, 'BEFORE FIRST', tbeg
    ierr = 1
    go to 77
else if ( tjd > tend ) then
    write ( *, 3 ) tjd, 'AFTER LAST', tend
    ierr = 2
    go to 77
end if

if ( ityp == 1 ) then

!         -------------------------------------------------------------
!         READ JULIAN DATES AND CIO RA VALUES FROM FORMATTED
!         SEQUENTIAL INPUT FILE

!         EACH RECORD OF THE FILE MUST CONTAIN A TDB JULIAN DATE
!         AND A CORRESPONDING CIO RA VALUE (WRT GCRS) IN ARCSECONDS

!         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL
!         -------------------------------------------------------------

!         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
    ntimes = ntimes + 1
    if ( ntimes == 1 ) then
        inquire ( file=filnam, exist=fileok )
        if ( .not. fileok ) then
            write ( *, 2 ) filnam
            ierr = 4
            go to 77
        end if
        open ( unit=lu, file=filnam, form='FORMATTED', &
               access='SEQUENTIAL', status='OLD' )
        read ( unit=lu, fmt=1 )
        do i = 1, nvals
            read ( unit=lu, fmt=4, end=50 ) t(i), r(i)
        end do
        tint = nint ( ( t(2) - t(1) ) * 1000.d0 ) / 1000.d0
        tbeg = t(middl)
        if ( tjd < tbeg ) go to 10
    end if

!         -------------------------------------------------------------

!         FILE READ SEQUENCE

20     dif = ( tjd - t(middl) ) / tint
    ndif = dif

!         BASIC DECISION ON HOW TO READ FILE
    if ( dif >= -0.1d0 .and. dif <= 1.1d0 ) then
!             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
        go to 70
    else if ( dif < 0.d0 ) then
!             REQUESTED JULIAN DATE BEFORE DATA PREVIOUSLY READ
        irec = ( t(nvals) - tbeg ) / tint
        nback = 3 * nvals
        if ( -dif <= 2 * nvals .and. irec > nback ) go to 34
        go to 25
    else if ( ndif > ( nvals + 1 ) ) then
!             REQUESTED JULIAN DATE FAR AHEAD OF DATA PREVIOUSLY READ
        nskip = ndif - nvals - 1
        go to 30
    else
!             REQUESTED JULIAN DATE A BIT AHEAD OF DATA PREVIOUSLY READ
        go to 40
    end if

!         REPOSITION FILE AT BEGINNING
25     rewind ( unit=lu )
    read ( unit=lu, fmt=1 )
    go to 36

!         FAST SKIP FORWARD
30     do i = 1, nskip
        read ( unit=lu, fmt=1, end=50 )
       end do
    go to 36

!         BACKSPACE FILE
34     do i = 1, nback
        backspace ( unit=lu )
       end do

!         FILL UP ARRAYS
36     do i = 1, nvals
        read ( unit=lu, fmt=4, end=50 ) t(i), r(i)
       end do
    go to 20

!         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
40     do i = 1, nvals - 1
        t(i) = t(i+1)
        r(i) = r(i+1)
       end do
    read ( unit=lu, fmt=4, end=50 ) t(nvals), r(nvals)

    go to 20

!         -------------------------------------------------------------

!         END OF FILE ENCOUNTERED
50  backspace ( unit=lu )
    backspace ( unit=lu )
    read ( unit=lu, fmt=4 ) tend
    tend = tend - ( nvals - middl - 1 ) * tint
    write ( *, 3 ) tjd, 'AFTER LAST', tend
    t(middl) = tend + tint
    ierr = 2
    go to 77

else if ( ityp == 2 ) then

!         -------------------------------------------------------------
!         READ JULIAN DATES AND CIO RA VALUES FROM UNFORMATTED
!         DIRECT ACCESS INPUT FILE

!         EACH RECORD OF THE FILE (EXCEPT THE FIRST) MUST CONTAIN A
!         TDB JULIAN DATE AND A CORRESPONDING CIO RA VALUE (WRT GCRS)
!         IN ARCSECONDS

!         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL

!         THE FIRST RECORD OF THE FILE IS SPECIAL AND MUST CONTAIN THE
!         TOTAL NUMBER OF RECORDS IN THE FILE
!         -------------------------------------------------------------

!         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
    ntimes = ntimes + 1
    if ( ntimes == 1 ) then
        inquire ( file=filnam, exist=fileok )
        if ( .not. fileok ) then
            write ( *, 2 ) filnam
            ierr = 4
            go to 77
        end if
        open ( unit=lu, file=filnam, form='UNFORMATTED', &
               access='DIRECT', recl=16, status='OLD' )
        read ( unit=lu, rec=1 ) nrecs
        do i = 1, nvals
            read ( unit=lu, rec=i+1 ) t(i), r(i)
        end do
        tint = nint ( ( t(2) - t(1) ) * 1000.d0 ) / 1000.d0
        tbeg = t(middl)
        tend = t(1) + ( nrecs - nvals + middl ) * tint
        t1 = t(1)
        lrec = 1
        maxrec = nrecs - nvals + 1
        if ( tjd < tbeg .or. tjd > tend ) go to 10
    end if

!         -------------------------------------------------------------

!         FILE READ SEQUENCE

60     dif = ( tjd - t(middl) ) / tint
!         IREC IS THE DATA RECORD NUMBER (PHYSICAL RECORD NUMBER - 1)
!         OF THE FIRST RECORD IN THE SEQUENCE OF NVALS RECORDS WITH
!         THE RELEVANT DATA TO BE RETURNED
    irec = ( ( tjd - t1 ) / tint ) - middl + 1.9d0
    if ( irec < 1      ) irec = 1
    if ( irec > maxrec ) irec = maxrec

!         BASIC DECISION ON HOW TO READ FILE
    if ( dif >= -0.1d0 .and. dif <= 1.1d0 ) then
!             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
        go to 70
    else if ( irec > lrec .and. irec - lrec <= middl ) then
!             REQUESTED JULIAN DATE JUST AHEAD OF DATA PREVIOUSLY READ
        go to 62
    else
!             REQUESTED JULIAN DATE IN DIFFERENT PART OF FILE
        go to 66
    end if

!         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
62     do i = 1, nvals - 1
        t(i) = t(i+1)
        r(i) = r(i+1)
       end do
    read ( unit=lu, rec=lrec+nvals+1 ) t(nvals), r(nvals)
    lrec = lrec + 1
    go to 60

!         GO DIRECTLY TO PROPER RECORD RANGE AND FILL UP ARRAYS
66     do i = 1, nvals
        read ( unit=lu, rec=irec+i ) t(i), r(i)
       end do
    lrec = irec

!         -------------------------------------------------------------

end if

!     GOT DATA, SO FILL OUTPUT ARRAYS
70 do i = 1, nvals
    tlist(i) = t(i)
    ralist(i) = r(i)
   end do
ierr = 0

77 return


entry ciofil ( lunit, fileid, itype )
!
!     THIS ENTRY ALLOWS SPECIFICATION OF THE LOGICAL UNIT NUMBER AND
!     FILE IDENTIFIER OF THE INPUT FILE CONTAINING A TIME SERIES OF CIO
!     RA VALUES.
!
!          LUNIT  = LOGIAL UNIT NUMBER TO BE USED FOR FILE (IN)
!          FILEID = FILE ID (IN)
!          ITYPE  = TYPE OF FILE (IN)
!                   SET ITYPE=1 FOR FORMATTED SEQUENTIAL FILE
!                   SET ITYPE=2 FOR UNFORMATTED BINARY FILE
!                   SET ITYPE=0 OR ANYTHING ELSE TO CLOSE THE CURRENT
!                               FILE (LUNIT AND FILEID IGNORED)
!
!     NOTE:  AFTER A CALL TO CIOFIL WITH ITYPE=0, THE ORIGINAL OR A
!     DIFFERENT FILE OF CIO RA VALUES CAN BE ACCESSED BY SUBSEQUENT
!     CALLS TO CIORD, BUT ONLY AFTER ANOTHER CALL TO CIOFIL WITH
!     ITYPE=1 OR 2.
!
!
if ( itype == 1 .or. itype == 2 ) then
    lu = lunit
    filnam = fileid
    ityp = itype
else
    close ( unit=lu )
    ntimes = 0
    tbeg = 0.d0
    tend = 1.d10
    fileok = .false.
end if

return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE RETURNS THE ORTHONORMAL BASIS VECTORS, WITH
!  RESPECT TO THE GCRS (GEOCENTRIC ICRS), OF THE CELESTIAL
!  INTERMEDIATE SYSTEM DEFINED BY THE CELESTIAL INTERMEDIATE POLE
!  (CIP) (IN THE Z DIRECTION) AND THE CELESTIAL INTERMEDIATE ORIGIN
!  (CIO) (IN THE X DIRECTION).  A TDB JULIAN DATE AND THE RIGHT
!  ASCENSION OF THE CIO AT THAT DATE IS REQUIRED AS INPUT.  THE
!  RIGHT ASCENSION OF THE CIO CAN BE WITH RESPECT TO EITHER THE
!  GCRS ORIGIN OR THE TRUE EQUINOX OF DATE -- DIFFERENT ALGORITHMS
!  ARE USED IN THE TWO CASES.
!
!       TJD    = TDB JULIAN DATE (IN)
!       RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (IN)
!       K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
!                EXPRESSED (IN)
!                SET K=1 FOR GCRS
!                SET K=2 FOR TRUE EQUATOR AND EQUINOX OF DATE
!       X      = UNIT VECTOR TOWARD THE CIO, EQUATORIAL RECTANGULAR
!                COORDINATES, REFERRED TO THE GCRS (OUT)
!       Y      = UNIT VECTOR TOWARD THE Y-DIRECTION, EQUATORIAL
!                RECTANGULAR COORDINATES, REFERRED TO THE GCRS (OUT)
!       Z      = UNIT VECTOR TOWARD NORTH CELESTIAL POLE (CIP),
!                EQUATORIAL RECTANGULAR COORDINATES, REFERRED TO
!                THE GCRS (OUT)

subroutine ciobas ( tjd, racio, k,   x, y, z )

double precision tjd,racio,x,y,z,xx,yy,zz,w0,w1,w2,z0, &
     pi,radcon,t0,tlast,sinra,cosra,xmag,dabs,dsin,dcos,dsqrt
dimension x(3), y(3), z(3), xx(3), yy(3), zz(3), z0(3), &
     w0(3), w1(3), w2(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( radcon = pi / 180.d0           )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data z0 / 0.d0, 0.d0, 1.d0 /,   tlast / 0.d0 /,   klast / 0 /

3 format ( ' CIOBAS ERROR: INVALID VALUE FOR K FOR JD ', &
     f10.1 )

!     USE LAST-COMPUTED BASIS VECTORS IF POSSIBLE
if ( dabs ( tjd - tlast ) <= 1.d-8 .and. k == klast ) &
   go to 60

!     COMPUTE UNIT VECTOR Z TOWARD CELESTIAL POLE (CIP)
call nutate ( -tjd, z0,   w1 )
call preces ( tjd, w1, t0,   w2 )
call frame ( w2, -1,   zz )

! --- RA OF CIO EXPRESSED IN GCRS -------------------------------------

if ( k == 1 ) then

!         COMPUTE VECTOR X TOWARD CIO IN GCRS
    sinra = dsin ( racio * 15.d0 * radcon )
    cosra = dcos ( racio * 15.d0 * radcon )
    xx(1) =  zz(3) * cosra
    xx(2) =  zz(3) * sinra
    xx(3) = -zz(1) * cosra - zz(2) * sinra

!         NORMALIZE VECTOR X
    xmag = dsqrt ( xx(1)**2 + xx(2)**2 + xx(3)**2 )
    xx(1) = xx(1) / xmag
    xx(2) = xx(2) / xmag
    xx(3) = xx(3) / xmag

!         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
    yy(1) = zz(2) * xx(3) - zz(3) * xx(2)
    yy(2) = zz(3) * xx(1) - zz(1) * xx(3)
    yy(3) = zz(1) * xx(2) - zz(2) * xx(1)

! --- RA OF CIO EXPRESSED IN EQUATOR-AND-EQUINOX OF DATE SYSTEM -------

else if ( k == 2 ) then

!         CONSTRUCT UNIT VECTOR TOWARD CIO
!         IN EQUATOR-AND-EQUINOX-OF-DATE SYSTEM
    w0(1) = dcos ( racio * 15.d0 * radcon )
    w0(2) = dsin ( racio * 15.d0 * radcon )
    w0(3) = 0.d0

!         ROTATE THE VECTOR INTO THE GCRS TO FORM UNIT VECTOR X
    call nutate ( -tjd, w0,   w1 )
    call preces ( tjd, w1, t0,   w2 )
    call frame ( w2, -1,   xx )

!         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
    yy(1) = zz(2) * xx(3) - zz(3) * xx(2)
    yy(2) = zz(3) * xx(1) - zz(1) * xx(3)
    yy(3) = zz(1) * xx(2) - zz(2) * xx(1)

! ---------------------------------------------------------------------

else

    write ( *, 3 ) tjd
    return

end if

! ---------------------------------------------------------------------

tlast = tjd
klast = k

60 do j = 1, 3
    x(j) = xx(j)
    y(j) = yy(j)
    z(j) = zz(j)
   end do

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE RETURNS THE VALUE OF THE EARTH ROTATION ANGLE
!  (THETA) FOR A GIVEN UT1 JULIAN DATE.  THE EXPRESSION USED IS
!  TAKEN FROM THE NOTE TO IAU RESOLUTION B1.8 OF 2000.
!
!      DATE1  = HIGH-ORDER PART OF UT1 JULIAN DATE (IN)
!      DATE2  = LOW-ORDER PART OF UT1 JULIAN DATE (IN)
!      THETA  = EARTH ROTATION ANGLE IN DEGREES (OUT)

subroutine erot (date1,date2,theta)

double precision date1, date2, theta, t0, thet1, thet2, thet3, &
     dmod

data t0 / 2451545.0d0 /

!     THE ALGORITHM USED BELOW IS EQUIVALENT TO THE CANNONICAL
!     THETA = 0.7790572732640D0 + 1.00273781191135448D0 * T,
!     WHERE T IS THE TIME IN UT1 DAYS FROM 2000.0 (T=DATE1+DATE2-T0),
!     BUT IT AVOIDS MANY TWO-PI 'WRAPS' THAT DECREASE PRECISION
!     (ADOPTED FROM SOFA ROUTINE IAU_ERA00 BY PAT WALLACE; SEE ALSO
!     EXPRESSION AT TOP OF PAGE 35 OF IERS CONVENTIONS (1996))

thet1 = 0.7790572732640d0 + 0.00273781191135448d0 * ( date1 - t0 )
thet2 =                     0.00273781191135448d0 *   date2
thet3 = dmod ( date1, 1.d0 ) + dmod ( date2, 1.d0 )
theta = dmod ( thet1 + thet2 + thet3, 1.d0 ) * 360.d0
if ( theta < 0.d0 ) theta = theta + 360.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES THE INTERMEDIATE RIGHT ASCENSION
!  OF THE EQUINOX AT JULIAN DATE TJD, USING AN ANALYTICAL EXPRESSION
!  FOR THE ACCUMULATED PRECESSION IN RIGHT ASCENSION.  FOR THE
!  TRUE EQUINOX THE RESULT IS THE EQUATION OF THE ORIGINS.
!
!       TJD    = TDB JULIAN DATE (IN)
!       K      = EQUINOX SELECTION CODE (IN)
!                SET K=0 FOR MEAN EQUINOX
!                SET K=1 FOR TRUE EQUINOX (EQUATION OF THE ORIGINS)
!       RADIF  = INTERMEDIATE RIGHT ASCENSION OF THE EQUINOX,
!                IN HOURS (+ OR -) (OUT)

subroutine eqxra ( tjd, k,    raeq )

double precision tjd,raeq,t0,tlast,ee,eqeq,t,a,precra,dabs
save

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /
data tlast / 0.d0 /,   ee / 0.d0 /

t = ( tjd - t0 ) / 36525.d0

!     FOR THE TRUE EQUINOX, OBTAIN THE EQUATION OF THE EQUINOXES IN
!     TIME SECONDS, WHICH INCLUDES THE 'COMPLIMENTARY TERMS'
if ( k == 1 ) then
    if ( dabs ( tjd - tlast ) > 1.d-8 ) then
        call etilt ( tjd,   a, a, ee, a, a )
        tlast = tjd
    end if
    eqeq = ee
else
    eqeq = 0.d0
end if

!     PRECESSION IN RA IN ARCSECONDS TAKEN FROM CAPITAINE ET AL. (2003),
!     ASTRONOMY AND ASTROPHYSICS 412, 567-586, EQ. (42)
precra = 0.014506d0 + &
         ( ( ( ( -    0.0000000368d0   * t &
                 -    0.000029956d0  ) * t &
                 -    0.00000044d0   ) * t &
                 +    1.3915817d0    ) * t &
                 + 4612.156534d0     ) * t

raeq = - ( precra / 15.d0 + eqeq ) / 3600.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE ALLOWS FOR THE SPECIFICATION OF THE DELTA-T VALUE
!  (DELTA-T = TT - UT1) TO BE USED IN THE CALCULATION OF SIDEREAL
!  TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS
!  THESE CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE
!  TIME ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME
!  ARGUMENT FOR THE PRECESSION AND NUTATION COMPONENTS.  THIS
!  SUBROUTINE, IF USED, SHOULD BE CALLED BEFORE ANY SUBROUTINE
!  RELATED TO EARTH ROTATION (E.G., SIDTIM OR TERCEL) FOR A GIVEN
!  DATE.  THE VALUE OF DELTA-T SPECIFIED HERE WILL BE USED UNTIL
!  EXPLICITLY CHANGED.
!
!       DELT   = VALUE OF DELTA-T (TT-UT1) IN SECONDS (IN)
!
!  NOTE 1:  THE COMPUTED VALUE OF SIDEREAL TIME, AND THE EQUIVALENT
!  EARTH ORIENTATION ANGLES, ARE RELATIVELY INSENSITIVE TO THE VALUE
!  OF DELTA-T: UP TO ONLY ~3 MICROARCSECONDS PER SECOND OF DELTA-T.
!  THEREFORE, FOR MANY APPLICATIONS, THIS SUBROUTINE EITHER NEED NOT
!  BE CALLED AT ALL, OR CAN BE CALLED JUST ONCE FOR A WIDE RANGE OF
!  DATES (E.G., A YEAR).  IF THIS CALL IS NOT USED, A DEFAULT
!  VALUE OF DELTA-T OF 64 SECONDS IS USED, WHICH IS APPROPRIATE TO
!  2000.0.
!
!  NOTE 2:  THE INPUT TIME ARGUMENTS TO SIDTIM AND TERCEL (TJDH AND
!  TJDL) ARE EXPRESSED IN UT1 REGARDLESS OF WHETHER THIS CALL IS
!  USED.

subroutine setdt ( delt )

double precision deltat, dt, delt
save dt

!     DEFAULT VALUE OF DELTA-T IN DAYS, EQUIVALENT TO 64 SECONDS,
!     THE APPROXIMATE VALUE AT 2000.0
data dt / 0.00074074d0 /

dt = delt / 86400.d0

return


entry getdt ( deltat )

!     THIS ENTRY RETURNS THE CURRENT VALUE OF DELTA-T
!     (DELTA-T = TT - UT1), PREVIOUSLY SET BY THE USER.  THE VALUE
!     RETURNED IS TO BE USED IN THE CALCULATION OF SIDEREAL TIME AND
!     THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS THESE
!     CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE TIME
!     ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME ARGUMENT
!     FOR THE PRECESSION AND NUTATION COMPONENTS.
!
!          DELTAT = VALUE OF DELTA-T (TT-UT1) IN DAYS (OUT)


deltat = dt

return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE ALLOWS THE USER TO SPECIFY THE 'MODE' VALUE,
!  WHICH DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
!  TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
!  ACCURACY OF NUTATION AND RELATED CALCULATIONS.
!
!       MODE   = SELECTION FOR METHOD AND ACCURACY (IN)
!                SET MODE=0 FOR CIO-BASED METHOD, FULL ACCURACY
!                SET MODE=1 FOR CIO-BASED METHOD, REDUCED ACCURACY
!                SET MODE=2 FOR EQUINOX-BASED METHOD, FULL ACCURACY
!                SET MODE=3 FOR EQUINOX-BASED METHOD, REDUCED
!                               ACCURACY
!
!  NOTE: OTHER ENTRY POINTS ARE PROVIDED TO ALLOW THE METHOD AND
!  ACCURACY TO BE SPECIFIED IN A MORE OBVIOUS WAY:
!  MODE=0 CAN BE SET BY CALL CIOTIO AND CALL HIACC
!  MODE=1 CAN BE SET BY CALL CIOTIO AND CALL LOACC
!  MODE=2 CAN BE SET BY CALL EQINOX AND CALL HIACC
!  MODE=3 CAN BE SET BY CALL EQINOX AND CALL LOACC

subroutine setmod ( mode )

save imode, lmode

data imode, lmode / 2, 2 /

lmode = imode
imode = mode

return


entry ciotio
lmode = imode
if ( imode >= 2 ) imode = imode - 2
return


entry eqinox
lmode = imode
if ( imode <= 1 ) imode = imode + 2
return


entry loacc
lmode = imode
if ( mod ( imode, 2 ) == 0 ) imode = imode + 1
return


entry hiacc
lmode = imode
if ( mod ( imode, 2 ) == 1 ) imode = imode - 1
return


entry resume
imode = lmode
return


entry getmod ( mode )
!
!     THIS SUBROUTINE RETURNS THE 'MODE' VALUE, WHICH
!     DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
!     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
!     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
!
!          MODE   = SELECTION FOR METHOD AND ACCURACY (OUT)
!                   MODE=0 MEANS CIO-BASED METHOD, FULL ACCURACY
!                   MODE=1 MEANS CIO-BASED METHOD, REDUCED ACCURACY
!                   MODE=2 MEANS EQUINOX-BASED METHOD, FULL ACCURACY
!                   MODE=3 MEANS EQUINOX-BASED METHOD, REDUCED ACCURACY
!
!
mode = imode

return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE ALLOWS THE USER TO RETRIEVE THE LAST COMPUTED
!  POSITION ON THE SKY AS A UNIT VECTOR.
!
!       UNITV  = UNIT VECTOR TOWARD LAST COMPUTED POSITION ON THE
!                SKY, IN THE COORDINATE SYSTEM USED FOR THAT
!                POSITION (OUT)

subroutine getvec ( unitv )

double precision unitv, p, pos, r, dsqrt
dimension unitv(3), p(3), pos(3)
save p

r = dsqrt ( p(1)**2 + p(2)**2 + p(3)**2 )

do j = 1, 3
    unitv(j) = p(j) / r
end do

return


entry setvec ( pos )
!
!     THIS ENTRY STORES THE LAST COMPUTED POSITION ON THE SKY.
!
!          POS    = VECTOR TOWARD LAST COMPUTED POSITION ON THE
!                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
!                   POSITION (IN)
!
!
do j = 1, 3
    p(j) = pos(j)
end do

return

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
!  TIME.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT TIME VALUE
!  CAN BE IN ANY UT-LIKE TIME SCALE (UTC, UT1, TT, ETC.) - OUTPUT
!  JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
!  VAN FLANDERN.
!
!       I      = YEAR (IN)
!       M      = MONTH NUMBER (IN)
!       K      = DAY OF MONTH (IN)
!       H      = UT HOURS (IN)
!       TJD    = JULIAN DATE (OUT)

subroutine juldat (i,m,k,h,tjd)

double precision h,tjd

!     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
jd = k-32075+1461*(i+4800+(m-14)/12)/4+367*(m-2-(m-14)/12*12)/12 &
     -3*((i+4900+(m-14)/12)/100)/4
tjd = jd - 0.5d0 + h/24.d0

end
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE COMPUTES CALENDAR DATE AND TIME, GIVEN JULIAN
!  DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE TIME SCALE
!  (UTC, UT1, TT, ETC.) - OUTPUT TIME VALUE WILL HAVE SAME BASIS.
!  OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
!  VAN FLANDERN.
!
!       TJD    = JULIAN DATE (IN)
!       I      = YEAR (OUT)
!       M      = MONTH NUMBER (OUT)
!       K      = DAY OF MONTH (OUT)
!       H      = UT HOURS (OUT)

subroutine caldat (tjd,i,m,k,h)

double precision tjd,h,djd,dmod

djd = tjd + 0.5d0
jd = djd
h = dmod (djd,1.d0) * 24.d0
!     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
l = jd + 68569
n = 4*l/146097
l = l - (146097*n+3)/4
!     I=YEAR, M=MONTH, K=DAY
i = 4000*(l+1)/1461001
l = l - 1461*i/4 + 31
m = 80*l/2447
k = l - 2447*m/80
l = m / 11
m = m + 2 - 12*l
i = 100*(n-49) + i + l

end subroutine caldat
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE SUPPLIES THE VALUES OF ASTRONOMICAL CONSTANTS.
!
!      NAME   = NAME OF CONSTANT WHOSE VALUE IS DESIRED (IN)
!               'C'         SPEED OF LIGHT IN METERS/SECOND
!               'C(AU/DAY)' SPEED OF LIGHT IN AU/DAY
!               'AU'        LENGTH OF ASTRONOMICAL UNIT IN METERS
!               'AU(SEC)'   LENGTH OF ASTRONOMICAL UNIT IN SECONDS
!               'GS'        HELIOCENTRIC GRAVITATIONAL CONSTANT
!                              IN METERS**3/SECOND**2
!               'GE'        GEOCENTRIC GRAVITATIONAL CONSTANT
!                              IN METERS**3/SECOND**2
!               'ERAD'      EQUATORIAL RADIUS OF EARTH IN METERS
!               'F'         FLATTENING FACTOR OF EARTH
!               'ANGVEL'    NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY
!                              OF EARTH IN RADIANS/SECOND
!               'MASS_SUN'  RECIPROCAL MASS OF THE SUN
!               'MASS_EAR'  RECIPROCAL MASS OF THE EARTH
!               'MASS_MOO'  RECIPROCAL MASS OF THE MOON
!               'MASS_MER'  RECIPROCAL MASS OF MERCURY
!                   :             :      :        :
!               'MASS_PLU'  RECIPROCAL MASS OF PLUTO
!      FACTOR = FACTOR BY WHICH CONSTANT VALUE IS TO BE MULTIPLIED
!               (IN)
!      CONST  = CONSTANT VALUE AFTER MULTIPLICATION BY FACTOR (OUT)

subroutine astcon (name,factor,const)

double precision factor,const,c,ausec
character name*(*)

!     NOTE:  THESE CONSTANT VALUES ARE BASED ON THE TDB SECOND WHERE
!     APPLICABLE.

!     SPEED OF LIGHT IN METERS/SECOND IS A DEFINING PHYSICAL CONSTANT
data c / 299792458.d0 /

!     LIGHT-TIME FOR ONE ASTRONOMICAL UNIT IN SECONDS, FROM DE-405
data ausec / 499.0047838061d0 /

!     SPEED OF LIGHT IN METERS/SECOND
if ( name == 'C' ) then
    const = c

!     SPEED OF LIGHT IN AU/DAY
else if ( name == 'C(AU/DAY)' ) then
    const = 86400.d0 / ausec

!     LENGTH OF ASTRONOMICAL UNIT IN METERS
else if ( name == 'AU' ) then
    const = ausec * c

!     LENGTH OF ASTRONOMICAL UNIT IN SECONDS
else if ( name == 'AU(SEC)' ) then
    const = ausec

!     HELIOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
!     DE-405
else if ( name == 'GS' ) then
    const = 1.32712440017987d20

!     GEOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
!     DE-405
else if ( name == 'GE' ) then
    const = 3.98600433d14

!     EQUATORIAL RADIUS OF EARTH IN METERS, FROM IERS CONVENTIONS (2003)
else if ( name == 'ERAD' ) then
    const = 6378136.6d0

!     FLATTENING FACTOR OF EARTH, FROM IERS CONVENTIONS (2003)
else if ( name == 'F' ) then
    const = 1.d0 / 298.25642d0

!     NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF EARTH
!     IN RADIANS/SECOND, FROM IERS CONVENTIONS (2003)
else if ( name == 'ANGVEL' ) then
    const = 7.2921150d-5

!     RECIPROCAL MASSES OF SOLAR SYSTEM BODIES, FROM DE-405
!     (SUN MASS / BODY MASS)
else if ( name(1:4) == 'MASS' ) then

    const = 1.d0
    if ( name(6:8) == 'SUN' ) const =         1.d0
    if ( name(6:8) == 'MOO' ) const =  27068700.387534d0
    if ( name(6:8) == 'MER' ) const =   6023600.d0
    if ( name(6:8) == 'VEN' ) const =    408523.71d0
    if ( name(6:8) == 'EAR' ) const =    332946.050895d0
    if ( name(6:8) == 'MAR' ) const =   3098708.d0
    if ( name(6:8) == 'JUP' ) const =      1047.3486d0
    if ( name(6:8) == 'SAT' ) const =      3497.898d0
    if ( name(6:8) == 'URA' ) const =     22902.98d0
    if ( name(6:8) == 'NEP' ) const =     19412.24d0
    if ( name(6:8) == 'PLU' ) const = 135200000.d0
    if ( name(6:8) == 'EMB' ) const =    328900.561400d0

end if

const = const * factor

end subroutine astcon
!***********************************************************************

!***********************************************************************
!>
!  THIS SUBROUTINE RETURNS THE VALUES FOR NUTATION IN LONGITUDE AND
!  NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
!
!       T     = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
!       DPSI  = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
!       DEPS  = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)

subroutine nod (t,dpsi,deps)

double precision t,dpsi,deps,pi,seccon,t0,t1,dp,de
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

!     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

!     GET METHOD/ACCURACY MODE
call getmod ( mode )

t1 = t * 36525.d0

!     =================================================================
!     EVALUATE NUTATION SERIES
!     RESULTING NUTATION IN LONGITUDE AND OBLIQUITY IN ARCSECONDS

!     CALL SUBROUTINE TO EVALUATE NUTATION SERIES
if ( mod ( mode, 2 ) == 0 ) then
!         HIGH ACCURACY MODE -- IERS SUBROUTINE
    call nu2000a ( t0, t1, dp, de )
else
!         LOW ACCURACY MODE -- MODIFICATION OF IERS SUBROUTINE
    call nu2000k ( t0, t1, dp, de )
end if
dpsi = dp * seccon
deps = de * seccon

end subroutine nod
!***********************************************************************

!***********************************************************************
!>
!  Nutation, IAU 2000A model (MHB_2000 without FCN).
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Given:
!     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
!
!  Returned:
!     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
!
!  This revision:  2002 November 25

subroutine nu2000a ( date1, date2, dpsi, deps )

implicit none

double precision date1, date2, dpsi, deps

!  Arcseconds to radians
double precision das2r
parameter ( das2r = 4.848136811095359935899141d-6 )

!  Milliarcseconds to radians
double precision dmas2r
parameter ( dmas2r = das2r / 1d3 )

!  Arc seconds in a full circle
double precision turnas
parameter ( turnas = 1296000d0 )

!  2Pi
double precision d2pi
parameter ( d2pi = 6.283185307179586476925287d0 )

!  Units of 0.1 microarcsecond to radians
double precision u2r
parameter ( u2r = das2r/1d7 )

!  Reference epoch (J2000), JD
double precision dj0
parameter ( dj0 = 2451545d0 )

!  Days per Julian century
double precision djc
parameter ( djc = 36525d0 )

!  Miscellaneous
double precision t, el, elp, f, d, om, arg, dp, de, sarg, carg, &
                 dpsils, depsls, &
                 al, alsu, af, ad, aom, alme, alve, alea, alma, &
                 alju, alsa, alur, alne, apa, dpsipl, depspl
integer i, j

!  -------------------------
!  Luni-Solar nutation model
!  -------------------------

!  Number of terms in the luni-solar nutation model
integer nls
parameter ( nls = 678 )

!  Coefficients for fundamental arguments
integer nals(5,nls)

!  Longitude and obliquity coefficients
double precision cls(6,nls)

!  ---------------
!  Planetary terms
!  ---------------

!  Number of terms in the planetary nutation model
integer npl
parameter ( npl = 687 )

!  Coefficients for fundamental arguments
integer napl(14,npl)

!  Longitude and obliquity coefficients
integer icpl(4,npl)

!  ----------------------------------------
!  Tables of argument and term coefficients
!  ----------------------------------------

!
!  Luni-Solar argument multipliers:
!               L     L'    F     D     Om

data ( ( nals(i,j), i=1,5 ), j=  1, 10 ) / &
          0,    0,    0,    0,    1, &
          0,    0,    2,   -2,    2, &
          0,    0,    2,    0,    2, &
          0,    0,    0,    0,    2, &
          0,    1,    0,    0,    0, &
          0,    1,    2,   -2,    2, &
          1,    0,    0,    0,    0, &
          0,    0,    2,    0,    1, &
          1,    0,    2,    0,    2, &
          0,   -1,    2,   -2,    2 /
data ( ( nals(i,j), i=1,5 ), j= 11, 20 ) / &
          0,    0,    2,   -2,    1, &
         -1,    0,    2,    0,    2, &
         -1,    0,    0,    2,    0, &
          1,    0,    0,    0,    1, &
         -1,    0,    0,    0,    1, &
         -1,    0,    2,    2,    2, &
          1,    0,    2,    0,    1, &
         -2,    0,    2,    0,    1, &
          0,    0,    0,    2,    0, &
          0,    0,    2,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j= 21, 30 ) / &
          0,   -2,    2,   -2,    2, &
         -2,    0,    0,    2,    0, &
          2,    0,    2,    0,    2, &
          1,    0,    2,   -2,    2, &
         -1,    0,    2,    0,    1, &
          2,    0,    0,    0,    0, &
          0,    0,    2,    0,    0, &
          0,    1,    0,    0,    1, &
         -1,    0,    0,    2,    1, &
          0,    2,    2,   -2,    2 /
data ( ( nals(i,j), i=1,5 ), j= 31, 40 ) / &
          0,    0,   -2,    2,    0, &
          1,    0,    0,   -2,    1, &
          0,   -1,    0,    0,    1, &
         -1,    0,    2,    2,    1, &
          0,    2,    0,    0,    0, &
          1,    0,    2,    2,    2, &
         -2,    0,    2,    0,    0, &
          0,    1,    2,    0,    2, &
          0,    0,    2,    2,    1, &
          0,   -1,    2,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j= 41, 50 ) / &
          0,    0,    0,    2,    1, &
          1,    0,    2,   -2,    1, &
          2,    0,    2,   -2,    2, &
         -2,    0,    0,    2,    1, &
          2,    0,    2,    0,    1, &
          0,   -1,    2,   -2,    1, &
          0,    0,    0,   -2,    1, &
         -1,   -1,    0,    2,    0, &
          2,    0,    0,   -2,    1, &
          1,    0,    0,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j= 51, 60 ) / &
          0,    1,    2,   -2,    1, &
          1,   -1,    0,    0,    0, &
         -2,    0,    2,    0,    2, &
          3,    0,    2,    0,    2, &
          0,   -1,    0,    2,    0, &
          1,   -1,    2,    0,    2, &
          0,    0,    0,    1,    0, &
         -1,   -1,    2,    2,    2, &
         -1,    0,    2,    0,    0, &
          0,   -1,    2,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j= 61, 70 ) / &
         -2,    0,    0,    0,    1, &
          1,    1,    2,    0,    2, &
          2,    0,    0,    0,    1, &
         -1,    1,    0,    1,    0, &
          1,    1,    0,    0,    0, &
          1,    0,    2,    0,    0, &
         -1,    0,    2,   -2,    1, &
          1,    0,    0,    0,    2, &
         -1,    0,    0,    1,    0, &
          0,    0,    2,    1,    2 /
data ( ( nals(i,j), i=1,5 ), j= 71, 80 ) / &
         -1,    0,    2,    4,    2, &
         -1,    1,    0,    1,    1, &
          0,   -2,    2,   -2,    1, &
          1,    0,    2,    2,    1, &
         -2,    0,    2,    2,    2, &
         -1,    0,    0,    0,    2, &
          1,    1,    2,   -2,    2, &
         -2,    0,    2,    4,    2, &
         -1,    0,    4,    0,    2, &
          2,    0,    2,   -2,    1 /
data ( ( nals(i,j), i=1,5 ), j= 81, 90 ) / &
          2,    0,    2,    2,    2, &
          1,    0,    0,    2,    1, &
          3,    0,    0,    0,    0, &
          3,    0,    2,   -2,    2, &
          0,    0,    4,   -2,    2, &
          0,    1,    2,    0,    1, &
          0,    0,   -2,    2,    1, &
          0,    0,    2,   -2,    3, &
         -1,    0,    0,    4,    0, &
          2,    0,   -2,    0,    1 /
data ( ( nals(i,j), i=1,5 ), j= 91,100 ) / &
         -2,    0,    0,    4,    0, &
         -1,   -1,    0,    2,    1, &
         -1,    0,    0,    1,    1, &
          0,    1,    0,    0,    2, &
          0,    0,   -2,    0,    1, &
          0,   -1,    2,    0,    1, &
          0,    0,    2,   -1,    2, &
          0,    0,    2,    4,    2, &
         -2,   -1,    0,    2,    0, &
          1,    1,    0,   -2,    1 /
data ( ( nals(i,j), i=1,5 ), j=101,110 ) / &
         -1,    1,    0,    2,    0, &
         -1,    1,    0,    1,    2, &
          1,   -1,    0,    0,    1, &
          1,   -1,    2,    2,    2, &
         -1,    1,    2,    2,    2, &
          3,    0,    2,    0,    1, &
          0,    1,   -2,    2,    0, &
         -1,    0,    0,   -2,    1, &
          0,    1,    2,    2,    2, &
         -1,   -1,    2,    2,    1 /
data ( ( nals(i,j), i=1,5 ), j=111,120 ) / &
          0,   -1,    0,    0,    2, &
          1,    0,    2,   -4,    1, &
         -1,    0,   -2,    2,    0, &
          0,   -1,    2,    2,    1, &
          2,   -1,    2,    0,    2, &
          0,    0,    0,    2,    2, &
          1,   -1,    2,    0,    1, &
         -1,    1,    2,    0,    2, &
          0,    1,    0,    2,    0, &
          0,   -1,   -2,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=121,130 ) / &
          0,    3,    2,   -2,    2, &
          0,    0,    0,    1,    1, &
         -1,    0,    2,    2,    0, &
          2,    1,    2,    0,    2, &
          1,    1,    0,    0,    1, &
          1,    1,    2,    0,    1, &
          2,    0,    0,    2,    0, &
          1,    0,   -2,    2,    0, &
         -1,    0,    0,    2,    2, &
          0,    1,    0,    1,    0 /
data ( ( nals(i,j), i=1,5 ), j=131,140 ) / &
          0,    1,    0,   -2,    1, &
         -1,    0,    2,   -2,    2, &
          0,    0,    0,   -1,    1, &
         -1,    1,    0,    0,    1, &
          1,    0,    2,   -1,    2, &
          1,   -1,    0,    2,    0, &
          0,    0,    0,    4,    0, &
          1,    0,    2,    1,    2, &
          0,    0,    2,    1,    1, &
          1,    0,    0,   -2,    2 /
data ( ( nals(i,j), i=1,5 ), j=141,150 ) / &
         -1,    0,    2,    4,    1, &
          1,    0,   -2,    0,    1, &
          1,    1,    2,   -2,    1, &
          0,    0,    2,    2,    0, &
         -1,    0,    2,   -1,    1, &
         -2,    0,    2,    2,    1, &
          4,    0,    2,    0,    2, &
          2,   -1,    0,    0,    0, &
          2,    1,    2,   -2,    2, &
          0,    1,    2,    1,    2 /
data ( ( nals(i,j), i=1,5 ), j=151,160 ) / &
          1,    0,    4,   -2,    2, &
         -1,   -1,    0,    0,    1, &
          0,    1,    0,    2,    1, &
         -2,    0,    2,    4,    1, &
          2,    0,    2,    0,    0, &
          1,    0,    0,    1,    0, &
         -1,    0,    0,    4,    1, &
         -1,    0,    4,    0,    1, &
          2,    0,    2,    2,    1, &
          0,    0,    2,   -3,    2 /
data ( ( nals(i,j), i=1,5 ), j=161,170 ) / &
         -1,   -2,    0,    2,    0, &
          2,    1,    0,    0,    0, &
          0,    0,    4,    0,    2, &
          0,    0,    0,    0,    3, &
          0,    3,    0,    0,    0, &
          0,    0,    2,   -4,    1, &
          0,   -1,    0,    2,    1, &
          0,    0,    0,    4,    1, &
         -1,   -1,    2,    4,    2, &
          1,    0,    2,    4,    2 /
data ( ( nals(i,j), i=1,5 ), j=171,180 ) / &
         -2,    2,    0,    2,    0, &
         -2,   -1,    2,    0,    1, &
         -2,    0,    0,    2,    2, &
         -1,   -1,    2,    0,    2, &
          0,    0,    4,   -2,    1, &
          3,    0,    2,   -2,    1, &
         -2,   -1,    0,    2,    1, &
          1,    0,    0,   -1,    1, &
          0,   -2,    0,    2,    0, &
         -2,    0,    0,    4,    1 /
data ( ( nals(i,j), i=1,5 ), j=181,190 ) / &
         -3,    0,    0,    0,    1, &
          1,    1,    2,    2,    2, &
          0,    0,    2,    4,    1, &
          3,    0,    2,    2,    2, &
         -1,    1,    2,   -2,    1, &
          2,    0,    0,   -4,    1, &
          0,    0,    0,   -2,    2, &
          2,    0,    2,   -4,    1, &
         -1,    1,    0,    2,    1, &
          0,    0,    2,   -1,    1 /
data ( ( nals(i,j), i=1,5 ), j=191,200 ) / &
          0,   -2,    2,    2,    2, &
          2,    0,    0,    2,    1, &
          4,    0,    2,   -2,    2, &
          2,    0,    0,   -2,    2, &
          0,    2,    0,    0,    1, &
          1,    0,    0,   -4,    1, &
          0,    2,    2,   -2,    1, &
         -3,    0,    0,    4,    0, &
         -1,    1,    2,    0,    1, &
         -1,   -1,    0,    4,    0 /
data ( ( nals(i,j), i=1,5 ), j=201,210 ) / &
         -1,   -2,    2,    2,    2, &
         -2,   -1,    2,    4,    2, &
          1,   -1,    2,    2,    1, &
         -2,    1,    0,    2,    0, &
         -2,    1,    2,    0,    1, &
          2,    1,    0,   -2,    1, &
         -3,    0,    2,    0,    1, &
         -2,    0,    2,   -2,    1, &
         -1,    1,    0,    2,    2, &
          0,   -1,    2,   -1,    2 /
data ( ( nals(i,j), i=1,5 ), j=211,220 ) / &
         -1,    0,    4,   -2,    2, &
          0,   -2,    2,    0,    2, &
         -1,    0,    2,    1,    2, &
          2,    0,    0,    0,    2, &
          0,    0,    2,    0,    3, &
         -2,    0,    4,    0,    2, &
         -1,    0,   -2,    0,    1, &
         -1,    1,    2,    2,    1, &
          3,    0,    0,    0,    1, &
         -1,    0,    2,    3,    2 /
data ( ( nals(i,j), i=1,5 ), j=221,230 ) / &
          2,   -1,    2,    0,    1, &
          0,    1,    2,    2,    1, &
          0,   -1,    2,    4,    2, &
          2,   -1,    2,    2,    2, &
          0,    2,   -2,    2,    0, &
         -1,   -1,    2,   -1,    1, &
          0,   -2,    0,    0,    1, &
          1,    0,    2,   -4,    2, &
          1,   -1,    0,   -2,    1, &
         -1,   -1,    2,    0,    1 /
data ( ( nals(i,j), i=1,5 ), j=231,240 ) / &
          1,   -1,    2,   -2,    2, &
         -2,   -1,    0,    4,    0, &
         -1,    0,    0,    3,    0, &
         -2,   -1,    2,    2,    2, &
          0,    2,    2,    0,    2, &
          1,    1,    0,    2,    0, &
          2,    0,    2,   -1,    2, &
          1,    0,    2,    1,    1, &
          4,    0,    0,    0,    0, &
          2,    1,    2,    0,    1 /
data ( ( nals(i,j), i=1,5 ), j=241,250 ) / &
          3,   -1,    2,    0,    2, &
         -2,    2,    0,    2,    1, &
          1,    0,    2,   -3,    1, &
          1,    1,    2,   -4,    1, &
         -1,   -1,    2,   -2,    1, &
          0,   -1,    0,   -1,    1, &
          0,   -1,    0,   -2,    1, &
         -2,    0,    0,    0,    2, &
         -2,    0,   -2,    2,    0, &
         -1,    0,   -2,    4,    0 /
data ( ( nals(i,j), i=1,5 ), j=251,260 ) / &
          1,   -2,    0,    0,    0, &
          0,    1,    0,    1,    1, &
         -1,    2,    0,    2,    0, &
          1,   -1,    2,   -2,    1, &
          1,    2,    2,   -2,    2, &
          2,   -1,    2,   -2,    2, &
          1,    0,    2,   -1,    1, &
          2,    1,    2,   -2,    1, &
         -2,    0,    0,   -2,    1, &
          1,   -2,    2,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=261,270 ) / &
          0,    1,    2,    1,    1, &
          1,    0,    4,   -2,    1, &
         -2,    0,    4,    2,    2, &
          1,    1,    2,    1,    2, &
          1,    0,    0,    4,    0, &
          1,    0,    2,    2,    0, &
          2,    0,    2,    1,    2, &
          3,    1,    2,    0,    2, &
          4,    0,    2,    0,    1, &
         -2,   -1,    2,    0,    0 /
data ( ( nals(i,j), i=1,5 ), j=271,280 ) / &
          0,    1,   -2,    2,    1, &
          1,    0,   -2,    1,    0, &
          0,   -1,   -2,    2,    1, &
          2,   -1,    0,   -2,    1, &
         -1,    0,    2,   -1,    2, &
          1,    0,    2,   -3,    2, &
          0,    1,    2,   -2,    3, &
          0,    0,    2,   -3,    1, &
         -1,    0,   -2,    2,    1, &
          0,    0,    2,   -4,    2 /
data ( ( nals(i,j), i=1,5 ), j=281,290 ) / &
         -2,    1,    0,    0,    1, &
         -1,    0,    0,   -1,    1, &
          2,    0,    2,   -4,    2, &
          0,    0,    4,   -4,    4, &
          0,    0,    4,   -4,    2, &
         -1,   -2,    0,    2,    1, &
         -2,    0,    0,    3,    0, &
          1,    0,   -2,    2,    1, &
         -3,    0,    2,    2,    2, &
         -3,    0,    2,    2,    1 /
data ( ( nals(i,j), i=1,5 ), j=291,300 ) / &
         -2,    0,    2,    2,    0, &
          2,   -1,    0,    0,    1, &
         -2,    1,    2,    2,    2, &
          1,    1,    0,    1,    0, &
          0,    1,    4,   -2,    2, &
         -1,    1,    0,   -2,    1, &
          0,    0,    0,   -4,    1, &
          1,   -1,    0,    2,    1, &
          1,    1,    0,    2,    1, &
         -1,    2,    2,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j=301,310 ) / &
          3,    1,    2,   -2,    2, &
          0,   -1,    0,    4,    0, &
          2,   -1,    0,    2,    0, &
          0,    0,    4,    0,    1, &
          2,    0,    4,   -2,    2, &
         -1,   -1,    2,    4,    1, &
          1,    0,    0,    4,    1, &
          1,   -2,    2,    2,    2, &
          0,    0,    2,    3,    2, &
         -1,    1,    2,    4,    2 /
data ( ( nals(i,j), i=1,5 ), j=311,320 ) / &
          3,    0,    0,    2,    0, &
         -1,    0,    4,    2,    2, &
          1,    1,    2,    2,    1, &
         -2,    0,    2,    6,    2, &
          2,    1,    2,    2,    2, &
         -1,    0,    2,    6,    2, &
          1,    0,    2,    4,    1, &
          2,    0,    2,    4,    2, &
          1,    1,   -2,    1,    0, &
         -3,    1,    2,    1,    2 /
data ( ( nals(i,j), i=1,5 ), j=321,330 ) / &
          2,    0,   -2,    0,    2, &
         -1,    0,    0,    1,    2, &
         -4,    0,    2,    2,    1, &
         -1,   -1,    0,    1,    0, &
          0,    0,   -2,    2,    2, &
          1,    0,    0,   -1,    2, &
          0,   -1,    2,   -2,    3, &
         -2,    1,    2,    0,    0, &
          0,    0,    2,   -2,    4, &
         -2,   -2,    0,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=331,340 ) / &
         -2,    0,   -2,    4,    0, &
          0,   -2,   -2,    2,    0, &
          1,    2,    0,   -2,    1, &
          3,    0,    0,   -4,    1, &
         -1,    1,    2,   -2,    2, &
          1,   -1,    2,   -4,    1, &
          1,    1,    0,   -2,    2, &
         -3,    0,    2,    0,    0, &
         -3,    0,    2,    0,    2, &
         -2,    0,    0,    1,    0 /
data ( ( nals(i,j), i=1,5 ), j=341,350 ) / &
          0,    0,   -2,    1,    0, &
         -3,    0,    0,    2,    1, &
         -1,   -1,   -2,    2,    0, &
          0,    1,    2,   -4,    1, &
          2,    1,    0,   -4,    1, &
          0,    2,    0,   -2,    1, &
          1,    0,    0,   -3,    1, &
         -2,    0,    2,   -2,    2, &
         -2,   -1,    0,    0,    1, &
         -4,    0,    0,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=351,360 ) / &
          1,    1,    0,   -4,    1, &
         -1,    0,    2,   -4,    1, &
          0,    0,    4,   -4,    1, &
          0,    3,    2,   -2,    2, &
         -3,   -1,    0,    4,    0, &
         -3,    0,    0,    4,    1, &
          1,   -1,   -2,    2,    0, &
         -1,   -1,    0,    2,    2, &
          1,   -2,    0,    0,    1, &
          1,   -1,    0,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=361,370 ) / &
          0,    0,    0,    1,    2, &
         -1,   -1,    2,    0,    0, &
          1,   -2,    2,   -2,    2, &
          0,   -1,    2,   -1,    1, &
         -1,    0,    2,    0,    3, &
          1,    1,    0,    0,    2, &
         -1,    1,    2,    0,    0, &
          1,    2,    0,    0,    0, &
         -1,    2,    2,    0,    2, &
         -1,    0,    4,   -2,    1 /
data ( ( nals(i,j), i=1,5 ), j=371,380 ) / &
          3,    0,    2,   -4,    2, &
          1,    2,    2,   -2,    1, &
          1,    0,    4,   -4,    2, &
         -2,   -1,    0,    4,    1, &
          0,   -1,    0,    2,    2, &
         -2,    1,    0,    4,    0, &
         -2,   -1,    2,    2,    1, &
          2,    0,   -2,    2,    0, &
          1,    0,    0,    1,    1, &
          0,    1,    0,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j=381,390 ) / &
          1,   -1,    2,   -1,    2, &
         -2,    0,    4,    0,    1, &
          2,    1,    0,    0,    1, &
          0,    1,    2,    0,    0, &
          0,   -1,    4,   -2,    2, &
          0,    0,    4,   -2,    4, &
          0,    2,    2,    0,    1, &
         -3,    0,    0,    6,    0, &
         -1,   -1,    0,    4,    1, &
          1,   -2,    0,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=391,400 ) / &
         -1,    0,    0,    4,    2, &
         -1,   -2,    2,    2,    1, &
         -1,    0,    0,   -2,    2, &
          1,    0,   -2,   -2,    1, &
          0,    0,   -2,   -2,    1, &
         -2,    0,   -2,    0,    1, &
          0,    0,    0,    3,    1, &
          0,    0,    0,    3,    0, &
         -1,    1,    0,    4,    0, &
         -1,   -1,    2,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=401,410 ) / &
         -2,    0,    2,    3,    2, &
          1,    0,    0,    2,    2, &
          0,   -1,    2,    1,    2, &
          3,   -1,    0,    0,    0, &
          2,    0,    0,    1,    0, &
          1,   -1,    2,    0,    0, &
          0,    0,    2,    1,    0, &
          1,    0,    2,    0,    3, &
          3,    1,    0,    0,    0, &
          3,   -1,    2,   -2,    2 /
data ( ( nals(i,j), i=1,5 ), j=411,420 ) / &
          2,    0,    2,   -1,    1, &
          1,    1,    2,    0,    0, &
          0,    0,    4,   -1,    2, &
          1,    2,    2,    0,    2, &
         -2,    0,    0,    6,    0, &
          0,   -1,    0,    4,    1, &
         -2,   -1,    2,    4,    1, &
          0,   -2,    2,    2,    1, &
          0,   -1,    2,    2,    0, &
         -1,    0,    2,    3,    1 /
data ( ( nals(i,j), i=1,5 ), j=421,430 ) / &
         -2,    1,    2,    4,    2, &
          2,    0,    0,    2,    2, &
          2,   -2,    2,    0,    2, &
         -1,    1,    2,    3,    2, &
          3,    0,    2,   -1,    2, &
          4,    0,    2,   -2,    1, &
         -1,    0,    0,    6,    0, &
         -1,   -2,    2,    4,    2, &
         -3,    0,    2,    6,    2, &
         -1,    0,    2,    4,    0 /
data ( ( nals(i,j), i=1,5 ), j=431,440 ) / &
          3,    0,    0,    2,    1, &
          3,   -1,    2,    0,    1, &
          3,    0,    2,    0,    0, &
          1,    0,    4,    0,    2, &
          5,    0,    2,   -2,    2, &
          0,   -1,    2,    4,    1, &
          2,   -1,    2,    2,    1, &
          0,    1,    2,    4,    2, &
          1,   -1,    2,    4,    2, &
          3,   -1,    2,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j=441,450 ) / &
          3,    0,    2,    2,    1, &
          5,    0,    2,    0,    2, &
          0,    0,    2,    6,    2, &
          4,    0,    2,    2,    2, &
          0,   -1,    1,   -1,    1, &
         -1,    0,    1,    0,    3, &
          0,   -2,    2,   -2,    3, &
          1,    0,   -1,    0,    1, &
          2,   -2,    0,   -2,    1, &
         -1,    0,    1,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=451,460 ) / &
         -1,    0,    1,    0,    1, &
         -1,   -1,    2,   -1,    2, &
         -2,    2,    0,    2,    2, &
         -1,    0,    1,    0,    0, &
         -4,    1,    2,    2,    2, &
         -3,    0,    2,    1,    1, &
         -2,   -1,    2,    0,    2, &
          1,    0,   -2,    1,    1, &
          2,   -1,   -2,    0,    1, &
         -4,    0,    2,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=461,470 ) / &
         -3,    1,    0,    3,    0, &
         -1,    0,   -1,    2,    0, &
          0,   -2,    0,    0,    2, &
          0,   -2,    0,    0,    2, &
         -3,    0,    0,    3,    0, &
         -2,   -1,    0,    2,    2, &
         -1,    0,   -2,    3,    0, &
         -4,    0,    0,    4,    0, &
          2,    1,   -2,    0,    1, &
          2,   -1,    0,   -2,    2 /
data ( ( nals(i,j), i=1,5 ), j=471,480 ) / &
          0,    0,    1,   -1,    0, &
         -1,    2,    0,    1,    0, &
         -2,    1,    2,    0,    2, &
          1,    1,    0,   -1,    1, &
          1,    0,    1,   -2,    1, &
          0,    2,    0,    0,    2, &
          1,   -1,    2,   -3,    1, &
         -1,    1,    2,   -1,    1, &
         -2,    0,    4,   -2,    2, &
         -2,    0,    4,   -2,    1 /
data ( ( nals(i,j), i=1,5 ), j=481,490 ) / &
         -2,   -2,    0,    2,    1, &
         -2,    0,   -2,    4,    0, &
          1,    2,    2,   -4,    1, &
          1,    1,    2,   -4,    2, &
         -1,    2,    2,   -2,    1, &
          2,    0,    0,   -3,    1, &
         -1,    2,    0,    0,    1, &
          0,    0,    0,   -2,    0, &
         -1,   -1,    2,   -2,    2, &
         -1,    1,    0,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=491,500 ) / &
          0,    0,    0,   -1,    2, &
         -2,    1,    0,    1,    0, &
          1,   -2,    0,   -2,    1, &
          1,    0,   -2,    0,    2, &
         -3,    1,    0,    2,    0, &
         -1,    1,   -2,    2,    0, &
         -1,   -1,    0,    0,    2, &
         -3,    0,    0,    2,    0, &
         -3,   -1,    0,    2,    0, &
          2,    0,    2,   -6,    1 /
data ( ( nals(i,j), i=1,5 ), j=501,510 ) / &
          0,    1,    2,   -4,    2, &
          2,    0,    0,   -4,    2, &
         -2,    1,    2,   -2,    1, &
          0,   -1,    2,   -4,    1, &
          0,    1,    0,   -2,    2, &
         -1,    0,    0,   -2,    0, &
          2,    0,   -2,   -2,    1, &
         -4,    0,    2,    0,    1, &
         -1,   -1,    0,   -1,    1, &
          0,    0,   -2,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=511,520 ) / &
         -3,    0,    0,    1,    0, &
         -1,    0,   -2,    1,    0, &
         -2,    0,   -2,    2,    1, &
          0,    0,   -4,    2,    0, &
         -2,   -1,   -2,    2,    0, &
          1,    0,    2,   -6,    1, &
         -1,    0,    2,   -4,    2, &
          1,    0,    0,   -4,    2, &
          2,    1,    2,   -4,    2, &
          2,    1,    2,   -4,    1 /
data ( ( nals(i,j), i=1,5 ), j=521,530 ) / &
          0,    1,    4,   -4,    4, &
          0,    1,    4,   -4,    2, &
         -1,   -1,   -2,    4,    0, &
         -1,   -3,    0,    2,    0, &
         -1,    0,   -2,    4,    1, &
         -2,   -1,    0,    3,    0, &
          0,    0,   -2,    3,    0, &
         -2,    0,    0,    3,    1, &
          0,   -1,    0,    1,    0, &
         -3,    0,    2,    2,    0 /
data ( ( nals(i,j), i=1,5 ), j=531,540 ) / &
          1,    1,   -2,    2,    0, &
         -1,    1,    0,    2,    2, &
          1,   -2,    2,   -2,    1, &
          0,    0,    1,    0,    2, &
          0,    0,    1,    0,    1, &
          0,    0,    1,    0,    0, &
         -1,    2,    0,    2,    1, &
          0,    0,    2,    0,    2, &
         -2,    0,    2,    0,    2, &
          2,    0,    0,   -1,    1 /
data ( ( nals(i,j), i=1,5 ), j=541,550 ) / &
          3,    0,    0,   -2,    1, &
          1,    0,    2,   -2,    3, &
          1,    2,    0,    0,    1, &
          2,    0,    2,   -3,    2, &
         -1,    1,    4,   -2,    2, &
         -2,   -2,    0,    4,    0, &
          0,   -3,    0,    2,    0, &
          0,    0,   -2,    4,    0, &
         -1,   -1,    0,    3,    0, &
         -2,    0,    0,    4,    2 /
data ( ( nals(i,j), i=1,5 ), j=551,560 ) / &
         -1,    0,    0,    3,    1, &
          2,   -2,    0,    0,    0, &
          1,   -1,    0,    1,    0, &
         -1,    0,    0,    2,    0, &
          0,   -2,    2,    0,    1, &
         -1,    0,    1,    2,    1, &
         -1,    1,    0,    3,    0, &
         -1,   -1,    2,    1,    2, &
          0,   -1,    2,    0,    0, &
         -2,    1,    2,    2,    1 /
data ( ( nals(i,j), i=1,5 ), j=561,570 ) / &
          2,   -2,    2,   -2,    2, &
          1,    1,    0,    1,    1, &
          1,    0,    1,    0,    1, &
          1,    0,    1,    0,    0, &
          0,    2,    0,    2,    0, &
          2,   -1,    2,   -2,    1, &
          0,   -1,    4,   -2,    1, &
          0,    0,    4,   -2,    3, &
          0,    1,    4,   -2,    1, &
          4,    0,    2,   -4,    2 /
data ( ( nals(i,j), i=1,5 ), j=571,580 ) / &
          2,    2,    2,   -2,    2, &
          2,    0,    4,   -4,    2, &
         -1,   -2,    0,    4,    0, &
         -1,   -3,    2,    2,    2, &
         -3,    0,    2,    4,    2, &
         -3,    0,    2,   -2,    1, &
         -1,   -1,    0,   -2,    1, &
         -3,    0,    0,    0,    2, &
         -3,    0,   -2,    2,    0, &
          0,    1,    0,   -4,    1 /
data ( ( nals(i,j), i=1,5 ), j=581,590 ) / &
         -2,    1,    0,   -2,    1, &
         -4,    0,    0,    0,    1, &
         -1,    0,    0,   -4,    1, &
         -3,    0,    0,   -2,    1, &
          0,    0,    0,    3,    2, &
         -1,    1,    0,    4,    1, &
          1,   -2,    2,    0,    1, &
          0,    1,    0,    3,    0, &
         -1,    0,    2,    2,    3, &
          0,    0,    2,    2,    2 /
data ( ( nals(i,j), i=1,5 ), j=591,600 ) / &
         -2,    0,    2,    2,    2, &
         -1,    1,    2,    2,    0, &
          3,    0,    0,    0,    2, &
          2,    1,    0,    1,    0, &
          2,   -1,    2,   -1,    2, &
          0,    0,    2,    0,    1, &
          0,    0,    3,    0,    3, &
          0,    0,    3,    0,    2, &
         -1,    2,    2,    2,    1, &
         -1,    0,    4,    0,    0 /
data ( ( nals(i,j), i=1,5 ), j=601,610 ) / &
          1,    2,    2,    0,    1, &
          3,    1,    2,   -2,    1, &
          1,    1,    4,   -2,    2, &
         -2,   -1,    0,    6,    0, &
          0,   -2,    0,    4,    0, &
         -2,    0,    0,    6,    1, &
         -2,   -2,    2,    4,    2, &
          0,   -3,    2,    2,    2, &
          0,    0,    0,    4,    2, &
         -1,   -1,    2,    3,    2 /
data ( ( nals(i,j), i=1,5 ), j=611,620 ) / &
         -2,    0,    2,    4,    0, &
          2,   -1,    0,    2,    1, &
          1,    0,    0,    3,    0, &
          0,    1,    0,    4,    1, &
          0,    1,    0,    4,    0, &
          1,   -1,    2,    1,    2, &
          0,    0,    2,    2,    3, &
          1,    0,    2,    2,    2, &
         -1,    0,    2,    2,    2, &
         -2,    0,    4,    2,    1 /
data ( ( nals(i,j), i=1,5 ), j=621,630 ) / &
          2,    1,    0,    2,    1, &
          2,    1,    0,    2,    0, &
          2,   -1,    2,    0,    0, &
          1,    0,    2,    1,    0, &
          0,    1,    2,    2,    0, &
          2,    0,    2,    0,    3, &
          3,    0,    2,    0,    2, &
          1,    0,    2,    0,    2, &
          1,    0,    3,    0,    3, &
          1,    1,    2,    1,    1 /
data ( ( nals(i,j), i=1,5 ), j=631,640 ) / &
          0,    2,    2,    2,    2, &
          2,    1,    2,    0,    0, &
          2,    0,    4,   -2,    1, &
          4,    1,    2,   -2,    2, &
         -1,   -1,    0,    6,    0, &
         -3,   -1,    2,    6,    2, &
         -1,    0,    0,    6,    1, &
         -3,    0,    2,    6,    1, &
          1,   -1,    0,    4,    1, &
          1,   -1,    0,    4,    0 /
data ( ( nals(i,j), i=1,5 ), j=641,650 ) / &
         -2,    0,    2,    5,    2, &
          1,   -2,    2,    2,    1, &
          3,   -1,    0,    2,    0, &
          1,   -1,    2,    2,    0, &
          0,    0,    2,    3,    1, &
         -1,    1,    2,    4,    1, &
          0,    1,    2,    3,    2, &
         -1,    0,    4,    2,    1, &
          2,    0,    2,    1,    1, &
          5,    0,    0,    0,    0 /
data ( ( nals(i,j), i=1,5 ), j=651,660 ) / &
          2,    1,    2,    1,    2, &
          1,    0,    4,    0,    1, &
          3,    1,    2,    0,    1, &
          3,    0,    4,   -2,    2, &
         -2,   -1,    2,    6,    2, &
          0,    0,    0,    6,    0, &
          0,   -2,    2,    4,    2, &
         -2,    0,    2,    6,    1, &
          2,    0,    0,    4,    1, &
          2,    0,    0,    4,    0 /
data ( ( nals(i,j), i=1,5 ), j=661,670 ) / &
          2,   -2,    2,    2,    2, &
          0,    0,    2,    4,    0, &
          1,    0,    2,    3,    2, &
          4,    0,    0,    2,    0, &
          2,    0,    2,    2,    0, &
          0,    0,    4,    2,    2, &
          4,   -1,    2,    0,    2, &
          3,    0,    2,    1,    2, &
          2,    1,    2,    2,    1, &
          4,    1,    2,    0,    2 /
data ( ( nals(i,j), i=1,5 ), j=671,678 ) / &
         -1,   -1,    2,    6,    2, &
         -1,    0,    2,    6,    1, &
          1,   -1,    2,    4,    1, &
          1,    1,    2,    4,    2, &
          3,    1,    2,    2,    2, &
          5,    0,    2,    0,    1, &
          2,   -1,    2,    4,    2, &
          2,    0,    2,    4,    1 /

!
!  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
!  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
!

data ( ( cls(i,j), i=1,6 ), j=  1, 10 ) / &
 -172064161d0, -174666d0,  33386d0, 92052331d0,  9086d0, 15377d0, &
  -13170906d0,   -1675d0, -13696d0,  5730336d0, -3015d0, -4587d0, &
   -2276413d0,    -234d0,   2796d0,   978459d0,  -485d0,  1374d0, &
    2074554d0,     207d0,   -698d0,  -897492d0,   470d0,  -291d0, &
    1475877d0,   -3633d0,  11817d0,    73871d0,  -184d0, -1924d0, &
    -516821d0,    1226d0,   -524d0,   224386d0,  -677d0,  -174d0, &
     711159d0,      73d0,   -872d0,    -6750d0,     0d0,   358d0, &
    -387298d0,    -367d0,    380d0,   200728d0,    18d0,   318d0, &
    -301461d0,     -36d0,    816d0,   129025d0,   -63d0,   367d0, &
     215829d0,    -494d0,    111d0,   -95929d0,   299d0,   132d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 11, 20 ) / &
     128227d0,     137d0,    181d0,   -68982d0,    -9d0,    39d0, &
     123457d0,      11d0,     19d0,   -53311d0,    32d0,    -4d0, &
     156994d0,      10d0,   -168d0,    -1235d0,     0d0,    82d0, &
      63110d0,      63d0,     27d0,   -33228d0,     0d0,    -9d0, &
     -57976d0,     -63d0,   -189d0,    31429d0,     0d0,   -75d0, &
     -59641d0,     -11d0,    149d0,    25543d0,   -11d0,    66d0, &
     -51613d0,     -42d0,    129d0,    26366d0,     0d0,    78d0, &
      45893d0,      50d0,     31d0,   -24236d0,   -10d0,    20d0, &
      63384d0,      11d0,   -150d0,    -1220d0,     0d0,    29d0, &
     -38571d0,      -1d0,    158d0,    16452d0,   -11d0,    68d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 21, 30 ) / &
      32481d0,       0d0,      0d0,   -13870d0,     0d0,     0d0, &
     -47722d0,       0d0,    -18d0,      477d0,     0d0,   -25d0, &
     -31046d0,      -1d0,    131d0,    13238d0,   -11d0,    59d0, &
      28593d0,       0d0,     -1d0,   -12338d0,    10d0,    -3d0, &
      20441d0,      21d0,     10d0,   -10758d0,     0d0,    -3d0, &
      29243d0,       0d0,    -74d0,     -609d0,     0d0,    13d0, &
      25887d0,       0d0,    -66d0,     -550d0,     0d0,    11d0, &
     -14053d0,     -25d0,     79d0,     8551d0,    -2d0,   -45d0, &
      15164d0,      10d0,     11d0,    -8001d0,     0d0,    -1d0, &
     -15794d0,      72d0,    -16d0,     6850d0,   -42d0,    -5d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 31, 40 ) / &
      21783d0,       0d0,     13d0,     -167d0,     0d0,    13d0, &
     -12873d0,     -10d0,    -37d0,     6953d0,     0d0,   -14d0, &
     -12654d0,      11d0,     63d0,     6415d0,     0d0,    26d0, &
     -10204d0,       0d0,     25d0,     5222d0,     0d0,    15d0, &
      16707d0,     -85d0,    -10d0,      168d0,    -1d0,    10d0, &
      -7691d0,       0d0,     44d0,     3268d0,     0d0,    19d0, &
     -11024d0,       0d0,    -14d0,      104d0,     0d0,     2d0, &
       7566d0,     -21d0,    -11d0,    -3250d0,     0d0,    -5d0, &
      -6637d0,     -11d0,     25d0,     3353d0,     0d0,    14d0, &
      -7141d0,      21d0,      8d0,     3070d0,     0d0,     4d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 41, 50 ) / &
      -6302d0,     -11d0,      2d0,     3272d0,     0d0,     4d0, &
       5800d0,      10d0,      2d0,    -3045d0,     0d0,    -1d0, &
       6443d0,       0d0,     -7d0,    -2768d0,     0d0,    -4d0, &
      -5774d0,     -11d0,    -15d0,     3041d0,     0d0,    -5d0, &
      -5350d0,       0d0,     21d0,     2695d0,     0d0,    12d0, &
      -4752d0,     -11d0,     -3d0,     2719d0,     0d0,    -3d0, &
      -4940d0,     -11d0,    -21d0,     2720d0,     0d0,    -9d0, &
       7350d0,       0d0,     -8d0,      -51d0,     0d0,     4d0, &
       4065d0,       0d0,      6d0,    -2206d0,     0d0,     1d0, &
       6579d0,       0d0,    -24d0,     -199d0,     0d0,     2d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 51, 60 ) / &
       3579d0,       0d0,      5d0,    -1900d0,     0d0,     1d0, &
       4725d0,       0d0,     -6d0,      -41d0,     0d0,     3d0, &
      -3075d0,       0d0,     -2d0,     1313d0,     0d0,    -1d0, &
      -2904d0,       0d0,     15d0,     1233d0,     0d0,     7d0, &
       4348d0,       0d0,    -10d0,      -81d0,     0d0,     2d0, &
      -2878d0,       0d0,      8d0,     1232d0,     0d0,     4d0, &
      -4230d0,       0d0,      5d0,      -20d0,     0d0,    -2d0, &
      -2819d0,       0d0,      7d0,     1207d0,     0d0,     3d0, &
      -4056d0,       0d0,      5d0,       40d0,     0d0,    -2d0, &
      -2647d0,       0d0,     11d0,     1129d0,     0d0,     5d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 61, 70 ) / &
      -2294d0,       0d0,    -10d0,     1266d0,     0d0,    -4d0, &
       2481d0,       0d0,     -7d0,    -1062d0,     0d0,    -3d0, &
       2179d0,       0d0,     -2d0,    -1129d0,     0d0,    -2d0, &
       3276d0,       0d0,      1d0,       -9d0,     0d0,     0d0, &
      -3389d0,       0d0,      5d0,       35d0,     0d0,    -2d0, &
       3339d0,       0d0,    -13d0,     -107d0,     0d0,     1d0, &
      -1987d0,       0d0,     -6d0,     1073d0,     0d0,    -2d0, &
      -1981d0,       0d0,      0d0,      854d0,     0d0,     0d0, &
       4026d0,       0d0,   -353d0,     -553d0,     0d0,  -139d0, &
       1660d0,       0d0,     -5d0,     -710d0,     0d0,    -2d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 71, 80 ) / &
      -1521d0,       0d0,      9d0,      647d0,     0d0,     4d0, &
       1314d0,       0d0,      0d0,     -700d0,     0d0,     0d0, &
      -1283d0,       0d0,      0d0,      672d0,     0d0,     0d0, &
      -1331d0,       0d0,      8d0,      663d0,     0d0,     4d0, &
       1383d0,       0d0,     -2d0,     -594d0,     0d0,    -2d0, &
       1405d0,       0d0,      4d0,     -610d0,     0d0,     2d0, &
       1290d0,       0d0,      0d0,     -556d0,     0d0,     0d0, &
      -1214d0,       0d0,      5d0,      518d0,     0d0,     2d0, &
       1146d0,       0d0,     -3d0,     -490d0,     0d0,    -1d0, &
       1019d0,       0d0,     -1d0,     -527d0,     0d0,    -1d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 81, 90 ) / &
      -1100d0,       0d0,      9d0,      465d0,     0d0,     4d0, &
       -970d0,       0d0,      2d0,      496d0,     0d0,     1d0, &
       1575d0,       0d0,     -6d0,      -50d0,     0d0,     0d0, &
        934d0,       0d0,     -3d0,     -399d0,     0d0,    -1d0, &
        922d0,       0d0,     -1d0,     -395d0,     0d0,    -1d0, &
        815d0,       0d0,     -1d0,     -422d0,     0d0,    -1d0, &
        834d0,       0d0,      2d0,     -440d0,     0d0,     1d0, &
       1248d0,       0d0,      0d0,     -170d0,     0d0,     1d0, &
       1338d0,       0d0,     -5d0,      -39d0,     0d0,     0d0, &
        716d0,       0d0,     -2d0,     -389d0,     0d0,    -1d0 /      !
data ( ( cls(i,j), i=1,6 ), j= 91,100 ) / &
       1282d0,       0d0,     -3d0,      -23d0,     0d0,     1d0, &
        742d0,       0d0,      1d0,     -391d0,     0d0,     0d0, &
       1020d0,       0d0,    -25d0,     -495d0,     0d0,   -10d0, &
        715d0,       0d0,     -4d0,     -326d0,     0d0,     2d0, &
       -666d0,       0d0,     -3d0,      369d0,     0d0,    -1d0, &
       -667d0,       0d0,      1d0,      346d0,     0d0,     1d0, &
       -704d0,       0d0,      0d0,      304d0,     0d0,     0d0, &
       -694d0,       0d0,      5d0,      294d0,     0d0,     2d0, &
      -1014d0,       0d0,     -1d0,        4d0,     0d0,    -1d0, &
       -585d0,       0d0,     -2d0,      316d0,     0d0,    -1d0 /      !
data ( ( cls(i,j), i=1,6 ), j=101,110 ) / &
       -949d0,       0d0,      1d0,        8d0,     0d0,    -1d0, &
       -595d0,       0d0,      0d0,      258d0,     0d0,     0d0, &
        528d0,       0d0,      0d0,     -279d0,     0d0,     0d0, &
       -590d0,       0d0,      4d0,      252d0,     0d0,     2d0, &
        570d0,       0d0,     -2d0,     -244d0,     0d0,    -1d0, &
       -502d0,       0d0,      3d0,      250d0,     0d0,     2d0, &
       -875d0,       0d0,      1d0,       29d0,     0d0,     0d0, &
       -492d0,       0d0,     -3d0,      275d0,     0d0,    -1d0, &
        535d0,       0d0,     -2d0,     -228d0,     0d0,    -1d0, &
       -467d0,       0d0,      1d0,      240d0,     0d0,     1d0 /      !
data ( ( cls(i,j), i=1,6 ), j=111,120 ) / &
        591d0,       0d0,      0d0,     -253d0,     0d0,     0d0, &
       -453d0,       0d0,     -1d0,      244d0,     0d0,    -1d0, &
        766d0,       0d0,      1d0,        9d0,     0d0,     0d0, &
       -446d0,       0d0,      2d0,      225d0,     0d0,     1d0, &
       -488d0,       0d0,      2d0,      207d0,     0d0,     1d0, &
       -468d0,       0d0,      0d0,      201d0,     0d0,     0d0, &
       -421d0,       0d0,      1d0,      216d0,     0d0,     1d0, &
        463d0,       0d0,      0d0,     -200d0,     0d0,     0d0, &
       -673d0,       0d0,      2d0,       14d0,     0d0,     0d0, &
        658d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=121,130 ) / &
       -438d0,       0d0,      0d0,      188d0,     0d0,     0d0, &
       -390d0,       0d0,      0d0,      205d0,     0d0,     0d0, &
        639d0,     -11d0,     -2d0,      -19d0,     0d0,     0d0, &
        412d0,       0d0,     -2d0,     -176d0,     0d0,    -1d0, &
       -361d0,       0d0,      0d0,      189d0,     0d0,     0d0, &
        360d0,       0d0,     -1d0,     -185d0,     0d0,    -1d0, &
        588d0,       0d0,     -3d0,      -24d0,     0d0,     0d0, &
       -578d0,       0d0,      1d0,        5d0,     0d0,     0d0, &
       -396d0,       0d0,      0d0,      171d0,     0d0,     0d0, &
        565d0,       0d0,     -1d0,       -6d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=131,140 ) / &
       -335d0,       0d0,     -1d0,      184d0,     0d0,    -1d0, &
        357d0,       0d0,      1d0,     -154d0,     0d0,     0d0, &
        321d0,       0d0,      1d0,     -174d0,     0d0,     0d0, &
       -301d0,       0d0,     -1d0,      162d0,     0d0,     0d0, &
       -334d0,       0d0,      0d0,      144d0,     0d0,     0d0, &
        493d0,       0d0,     -2d0,      -15d0,     0d0,     0d0, &
        494d0,       0d0,     -2d0,      -19d0,     0d0,     0d0, &
        337d0,       0d0,     -1d0,     -143d0,     0d0,    -1d0, &
        280d0,       0d0,     -1d0,     -144d0,     0d0,     0d0, &
        309d0,       0d0,      1d0,     -134d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=141,150 ) / &
       -263d0,       0d0,      2d0,      131d0,     0d0,     1d0, &
        253d0,       0d0,      1d0,     -138d0,     0d0,     0d0, &
        245d0,       0d0,      0d0,     -128d0,     0d0,     0d0, &
        416d0,       0d0,     -2d0,      -17d0,     0d0,     0d0, &
       -229d0,       0d0,      0d0,      128d0,     0d0,     0d0, &
        231d0,       0d0,      0d0,     -120d0,     0d0,     0d0, &
       -259d0,       0d0,      2d0,      109d0,     0d0,     1d0, &
        375d0,       0d0,     -1d0,       -8d0,     0d0,     0d0, &
        252d0,       0d0,      0d0,     -108d0,     0d0,     0d0, &
       -245d0,       0d0,      1d0,      104d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=151,160 ) / &
        243d0,       0d0,     -1d0,     -104d0,     0d0,     0d0, &
        208d0,       0d0,      1d0,     -112d0,     0d0,     0d0, &
        199d0,       0d0,      0d0,     -102d0,     0d0,     0d0, &
       -208d0,       0d0,      1d0,      105d0,     0d0,     0d0, &
        335d0,       0d0,     -2d0,      -14d0,     0d0,     0d0, &
       -325d0,       0d0,      1d0,        7d0,     0d0,     0d0, &
       -187d0,       0d0,      0d0,       96d0,     0d0,     0d0, &
        197d0,       0d0,     -1d0,     -100d0,     0d0,     0d0, &
       -192d0,       0d0,      2d0,       94d0,     0d0,     1d0, &
       -188d0,       0d0,      0d0,       83d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=161,170 ) / &
        276d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
       -286d0,       0d0,      1d0,        6d0,     0d0,     0d0, &
        186d0,       0d0,     -1d0,      -79d0,     0d0,     0d0, &
       -219d0,       0d0,      0d0,       43d0,     0d0,     0d0, &
        276d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
       -153d0,       0d0,     -1d0,       84d0,     0d0,     0d0, &
       -156d0,       0d0,      0d0,       81d0,     0d0,     0d0, &
       -154d0,       0d0,      1d0,       78d0,     0d0,     0d0, &
       -174d0,       0d0,      1d0,       75d0,     0d0,     0d0, &
       -163d0,       0d0,      2d0,       69d0,     0d0,     1d0 /      !
data ( ( cls(i,j), i=1,6 ), j=171,180 ) / &
       -228d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         91d0,       0d0,     -4d0,      -54d0,     0d0,    -2d0, &
        175d0,       0d0,      0d0,      -75d0,     0d0,     0d0, &
       -159d0,       0d0,      0d0,       69d0,     0d0,     0d0, &
        141d0,       0d0,      0d0,      -72d0,     0d0,     0d0, &
        147d0,       0d0,      0d0,      -75d0,     0d0,     0d0, &
       -132d0,       0d0,      0d0,       69d0,     0d0,     0d0, &
        159d0,       0d0,    -28d0,      -54d0,     0d0,    11d0, &
        213d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
        123d0,       0d0,      0d0,      -64d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=181,190 ) / &
       -118d0,       0d0,     -1d0,       66d0,     0d0,     0d0, &
        144d0,       0d0,     -1d0,      -61d0,     0d0,     0d0, &
       -121d0,       0d0,      1d0,       60d0,     0d0,     0d0, &
       -134d0,       0d0,      1d0,       56d0,     0d0,     1d0, &
       -105d0,       0d0,      0d0,       57d0,     0d0,     0d0, &
       -102d0,       0d0,      0d0,       56d0,     0d0,     0d0, &
        120d0,       0d0,      0d0,      -52d0,     0d0,     0d0, &
        101d0,       0d0,      0d0,      -54d0,     0d0,     0d0, &
       -113d0,       0d0,      0d0,       59d0,     0d0,     0d0, &
       -106d0,       0d0,      0d0,       61d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=191,200 ) / &
       -129d0,       0d0,      1d0,       55d0,     0d0,     0d0, &
       -114d0,       0d0,      0d0,       57d0,     0d0,     0d0, &
        113d0,       0d0,     -1d0,      -49d0,     0d0,     0d0, &
       -102d0,       0d0,      0d0,       44d0,     0d0,     0d0, &
        -94d0,       0d0,      0d0,       51d0,     0d0,     0d0, &
       -100d0,       0d0,     -1d0,       56d0,     0d0,     0d0, &
         87d0,       0d0,      0d0,      -47d0,     0d0,     0d0, &
        161d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         96d0,       0d0,      0d0,      -50d0,     0d0,     0d0, &
        151d0,       0d0,     -1d0,       -5d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=201,210 ) / &
       -104d0,       0d0,      0d0,       44d0,     0d0,     0d0, &
       -110d0,       0d0,      0d0,       48d0,     0d0,     0d0, &
       -100d0,       0d0,      1d0,       50d0,     0d0,     0d0, &
         92d0,       0d0,     -5d0,       12d0,     0d0,    -2d0, &
         82d0,       0d0,      0d0,      -45d0,     0d0,     0d0, &
         82d0,       0d0,      0d0,      -45d0,     0d0,     0d0, &
        -78d0,       0d0,      0d0,       41d0,     0d0,     0d0, &
        -77d0,       0d0,      0d0,       43d0,     0d0,     0d0, &
          2d0,       0d0,      0d0,       54d0,     0d0,     0d0, &
         94d0,       0d0,      0d0,      -40d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=211,220 ) / &
        -93d0,       0d0,      0d0,       40d0,     0d0,     0d0, &
        -83d0,       0d0,     10d0,       40d0,     0d0,    -2d0, &
         83d0,       0d0,      0d0,      -36d0,     0d0,     0d0, &
        -91d0,       0d0,      0d0,       39d0,     0d0,     0d0, &
        128d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
        -79d0,       0d0,      0d0,       34d0,     0d0,     0d0, &
        -83d0,       0d0,      0d0,       47d0,     0d0,     0d0, &
         84d0,       0d0,      0d0,      -44d0,     0d0,     0d0, &
         83d0,       0d0,      0d0,      -43d0,     0d0,     0d0, &
         91d0,       0d0,      0d0,      -39d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=221,230 ) / &
        -77d0,       0d0,      0d0,       39d0,     0d0,     0d0, &
         84d0,       0d0,      0d0,      -43d0,     0d0,     0d0, &
        -92d0,       0d0,      1d0,       39d0,     0d0,     0d0, &
        -92d0,       0d0,      1d0,       39d0,     0d0,     0d0, &
        -94d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         68d0,       0d0,      0d0,      -36d0,     0d0,     0d0, &
        -61d0,       0d0,      0d0,       32d0,     0d0,     0d0, &
         71d0,       0d0,      0d0,      -31d0,     0d0,     0d0, &
         62d0,       0d0,      0d0,      -34d0,     0d0,     0d0, &
        -63d0,       0d0,      0d0,       33d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=231,240 ) / &
        -73d0,       0d0,      0d0,       32d0,     0d0,     0d0, &
        115d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
       -103d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         63d0,       0d0,      0d0,      -28d0,     0d0,     0d0, &
         74d0,       0d0,      0d0,      -32d0,     0d0,     0d0, &
       -103d0,       0d0,     -3d0,        3d0,     0d0,    -1d0, &
        -69d0,       0d0,      0d0,       30d0,     0d0,     0d0, &
         57d0,       0d0,      0d0,      -29d0,     0d0,     0d0, &
         94d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
         64d0,       0d0,      0d0,      -33d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=241,250 ) / &
        -63d0,       0d0,      0d0,       26d0,     0d0,     0d0, &
        -38d0,       0d0,      0d0,       20d0,     0d0,     0d0, &
        -43d0,       0d0,      0d0,       24d0,     0d0,     0d0, &
        -45d0,       0d0,      0d0,       23d0,     0d0,     0d0, &
         47d0,       0d0,      0d0,      -24d0,     0d0,     0d0, &
        -48d0,       0d0,      0d0,       25d0,     0d0,     0d0, &
         45d0,       0d0,      0d0,      -26d0,     0d0,     0d0, &
         56d0,       0d0,      0d0,      -25d0,     0d0,     0d0, &
         88d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
        -75d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=251,260 ) / &
         85d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         49d0,       0d0,      0d0,      -26d0,     0d0,     0d0, &
        -74d0,       0d0,     -3d0,       -1d0,     0d0,    -1d0, &
        -39d0,       0d0,      0d0,       21d0,     0d0,     0d0, &
         45d0,       0d0,      0d0,      -20d0,     0d0,     0d0, &
         51d0,       0d0,      0d0,      -22d0,     0d0,     0d0, &
        -40d0,       0d0,      0d0,       21d0,     0d0,     0d0, &
         41d0,       0d0,      0d0,      -21d0,     0d0,     0d0, &
        -42d0,       0d0,      0d0,       24d0,     0d0,     0d0, &
        -51d0,       0d0,      0d0,       22d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=261,270 ) / &
        -42d0,       0d0,      0d0,       22d0,     0d0,     0d0, &
         39d0,       0d0,      0d0,      -21d0,     0d0,     0d0, &
         46d0,       0d0,      0d0,      -18d0,     0d0,     0d0, &
        -53d0,       0d0,      0d0,       22d0,     0d0,     0d0, &
         82d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
         81d0,       0d0,     -1d0,       -4d0,     0d0,     0d0, &
         47d0,       0d0,      0d0,      -19d0,     0d0,     0d0, &
         53d0,       0d0,      0d0,      -23d0,     0d0,     0d0, &
        -45d0,       0d0,      0d0,       22d0,     0d0,     0d0, &
        -44d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=271,280 ) / &
        -33d0,       0d0,      0d0,       16d0,     0d0,     0d0, &
        -61d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         28d0,       0d0,      0d0,      -15d0,     0d0,     0d0, &
        -38d0,       0d0,      0d0,       19d0,     0d0,     0d0, &
        -33d0,       0d0,      0d0,       21d0,     0d0,     0d0, &
        -60d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         48d0,       0d0,      0d0,      -10d0,     0d0,     0d0, &
         27d0,       0d0,      0d0,      -14d0,     0d0,     0d0, &
         38d0,       0d0,      0d0,      -20d0,     0d0,     0d0, &
         31d0,       0d0,      0d0,      -13d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=281,290 ) / &
        -29d0,       0d0,      0d0,       15d0,     0d0,     0d0, &
         28d0,       0d0,      0d0,      -15d0,     0d0,     0d0, &
        -32d0,       0d0,      0d0,       15d0,     0d0,     0d0, &
         45d0,       0d0,      0d0,       -8d0,     0d0,     0d0, &
        -44d0,       0d0,      0d0,       19d0,     0d0,     0d0, &
         28d0,       0d0,      0d0,      -15d0,     0d0,     0d0, &
        -51d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -36d0,       0d0,      0d0,       20d0,     0d0,     0d0, &
         44d0,       0d0,      0d0,      -19d0,     0d0,     0d0, &
         26d0,       0d0,      0d0,      -14d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=291,300 ) / &
        -60d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         35d0,       0d0,      0d0,      -18d0,     0d0,     0d0, &
        -27d0,       0d0,      0d0,       11d0,     0d0,     0d0, &
         47d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         36d0,       0d0,      0d0,      -15d0,     0d0,     0d0, &
        -36d0,       0d0,      0d0,       20d0,     0d0,     0d0, &
        -35d0,       0d0,      0d0,       19d0,     0d0,     0d0, &
        -37d0,       0d0,      0d0,       19d0,     0d0,     0d0, &
         32d0,       0d0,      0d0,      -16d0,     0d0,     0d0, &
         35d0,       0d0,      0d0,      -14d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=301,310 ) / &
         32d0,       0d0,      0d0,      -13d0,     0d0,     0d0, &
         65d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         47d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         32d0,       0d0,      0d0,      -16d0,     0d0,     0d0, &
         37d0,       0d0,      0d0,      -16d0,     0d0,     0d0, &
        -30d0,       0d0,      0d0,       15d0,     0d0,     0d0, &
        -32d0,       0d0,      0d0,       16d0,     0d0,     0d0, &
        -31d0,       0d0,      0d0,       13d0,     0d0,     0d0, &
         37d0,       0d0,      0d0,      -16d0,     0d0,     0d0, &
         31d0,       0d0,      0d0,      -13d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=311,320 ) / &
         49d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         32d0,       0d0,      0d0,      -13d0,     0d0,     0d0, &
         23d0,       0d0,      0d0,      -12d0,     0d0,     0d0, &
        -43d0,       0d0,      0d0,       18d0,     0d0,     0d0, &
         26d0,       0d0,      0d0,      -11d0,     0d0,     0d0, &
        -32d0,       0d0,      0d0,       14d0,     0d0,     0d0, &
        -29d0,       0d0,      0d0,       14d0,     0d0,     0d0, &
        -27d0,       0d0,      0d0,       12d0,     0d0,     0d0, &
         30d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -11d0,       0d0,      0d0,        5d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=321,330 ) / &
        -21d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
        -34d0,       0d0,      0d0,       15d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
        -36d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -9d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
        -21d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
        -29d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
        -15d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
        -20d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=331,340 ) / &
         28d0,       0d0,      0d0,        0d0,     0d0,    -2d0, &
         17d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -22d0,       0d0,      0d0,       12d0,     0d0,     0d0, &
        -14d0,       0d0,      0d0,        7d0,     0d0,     0d0, &
         24d0,       0d0,      0d0,      -11d0,     0d0,     0d0, &
         11d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
         14d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
         24d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         18d0,       0d0,      0d0,       -8d0,     0d0,     0d0, &
        -38d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=341,350 ) / &
        -31d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -16d0,       0d0,      0d0,        8d0,     0d0,     0d0, &
         29d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -18d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
        -17d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
         16d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
         22d0,       0d0,      0d0,      -12d0,     0d0,     0d0, &
         20d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=351,360 ) / &
        -13d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
        -17d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
        -14d0,       0d0,      0d0,        8d0,     0d0,     0d0, &
          0d0,       0d0,      0d0,       -7d0,     0d0,     0d0, &
         14d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         19d0,       0d0,      0d0,      -10d0,     0d0,     0d0, &
        -34d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -20d0,       0d0,      0d0,        8d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -5d0,     0d0,     0d0, &
        -18d0,       0d0,      0d0,        7d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=361,370 ) / &
         13d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
         17d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
         15d0,       0d0,      0d0,       -8d0,     0d0,     0d0, &
        -11d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         13d0,       0d0,      0d0,       -5d0,     0d0,     0d0, &
        -18d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -35d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
        -19d0,       0d0,      0d0,       10d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=371,380 ) / &
        -26d0,       0d0,      0d0,       11d0,     0d0,     0d0, &
          8d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         10d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
        -21d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
        -15d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -5d0,     0d0,     0d0, &
        -29d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -19d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
         12d0,       0d0,      0d0,       -5d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=381,390 ) / &
         22d0,       0d0,      0d0,       -9d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
        -20d0,       0d0,      0d0,       11d0,     0d0,     0d0, &
        -20d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -17d0,       0d0,      0d0,        7d0,     0d0,     0d0, &
         15d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          8d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
         14d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
         25d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=391,400 ) / &
        -13d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
        -14d0,       0d0,      0d0,        8d0,     0d0,     0d0, &
         13d0,       0d0,      0d0,       -5d0,     0d0,     0d0, &
        -17d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
         10d0,       0d0,      0d0,       -6d0,     0d0,     0d0, &
        -15d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -22d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         28d0,       0d0,      0d0,       -1d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=401,410 ) / &
         15d0,       0d0,      0d0,       -7d0,     0d0,     0d0, &
         23d0,       0d0,      0d0,      -10d0,     0d0,     0d0, &
         12d0,       0d0,      0d0,       -5d0,     0d0,     0d0, &
         29d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
        -25d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         22d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -18d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         15d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
        -23d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         12d0,       0d0,      0d0,       -5d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=411,420 ) / &
         -8d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
        -19d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         21d0,       0d0,      0d0,       -9d0,     0d0,     0d0, &
         23d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
        -16d0,       0d0,      0d0,        8d0,     0d0,     0d0, &
        -19d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
        -22d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
         27d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         16d0,       0d0,      0d0,       -8d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=421,430 ) / &
         19d0,       0d0,      0d0,       -8d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
         -9d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         -9d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         -8d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         18d0,       0d0,      0d0,       -9d0,     0d0,     0d0, &
         16d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
        -10d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
        -23d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
         16d0,       0d0,      0d0,       -1d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=431,440 ) / &
        -12d0,       0d0,      0d0,        6d0,     0d0,     0d0, &
         -8d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         30d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         24d0,       0d0,      0d0,      -10d0,     0d0,     0d0, &
         10d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
        -16d0,       0d0,      0d0,        7d0,     0d0,     0d0, &
        -16d0,       0d0,      0d0,        7d0,     0d0,     0d0, &
         17d0,       0d0,      0d0,       -7d0,     0d0,     0d0, &
        -24d0,       0d0,      0d0,       10d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        5d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=441,450 ) / &
        -24d0,       0d0,      0d0,       11d0,     0d0,     0d0, &
        -23d0,       0d0,      0d0,        9d0,     0d0,     0d0, &
        -13d0,       0d0,      0d0,        5d0,     0d0,     0d0, &
        -15d0,       0d0,      0d0,        7d0,     0d0,     0d0, &
          0d0,       0d0,  -1988d0,        0d0,     0d0, -1679d0, &
          0d0,       0d0,    -63d0,        0d0,     0d0,   -27d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          0d0,       0d0,      5d0,        0d0,     0d0,     4d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          0d0,       0d0,    364d0,        0d0,     0d0,   176d0 /      !
data ( ( cls(i,j), i=1,6 ), j=451,460 ) / &
          0d0,       0d0,  -1044d0,        0d0,     0d0,  -891d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          0d0,       0d0,    330d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=461,470 ) / &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          0d0,       0d0,      5d0,        0d0,     0d0,     0d0, &
          0d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=471,480 ) / &
         -5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
          0d0,       0d0,    -12d0,        0d0,     0d0,   -10d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=481,490 ) / &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          0d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=491,500 ) / &
         -8d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=501,510 ) / &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,       -3d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=511,520 ) / &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          8d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=521,530 ) / &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
        -13d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=531,540 ) / &
         10d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         10d0,       0d0,     13d0,        6d0,     0d0,    -5d0, &
          0d0,       0d0,     30d0,        0d0,     0d0,    14d0, &
          0d0,       0d0,   -162d0,        0d0,     0d0,  -138d0, &
          0d0,       0d0,     75d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=541,550 ) / &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          9d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=551,560 ) / &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -6d0,       0d0,     -3d0,        3d0,     0d0,     1d0, &
          0d0,       0d0,     -3d0,        0d0,     0d0,    -2d0, &
         11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=561,570 ) / &
         -1d0,       0d0,      3d0,        3d0,     0d0,    -1d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          0d0,       0d0,    -13d0,        0d0,     0d0,   -11d0, &
          3d0,       0d0,      6d0,        0d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        3d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=571,580 ) / &
          8d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          8d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=581,590 ) / &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -8d0,       0d0,      0d0,        4d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=591,600 ) / &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          0d0,       0d0,    -26d0,        0d0,     0d0,   -11d0, &
          0d0,       0d0,    -10d0,        0d0,     0d0,    -5d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
        -13d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=601,610 ) / &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -7d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=611,620 ) / &
         13d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
        -11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=621,630 ) / &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
        -12d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
          0d0,       0d0,     -5d0,        0d0,     0d0,    -2d0, &
         -7d0,       0d0,      0d0,        4d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=631,640 ) / &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         12d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=641,650 ) / &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          6d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=651,660 ) / &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,       -4d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -5d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        3d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         10d0,       0d0,      0d0,        0d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=661,670 ) / &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          7d0,       0d0,      0d0,       -3d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
         11d0,       0d0,      0d0,        0d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
         -6d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          5d0,       0d0,      0d0,       -2d0,     0d0,     0d0 /      !
data ( ( cls(i,j), i=1,6 ), j=671,678 ) / &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -4d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0, &
          4d0,       0d0,      0d0,       -2d0,     0d0,     0d0, &
          3d0,       0d0,      0d0,       -1d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        1d0,     0d0,     0d0, &
         -3d0,       0d0,      0d0,        2d0,     0d0,     0d0 /      !

!
!  Planetary argument multipliers:
!              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre

data ( ( napl(i,j), i=1,14 ), j=  1, 10 ) / &
         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2, &
         0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 11, 20 ) / &
         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j= 21, 30 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2, &
         2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0, &
         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 31, 40 ) / &
        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0, &
        -2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
        -1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2, &
         0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 41, 50 ) / &
         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0, &
        -2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0, &
        -1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 51, 60 ) / &
        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0, &
        -2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2, &
        -1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 61, 70 ) / &
         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0, &
         1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0, &
         0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 71, 80 ) / &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0, &
         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0, &
         2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0, &
         0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j= 81, 90 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0, &
         1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j= 91,100 ) / &
        -2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0, &
        -2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0, &
        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
        -1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=101,110 ) / &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0, &
         2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0, &
         0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=111,120 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2, &
         0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2, &
         1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0, &
         1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=121,130 ) / &
         0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0, &
         0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0, &
        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=131,140 ) / &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0, &
        -2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=141,150 ) / &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0, &
         2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0, &
        -2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=151,160 ) / &
         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1, &
        -1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=161,170 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2, &
         0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=171,180 ) / &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0, &
         1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=181,190 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2, &
         0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=191,200 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2, &
        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=201,210 ) / &
         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=211,220 ) / &
         0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=221,230 ) / &
         0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=231,240 ) / &
         0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=241,250 ) / &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2, &
         0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2, &
         0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=251,260 ) / &
         0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2, &
         0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=261,270 ) / &
         0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
        -2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=271,280 ) / &
         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2, &
         0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=281,290 ) / &
         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=291,300 ) / &
         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1, &
         0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=301,310 ) / &
         0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=311,320 ) / &
        -2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=321,330 ) / &
         0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=331,340 ) / &
         0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2, &
        -2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=341,350 ) / &
         0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=351,360 ) / &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=361,370 ) / &
         0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=371,380 ) / &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=381,390 ) / &
         0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2, &
         0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=391,400 ) / &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=401,410 ) / &
         0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=411,420 ) / &
         0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=421,430 ) / &
         0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=431,440 ) / &
         0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=441,450 ) / &
         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=451,460 ) / &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=461,470 ) / &
         0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=471,480 ) / &
         0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=481,490 ) / &
         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=491,500 ) / &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2, &
         0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=501,510 ) / &
         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=511,520 ) / &
         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=521,530 ) / &
         0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=531,540 ) / &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=541,550 ) / &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2, &
         0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=551,560 ) / &
         0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=561,570 ) / &
         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=571,580 ) / &
         0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=581,590 ) / &
         0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=591,600 ) / &
         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=601,610 ) / &
         0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1 /
data ( ( napl(i,j), i=1,14 ), j=611,620 ) / &
         0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 /
data ( ( napl(i,j), i=1,14 ), j=621,630 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, &
         1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
        -1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=631,640 ) / &
        -1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
        -2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
         1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
        -1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0, &
         1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=641,650 ) / &
        -1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
        -1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
        -1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0, &
         1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=651,660 ) / &
         0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=661,670 ) / &
         0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0, &
         1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0, &
         0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=671,680 ) / &
         0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
         2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
data ( ( napl(i,j), i=1,14 ), j=681,687 ) / &
         1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
         2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
        -1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /

!
!  Planetary nutation coefficients, unit 1e-7 arcsec:
!  longitude (sin, cos), obliquity (sin, cos)
!

data ( ( icpl(i,j), i=1,4 ), j=  1, 10 ) / &
       1440,          0,          0,          0, &
         56,       -117,        -42,        -40, &
        125,        -43,          0,        -54, &
          0,          5,          0,          0, &
          3,         -7,         -3,          0, &
          3,          0,          0,         -2, &
       -114,          0,          0,         61, &
       -219,         89,          0,          0, &
         -3,          0,          0,          0, &
       -462,       1604,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j= 11, 20 ) / &
         99,          0,          0,        -53, &
         -3,          0,          0,          2, &
          0,          6,          2,          0, &
          3,          0,          0,          0, &
        -12,          0,          0,          0, &
         14,       -218,        117,          8, &
         31,       -481,       -257,        -17, &
       -491,        128,          0,          0, &
      -3084,       5123,       2735,       1647, &
      -1444,       2409,      -1286,       -771 /
data ( ( icpl(i,j), i=1,4 ), j= 21, 30 ) / &
         11,        -24,        -11,         -9, &
         26,         -9,          0,          0, &
        103,        -60,          0,          0, &
          0,        -13,         -7,          0, &
        -26,        -29,        -16,         14, &
          9,        -27,        -14,         -5, &
         12,          0,          0,         -6, &
         -7,          0,          0,          0, &
          0,         24,          0,          0, &
        284,          0,          0,       -151 /
data ( ( icpl(i,j), i=1,4 ), j= 31, 40 ) / &
        226,        101,          0,          0, &
          0,         -8,         -2,          0, &
          0,         -6,         -3,          0, &
          5,          0,          0,         -3, &
        -41,        175,         76,         17, &
          0,         15,          6,          0, &
        425,        212,       -133,        269, &
       1200,        598,        319,       -641, &
        235,        334,          0,          0, &
         11,        -12,         -7,         -6 /
data ( ( icpl(i,j), i=1,4 ), j= 41, 50 ) / &
          5,         -6,          3,          3, &
         -5,          0,          0,          3, &
          6,          0,          0,         -3, &
         15,          0,          0,          0, &
         13,          0,          0,         -7, &
         -6,         -9,          0,          0, &
        266,        -78,          0,          0, &
       -460,       -435,       -232,        246, &
          0,         15,          7,          0, &
         -3,          0,          0,          2 /
data ( ( icpl(i,j), i=1,4 ), j= 51, 60 ) / &
          0,        131,          0,          0, &
          4,          0,          0,          0, &
          0,          3,          0,          0, &
          0,          4,          2,          0, &
          0,          3,          0,          0, &
        -17,        -19,        -10,          9, &
         -9,        -11,          6,         -5, &
         -6,          0,          0,          3, &
        -16,          8,          0,          0, &
          0,          3,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j= 61, 70 ) / &
         11,         24,         11,         -5, &
         -3,         -4,         -2,          1, &
          3,          0,          0,         -1, &
          0,         -8,         -4,          0, &
          0,          3,          0,          0, &
          0,          5,          0,          0, &
          0,          3,          2,          0, &
         -6,          4,          2,          3, &
         -3,         -5,          0,          0, &
         -5,          0,          0,          2 /
data ( ( icpl(i,j), i=1,4 ), j= 71, 80 ) / &
          4,         24,         13,         -2, &
        -42,         20,          0,          0, &
        -10,        233,          0,          0, &
         -3,          0,          0,          1, &
         78,        -18,          0,          0, &
          0,          3,          1,          0, &
          0,         -3,         -1,          0, &
          0,         -4,         -2,          1, &
          0,         -8,         -4,         -1, &
          0,         -5,          3,          0 /
data ( ( icpl(i,j), i=1,4 ), j= 81, 90 ) / &
         -7,          0,          0,          3, &
        -14,          8,          3,          6, &
          0,          8,         -4,          0, &
          0,         19,         10,          0, &
         45,        -22,          0,          0, &
         -3,          0,          0,          0, &
          0,         -3,          0,          0, &
          0,          3,          0,          0, &
          3,          5,          3,         -2, &
         89,        -16,         -9,        -48 /
data ( ( icpl(i,j), i=1,4 ), j= 91,100 ) / &
          0,          3,          0,          0, &
         -3,          7,          4,          2, &
       -349,        -62,          0,          0, &
        -15,         22,          0,          0, &
         -3,          0,          0,          0, &
        -53,          0,          0,          0, &
          5,          0,          0,         -3, &
          0,         -8,          0,          0, &
         15,         -7,         -4,         -8, &
         -3,          0,          0,          1 /
data ( ( icpl(i,j), i=1,4 ), j=101,110 ) / &
        -21,        -78,          0,          0, &
         20,        -70,        -37,        -11, &
          0,          6,          3,          0, &
          5,          3,          2,         -2, &
        -17,         -4,         -2,          9, &
          0,          6,          3,          0, &
         32,         15,         -8,         17, &
        174,         84,         45,        -93, &
         11,         56,          0,          0, &
        -66,        -12,         -6,         35 /
data ( ( icpl(i,j), i=1,4 ), j=111,120 ) / &
         47,          8,          4,        -25, &
          0,          8,          4,          0, &
         10,        -22,        -12,         -5, &
         -3,          0,          0,          2, &
        -24,         12,          0,          0, &
          5,         -6,          0,          0, &
          3,          0,          0,         -2, &
          4,          3,          1,         -2, &
          0,         29,         15,          0, &
         -5,         -4,         -2,          2 /
data ( ( icpl(i,j), i=1,4 ), j=121,130 ) / &
          8,         -3,         -1,         -5, &
          0,         -3,          0,          0, &
         10,          0,          0,          0, &
          3,          0,          0,         -2, &
         -5,          0,          0,          3, &
         46,         66,         35,        -25, &
        -14,          7,          0,          0, &
          0,          3,          2,          0, &
         -5,          0,          0,          0, &
        -68,        -34,        -18,         36 /
data ( ( icpl(i,j), i=1,4 ), j=131,140 ) / &
          0,         14,          7,          0, &
         10,         -6,         -3,         -5, &
         -5,         -4,         -2,          3, &
         -3,          5,          2,          1, &
         76,         17,          9,        -41, &
         84,        298,        159,        -45, &
          3,          0,          0,         -1, &
         -3,          0,          0,          2, &
         -3,          0,          0,          1, &
        -82,        292,        156,         44 /
data ( ( icpl(i,j), i=1,4 ), j=141,150 ) / &
        -73,         17,          9,         39, &
         -9,        -16,          0,          0, &
          3,          0,         -1,         -2, &
         -3,          0,          0,          0, &
         -9,         -5,         -3,          5, &
       -439,          0,          0,          0, &
         57,        -28,        -15,        -30, &
          0,         -6,         -3,          0, &
         -4,          0,          0,          2, &
        -40,         57,         30,         21 /
data ( ( icpl(i,j), i=1,4 ), j=151,160 ) / &
         23,          7,          3,        -13, &
        273,         80,         43,       -146, &
       -449,        430,          0,          0, &
         -8,        -47,        -25,          4, &
          6,         47,         25,         -3, &
          0,         23,         13,          0, &
         -3,          0,          0,          2, &
          3,         -4,         -2,         -2, &
        -48,       -110,        -59,         26, &
         51,        114,         61,        -27 /
data ( ( icpl(i,j), i=1,4 ), j=161,170 ) / &
       -133,          0,          0,         57, &
          0,          4,          0,          0, &
        -21,         -6,         -3,         11, &
          0,         -3,         -1,          0, &
        -11,        -21,        -11,          6, &
        -18,       -436,       -233,          9, &
         35,         -7,          0,          0, &
          0,          5,          3,          0, &
         11,         -3,         -1,         -6, &
         -5,         -3,         -1,          3 /
data ( ( icpl(i,j), i=1,4 ), j=171,180 ) / &
        -53,         -9,         -5,         28, &
          0,          3,          2,          1, &
          4,          0,          0,         -2, &
          0,         -4,          0,          0, &
        -50,        194,        103,         27, &
        -13,         52,         28,          7, &
        -91,        248,          0,          0, &
          6,         49,         26,         -3, &
         -6,        -47,        -25,          3, &
          0,          5,          3,          0 /
data ( ( icpl(i,j), i=1,4 ), j=181,190 ) / &
         52,         23,         10,        -23, &
         -3,          0,          0,          1, &
          0,          5,          3,          0, &
         -4,          0,          0,          0, &
         -4,          8,          3,          2, &
         10,          0,          0,          0, &
          3,          0,          0,         -2, &
          0,          8,          4,          0, &
          0,          8,          4,          1, &
         -4,          0,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=191,200 ) / &
         -4,          0,          0,          0, &
         -8,          4,          2,          4, &
          8,         -4,         -2,         -4, &
          0,         15,          7,          0, &
       -138,          0,          0,          0, &
          0,         -7,         -3,          0, &
          0,         -7,         -3,          0, &
         54,          0,          0,        -29, &
          0,         10,          4,          0, &
         -7,          0,          0,          3 /
data ( ( icpl(i,j), i=1,4 ), j=201,210 ) / &
        -37,         35,         19,         20, &
          0,          4,          0,          0, &
         -4,          9,          0,          0, &
          8,          0,          0,         -4, &
         -9,        -14,         -8,          5, &
         -3,         -9,         -5,          3, &
       -145,         47,          0,          0, &
        -10,         40,         21,          5, &
         11,        -49,        -26,         -7, &
      -2150,          0,          0,        932 /
data ( ( icpl(i,j), i=1,4 ), j=211,220 ) / &
        -12,          0,          0,          5, &
         85,          0,          0,        -37, &
          4,          0,          0,         -2, &
          3,          0,          0,         -2, &
        -86,        153,          0,          0, &
         -6,          9,          5,          3, &
          9,        -13,         -7,         -5, &
         -8,         12,          6,          4, &
        -51,          0,          0,         22, &
        -11,       -268,       -116,          5 /
data ( ( icpl(i,j), i=1,4 ), j=221,230 ) / &
          0,         12,          5,          0, &
          0,          7,          3,          0, &
         31,          6,          3,        -17, &
        140,         27,         14,        -75, &
         57,         11,          6,        -30, &
        -14,        -39,          0,          0, &
          0,         -6,         -2,          0, &
          4,         15,          8,         -2, &
          0,          4,          0,          0, &
         -3,          0,          0,          1 /
data ( ( icpl(i,j), i=1,4 ), j=231,240 ) / &
          0,         11,          5,          0, &
          9,          6,          0,          0, &
         -4,         10,          4,          2, &
          5,          3,          0,          0, &
         16,          0,          0,         -9, &
         -3,          0,          0,          0, &
          0,          3,          2,         -1, &
          7,          0,          0,         -3, &
        -25,         22,          0,          0, &
         42,        223,        119,        -22 /
data ( ( icpl(i,j), i=1,4 ), j=241,250 ) / &
        -27,       -143,        -77,         14, &
          9,         49,         26,         -5, &
      -1166,          0,          0,        505, &
         -5,          0,          0,          2, &
         -6,          0,          0,          3, &
         -8,          0,          1,          4, &
          0,         -4,          0,          0, &
        117,          0,          0,        -63, &
         -4,          8,          4,          2, &
          3,          0,          0,         -2 /
data ( ( icpl(i,j), i=1,4 ), j=251,260 ) / &
         -5,          0,          0,          2, &
          0,         31,          0,          0, &
         -5,          0,          1,          3, &
          4,          0,          0,         -2, &
         -4,          0,          0,          2, &
        -24,        -13,         -6,         10, &
          3,          0,          0,          0, &
          0,        -32,        -17,          0, &
          8,         12,          5,         -3, &
          3,          0,          0,         -1 /
data ( ( icpl(i,j), i=1,4 ), j=261,270 ) / &
          7,         13,          0,          0, &
         -3,         16,          0,          0, &
         50,          0,          0,        -27, &
          0,         -5,         -3,          0, &
         13,          0,          0,          0, &
          0,          5,          3,          1, &
         24,          5,          2,        -11, &
          5,        -11,         -5,         -2, &
         30,         -3,         -2,        -16, &
         18,          0,          0,         -9 /
data ( ( icpl(i,j), i=1,4 ), j=271,280 ) / &
          8,        614,          0,          0, &
          3,         -3,         -1,         -2, &
          6,         17,          9,         -3, &
         -3,         -9,         -5,          2, &
          0,          6,          3,         -1, &
       -127,         21,          9,         55, &
          3,          5,          0,          0, &
         -6,        -10,         -4,          3, &
          5,          0,          0,          0, &
         16,          9,          4,         -7 /
data ( ( icpl(i,j), i=1,4 ), j=281,290 ) / &
          3,          0,          0,         -2, &
          0,         22,          0,          0, &
          0,         19,         10,          0, &
          7,          0,          0,         -4, &
          0,         -5,         -2,          0, &
          0,          3,          1,          0, &
         -9,          3,          1,          4, &
         17,          0,          0,         -7, &
          0,         -3,         -2,         -1, &
        -20,         34,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=291,300 ) / &
        -10,          0,          1,          5, &
         -4,          0,          0,          2, &
         22,        -87,          0,          0, &
         -4,          0,          0,          2, &
         -3,         -6,         -2,          1, &
        -16,         -3,         -1,          7, &
          0,         -3,         -2,          0, &
          4,          0,          0,          0, &
        -68,         39,          0,          0, &
         27,          0,          0,        -14 /
data ( ( icpl(i,j), i=1,4 ), j=301,310 ) / &
          0,         -4,          0,          0, &
        -25,          0,          0,          0, &
        -12,         -3,         -2,          6, &
          3,          0,          0,         -1, &
          3,         66,         29,         -1, &
        490,          0,          0,       -213, &
        -22,         93,         49,         12, &
         -7,         28,         15,          4, &
         -3,         13,          7,          2, &
        -46,         14,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=311,320 ) / &
         -5,          0,          0,          0, &
          2,          1,          0,          0, &
          0,         -3,          0,          0, &
        -28,          0,          0,         15, &
          5,          0,          0,         -2, &
          0,          3,          0,          0, &
        -11,          0,          0,          5, &
          0,          3,          1,          0, &
         -3,          0,          0,          1, &
         25,        106,         57,        -13 /
data ( ( icpl(i,j), i=1,4 ), j=321,330 ) / &
          5,         21,         11,         -3, &
       1485,          0,          0,          0, &
         -7,        -32,        -17,          4, &
          0,          5,          3,          0, &
         -6,         -3,         -2,          3, &
         30,         -6,         -2,        -13, &
         -4,          4,          0,          0, &
        -19,          0,          0,         10, &
          0,          4,          2,         -1, &
          0,          3,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=331,340 ) / &
          4,          0,          0,         -2, &
          0,         -3,         -1,          0, &
         -3,          0,          0,          0, &
          5,          3,          1,         -2, &
          0,         11,          0,          0, &
        118,          0,          0,        -52, &
          0,         -5,         -3,          0, &
        -28,         36,          0,          0, &
          5,         -5,          0,          0, &
         14,        -59,        -31,         -8 /
data ( ( icpl(i,j), i=1,4 ), j=341,350 ) / &
          0,          9,          5,          1, &
       -458,          0,          0,        198, &
          0,        -45,        -20,          0, &
          9,          0,          0,         -5, &
          0,         -3,          0,          0, &
          0,         -4,         -2,         -1, &
         11,          0,          0,         -6, &
          6,          0,          0,         -2, &
        -16,         23,          0,          0, &
          0,         -4,         -2,          0 /
data ( ( icpl(i,j), i=1,4 ), j=351,360 ) / &
         -5,          0,          0,          2, &
       -166,        269,          0,          0, &
         15,          0,          0,         -8, &
         10,          0,          0,         -4, &
        -78,         45,          0,          0, &
          0,         -5,         -2,          0, &
          7,          0,          0,         -4, &
         -5,        328,          0,          0, &
          3,          0,          0,         -2, &
          5,          0,          0,         -2 /
data ( ( icpl(i,j), i=1,4 ), j=361,370 ) / &
          0,          3,          1,          0, &
         -3,          0,          0,          0, &
         -3,          0,          0,          0, &
          0,         -4,         -2,          0, &
      -1223,        -26,          0,          0, &
          0,          7,          3,          0, &
          3,          0,          0,          0, &
          0,          3,          2,          0, &
         -6,         20,          0,          0, &
       -368,          0,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=371,380 ) / &
        -75,          0,          0,          0, &
         11,          0,          0,         -6, &
          3,          0,          0,         -2, &
         -3,          0,          0,          1, &
        -13,        -30,          0,          0, &
         21,          3,          0,          0, &
         -3,          0,          0,          1, &
         -4,          0,          0,          2, &
          8,        -27,          0,          0, &
        -19,        -11,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=381,390 ) / &
         -4,          0,          0,          2, &
          0,          5,          2,          0, &
         -6,          0,          0,          2, &
         -8,          0,          0,          0, &
         -1,          0,          0,          0, &
        -14,          0,          0,          6, &
          6,          0,          0,          0, &
        -74,          0,          0,         32, &
          0,         -3,         -1,          0, &
          4,          0,          0,         -2 /
data ( ( icpl(i,j), i=1,4 ), j=391,400 ) / &
          8,         11,          0,          0, &
          0,          3,          2,          0, &
       -262,          0,          0,        114, &
          0,         -4,          0,          0, &
         -7,          0,          0,          4, &
          0,        -27,        -12,          0, &
        -19,         -8,         -4,          8, &
        202,          0,          0,        -87, &
         -8,         35,         19,          5, &
          0,          4,          2,          0 /
data ( ( icpl(i,j), i=1,4 ), j=401,410 ) / &
         16,         -5,          0,          0, &
          5,          0,          0,         -3, &
          0,         -3,          0,          0, &
          1,          0,          0,          0, &
        -35,        -48,        -21,         15, &
         -3,         -5,         -2,          1, &
          6,          0,          0,         -3, &
          3,          0,          0,         -1, &
          0,         -5,          0,          0, &
         12,         55,         29,         -6 /
data ( ( icpl(i,j), i=1,4 ), j=411,420 ) / &
          0,          5,          3,          0, &
       -598,          0,          0,          0, &
         -3,        -13,         -7,          1, &
         -5,         -7,         -3,          2, &
          3,          0,          0,         -1, &
          5,         -7,          0,          0, &
          4,          0,          0,         -2, &
         16,         -6,          0,          0, &
          8,         -3,          0,          0, &
          8,        -31,        -16,         -4 /
data ( ( icpl(i,j), i=1,4 ), j=421,430 ) / &
          0,          3,          1,          0, &
        113,          0,          0,        -49, &
          0,        -24,        -10,          0, &
          4,          0,          0,         -2, &
         27,          0,          0,          0, &
         -3,          0,          0,          1, &
          0,         -4,         -2,          0, &
          5,          0,          0,         -2, &
          0,         -3,          0,          0, &
        -13,          0,          0,          6 /
data ( ( icpl(i,j), i=1,4 ), j=431,440 ) / &
          5,          0,          0,         -2, &
        -18,        -10,         -4,          8, &
         -4,        -28,          0,          0, &
         -5,          6,          3,          2, &
         -3,          0,          0,          1, &
         -5,         -9,         -4,          2, &
         17,          0,          0,         -7, &
         11,          4,          0,          0, &
          0,         -6,         -2,          0, &
         83,         15,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=441,450 ) / &
         -4,          0,          0,          2, &
          0,       -114,        -49,          0, &
        117,          0,          0,        -51, &
         -5,         19,         10,          2, &
         -3,          0,          0,          0, &
         -3,          0,          0,          2, &
          0,         -3,         -1,          0, &
          3,          0,          0,          0, &
          0,         -6,         -2,          0, &
        393,          3,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=451,460 ) / &
         -4,         21,         11,          2, &
         -6,          0,         -1,          3, &
         -3,          8,          4,          1, &
          8,          0,          0,          0, &
         18,        -29,        -13,         -8, &
          8,         34,         18,         -4, &
         89,          0,          0,          0, &
          3,         12,          6,         -1, &
         54,        -15,         -7,        -24, &
          0,          3,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=461,470 ) / &
          3,          0,          0,         -1, &
          0,         35,          0,          0, &
       -154,        -30,        -13,         67, &
         15,          0,          0,          0, &
          0,          4,          2,          0, &
          0,          9,          0,          0, &
         80,        -71,        -31,        -35, &
          0,        -20,         -9,          0, &
         11,          5,          2,         -5, &
         61,        -96,        -42,        -27 /
data ( ( icpl(i,j), i=1,4 ), j=471,480 ) / &
         14,          9,          4,         -6, &
        -11,         -6,         -3,          5, &
          0,         -3,         -1,          0, &
        123,       -415,       -180,        -53, &
          0,          0,          0,        -35, &
         -5,          0,          0,          0, &
          7,        -32,        -17,         -4, &
          0,         -9,         -5,          0, &
          0,         -4,          2,          0, &
        -89,          0,          0,         38 /
data ( ( icpl(i,j), i=1,4 ), j=481,490 ) / &
          0,        -86,        -19,         -6, &
          0,          0,        -19,          6, &
       -123,       -416,       -180,         53, &
          0,         -3,         -1,          0, &
         12,         -6,         -3,         -5, &
        -13,          9,          4,          6, &
          0,        -15,         -7,          0, &
          3,          0,          0,         -1, &
        -62,        -97,        -42,         27, &
        -11,          5,          2,          5 /
data ( ( icpl(i,j), i=1,4 ), j=491,500 ) / &
          0,        -19,         -8,          0, &
         -3,          0,          0,          1, &
          0,          4,          2,          0, &
          0,          3,          0,          0, &
          0,          4,          2,          0, &
        -85,        -70,        -31,         37, &
        163,        -12,         -5,        -72, &
        -63,        -16,         -7,         28, &
        -21,        -32,        -14,          9, &
          0,         -3,         -1,          0 /
data ( ( icpl(i,j), i=1,4 ), j=501,510 ) / &
          3,          0,          0,         -2, &
          0,          8,          0,          0, &
          3,         10,          4,         -1, &
          3,          0,          0,         -1, &
          0,         -7,         -3,          0, &
          0,         -4,         -2,          0, &
          6,         19,          0,          0, &
          5,       -173,        -75,         -2, &
          0,         -7,         -3,          0, &
          7,        -12,         -5,         -3 /
data ( ( icpl(i,j), i=1,4 ), j=511,520 ) / &
         -3,          0,          0,          2, &
          3,         -4,         -2,         -1, &
         74,          0,          0,        -32, &
         -3,         12,          6,          2, &
         26,        -14,         -6,        -11, &
         19,          0,          0,         -8, &
          6,         24,         13,         -3, &
         83,          0,          0,          0, &
          0,        -10,         -5,          0, &
         11,         -3,         -1,         -5 /
data ( ( icpl(i,j), i=1,4 ), j=521,530 ) / &
          3,          0,          1,         -1, &
          3,          0,          0,         -1, &
         -4,          0,          0,          0, &
          5,        -23,        -12,         -3, &
       -339,          0,          0,        147, &
          0,        -10,         -5,          0, &
          5,          0,          0,          0, &
          3,          0,          0,         -1, &
          0,         -4,         -2,          0, &
         18,         -3,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=531,540 ) / &
          9,        -11,         -5,         -4, &
         -8,          0,          0,          4, &
          3,          0,          0,         -1, &
          0,          9,          0,          0, &
          6,         -9,         -4,         -2, &
         -4,        -12,          0,          0, &
         67,        -91,        -39,        -29, &
         30,        -18,         -8,        -13, &
          0,          0,          0,          0, &
          0,       -114,        -50,          0 /
data ( ( icpl(i,j), i=1,4 ), j=541,550 ) / &
          0,          0,          0,         23, &
        517,         16,          7,       -224, &
          0,         -7,         -3,          0, &
        143,         -3,         -1,        -62, &
         29,          0,          0,        -13, &
         -4,          0,          0,          2, &
         -6,          0,          0,          3, &
          5,         12,          5,         -2, &
        -25,          0,          0,         11, &
         -3,          0,          0,          1 /
data ( ( icpl(i,j), i=1,4 ), j=551,560 ) / &
          0,          4,          2,          0, &
        -22,         12,          5,         10, &
         50,          0,          0,        -22, &
          0,          7,          4,          0, &
          0,          3,          1,          0, &
         -4,          4,          2,          2, &
         -5,        -11,         -5,          2, &
          0,          4,          2,          0, &
          4,         17,          9,         -2, &
         59,          0,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=561,570 ) / &
          0,         -4,         -2,          0, &
         -8,          0,          0,          4, &
         -3,          0,          0,          0, &
          4,        -15,         -8,         -2, &
        370,         -8,          0,       -160, &
          0,          0,         -3,          0, &
          0,          3,          1,          0, &
         -6,          3,          1,          3, &
          0,          6,          0,          0, &
        -10,          0,          0,          4 /
data ( ( icpl(i,j), i=1,4 ), j=571,580 ) / &
          0,          9,          4,          0, &
          4,         17,          7,         -2, &
         34,          0,          0,        -15, &
          0,          5,          3,          0, &
         -5,          0,          0,          2, &
        -37,         -7,         -3,         16, &
          3,         13,          7,         -2, &
         40,          0,          0,          0, &
          0,         -3,         -2,          0, &
       -184,         -3,         -1,         80 /
data ( ( icpl(i,j), i=1,4 ), j=581,590 ) / &
         -3,          0,          0,          1, &
         -3,          0,          0,          0, &
          0,        -10,         -6,         -1, &
         31,         -6,          0,        -13, &
         -3,        -32,        -14,          1, &
         -7,          0,          0,          3, &
          0,         -8,         -4,          0, &
          3,         -4,          0,          0, &
          0,          4,          0,          0, &
          0,          3,          1,          0 /
data ( ( icpl(i,j), i=1,4 ), j=591,600 ) / &
         19,        -23,        -10,          2, &
          0,          0,          0,        -10, &
          0,          3,          2,          0, &
          0,          9,          5,         -1, &
         28,          0,          0,          0, &
          0,         -7,         -4,          0, &
          8,         -4,          0,         -4, &
          0,          0,         -2,          0, &
          0,          3,          0,          0, &
         -3,          0,          0,          1 /
data ( ( icpl(i,j), i=1,4 ), j=601,610 ) / &
         -9,          0,          1,          4, &
          3,         12,          5,         -1, &
         17,         -3,         -1,          0, &
          0,          7,          4,          0, &
         19,          0,          0,          0, &
          0,         -5,         -3,          0, &
         14,         -3,          0,         -1, &
          0,          0,         -1,          0, &
          0,          0,          0,         -5, &
          0,          5,          3,          0 /
data ( ( icpl(i,j), i=1,4 ), j=611,620 ) / &
         13,          0,          0,          0, &
          0,         -3,         -2,          0, &
          2,          9,          4,          3, &
          0,          0,          0,         -4, &
          8,          0,          0,          0, &
          0,          4,          2,          0, &
          6,          0,          0,         -3, &
          6,          0,          0,          0, &
          0,          3,          1,          0, &
          5,          0,          0,         -2 /
data ( ( icpl(i,j), i=1,4 ), j=621,630 ) / &
          3,          0,          0,         -1, &
         -3,          0,          0,          0, &
          6,          0,          0,          0, &
          7,          0,          0,          0, &
         -4,          0,          0,          0, &
          4,          0,          0,          0, &
          6,          0,          0,          0, &
          0,         -4,          0,          0, &
          0,         -4,          0,          0, &
          5,          0,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=631,640 ) / &
         -3,          0,          0,          0, &
          4,          0,          0,          0, &
         -5,          0,          0,          0, &
          4,          0,          0,          0, &
          0,          3,          0,          0, &
         13,          0,          0,          0, &
         21,         11,          0,          0, &
          0,         -5,          0,          0, &
          0,         -5,         -2,          0, &
          0,          5,          3,          0 /
data ( ( icpl(i,j), i=1,4 ), j=641,650 ) / &
          0,         -5,          0,          0, &
         -3,          0,          0,          2, &
         20,         10,          0,          0, &
        -34,          0,          0,          0, &
        -19,          0,          0,          0, &
          3,          0,          0,         -2, &
         -3,          0,          0,          1, &
         -6,          0,          0,          3, &
         -4,          0,          0,          0, &
          3,          0,          0,          0 /
data ( ( icpl(i,j), i=1,4 ), j=651,660 ) / &
          3,          0,          0,          0, &
          4,          0,          0,          0, &
          3,          0,          0,         -1, &
          6,          0,          0,         -3, &
         -8,          0,          0,          3, &
          0,          3,          1,          0, &
         -3,          0,          0,          0, &
          0,         -3,         -2,          0, &
        126,        -63,        -27,        -55, &
         -5,          0,          1,          2 /
data ( ( icpl(i,j), i=1,4 ), j=661,670 ) / &
         -3,         28,         15,          2, &
          5,          0,          1,         -2, &
          0,          9,          4,          1, &
          0,          9,          4,         -1, &
       -126,        -63,        -27,         55, &
          3,          0,          0,         -1, &
         21,        -11,         -6,        -11, &
          0,         -4,          0,          0, &
        -21,        -11,         -6,         11, &
         -3,          0,          0,          1 /
data ( ( icpl(i,j), i=1,4 ), j=671,680 ) / &
          0,          3,          1,          0, &
          8,          0,          0,         -4, &
         -6,          0,          0,          3, &
         -3,          0,          0,          1, &
          3,          0,          0,         -1, &
         -3,          0,          0,          1, &
         -5,          0,          0,          2, &
         24,        -12,         -5,        -11, &
          0,          3,          1,          0, &
          0,          3,          1,          0 /
data ( ( icpl(i,j), i=1,4 ), j=681,687 ) / &
          0,          3,          2,          0, &
        -24,        -12,         -5,         10, &
          4,          0,         -1,         -2, &
         13,          0,          0,         -6, &
          7,          0,          0,         -3, &
          3,          0,          0,         -1, &
          3,          0,          0,         -1 /

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and given date (JC).
t = ( ( date1-dj0 ) + date2 ) / djc

!  -------------------
!  LUNI-SOLAR NUTATION
!  -------------------

!
!  Fundamental (Delaunay) arguments from Simon et al. (1994)
!

!  Mean anomaly of the Moon.
el  = mod (         485868.249036d0 + &
            t*( 1717915923.2178d0 + &
            t*(         31.8792d0 + &
            t*(          0.051635d0 + &
            t*(        - 0.00024470d0 )))), turnas ) * das2r

!  Mean anomaly of the Sun.
elp = mod (        1287104.79305d0 + &
            t*(  129596581.0481d0 + &
            t*(        - 0.5532d0 + &
            t*(          0.000136d0 + &
            t*(        - 0.00001149d0 )))), turnas ) * das2r

!  Mean argument of the latitude of the Moon.
f   = mod (         335779.526232d0 + &
            t*( 1739527262.8478d0 + &
            t*(       - 12.7512d0 + &
            t*(       -  0.001037d0 + &
            t*(          0.00000417d0 )))), turnas ) * das2r

!  Mean elongation of the Moon from the Sun.
d   = mod (        1072260.70369d0 + &
            t*( 1602961601.2090d0 + &
            t*(        - 6.3706d0 + &
            t*(          0.006593d0 + &
            t*(        - 0.00003169d0 )))), turnas ) * das2r

!  Mean longitude of the ascending node of the Moon.
om  = mod (         450160.398036d0 + &
            t*(  - 6962890.5431d0 + &
            t*(          7.4722d0 + &
            t*(          0.007702d0 + &
            t*(        - 0.00005939d0 )))), turnas ) * das2r

!  Initialize the nutation values.
dp = 0d0
de = 0d0

!  Summation of luni-solar nutation series (in reverse order).
do i = nls, 1, -1

!     Argument and functions.
   arg = mod ( dble ( nals(1,i) ) * el  + &
               dble ( nals(2,i) ) * elp + &
               dble ( nals(3,i) ) * f   + &
               dble ( nals(4,i) ) * d   + &
               dble ( nals(5,i) ) * om, d2pi )
   sarg = sin(arg)
   carg = cos(arg)

!     Term.
   dp = dp + ( cls(1,i) + cls(2,i) * t ) * sarg &
           +   cls(3,i)                  * carg
   de = de + ( cls(4,i) + cls(5,i) * t ) * carg &
           +   cls(6,i)                  * sarg

end do

!  Convert from 0.1 microarcsec units to radians.
dpsils = dp * u2r
depsls = de * u2r

!  ------------------
!  PLANETARY NUTATION
!  ------------------

!  Mean anomaly of the Moon.
al   = mod ( 2.35555598d0 + 8328.6914269554d0 * t, d2pi )

!  Mean anomaly of the Sun.
alsu = mod ( 6.24006013d0 + 628.301955d0 * t, d2pi )

!  Mean argument of the latitude of the Moon.
af   = mod ( 1.627905234d0 + 8433.466158131d0 * t, d2pi )

!  Mean elongation of the Moon from the Sun.
ad   = mod ( 5.198466741d0 + 7771.3771468121d0 * t, d2pi )

!  Mean longitude of the ascending node of the Moon.
aom  = mod ( 2.18243920d0 - 33.757045d0 * t, d2pi )

!  General accumulated precession in longitude.
apa  = ( 0.02438175d0 + 0.00000538691d0 * t ) * t

!  Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
alme = mod ( 4.402608842d0 + 2608.7903141574d0 * t, d2pi )
alve = mod ( 3.176146697d0 + 1021.3285546211d0 * t, d2pi )
alea = mod ( 1.753470314d0 +  628.3075849991d0 * t, d2pi )
alma = mod ( 6.203480913d0 +  334.0612426700d0 * t, d2pi )
alju = mod ( 0.599546497d0 +   52.9690962641d0 * t, d2pi )
alsa = mod ( 0.874016757d0 +   21.3299104960d0 * t, d2pi )
alur = mod ( 5.481293871d0 +    7.4781598567d0 * t, d2pi )
alne = mod ( 5.321159000d0 +    3.8127774000d0 * t, d2pi )

!  Initialize the nutation values.
dp = 0d0
de = 0d0

!  Summation of planetary nutation series (in reverse order).
do i = npl, 1, -1

!     Argument and functions.
   arg = mod ( dble ( napl( 1,i) ) * al   + &
               dble ( napl( 2,i) ) * alsu + &
               dble ( napl( 3,i) ) * af   + &
               dble ( napl( 4,i) ) * ad   + &
               dble ( napl( 5,i) ) * aom  + &
               dble ( napl( 6,i) ) * alme + &
               dble ( napl( 7,i) ) * alve + &
               dble ( napl( 8,i) ) * alea + &
               dble ( napl( 9,i) ) * alma + &
               dble ( napl(10,i) ) * alju + &
               dble ( napl(11,i) ) * alsa + &
               dble ( napl(12,i) ) * alur + &
               dble ( napl(13,i) ) * alne + &
               dble ( napl(14,i) ) * apa, d2pi )
   sarg = sin(arg)
   carg = cos(arg)

!     Term.
   dp = dp + dble( icpl(1,i)) * sarg + dble( icpl(2,i)) * carg
   de = de + dble( icpl(3,i)) * sarg + dble( icpl(4,i)) * carg

end do

!  Convert from 0.1 microarcsec units to radians.
dpsipl = dp * u2r
depspl = de * u2r

!  -----
!  TOTAL
!  -----

!  Add planetary and luni-solar components.
dpsi = dpsipl + dpsils
deps = depspl + depsls

end subroutine nu2000a
!***********************************************************************

!***********************************************************************
!>
!  Nutation, IAU 2000A model (MHB_2000 without FCN) MODIFIED.  Series
!     truncated for speed of execution, and using Simon et al. (1994)
!     fundamental arguments throughout.  Accuracy, compared to
!     IAU 2000 A series, is 0.1 mas in delta psi and 0.04 mas in
!     delta epsilon and delta psi sin(epsilon) over 6 centuries
!     centered at year 2000 (99% of errors less than these values).
!
!  Modified form of NU2000A, by Pat Wallace, given in subroutine annex
!  to Chapter 5 of IERS Conventions (2003).
!
!  Given:
!     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
!
!  Returned:
!     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
!
!  This revision:  2002 November 25
!                  2004 March 1     (by G. Kaplan)

subroutine nu2000k ( date1, date2, dpsi, deps )

implicit none

double precision date1, date2, dpsi, deps

!  Arcseconds to radians
double precision das2r
parameter ( das2r = 4.848136811095359935899141d-6 )

!  Milliarcseconds to radians
double precision dmas2r
parameter ( dmas2r = das2r / 1d3 )

!  Arc seconds in a full circle
double precision turnas
parameter ( turnas = 1296000d0 )

!  2Pi
double precision d2pi
parameter ( d2pi = 6.283185307179586476925287d0 )

!  Units of 0.1 microarcsecond to radians
double precision u2r
parameter ( u2r = das2r/1d7 )

!  Reference epoch (J2000), JD
double precision dj0
parameter ( dj0 = 2451545d0 )

!  Days per Julian century
double precision djc
parameter ( djc = 36525d0 )

!  Miscellaneous
double precision t, el, elp, f, d, om, arg, dp, de, sarg, carg, &
                 dpsils, depsls, &
                 alme, alve, alea, alma, alju, alsa, alur, alne, &
                 apa, &
                 dpsipl, depspl
integer i, j

!  -------------------------
!  Luni-Solar nutation model
!  -------------------------

!  Number of terms in the luni-solar nutation model
integer nls
parameter ( nls = 323 )

!  Coefficients for fundamental arguments
integer nals(5,nls)

!  Longitude and obliquity coefficients
double precision cls(6,nls)

!  ---------------
!  Planetary terms
!  ---------------

!  Number of terms in the planetary nutation model
integer npl
parameter ( npl = 165 )

!  Coefficients for fundamental arguments
integer napl(14,npl)

!  Longitude and obliquity coefficients
double precision cpl(4,npl)

!  ----------------------------------------
!  Tables of argument and term coefficients
!  ----------------------------------------

!
!  Luni-Solar argument multipliers:
!               L     L'    F     D     Om
!
data ( ( nals(i,j), i=1,5 ), j=  1, 20 ) / &
          0,    0,    0,    0,    1, &
          0,    0,    2,   -2,    2, &
          0,    0,    2,    0,    2, &
          0,    0,    0,    0,    2, &
          0,    1,    0,    0,    0, &
          0,    1,    2,   -2,    2, &
          1,    0,    0,    0,    0, &
          0,    0,    2,    0,    1, &
          1,    0,    2,    0,    2, &
          0,   -1,    2,   -2,    2, &
          0,    0,    2,   -2,    1, &
         -1,    0,    2,    0,    2, &
         -1,    0,    0,    2,    0, &
          1,    0,    0,    0,    1, &
         -1,    0,    0,    0,    1, &
         -1,    0,    2,    2,    2, &
          1,    0,    2,    0,    1, &
         -2,    0,    2,    0,    1, &
          0,    0,    0,    2,    0, &
          0,    0,    2,    2,    2/
data ( ( nals(i,j), i=1,5 ), j= 21, 40 ) / &
          0,   -2,    2,   -2,    2, &
         -2,    0,    0,    2,    0, &
          2,    0,    2,    0,    2, &
          1,    0,    2,   -2,    2, &
         -1,    0,    2,    0,    1, &
          2,    0,    0,    0,    0, &
          0,    0,    2,    0,    0, &
          0,    1,    0,    0,    1, &
         -1,    0,    0,    2,    1, &
          0,    2,    2,   -2,    2, &
          0,    0,   -2,    2,    0, &
          1,    0,    0,   -2,    1, &
          0,   -1,    0,    0,    1, &
         -1,    0,    2,    2,    1, &
          0,    2,    0,    0,    0, &
          1,    0,    2,    2,    2, &
         -2,    0,    2,    0,    0, &
          0,    1,    2,    0,    2, &
          0,    0,    2,    2,    1, &
          0,   -1,    2,    0,    2/
data ( ( nals(i,j), i=1,5 ), j= 41, 60 ) / &
          0,    0,    0,    2,    1, &
          1,    0,    2,   -2,    1, &
          2,    0,    2,   -2,    2, &
         -2,    0,    0,    2,    1, &
          2,    0,    2,    0,    1, &
          0,   -1,    2,   -2,    1, &
          0,    0,    0,   -2,    1, &
         -1,   -1,    0,    2,    0, &
          2,    0,    0,   -2,    1, &
          1,    0,    0,    2,    0, &
          0,    1,    2,   -2,    1, &
          1,   -1,    0,    0,    0, &
         -2,    0,    2,    0,    2, &
          3,    0,    2,    0,    2, &
          0,   -1,    0,    2,    0, &
          1,   -1,    2,    0,    2, &
          0,    0,    0,    1,    0, &
         -1,   -1,    2,    2,    2, &
         -1,    0,    2,    0,    0, &
          0,   -1,    2,    2,    2/
data ( ( nals(i,j), i=1,5 ), j= 61, 80 ) / &
         -2,    0,    0,    0,    1, &
          1,    1,    2,    0,    2, &
          2,    0,    0,    0,    1, &
         -1,    1,    0,    1,    0, &
          1,    1,    0,    0,    0, &
          1,    0,    2,    0,    0, &
         -1,    0,    2,   -2,    1, &
          1,    0,    0,    0,    2, &
         -1,    0,    0,    1,    0, &
          0,    0,    2,    1,    2, &
         -1,    0,    2,    4,    2, &
         -1,    1,    0,    1,    1, &
          0,   -2,    2,   -2,    1, &
          1,    0,    2,    2,    1, &
         -2,    0,    2,    2,    2, &
         -1,    0,    0,    0,    2, &
          1,    1,    2,   -2,    2, &
         -2,    0,    2,    4,    2, &
         -1,    0,    4,    0,    2, &
          2,    0,    2,   -2,    1/
data ( ( nals(i,j), i=1,5 ), j= 81,100 ) / &
          2,    0,    2,    2,    2, &
          1,    0,    0,    2,    1, &
          3,    0,    0,    0,    0, &
          3,    0,    2,   -2,    2, &
          0,    0,    4,   -2,    2, &
          0,    1,    2,    0,    1, &
          0,    0,   -2,    2,    1, &
          0,    0,    2,   -2,    3, &
         -1,    0,    0,    4,    0, &
          2,    0,   -2,    0,    1, &
         -2,    0,    0,    4,    0, &
         -1,   -1,    0,    2,    1, &
         -1,    0,    0,    1,    1, &
          0,    1,    0,    0,    2, &
          0,    0,   -2,    0,    1, &
          0,   -1,    2,    0,    1, &
          0,    0,    2,   -1,    2, &
          0,    0,    2,    4,    2, &
         -2,   -1,    0,    2,    0, &
          1,    1,    0,   -2,    1/
data ( ( nals(i,j), i=1,5 ), j=101,120 ) / &
         -1,    1,    0,    2,    0, &
         -1,    1,    0,    1,    2, &
          1,   -1,    0,    0,    1, &
          1,   -1,    2,    2,    2, &
         -1,    1,    2,    2,    2, &
          3,    0,    2,    0,    1, &
          0,    1,   -2,    2,    0, &
         -1,    0,    0,   -2,    1, &
          0,    1,    2,    2,    2, &
         -1,   -1,    2,    2,    1, &
          0,   -1,    0,    0,    2, &
          1,    0,    2,   -4,    1, &
         -1,    0,   -2,    2,    0, &
          0,   -1,    2,    2,    1, &
          2,   -1,    2,    0,    2, &
          0,    0,    0,    2,    2, &
          1,   -1,    2,    0,    1, &
         -1,    1,    2,    0,    2, &
          0,    1,    0,    2,    0, &
          0,   -1,   -2,    2,    0/
data ( ( nals(i,j), i=1,5 ), j=121,140 ) / &
          0,    3,    2,   -2,    2, &
          0,    0,    0,    1,    1, &
         -1,    0,    2,    2,    0, &
          2,    1,    2,    0,    2, &
          1,    1,    0,    0,    1, &
          1,    1,    2,    0,    1, &
          2,    0,    0,    2,    0, &
          1,    0,   -2,    2,    0, &
         -1,    0,    0,    2,    2, &
          0,    1,    0,    1,    0, &
          0,    1,    0,   -2,    1, &
         -1,    0,    2,   -2,    2, &
          0,    0,    0,   -1,    1, &
         -1,    1,    0,    0,    1, &
          1,    0,    2,   -1,    2, &
          1,   -1,    0,    2,    0, &
          0,    0,    0,    4,    0, &
          1,    0,    2,    1,    2, &
          0,    0,    2,    1,    1, &
          1,    0,    0,   -2,    2/
data ( ( nals(i,j), i=1,5 ), j=141,160 ) / &
         -1,    0,    2,    4,    1, &
          1,    0,   -2,    0,    1, &
          1,    1,    2,   -2,    1, &
          0,    0,    2,    2,    0, &
         -1,    0,    2,   -1,    1, &
         -2,    0,    2,    2,    1, &
          4,    0,    2,    0,    2, &
          2,   -1,    0,    0,    0, &
          2,    1,    2,   -2,    2, &
          0,    1,    2,    1,    2, &
          1,    0,    4,   -2,    2, &
         -1,   -1,    0,    0,    1, &
          0,    1,    0,    2,    1, &
         -2,    0,    2,    4,    1, &
          2,    0,    2,    0,    0, &
          1,    0,    0,    1,    0, &
         -1,    0,    0,    4,    1, &
         -1,    0,    4,    0,    1, &
          2,    0,    2,    2,    1, &
          0,    0,    2,   -3,    2/
data ( ( nals(i,j), i=1,5 ), j=161,180 ) / &
         -1,   -2,    0,    2,    0, &
          2,    1,    0,    0,    0, &
          0,    0,    4,    0,    2, &
          0,    0,    0,    0,    3, &
          0,    3,    0,    0,    0, &
          0,    0,    2,   -4,    1, &
          0,   -1,    0,    2,    1, &
          0,    0,    0,    4,    1, &
         -1,   -1,    2,    4,    2, &
          1,    0,    2,    4,    2, &
         -2,    2,    0,    2,    0, &
         -2,   -1,    2,    0,    1, &
         -2,    0,    0,    2,    2, &
         -1,   -1,    2,    0,    2, &
          0,    0,    4,   -2,    1, &
          3,    0,    2,   -2,    1, &
         -2,   -1,    0,    2,    1, &
          1,    0,    0,   -1,    1, &
          0,   -2,    0,    2,    0, &
         -2,    0,    0,    4,    1/
data ( ( nals(i,j), i=1,5 ), j=181,200 ) / &
         -3,    0,    0,    0,    1, &
          1,    1,    2,    2,    2, &
          0,    0,    2,    4,    1, &
          3,    0,    2,    2,    2, &
         -1,    1,    2,   -2,    1, &
          2,    0,    0,   -4,    1, &
          0,    0,    0,   -2,    2, &
          2,    0,    2,   -4,    1, &
         -1,    1,    0,    2,    1, &
          0,    0,    2,   -1,    1, &
          0,   -2,    2,    2,    2, &
          2,    0,    0,    2,    1, &
          4,    0,    2,   -2,    2, &
          2,    0,    0,   -2,    2, &
          0,    2,    0,    0,    1, &
          1,    0,    0,   -4,    1, &
          0,    2,    2,   -2,    1, &
         -3,    0,    0,    4,    0, &
         -1,    1,    2,    0,    1, &
         -1,   -1,    0,    4,    0/
data ( ( nals(i,j), i=1,5 ), j=201,220 ) / &
         -1,   -2,    2,    2,    2, &
         -2,   -1,    2,    4,    2, &
          1,   -1,    2,    2,    1, &
         -2,    1,    0,    2,    0, &
         -2,    1,    2,    0,    1, &
          2,    1,    0,   -2,    1, &
         -3,    0,    2,    0,    1, &
         -2,    0,    2,   -2,    1, &
         -1,    1,    0,    2,    2, &
          0,   -1,    2,   -1,    2, &
         -1,    0,    4,   -2,    2, &
          0,   -2,    2,    0,    2, &
         -1,    0,    2,    1,    2, &
          2,    0,    0,    0,    2, &
          0,    0,    2,    0,    3, &
         -2,    0,    4,    0,    2, &
         -1,    0,   -2,    0,    1, &
         -1,    1,    2,    2,    1, &
          3,    0,    0,    0,    1, &
         -1,    0,    2,    3,    2/
data ( ( nals(i,j), i=1,5 ), j=221,240 ) / &
          2,   -1,    2,    0,    1, &
          0,    1,    2,    2,    1, &
          0,   -1,    2,    4,    2, &
          2,   -1,    2,    2,    2, &
          0,    2,   -2,    2,    0, &
         -1,   -1,    2,   -1,    1, &
          0,   -2,    0,    0,    1, &
          1,    0,    2,   -4,    2, &
          1,   -1,    0,   -2,    1, &
         -1,   -1,    2,    0,    1, &
          1,   -1,    2,   -2,    2, &
         -2,   -1,    0,    4,    0, &
         -1,    0,    0,    3,    0, &
         -2,   -1,    2,    2,    2, &
          0,    2,    2,    0,    2, &
          1,    1,    0,    2,    0, &
          2,    0,    2,   -1,    2, &
          1,    0,    2,    1,    1, &
          4,    0,    0,    0,    0, &
          2,    1,    2,    0,    1/
data ( ( nals(i,j), i=1,5 ), j=241,260 ) / &
          3,   -1,    2,    0,    2, &
         -2,    2,    0,    2,    1, &
          1,    0,    2,   -3,    1, &
          1,    1,    2,   -4,    1, &
         -1,   -1,    2,   -2,    1, &
          0,   -1,    0,   -1,    1, &
          0,   -1,    0,   -2,    1, &
         -2,    0,    0,    0,    2, &
         -2,    0,   -2,    2,    0, &
         -1,    0,   -2,    4,    0, &
          1,   -2,    0,    0,    0, &
          0,    1,    0,    1,    1, &
         -1,    2,    0,    2,    0, &
          1,   -1,    2,   -2,    1, &
          1,    2,    2,   -2,    2, &
          2,   -1,    2,   -2,    2, &
          1,    0,    2,   -1,    1, &
          2,    1,    2,   -2,    1, &
         -2,    0,    0,   -2,    1, &
          1,   -2,    2,    0,    2/
data ( ( nals(i,j), i=1,5 ), j=261,280 ) / &
          0,    1,    2,    1,    1, &
          1,    0,    4,   -2,    1, &
         -2,    0,    4,    2,    2, &
          1,    1,    2,    1,    2, &
          1,    0,    0,    4,    0, &
          1,    0,    2,    2,    0, &
          2,    0,    2,    1,    2, &
          3,    1,    2,    0,    2, &
          4,    0,    2,    0,    1, &
         -2,   -1,    2,    0,    0, &
          0,    1,   -2,    2,    1, &
          1,    0,   -2,    1,    0, &
          2,   -1,    0,   -2,    1, &
         -1,    0,    2,   -1,    2, &
          1,    0,    2,   -3,    2, &
          0,    1,    2,   -2,    3, &
         -1,    0,   -2,    2,    1, &
          0,    0,    2,   -4,    2, &
          2,    0,    2,   -4,    2, &
          0,    0,    4,   -4,    4/
data ( ( nals(i,j), i=1,5 ), j=281,300 ) / &
          0,    0,    4,   -4,    2, &
         -2,    0,    0,    3,    0, &
          1,    0,   -2,    2,    1, &
         -3,    0,    2,    2,    2, &
         -2,    0,    2,    2,    0, &
          2,   -1,    0,    0,    1, &
          1,    1,    0,    1,    0, &
          0,    1,    4,   -2,    2, &
         -1,    1,    0,   -2,    1, &
          0,    0,    0,   -4,    1, &
          1,   -1,    0,    2,    1, &
          1,    1,    0,    2,    1, &
         -1,    2,    2,    2,    2, &
          3,    1,    2,   -2,    2, &
          0,   -1,    0,    4,    0, &
          2,   -1,    0,    2,    0, &
          0,    0,    4,    0,    1, &
          2,    0,    4,   -2,    2, &
         -1,   -1,    2,    4,    1, &
          1,    0,    0,    4,    1/
data ( ( nals(i,j), i=1,5 ), j=301,320 ) / &
          1,   -2,    2,    2,    2, &
          0,    0,    2,    3,    2, &
         -1,    1,    2,    4,    2, &
          3,    0,    0,    2,    0, &
         -1,    0,    4,    2,    2, &
         -2,    0,    2,    6,    2, &
         -1,    0,    2,    6,    2, &
          1,    1,   -2,    1,    0, &
         -1,    0,    0,    1,    2, &
         -1,   -1,    0,    1,    0, &
         -2,    0,    0,    1,    0, &
          0,    0,   -2,    1,    0, &
          1,   -1,   -2,    2,    0, &
          1,    2,    0,    0,    0, &
          3,    0,    2,    0,    0, &
          0,   -1,    1,   -1,    1, &
         -1,    0,    1,    0,    3, &
         -1,    0,    1,    0,    2, &
         -1,    0,    1,    0,    1, &
         -1,    0,    1,    0,    0/
data ( ( nals(i,j), i=1,5 ), j=321,323 ) / &
          0,    0,    1,    0,    2, &
          0,    0,    1,    0,    1, &
          0,    0,    1,    0,    0/

!
!  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
!  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
!
data ( ( cls(i,j), i=1,6 ), j=  1, 20 ) / &
-172064161.d0,-174666.d0, 33386.d0, 92052331.d0, 9086.d0,15377.d0, &    !
 -13170906.d0,  -1675.d0,-13696.d0,  5730336.d0,-3015.d0,-4587.d0, &    !
  -2276413.d0,   -234.d0,  2796.d0,   978459.d0, -485.d0, 1374.d0, &    !
   2074554.d0,    207.d0,  -698.d0,  -897492.d0,  470.d0, -291.d0, &    !
   1475877.d0,  -3633.d0, 11817.d0,    73871.d0, -184.d0,-1924.d0, &    !
   -516821.d0,   1226.d0,  -524.d0,   224386.d0, -677.d0, -174.d0, &    !
    711159.d0,     73.d0,  -872.d0,    -6750.d0,    0.d0,  358.d0, &    !
   -387298.d0,   -367.d0,   380.d0,   200728.d0,   18.d0,  318.d0, &    !
   -301461.d0,    -36.d0,   816.d0,   129025.d0,  -63.d0,  367.d0, &    !
    215829.d0,   -494.d0,   111.d0,   -95929.d0,  299.d0,  132.d0, &    !
    128227.d0,    137.d0,   181.d0,   -68982.d0,   -9.d0,   39.d0, &    !
    123457.d0,     11.d0,    19.d0,   -53311.d0,   32.d0,   -4.d0, &    !
    156994.d0,     10.d0,  -168.d0,    -1235.d0,    0.d0,   82.d0, &    !
     63110.d0,     63.d0,    27.d0,   -33228.d0,    0.d0,   -9.d0, &    !
    -57976.d0,    -63.d0,  -189.d0,    31429.d0,    0.d0,  -75.d0, &    !
    -59641.d0,    -11.d0,   149.d0,    25543.d0,  -11.d0,   66.d0, &    !
    -51613.d0,    -42.d0,   129.d0,    26366.d0,    0.d0,   78.d0, &    !
     45893.d0,     50.d0,    31.d0,   -24236.d0,  -10.d0,   20.d0, &    !
     63384.d0,     11.d0,  -150.d0,    -1220.d0,    0.d0,   29.d0, &    !
    -38571.d0,     -1.d0,   158.d0,    16452.d0,  -11.d0,   68.d0/      !
data ( ( cls(i,j), i=1,6 ), j= 21, 40 ) / &
     32481.d0,      0.d0,     0.d0,   -13870.d0,    0.d0,    0.d0, &    !
    -47722.d0,      0.d0,   -18.d0,      477.d0,    0.d0,  -25.d0, &    !
    -31046.d0,     -1.d0,   131.d0,    13238.d0,  -11.d0,   59.d0, &    !
     28593.d0,      0.d0,    -1.d0,   -12338.d0,   10.d0,   -3.d0, &    !
     20441.d0,     21.d0,    10.d0,   -10758.d0,    0.d0,   -3.d0, &    !
     29243.d0,      0.d0,   -74.d0,     -609.d0,    0.d0,   13.d0, &    !
     25887.d0,      0.d0,   -66.d0,     -550.d0,    0.d0,   11.d0, &    !
    -14053.d0,    -25.d0,    79.d0,     8551.d0,   -2.d0,  -45.d0, &    !
     15164.d0,     10.d0,    11.d0,    -8001.d0,    0.d0,   -1.d0, &    !
    -15794.d0,     72.d0,   -16.d0,     6850.d0,  -42.d0,   -5.d0, &    !
     21783.d0,      0.d0,    13.d0,     -167.d0,    0.d0,   13.d0, &    !
    -12873.d0,    -10.d0,   -37.d0,     6953.d0,    0.d0,  -14.d0, &    !
    -12654.d0,     11.d0,    63.d0,     6415.d0,    0.d0,   26.d0, &    !
    -10204.d0,      0.d0,    25.d0,     5222.d0,    0.d0,   15.d0, &    !
     16707.d0,    -85.d0,   -10.d0,      168.d0,   -1.d0,   10.d0, &    !
     -7691.d0,      0.d0,    44.d0,     3268.d0,    0.d0,   19.d0, &    !
    -11024.d0,      0.d0,   -14.d0,      104.d0,    0.d0,    2.d0, &    !
      7566.d0,    -21.d0,   -11.d0,    -3250.d0,    0.d0,   -5.d0, &    !
     -6637.d0,    -11.d0,    25.d0,     3353.d0,    0.d0,   14.d0, &    !
     -7141.d0,     21.d0,     8.d0,     3070.d0,    0.d0,    4.d0/      !
data ( ( cls(i,j), i=1,6 ), j= 41, 60 ) / &
     -6302.d0,    -11.d0,     2.d0,     3272.d0,    0.d0,    4.d0, &    !
      5800.d0,     10.d0,     2.d0,    -3045.d0,    0.d0,   -1.d0, &    !
      6443.d0,      0.d0,    -7.d0,    -2768.d0,    0.d0,   -4.d0, &    !
     -5774.d0,    -11.d0,   -15.d0,     3041.d0,    0.d0,   -5.d0, &    !
     -5350.d0,      0.d0,    21.d0,     2695.d0,    0.d0,   12.d0, &    !
     -4752.d0,    -11.d0,    -3.d0,     2719.d0,    0.d0,   -3.d0, &    !
     -4940.d0,    -11.d0,   -21.d0,     2720.d0,    0.d0,   -9.d0, &    !
      7350.d0,      0.d0,    -8.d0,      -51.d0,    0.d0,    4.d0, &    !
      4065.d0,      0.d0,     6.d0,    -2206.d0,    0.d0,    1.d0, &    !
      6579.d0,      0.d0,   -24.d0,     -199.d0,    0.d0,    2.d0, &    !
      3579.d0,      0.d0,     5.d0,    -1900.d0,    0.d0,    1.d0, &    !
      4725.d0,      0.d0,    -6.d0,      -41.d0,    0.d0,    3.d0, &    !
     -3075.d0,      0.d0,    -2.d0,     1313.d0,    0.d0,   -1.d0, &    !
     -2904.d0,      0.d0,    15.d0,     1233.d0,    0.d0,    7.d0, &    !
      4348.d0,      0.d0,   -10.d0,      -81.d0,    0.d0,    2.d0, &    !
     -2878.d0,      0.d0,     8.d0,     1232.d0,    0.d0,    4.d0, &    !
     -4230.d0,      0.d0,     5.d0,      -20.d0,    0.d0,   -2.d0, &    !
     -2819.d0,      0.d0,     7.d0,     1207.d0,    0.d0,    3.d0, &    !
     -4056.d0,      0.d0,     5.d0,       40.d0,    0.d0,   -2.d0, &    !
     -2647.d0,      0.d0,    11.d0,     1129.d0,    0.d0,    5.d0/      !
data ( ( cls(i,j), i=1,6 ), j= 61, 80 ) / &
     -2294.d0,      0.d0,   -10.d0,     1266.d0,    0.d0,   -4.d0, &    !
      2481.d0,      0.d0,    -7.d0,    -1062.d0,    0.d0,   -3.d0, &    !
      2179.d0,      0.d0,    -2.d0,    -1129.d0,    0.d0,   -2.d0, &    !
      3276.d0,      0.d0,     1.d0,       -9.d0,    0.d0,    0.d0, &    !
     -3389.d0,      0.d0,     5.d0,       35.d0,    0.d0,   -2.d0, &    !
      3339.d0,      0.d0,   -13.d0,     -107.d0,    0.d0,    1.d0, &    !
     -1987.d0,      0.d0,    -6.d0,     1073.d0,    0.d0,   -2.d0, &    !
     -1981.d0,      0.d0,     0.d0,      854.d0,    0.d0,    0.d0, &    !
      4026.d0,      0.d0,  -353.d0,     -553.d0,    0.d0, -139.d0, &    !
      1660.d0,      0.d0,    -5.d0,     -710.d0,    0.d0,   -2.d0, &    !
     -1521.d0,      0.d0,     9.d0,      647.d0,    0.d0,    4.d0, &    !
      1314.d0,      0.d0,     0.d0,     -700.d0,    0.d0,    0.d0, &    !
     -1283.d0,      0.d0,     0.d0,      672.d0,    0.d0,    0.d0, &    !
     -1331.d0,      0.d0,     8.d0,      663.d0,    0.d0,    4.d0, &    !
      1383.d0,      0.d0,    -2.d0,     -594.d0,    0.d0,   -2.d0, &    !
      1405.d0,      0.d0,     4.d0,     -610.d0,    0.d0,    2.d0, &    !
      1290.d0,      0.d0,     0.d0,     -556.d0,    0.d0,    0.d0, &    !
     -1214.d0,      0.d0,     5.d0,      518.d0,    0.d0,    2.d0, &    !
      1146.d0,      0.d0,    -3.d0,     -490.d0,    0.d0,   -1.d0, &    !
      1019.d0,      0.d0,    -1.d0,     -527.d0,    0.d0,   -1.d0/      !
data ( ( cls(i,j), i=1,6 ), j= 81,100 ) / &
     -1100.d0,      0.d0,     9.d0,      465.d0,    0.d0,    4.d0, &    !
      -970.d0,      0.d0,     2.d0,      496.d0,    0.d0,    1.d0, &    !
      1575.d0,      0.d0,    -6.d0,      -50.d0,    0.d0,    0.d0, &    !
       934.d0,      0.d0,    -3.d0,     -399.d0,    0.d0,   -1.d0, &    !
       922.d0,      0.d0,    -1.d0,     -395.d0,    0.d0,   -1.d0, &    !
       815.d0,      0.d0,    -1.d0,     -422.d0,    0.d0,   -1.d0, &    !
       834.d0,      0.d0,     2.d0,     -440.d0,    0.d0,    1.d0, &    !
      1248.d0,      0.d0,     0.d0,     -170.d0,    0.d0,    1.d0, &    !
      1338.d0,      0.d0,    -5.d0,      -39.d0,    0.d0,    0.d0, &    !
       716.d0,      0.d0,    -2.d0,     -389.d0,    0.d0,   -1.d0, &    !
      1282.d0,      0.d0,    -3.d0,      -23.d0,    0.d0,    1.d0, &    !
       742.d0,      0.d0,     1.d0,     -391.d0,    0.d0,    0.d0, &    !
      1020.d0,      0.d0,   -25.d0,     -495.d0,    0.d0,  -10.d0, &    !
       715.d0,      0.d0,    -4.d0,     -326.d0,    0.d0,    2.d0, &    !
      -666.d0,      0.d0,    -3.d0,      369.d0,    0.d0,   -1.d0, &    !
      -667.d0,      0.d0,     1.d0,      346.d0,    0.d0,    1.d0, &    !
      -704.d0,      0.d0,     0.d0,      304.d0,    0.d0,    0.d0, &    !
      -694.d0,      0.d0,     5.d0,      294.d0,    0.d0,    2.d0, &    !
     -1014.d0,      0.d0,    -1.d0,        4.d0,    0.d0,   -1.d0, &    !
      -585.d0,      0.d0,    -2.d0,      316.d0,    0.d0,   -1.d0/      !
data ( ( cls(i,j), i=1,6 ), j=101,120 ) / &
      -949.d0,      0.d0,     1.d0,        8.d0,    0.d0,   -1.d0, &    !
      -595.d0,      0.d0,     0.d0,      258.d0,    0.d0,    0.d0, &    !
       528.d0,      0.d0,     0.d0,     -279.d0,    0.d0,    0.d0, &    !
      -590.d0,      0.d0,     4.d0,      252.d0,    0.d0,    2.d0, &    !
       570.d0,      0.d0,    -2.d0,     -244.d0,    0.d0,   -1.d0, &    !
      -502.d0,      0.d0,     3.d0,      250.d0,    0.d0,    2.d0, &    !
      -875.d0,      0.d0,     1.d0,       29.d0,    0.d0,    0.d0, &    !
      -492.d0,      0.d0,    -3.d0,      275.d0,    0.d0,   -1.d0, &    !
       535.d0,      0.d0,    -2.d0,     -228.d0,    0.d0,   -1.d0, &    !
      -467.d0,      0.d0,     1.d0,      240.d0,    0.d0,    1.d0, &    !
       591.d0,      0.d0,     0.d0,     -253.d0,    0.d0,    0.d0, &    !
      -453.d0,      0.d0,    -1.d0,      244.d0,    0.d0,   -1.d0, &    !
       766.d0,      0.d0,     1.d0,        9.d0,    0.d0,    0.d0, &    !
      -446.d0,      0.d0,     2.d0,      225.d0,    0.d0,    1.d0, &    !
      -488.d0,      0.d0,     2.d0,      207.d0,    0.d0,    1.d0, &    !
      -468.d0,      0.d0,     0.d0,      201.d0,    0.d0,    0.d0, &    !
      -421.d0,      0.d0,     1.d0,      216.d0,    0.d0,    1.d0, &    !
       463.d0,      0.d0,     0.d0,     -200.d0,    0.d0,    0.d0, &    !
      -673.d0,      0.d0,     2.d0,       14.d0,    0.d0,    0.d0, &    !
       658.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=121,140 ) / &
      -438.d0,      0.d0,     0.d0,      188.d0,    0.d0,    0.d0, &    !
      -390.d0,      0.d0,     0.d0,      205.d0,    0.d0,    0.d0, &    !
       639.d0,    -11.d0,    -2.d0,      -19.d0,    0.d0,    0.d0, &    !
       412.d0,      0.d0,    -2.d0,     -176.d0,    0.d0,   -1.d0, &    !
      -361.d0,      0.d0,     0.d0,      189.d0,    0.d0,    0.d0, &    !
       360.d0,      0.d0,    -1.d0,     -185.d0,    0.d0,   -1.d0, &    !
       588.d0,      0.d0,    -3.d0,      -24.d0,    0.d0,    0.d0, &    !
      -578.d0,      0.d0,     1.d0,        5.d0,    0.d0,    0.d0, &    !
      -396.d0,      0.d0,     0.d0,      171.d0,    0.d0,    0.d0, &    !
       565.d0,      0.d0,    -1.d0,       -6.d0,    0.d0,    0.d0, &    !
      -335.d0,      0.d0,    -1.d0,      184.d0,    0.d0,   -1.d0, &    !
       357.d0,      0.d0,     1.d0,     -154.d0,    0.d0,    0.d0, &    !
       321.d0,      0.d0,     1.d0,     -174.d0,    0.d0,    0.d0, &    !
      -301.d0,      0.d0,    -1.d0,      162.d0,    0.d0,    0.d0, &    !
      -334.d0,      0.d0,     0.d0,      144.d0,    0.d0,    0.d0, &    !
       493.d0,      0.d0,    -2.d0,      -15.d0,    0.d0,    0.d0, &    !
       494.d0,      0.d0,    -2.d0,      -19.d0,    0.d0,    0.d0, &    !
       337.d0,      0.d0,    -1.d0,     -143.d0,    0.d0,   -1.d0, &    !
       280.d0,      0.d0,    -1.d0,     -144.d0,    0.d0,    0.d0, &    !
       309.d0,      0.d0,     1.d0,     -134.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=141,160 ) / &
      -263.d0,      0.d0,     2.d0,      131.d0,    0.d0,    1.d0, &    !
       253.d0,      0.d0,     1.d0,     -138.d0,    0.d0,    0.d0, &    !
       245.d0,      0.d0,     0.d0,     -128.d0,    0.d0,    0.d0, &    !
       416.d0,      0.d0,    -2.d0,      -17.d0,    0.d0,    0.d0, &    !
      -229.d0,      0.d0,     0.d0,      128.d0,    0.d0,    0.d0, &    !
       231.d0,      0.d0,     0.d0,     -120.d0,    0.d0,    0.d0, &    !
      -259.d0,      0.d0,     2.d0,      109.d0,    0.d0,    1.d0, &    !
       375.d0,      0.d0,    -1.d0,       -8.d0,    0.d0,    0.d0, &    !
       252.d0,      0.d0,     0.d0,     -108.d0,    0.d0,    0.d0, &    !
      -245.d0,      0.d0,     1.d0,      104.d0,    0.d0,    0.d0, &    !
       243.d0,      0.d0,    -1.d0,     -104.d0,    0.d0,    0.d0, &    !
       208.d0,      0.d0,     1.d0,     -112.d0,    0.d0,    0.d0, &    !
       199.d0,      0.d0,     0.d0,     -102.d0,    0.d0,    0.d0, &    !
      -208.d0,      0.d0,     1.d0,      105.d0,    0.d0,    0.d0, &    !
       335.d0,      0.d0,    -2.d0,      -14.d0,    0.d0,    0.d0, &    !
      -325.d0,      0.d0,     1.d0,        7.d0,    0.d0,    0.d0, &    !
      -187.d0,      0.d0,     0.d0,       96.d0,    0.d0,    0.d0, &    !
       197.d0,      0.d0,    -1.d0,     -100.d0,    0.d0,    0.d0, &    !
      -192.d0,      0.d0,     2.d0,       94.d0,    0.d0,    1.d0, &    !
      -188.d0,      0.d0,     0.d0,       83.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=161,180 ) / &
       276.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
      -286.d0,      0.d0,     1.d0,        6.d0,    0.d0,    0.d0, &    !
       186.d0,      0.d0,    -1.d0,      -79.d0,    0.d0,    0.d0, &    !
      -219.d0,      0.d0,     0.d0,       43.d0,    0.d0,    0.d0, &    !
       276.d0,      0.d0,     0.d0,        2.d0,    0.d0,    0.d0, &    !
      -153.d0,      0.d0,    -1.d0,       84.d0,    0.d0,    0.d0, &    !
      -156.d0,      0.d0,     0.d0,       81.d0,    0.d0,    0.d0, &    !
      -154.d0,      0.d0,     1.d0,       78.d0,    0.d0,    0.d0, &    !
      -174.d0,      0.d0,     1.d0,       75.d0,    0.d0,    0.d0, &    !
      -163.d0,      0.d0,     2.d0,       69.d0,    0.d0,    1.d0, &    !
      -228.d0,      0.d0,     0.d0,        1.d0,    0.d0,    0.d0, &    !
        91.d0,      0.d0,    -4.d0,      -54.d0,    0.d0,   -2.d0, &    !
       175.d0,      0.d0,     0.d0,      -75.d0,    0.d0,    0.d0, &    !
      -159.d0,      0.d0,     0.d0,       69.d0,    0.d0,    0.d0, &    !
       141.d0,      0.d0,     0.d0,      -72.d0,    0.d0,    0.d0, &    !
       147.d0,      0.d0,     0.d0,      -75.d0,    0.d0,    0.d0, &    !
      -132.d0,      0.d0,     0.d0,       69.d0,    0.d0,    0.d0, &    !
       159.d0,      0.d0,   -28.d0,      -54.d0,    0.d0,   11.d0, &    !
       213.d0,      0.d0,     0.d0,       -4.d0,    0.d0,    0.d0, &    !
       123.d0,      0.d0,     0.d0,      -64.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=181,200 ) / &
      -118.d0,      0.d0,    -1.d0,       66.d0,    0.d0,    0.d0, &    !
       144.d0,      0.d0,    -1.d0,      -61.d0,    0.d0,    0.d0, &    !
      -121.d0,      0.d0,     1.d0,       60.d0,    0.d0,    0.d0, &    !
      -134.d0,      0.d0,     1.d0,       56.d0,    0.d0,    1.d0, &    !
      -105.d0,      0.d0,     0.d0,       57.d0,    0.d0,    0.d0, &    !
      -102.d0,      0.d0,     0.d0,       56.d0,    0.d0,    0.d0, &    !
       120.d0,      0.d0,     0.d0,      -52.d0,    0.d0,    0.d0, &    !
       101.d0,      0.d0,     0.d0,      -54.d0,    0.d0,    0.d0, &    !
      -113.d0,      0.d0,     0.d0,       59.d0,    0.d0,    0.d0, &    !
      -106.d0,      0.d0,     0.d0,       61.d0,    0.d0,    0.d0, &    !
      -129.d0,      0.d0,     1.d0,       55.d0,    0.d0,    0.d0, &    !
      -114.d0,      0.d0,     0.d0,       57.d0,    0.d0,    0.d0, &    !
       113.d0,      0.d0,    -1.d0,      -49.d0,    0.d0,    0.d0, &    !
      -102.d0,      0.d0,     0.d0,       44.d0,    0.d0,    0.d0, &    !
       -94.d0,      0.d0,     0.d0,       51.d0,    0.d0,    0.d0, &    !
      -100.d0,      0.d0,    -1.d0,       56.d0,    0.d0,    0.d0, &    !
        87.d0,      0.d0,     0.d0,      -47.d0,    0.d0,    0.d0, &    !
       161.d0,      0.d0,     0.d0,       -1.d0,    0.d0,    0.d0, &    !
        96.d0,      0.d0,     0.d0,      -50.d0,    0.d0,    0.d0, &    !
       151.d0,      0.d0,    -1.d0,       -5.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=201,220 ) / &
      -104.d0,      0.d0,     0.d0,       44.d0,    0.d0,    0.d0, &    !
      -110.d0,      0.d0,     0.d0,       48.d0,    0.d0,    0.d0, &    !
      -100.d0,      0.d0,     1.d0,       50.d0,    0.d0,    0.d0, &    !
        92.d0,      0.d0,    -5.d0,       12.d0,    0.d0,   -2.d0, &    !
        82.d0,      0.d0,     0.d0,      -45.d0,    0.d0,    0.d0, &    !
        82.d0,      0.d0,     0.d0,      -45.d0,    0.d0,    0.d0, &    !
       -78.d0,      0.d0,     0.d0,       41.d0,    0.d0,    0.d0, &    !
       -77.d0,      0.d0,     0.d0,       43.d0,    0.d0,    0.d0, &    !
         2.d0,      0.d0,     0.d0,       54.d0,    0.d0,    0.d0, &    !
        94.d0,      0.d0,     0.d0,      -40.d0,    0.d0,    0.d0, &    !
       -93.d0,      0.d0,     0.d0,       40.d0,    0.d0,    0.d0, &    !
       -83.d0,      0.d0,    10.d0,       40.d0,    0.d0,   -2.d0, &    !
        83.d0,      0.d0,     0.d0,      -36.d0,    0.d0,    0.d0, &    !
       -91.d0,      0.d0,     0.d0,       39.d0,    0.d0,    0.d0, &    !
       128.d0,      0.d0,     0.d0,       -1.d0,    0.d0,    0.d0, &    !
       -79.d0,      0.d0,     0.d0,       34.d0,    0.d0,    0.d0, &    !
       -83.d0,      0.d0,     0.d0,       47.d0,    0.d0,    0.d0, &    !
        84.d0,      0.d0,     0.d0,      -44.d0,    0.d0,    0.d0, &    !
        83.d0,      0.d0,     0.d0,      -43.d0,    0.d0,    0.d0, &    !
        91.d0,      0.d0,     0.d0,      -39.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=221,240 ) / &
       -77.d0,      0.d0,     0.d0,       39.d0,    0.d0,    0.d0, &    !
        84.d0,      0.d0,     0.d0,      -43.d0,    0.d0,    0.d0, &    !
       -92.d0,      0.d0,     1.d0,       39.d0,    0.d0,    0.d0, &    !
       -92.d0,      0.d0,     1.d0,       39.d0,    0.d0,    0.d0, &    !
       -94.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
        68.d0,      0.d0,     0.d0,      -36.d0,    0.d0,    0.d0, &    !
       -61.d0,      0.d0,     0.d0,       32.d0,    0.d0,    0.d0, &    !
        71.d0,      0.d0,     0.d0,      -31.d0,    0.d0,    0.d0, &    !
        62.d0,      0.d0,     0.d0,      -34.d0,    0.d0,    0.d0, &    !
       -63.d0,      0.d0,     0.d0,       33.d0,    0.d0,    0.d0, &    !
       -73.d0,      0.d0,     0.d0,       32.d0,    0.d0,    0.d0, &    !
       115.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
      -103.d0,      0.d0,     0.d0,        2.d0,    0.d0,    0.d0, &    !
        63.d0,      0.d0,     0.d0,      -28.d0,    0.d0,    0.d0, &    !
        74.d0,      0.d0,     0.d0,      -32.d0,    0.d0,    0.d0, &    !
      -103.d0,      0.d0,    -3.d0,        3.d0,    0.d0,   -1.d0, &    !
       -69.d0,      0.d0,     0.d0,       30.d0,    0.d0,    0.d0, &    !
        57.d0,      0.d0,     0.d0,      -29.d0,    0.d0,    0.d0, &    !
        94.d0,      0.d0,     0.d0,       -4.d0,    0.d0,    0.d0, &    !
        64.d0,      0.d0,     0.d0,      -33.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=241,260 ) / &
       -63.d0,      0.d0,     0.d0,       26.d0,    0.d0,    0.d0, &    !
       -38.d0,      0.d0,     0.d0,       20.d0,    0.d0,    0.d0, &    !
       -43.d0,      0.d0,     0.d0,       24.d0,    0.d0,    0.d0, &    !
       -45.d0,      0.d0,     0.d0,       23.d0,    0.d0,    0.d0, &    !
        47.d0,      0.d0,     0.d0,      -24.d0,    0.d0,    0.d0, &    !
       -48.d0,      0.d0,     0.d0,       25.d0,    0.d0,    0.d0, &    !
        45.d0,      0.d0,     0.d0,      -26.d0,    0.d0,    0.d0, &    !
        56.d0,      0.d0,     0.d0,      -25.d0,    0.d0,    0.d0, &    !
        88.d0,      0.d0,     0.d0,        2.d0,    0.d0,    0.d0, &    !
       -75.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
        85.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
        49.d0,      0.d0,     0.d0,      -26.d0,    0.d0,    0.d0, &    !
       -74.d0,      0.d0,    -3.d0,       -1.d0,    0.d0,   -1.d0, &    !
       -39.d0,      0.d0,     0.d0,       21.d0,    0.d0,    0.d0, &    !
        45.d0,      0.d0,     0.d0,      -20.d0,    0.d0,    0.d0, &    !
        51.d0,      0.d0,     0.d0,      -22.d0,    0.d0,    0.d0, &    !
       -40.d0,      0.d0,     0.d0,       21.d0,    0.d0,    0.d0, &    !
        41.d0,      0.d0,     0.d0,      -21.d0,    0.d0,    0.d0, &    !
       -42.d0,      0.d0,     0.d0,       24.d0,    0.d0,    0.d0, &    !
       -51.d0,      0.d0,     0.d0,       22.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=261,280 ) / &
       -42.d0,      0.d0,     0.d0,       22.d0,    0.d0,    0.d0, &    !
        39.d0,      0.d0,     0.d0,      -21.d0,    0.d0,    0.d0, &    !
        46.d0,      0.d0,     0.d0,      -18.d0,    0.d0,    0.d0, &    !
       -53.d0,      0.d0,     0.d0,       22.d0,    0.d0,    0.d0, &    !
        82.d0,      0.d0,     0.d0,       -4.d0,    0.d0,    0.d0, &    !
        81.d0,      0.d0,    -1.d0,       -4.d0,    0.d0,    0.d0, &    !
        47.d0,      0.d0,     0.d0,      -19.d0,    0.d0,    0.d0, &    !
        53.d0,      0.d0,     0.d0,      -23.d0,    0.d0,    0.d0, &    !
       -45.d0,      0.d0,     0.d0,       22.d0,    0.d0,    0.d0, &    !
       -44.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
       -33.d0,      0.d0,     0.d0,       16.d0,    0.d0,    0.d0, &    !
       -61.d0,      0.d0,     0.d0,        1.d0,    0.d0,    0.d0, &    !
       -38.d0,      0.d0,     0.d0,       19.d0,    0.d0,    0.d0, &    !
       -33.d0,      0.d0,     0.d0,       21.d0,    0.d0,    0.d0, &    !
       -60.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
        48.d0,      0.d0,     0.d0,      -10.d0,    0.d0,    0.d0, &    !
        38.d0,      0.d0,     0.d0,      -20.d0,    0.d0,    0.d0, &    !
        31.d0,      0.d0,     0.d0,      -13.d0,    0.d0,    0.d0, &    !
       -32.d0,      0.d0,     0.d0,       15.d0,    0.d0,    0.d0, &    !
        45.d0,      0.d0,     0.d0,       -8.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=281,300 ) / &
       -44.d0,      0.d0,     0.d0,       19.d0,    0.d0,    0.d0, &    !
       -51.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -36.d0,      0.d0,     0.d0,       20.d0,    0.d0,    0.d0, &    !
        44.d0,      0.d0,     0.d0,      -19.d0,    0.d0,    0.d0, &    !
       -60.d0,      0.d0,     0.d0,        2.d0,    0.d0,    0.d0, &    !
        35.d0,      0.d0,     0.d0,      -18.d0,    0.d0,    0.d0, &    !
        47.d0,      0.d0,     0.d0,       -1.d0,    0.d0,    0.d0, &    !
        36.d0,      0.d0,     0.d0,      -15.d0,    0.d0,    0.d0, &    !
       -36.d0,      0.d0,     0.d0,       20.d0,    0.d0,    0.d0, &    !
       -35.d0,      0.d0,     0.d0,       19.d0,    0.d0,    0.d0, &    !
       -37.d0,      0.d0,     0.d0,       19.d0,    0.d0,    0.d0, &    !
        32.d0,      0.d0,     0.d0,      -16.d0,    0.d0,    0.d0, &    !
        35.d0,      0.d0,     0.d0,      -14.d0,    0.d0,    0.d0, &    !
        32.d0,      0.d0,     0.d0,      -13.d0,    0.d0,    0.d0, &    !
        65.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
        47.d0,      0.d0,     0.d0,       -1.d0,    0.d0,    0.d0, &    !
        32.d0,      0.d0,     0.d0,      -16.d0,    0.d0,    0.d0, &    !
        37.d0,      0.d0,     0.d0,      -16.d0,    0.d0,    0.d0, &    !
       -30.d0,      0.d0,     0.d0,       15.d0,    0.d0,    0.d0, &    !
       -32.d0,      0.d0,     0.d0,       16.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=301,320 ) / &
       -31.d0,      0.d0,     0.d0,       13.d0,    0.d0,    0.d0, &    !
        37.d0,      0.d0,     0.d0,      -16.d0,    0.d0,    0.d0, &    !
        31.d0,      0.d0,     0.d0,      -13.d0,    0.d0,    0.d0, &    !
        49.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
        32.d0,      0.d0,     0.d0,      -13.d0,    0.d0,    0.d0, &    !
       -43.d0,      0.d0,     0.d0,       18.d0,    0.d0,    0.d0, &    !
       -32.d0,      0.d0,     0.d0,       14.d0,    0.d0,    0.d0, &    !
        30.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -34.d0,      0.d0,     0.d0,       15.d0,    0.d0,    0.d0, &    !
       -36.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -38.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -31.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -34.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
       -35.d0,      0.d0,     0.d0,        0.d0,    0.d0,    0.d0, &    !
        30.d0,      0.d0,     0.d0,       -2.d0,    0.d0,    0.d0, &    !
         0.d0,      0.d0, -1988.d0,        0.d0,    0.d0,-1679.d0, &    !
         0.d0,      0.d0,   -63.d0,        0.d0,    0.d0,  -27.d0, &    !
         0.d0,      0.d0,   364.d0,        0.d0,    0.d0,  176.d0, &    !
         0.d0,      0.d0, -1044.d0,        0.d0,    0.d0, -891.d0, &    !
         0.d0,      0.d0,   330.d0,        0.d0,    0.d0,    0.d0/      !
data ( ( cls(i,j), i=1,6 ), j=321,323 ) / &
         0.d0,      0.d0,    30.d0,        0.d0,    0.d0,   14.d0, &    !
         0.d0,      0.d0,  -162.d0,        0.d0,    0.d0, -138.d0, &    !
         0.d0,      0.d0,    75.d0,        0.d0,    0.d0,    0.d0/      !

!
!  Planetary argument multipliers:
!              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
!
data ( ( napl(i,j), i=1,14), j=  1, 20 ) / &
         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2, &
         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1, &
         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0, &
         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0/
data ( ( napl(i,j), i=1,14), j= 21, 40 ) / &
         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0, &
        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0, &
         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1, &
         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0, &
        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0/
data ( ( napl(i,j), i=1,14), j= 41, 60 ) / &
         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0, &
        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2, &
         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0, &
         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0, &
        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1/
data ( ( napl(i,j), i=1,14), j= 61, 80 ) / &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2, &
        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0, &
         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1, &
         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0/
data ( ( napl(i,j), i=1,14), j= 81,100 ) / &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2, &
         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2/
data ( ( napl(i,j), i=1,14), j=101,120 ) / &
         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1/
data ( ( napl(i,j), i=1,14), j=121,140 ) / &
         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2/
data ( ( napl(i,j), i=1,14), j=141,160 ) / &
         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2/
data ( ( napl(i,j), i=1,14), j=161,165 ) / &
         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2, &
         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2, &
        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0, &
         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0, &
        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0/

!
!  Planetary nutation coefficients, unit 1e-7 arcsec:
!  longitude (sin, cos), obliquity (sin, cos)
!
data ( ( cpl(i,j), i=1,4 ), j=  1, 20 ) / &
    1440.d0,       0.d0,       0.d0,       0.d0, &
      56.d0,    -117.d0,     -42.d0,     -40.d0, &
     125.d0,     -43.d0,       0.d0,     -54.d0, &
    -114.d0,       0.d0,       0.d0,      61.d0, &
    -219.d0,      89.d0,       0.d0,       0.d0, &
    -462.d0,    1604.d0,       0.d0,       0.d0, &
      99.d0,       0.d0,       0.d0,     -53.d0, &
      14.d0,    -218.d0,     117.d0,       8.d0, &
      31.d0,    -481.d0,    -257.d0,     -17.d0, &
    -491.d0,     128.d0,       0.d0,       0.d0, &
   -3084.d0,    5123.d0,    2735.d0,    1647.d0, &
   -1444.d0,    2409.d0,   -1286.d0,    -771.d0, &
     103.d0,     -60.d0,       0.d0,       0.d0, &
     -26.d0,     -29.d0,     -16.d0,      14.d0, &
     284.d0,       0.d0,       0.d0,    -151.d0, &
     226.d0,     101.d0,       0.d0,       0.d0, &
     -41.d0,     175.d0,      76.d0,      17.d0, &
     425.d0,     212.d0,    -133.d0,     269.d0, &
    1200.d0,     598.d0,     319.d0,    -641.d0, &
     235.d0,     334.d0,       0.d0,       0.d0/
data ( ( cpl(i,j), i=1,4 ), j= 21, 40 ) / &
     266.d0,     -78.d0,       0.d0,       0.d0, &
    -460.d0,    -435.d0,    -232.d0,     246.d0, &
       0.d0,     131.d0,       0.d0,       0.d0, &
     -42.d0,      20.d0,       0.d0,       0.d0, &
     -10.d0,     233.d0,       0.d0,       0.d0, &
      78.d0,     -18.d0,       0.d0,       0.d0, &
      45.d0,     -22.d0,       0.d0,       0.d0, &
      89.d0,     -16.d0,      -9.d0,     -48.d0, &
    -349.d0,     -62.d0,       0.d0,       0.d0, &
     -53.d0,       0.d0,       0.d0,       0.d0, &
     -21.d0,     -78.d0,       0.d0,       0.d0, &
      20.d0,     -70.d0,     -37.d0,     -11.d0, &
      32.d0,      15.d0,      -8.d0,      17.d0, &
     174.d0,      84.d0,      45.d0,     -93.d0, &
      11.d0,      56.d0,       0.d0,       0.d0, &
     -66.d0,     -12.d0,      -6.d0,      35.d0, &
      47.d0,       8.d0,       4.d0,     -25.d0, &
      46.d0,      66.d0,      35.d0,     -25.d0, &
     -68.d0,     -34.d0,     -18.d0,      36.d0, &
      76.d0,      17.d0,       9.d0,     -41.d0/
data ( ( cpl(i,j), i=1,4 ), j= 41, 60 ) / &
      84.d0,     298.d0,     159.d0,     -45.d0, &
     -82.d0,     292.d0,     156.d0,      44.d0, &
     -73.d0,      17.d0,       9.d0,      39.d0, &
    -439.d0,       0.d0,       0.d0,       0.d0, &
      57.d0,     -28.d0,     -15.d0,     -30.d0, &
     -40.d0,      57.d0,      30.d0,      21.d0, &
     273.d0,      80.d0,      43.d0,    -146.d0, &
    -449.d0,     430.d0,       0.d0,       0.d0, &
      -8.d0,     -47.d0,     -25.d0,       4.d0, &
       6.d0,      47.d0,      25.d0,      -3.d0, &
     -48.d0,    -110.d0,     -59.d0,      26.d0, &
      51.d0,     114.d0,      61.d0,     -27.d0, &
    -133.d0,       0.d0,       0.d0,      57.d0, &
     -18.d0,    -436.d0,    -233.d0,       9.d0, &
      35.d0,      -7.d0,       0.d0,       0.d0, &
     -53.d0,      -9.d0,      -5.d0,      28.d0, &
     -50.d0,     194.d0,     103.d0,      27.d0, &
     -13.d0,      52.d0,      28.d0,       7.d0, &
     -91.d0,     248.d0,       0.d0,       0.d0, &
       6.d0,      49.d0,      26.d0,      -3.d0/
data ( ( cpl(i,j), i=1,4 ), j= 61, 80 ) / &
      -6.d0,     -47.d0,     -25.d0,       3.d0, &
      52.d0,      23.d0,      10.d0,     -23.d0, &
    -138.d0,       0.d0,       0.d0,       0.d0, &
      54.d0,       0.d0,       0.d0,     -29.d0, &
     -37.d0,      35.d0,      19.d0,      20.d0, &
    -145.d0,      47.d0,       0.d0,       0.d0, &
     -10.d0,      40.d0,      21.d0,       5.d0, &
      11.d0,     -49.d0,     -26.d0,      -7.d0, &
   -2150.d0,       0.d0,       0.d0,     932.d0, &
      85.d0,       0.d0,       0.d0,     -37.d0, &
     -86.d0,     153.d0,       0.d0,       0.d0, &
     -51.d0,       0.d0,       0.d0,      22.d0, &
     -11.d0,    -268.d0,    -116.d0,       5.d0, &
      31.d0,       6.d0,       3.d0,     -17.d0, &
     140.d0,      27.d0,      14.d0,     -75.d0, &
      57.d0,      11.d0,       6.d0,     -30.d0, &
     -14.d0,     -39.d0,       0.d0,       0.d0, &
     -25.d0,      22.d0,       0.d0,       0.d0, &
      42.d0,     223.d0,     119.d0,     -22.d0, &
     -27.d0,    -143.d0,     -77.d0,      14.d0/
data ( ( cpl(i,j), i=1,4 ), j= 81,100 ) / &
       9.d0,      49.d0,      26.d0,      -5.d0, &
   -1166.d0,       0.d0,       0.d0,     505.d0, &
     117.d0,       0.d0,       0.d0,     -63.d0, &
       0.d0,      31.d0,       0.d0,       0.d0, &
       0.d0,     -32.d0,     -17.d0,       0.d0, &
      50.d0,       0.d0,       0.d0,     -27.d0, &
      30.d0,      -3.d0,      -2.d0,     -16.d0, &
       8.d0,     614.d0,       0.d0,       0.d0, &
    -127.d0,      21.d0,       9.d0,      55.d0, &
     -20.d0,      34.d0,       0.d0,       0.d0, &
      22.d0,     -87.d0,       0.d0,       0.d0, &
     -68.d0,      39.d0,       0.d0,       0.d0, &
       3.d0,      66.d0,      29.d0,      -1.d0, &
     490.d0,       0.d0,       0.d0,    -213.d0, &
     -22.d0,      93.d0,      49.d0,      12.d0, &
     -46.d0,      14.d0,       0.d0,       0.d0, &
      25.d0,     106.d0,      57.d0,     -13.d0, &
    1485.d0,       0.d0,       0.d0,       0.d0, &
      -7.d0,     -32.d0,     -17.d0,       4.d0, &
      30.d0,      -6.d0,      -2.d0,     -13.d0/
data ( ( cpl(i,j), i=1,4 ), j=101,120 ) / &
     118.d0,       0.d0,       0.d0,     -52.d0, &
     -28.d0,      36.d0,       0.d0,       0.d0, &
      14.d0,     -59.d0,     -31.d0,      -8.d0, &
    -458.d0,       0.d0,       0.d0,     198.d0, &
       0.d0,     -45.d0,     -20.d0,       0.d0, &
    -166.d0,     269.d0,       0.d0,       0.d0, &
     -78.d0,      45.d0,       0.d0,       0.d0, &
      -5.d0,     328.d0,       0.d0,       0.d0, &
   -1223.d0,     -26.d0,       0.d0,       0.d0, &
    -368.d0,       0.d0,       0.d0,       0.d0, &
     -75.d0,       0.d0,       0.d0,       0.d0, &
     -13.d0,     -30.d0,       0.d0,       0.d0, &
     -74.d0,       0.d0,       0.d0,      32.d0, &
    -262.d0,       0.d0,       0.d0,     114.d0, &
     202.d0,       0.d0,       0.d0,     -87.d0, &
      -8.d0,      35.d0,      19.d0,       5.d0, &
     -35.d0,     -48.d0,     -21.d0,      15.d0, &
      12.d0,      55.d0,      29.d0,      -6.d0, &
    -598.d0,       0.d0,       0.d0,       0.d0, &
       8.d0,     -31.d0,     -16.d0,      -4.d0/
data ( ( cpl(i,j), i=1,4 ), j=121,140 ) / &
     113.d0,       0.d0,       0.d0,     -49.d0, &
      83.d0,      15.d0,       0.d0,       0.d0, &
       0.d0,    -114.d0,     -49.d0,       0.d0, &
     117.d0,       0.d0,       0.d0,     -51.d0, &
     393.d0,       3.d0,       0.d0,       0.d0, &
      18.d0,     -29.d0,     -13.d0,      -8.d0, &
       8.d0,      34.d0,      18.d0,      -4.d0, &
      89.d0,       0.d0,       0.d0,       0.d0, &
      54.d0,     -15.d0,      -7.d0,     -24.d0, &
       0.d0,      35.d0,       0.d0,       0.d0, &
    -154.d0,     -30.d0,     -13.d0,      67.d0, &
      80.d0,     -71.d0,     -31.d0,     -35.d0, &
      61.d0,     -96.d0,     -42.d0,     -27.d0, &
     123.d0,    -415.d0,    -180.d0,     -53.d0, &
       0.d0,       0.d0,       0.d0,     -35.d0, &
       7.d0,     -32.d0,     -17.d0,      -4.d0, &
     -89.d0,       0.d0,       0.d0,      38.d0, &
       0.d0,     -86.d0,     -19.d0,      -6.d0, &
    -123.d0,    -416.d0,    -180.d0,      53.d0, &
     -62.d0,     -97.d0,     -42.d0,      27.d0/
data ( ( cpl(i,j), i=1,4 ), j=141,160 ) / &
     -85.d0,     -70.d0,     -31.d0,      37.d0, &
     163.d0,     -12.d0,      -5.d0,     -72.d0, &
     -63.d0,     -16.d0,      -7.d0,      28.d0, &
     -21.d0,     -32.d0,     -14.d0,       9.d0, &
       5.d0,    -173.d0,     -75.d0,      -2.d0, &
      74.d0,       0.d0,       0.d0,     -32.d0, &
      83.d0,       0.d0,       0.d0,       0.d0, &
    -339.d0,       0.d0,       0.d0,     147.d0, &
      67.d0,     -91.d0,     -39.d0,     -29.d0, &
      30.d0,     -18.d0,      -8.d0,     -13.d0, &
       0.d0,    -114.d0,     -50.d0,       0.d0, &
     517.d0,      16.d0,       7.d0,    -224.d0, &
     143.d0,      -3.d0,      -1.d0,     -62.d0, &
      50.d0,       0.d0,       0.d0,     -22.d0, &
      59.d0,       0.d0,       0.d0,       0.d0, &
     370.d0,      -8.d0,       0.d0,    -160.d0, &
      34.d0,       0.d0,       0.d0,     -15.d0, &
     -37.d0,      -7.d0,      -3.d0,      16.d0, &
      40.d0,       0.d0,       0.d0,       0.d0, &
    -184.d0,      -3.d0,      -1.d0,      80.d0/
data ( ( cpl(i,j), i=1,4 ), j=161,165 ) / &
      31.d0,      -6.d0,       0.d0,     -13.d0, &
      -3.d0,     -32.d0,     -14.d0,       1.d0, &
     -34.d0,       0.d0,       0.d0,       0.d0, &
     126.d0,     -63.d0,     -27.d0,     -55.d0, &
    -126.d0,     -63.d0,     -27.d0,      55.d0/

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and given date (JC).
t = ( ( date1-dj0 ) + date2 ) / djc

!  -------------------
!  LUNI-SOLAR NUTATION
!  -------------------

!
!  Fundamental (Delaunay) arguments from Simon et al. (1994)
!
call funarg ( t,   el, elp, f, d, om )

!  Initialize the nutation values.
dp = 0.d0
de = 0.d0

!  Summation of luni-solar nutation series (in reverse order).
do i = nls, 1, -1

!     Argument and functions.
   arg = mod ( dble ( nals(1,i) ) * el  + &
               dble ( nals(2,i) ) * elp + &
               dble ( nals(3,i) ) * f   + &
               dble ( nals(4,i) ) * d   + &
               dble ( nals(5,i) ) * om, d2pi )
   sarg = sin(arg)
   carg = cos(arg)

!     Term.
   dp = dp + ( cls(1,i) + cls(2,i) * t ) * sarg &
           +   cls(3,i)                  * carg
   de = de + ( cls(4,i) + cls(5,i) * t ) * carg &
           +   cls(6,i)                  * sarg

end do

!  Convert from 0.1 microarcsec units to radians.
dpsils = dp * u2r
depsls = de * u2r

!  ------------------
!  PLANETARY NUTATION
!  ------------------

!  Planetary longitudes, Mercury through Neptune, wrt mean dynamical
!  ecliptic and equinox of J2000, with high order terms omitted
!  (Simon et al. 1994, 5.8.1-5.8.8).
alme = mod ( 4.402608842461d0 + 2608.790314157421d0 * t, d2pi )
alve = mod ( 3.176146696956d0 + 1021.328554621099d0 * t, d2pi )
alea = mod ( 1.753470459496d0 +  628.307584999142d0 * t, d2pi )
alma = mod ( 6.203476112911d0 +  334.061242669982d0 * t, d2pi )
alju = mod ( 0.599547105074d0 +   52.969096264064d0 * t, d2pi )
alsa = mod ( 0.874016284019d0 +   21.329910496032d0 * t, d2pi )
alur = mod ( 5.481293871537d0 +    7.478159856729d0 * t, d2pi )
alne = mod ( 5.311886286677d0 +    3.813303563778d0 * t, d2pi )

!  General precession in longitude (Simon et al. 1994), equivalent
!  to 5028.8200 arcsec/cy at J2000.
apa = ( 0.024380407358d0 + 0.000005391235d0 * t ) * t

!  Initialize the nutation values.
dp = 0.d0
de = 0.d0

!  Summation of planetary nutation series (in reverse order).
do i = npl, 1, -1

!     Argument and functions.
   arg = mod ( dble ( napl( 1,i) ) * el   + &
               dble ( napl( 2,i) ) * elp  + &
               dble ( napl( 3,i) ) * f    + &
               dble ( napl( 4,i) ) * d    + &
               dble ( napl( 5,i) ) * om   + &
               dble ( napl( 6,i) ) * alme + &
               dble ( napl( 7,i) ) * alve + &
               dble ( napl( 8,i) ) * alea + &
               dble ( napl( 9,i) ) * alma + &
               dble ( napl(10,i) ) * alju + &
               dble ( napl(11,i) ) * alsa + &
               dble ( napl(12,i) ) * alur + &
               dble ( napl(13,i) ) * alne + &
               dble ( napl(14,i) ) * apa, d2pi )
   sarg = sin(arg)
   carg = cos(arg)

!     Term.
   dp = dp + cpl(1,i) * sarg + cpl(2,i) * carg
   de = de + cpl(3,i) * sarg + cpl(4,i) * carg

end do

!  Convert from 0.1 microarcsec units to radians.
dpsipl = dp * u2r
depspl = de * u2r

!  -----
!  TOTAL
!  -----

!  Add planetary and luni-solar components.
dpsi = dpsipl + dpsils
deps = depspl + depsls

end
!***********************************************************************

!***********************************************************************
!>
!  Equation of the equinoxes complementary terms, consistent with
!  IAU 2000 resolutions.
!
!  Annexe to IERS Conventions 2000, Chapter 5
!
!  Capitaine, N., Wallace, P.T., & McCarthy, D.D. (2003). Astron. &
!    Astrophys. 406, pp. 1135-1149, Table 3.
!  IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
!    (Table 5.2e presented in the printed publication is a truncated
!    series. The full series, which is used in NOVAS, is available on
!    the IERS Conventions Center website in file tab5.2e.txt.)
!    ftp://tai.bipm.org/iers/conv2010/chapter5/
!
!  Given:
!     DATE1,DATE2   d    TT date (JD = DATE1+DATE2)
!
!  Returned:
!     EECT00        d    complementary terms (radians)
!
!  This revision:  2002 November 13
!                  References updated 2010 November 26

double precision function eect2000 ( date1, date2 )

implicit none

double precision date1, date2

!  2Pi
double precision d2pi
parameter ( d2pi = 6.283185307179586476925287d0 )

!  Arcseconds to radians
double precision das2r
parameter ( das2r = 4.848136811095359935899141d-6 )

!  Reference epoch (J2000), JD
double precision dj0
parameter ( dj0 = 2451545d0 )

!  Days per Julian century
double precision djc
parameter ( djc = 36525d0 )

!  Time since J2000, in Julian centuries
double precision t

!  Miscellaneous
integer i, j
double precision a, s0, s1
!double precision anmp

!  Fundamental arguments
double precision fa(14)

!  -----------------------------------------
!  The series for the EE complementary terms
!  -----------------------------------------

!  Number of terms in the series
integer ne0, ne1
parameter ( ne0=  33, ne1=  1 )

!  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
integer ke0 ( 14, ne0 ), &
        ke1 ( 14, ne1 )

!  Sine and cosine coefficients
double precision se0 ( 2, ne0 ), &
                 se1 ( 2, ne1 )

!  Argument coefficients for t^0
data ( ( ke0(i,j), i=1,14), j =    1,   10 ) / &
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
data ( ( ke0(i,j), i=1,14), j =   11,   20 ) / &
  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0, &
  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
data ( ( ke0(i,j), i=1,14), j =   21,   30 ) / &
  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1, &
  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
data ( ( ke0(i,j), i=1,14), j =   31,  ne0 ) / &
  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

!  Argument coefficients for t^1
data ( ( ke1(i,j), i=1,14), j =    1,  ne1 ) / &
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

!  Sine and cosine coefficients for t^0
data ( ( se0(i,j), i=1,2), j =    1,   10 ) / &
            +2640.96d-6,          -0.39d-6, &
              +63.52d-6,          -0.02d-6, &
              +11.75d-6,          +0.01d-6, &
              +11.21d-6,          +0.01d-6, &
               -4.55d-6,          +0.00d-6, &
               +2.02d-6,          +0.00d-6, &
               +1.98d-6,          +0.00d-6, &
               -1.72d-6,          +0.00d-6, &
               -1.41d-6,          -0.01d-6, &
               -1.26d-6,          -0.01d-6 /
data ( ( se0(i,j), i=1,2), j =   11,   20 ) / &
               -0.63d-6,          +0.00d-6, &
               -0.63d-6,          +0.00d-6, &
               +0.46d-6,          +0.00d-6, &
               +0.45d-6,          +0.00d-6, &
               +0.36d-6,          +0.00d-6, &
               -0.24d-6,          -0.12d-6, &
               +0.32d-6,          +0.00d-6, &
               +0.28d-6,          +0.00d-6, &
               +0.27d-6,          +0.00d-6, &
               +0.26d-6,          +0.00d-6 /
data ( ( se0(i,j), i=1,2), j =   21,   30 ) / &
               -0.21d-6,          +0.00d-6, &
               +0.19d-6,          +0.00d-6, &
               +0.18d-6,          +0.00d-6, &
               -0.10d-6,          +0.05d-6, &
               +0.15d-6,          +0.00d-6, &
               -0.14d-6,          +0.00d-6, &
               +0.14d-6,          +0.00d-6, &
               -0.14d-6,          +0.00d-6, &
               +0.14d-6,          +0.00d-6, &
               +0.13d-6,          +0.00d-6 /
data ( ( se0(i,j), i=1,2), j =   31,  ne0 ) / &
               -0.11d-6,          +0.00d-6, &
               +0.11d-6,          +0.00d-6, &
               +0.11d-6,          +0.00d-6 /

!  Sine and cosine coefficients for t^1
data ( ( se1(i,j), i=1,2), j =    1,  ne1 ) / &
               -0.87d-6,          +0.00d-6 /

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Interval between fundamental epoch J2000.0 and current date (JC).
t = ( ( date1-dj0 ) + date2 ) / djc

!  Fundamental Arguments (from IERS Conventions 2000)

!  Mean Anomaly of the Moon.
fa(1) = anmp ( ( 485868.249036d0 + &
               ( 715923.2178d0 + &
               (     31.8792d0 + &
               (      0.051635d0 + &
               (     -0.00024470d0 ) &
               * t ) * t ) * t ) * t ) * das2r &
               + mod ( 1325d0*t, 1d0 ) * d2pi )

!  Mean Anomaly of the Sun.
fa(2) = anmp ( ( 1287104.793048d0 + &
               ( 1292581.0481d0 + &
               (      -0.5532d0 + &
               (      +0.000136d0 + &
               (      -0.00001149d0 ) &
               * t ) * t ) * t ) * t ) * das2r &
               + mod ( 99d0*t, 1d0 ) * d2pi )

!  Mean Longitude of the Moon minus Mean Longitude of the Ascending
!  Node of the Moon.
fa(3) = anmp ( (  335779.526232d0 + &
               (  295262.8478d0 + &
               (     -12.7512d0 + &
               (      -0.001037d0 + &
               (       0.00000417d0 ) &
               * t ) * t ) * t ) * t ) * das2r &
               + mod ( 1342d0*t, 1d0 ) * d2pi )

!  Mean Elongation of the Moon from the Sun.
fa(4) = anmp ( ( 1072260.703692d0 + &
               ( 1105601.2090d0 + &
               (      -6.3706d0 + &
               (       0.006593d0 + &
               (      -0.00003169d0 ) &
               * t ) * t ) * t ) * t ) * das2r &
               + mod ( 1236d0*t, 1d0 ) * d2pi )

!  Mean Longitude of the Ascending Node of the Moon.
fa(5) = anmp ( (  450160.398036d0 + &
               ( -482890.5431d0 + &
               (       7.4722d0 + &
               (       0.007702d0 + &
               (      -0.00005939d0 ) &
               * t ) * t ) * t ) * t ) * das2r &
               + mod ( -5d0*t, 1d0 ) * d2pi )

fa( 6) = anmp ( 4.402608842d0 + 2608.7903141574d0 * t )
fa( 7) = anmp ( 3.176146697d0 + 1021.3285546211d0 * t )
fa( 8) = anmp ( 1.753470314d0 +  628.3075849991d0 * t )
fa( 9) = anmp ( 6.203480913d0 +  334.0612426700d0 * t )
fa(10) = anmp ( 0.599546497d0 +   52.9690962641d0 * t )
fa(11) = anmp ( 0.874016757d0 +   21.3299104960d0 * t )
fa(12) = anmp ( 5.481293872d0 +    7.4781598567d0 * t )
fa(13) = anmp ( 5.311886287d0 +    3.8133035638d0 * t )
fa(14) =      ( 0.024381750d0 +    0.00000538691d0 * t ) * t

!  Evaluate the EE complementary terms.
s0 = 0d0
s1 = 0d0

do i = ne0,1,-1
   a = 0d0
   do j=1,14
      a = a + dble(ke0(j,i))*fa(j)
   end do
   s0 = s0 + ( se0(1,i)*sin(a) + se0(2,i)*cos(a) )
end do
do i = ne1,1,-1
   a = 0d0
   do j=1,14
      a = a + dble(ke1(j,i))*fa(j)
   end do
   s1 = s1 + ( se1(1,i)*sin(a) + se1(2,i)*cos(a) )
end do
eect2000 = ( s0 + s1 * t ) * das2r

end function eect2000
!***********************************************************************

!***********************************************************************
!>
!  Normalize angle into the range -pi <= A < +pi.

    double precision function anmp ( a )

    implicit none

    double precision a

    double precision dpi, d2pi
    parameter ( dpi = 3.141592653589793238462643d0, &
                d2pi = 6.283185307179586476925287d0 )

    double precision w

    w = mod(a,d2pi)
    if ( abs(w) >= dpi ) w = w - sign(d2pi,a)
    anmp = w

    end function anmp
!***********************************************************************

!***********************************************************************
  end module novas_module
!***********************************************************************