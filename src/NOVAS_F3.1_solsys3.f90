!** SOLSYS VERSION 3 PACKAGE: SOLSYS, SUN, IDSS ***

!***********************************************************************
!>
!  SUBROUTINE SOLSYS VERSION 3.
!  THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF THE
!  EARTH AT EPOCH TJD BY EVALUATING A CLOSED-FORM THEORY WITHOUT
!  REFERENCE TO AN EXTERNAL FILE.  THIS ROUTINE CAN ALSO PROVIDE
!  THE POSITION AND VELOCITY OF THE SUN.
!
!       TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
!       M    = BODY IDENTIFICATION NUMBER (IN)
!              SET M=0 OR M=1 FOR THE SUN
!              SET M=2 OR M=3 FOR THE EARTH
!       K    = ORIGIN SELECTION CODE (IN)
!              SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
!              SET K=1 FOR ORIGIN AT CENTER OF SUN
!       POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
!              OF J2000.0, COMPONENTS IN AU (OUT)
!       VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
!              OF J2000.0, COMPONENTS IN AU/DAY (OUT)
!       IERR = ERROR INDICATOR (OUT)
!              IERR=0 MEANS EVERYTHING OK
!              IERR=1 MEANS TJD BEFORE FIRST ALLOWED DATE
!              IERR=2 MEANS TJD AFTER LAST ALLOWED DATE

subroutine solsys (tjd,m,k,pos,vel,ierr)

 use novas_module, only: preces

 double precision tjd,pos,vel,pi,twopi,t0,obl,el,c,p,tlast, &
     pm,pa,pe,pj,po,pw,pl,pn, &
     tmass,se,ce,si,ci,sn,cn,sw,cw,p1,p2,p3,q1,q2,q3,roote,a,b, &
     qjd,e,mlon,ma,u,sinu,cosu,anr,pplan,vplan,f,pbary,vbary, &
     dfloat,dabs,dmod,dsin,dcos,dsqrt
dimension pos(3), vel(3), el(21), c(13), p(3,3), &
     pm(4), pa(4), pe(4), pj(4), po(4), pw(4), pl(4), pn(4), &
     a(3,4), b(3,4), pplan(3), vplan(3), pbary(3), vbary(3)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( twopi  = 2.d0 * pi             )
parameter ( t0     = 2451545.0d0           )
parameter ( obl    = 23.43927944d0         )
! T0 = TDB JULIAN DATE OF EPOCH J2000.0
! OBL = OBLIQUITY OF ECLIPTIC AT EPOCH J2000.0

data el, c, p / 43*0.d0 /,   tlast / 0.d0 /

! ARRAYS BELOW CONTAIN MASSES AND ORBITAL ELEMENTS OF THE FOUR
! LARGEST PLANETS (SEE EXPLANATORY SUPPLEMENT (1992), P. 316)
! WITH ANGLES IN RADIANS
! THIS DATA USED FOR BARYCENTER COMPUTATIONS ONLY
!             JUPITER        SATURN        URANUS       NEPTUNE
data pm /  1047.349d+0,  3497.898d+0,   22903.0d+0,   19412.2d+0 /      !
data pa /  5.203363d+0,  9.537070d+0, 19.191264d+0, 30.068963d+0 /      !
data pe /  0.048393d+0,  0.054151d+0,  0.047168d+0,  0.008586d+0 /      !
data pj /  0.022782d+0,  0.043362d+0,  0.013437d+0,  0.030878d+0 /      !
data po /  1.755036d+0,  1.984702d+0,  1.295556d+0,  2.298977d+0 /      !
data pw /  0.257503d+0,  1.613242d+0,  2.983889d+0,  0.784898d+0 /      !
data pl /  0.600470d+0,  0.871693d+0,  5.466933d+0,  5.321160d+0 /      !
data pn /  1.450138d-3,  5.841727d-4,  2.047497d-4,  1.043891d-4 /      !

if ( tlast < 1.d0 ) then
    ! first time computations
    ! mass of sun plus four inner planets
    tmass = 1.d0 + 5.977d-6
    se = dsin ( obl * pi / 180.d0 )
    ce = dcos ( obl * pi / 180.d0 )
    do i = 1, 4
        tmass = tmass + 1.d0 / pm(i)
        ! compute sine and cosine of orbital angles
        si = dsin ( pj(i) )
        ci = dcos ( pj(i) )
        sn = dsin ( po(i) )
        cn = dcos ( po(i) )
        sw = dsin ( pw(i) - po(i) )
        cw = dcos ( pw(i) - po(i) )
        ! compute p and q vectors (see brouwer & clemence (1961),
        ! methods of celestial mechanics, pp. 35-36.)
        p1 =    cw * cn - sw * sn * ci
        p2 = (  cw * sn + sw * cn * ci ) * ce - sw * si * se
        p3 = (  cw * sn + sw * cn * ci ) * se + sw * si * ce
        q1 =   -sw * cn - cw * sn * ci
        q2 = ( -sw * sn + cw * cn * ci ) * ce - cw * si * se
        q3 = ( -sw * sn + cw * cn * ci ) * se + cw * si * ce
        roote = dsqrt ( 1.d0 - pe(i)**2 )
        a(1,i) = pa(i) * p1
        a(2,i) = pa(i) * p2
        a(3,i) = pa(i) * p3
        b(1,i) = pa(i) * roote * q1
        b(2,i) = pa(i) * roote * q2
        b(3,i) = pa(i) * roote * q3
    end do
    tlast = 1.d0
end if

ierr = 0
! VALID DATES ARE WITHIN 3 CENTURIES OF J2000, ALTHOUGH RESULTS
! DETERIORATE GRADUALLY
if ( tjd < 2340000.5d0 ) ierr = 1
if ( tjd > 2560000.5d0 ) ierr = 2
if ( ierr /= 0 ) return
if ( m >= 2 ) then
    ! FORM HELIOCENTRIC COORDINATES OF EARTH
    ! VELOCITIES ARE OBTAINED FROM CRUDE NUMERICAL DIFFERENTIATION
    do i = 1, 3
        qjd = tjd + dfloat(i-2) * 0.1d0
        ! SUBROUTINE SUN COMPUTES EARTH-SUN VECTOR
        call sun ( qjd, el, c )
        call preces ( qjd, c(11), t0, pos )
        p(i,1) = -pos(1)
        p(i,2) = -pos(2)
        p(i,3) = -pos(3)
    end do
    do j=1,3
        pos(j) =   p(2,j)
        vel(j) = ( p(3,j) - p(1,j) ) / 0.2d0
    end do
    if ( k >= 1 ) return
else
    ! FORM HELIOCENTRIC COORDINATES OF SUN
    do j=1,3
        pos(j) = 0.d0
        vel(j) = 0.d0
    end do
    if ( k >= 1 ) return
end if

! IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
! SOLAR SYSTEM BARYCENTER COORDINATES COMPUTED FROM KEPLERIAN
! APPROXIMATIONS OF THE COORDINATES OF THE FOUR LARGEST PLANETS
if ( dabs ( tjd - tlast ) >= 1.d-6 ) then
    do j = 1, 3
        pbary(j) = 0.d0
        vbary(j) = 0.d0
    end do
    ! THE FOLLOWING LOOP CYCLES ONCE FOR EACH OF THE FOUR LARGE PLANETS
    do i = 1, 4
        ! COMPUTE MEAN LONGITUDE, MEAN ANOMALY, AND ECCENTRIC ANOMOLY
        e = pe(i)
        mlon = pl(i) + pn(i) * ( tjd - t0 )
        ma = dmod ( mlon - pw(i), twopi )
        u = ma + e * dsin ( ma ) + 0.5d0 * e * e * dsin ( 2.d0 * ma )
        sinu = dsin ( u )
        cosu = dcos ( u )
        ! COMPUTE VELOCITY FACTOR
        anr = pn(i) / ( 1.d0 - e * cosu )
        ! COMPUTE PLANET'S POSITION AND VELOCITY WRT EQ & EQ J2000
        pplan(1) = a(1,i) * ( cosu - e ) + b(1,i) * sinu
        pplan(2) = a(2,i) * ( cosu - e ) + b(2,i) * sinu
        pplan(3) = a(3,i) * ( cosu - e ) + b(3,i) * sinu
        vplan(1) = anr * ( -a(1,i) * sinu + b(1,i) * cosu )
        vplan(2) = anr * ( -a(2,i) * sinu + b(2,i) * cosu )
        vplan(3) = anr * ( -a(3,i) * sinu + b(3,i) * cosu )
        ! COMPUTE MASS FACTOR AND ADD IN TO TOTAL DISPLACEMENT
        f = 1.d0 / ( pm(i) * tmass )
        pbary(1) = pbary(1) + pplan(1) * f
        pbary(2) = pbary(2) + pplan(2) * f
        pbary(3) = pbary(3) + pplan(3) * f
        vbary(1) = vbary(1) + vplan(1) * f
        vbary(2) = vbary(2) + vplan(2) * f
        vbary(3) = vbary(3) + vplan(3) * f
    end do
    tlast = tjd
end if
do j=1,3
    pos(j) = pos(j) - pbary(j)
    vel(j) = vel(j) - vbary(j)
end do

end subroutine solsys
!***********************************************************************

!***********************************************************************
!>
!  FOR USE WITH SUBROUTINE SOLSYS VERSION 3.
!  THIS SUBROUTINE COMPUTES THE COORDINATES OF THE EARTH-SUN
!  POSITION VECTOR WITH RESPECT TO THE ECLIPTIC AND EQUATOR
!  OF DATE.  A MODIFIED FORM OF NEWCOMB'S THEORY ('TABLES OF THE
!  SUN', 1898) IS USED.  ONLY THE LARGEST PERIODIC PERTURBATIONS
!  ARE EVALUATED, AND VAN FLANDERN'S EXPRESSIONS FOR THE FUNDAMENTAL
!  ARGUMENTS ('IMPROVED MEAN ELEMENTS FOR THE EARTH AND MOON', 1981)
!  ARE USED.  THE ABSOLUTE ACCURACY IS NO WORSE THAN 1 ARCSECOND
!  (AVERAGE ERROR ABOUT 0.2 ARCSECOND) OVER 1800-2200.
!  (ADAPTED FROM SUBROUTINE IAUSUN BY P. M. JANICZEK, USNO.)
!
!       DJ   = TDB JULIAN DATE OF DESIRED EPOCH (IN)
!       EL   = ARRAY OF ORBITAL ELEMENTS (SEE BELOW) FOR
!              EPOCH DJ (OUT)
!       C    = ARRAY OF COORDINATES (SEE BELOW) FOR
!              EPOCH DJ (OUT)

subroutine sun (dj,el,c)

double precision dj,el,c,t,tp,t20,ro,gv,gm,gj,gs,dl,dr,db,dg, &
 dblarg,d,twopi,str,rtd,r,tr, &
 sino,coso,sinl,cosl,sinb,cosb, &
 dsin,dcos,dmod
!
dimension el(21)
!
!     EL( 1)= SEMI-MAJOR AXIS, AU
!     EL( 2)= ORBITAL ECCENTRICITY
!     EL( 5)= LONGITUDE OF PERIGEE, RADIANS
!     EL( 9)= UNPERTURBED MEAN LONGITUDE, RADIANS
!     EL(10)= MEAN ANOMALY, AFFECTED BY LONG-PD PERTURBATIONS, RADIANS
!     EL(11)= UNPERTURBED RADIUS, AU
!     EL(12)= EQUATION OF THE CENTER, RADIANS
!     EL(13)= MEAN OBLIQUITY OF ECLIPTIC, RADIANS
!     EL(14)= MEAN LONGITUDE OF MOON, RADIANS
!     EL(15)= MEAN ANOMALY OF MOON, RADIANS
!     EL(16)= LUNAR MEAN ARGUMENT OF LATITUDE, RADIANS
!     EL(17)= MEAN LONGITUDE OF LUNAR ASCENDING NODE, RADIANS
!     EL(21)= MEAN LONGITUDE OF MOON'S PERIGEE, RADIANS
!             (REMAINING ELEMENTS OF ARRAY EL NOT USED)
!
dimension c(13)
!
!     C( 1) = PERTURBED RADIUS VECTOR, AU
!     C( 2) = SAME AS C(4), DEGREES
!     C( 3) = SAME AS C(5), DEGREES
!     C( 4) = ECLIPTIC LONGITUDE WRT MEAN ECL & EQUX OF DATE, RADIANS
!     C( 5) = ECLIPTIC LATITUDE  WRT MEAN ECL        OF DATE, RADIANS
!     C(11) = EQUATORIAL X WRT MEAN EQU & EQUX OF DATE, AU
!     C(12) = EQUATORIAL Y WRT MEAN EQU & EQUX OF DATE, AU
!     C(13) = EQUATORIAL Z WRT MEAN EQU & EQUX OF DATE, AU
!             (REMAINING ELEMENTS OF ARRAY C NOT USED)
!
!
!***********************************************************************
!
!     PART I    TABLES OF THE PERTURBATIONS
!
dimension x(8,46), x1(80), x2(80), x3(80), x4(80), x5(48)
equivalence (x(1, 1),x1(1))
equivalence (x(1,11),x2(1))
equivalence (x(1,21),x3(1))
equivalence (x(1,31),x4(1))
equivalence (x(1,41),x5(1))
!
!     PERTURBATIONS BY VENUS
!                  J    I     VC      VS    RHOC    RHOS      BC     BS
data x1 /  - 1.,  0., +  33.,-  67., -  85.,-  39., +  24.,-  17., &    !
           - 1.,+ 1., +2353.,-4228., -2062.,-1146., -   4.,+   3., &    !
           - 1.,+ 2., -  65.,-  34., +  68.,-  14., +   6.,-  92., &    !
           - 2.,+ 1., -  99.,+  60., +  84.,+ 136., +  23.,-   3., &    !
           - 2.,+ 2., -4702.,+2903., +3593.,+5822., +  10.,-   6., &    !
           - 2.,+ 3., +1795.,-1737., - 596.,- 632., +  37.,-  56., &    !
           - 3.,+ 3., - 666.,+  27., +  44.,+1044., +   8.,+   1., &    !
           - 3.,+ 4., +1508.,- 397., - 381.,-1448., + 185.,- 100., &    !
           - 3.,+ 5., + 763.,- 684., + 126.,+ 148., +   6.,-   3., &    !
           - 4.,+ 4., - 188.,-  93., - 166.,+ 337.,     0.,    0./      !
data x2 /  - 4.,+ 5., - 139.,-  38., -  51.,+ 189., -  31.,-   1., &    !
           - 4.,+ 6., + 146.,-  42., -  25.,-  91., +  12.,    0., &    !
           - 5.,+ 5., -  47.,-  69., - 134.,+  93.,     0.,    0., &    !
           - 5.,+ 7., - 119.,-  33., -  37.,+ 136., -  18.,-   6., &    !
           - 5.,+ 8., + 154.,    0.,     0.,-  26.,     0.,    0., &    !
           - 6.,+ 6., -   4.,-  38., -  80.,+   8.,     0.,    0., &    !
!
!     PERTURBATIONS BY MARS
!                  J    I     VC      VS    RHOC    RHOS      BC     BS
           + 1.,- 1., - 216.,- 167., -  92.,+ 119.,     0.,    0., &    !
           + 2.,- 2., +1963.,- 567., - 573.,-1976.,     0.,-   8., &    !
           + 2.,- 1., -1659.,- 617., +  64.,- 137.,     0.,    0., &    !
           + 3.,- 3., +  53.,- 118., - 154.,-  67.,     0.,    0./      !
data x3 /  + 3.,- 2., + 396.,- 153., -  77.,- 201.,     0.,    0., &    !
           + 4.,- 3., - 131.,+ 483., + 461.,+ 125., +   7.,+   1., &    !
           + 4.,- 2., + 526.,- 256., +  43.,+  96.,     0.,    0., &    !
           + 5.,- 4., +  49.,+  69., +  87.,-  62.,     0.,    0., &    !
           + 5.,- 3., -  38.,+ 200., +  87.,+  17.,     0.,    0., &    !
           + 6.,- 4., - 104.,- 113., - 102.,+  94.,     0.,    0., &    !
           + 6.,- 3., -  11.,+ 100., -  27.,-   4.,     0.,    0., &    !
           + 7.,- 4., -  78.,-  72., -  26.,+  28.,     0.,    0., &    !
           + 9.,- 5., +  60.,-  15., -   4.,-  17.,     0.,    0., &    !
           +15.,- 8., + 200.,-  30., -   1.,-   6.,     0.,    0./      !
!
!     PERTURBATIONS BY JUPITER
!                  J    I     VC      VS    RHOC    RHOS      BC     BS
data x4 /  + 1.,- 2., - 155.,-  52., -  78.,+ 193., +   7.,    0., &    !
           + 1.,- 1., -7208.,+  59., +  56.,+7067., -   1.,+  17., &    !
           + 1.,  0., - 307.,-2582., + 227.,-  89., +  16.,    0., &    !
           + 1.,+ 1., +   8.,-  73., +  79.,+   9., +   1.,+  23., &    !
           + 2.,- 3., +  11.,+  68., + 102.,-  17.,     0.,    0., &    !
           + 2.,- 2., + 136.,+2728., +4021.,- 203.,     0.,    0., &    !
           + 2.,- 1., - 537.,+1518., +1376.,+ 486., +  13.,+ 166., &    !
           + 3.,- 3., - 162.,+  27., +  43.,+ 278.,     0.,    0., &    !
           + 3.,- 2., +  71.,+ 551., + 796.,- 104., +   6.,-   1., &    !
           + 3.,- 1., -  31.,+ 208., + 172.,+  26., +   1.,+  18./      !
data x5 /  + 4.,- 3., -  43.,+   9., +  13.,+  73.,     0.,    0., &    !
           + 4.,- 2., +  17.,+  78., + 110.,-  24.,     0.,    0., &    !
!
!     PERTURBATIONS BY SATURN
!                  J    I     VC      VS    RHOC    RHOS      BC     BS
           + 1.,- 1., -  77.,+ 412., + 422.,+  79., +   1.,+   6., &    !
           + 1.,  0., -   3.,- 320., +   8.,-   1.,     0.,    0., &    !
           + 2.,- 2., +  38.,- 101., - 152.,-  57.,     0.,    0., &    !
           + 2.,- 1., +  45.,- 103., - 103.,-  44.,     0.,    0./      !
!
!
!***********************************************************************
!
!     PART II   NECESSARY PRELIMINARIES
!
data twopi /6.283185307179586d0/
data str   /206264806.2470964d0/
data rtd   /57.295779513082321d0/
data r     /1296000.0d0/
tr = 1000.0d0 / str
!
!     T  = TIME IN JULIAN CENTURIES FROM 1900 JANUARY 0
t  = (dj - 2415020.d0)/36525.d0
!
!     TP = TIME IN JULIAN YEARS     FROM 1850 JANUARY 0
tp = (dj - 2396758.d0)/365.25d0
!
!     T20= TIME IN JULIAN CENTURIES FROM J2000.0
t20= (dj - 2451545.d0)/36525.d0
!
!
!***********************************************************************
!
!     PART III  COMPUTATION OF ELLIPTIC ELEMENTS AND SECULAR TERMS
!
!     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
el( 1) = 1.00000030007166d0
el( 2) = 0.016708320d0 + (-0.42229d-04 - 0.126d-06 * t20) * t20
el( 5) = 1018578.046d0 + (6190.046d0 + &
                (1.666d0 + 0.012d0 * t20) * t20) * t20
el( 5) = el( 5) * tr
el( 9) = 1009677.850d0 + (100.0d0 * r + 2771.27d0 + &
                1.089d0 * t20) * t20
el( 9) = dmod (el( 9) * tr, twopi)
el(10) = 1287099.804d0 + (99.0d0 * r + 1292581.224d0 + &
                (-0.577d0 - 0.012d0 * t20) * t20) * t20
el(10) = dmod (el(10) * tr, twopi)

!     EXPRESSION FOR OBLIQUITY FROM P03 (IAU 2006) PRECESSION
el(13) = 84381.406d0 + (-46.836769d0 + &
               (-0.0001831d0 + 0.00200340d0 * t20) * t20) * t20
el(13) = el(13) * tr

!     KAPLAN CORRECTION TO SUN'S MEAN LONGITUDE TO FIT DE405 OVER
!     INTERVAL 1800-2200, USING P03 (IAU 2006) PRECESSION
el(9) = el(9) &
      + ( 0.1320d0 - 0.1355d0 * t20 ) * tr

!
!***********************************************************************
!
!     PART IV   LUNAR TERMS
!
!     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
el(14) = 785939.157d0 + (1336.0d0 * r + 1108372.598d0 &
                + (-5.802d0 + 0.019d0 * t20) * t20) * t20
el(14) = dmod (el(14) * tr, twopi)
el(17) = 450160.280d0 + (-5.0d0 * r - 482890.539d0 + &
                (7.455d0 + 0.008d0 * t20) * t20) * t20
el(17) = dmod (el(17) * tr, twopi)
el(21) = 300072.424d0 + (11.0d0 * r + 392449.965d0 + &
                (-37.112d0 - 0.045d0 * t20) * t20) * t20
el(21) = dmod (el(21) * tr, twopi)
!
!     DERIVED ARGUMENTS
el(15) = el(14) - el(21)
el(16) = el(14) - el(17)
el(15) = dmod (el(15),twopi)
el(16) = dmod (el(16),twopi)
!     MEAN ELONGATION
d      = el(14) - el(9)
!
!     COMBINATIONS OF ARGUMENTS AND THE PERTURBATIONS
d = dmod (d,twopi)
arg = d
dl =    +  6469.*sin(arg) +  13.*sin(3.*arg)
dr =    + 13390.*cos(arg) +  30.*cos(3.*arg)
!
dblarg = d + el(15)
dblarg = dmod (dblarg,twopi)
arg = dblarg
dl = dl +  177.*sin(arg)
dr = dr +  370.*cos(arg)
!
dblarg = d - el(15)
dblarg = dmod (dblarg,twopi)
arg = dblarg
dl = dl -  424.*sin(arg)
dr = dr - 1330.*cos(arg)
!
dblarg = 3.d0*d - el(15)
dblarg = dmod (dblarg,twopi)
arg = dblarg
dl = dl +   39.*sin(arg)
dr = dr +   80.*cos(arg)
!
dblarg = d + el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
dl = dl -   64.*sin(arg)
dr = dr -  140.*cos(arg)
!
dblarg = d - el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
dl = dl +  172.*sin(arg)
dr = dr +  360.*cos(arg)
!
el(16) = dmod (el(16),twopi)
arg = el(16)
db =    + 576.*sin(arg)
!
!
!***********************************************************************
!
!     PART V    COMPUTATION OF PERIODIC PERTURBATIONS
!
!     THE PERTURBING MEAN ANOMALIES
!
gv  = 0.19984020d+01 + .1021322923d+02*tp
gm  = 0.19173489d+01 + .3340556174d+01*tp
gj  = 0.25836283d+01 + .5296346478d+00*tp
gs  = 0.49692316d+01 + .2132432808d+00*tp
gv  = dmod (gv,twopi)
gm  = dmod (gm,twopi)
gj  = dmod (gj,twopi)
gs  = dmod (gs,twopi)
!
!
!     MODIFICATION OF FUNDAMENTAL ARGUMENTS
!
!     APPLICATION OF THE JUPITER-SATURN GREAT INEQUALITY
!     TO JUPITER'S MEAN ANOMALY
!
gj = gj + 0.579904067d-02 * dsin (5.d0*gs - 2.d0*gj &
                 + 1.1719644977d0 - 0.397401726d-03*tp)
gj = dmod (gj,twopi)
!
!     LONG PERIOD PERTURBATIONS OF MEAN ANOMALY
!
st = t
!                ARGUMENT IS ( 4 MARS - 7 EARTH + 3 VENUS )
dg = 266.* sin (0.555015 + 2.076942*st) &
!                ARGUMENT IS ( 3 JUPITER - 8 MARS + 4 EARTH )
    + 6400.* sin (4.035027 + 0.3525565*st) &
!                ARGUMENT IS ( 13 EARTH - 8 VENUS )
    + (1882.-16.*st) * sin (0.9990265 + 2.622706*st)
!
!
!     COMPUTATION OF THE EQUATION OF THE CENTER
!
!     FORM PERTURBED MEAN ANOMALY
el(10) = dg/str + el(10)
el(10) = dmod (el(10),twopi)
el(12) =   dsin(     el(10)) * (6910057.d0 -(17240.d0+52.d0*t)*t) &
         + dsin(2.d0*el(10)) * (  72338.d0 -    361.d0*t) &
         + dsin(3.d0*el(10)) * (   1054.d0 -      1.d0*t)
!
!     THE UNPERTURBED RADIUS VECTOR
ro     =                          30570.d0 -    150.d0*t &
         - dcos(     el(10)) * (7274120.d0 - (18140.d0+50.d0*t)*t) &    !
         - dcos(2.d0*el(10)) * (  91380.d0 -    460.d0*t) &
         - dcos(3.d0*el(10)) * (   1450.d0 -     10.d0*t)
el(11) = 10.d0**(ro*1.d-09)
!
!
!     SELECTED PLANETARY PERTURBATIONS FROM NEWCOMB'S THEORY FOLLOW
!
!     PERTURBATIONS BY VENUS
do 20 k=1,16
!     ARGUMENT J * VENUS +   I * EARTH
dblarg = x(1,k)*gv + x(2,k)*el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
cs  = cos(arg)
ss  = sin(arg)
dl  =(x(3,k)*cs  + x(4,k)*ss )+ dl
dr  =(x(5,k)*cs  + x(6,k)*ss )+ dr
db  =(x(7,k)*cs  + x(8,k)*ss )+ db
20 continue
!
!     PERTURBATIONS BY MARS
do 30 k=17,30
!     ARGUMENT  J * MARS +   I * EARTH
dblarg = x(1,k)*gm + x(2,k)*el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
cs  = cos(arg)
ss  = sin(arg)
dl  =(x(3,k)*cs  + x(4,k)*ss )+ dl
dr  =(x(5,k)*cs  + x(6,k)*ss )+ dr
db  =(x(7,k)*cs  + x(8,k)*ss )+ db
30 continue
!
!     PERTURBATIONS BY JUPITER
do 40 k=31,42
!     ARGUMENT J*JUPITER +   I * EARTH
dblarg = x(1,k)*gj + x(2,k)*el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
cs  = cos(arg)
ss  = sin(arg)
dl  =(x(3,k)*cs  + x(4,k)*ss )+ dl
dr  =(x(5,k)*cs  + x(6,k)*ss )+ dr
db  =(x(7,k)*cs  + x(8,k)*ss )+ db
40 continue
!
!     PERTURBATIONS BY SATURN
do 50 k=43,46
!     ARGUMENT J*SATURN  +   I * EARTH
dblarg = x(1,k)*gs + x(2,k)*el(10)
dblarg = dmod (dblarg,twopi)
arg = dblarg
cs  = cos(arg)
ss  = sin(arg)
dl  =(x(3,k)*cs  + x(4,k)*ss )+ dl
dr  =(x(5,k)*cs  + x(6,k)*ss )+ dr
db  =(x(7,k)*cs  + x(8,k)*ss )+ db
50 continue
!
!
!***********************************************************************
!
!     PART VI   COMPUTATION OF ECLIPTIC AND EQUATORIAL COORDINATES
!
c(1) = el(11)*10.d0**(dr*1.d-09)
c(4) = (dl + dg + el(12))/str + el(9)
c(4) = dmod (c(4),twopi)
c(5) = db/str
c(2) = c(4)*rtd
c(3) = c(5)*rtd
sino = dsin (el(13))
coso = dcos (el(13))
sinl = dsin (c(4))
cosl = dcos (c(4))
sinb = dsin (c(5))
cosb = dcos (c(5))
c(11) = c(1) * (cosb * cosl)
c(12) = c(1) * (cosb * sinl * coso - sinb * sino)
c(13) = c(1) * (cosb * sinl * sino + sinb * coso)

end subroutine sun
!***********************************************************************

!***********************************************************************
!>
!  THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
!  FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
!
!      NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
!               'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
!               UPPER-CASE LETTERS (IN)
!      IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
!               (FUNCTION VALUE RETURNED)
!
!  NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
!  BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
!  MIGHT USE MORE LETTERS.
!
!  NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=1, SINCE SOLSYS
!  VERSION 3 DOES NOT PROCESS SPLIT JULIAN DATES.
!
!  NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
!  THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.

integer function idss ( name )

character name*(*), namein*3, names*3
dimension names(35), ids(35)

data names / 'SUN', 'EAR', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---'  /
data ids   /     0,     3,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0  /
data num   / 2 /

3 format ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', a )

idss = -9999
namein = name

!     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
do i = 1, num
    if ( namein == names(i) ) then
        idss = ids(i)
        return
    end if
end do

!     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES
if ( namein == 'JD ' ) then
!         IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
!         JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE
    idss = 1
    return
end if

write ( *, 3 ) name

end function idss
!***********************************************************************