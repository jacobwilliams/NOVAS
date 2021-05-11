!  NOVAS FORTRAN VERS F3.1
!  ALTERNATIVE VERSIONS OF SOME SUBROUTINES

!***********************************************************************
!                                                                      *
!                              N O V A S                               *
!           NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE               *
!                                                                      *
!                            G. H. KAPLAN                              *
!                        U.S. NAVAL OBSERVATORY                        *
!                                                                      *
!***********************************************************************


!***********************************************************************
!>
!  SUBROUTINE GRVDEF VERSION 2.
!  THIS SUBROUTINE COMPUTES THE TOTAL GRAVITATIONAL DEFLECTION OF
!  LIGHT FOR THE OBSERVED OBJECT DUE TO THE MAJOR GRAVITATING BODIES
!  IN THE SOLAR SYSTEM.
!  THIS VERSION IS A DUMMY.  NO CORRECTION IS APPLIED.
!
!       TJD    = (NOT USED)
!       LOC    = (NOT USED)
!       POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
!                TO ICRS AXES, COMPONENTS IN AU (IN)
!       POBS   = (NOT USED)
!       POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
!                ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
!                TO ICRS AXES, CORRECTED FOR GRAVITATIONAL
!                DEFLECTION, COMPONENTS IN AU (OUT)
!
!  NOTE:  IN THIS VERSION, POS2 = POS1.  THE GRAVITATIONAL DEFLECTION
!  THAT IS NEGLECTED HERE CAN REACH 1.8 ARCSECONDS AT THE LIMB OF
!  THE SUN, BUT IS LESS THAN 10 MILLIARCSECONDS OVER THE AREA OF SKY
!  MORE THAN 45 DEGREES FROM THE SUN.  SEE TABLE 3.26.1 ON PAGE 138
!  OF THE EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC (1992).

subroutine grvdef (tjd,loc,pos1,pobs,pos2)

double precision tjd,pos1,pobs,pos2
dimension pos1(3), pobs(3), pos2(3)

do j = 1, 3
    pos2(j) = pos1(j)
end do

end subroutine grvdef
!***********************************************************************

!***********************************************************************
!>
!  SUBROUTINE NOD VERSION 2.
!  IN LOW-ACCURACY MODE, THIS SUBROUTINE EVALUATES A SHORT
!  NUTATION SERIES AND RETURNS APPROXIMATE VALUES FOR NUTATION IN
!  LONGITUDE AND NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
!  IN THIS MODE, ONLY THE LARGEST 13 TERMS OF THE IAU 2000A NUTATION
!  SERIES ARE EVALUATED.  IN HIGH-ACCURACY MODE, THE STANDARD IERS
!  SUBROUINE IS CALLED TO EVALUATE THE FULL IAU 2000A NUTATION
!  SERIES.
!
!       T    = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
!       DPSI = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
!       DEPS = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
!
!  NOTE:  IN LOW-ACCURACY MODE, MAX ERROR IN DPSI < 0.05 ARCSEC,
!  MAX ERROR IN DEPS < 0.02 ARCSEC, AVERAGE ERROR ABOUT 1/4 OF MAX.

subroutine nod (t,dpsi,deps)

double precision t,dpsi,deps,pi,seccon,t0,t1,dp,de, &
     x,el,elp,f,d,om,arg,dsin,dcos
dimension x(9,13)
save

parameter ( pi     = 3.14159265358979324d0 )
parameter ( seccon = 180.d0 * 3600.d0 / pi )

! T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
data t0 / 2451545.00000000d0 /

! LARGEST 13 TERMS OF IAU 2000A NUTATION SERIES, WITH PRECISION
! OF COEFFICIENTS TRUNCATED
data x / 0., 0., 0., 0., 1., -17.2064,-0.01747, 9.2052, 0.00091, &
         0., 0., 2.,-2., 2.,  -1.3171,-0.00017, 0.5730,-0.00030, &
         0., 0., 2., 0., 2.,  -0.2276,-0.00002, 0.0978,-0.00005, &
         0., 0., 0., 0., 2.,   0.2075, 0.00002,-0.0897, 0.00005, &
         0., 1., 0., 0., 0.,   0.1476,-0.00036, 0.0074,-0.00002, &
         0., 1., 2.,-2., 2.,  -0.0517, 0.00012, 0.0224,-0.00007, &
         1., 0., 0., 0., 0.,   0.0711, 0.00001,-0.0007, 0.00000, &
         0., 0., 2., 0., 1.,  -0.0387,-0.00004, 0.0201, 0.00000, &
         1., 0., 2., 0., 2.,  -0.0301, 0.00000, 0.0129,-0.00001, &
         0.,-1., 2.,-2., 2.,   0.0216,-0.00005,-0.0096, 0.00003, &
         0., 0., 2.,-2., 1.,   0.0128, 0.00001,-0.0069,-0.00000, &
        -1., 0., 2., 0., 2.,   0.0123, 0.00000,-0.0053, 0.00000, &
        -1., 0., 0., 2., 0.,   0.0157, 0.00000,-0.0001, 0.00000 /
! REMAINING TERMS ALL HAVE AMPLITUDES < 0.01 ARCSECOND

! GET METHOD/ACCURACY MODE
call getmod ( mode )

if ( mod ( mode, 2 ) == 0 ) then

    ! HIGH ACCURACY MODE -- CALL IERS SUBROUTINE

    t1 = t * 36525.d0
    call nu2000a ( t0, t1, dp, de )
    dpsi = dp * seccon
    deps = de * seccon

else

    ! LOW ACCURACY MODE -- EVALUATE SHORT NUTATION SERIES ABOVE

    ! COMPUTATION OF FUNDAMENTAL ARGUMENTS
    call funarg ( t,   el, elp, f, d, om )

    dpsi = 0.d0
    deps = 0.d0

    ! SUM NUTATION SERIES TERMS
    do i = 13, 1, -1
        arg = x(1,i) * el &
            + x(2,i) * elp &
            + x(3,i) * f &
            + x(4,i) * d &
            + x(5,i) * om
        dpsi = ( x(6,i) + x(7,i) * t ) * dsin ( arg ) + dpsi
        deps = ( x(8,i) + x(9,i) * t ) * dcos ( arg ) + deps
    end do

    ! ADD IN OUT-OF-PHASE COMPONENT OF PRINCIPAL (18.6-YEAR) TERM
    ! (TO AVOID SMALL BUT LONG-TERM BIAS IN RESULTS)
    dpsi = dpsi + 0.0033d0 * dcos ( om )
    deps = deps + 0.0015d0 * dsin ( om )

end if

end
!***********************************************************************