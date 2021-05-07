!
!     C I O _ F I L E
!
!     PROGRAM TO PRODUCE A BINARY DIRECT ACCESS FILE OF RIGHT ASCENSION
!     VALUES OF THE CELESTIAL INTERMEDIATE ORIGIN (CIO), GIVEN A
!     FORMATTED SEQUENTIAL FILE OF THE SAME DATA.  EACH INPUT AND
!     OUTPUT DATA RECORD CONTAINS A TDB JULIAN DATE AND A
!     RIGHT ASCENSION VALUE (WRT GCRS) IN ARCSECONDS.
!
!
program cio_file 

double precision tdbjd, ciora
character infil*24, outfil*24, iden*40

1 format ( a )
2 format ( f16.6,  f24.14 )

!     GET FILE IDENTIFIERS 
infil  = 'CIO_RA.TXT'
outfil = 'CIO_RA.DA'
!      WRITE ( *, * ) 'ENTER INPUT FILENAME: '
!      READ ( *, * ) INFIL
!      WRITE ( *, * ) 'ENTER OUTPUT FILENAME: '
!      READ ( *, * ) OUTFIL

!     OPEN INPUT FILE
open ( unit=18, file=infil, form='FORMATTED', &
     access='SEQUENTIAL', status='OLD' )

!     OPEN OUTPUT FILE
open ( unit=19, file=outfil, form='UNFORMATTED', &
     access='DIRECT', recl=16, status='UNKNOWN' )

!     READ INPUT FILE IDENTIFIER
read ( unit=18, fmt=1 ) iden

n = 1

!     MAIN READ-WRITE LOOP
50 read ( unit=18, fmt=2, end=70 ) tdbjd, ciora
n = n + 1
write ( unit=19, rec=n ) tdbjd, ciora

if ( mod ( n, 1000 ) .eq. 0 ) then
    iyear = nint ( ( tdbjd - 2451545.0 ) / 365.25d0 + 2000.d0 )
    write ( *, * ) 'DONE THROUGH RECORD ', n, '      YEAR ', iyear      !
end if

go to 50

!     WRITE NUMBER OF DATA RECORDS IN FIRST RECORD OF OUTPUT FILE
!     ALONG WITH FIRST 12 CHARACTERS OF INPUT FILE IDENTIFIER
70 write ( unit=19, rec=1 ) n-1, iden(1:12)

!     FINISH UP
write ( *, * ) n, ' TOTAL RECORDS WRITTEN'
close ( unit=18 )
close ( unit=19 )

stop

end program cio_file

