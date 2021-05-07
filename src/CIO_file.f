*
*     C I O _ F I L E
*
*     PROGRAM TO PRODUCE A BINARY DIRECT ACCESS FILE OF RIGHT ASCENSION
*     VALUES OF THE CELESTIAL INTERMEDIATE ORIGIN (CIO), GIVEN A
*     FORMATTED SEQUENTIAL FILE OF THE SAME DATA.  EACH INPUT AND
*     OUTPUT DATA RECORD CONTAINS A TDB JULIAN DATE AND A
*     RIGHT ASCENSION VALUE (WRT GCRS) IN ARCSECONDS.
*
*
      PROGRAM CIO_FILE 

      DOUBLE PRECISION TDBJD, CIORA
      CHARACTER INFIL*24, OUTFIL*24, IDEN*40

   1  FORMAT ( A )
   2  FORMAT ( F16.6,  F24.14 )

*     GET FILE IDENTIFIERS 
      INFIL  = 'CIO_RA.TXT'
      OUTFIL = 'CIO_RA.DA'
C      WRITE ( *, * ) 'ENTER INPUT FILENAME: '
C      READ ( *, * ) INFIL
C      WRITE ( *, * ) 'ENTER OUTPUT FILENAME: '
C      READ ( *, * ) OUTFIL

*     OPEN INPUT FILE
      OPEN ( UNIT=18, FILE=INFIL, FORM='FORMATTED',
     .     ACCESS='SEQUENTIAL', STATUS='OLD' )

*     OPEN OUTPUT FILE
      OPEN ( UNIT=19, FILE=OUTFIL, FORM='UNFORMATTED',
     .     ACCESS='DIRECT', RECL=16, STATUS='UNKNOWN' )

*     READ INPUT FILE IDENTIFIER
      READ ( UNIT=18, FMT=1 ) IDEN

      N = 1

*     MAIN READ-WRITE LOOP
  50  READ ( UNIT=18, FMT=2, END=70 ) TDBJD, CIORA
      N = N + 1
      WRITE ( UNIT=19, REC=N ) TDBJD, CIORA

      IF ( MOD ( N, 1000 ) .EQ. 0 ) THEN
          IYEAR = NINT ( ( TDBJD - 2451545.0 ) / 365.25D0 + 2000.D0 )
          WRITE ( *, * ) 'DONE THROUGH RECORD ', N, '      YEAR ', IYEAR
      END IF

      GO TO 50

*     WRITE NUMBER OF DATA RECORDS IN FIRST RECORD OF OUTPUT FILE
*     ALONG WITH FIRST 12 CHARACTERS OF INPUT FILE IDENTIFIER
  70  WRITE ( UNIT=19, REC=1 ) N-1, IDEN(1:12)

*     FINISH UP
      WRITE ( *, * ) N, ' TOTAL RECORDS WRITTEN'
      CLOSE ( UNIT=18 )
      CLOSE ( UNIT=19 )

      STOP

      END PROGRAM CIO_FILE
