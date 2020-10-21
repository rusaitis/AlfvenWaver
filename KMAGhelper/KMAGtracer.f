      program KMAG_PROGRAM
      ! ---------------------------------------------------------------
      !       _  ____  __    _    ____ _
      !      | |/ /  \/  |  / \  / ___| |_ _ __ __ _  ___ ___ _ __
      !      | ' /| |\/| | / _ \| |  _| __| '__/ _` |/ __/ _ \ '__|
      !      | . \| |  | |/ ___ \ |_| | |_| | | (_| | (_|  __/ |
      !      |_|\_\_|  |_/_/   \_\____|\__|_|  \__,_|\___\___|_|
      !
      ! Written by:               Liutauras (Leo) Rusaitis
      !                           10-13-2020
      !
      !                           Space Physics PhD Student,
      !                           Earth, Planetary, and Space Sciences,
      !                           University of California, Los Angeles
      !                           GitHub: https://github.com/rusaitis
      !                           Contact: rusaitis@ucla.edu
      ! ---------------------------------------------------------------
      !
      ! Tracing procedures interafacing with the KMAG magnetic field
      ! model.
      !
      ! ---------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 BY_IMF, BZ_IMF, Dp
      real*8 TIME, R, TH, PH
      real*8 TIME0, R0, TH0, PH0
      real*8 BR, BTH, BPH
      real*8 dX3C, dY3C, dZ3C
      real*8 X, Y, Z
      real*8 X0, Y0, Z0
      real*8 dX, dY, dZ
      real*8 dX0, dY0, dZ0
      real*8 RLT, BRM, BTM, BPM
      real*8 STEP
      real*8 PLANET_SURFACE, MAX_R_DISTANCE
      integer n
      REAL*8 h
      REAL*8, DIMENSION(3) :: YY, XX, XOUT, XX0
      REAL*8, DIMENSION(5) :: P
      REAL*8 SQRT_VEC
      REAL*8 INVL
      character*5 epoch
      integer i, j, ITER_COUNT 
      integer STEP_COUNT, MAX_ITER_COUNT
      integer DIR
      integer KEEP_TRACING
      integer DIALOG_SWITCH
      integer LINE_NUM
      integer FLAG
      integer TRACE_NUMBER_MAX
      integer FILE_ID
      integer TRACE_OUTCOME
      Character*3 From, To
      Character*3 IN_COORD, OUT_COORD, COORD(2)
      Character*3 UNITS
      Character*30 COMMENT
      Character*1 TRACE
      Character*3 TRACE_COORD
      character(len=100) :: filename
      character*3 :: run_number
      character*30 :: OUTDIR
      character*19 :: date_str
      character(len=12) :: DESCRIPTION
      external derivs
      PARAMETER (PI = 3.1415927)

      !/////////////////////// OPEN INPUT FILE /////////////////////////
      OPEN (1, file = './KMAGhelper/kmag_params.txt')
      READ(1,*)
      READ(UNIT=1, FMT=*, iostat=io) time,BY_IMF,BZ_IMF,Dp,IN_COORD,
     +COORD(2),COMMENT,TRACE,COORD(1),STEP
      CLOSE (1)
      !/////////////////////////////////////////////////////////////////
      
      n = 3 ! NUMBER OF DIMENSIONS
      epoch = 'j2000' ! EPOCH TYPE ('j2000' | 'ctime')
      OUTDIR = 'Output'
      PLANET_SURFACE = 1. ! R_planet
      MAX_ITER_COUNT = 100000
      TRACE_COORD = 'CAR' ! ('SPH' | 'CAR')
      MAX_R_DISTANCE = 300. ! IN PLANET RADII
      DIR = 1 ! DIRECTION OF TRACING
      KEEP_TRACING = 1
      LINE_NUM = 0
      DIALOG_SWITCH = 0
      TRACE_NUMBER_MAX = 2 ! 2:TRACE FORWARD & BACK, 1: TRACE FORWARD 
      DESCRIPTION = 'EVENT'
      FILE_ID = 0
      !/////////////////////// DISCONTINUITIES /////////////////////////
      OPEN(5,FILE='./KMAGhelper/discontinuities.txt', Access = 'APPEND',
     +STATUS='UNKNOWN')
      !/////////////////////////////////////////////////////////////////

      IF (TRACE.EQ.'N') THEN
         WRITE(filename, '(I6.6)') FILE_ID
         CALL OPEN_FILE(4, OUTDIR, filename)
         CALL PRINT_HEADERS(4, TRACE_COORD)
      END IF 

      OPEN(1, FILE = './KMAGhelper/kmag_input.txt')
      READ(1,*)

      !////////////////// READ INPUT LINE AT A TIME ////////////////////
      DO
         READ(UNIT=1, FMT=*, iostat=io) TIME,XX(1),XX(2),XX(3)
         IF (io/=0) EXIT
         P(1)=TIME; P(2)=BY_IMF; P(3)=BZ_IMF; P(4)=DP; P(5)=0
   
         !////////////////////// INITIATE TRACING //////////////////////
         IF (TRACE.EQ.'Y') THEN
            IF ((COORD(1).EQ.'SPH').AND.(TRACE_COORD.EQ.'CAR')) THEN
               CALL SPH2CAR(XX, XX0)
               R=XX(1)
               CALL KROT(IN_COORD,COORD(2),XX0,XX,P(1),epoch)
               XX0 = XX
            END IF
            IF ((COORD(1).EQ.'SPH').AND.(TRACE_COORD.EQ.'SPH')) THEN
               CALL SPH2CAR(XX, XX0)
               CALL KROT(IN_COORD,COORD(2),XX0,XX,P(1),epoch)
               CALL CAR2SPH(XX, XX0)
               XX = XX0
               R=XX(1)
            END IF
            IF ((COORD(1).EQ.'CAR').AND.(TRACE_COORD.EQ.'SPH')) THEN
               R=SQRT_VEC(XX)
            END IF
            XX0 = XX
            TIME0 = TIME
            TRACE_NUMBER = 0
            KEEP_TRACING = 1
            STEP_COUNT = 0
            WRITE(filename, '(I6.6)') FILE_ID
            CALL OPEN_FILE(4, OUTDIR, filename)
            CALL PRINT_HEADERS(4, TRACE_COORD)
   
            !///////////////////////// TRACING /////////////////////////
   
            DO WHILE (KEEP_TRACING .EQ. 1)
               CALL rk4(XX,YY,n,STEP,dir,XOUT,derivs,P,COORD, 
     +         TRACE_COORD)
C                CALL euler(XX,YY,n,STEP,dir,XOUT,derivs,P,COORD,
C      +         TRACE_COORD)
               CALL CHECK_POSITION(XX0, XX, XOUT, COORD, TRACE_COORD, 
     +         DESCRIPTION)
               XX = XOUT
               !--------------- CHECK FOR DISCONTINUITY (FLAG) ---------
               !IF (FLAG.EQ.1) THEN
               !      CALL WRITE_LINE(5, P(1), XX, YY)
               !END IF
               !--------------------------------------------------------
               !                  * TRACING CONDITIONALS *
               IF (TRACE_COORD.EQ.'SPH') R=XX(1)
               IF (TRACE_COORD.EQ.'CAR') THEN
                  R=DSQRT( XX(1)**2 + XX(2)**2 +1.0* XX(3)**2 )
               END IF
               STEP_COUNT = STEP_COUNT + 1
               CALL CHECK_TRACING(KEEP_TRACING,R,STEP_COUNT,
     +         MAX_ITER_COUNT,MAX_R_DISTANCE,PLANET_SURFACE,STEP,
     +         DIALOG_SWITCH,TRACE_OUTCOME)
               IF (TRACE_OUTCOME .EQ. 1) DESCRIPTION = 'R=1'
               IF (TRACE_OUTCOME .EQ. 2) DESCRIPTION = 'R=i'
               IF (TRACE_OUTCOME .EQ. 3) DESCRIPTION = 'MAX'
               CALL WRITE_LINE(4, P(1), XX, YY, DESCRIPTION)
               IF (KEEP_TRACING .EQ. 0) THEN  
                     TRACE_NUMBER = TRACE_NUMBER + 1
                     IF (TRACE_NUMBER .LT. TRACE_NUMBER_MAX) THEN
                           TIME=TIME0; XX=XX0  !GO BACK TO START POS
                           DIR = DIR * (-1)    !CHANGE TRACE DIRECTION
                           KEEP_TRACING = 1    !CONTINUE TRACING AGAIN
                           STEP_COUNT = 0
                           DESCRIPTION = 'EVENT'
                     END IF
               END IF
            END DO ! END MAG FIELD LINE TRACE
            !///////////////////////////////////////////////////////////
   
         ELSE
            !//////////// IF NO TRACING - JUST GET MAG FIELDS //////////
            CALL EULER(XX,YY,n,STEP,dir,XOUT,derivs,P,COORD,TRACE_COORD)
            CALL WRITE_LINE(4, P(1), XX, YY, '-')
            !///////////////////////////////////////////////////////////
         END IF
         FILE_ID = FILE_ID + 1
      END DO
C******************************CLOSE THE PROGRAM************************
      CLOSE (4)
      CLOSE (1)
      END
C****************************** FLAG CHECK *****************************
      SUBROUTINE CHECK_POSITION(X0, XIN, XOUT, COORD, TRACE_COORD, 
     +DESCRIPTION)
      double precision X0(3), XIN(3), XOUT(3)
      character*3 COORD(2), TRACE_COORD
      character(len=*) :: DESCRIPTION
      PARAMETER (PI=3.141593)
      DESCRIPTION = '-'
      IF (TRACE_COORD.EQ.'SPH') THEN
         IF ((PI/2-XOUT(2))*(PI/2-XIN(2)).LE.0) THEN
            DESCRIPTION = 'PSH'
         END IF
      ELSE IF (TRACE_COORD.EQ.'CAR') THEN
         IF (XOUT(3)*XIN(3).LE.0) THEN
            DESCRIPTION = 'PSH'
         END IF
      END IF
      IF ((XIN(1).EQ.X0(1)).AND.
     +    (XIN(2).EQ.X0(2)).AND.
     +    (XIN(3).EQ.X0(3))) THEN
         DESCRIPTION = 'EVENT'
      END IF
      RETURN
      END
C************************ PRINT OUTPUT COLUMN NAMES ********************
      SUBROUTINE OPEN_FILE(FILE_NUM, OUTDIR, FILENAME)
      integer FILE_NUM, LINE_NUM
      character*30 OUTDIR
      character*30 FILENAME
      character*120 PATH
      WRITE(PATH,'(a,a,a,a)') 
     +trim(adjustl(OUTDIR)), "/", trim(adjustl(FILENAME)), ".txt"
      OPEN(FILE_NUM, FILE=trim(adjustl(PATH)),STATUS='UNKNOWN')
      RETURN
      END
C************************ PRINT OUTPUT COLUMN NAMES ********************
      SUBROUTINE PRINT_HEADERS(FILE_NUM, COORD)
      integer FILE_NUM
      character*3 COORD
      IF (COORD.EQ.'SPH') THEN
         Write(FILE_NUM,'(9A16)') 'TIME','R','TH','PHI',
     +                  'BRM','BTM','BPM','BMM','COMMENT'
      ELSE IF (COORD.EQ.'CAR') THEN
         Write(FILE_NUM,'(9A16)') 'TIME','X','Y','Z',
     +                  'BX','BY','BZ','BMM','COMMENT'
      END IF
      RETURN
      END
C****************************** WRITE A LINE ***************************
      SUBROUTINE WRITE_LINE(FILE_NUM, TIME, X, B, COMMENT)
      integer FILE_NUM
      real*8 TIME, BMM
      character*3 c
      character(len=*) :: COMMENT
C       character*3,optional COMMENT
      double precision X(3),B(3)
      BMM = DSQRT(B(1)**2+B(2)**2+B(3)**2)
      write(FILE_NUM,'(8F16.3, a16)') TIME,X(1),X(2),X(3),B(1),B(2),
     +B(3),BMM,COMMENT
      RETURN
      END
C***********************************************************************
      SUBROUTINE CHECK_TRACING(KEEP_TRACING,R,STEP_COUNT,MAX_ITER_COUNT,
     + MAX_R_DISTANCE,PLANET_SURFACE,STEP, DIALOG_SWITCH, OUTCOME)
      real*8 PLANET_SURFACE, MAX_R_DISTANCE, STEP
      real*8 R
      integer STEP_COUNT, MAX_ITER_COUNT
      integer DIR
      integer OUTCOME
      integer DIALOG_SWITCH, KEEP_TRACING
      character(len=70) :: DIALOG_TEXT
      KEEP_TRACING = 1
      OUTCOME = 0
      IF ( STEP_COUNT .GT. MAX_ITER_COUNT ) THEN
        DIALOG_TEXT='Program terminated. Max iteration count reached.'
        KEEP_TRACING = 0
        OUTCOME = 3
      ELSE IF ( R .GT. MAX_R_DISTANCE ) THEN
        DIALOG_TEXT='Program terminated. Beyond 80RS.'
        KEEP_TRACING = 0
        OUTCOME = 2 
      ELSE IF ( R .LT. PLANET_SURFACE ) THEN
        DIALOG_TEXT='Program terminated. Reached Saturn.'
        KEEP_TRACING = 0
        OUTCOME = 1 
      END IF
      IF ((DIALOG_SWITCH .EQ. 1).AND.(KEEP_TRACING .EQ. 0 )) THEN
            print*, DIALOG_TEXT
      END IF
      RETURN 
      END

C************************CONVERT INPUT KSO COORD TO S3C*****************
      FUNCTION SQRT_VEC(X)
      double precision X(3)
      real*8 SQRT_VEC
      SQRT_VEC = DSQRT( X(1)**2 + X(2)**2 + X(3)**2 )
      RETURN
      END FUNCTION
C************************CONVERT INPUT KSO COORD TO S3C*****************
      SUBROUTINE CAR2SPH(X,R)
      double precision X(3),R(3), RHO
      PARAMETER (PI=3.141593)
      R(3) = datan2(X(2),X(1))
      RHO = dsqrt(X(1)**2+X(2)**2)
      R(1) = dsqrt(X(1)**2 + X(2)**2 + X(3)**2)
      R(2) = datan2(RHO, X(3)) !COLATITUDE
      RETURN
      END
C************************CONVERT INPUT KSO COORD TO S3C*****************
      SUBROUTINE SPH2CAR(R,X)
      double precision X(3),R(3)
      PARAMETER (PI=3.141593)
      X(1)=R(1)*DSIN(R(2))*DCOS(R(3))
      X(2)=R(1)*DSIN(R(2))*DSIN(R(3))
      X(3)=R(1)*DCOS(R(2))
      RETURN
      END
C************************ 4TH ORDER RUNGE KUTTA ************************

      subroutine rk4(X,Y,N,H,DIR,XOUT,DERIVS,P,COORD, TRACE_COORD)
      implicit none
      integer N,NMAX,I
      integer dir
      double precision h,X(N),DX(N),Y(N),XOUT(N)
      double precision P(5)
      Character*3 COORD(2), TRACE_COORD
      external derivs
      parameter(nmax=50)
      double precision h6,hh,dxm(nmax),dxt(nmax),xt(nmax)
      hh = dir * h * 0.5d0
      h6 = dir * h / 6d0
      CALL DERIVS(X,  DX, Y,  N, P, COORD, TRACE_COORD)
      DO i=1,N
         XT(i) = X(i)+hh*DX(i)
      ENDDO
      CALL DERIVS(XT, DXT, Y, N, P, COORD, TRACE_COORD)
      DO i=1,N
         XT(i) = X(i)+hh*DXT(i)
      ENDDO
      CALL DERIVS(XT, DXM, Y, N, P, COORD, TRACE_COORD)
      DO i=1,N
         XT(i) = X(i)+h*DXM(i)
         DXM(i) = DXT(i)+DXM(i)
      ENDDO
      CALL DERIVS(XT, DXT, Y, N, P, COORD, TRACE_COORD)
      DO i=1,N
         XOUT(i)=X(i)+h6*(DX(i)+DXT(i)+2D0*DXM(i))
      ENDDO
      RETURN
      END

C*************************** EULER (1'ST ORDER) ************************

      subroutine euler(X,Y,N,H,DIR,XOUT,DERIVS,P,COORD, TRACE_COORD)
      implicit none
      integer N,NMAX,I
      integer dir
      double precision h,X(N),DX(N),Y(N),XOUT(N), DX0(N)
      double precision P(5)
      external derivs
      parameter(nmax=50)
      Character*3 COORD(2), TRACE_COORD
      double precision h6,hh,dxm(nmax),dxt(nmax),xt(nmax)
      CALL DERIVS(X,  DX, Y,  N, P, COORD, TRACE_COORD)
      DO i=1,N
         XOUT(i)=X(i)+dir*h*DX(i)
      ENDDO
      RETURN
      END

C SUBROUTINE DERIVS(X,Y,DYDX), WHICH RETURNS DERIVATIVES DYDX AT X.
      SUBROUTINE DERIVS(X,DX,Y,N,P,COORD, TRACE_COORD)
      IMPLICIT REAL*8(A-H,O-Z)
      integer N
      double precision H,DX(N),X(N),X2(N),X0(N),Y(N),Y2(N), XC(N)
      double precision P(5), YT
      double precision XT
      double precision M
      Character*5 E
      Character*3 COORD(2), FROM, TO, TRACE_COORD
      Dimension Vecin(3),Vecout(3),Vecout2(3)
      Parameter (PI = 3.1415927)
C       TIME = P(1); BY_IMF = P(2); BZ_IMF = P(3); DP = P(4); RLT = P(5)
      E = 'j2000'
      X0 = X !SAVE POSITION IN ORIGINAL COORDS 

      !/////////////////// CONVERT TO CARTESIAN COORDS /////////////////
      IF (TRACE_COORD.EQ.'SPH') CALL SPH2CAR(X0,X)
      !/////////////////////////////////////////////////////////////////

      !///////////// CONVERT CARTESIAN INPUT COORD TO S3C //////////////
      CALL KROT(COORD(2),'S3C',X,Vecout,P(1),E)
      CALL CAR2SPH(Vecout, X)
      !/////////////////////////////////////////////////////////////////

      !CALCULATE MAGNETIC FIELD AT THE GIVEN TIME & POSITION
      CALL KMAG(P(1),E,X(1),X(2),X(3),
     +P(5),P(2),P(3),Y(1),Y(2),Y(3),P(4))

      !CALCULATE MAGNETIC FIELD FOR DIPOPE
C       XT = SQRT(X(1)**2+X(2)**2+X(3)**2)
C       M =  20E-6 !T
C       Y(1) = 2*M*COS(X(2))/X(1)**3
C       Y(2) = M*SIN(X(2))/X(1)**3
C       Y(3) = 0.
      !AN EQUIVALENT EXPRESSION FOR A DIPOPE
C       Y(1) = 3*M*X(1)*X(3)/XT**5
C       Y(2) = 3*M*X(2)*X(3)/XT**5
C       Y(3) = M*(3*X(3)**2-XT**2)/XT**5

      !///////////// CALCULATE INCREMENT AND SPERICAL CORDS ////////////
      IF (TRACE_COORD.EQ.'SPH') THEN
         YT = SQRT(Y(1)**2+Y(2)**2+Y(3)**2);
         DX = (/Y(1)/YT, Y(2)/YT/X0(1), Y(3)/YT/X0(1)/DSIN(X0(2))/)
      END IF
      IF (TRACE_COORD.EQ.'CAR') THEN
         CALL SPH2CAR_MAG(Y(1),Y(2),Y(3),Y2(1),Y2(2),Y2(3),X(2),X(3))
         CALL KROT('S3C',COORD(2),Y2,Y,P(1),E)
         YT = SQRT(Y(1)**2+Y(2)**2+Y(3)**2);
         DX = (/ Y(1)/YT, Y(2)/YT, Y(3)/YT /)
      END IF
      !/////////////////////////////////////////////////////////////////

      X = X0 !RESET X TO ORIGINAL CHOSEN COORDS
      RETURN
      END