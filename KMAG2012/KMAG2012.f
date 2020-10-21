C*********************************************************************************************
C*********************************************************************************************
	SUBROUTINE 	KMAG(TIME,EPOCH,R,THETA,PHI,RLT,BY_IMF,BZ_IMF,
     ~	BRM,BTM,BPM,Dp)
c	R,	THETA and PHI are radial distance, colatitude (in radians) and longitude (in radians)
c	of the spacecraft. The program returns, local time (LT, in hours) and the model 
c	magnetic field in system 3 shperical coordinates.


	IMPLICIT REAL*8(A-H,O-Z)
	Dimension pvecin(3),pvecou(3)
	Dimension bvecin(3),bvecou(3)
	Real*8 time
	character*5 Epoch
	DIMENSION C(5)
	DIMENSION Brscsarr(5),Btscsarr(5),Bpscsarr(5)
	!DATA NL/4/ !number of current sheet modes
	DATA dql /1.41/

      NL=5

C	Strength of tail modes
      ics_on=1
      NL=5
      DATA C /
     ~   2.023453,  
     ~   1.500525, 
     ~   0.441146,  
     ~   0.286629,  
     ~   0.152069/ 
      

C	DATA C/
C     ~  -0.553754,
C     ~   0.714655,
C     ~   0.119748,
C     ~   1.237026/




	DATA Dp0/0.017/ !nominal dynamic pressure


	DATA a1/9.696/
	DATA a2/0.2351/
	DATA dK/0.71516/ !for selfsimilar mp

	DATA dhD/2.0/ !nominal half thickness of the cs


	r0=a1*Dp0**(-a2)


	r0or=r0/(a1*Dp**(-a2))


	D=dhD/r0or !Half thickness of the CS


C
c	FIRST CALCULATE THE INTERNAL FIELD HIGHER THAN DIPOLE 

      CALL  KRONIAN_HIGHER(3,R,THETA,PHI,BRI,BTI,BPI)

C
C	NOW CALCULATE THE AXISYMMETRIC DIPOLE AND ITS SHIELD FIELD
C
      XS3=R*DSIN(THETA)*DCOS(PHI)
      YS3=R*DSIN(THETA)*DSIN(PHI)
      ZS3=R*DCOS(THETA)
C
C	ROTATE TRAJECTORY INTO DIS COORDINATES	
C
	pvecin(1)=xs3
	pvecin(2)=ys3
	pvecin(3)=zs3	
	CALL KROT('S3C','DIS',pvecin,pvecou,time,Epoch)
	XDIS=pvecou(1)
	YDIS=pvecou(2)
	ZDIS=pvecou(3)
c	Calculate the mapped location
	call mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)

	RMAP=dsqrt(XMAP**2+YMAP**2+ZMAP**2)
	

	CALL checkIfInsideMpSATURN(RMAP,XMAP,is_inside_mp,Dp,a1,a2,dK)


	CALL dipole_shielded(time,XMAP,YMAP,ZMAP,BXD,BYD,BZD,r0or)
C
c
c	NOW MAP IT AND ROTATE INTO SYSTEM III
C
	CALL Mapped_field(time,XMAP,YMAP,ZMAP,BXD,BYD,BZD,BXDMAP,BYDMAP,
     ~BZDMAP,r0or,Epoch)
	bvecin(1)=BXDMAP
	bvecin(2)=BYDMAP
	bvecin(3)=BZDMAP
	CALL KROT('DIS','S3C',bvecin,bvecou,time,Epoch)
	BXDS3=bvecou(1)
	BYDS3=bvecou(2)
	BZDS3=bvecou(3)
      BRD=BXDS3*DSIN(THETA)*DCOS(PHI)+BYDS3*DSIN(THETA)*
     .DSIN(PHI)+BZDS3*DCOS(THETA)
      BTD=BXDS3*DCOS(THETA)*DCOS(PHI)+BYDS3*DCOS(THETA)*DSIN(PHI)-
     .BZDS3*DSIN(THETA)
      BPD=-BXDS3*DSIN(PHI)+BYDS3*DCOS(PHI)
	  
C	NOW CALCULATE THE FIELD FROM CURRENT SHEET AND ITS SHIELDING
	call shielded_csheetfield(time,R,THETA,PHI,RLT,Brscs,Btscs,Bpscs,
     ~	NL,C,D,r0or,dql,Epoch,Brscsarr,Btscsarr,Bpscsarr)


C     NOW get the influence of the IMF
	
	CALL getIMF_penetration(BY_IMF,BZ_IMF,BY_p,BZ_p)

	bvecin(1)=0.0
	bvecin(2)=BY_p
	bvecin(3)=BZ_p
	CALL KROT('KSM','S3C',bvecin,bvecou,time,Epoch)
	BX_ps3=bvecou(1)
	BY_ps3=bvecou(2)
	BZ_ps3=bvecou(3)
	
	CALL CAR2SPH_MAG(BX_ps3,BY_ps3,BZ_ps3,BR_p,BT_p,BP_p,THETA,PHI)

	
c	Now add all the fields

	IF (is_inside_mp .eq. 1) THEN

        
        !higher order internal field + shieded dipole + shielded current sheet field + IMF penetration
        IF (ics_on .eq. 1) THEN
		    Brm=BRI+BRD + Brscs + Br_p
		    Btm=BTI+BTD + Btscs + Bt_p
		    Bpm=BPI+BPD + Bpscs + Bp_p
        ELSE
            Brm=BRI+BRD + Br_p
		    Btm=BTI+BTD + Bt_p
		    Bpm=BPI+BPD + Bp_p
        ENDIF

	!Brm=BRD + Brscs + Br_p
	!Btm=BTD + Btscs + Bt_p
	!Bpm=BPD + Bpscs + Bp_p


	!Brm=BRI
	!Btm=BTI
	!Bpm=BPI


!	Brm=Brscs 
!	Btm=Btscs 
!	Bpm=Bpscs 
	ELSE

		Brm=Br_p
		Btm=Bt_p
		Bpm=Bp_p
		
		!set current sheet field to zero if outside mp
		Brscs=0.0
		Btscs=0.0
		Bpscs=0.0
		
		DO i=1,NL
		    Brscsarr(i)=0.0
		    Btscsarr(i)=0.0
		    Bpscsarr(i)=0.0
		ENDDO

	ENDIF

	Return
	End
C*********************************************************************************
C
      SUBROUTINE KRONIAN_HIGHER (NM,R,T,F,BR,BT,BF)
C
C     BASED ON THE SUBROUTINE IGRF WRITTEN BY N.A. TSYGANENKO (1979)
C     MODIFIED BY KRISHAN KHURANA, JULY, 1996. AND   NOV. 2004.
C
C     CALCULATES COMPONENTS OF MAIN KRONIAN FIELD IN RIGHT HANDED SPHERICAL
C     COORD SYSTEM. BASED ON THE  SPHERICAL HARMONIC COEFFICIENTS GIVEN BY 
C     ACUNA ET AL. [1983] (Z3 MODEL,  V1+V2 DATA) 
C     MAXIMUM ORDER OF HARMONICS TAKEN INTO ACCOUNT (NOT MORE THAN 0RDER 3)
C
C
C     IT IS ASSUMED THAT THE TRAJECTORY IS IS IN RIGHT HANDED S III COORDINATES.
C     THE OUTPUT IS ALSO IN RTP (RH) COORDINATES.
C
C            INPUT:  NM (INTEGER)- MAXIMUM ORDER OF HARMONICS TAKEN
C                                  INTO ACCOUNT (NM.LE.12)
C
C                    R,T,F (REAL)- POSITION OF DESIRED FIELD VALUE IN
C                                  SPHERICAL JOVIGRAPHIC COORDINATE SYSTEM
C                                  (R IN PLANET RADII, COLATITUDE T AND
C                                   RH LONGITUDE F IN RADIANS)
C
C           OUTPUT: BR,BT,BF (REAL)- COMPONENTS OF THE INTERNAL PORTION
C                                    OF THE MAIN MAGNETIC FIELD IN
C                                    SPHERICAL S III COORD SYSTEM
C                                    (VALUES GIVEN IN GAMMA)
C
C***************************************************************************
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BK,BM
      DIMENSION A(13),B(13),G(91),H(91),REC(91)
	DATA REC/ 91 * 1.0/


C	DATA G/0.,0.,24.,1449.,-179.,-35.,2012.,
C     ~       53.,92.,-40.,
C     ~       81*0.0/

C	DATA H/0.,0.,-0.9110,0.,43.5680,-44.7762,0.,
C     ~       140.9273,-92.8876,68.2759,
C     ~       81*0.0/

      ! newer coefficients 2012
      DATA G/0.,0.,0.,1586.,0.,0.,2374.,
     ~       0.,0.,0.,65.,0.,0.,0.,0.,185.,
     ~       75*0.0/

      DATA H/91*0.0/
C      1.Cao, H., C.T. Russell, U.R. Christensen, M.K. Dougherty, and M.E. Burton(2011), Saturn’s very axisymmetric magnetic field: No detectable secular variation or tilt, Earth and Planet. Sci. Lett., 304, 22-28 
C      2.Cao, H., Russell, C.T., Christensen, U.R., Wicht, J., Dougherty, M.K., 2012. Saturn’s High Degree Magnetic Moments: Evidence for a Unique Planetary Dynamo. Submitted to Icarus)


      DATA FIRSTI/0.0/
C     WRITE(1,'(5F15.3)')(G(I),I=1,10)
C     WRITE(1,'(5F14.3)')(H(I),I=1,10)
      IF(FIRSTI.EQ.0.0) GO TO 1
      GO TO 6
 1    FIRSTI=1.0
      G(1)=0.
      H(1)=0.
      KNM=15
      DO 2 N=1,13
      N2=2*N-1
      N2=N2*(N2-2)
      DO 2 M=1,N
      MN=N*(N-1)/2+M
  2   REC(MN)=DBLE((N-M)*(N+M-2))/DBLE(N2)
      S=1.
      DO 5 N=2,13
      MN=N*(N-1)/2+1
      S=S*DBLE(2*N-3)/DBLE(N-1)
      G(MN)=G(MN)*S
      H(MN)=H(MN)*S
      P=S
      DO 5 M=2,N
      AA=1.
      IF (M.EQ.2) AA=2.
      P=P*DSQRT(AA*DBLE(N-M+1)/DBLE(N+M-2))
      MNN=MN+M-1
      G(MNN)=G(MNN)*P
  5   H(MNN)=H(MNN)*P
  6   IF(KNM.EQ.NM) GOTO 61
      KNM=NM
      K=KNM+1
 61   PP=1./R
      P=PP
      DO 7 N=1,K
      P=P*PP
      A(N)=P
  7   B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=DCOS(F)
      SF=DSIN(F)
      C=DCOS(U)
      S=DSIN(U)
      BK=(S.LT.1.D-5)
      DO 12 M=1,K
      BM=(M.EQ.1)
      IF(BM) GOTO 8
      MM=M-1
      W=X
      X=W*CF+Y*SF
      Y=Y*CF-W*SF
      GOTO 9
  8   X=0.
      Y=1.
  9   Q=P
      Z=D
      BI=0.
      P2=0.
      D2=0.
      DO 11 N=M,K
      AN=A(N)
      MN=N*(N-1)/2+M
      E=G(MN)
      HH=H(MN)
      W=E*Y+HH*X
      IF (DABS(P2).LT.1.D-38) P2=0.0
      IF (DABS(Q).LT.1.D-38) Q=0.0
      BBR=BBR+B(N)*W*Q
      BBT=BBT-AN*W*Z
      IF(BM) GOTO 10
      QQ=Q
      IF(BK) QQ=Z
      BI=BI+AN*(E*X-HH*Y)*QQ
  10  XK=REC(MN)
      DP=C*Z-S*Q-XK*D2
      PM=C*Q-XK*P2
      D2=Z
      P2=Q
      Z=DP
      Q=PM
   11 CONTINUE
      D=S*D+C*P
      P=S*P
      IF(BM) GOTO 12
      BI=BI*MM
      BBF=BBF+BI
  12  CONTINUE
      BR=BBR
      BT=BBT
      IF(BK) GOTO 13
      BF=BBF/S
      GOTO 14
  13  IF(C.LT.0.) BBF=-BBF
      BF=BBF
  14  CONTINUE
      RETURN
      END
C
C
  
C*************************************************************************

	Subroutine dipole_shielded(time,x,y,z,Bx,By,Bz,r0or)
C           WRITTEN BY K. K. KHURANA     11/2004

C           x,y,z input position
C		  OUTPUT: mag .filed Bx,By,Bz at x,y,z
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION pvecin(3),pvecou(3)

	DATA B0x/0.0/
	DATA B0y/0.0/
	!DATA B0z/21160.0/
	!new dipole strength 2012
C       DATA B0z/26191.0/  !for references see KRONIAN_HIGHER
	DATA B0z/21191.0/  !for references see KRONIAN_HIGHER


c	We will  first calculate the field of the dipole



	call dipole(B0x,B0y,B0z,x,y,z,Bxd,Byd,Bzd)
c
c	Now calculate the perpendicular dipole shielding field
c
	if( (z .eq. 0) .and. (y .eq. 0)) Then
	phi=0.0
	go to 1
	end if
	phi=datan2(z,y)
    1	rho=y*dcos(phi)+z*dsin(phi)



	Call B_mp_perp(rho,phi,x,Brho1,Bphi1,Bx1,8,r0or)


	By1=Brho1*dcos(phi)-Bphi1*dsin(phi)
	Bz1=Brho1*dsin(phi)+Bphi1*dcos(phi)



	Bx=Bxd -Bx1
	By=Byd -By1
	Bz=Bzd -Bz1



	Return
	End
c******************************************************************
c
	Subroutine dipole(B0x,B0y,B0z,x,y,z,Bx,By,Bz)
	IMPLICIT REAL*8(A-H,O-Z)

c	Calculates the field of Jupiter's dipole for shielding in the magnetopause
c	(B0x, B0y, B0z) is the dipole moment
c	x, y, z is the position vector
c     Bx, By, Bz is the output field vector
c
	r=dsqrt(x**2+y**2+z**2)
	a11=(3*x**2-r**2)/r**5
	a12=(3*x*y)/r**5
	a13=(3*x*z)/r**5
	a21=(3*x*y)/r**5
	a22=(3*y**2-r**2)/r**5
	a23=(3*y*z)/r**5
	a31=(3*x*z)/r**5
	a32=(3*y*z)/r**5
	a33=(3*z**2-r**2)/r**5
	Bx=a11*B0x+a12*B0y+a13*B0z
	By=a21*B0x+a22*B0y+a23*B0z
	Bz=a31*B0x+a32*B0y+a33*B0z
	Return
	End
C*************************************************************************
	Subroutine B_mp_perp(rho,phi,x,Bperpr,Bperpf,Bperpx,M,r0or)
	IMPLICIT REAL*8(A-H,O-Z)
	Dimension a(8),b(8)
C

	Data a/3.32125816385087D-004, -0.784083155434928,     
     ~	-0.416208789540178, -3.60127324820496D-002, -0.273188407067209,     
     ~   0.736497845500708, -0.939713175408542,0.433616725727916/     
      Data b/9.24095494266823, 30.4656008995161,     
     ~	   18.1008047927882, 60.6014301282221,  125.347578447821,     
     ~	   253.707136094045, 509.877831023220, 1021.96079781869/ 


	Bperpr=0.
	Bperpf=0.
	Bperpx=0.
	Do 4 K=1,M
	Term1=Dsin(phi)*Dexp(x/(b(K)/r0or))
	Term2=bessj1(rho/(b(K)/r0or))/(rho/(b(K)/r0or))
     ~	-bessj0(rho/(b(K)/r0or))
	Term3=-Dcos(phi)*Dexp(x/(b(K)/r0or))
	Term4=bessj1(rho/(b(K)/r0or))/(rho/(b(K)/r0or))
	Term5=-Dsin(phi)*Dexp(x/(b(K)/r0or))
	Term6=bessj1(rho/(b(K)/r0or))
    	Bperpr=Bperpr+a(K)*r0or**3*Term1*Term2
	Bperpf=Bperpf+a(K)*r0or**3*Term3*Term4
      Bperpx=Bperpx+a(K)*r0or**3*Term5*Term6
    4 continue

      Return
	End
C***************************************************************
      Function DBSJ2(x)
	IMPLICIT REAL*8(A-H,O-Z)
	DBSJ2=(2.0/x)*bessj1(x)-bessj0(x)
	Return
	End
C***********************************************************************

      Function DBSJ3(x)
	IMPLICIT REAL*8(A-H,O-Z)
	DBSJ3=(4.0/x)*DBSJ2(x)-bessj1(x)
	Return
	End




C***********************************************************************

C     Hannes needs his own Bessl functions

C***********************************************************************


C     Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
	function bessj1(x)
      IMPLICIT REAL*8(A-H,O-Z)
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     #242396853.1d0,
     #-2972611.439d0,15704.48260d0,-30.160366606d0/,
     #s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     #18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     #.2457520174d-5,
     #-.240337019d-6/, q1,q2,q3,q4,q5
     #/.04687499995d0,-.2002690873d-3,
     #.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(dabs(x).lt.8) then 
	    y=x**2
	    bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     #      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
	    ax=dabs(x)
	    z=8./ax
	    y=z**2
	    xx=ax-2.356194491
	    bessj1=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*
     #      (p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     #      *dsign(dble(1.),x)
      endif
      return
      end

C***********************************************************************

C     Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
      function bessj0(x)
      IMPLICIT REAL*8(A-H,O-Z)
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     #-.2073370639d-5,.2093887211d-6/
     #,q1,q2,q3,q4,q5/-.1562499995d-1,
     #.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     #651619640.7d0,
     #-11214424.18d0,77392.33017d0,-184.9052456d0/,
     #s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     #9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8) then 
	    y=x**2
	    bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     #      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
	    ax=dabs(x)
	    z=8./ax
	    y=z**2
	    xx=ax-.785398164
	    bessj0=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+
     #      y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end

C**********************************************************************
c
	Subroutine KROT(From,To,Vecin,Vecout,time,epoch)
c
c	INPUTS: From,	a chracter*3 variable denoting the incoming coordinate system
c			To:		a chracter*3 variable denoting the outgoing coordinate system
c			Vecin	a variable of dimension 3 containing the incoming vector
c			Vecout:	a variable of dimension 3 containing the outgoing vector
c			time	A double precision variable denoting number of seconds from an epoch 
c			epoch   A five letter character which can be either 'ctime', 'etime' or 'J2000'/upper or lower case
c
c
C	Rotates the vector vecin from the coordinate system "From" to the
c	system "To". The rotated vector is placed in the vector Vecout.
c
c	We always first rotate the Vecin in "from" system to system III. Then we rotate 
c	it to the "To"coordinate system
c
c	The supported coordinate systems are:
c			S3C System III Cartesian (right-handed)
c			KSO Kronian  -Sun-Orbital
c			KSM Kronian-Sun-Magnetic (So far same as KSO)
c	        DIP Dipole (cartesian)
	Implicit Real*8(A-H,O-Z)
	Dimension Vecin(3),Vecout(3),vector(3)
	Character*3 From, To
	Character*5 epoch
	Real*8 dipole(3,3),first(3,3),second(3,3),dummy(3,3)
c	Real*4 stheta,sphi
	Real*8 time,jtime,eetime
	Parameter (PI=3.1415927,twopi=2.*PI,radian=PI/180.,degree=180./PI)
c      Parameter (TH_DIP=0.0,PH_DIP=0.0) !The dipole parameters
	Data dipole(3,1),dipole(3,2)/0.,0./
	Data dipole(3,3)/1.0/														   
	Data dipole(2,1),dipole(2,2),dipole(2,3)/0.,1.0,0/
	data dipole(1,1),dipole(1,2)/1.0,0./
	data dipole(1,3)/0.0/


!	write(6,*)"dipole: ",dipole

c	First ascertain the time epoch

	if((epoch(1:5).eq."j2000").or.(epoch(1:5).eq."J2000")) then
	jtime=time
	Go to 101
	end if


	if((epoch(1:5).eq."i2000").or.(epoch(1:5).eq."I2000")) then
	jtime=time+32.184D0
	Go to 101
	end if


	  if(((epoch(1:1) .eq. "c") .or. epoch(1:1) .eq. "C") .and.
     ~     ((epoch(2:2) .eq. "t") .or. epoch(2:2) .eq. "T") .and.
     ~     ((epoch(3:3) .eq. "i") .or. epoch(3:3) .eq. "I") .and.
     ~     ((epoch(4:4) .eq. "m") .or. epoch(4:4) .eq. "M") .and.
     ~     ((epoch(5:5) .eq. "e") .or. epoch(5:5) .eq. "E")) then
           jtime=eetime(time)
	      go to 101
	  end if
c
	  if(((epoch(1:1) .eq. "e") .or. epoch(1:1) .eq. "E") .and.
     ~     ((epoch(2:2) .eq. "t") .or. epoch(2:2) .eq. "T") .and.
     ~     ((epoch(3:3) .eq. "i") .or. epoch(3:3) .eq. "I") .and.
     ~     ((epoch(4:4) .eq. "m") .or. epoch(4:4) .eq. "M") .and.
     ~     ((epoch(5:5) .eq. "e") .or. epoch(5:5) .eq. "E")) then
		 jtime=time
	     go to 101
	  end if
	Write(*,*)" I do not understand the time epoch."
	Go to 100
  
c
c	Next ascertain the incoming coordinate system
c	Is it system 3?
  101	If(((From(1:1) .eq. "s") .or. from(1:1) .eq. "S") .and.
     ~   ((From(2:2) .eq. "3") .or. from(2:2) .eq. "3") .and.
     ~   ((From(3:3) .eq. "c") .or. from(3:3) .eq. "C")) Go to 1
c	Is it JSO?
	If(((From(1:1) .eq. "k") .or. from(1:1) .eq. "K") .and.
     ~   ((From(2:2) .eq. "s") .or. from(2:2) .eq. "S") .and.
     ~   ((From(3:3) .eq. "o") .or. from(3:3) .eq. "O")) Go to 2
c	Is it JSM?
	If(((From(1:1) .eq. "k") .or. from(1:1) .eq. "K") .and.
     ~   ((From(2:2) .eq. "s") .or. from(2:2) .eq. "S") .and.
     ~   ((From(3:3) .eq. "m") .or. from(3:3) .eq. "M")) Go to 3
c	Is it Dipole?
	If(((From(1:1) .eq. "d") .or. from(1:1) .eq. "D") .and.
     ~   ((From(2:2) .eq. "i") .or. from(2:2) .eq. "I") .and.
     ~   ((From(3:3) .eq. "p") .or. from(3:3) .eq. "P")) Go to 4
c	Is it Dipole Sun?
	If(((From(1:1) .eq. "d") .or. from(1:1) .eq. "D") .and.
     ~   ((From(2:2) .eq. "i") .or. from(2:2) .eq. "I") .and.
     ~   ((From(3:3) .eq. "s") .or. from(3:3) .eq. "S")) Go to 5

	Write(*,*)" I do not understand the incoming coordinate system"
	Go to 100
    1 continue
c	Write(*,*) " Incoming coordinate system is S3 Cartesian"
	First(1,1)=1.0
	First(1,2)=0.0
	First(1,3)=0.0
	First(2,1)=0.0
	First(2,2)=1.0
	First(2,3)=0.0
	First(3,1)=0.0
	First(3,2)=0.0
	First(3,3)=1.0
	Go to 60
    2 continue
c	Write(*,*) " Incoming coordinate system is KSO"
c	Calculate the rotation matrix to go from KSO to S3R
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSO unit vector etc.
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
C	Calculate the Z axis of the KSO coordinate system.
c
c	Define X component in system III of ZKSO unit vector etc.
	dummy(3,1)=cos(Ztheta)*cos(Zphi)
	dummy(3,2)=cos(Ztheta)*sin(Zphi)
	dummy(3,3)=sin(Ztheta)
c
c	Define X component in system III of YKSO unit vector etc.

	dummy(2,1)=dummy(3,2)*dummy(1,3)-dummy(3,3)*dummy(1,2)
	dummy(2,2)=dummy(3,3)*dummy(1,1)-dummy(3,1)*dummy(1,3)
	dummy(2,3)=dummy(3,1)*dummy(1,2)-dummy(3,2)*dummy(1,1)
c	The transpose of the dummy matrix takes KSO to System 3.
	Do 21 i=1,3
	Do 21 j=1,3
   21	first(i,j)=dummy(j,i)
	go to 60

    3 continue
c	Write(*,*) " Incoming coordinate system is KSM"

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)

c	Now define the KSM transpose vector
c	Define X component in system III of XJSM unit vector etc.

	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c	Now define the Y vector so that it is perpendicular to the dipole vector and X 
	dummy(2,1)=dummy(1,3)*dipole(3,2)-dummy(1,2)*dipole(3,3)
	dummy(2,2)=dummy(1,1)*dipole(3,3)-dummy(1,3)*dipole(3,1)
	dummy(2,3)=dummy(1,2)*dipole(3,1)-dummy(1,1)*dipole(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom
c	Now define the z vector
	dummy(3,1)=dummy(1,2)*dummy(2,3)-dummy(1,3)*dummy(2,2)
	dummy(3,2)=dummy(1,3)*dummy(2,1)-dummy(1,1)*dummy(2,3)
	dummy(3,3)=dummy(1,1)*dummy(2,2)-dummy(1,2)*dummy(2,1)
	Do 30 i=1,3
	Do 30 j=1,3
   30 first(i,j)=dummy(j,i)
	Go to 60
c

    4 continue
c	Write(*,*) " Incoming coordinate system is DIP"
	Do 40 i=1,3
	Do 40 j=1,3
   40 first(i,j)=dipole(j,i)
	go to 60
c
    5 continue
c	Write(*,*) " Incoming coordinate system is DIS"

c	Define  ZDIS unit vector 
	dummy(3,1)=	dipole(3,1)
	dummy(3,2)=	dipole(3,2)
	dummy(3,3)=	dipole(3,3)

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)


c	Define  XKSO unit vector 
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c
C	Define YDIS unit vector	(as normal to ZDIS and XKSO)
	dummy(2,1)=dummy(1,3)*dummy(3,2)-dummy(1,2)*dummy(3,3)
	dummy(2,2)=dummy(1,1)*dummy(3,3)-dummy(1,3)*dummy(3,1)
	dummy(2,3)=dummy(1,2)*dummy(3,1)-dummy(1,1)*dummy(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom

c	Define XDIS unit vector as normal to YDIS and ZDIS
c
	dummy(1,1)=dummy(2,2)*dummy(3,3)-dummy(2,3)*dummy(3,2)
	dummy(1,2)=dummy(2,3)*dummy(3,1)-dummy(2,1)*dummy(3,3)
	dummy(1,3)=dummy(2,1)*dummy(3,2)-dummy(2,2)*dummy(3,1)


c	The transpose of the dummy matrix takes DIS to System 3.
	Do 22 i=1,3
	Do 22 j=1,3
   22	first(i,j)=dummy(j,i)


c
c
   60 continue
c
c	Now Figure out the outgoing coordinate system

	If(((To(1:1) .eq. "s") .or. To(1:1) .eq. "S") .and.
     ~   ((To(2:2) .eq. "3") .or. To(2:2) .eq. "3") .and.
     ~   ((To(3:3) .eq. "c") .or. To(3:3) .eq. "C")) Go to 6
c	Is it JSO?
	If(((To(1:1) .eq. "k") .or. To(1:1) .eq. "K") .and.
     ~   ((To(2:2) .eq. "s") .or. To(2:2) .eq. "S") .and.
     ~   ((To(3:3) .eq. "o") .or. To(3:3) .eq. "O")) Go to 7
c	Is it JSM?
	If(((To(1:1) .eq. "k") .or. To(1:1) .eq. "K") .and.
     ~   ((To(2:2) .eq. "s") .or. To(2:2) .eq. "S") .and.
     ~   ((To(3:3) .eq. "m") .or. To(3:3) .eq. "M")) Go to 8
c	Is it Dipole?
	If(((To(1:1) .eq. "d") .or. To(1:1) .eq. "D") .and.
     ~   ((To(2:2) .eq. "i") .or. To(2:2) .eq. "I") .and.
     ~   ((To(3:3) .eq. "p") .or. To(3:3) .eq. "P")) Go to 9
c	Is it DIS?
	If(((To(1:1) .eq. "d") .or. To(1:1) .eq. "D") .and.
     ~   ((To(2:2) .eq. "i") .or. To(2:2) .eq. "I") .and.
     ~   ((To(3:3) .eq. "s") .or. To(3:3) .eq. "S")) Go to 10
c
	Write(*,*)" I do not understand the outgoing coordinate system"
	Go to 100

    6	Continue
c	Write(*,*) " The outgoing coordinate system is S3 cartesian"
	Second(1,1)=1.0
	Second(1,2)=0.0
	Second(1,3)=0.0
	Second(2,1)=0.0
	Second(2,2)=1.0
	Second(2,3)=0.0
	Second(3,1)=0.0
	Second(3,2)=0.0
	Second(3,3)=1.0
	Go to 70
c
    7 continue
c	Write(*,*) " The outgoing coordinate system is KSO"
c	Get the matrix that Rotates from System III cartesian into KSO.	
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSO unit vector etc.
	Second(1,1)=cos(stheta)*cos(sphi)
	Second(1,2)=cos(stheta)*sin(sphi)
	Second(1,3)=sin(stheta)
c
c	Define X component in system III of ZKSO unit vector etc.
	Second(3,1)=cos(Ztheta)*cos(Zphi)
	Second(3,2)=cos(Ztheta)*sin(Zphi)
	Second(3,3)=sin(Ztheta)
c
c	Check orthogonality of the new coordinate system
	dotprod=second(3,1)*second(1,1)+second(3,2)*second(1,2)+
     .second(3,3)*second(1,3)
C	Write(*,*) dotprod
c	Define X component in system III of YKSO unit vector etc.

	Second(2,1)=Second(3,2)*Second(1,3)-Second(3,3)*Second(1,2)
	Second(2,2)=Second(3,3)*Second(1,1)-Second(3,1)*Second(1,3)
	Second(2,3)=Second(3,1)*Second(1,2)-Second(3,2)*Second(1,1)
c	seconm=	sqrt(second(2,1)**2+second(2,2)**2+second(2,3)**2)
c	Write(*,*)seconm
	Go to 70
    8 continue
c	Write(*,*) " The outgoing coordinate system is KSM"

c	Get the matrix that Rotates from System III cartesian into KSM.	
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSM unit vector etc.
	Second(1,1)=cos(stheta)*cos(sphi)
	Second(1,2)=cos(stheta)*sin(sphi)
	Second(1,3)=sin(stheta)
c	Now define the Y vector so that it is perpendicular to  the dipole vector and X 
	Second(2,1)=second(1,3)*dipole(3,2)-second(1,2)*dipole(3,3)
	Second(2,2)=second(1,1)*dipole(3,3)-second(1,3)*dipole(3,1)
	Second(2,3)=second(1,2)*dipole(3,1)-second(1,1)*dipole(3,2)
	denom=sqrt(second(2,1)**2+second(2,2)**2+second(2,3)**2)
	Second(2,1)=second(2,1)/denom
	Second(2,2)=second(2,2)/denom
	Second(2,3)=second(2,3)/denom
c	Now define the z vector
	Second(3,1)=second(1,2)*second(2,3)-second(1,3)*second(2,2)
	Second(3,2)=second(1,3)*second(2,1)-second(1,1)*second(2,3)
	Second(3,3)=second(1,1)*second(2,2)-second(1,2)*second(2,1)
	Go to 70
c
    9 continue
c	Write(*,*) " The outgoing coordinate system is Dipole"
	Do 80 i=1,3
	Do 80 j=1,3
   80	Second(i,j)=dipole(i,j)  
	go to 70
   10 Continue
c	Write(*,*) " Outgoing coordinate system is DIS"

c	Define  ZDIS unit vector 
	dummy(3,1)=	dipole(3,1)
	dummy(3,2)=	dipole(3,2)
	dummy(3,3)=	dipole(3,3)


	!write(6,*)"DUMMY1: ",dummy

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define  XKSO unit vector 
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c

	!write(6,*)"stheta: ",stheta
	!write(6,*)"sphi: ",sphi

	!write(6,*)"DUMMY2: ",dummy

C	Define YDIS unit vector	(as normal to ZDIS and XKSO)
	dummy(2,1)=dummy(1,3)*dummy(3,2)-dummy(1,2)*dummy(3,3)
	dummy(2,2)=dummy(1,1)*dummy(3,3)-dummy(1,3)*dummy(3,1)
	dummy(2,3)=dummy(1,2)*dummy(3,1)-dummy(1,1)*dummy(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom

      !write(6,*)"DUMMY3: ",dummy


c	Define XDIS unit vector as normal to YDIS and ZDIS
c
	dummy(1,1)=dummy(2,2)*dummy(3,3)-dummy(2,3)*dummy(3,2)
	dummy(1,2)=dummy(2,3)*dummy(3,1)-dummy(2,1)*dummy(3,3)
	dummy(1,3)=dummy(2,1)*dummy(3,2)-dummy(2,2)*dummy(3,1)


c	The  dummy matrix takes  System 3 to DIS.
	Do 23 i=1,3
	Do 23 j=1,3
   23	second(i,j)=dummy(i,j)
c
   70	continue
c	Now multimply vecin with first and second matrices to get the vecout
c	Write(*,*)" calling matmult 1"
	call matmult(first,vecin,vector,3,3,1)
c	Write(*,*)" calling matmult 2"
	call matmult(second,vector,vecout,3,3,1)

	!write(6,*)"dipole: ",dipole
	!write(6,*)"DUMMY: ",dummy
	!write(6,*)"FIRST: ",first
	!write(6,*)"sECOND: ",second



c
  100 Continue	
	Return
	End
c***********************************************************************

c***********************************************************************
	Subroutine KSun(time,stheta,sphi,Ztheta,Zphi)
c	INPUT:  J2000 time of the data point 
C	OUTPUTS: stheta, sphi, latitude and longitude  (in radians) of the Sun in system III (RH).  
C	OUTPUTS: Ztheta, Zphi, latitude and longitude  (in radians) of the Zkso in system III (RH).  
c	The equations are written in J2000 convention followed by PDS. 
c

c	Last updated August 26, 2004.
c
c	theta=a1*cos(omegay*t)+a2*sin(omegay*t)+
c	a3*cos(2.*omegay*t)+a4*sin(2.*omegay*t)+
c	a5*cos(3.*omegay*t)+a6*sin(3.*omegay*t)+
c	a7*cos(4.*omegay*t)+a8*sin(4.*omegay*t)+
c     a9*t**2+a10*t+a11
c	fphi=b1*cos(omegay*t)+b2*sin(omegay*t)+
c     b3*cos(2.*omegay*t)+b4*sin(2.*omegay*t)+
c     b5*cos(3.*omegay*t)+b6*sin(3.*omegay*t)+
c     b7*cos(4.*omegay*t)+b8*sin(4.*omegay*t)+
c	b9*t**2+b10*t+b11	(I assume I have despun Saturn first.)

c	Then we rotate into the System III coordinates

c	time is a double precision variable and is assumed to be in J2000. theta and phi are single precision 
c	variables.
c

c	
c	omega is Saturn's rotation rate
c	omegaz is saturn's rotation rate for use of Zkso lat and lon
c	omegay is Saturn's yearly orbital rate

c	Initialize variables
	Implicit Real*8(A-H,O-Z)
	Real*8 aa(11),bb(11), cc(11),x(11)
	Real*8 time,jtime1,fphi,yrsat,omega,D360,omegay,t,year,omegaz
c	Real*4 stheta,sphi,ztheta,Zphi
	Parameter (PI=3.141592654d0)
	PARAMETER (twopi=2.*PI,radian=PI/180.,degree=180./PI)
	Parameter (yrsat=10759.22d0*86400.0d0)
      Parameter (omega=360.D0/(10.65622222D0*3600.D0))
C       Parameter (omega=360.D0/(10.53333333D0*3600.D0)) ![Helled,2015]
	PARAMETER (omegay=2.*PI/yrsat)
      Parameter (omegaz=360.D0/(10.65622222D0*3600.D0)+3.880688963D-7)
C 	Parameter (omegaz=360.D0/(10.53333333D0*3600.D0)+3.880688963D-7)
	Parameter (year=86400D0*.36525D3,D360=360.D0)
	Parameter (jtime1=-.825767955816D+09)
	Parameter (Zthetd=63.271d0)
c	Parameter (Zthetd=-26.729)
	Data aa  /-.26237934D+02,.30525104D+01,-.13686733D+01,
     +.10182726D+00,-.30897805D+00,.89277949D-01,-.39638771D-01,
     +.97499653D-02,.40974159D-04,.14075684D-02,.13505775D+01/
c
	Data bb /-.32108124D+00,.58569866D+01,.71266272D+00,.33244539D+01,
     +.47529267D-01,.32770362D+00,.46935622D-01,.10720536D+00,
     + -.50594764D-03,.29439612D-01,.26423581D+03/
c
	Data cc/-.13200649D-02,-.18358429D-02,.64658927D-03,.12800487D-02,
     +.17618936D-04,-.58790898D-03,.49804081D-04,.42372137D-03,
     +.14744891D-04,.25369744D-01,.77943328D+02/

c	First calculate the latitude and longitude in non-rotating Saturn coordinates.

c	Calculate the best fit theta and fphi
	t=time-jtime1
	x(1)=cos(omegay*t)
	x(2)=sin(omegay*t)
	x(3)=cos(2.*omegay*t)
	x(4)=sin(2.*omegay*t)
	x(5)=cos(3.*omegay*t)
	x(6)=sin(3.*omegay*t)
	x(7)=cos(4.*omegay*t)
	x(8)=sin(4.*omegay*t)
	x(9)=(t/year)**2
	x(10)=t/year
	x(11)=1.0
	stheta=0.
c	fphi is phi in Saturn fixed (non-rotating) coordinate
	fphi=0.0
	Zfphi=0.0
	Do 1 j=1,11
	fphi=fphi+bb(j)*x(j)
	Zfphi=Zfphi+cc(j)*x(j)
    1 stheta=stheta+aa(j)*x(j)
C	Now rotate the longitude to Saturn System III
c	First Add the rotation of Saturn around the Sun.
c	fphi is the phi of the Sun as unspinning Saturn goes around the sun
	fphi=DMod(fphi+t/yrsat*360.d0, D360)
c	Next add the rotation of Saturn around its axis.
	sphi=DMod(fphi-t*omega, D360)
	if (sphi .lt. 0.0) sphi = sphi+360.0
	sphi=sphi*radian
	stheta=stheta*radian
c
C	Next rotate the longitude of Zfphi to Saturn System III
c	First Add the rotation of Saturn around the Sun.
c	Zfphi is the phi of Zkso as unspinning Saturn goes around the sun
	 Zfphi=DMod(Zfphi+t/yrsat*360.d0+180.D0, D360)
c	Next add the rotation of Saturn around its axis.
	Zphi=DMod(Zfphi-t*omegaz, D360)
	if (Zphi .lt. 0.0) Zphi = Zphi+360.0
	Zphi=Zphi*radian
	Ztheta=Zthetd*radian
	Return
      end 
c*******************************************************************
	Function eetime(ctime)
	IMPLICIT real*8(A-H,O-z)
	Real*8 ctime, eetime
	Tcor=0.
	if(ctime .ge. 189302400.000) Tcor = Tcor+10	
	if(ctime .ge. 205027200.000) Tcor = Tcor+1	
	if(ctime .ge. 220924800.000) Tcor = Tcor+1
	if(ctime .ge. 252460800.000) Tcor = Tcor+1
	if(ctime .ge. 283996800.000) Tcor = Tcor+1
	if(ctime .ge. 315532800.000) Tcor = Tcor+1
	if(ctime .ge. 347155200.000) Tcor = Tcor+1
	if(ctime .ge. 378691200.000) Tcor = Tcor+1
	if(ctime .ge. 410227200.000) Tcor = Tcor+1
	if(ctime .ge. 441763200.000) Tcor = Tcor+1
	if(ctime .ge. 489024000.000) Tcor = Tcor+1
	if(ctime .ge. 520560000.000) Tcor = Tcor+1
	if(ctime .ge. 552096000.000) Tcor = Tcor+1
	if(ctime .ge. 615254400.000) Tcor = Tcor+1
	if(ctime .ge. 694224000.000) Tcor = Tcor+1
	if(ctime .ge. 757382400.000) Tcor = Tcor+1
	if(ctime .ge. 788918400.000) Tcor = Tcor+1
	if(ctime .ge. 836179200.000) Tcor = Tcor+1
	if(ctime .ge. 867715200.000) Tcor = Tcor+1
	if(ctime .ge. 899251200.000) Tcor = Tcor+1
	if(ctime .ge. 946684800.000) Tcor = Tcor+1
	if(ctime .ge. 993945600.000) Tcor = Tcor+1
	if(ctime .ge. 1041379200.000) Tcor = Tcor+1
	eetime=ctime+dble(Tcor)-.1072958367816D10
	Return
	End
C*************************************************************************


      SUBROUTINE matmult(xx,yy,zz,nxrow,nxcol,nycol)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xx(nxrow,nxcol),yy(nxcol,nycol),zz(nxrow,nycol) 
	Do 1 i=1,nxrow
	Do 1 j=1,nycol
	zz(i,j)=0.
	Do 1 k=1,nxcol
	zz(i,j)=zz(i,j)+xx(i,k)*yy(k,j)
    1 continue
      RETURN
      END 
C******************************************************************
	Function ctimer(iyr,imon,iday,ihr,imin,sec)
	integer iyr,imon,iday,ihr,imin
	integer*4 ndays,doy
	Real*8 sec,ctimer

c	First calculate the number of days from Jan 1, 1966
		ndays=0
		If	(iyr .ge. 1966) then
		Do 1 i=1966,iyr-1
		ndays=ndays+365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or. 
     +	(mod(i,400)) .eq. 0) ndays=ndays+1
   1		continue
c	Now add the number of days of the current year
    	ndays=ndays+doy(iyr,imon,iday)-1	
    	ctimer = dble(ndays)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
     +dble(sec)
	go to 4
	end if

c	Calculate the seconds for the negative years
	If	(iyr .lt. 1966) then
	Do 2 i=iyr, 1965
	ndays=ndays-365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or. 
     +	(mod(i,400)) .eq. 0) ndays=ndays-1
    2 continue
c	Now subtract the number of days of the current year	
	ndays=ndays+doy(iyr,imon,iday)
	end if
    	ctimer = dble(ndays-1)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
     +dble(sec)
    4 continue
	idyr=doy(iyr,imon,iday)
	Write(*,*) " doy is ",idyr
	Write(*,*)" ndays is ",ndays
	Return 
	End
c*******************************************************************
	Function doy(iyr,imon,iday)
	Integer mon(12),doy
	data mon /31,28,31,30,31,30,31,31,30,31,30,31/
	doy=0
	Do 1 i = 2,imon
    	doy=doy+mon(i-1)
c	Add an extra day for February
	if(i .eq. 3) then
		if(((mod(iyr,4) .eq. 0) .and. (mod(iyr,100) .ne. 0)) .or. 
     +	(mod(iyr,400)) .eq. 0) doy=doy+1
	end if

    1 continue
	doy=doy+iday
	Return
	End
C********************************************************************************
C********************************************************************************
	Subroutine mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)
	IMPLICIT REAL*8(A-H,O-Z)
	data RH0/23.0/
	
	Dimension pvecin(3),pvecou(3)
	character*5 Epoch
C	r=sqrt(XKSM**2+YKSM**2+ZKSM**2)

C	r0or=1.0

      !write(*,*)r0or

	RH=RH0/r0or !hinging as function of pressure
c	First calculate the dipole tilt (si) in the KSM coordinates
	pvecin(1)=0.
	pvecin(2)=0.
	pvecin(3)=1.
	CALL KROT('S3C','KSM',pvecin,pvecou,time,epoch)
	si=	dacos(pvecou(3))
	if(pvecou(1) .lt. 0) si=-si
c	Now calculate the height of the current sheet
	rhoDIS=dsqrt(XDIS**2+YDIS**2)
	phiLT=datan2(YDIS,XDIS)
	theta=datan2(ZDIS,rhoDIS)
	Zcs=(rhoDIS-RH*dtanh(rhoDIS/RH))*dtan(-si) !sign of si?
	Thcs=datan2(Zcs,rhoDIS) !sign of Thcs?
	
	if (XDIS .gt. 0.0) Thcs=-Thcs ! change 11/09/2012
	
	rhomap= rhoDIS*dcos(thcs)+ZDIS*dsin(thcs) 
	Zmap=  -rhoDIS*dsin(thcs)+ZDIS*dcos(thcs)
	Xmap=rhomap*dcos(phiLT)
	Ymap=rhomap*dsin(phiLT)
	Return
	End
C********************************************************************************
C********************************************************************************
	Subroutine mapitMP(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)
	IMPLICIT REAL*8(A-H,O-Z)
	data RH0/23.0/
	Dimension pvecin(3),pvecou(3)
	character*5 Epoch
C	r=sqrt(XKSM**2+YKSM**2+ZKSM**2)

C	r0or=1.0

      !write(*,*)r0or

	RH=RH0/r0or !hinging as function of pressure
	!RH=0.0
	
c	First calculate the dipole tilt (si) in the KSM coordinates
	pvecin(1)=0.
	pvecin(2)=0.
	pvecin(3)=1.
	CALL KROT('S3C','KSM',pvecin,pvecou,time,epoch)
	si=	dacos(pvecou(3))
	if(pvecou(1) .lt. 0) si=-si
c	Now calculate the height of the current sheet
	rhoDIS=dsqrt(XDIS**2+YDIS**2)
	phiLT=datan2(YDIS,XDIS)
	theta=datan2(ZDIS,rhoDIS)
	Zcs=(rhoDIS-RH*dtanh(rhoDIS/RH))*dtan(+si) !sign of si?
	Thcs=datan2(Zcs,rhoDIS) !sign of Thcs?
	
	if (XDIS .gt. 0.0) Thcs=-Thcs ! change 11/09/2012
	
	
	rhomap= rhoDIS*dcos(thcs)+ZDIS*dsin(thcs) 
	Zmap=  -rhoDIS*dsin(thcs)+ZDIS*dcos(thcs)
	
	
	
	Xmap=rhomap*dcos(phiLT)
	Ymap=rhomap*dsin(phiLT)
	Return
	End
C********************************************************************************
C********************************************************************************
	SUBROUTINE Mapped_field(
     ~	time,X,Y,Z,BX,BY,BZ,BXMAP,BYMAP,BZMAP,r0or,Epoch)
	IMPLICIT REAL*8(A-H,O-Z)
	character*5 Epoch
C	DIMENSION pvecin(3),pvecou(3)
c
C	r0or=1.0
c
c
c	Now define the nine derivatives
c	These are calculated at the original location
c
	dx=0.01
	dy=0.01
	dz=0.01
	xp=x+dx
	xm=x-dx
	yp=y+dy
	ym=y-dy
	zp=z+dz
	zm=z-dz
c	
c	We begin with the x derivatives
	
	call mapit(time,xp,y,z,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,xm,y,z,xpm,ypm,zpm,r0or,Epoch)

	dxpdx=(xpp-xpm)/(2.*dx)
	dypdx=(ypp-ypm)/(2.*dx)
	dzpdx=(zpp-zpm)/(2.*dx)

c	Next the y derivatives
c

C	dxpdy=0.
C	dypdy=1.
C	dzpdy=0.


	call mapit(time,x,yp,z,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,x,ym,z,xpm,ypm,zpm,r0or,Epoch)


	dxpdy=(xpp-xpm)/(2.*dy)
	dypdy=(ypp-ypm)/(2.*dy)
	dzpdy=(zpp-zpm)/(2.*dy)


c	Next the z gradients 

	call mapit(time,x,y,zp,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,x,y,zm,xpm,ypm,zpm,r0or,Epoch)

	dxpdz=(xpp-xpm)/(2.*dz)
	dypdz=(ypp-ypm)/(2.*dz)
	dzpdz=(zpp-zpm)/(2.*dz)


c
C	Now calculate the T matrix
c  	Calculate the mapped location
c
	Txx=dypdy*dzpdz-dypdz*dzpdy
	Txy=dxpdz*dzpdy-dxpdy*dzpdz
	Txz=dxpdy*dypdz-dxpdz*dypdy
	Tyx=dypdz*dzpdx-dypdx*dzpdz
	Tyy=dxpdx*dzpdz-dxpdz*dzpdx
	Tyz=dxpdz*dypdx-dxpdx*dypdz
	Tzx=dypdx*dzpdy-dypdy*dzpdx
	Tzy=dxpdy*dzpdx-dxpdx*dzpdy
	Tzz=dxpdx*dypdy-dxpdy*dypdx
		
C	Now calculate the field at the mapped location 
c
      Bxmap=Txx*Bx+Txy*By+Txz*Bz	
      Bymap=Tyx*Bx+Tyy*By+Tyz*Bz	
      Bzmap=Tzx*Bx+Tzy*By+Tzz*Bz	
	return 
	End

C******************************************************************
	Subroutine shielded_csheetfield(time,R1,T1,P1,RLT,Brm1,Btm1,
     ~Bpm1,NL,C,D,r0or,dql,epoch,Brm1arr,Btm1arr,Bpm1arr)
      IMPLICIT REAL*8(A-H,O-Z)
c	This subroutine calculates the best fit current sheet by modeling only
c	Brho and Bz components. Once the coefficients for the 7 modes are known, 
c	the other subroutines MAPIT and Mapped_field introduce the stretch transformations 
c	and calculate the new values of Brho,Bphi and Bz. 
c	Real*4 stheta,sphi,ztheta,Zphi

	Dimension pvecin(3),pvecou(3)
	Dimension bvecin(3),bvecou(3)
	Real*8 time
	character*5 Epoch
	Dimension f(6,6),beta(6,6),C(NL)
	Dimension Xrho(NL),Xz(NL),Brm1arr(NL),Btm1arr(NL),Bpm1arr(NL)


C     RING CURRENT

C     ============

C     Peak AT 12

c	f(1,1)=    -68706.3686478971

c	f(1,2)=       -2.9523400227

c	f(1,3)=   100047.0971381391

c	f(1,4)=   -30783.7506416514

c	f(1,5)=       -0.0000642075

c	f(1,6)=       -0.0031683487



c	beta(1,1)=     9.2005569531

c	beta(1,2)=     2.1784062411

c	beta(1,3)=     8.8767954666

c     beta(1,4)=     8.2680804765

c	beta(1,5)=     0.1100281701

c	beta(1,6)=     0.4863609871



c     PEAK AT 8

C	f(1,1)=    99323.9271557975

C	f(1,2)=       22.6189468473

C	f(1,3)=      720.7929798854

C	f(1,4)=     -485.4538206945

C	f(1,5)=   -37918.4097869931

C	f(1,6)=   -61520.1671330672



c bi

C	beta(1,1)=   5.4851103313

C	beta(1,2)=   1.9438277867

C	beta(1,3)=   2.9025159249

C      beta(1,4)=   2.6449016349

C	beta(1,5)=   5.1921884427

C	beta(1,6)=   5.6729005214



c     PEAK AT 6

	!f(1,1)=	   -96480.9970060387
	!f(1,2)=		 -733.2997300981
	!f(1,3)=	   -17474.3949850542
	!f(1,4)=		13878.3552501633
	!f(1,5)=		99937.8587320000
	!f(1,6)=		  988.4579332559

	!beta(1,1)=		4.5216256036
	!beta(1,2)=		2.5132995869
	!beta(1,3)=		3.2135379771
	!beta(1,4)=		3.0509806322
	!beta(1,5)=		4.4629794036
	!beta(1,6)=		6.6822375235

c     PEAK AT 4

      !f(1,1)=	   69605.0534869644
      !f(1,2)=	   149.0772338679
      !f(1,3)=	   78909.3792557924
      !f(1,4)=	-54703.2221035385
      !f(1,5)=	-93857.8599750019
      !f(1,6)=	   -59.8779082335

C bi
      !beta(1,1)=    3.7051342090
      !beta(1,2)=    2.3672348371
      !beta(1,3)=    3.2236662867
      !beta(1,4)=    3.1391528825
      !beta(1,5)=    3.6275112423
      !beta(1,6)=    7.8003148590


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! fi: new mode Aug 2012 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !f(1,1)=   828.7022879010
      !f(1,2)= -2953.8174466019
      !f(1,3)=   587.3538467185
      !f(1,4)=  -237.6901233573
      !f(1,5)=  -303.3940846623
      !f(1,6)=  2500.5187349814
! bi
      !beta(1,1)=  59.2734758066
      !beta(1,2)=  34.4268323961
      !beta(1,3)=   4.0855094297
      !beta(1,4)=   3.1278017056
      !beta(1,5)=   8.7602255391
      !beta(1,6)=  30.9667525321
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! fi: new mode Aug 2012 start at 1.2 peak at 5
      !f(1,1)=     559.6197083702
      !f(1,2)=   -1968.8800437331
      !f(1,3)=     584.4840512599
      !f(1,4)=    -364.3594622078
      !f(1,5)=    -166.6634660654
      !f(1,6)=    1659.0247186932
!bi
       !beta(1,1)=  59.2711630145
       !beta(1,2)=  34.4177680800
       !beta(1,3)=   3.1524952826
       !beta(1,4)=   2.7593883598
       !beta(1,5)=   8.4957542421
       !beta(1,6)=  30.9550941820
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! fi: new mode Aug 2012 start at 2 peak at 5
      !f(1,1)=       9.7672928003
      !f(1,2)=      63.5144165850
      !f(1,3)=   -2975.8264925754
      !f(1,4)=    3960.8018682005
      !f(1,5)=     -65.1828634143
      !f(1,6)=    -863.4772958890
!bi
      !beta(1,1)=    1.5894275795
      !beta(1,2)=   34.3759030779
      !beta(1,3)=    3.6428139222
      !beta(1,4)=    3.4574753261
      !beta(1,5)=    8.2962286103
      !beta(1,6)=    2.8921770952
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! fi: new mode Aug 2012 start at 5 peak at 8
 !fi
!! fi: new mode Aug 2012 start at 5 peak at 8
 !fi
C      f(1,1)=   58371.7584883787
C      f(1,2)=      67.2696628374
C      f(1,3)=    4829.5415810648
C      f(1,4)=   -3247.0372511400
C      f(1,5)=  -59240.7786762627
C      f(1,6)=    -748.9828504765
 !bi
C       beta(1,1)=   4.3882059995
C       beta(1,2)=   2.5209029894
C       beta(1,3)=   3.2052095595
C       beta(1,4)=   3.0487351464
C       beta(1,5)=   4.3379025266
C       beta(1,6)=   8.2846806078

C 
 !fi
       f(1,1)=     -100.5052977296
       f(1,2)=     5366.2712701130
       f(1,3)=   -15885.9490677644
       f(1,4)=    40732.9959830497
       f(1,5)=   -32129.0816517276
       f(1,6)=     2215.9395715640
 !bi
       beta(1,1)=     2.7117925416
       beta(1,2)=     4.0941642385
       beta(1,3)=     4.6702786290
       beta(1,4)=     5.7582130292
       beta(1,5)=     6.2062614171
       beta(1,6)=     8.4355483099
       
C
       f(2,1)=   -78718.3907963260
       f(2,2)=      -54.8132316673
       f(2,3)=    -2414.9302620727
       f(2,4)=     1616.7325480546
       f(2,5)=    74549.4133726164
       f(2,6)=     5883.6406782211
C
       beta(2,1)=    5.0017041440
       beta(2,2)=    2.6863215815
       beta(2,3)=    3.5010631982
       beta(2,4)=    3.3102786530
       beta(2,5)=    4.9394705994
       beta(2,6)=    6.1620526042



        f(3,1)=     828.2831671297
        f(3,2)=  -61442.3855659502
        f(3,3)=  -22900.6452859245
        f(3,4)=   40528.1409725691
        f(3,5)=  -42635.0448098891
        f(3,6)=   87296.6633326159
     
        beta(3,1)=    72.3455716311
        beta(3,2)=    10.1486476486
        beta(3,3)=     8.0002370990
        beta(3,4)=     8.4224589108
        beta(3,5)=    12.8834880849
        beta(3,6)=    11.6138431830

        f(4,1)=        -52.9340216641
        f(4,2)=          4.9006893928
        f(4,3)=        103.5467051654
        f(4,4)=     -16743.1838794522
        f(4,5)=      26132.5659128414
        f(4,6)=      -7214.7617456977
 
        beta(4,1)=      4.7199163645
        beta(4,2)=      2.5733685428
        beta(4,3)=      9.5232221496
        beta(4,4)=     18.2512781842
        beta(4,5)=     20.4896316537
        beta(4,6)=     59.5697679599


        !f(5,1)=             9.6382346663
        !f(5,2)=            -1.9403694819
        !f(5,3)=          -110.7011194690
        !f(5,4)=          -212.0255553191
        !f(5,5)=        -15859.9202979478
        !f(5,6)=         55821.0012327883
     
        !beta(5,1)=          4.5856791348
        !beta(5,2)=          2.5276198276
        !beta(5,3)=          8.3627067607
        !beta(5,4)=         17.9703567386
        !beta(5,5)=         36.2969029093
        !beta(5,6)=         46.4652858681

!fi
        f(5,1)=     -23.3052813129
        f(5,2)=       2.2671216463
        f(5,3)=      -4.0457878393
        f(5,4)=   -2345.2488224382
        f(5,5)=   89556.6128080047
        f(5,6)=  -70560.3264679953
 !bi
        beta(5,1)=    4.6355942963
        beta(5,2)=    2.5383712259
        beta(5,3)=    8.4531654549
        beta(5,4)=   20.0059847755
        beta(5,5)=   42.8693646158
        beta(5,6)=   47.3861513994


       

C L=4
      ! taken: use old mode again 
C      f(2,1)=     -60.6568564967
C      f(2,2)=  -73849.9534670401
C      f(2,3)=   97442.8437444251
C      f(2,4)=  -23178.2043425704
C      f(2,5)=      66.3747390965
C      f(2,6)=      -1.4848752731

C bi
C      beta(2,1)=   48.1885954914
C      beta(2,2)=    5.7107159217
C      beta(2,3)=    5.6345236716
C      beta(2,4)=    5.3879744819
C      beta(2,5)=   89.0583961619
C      beta(2,6)=   22.6415533854
      
      
C L=4
C new mode Aug 2012 start at 2 peak at 5
      !f(2,1)=      1.0753834757
      !f(2,2)=    -87.8886394188
      !f(2,3)=    213.0301443067
      !f(2,4)=    180.5413030125
      !f(2,5)=    114.7177473372
      !f(2,6)=   -196.5847020514

C bi
      !beta(2,1)=     1.5383102306
      !beta(2,2)=    34.3790120138
      !beta(2,3)=     3.4830213559
      !beta(2,4)=     3.2635807863
      !beta(2,5)=     8.3908399544
      !beta(2,6)=     2.6750562404     


C L=8            
      ! taken: use old mode again
C      f(3,1)=      310.2640136545
C      f(3,2)=    24690.5295574541
C      f(3,3)=      786.0977568405
C      f(3,4)=   -24122.7960738661
C      f(3,5)=     -259.2663294173
C      f(3,6)=     -598.5252855461

C      beta(3,1)=      48.9404306370
C      beta(3,2)=       5.2832855659
C      beta(3,3)=      16.8568900868
C      beta(3,4)=       5.2434883472
C      beta(3,5)=      89.2258188514
C      beta(3,6)=      25.7902108924
      
      
C new mode Aug 2012 start at 2 peak at 5
      !f(3,1)=    34039.2023345358
      !f(3,2)=       60.7894247475
      !f(3,3)=    36273.8783023426
      !f(3,4)=   -25286.2224420778
      !f(3,5)=   -44812.4245264180
      !f(3,6)=      460.6339790425

      !beta(3,1)=       3.8750997980
      !beta(3,2)=       2.4309215256
      !beta(3,3)=       3.3669595839
      !beta(3,4)=       3.2782539511
      !beta(3,5)=       3.7923429792
      !beta(3,6)=      10.7125633476      


C L=15

      !f(4,1)=           -410.3260290054
      !f(4,2)=           1168.8579712825
      !f(4,3)=         -98495.0510431024
      !f(4,4)=          98918.3033927164
      !f(4,5)=             18.2392467328
      !f(4,6)=           1105.3849889265

      !beta(4,1) =          48.1078229991
      !beta(4,2)  =         13.7459627225
      !beta(4,3)   =         4.9997421200
      !beta(4,4)    =        5.0077621877
      !beta(4,5)     =      89.3921544487
      !beta(4,6)      =     27.4857443437

C  L = 16
C  new mode Aug 2012 start at 2 peak at 5      
      !f(4,1)=           812.9918681787
      !f(4,2)=           441.2821316145
      !f(4,3)=          1301.6739693211
      !f(4,4)=          -143.4003781212
      !f(4,5)=            69.1750175418
      !f(4,6)=           175.9294532933

      !beta(4,1) =           48.1662482305
      !beta(4,2)  =           7.4414931962
      !beta(4,3)   =         15.8508455321
      !beta(4,4)    =         1.9564847414
      !beta(4,5)     =        1.7260534657
      !beta(4,6)      =       3.2147164196      



C L=14

      !f(4,1)=             -93.7744100883
      !f(4,2)=            1229.3279885934
      !f(4,3)=          -22458.0945711871
      !f(4,4)=           91059.8513053717
      !f(4,5)=          -70480.8918194039
      !f(4,6)=            2864.4839104274


      !beta(4,1) =          45.0729200914
      !beta(4,2)  =         20.3745381016
      !beta(4,3)   =         6.4705187514
      !beta(4,4)    =        5.9188525245
      !beta(4,5)     =       5.8018896209
      !beta(4,6)      =      8.5158398479

C L=16
      ! taken: use old mode again
C      f(4,1)=     -359.4857268465
C      f(4,2)=  -718318.9196872473
C      f(4,3)=      952.9729049448
C      f(4,4)=   718848.3348118553
C      f(4,5)=       29.9449422503
C      f(4,6)=     1422.1605650097

C      beta(4,1)=    48.9404373087
C      beta(4,2)=     5.4911324770
C      beta(4,3)=    13.9191935499
C      beta(4,4)=     5.4924655928
C      beta(4,5)=    89.2258188514
C      beta(4,6)=    25.7560973761

C L=32
C      ! taken: use old mode again
C      f(5,1)=     4285.5052331319
C      f(5,2)=  -600859.2054973751
C      f(5,3)=     1283.3234583067
C      f(5,4)=   601269.3489620079
C      f(5,5)=     -421.4303718400
C      f(5,6)=     3031.9399845112


C      beta(5,1)=  48.9404306370  
C      beta(5,2)=   5.1544511385  
C      beta(5,3)=  13.9220169737  
C      beta(5,4)=   5.1557936291  
C      beta(5,5)=  89.2258188514  
C      beta(5,6)=  25.7556569433  
      
C L=32 
C   new mode Aug 2012 start at 2 peak at 5 

      !f(5,1)=   5336.9276606427
      !f(5,2)=    239.4345954380
      !f(5,3)=    696.6336721901
      !f(5,4)=    -55.6358364446
      !f(5,5)=  -2581.4422206295
      !f(5,6)=   2848.1479327744


      !beta(5,1)=   48.2614188891 
      !beta(5,2)=    4.1173201057 
      !beta(5,3)=    9.4612155568 
      !beta(5,4)=    2.4503125863 
      !beta(5,5)=   89.1104591707 
      !beta(5,6)=   21.1589243453       

C L=64
      ! taken: use old mode again
C      f(6,1)=      9938.2744323159
C      f(6,2)=   -500764.5439833679
C      f(6,3)=      1287.6877590235
C      f(6,4)=    501147.4432669706
C      f(6,5)=     14912.2836873664
C      f(6,6)=      3963.5336376826

C      beta(6,1)=   48.9368813776
C      beta(6,2)=    5.1508684107
c      beta(6,3)=   13.9644917853
C      beta(6,4)=    5.1524290713
C      beta(6,5)=   89.2282691400
C      beta(6,6)=   25.7556037045
      
C L=64
C   new mode Aug 2012 start at 2 peak at 5 

      !f(6,1)=     9145.1487257103
      !f(6,2)=   -99402.7338129873
      !f(6,3)=     1370.4350614189
      !f(6,4)=    99735.9196457199
      !f(6,5)=    16595.1966646073
      !f(6,6)=     3784.6648254325

      !beta(6,1)=     46.9580300712
      !beta(6,2)=      4.9433214823
      !beta(6,3)=     13.6979549612
      !beta(6,4)=      4.9502924999
      !beta(6,5)=     88.8751266506
      !beta(6,6)=     25.9235159381     

C L=12

C       f(6,1)=   -171.2655001298

C       f(6,2)=    856.1629846247

C       f(6,3)= -41925.1558800406

C       f(6,4)=  94908.0078132388

C       f(6,5)= -54688.5441954559

C       f(6,6)=   2646.2528663652

C bi

C        beta(6,1)=  45.1551416237

C        beta(6,2)=  19.3395079889

C        beta(6,3)=   6.3257807905

C        beta(6,4)=   5.9930084444

C        beta(6,5)=   5.8046579879

C        beta(6,6)=   8.4762240472







	DATA NI/6/
C     NL number of modes
C     NI mumber of coefficients per mode


c	We work in the KSM coordinate system

c	In this subroutine, index I=1-6 corresponds to each component of a mode
c	Index L = 1,7 denotes the seven L modes
c	Index j sums over observations (Br and Btheta) and is twice the number of
c	observations
c	index M is used to obtain least squares equations (same dimension as L)
      PI=4.0*datan(1.D0)
	twopi=2.*PI
	PI2=PI/2.
      DEGREE=180./PI
	RADIAN=PI/180.
	drho=0.05
	dz=0.05
c	D is the half thickness of current sheet

	M=5
	M1=25
	M2=2*M
C
C     END DEFINITIONS    
C..........................................................................
c......................................................................
	R_S3=R1
	THETA=T1
	PHI=P1
	rho=R_S3*DSIN(THETA)
	time=time
c	RLT=(posLT(j)-12.0)*15.0
c	If (RLT .LT. 0.0) RLT=RLT+360.0
c	RLT=RLT*Radian
      XS3=R_S3*DSIN(THETA)*DCOS(PHI)
      YS3=R_S3*DSIN(THETA)*DSIN(PHI)
      ZS3=R_S3*DCOS(THETA)
C
C   NOW ROTATE the trajectory INTO KSM COORDINATES
C 

	!write(6,*)"Xs3: ",Xs3
	!write(6,*)"Ys3: ",Ys3
	!write(6,*)"Zs3: ",Zs3

      !write(6,*)"time: ",time,epoch


	pvecin(1)=xs3
	pvecin(2)=ys3
	pvecin(3)=zs3	
	CALL KROT('S3C','DIS',pvecin,pvecou,time,epoch)
	XDIS=pvecou(1)
	YDIS=pvecou(2)
	ZDIS=pvecou(3)

	RLT=datan2(YDIS,XDIS)


c	Start calculating the field for each unit mode	at the mapped location
c	First map the trajectory

	!write(6,*)"XDIS: ",XDIS
	!write(6,*)"YDIS: ",YDIS
	!write(6,*)"ZDIS: ",ZDIS



	call mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)


	Z=ZMAP
	RHOMAG=sqrt(XMAP**2+YMAP**2)

	Do 2 L = 1, NL

	ZM=dabs(Z-dz)
	If (ZM < D) ZM= 0.5*(ZM**2/D+D)
	ZP=dabs(Z+dz)
	If (ZP < D) ZP= 0.5*(ZP**2/D+D)
	xlpp=0.
	xlpm=0.
	Do 3 i=1,NI
	S1p=dsqrt((beta(L,i)/r0or+ZP)**2+(RHOMAG+Beta(L,i)/r0or)**2)
	S2p=dsqrt((beta(L,i)/r0or+ZP)**2+(RHOMAG-Beta(L,i)/r0or)**2)
	S1m=dsqrt((beta(L,i)/r0or+ZM)**2+(RHOMAG+Beta(L,i)/r0or)**2)
	S2m=dsqrt((beta(L,i)/r0or+ZM)**2+(RHOMAG-Beta(L,i)/r0or)**2)
	tp=2*(Beta(L,i)/r0or)/(S1p+S2p)
	tm=2*(Beta(L,i)/r0or)/(S1m+S2m)
	AAp=tp*dsqrt(1.-tp**2)/(S1p*S2p)
	AAm=tm*dsqrt(1.-tm**2)/(S1m*S2m)
	xlpp=xlpp+f(L,i)*AAp*rhomag
	xlpm=xlpm+f(L,i)*AAm*rhomag
    3 continue
	dxpldz=(xlpp-xlpm)/(2.*dz)
	!write(6,*)"dxpldz: ",dxpldz
	Xrho(L)=-dxpldz
    2	continue

c	Now CALCULATE THE BZ COMPONENT
	Do 20 L = 1, NL
	rhom=RHOMAG-drho
	rhop=RHOMAG+drho

	xi=dabs(Z)
	If (dabs(Z) .le. D) xi= 0.5*(Z**2/D+D)

c
	xlpp=0.
	xlpm=0.
	Do 30 i=1,NI
	S1p=dsqrt((beta(L,i)/r0or+xi)**2+(rhop+Beta(L,i)/r0or)**2)
	S2p=dsqrt((beta(L,i)/r0or+xi)**2+(rhop-Beta(L,i)/r0or)**2)
	S1m=dsqrt((beta(L,i)/r0or+xi)**2+(rhom+Beta(L,i)/r0or)**2)
	S2m=dsqrt((beta(L,i)/r0or+xi)**2+(rhom-Beta(L,i)/r0or)**2)
	tp=2*(Beta(L,i)/r0or)/(S1p+S2p)
	tm=2*(Beta(L,i)/r0or)/(S1m+S2m)
	AAp=tp*dsqrt(1.-tp**2)/(S1p*S2p)
	AAm=tm*dsqrt(1.-tm**2)/(S1m*S2m)
	xlpp=xlpp+f(L,i)*AAp*rhop
	xlpm=xlpm+f(L,i)*AAm*rhom
   30 continue
	dxpldr=(rhop*xlpp-rhom*xlpm)/(2.*drho)
	!write(6,*)"dxpldr: ",dxpldr
	Xz(L)=dxpldr/RHOMAG
	!write(6,*)"RHOMAG: ",RHOMAG
   20	continue


c
c	Now calculate the new mapped field for each  UNIT mode
	Bx=0.0
	By=0.0
	Bz=0.0

C	Do 60 L=5,5
	Do 60 L=1,NL
	if(( ymap .eq. 0.) .and. (xmap .eq. 0.)) then
	phimap=0.
	go to 69
	end if
	phimap=atan2(YMAP,XMAP)
   69	Bx1=Xrho(L)*Dcos(phimap)
	By1=Xrho(L)*Dsin(phimap)
	Bz1=Xz(L)

	!write(6,*)"Bx1: ",Bx1
	!write(6,*)"By1: ",By1
	!write(6,*)"Bz1: ",Bz1



c	Now Add the ACTUAL shield field to the current sheet field
	XMAP1=XMAP 
	YMAP1=YMAP 
	ZMAP1=ZMAP 
	call shieldfield_csheet(XMAP1,YMAP1,ZMAP1,L,M,M1,M2,
     ~BX2,BY2,BZ2,r0or)
     

!      BX2=0.0
!	BY2=0.0
!	BZ2=0.0




	Bx=Bx+C(L)*r0or**dql*(Bx1*(1/(r0or))**2.0-Bx2/r0or)
	By=By+C(L)*r0or**dql*(By1*(1/(r0or))**2.0-By2/r0or)
	Bz=Bz+C(L)*r0or**dql*(Bz1*(1/(r0or))**2.0-Bz2/r0or)

      !for writing out single shielded tail modes
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	Bxi=C(L)*r0or**dql*(Bx1*(1/(r0or))**2.0-Bx2/r0or)
	Byi=C(L)*r0or**dql*(By1*(1/(r0or))**2.0-By2/r0or)
	Bzi=C(L)*r0or**dql*(Bz1*(1/(r0or))**2.0-Bz2/r0or)
	
	call Mapped_field(
     ~	time,XMAP,YMAP,ZMAP,Bxi,Byi,Bzi,Bxdisi,Bydisi,Bzdisi,r0or,Epoch)
c
	bvecin(1)=Bxdisi
	bvecin(2)=Bydisi
	bvecin(3)=Bzdisi

	!write(6,*)"Bxdis: ",Bxdis
	!write(6,*)"Bydis: ",Bydis
	!write(6,*)"Bzdis: ",Bzdis

c	Now rotate back into spherical coordinate system

	
	CALL KROT('DIS','S3C',bvecin,bvecou,time,epoch)

	Bxs3i=bvecou(1)
	Bys3i=bvecou(2)
	Bzs3i=bvecou(3)

	!write(6,*)"Bxs3: ",Bxs3
	!write(6,*)"Bys3: ",Bys3
	!write(6,*)"Bzs3: ",Bzs3

	Brm1arr(L)=Bxs3i*sin(theta)*cos(phi)+Bys3i*sin(theta)*sin(phi)+
     ~Bzs3i*cos(theta)
	Btm1arr(L)=Bxs3i*cos(theta)*cos(phi)+Bys3i*cos(theta)*sin(phi)-
     ~Bzs3i*sin(theta)
	Bpm1arr(L)=-Bxs3i*sin(phi)+Bys3i*cos(phi)
	
	!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	!end for writing out single tail modes
	
	
	
	!write(6,*)"Bx: ",Bx
	!write(6,*)"By: ",By
	!write(6,*)"Bz: ",Bz


   60 Continue
c
	!write(6,*)"XMAP: ",XMAP
	!write(6,*)"YMAP: ",YMAP
	!write(6,*)"ZMAP: ",ZMAP


	call Mapped_field(
     ~	time,XMAP,YMAP,ZMAP,Bx,By,Bz,Bxdis,Bydis,Bzdis,r0or,Epoch)
c
	bvecin(1)=Bxdis
	bvecin(2)=Bydis
	bvecin(3)=Bzdis

	!write(6,*)"Bxdis: ",Bxdis
	!write(6,*)"Bydis: ",Bydis
	!write(6,*)"Bzdis: ",Bzdis

c	Now rotate back into spherical coordinate system

	
	CALL KROT('DIS','S3C',bvecin,bvecou,time,epoch)

	Bxs3=bvecou(1)
	Bys3=bvecou(2)
	Bzs3=bvecou(3)

	!write(6,*)"Bxs3: ",Bxs3
	!write(6,*)"Bys3: ",Bys3
	!write(6,*)"Bzs3: ",Bzs3

	Brm1=Bxs3*sin(theta)*cos(phi)+Bys3*sin(theta)*sin(phi)+
     ~Bzs3*cos(theta)
	Btm1=Bxs3*cos(theta)*cos(phi)+Bys3*cos(theta)*sin(phi)-
     ~Bzs3*sin(theta)
	Bpm1=-Bxs3*sin(phi)+Bys3*cos(phi)

!	!write(6,*)"Brm1: ",Brm1
!	!write(6,*)"Btm1: ",Btm1
!	!write(6,*)"Bpm1: ",Bpm1


	Return
	End

C********************************************************************
	Subroutine shieldfield_csheet(X,Y,Z,L,M,M1,M2,
     ~Bxm,Bym,Bzm,r0or)
	IMPLICIT REAL*8(A-H,O-Z)
c

c
	Dimension ABCD(6,25),PQ(6,10)
	Dimension a(M,M),p(M2)

	Dimension b(M,M),c(M,M),d(M,M)
	Dimension q(M2),r(M2),s(M2)
	Dimension XXP(M1),XXM(M1)

C     MODE 1 Aug 2012 peak at 5
C	Data (abcd (1, j), j = 1, 25) /
C     ~    -0.35153006639006954571D-05,
C     ~     0.23575800697128829436D-03,
C     ~    -0.68523032683502002604D-02,
C     ~     0.50584160772160019492D-01,
C     ~    -0.75694100518610099582D-01,
C     ~     0.22456816943470985577D-03,
C     ~    -0.10445488138626792506D-02,
C     ~     0.80112778867057421461D-01,
C     ~    -0.64183468172852997657D+00,
C     ~     0.98809223559946512960D+00,
C     ~    -0.52762855732372444806D-02,
C     ~    -0.30779140045062824971D-01,
C     ~    -0.16265439828762690011D-01,
C     ~     0.14316815615149347529D+00,
C     ~    -0.43849234995754932109D+00,
C     ~     0.22216049621179477513D-01,
C     ~     0.13161835247679370964D+00,
C     ~     0.99049377626076032244D-01,
C     ~     0.26013968527285968868D+00,
C     ~     0.32873575249230019280D+00,
C     ~    -0.17098828184387373369D-01,
C     ~    -0.10079006634277780374D+00,
C     ~    -0.11356834871457452407D+00,
C     ~    -0.29911921101585359217D+00,
C     ~    -0.48562354745533009570D+00/
C      	Data (pq (1, j), j = 1, 10) /
C     ~        0.40000081062316894531D+01,
C     ~        0.80000104904174804688D+01,
C     ~        0.16000005722045898438D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02,
C     ~        0.39998939037322998047D+01,
C     ~        0.79999985694885253906D+01,
C     ~        0.16000000000000000000D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02/


C     MODE 2 Aug 2012 peak at 5
C	Data (abcd (2, j), j = 1, 25) /
C     ~     0.10678961640307714712D-04,
C     ~    -0.79277805987504465154D-03,
C     ~     0.20251500952258305599D-01,
C     ~    -0.14263622699416847461D+00,
C     ~     0.21053465945964125461D+00,
C     ~    -0.18532209610370927204D-03,
C     ~     0.53288622109895090323D-02,
C     ~    -0.23662576366437076958D+00,
C     ~     0.18459082306133383078D+01,
C     ~    -0.28062869665645946426D+01,
C     ~     0.36776835832946536305D-02,
C     ~     0.34357790313751211675D-01,
C     ~    -0.23388070299514979178D-01,
C     ~     0.33422041691935078234D-02,
C     ~     0.33369653661920106735D+00,
C     ~    -0.14215363077573411194D-01,
C     ~    -0.11038974961262414121D+00,
C     ~    -0.98945745397397172338D-01,
C     ~    -0.52453160708230328613D+00,
C     ~    -0.27580417633032239255D+00,
C     ~     0.10387043558141385857D-01,
C     ~     0.68529221559159858579D-01,
C     ~     0.26274248425618644542D-01,
C     ~    -0.35336925055235512971D+00,
C     ~     0.22982151697208269558D+00/
C      	Data (pq (2, j), j = 1, 10) /
C     ~        0.39999945163726806641D+01,
C     ~        0.79999995231628417969D+01,
C     ~        0.16000001907348632812D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02,
C     ~        0.39993464946746826172D+01,
C     ~        0.80000047683715820313D+01,
C     ~        0.15999997138977050781D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02/

C     MODE 3 Aug 2012 peak at 5
C	Data (abcd (3, j), j = 1, 25) /
C     ~     0.56805659376579659158D-04,
C     ~    -0.49826416171172595304D-02,
C     ~     0.12526839464659417223D+00,
C     ~    -0.87409583160880188757D+00,
C     ~     0.12856485167171063377D+01,
C     ~     0.11216145207766179542D-02,
C     ~     0.35698754430595128206D-01,
C     ~    -0.11826555719375608788D+01,
C     ~     0.91620587898341590005D+01,
C     ~    -0.13838867662866388031D+02,
C     ~    -0.33466353340866311639D-01,
C     ~    -0.94500618941237457538D-01,
C     ~    -0.48684028075674673497D+00,
C     ~     0.14607883360870395073D+01,
C     ~    -0.13643921002377688012D+01,
C     ~     0.15327105053855760719D+00,
C     ~     0.77253223257713465877D+00,
C     ~     0.66066721067945421986D+00,
C     ~    -0.86601218815076608237D+00,
C     ~     0.10860893867610756924D+01,
C     ~    -0.12297786782076687573D+00,
C     ~    -0.72241771816889999158D+00,
C     ~    -0.77546966054031540860D+00,
C     ~    -0.22509195030543880378D+01,
C     ~    -0.19436733171165148093D+01/
C      	Data (pq (3, j), j = 1, 10) /
C     ~        0.39999673366546630859D+01,
C     ~        0.79999799728393554688D+01,
C     ~        0.16000000000000000000D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02,
C     ~        0.39958155155181884766D+01,
C     ~        0.80000371932983398438D+01,
C     ~        0.15999984741210937500D+02,
C     ~        0.32000000000000000000D+02,
C     ~        0.64000000000000000000D+02/


C     MODE 4 Aug 2012 peak at 5
C	Data (abcd (4, j), j = 1, 25) /
C     ~     0.19638341783294487237D-03,
C     ~    -0.19807955436661750515D-01,
C     ~     0.50533167084825092452D+00,
C     ~    -0.35380207109579977498D+01,
C     ~     0.52069369271108874031D+01,
C     ~     0.83297371091017383127D-02,
C     ~     0.12861211658879426190D+00,
C     ~    -0.37415022082172026963D+01,
C     ~     0.28872340578355697005D+02,
C     ~    -0.43434963349345771633D+02,
C     ~    -0.22717934357869418172D+00,
C     ~    -0.81333065020493922059D+00,
C     ~    -0.26455667944533032276D+01,
C     ~     0.11153823484896403073D+02,
C     ~    -0.16209356038027308955D+02,
C     ~     0.10139287180009441958D+01,
C     ~     0.51585299887162561205D+01,
C     ~     0.42923466779422323469D+01,
C     ~    -0.21809683187273272154D+01,
C     ~     0.14371299951134901107D+02,
C     ~    -0.80264576829420197779D+00,
C     ~    -0.44983967922243719428D+01,
C     ~    -0.42669104398308794757D+01,
C     ~    -0.83595341151041111516D+01,
C     ~    -0.16784350456480702007D+02/
C      	Data (pq (4, j), j = 1, 10) /
C     ~        0.39999535083770751953D+01,
C     ~        0.79999637603759765625D+01,
C     ~        0.16000024795532226562D+02,
C     ~        0.32000011444091796875D+02,
C     ~        0.64000007629394531250D+02,
C     ~        0.39839067459106445312D+01,
C     ~        0.80001392364501953125D+01,
C     ~        0.15999942779541015625D+02,
C     ~        0.31999998092651367188D+02,
C     ~        0.64000000000000000000D+02/




C     MODE 5 Aug 2012 peak at 5
C	Data (abcd (5, j), j = 1, 25) /
C     ~     0.57602181502077927514D-08,
C     ~    -0.27198178058757190846D-06,
C     ~     0.14413677102796104145D-04,
C     ~    -0.11014990050355253606D-03,
C     ~     0.16067283476392421300D-03,
C     ~     0.13193945273706952556D-05,
C     ~    -0.89828918286282105447D-04,
C     ~     0.21515136058656212509D-01,
C     ~    -0.19682427702020968696D+00,
C     ~     0.32279192414489094043D+00,
C     ~    -0.67184429268400664261D-04,
C     ~    -0.14284490361830121802D-01,
C     ~    -0.35629192133510734219D+00,
C     ~     0.30705765682217287527D+01,
C     ~    -0.61419772176518296902D+01,
C     ~     0.38732153569629476468D-03,
C     ~     0.11458722673551952065D+00,
C     ~     0.17929207140451239866D+01,
C     ~    -0.89188944523992432778D+01,
C     ~     0.32228805685562107897D+02,
C     ~    -0.32942021863423845375D-03,
C     ~    -0.11155626404876511915D+00,
C     ~    -0.19432724644718772922D+01,
C     ~     0.46048240963900362033D+01,
C     ~    -0.52986016965987758454D+02/
C      	Data (pq (5, j), j = 1, 10) /
C     ~        0.20007674694061279297D+01,
C     ~        0.40067391395568847656D+01,
C     ~        0.80072870254516601563D+01,
C     ~        0.16005931854248046875D+02,
C     ~        0.32003528594970703125D+02,
C     ~        0.19872297048568725586D+01,
C     ~        0.39964907169342041016D+01,
C     ~        0.80006351470947265625D+01,
C     ~        0.16000080108642578125D+02,
C     ~        0.32000152587890625000D+02/

C     MODE 6 Aug 2012 peak at 5
C	Data (abcd (6, j), j = 1, 25) /
C     ~     0.12600037039950878022D-03,
C     ~    -0.10286886633016696646D-01,
C     ~    -0.12772309288678795891D+00,
C     ~     0.19989257632239989526D+01,
C     ~    -0.34321939046075433666D+01,
C     ~     0.55499013547815743974D-01,
C     ~     0.31173553004160692304D-01,
C     ~     0.61195985306330076270D+01,
C     ~    -0.52035160741602140888D+02,
C     ~     0.82976563376006623685D+02,
C     ~    -0.13749201131347199567D+01,
C     ~    -0.53997441327526773946D+01,
C     ~    -0.33885918241222711345D+02,
C     ~     0.22880910222910569018D+03,
C     ~    -0.38340140380148346821D+03,
C     ~     0.59693886907215576443D+01,
C     ~     0.31333698356997615520D+02,
C     ~     0.59768185557434726718D+02,
C     ~    -0.17342429700241308410D+03,
C     ~     0.45716000890957502634D+03,
C     ~    -0.46705148472034592189D+01,
C     ~    -0.26291656938929079956D+02,
C     ~    -0.47428705936755726214D+02,
C     ~     0.50799236585669959254D+02,
C     ~    -0.36579301953188536345D+03/
C     	 Data (pq (6, j), j = 1, 10) /
C     ~        0.40019817352294921875D+01,
C     ~        0.80021457672119140625D+01,
C     ~        0.16002088546752929688D+02,
C     ~        0.32001159667968750000D+02,
C     ~        0.64000549316406250000D+02,
C     ~        0.39563248157501220703D+01,
C     ~        0.80006628036499023438D+01,
C     ~        0.15999974250793457031D+02,
C     ~        0.32000030517578125000D+02,
C     ~        0.64000007629394531250D+02/



C	MODE 1 Peak at 6
C	Data (abcd (1, j), j = 1, 25) / 
C     ~    0.57900389912988492469D-05,
C     ~   -0.33867095260193265104D-02,
C     ~   -0.17885402929184579079D+00,
C     ~   -0.76105861802226257850D-01,
C     ~   -0.21707320837598329532D-01,
C     ~   -0.22806696336541141256D-04,
C     ~    0.13605855599644953723D-01,
C     ~    0.72655745785709626716D+00,
C     ~    0.33979077925435658968D+00,
C     ~    0.84015518189067677212D-01,
C     ~    0.38251019250246102387D-04,
C     ~   -0.23331745620323913747D-01,
C     ~   -0.12632616230609083896D+01,
C     ~   -0.66501059469106635901D+00,
C     ~   -0.11557383129450049530D+00,
C     ~   -0.46761660834837321942D-04,
C     ~    0.29188863483653651798D-01,
C     ~    0.16062192067320260058D+01,
C     ~    0.10008527972927729887D+01,
C     ~    0.96874076440719001368D-01,
C     ~    0.25527293306380216542D-04,
C     ~   -0.16076374426793411664D-01,
C     ~   -0.89077972230583615242D+00,
C     ~   -0.61419050551211071820D+00,
C     ~   -0.48351231840305022302D+00/

C     MODE 1 Peak at 4 (valid before Aug 2012)
C	Data (abcd (1, j), j = 1, 25) /
C     ~   -0.14004329253036342350D-07,
C     ~    0.37688393599679170797D-05,
C     ~    0.49241640240943258532D-04,
C     ~    0.27365751661153279172D-03,
C     ~   -0.29637572764194883845D-03,
C     ~    0.23608493437412536586D-06,
C     ~   -0.99590524882565070186D-04,
C     ~   -0.13716499714981711388D-02,
C     ~   -0.68651650806999411358D-02,
C     ~    0.13086363094019834996D-01,
C     ~   -0.17741658850262338020D-05,
C     ~    0.95936069671529971003D-03,
C     ~    0.14635178751353254966D-01,
C     ~    0.13851244022798143618D-01,
C     ~   -0.50732363273520704183D-01,
C     ~    0.59728889093050563374D-05,
C     ~   -0.35847147114331114892D-02,
C     ~   -0.60417634422173911445D-01,
C     ~   -0.24014471471244385192D-02,
C     ~   -0.87454103956755879778D-02,
C     ~   -0.44206961724389151058D-05,
C     ~    0.27227339013033198256D-02,
C     ~    0.47812648728699302935D-01,
C     ~   -0.49208189730297080544D-02,
C     ~    0.84257347502468000755D-02/

C     PEAK at 8 start at 5 (new mode Aug 2012)     
C      Data (abcd (1, j), j = 1, 25) /
C     ~   -0.22806650633874154589D-04,
C     ~   0.19981441572494572541D-02 ,
C     ~  -0.50917717624770160445D-01 ,
C     ~   0.35714413408705841801D+00 ,
C     ~  -0.52608610022227830783D+00 ,
C     ~  -0.39669626460951413474D-03 ,
C     ~  -0.14042301137871120859D-01 ,
C     ~   0.47028862177184993509D+00 ,
C     ~  -0.36578470766660089453D+01 ,
C     ~   0.55325325060575014646D+01 ,
C     ~   0.12202297433719602737D-01 ,
C     ~   0.35921039006745036914D-01 ,
C     ~   0.19367172629944376383D+00 ,
C     ~  -0.47038727079262165898D+00 ,
C     ~   0.31050626588489294022D+00 ,
C     ~  -0.56170941159746065607D-01 ,
C     ~  -0.29912370171159619048D+00 ,
C     ~  -0.30322589282814016443D+00 ,
C     ~   0.27108065541804371446D+00 ,
C     ~  -0.34160605247767883430D+00 ,
C     ~   0.45136745871513728190D-01 ,
C     ~   0.27576166498336385224D+00 ,
C     ~   0.23869604996261831187D+00 ,
C     ~   0.18925998312564620840D+00 ,
C     ~   0.44434492288249283209D+00 /
C       
       Data (abcd (1, j), j = 1, 25) /
     ~   0.52034884821428676664D-08,
     ~  -0.68983393140626007927D-06 ,
     ~   0.26013060818085627384D-04 ,
     ~  -0.20220185555762154508D-03 ,
     ~   0.30564180727029997440D-03 ,
     ~   0.19930189328513339478D-06 ,
     ~  -0.12444044081689495742D-04 ,
     ~   0.14540867652184175443D-02 ,
     ~  -0.12059368770663559875D-01 ,
     ~   0.18972457961433673856D-01 ,
     ~  -0.12712315070435666379D-04 ,
     ~  -0.10948713314089635051D-02 ,
     ~  -0.11325779706928780499D-01 ,
     ~   0.10236958766108690544D+00 ,
     ~  -0.19081140235653648585D+00 ,
     ~   0.78157561366167476293D-04 ,
     ~   0.81610134021544745997D-02 ,
     ~   0.34775577517361599966D-01 ,
     ~  -0.13775781931913483325D+00 ,
     ~   0.59210130685600637079D+00 ,
     ~  -0.67653598302448326595D-04 ,
     ~  -0.77757671746549184011D-02 ,
     ~  -0.43300045770593376115D-01 ,
     ~  -0.72121864531289603395D-01 ,
     ~  -0.10966408811380061739D+01 /
        
                                  

      Data (abcd (2, j), j = 1, 25) /
     ~   0.12545399310341492571D-08,
     ~  -0.11363771161960078800D-06,
     ~   0.55343644774000541920D-05,
     ~  -0.44032796493657441067D-04,
     ~   0.66582238656644651925D-04,
     ~   0.55943482168823719799D-07,
     ~  -0.97499376833970518089D-05,
     ~   0.18348839506297323111D-02,
     ~  -0.17055331333206070338D-01,
     ~   0.27862314450439854691D-01,
     ~  -0.49597135329461795178D-05,
     ~  -0.11392557185046436278D-02,
     ~  -0.24260911596492629183D-01,
     ~   0.24245633528367632747D+00,
     ~  -0.49494927634042140685D+00,
     ~   0.33779783670884271877D-04,
     ~   0.10006602582479727503D-01,
     ~   0.12919577290627834443D+00,
     ~  -0.61024568956482516668D+00,
     ~   0.27989901087956345904D+01,
     ~  -0.30128015057845978483D-04,
     ~  -0.10625414457624882125D-01,
     ~  -0.18652399199560046328D+00,
     ~  -0.47424636188494684319D-01,
     ~  -0.61415416903532733528D+01/
     

     
C
C     ###################################
C     New mode Aug 2012
C     MODE 1 Start at 2; Peak at 6
C     ###################################
C	Data (abcd (1, j), j = 1, 25) /
C     ~  -0.13774757805365866521D-05,
C     ~  -0.70971199764301553040D-03,
C     ~   0.27634465364814016219D-01,
C     ~  -0.22098384462359346725D+00,
C     ~   0.33708836353742760439D+00,
C     ~   0.10832025970333686131D-02,
C     ~  -0.53530485360866290889D-02,
C     ~   0.18944184087282156526D+00,
C     ~  -0.12716464964495621803D+01,
C     ~   0.18215778704469625104D+01,
C     ~  -0.26286979174857352803D-01,
C     ~  -0.48444880496650646884D-01,
C     ~  -0.52660660932882952245D+00,
C     ~   0.45082237972458321096D+01,
C     ~  -0.65842130925306410205D+01,
C     ~   0.11411754551073770292D+00,
C     ~   0.31979110223247919631D+00,
C     ~   0.19856752256512044674D+00,
C     ~  -0.49506264985306094317D+01,
C     ~   0.77412277628375987248D+01,
C     ~  -0.89074049080843556436D-01,
C     ~  -0.26584905850947632144D+00,
C     ~   0.29268786738959606808D-02,
C     ~   0.22086426652376278668D+01,
C     ~  -0.47584288669132392968D+01/

C     ###################################
C     New mode Aug 2012
C     MODE 1 Start at 1.2; Peak at 5
C     ###################################
C	Data (abcd (1, j), j = 1, 25) /
C     ~    0.27520534657227329261D-06,
C     ~   -0.56479711896642841723D-03,
C     ~    0.20710029128770408857D-01,
C     ~   -0.16323345825271243226D+00,
C     ~    0.24811642320274354678D+00,
C     ~    0.71681996468643278259D-03,
C     ~   -0.30471710664889337106D-02,
C     ~    0.10516844777484216089D+00,
C     ~   -0.68502745701355149954D+00,
C     ~    0.96848183576033330944D+00,
C     ~   -0.17509138573501854125D-01,
C     ~   -0.30505093121171573262D-01,
C     ~   -0.36239860918944366563D+00,
C     ~    0.30638637104303150238D+01,
C     ~   -0.44606235835225174213D+01,
C     ~    0.76215353970331842226D-01,
C     ~    0.21157611817542446020D+00,
C     ~    0.14114391054107056078D+00,
C     ~   -0.33790814687815076134D+01,
C     ~    0.52659027560436868498D+01,
C     ~   -0.59563477929664765398D-01,
C     ~   -0.17737804165846532412D+00,
C     ~   -0.45599955713709894312D-03,
C     ~    0.15056277265595638948D+01,
C     ~   -0.32357335890645630982D+01/

C	MODE 2
C	Data (abcd (2, j), j = 1, 25) /
C     ~    0.77668182761180606377D-08,
C     ~   -0.20196291280988356575D-05,
C     ~   -0.49273322968357806672D-04,
C     ~    0.21539440444312103473D-03,
C     ~   -0.83705446882959897436D-04,
C     ~    0.33536606594698534777D-04,
C     ~   -0.12866462057380196881D-01,
C     ~   -0.31621223228494135248D+00,
C     ~    0.61757501135689674143D+00,
C     ~   -0.11476830385717948779D+00,
C     ~   -0.17431745682254278229D-04,
C     ~    0.65663659425452749474D-02,
C     ~    0.15888921816088288352D+00,
C     ~   -0.33062421999133104755D+00,
C     ~    0.13245199895207815377D+00,
C     ~   -0.26812881293819263994D-04,
C     ~    0.10672183150391700845D-01,
C     ~    0.27287667133660016283D+00,
C     ~   -0.42942834594719494489D+00,
C     ~   -0.32044936014889375819D+00,
C     ~    0.10700204255233096706D-04,
C     ~   -0.43707712411099670646D-02,
C     ~   -0.11613927454426622443D+00,
C     ~    0.71015669683247804044D-01,
C     ~   -0.13462690359270563789D+01/
 
C	MODE 2
	Data (abcd (3, j), j = 1, 25) /
     ~   0.24824112592626616873D-08,
     ~  -0.22636768136616835943D-06,
     ~   0.11627747929157943881D-04,
     ~  -0.93185828541469004570D-04,
     ~   0.14182203142413209369D-03,
     ~  -0.22688551303297704850D-06,
     ~  -0.46061637003510315412D-04,
     ~   0.22179188550876874353D-02,
     ~  -0.21306193978350578683D-01,
     ~   0.34153314431143159957D-01,
     ~   0.57496428715017037080D-05,
     ~   0.97441205991077970046D-03,
     ~  -0.18779683259496859904D-01,
     ~   0.32105889582564206552D+00,
     ~  -0.60933710405131935595D+00,
     ~  -0.15042258847903980102D-04,
     ~  -0.20826093187111503757D-02,
     ~   0.73347474035392712377D-01,
     ~  -0.10606675688097344512D+01,
     ~   0.37511030714926558716D+01,
     ~   0.76517215572531172669D-05,
     ~  -0.16417219343249312910D-02,
     ~  -0.18174979954668427085D+00,
     ~   0.16023634510000572817D+00,
     ~  -0.87233263461601922018D+01/

   
   
   
   

C	MODE 3
C	Data (abcd (3, j), j = 1, 25) /
C     ~    0.92353191216470111868D-02,
C     ~   -0.47065512098254735917D+01,
C     ~   -0.22554133829293983026D+03,
C     ~   -0.23457137432873151894D+02,
C     ~    0.46111664280590520803D+02,
C     ~   -0.12239552063973810902D-01,
C     ~    0.62483152611269074938D+01,
C     ~    0.29997204397613383974D+03,
C     ~    0.31604341484015141539D+02,
C     ~   -0.57823703323961463951D+02,
C     ~    0.41170674967848803760D-02,
C     ~   -0.21219401709444536408D+01,
C     ~   -0.10293144491195860279D+03,
C     ~   -0.11947843043877295343D+02,
C     ~    0.14488934852695538602D+02,
C     ~   -0.18939743624590780868D-02,
C     ~    0.99324763584517850034D+00,
C     ~    0.49141758650519626883D+02,
C     ~    0.77481048808290280405D+01,
C     ~   -0.39887477569330695104D+01,
C     ~    0.78113976972219205663D-03,
C     ~   -0.41307234546100435323D+00,
C     ~   -0.20642069925114081563D+02,
C     ~   -0.40978583394053940125D+01,
C     ~   -0.35575024161355539575D+01/

C	MODE 3
	Data (abcd (4, j), j = 1, 25) /
     ~   0.27771830099459217676D-06,
     ~  -0.81341899879401334040D-04,
     ~   0.29282572672455508579D-02,
     ~  -0.22583395250099128071D-01,
     ~   0.34010285131856290985D-01,
     ~   0.99291569751668582507D-04,
     ~   0.26052436207605949789D-02,
     ~   0.10471966853985351098D+00,
     ~  -0.83263064611646853130D+00,
     ~   0.13428602164692611254D+01,
     ~  -0.36787553054767226091D-02,
     ~  -0.26316046069385012895D+00,
     ~  -0.93764343790019077929D+00,
     ~   0.50555823544919160639D+01,
     ~  -0.11196287302734999614D+02,
     ~   0.18441461144297600916D-01,
     ~   0.16391989138877585130D+01,
     ~   0.46270302869078419050D+01,
     ~   0.18359289360326336293D+01,
     ~   0.32155291790796020734D+02,
     ~  -0.15031144448049583301D-01,
     ~  -0.14580925618602209592D+01,
     ~  -0.54630082835657294282D+01,
     ~  -0.12099960332237163385D+02,
     ~  -0.53861339558258123361D+02/
     
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

	 
C	MODE 4
C	Data (abcd (4, j), j = 1, 25) /
C     ~    0.21992617916451395743D+02,
C     ~   -0.11399307544820429516D+05,
C     ~   -0.76667036088618019107D+06,
C     ~   -0.39185539354281475610D+06,
C     ~   -0.59239538019716295735D+03,
C     ~    0.47524457523937337910D+00,
C     ~   -0.24323272755420672908D+03,
C     ~   -0.16297746944208286734D+05,
C     ~   -0.78974933554931041612D+04,
C     ~   -0.34576773059694034806D+02,
C     ~   -0.67509882718360936237D+00,
C     ~    0.34607166925607515395D+03,
C     ~    0.23198896555695553978D+05,
C     ~    0.11312410520290614446D+05,
C     ~    0.45910310222559074233D+02,
C     ~   -0.21888994250814901576D+02,
C     ~    0.11346764342961157545D+05,
C     ~    0.76316098499373987706D+06,
C     ~    0.39024560870819433588D+06,
C     ~    0.58817978798588868016D+03,
C     ~    0.96230585996413289251D-01,
C     ~   -0.50295745279677515071D+02,
C     ~   -0.33917785686196171290D+04,
C     ~   -0.18056934409611232084D+04,
C     ~   -0.22829873613259561437D+02/

C	MODE 4
C	Data (abcd (5, j), j = 1, 25) /
C     ~  -0.45915862230592986108D-03,
C     ~   0.93854011964554895364D-01,
C     ~  -0.38314667487779807686D+01,
C     ~   0.30828311901267987594D+02,
C     ~  -0.47112777682193481610D+02,
C     ~   0.53309485417775998028D-01,
C     ~  -0.18620821471067905417D+00,
C     ~   0.24001863168454320885D+02,
C     ~  -0.21574676433475329418D+03,
C     ~   0.34315573083493438844D+03,
C     ~  -0.13148027970186495494D+01,
C     ~  -0.76415371208712992157D+01,
C     ~  -0.21148444432406684257D+02,
C     ~   0.23239199306603961759D+03,
C     ~  -0.48969517910810122885D+03,
C     ~   0.56621555612493432719D+01,
C     ~   0.35897685020618332885D+02,
C     ~   0.40435286626670226440D+02,
C     ~  -0.91416647341691259498D+01,
C     ~   0.50909741444276653510D+03,
C     ~  -0.44092277778652038123D+01,
C     ~  -0.28613448638516562283D+02,
C     ~  -0.44306202263644621553D+02,
C     ~  -0.53377514196322152884D+02,
C     ~  -0.47645787621277645485D+03/
 
C 	Data (abcd (5, j), j = 1, 25) /    
C     ~   0.29184501829407924565D-03,
C     ~  -0.35174081943346383050D-01,
C     ~   0.82480246693769743160D+00,
C     ~  -0.55640837422291227199D+01,
C     ~   0.80920395915216154492D+01,
C     ~   0.33524297003114426896D-01,
C     ~   0.27659122634700405996D+00,
C     ~  -0.51499706442952035701D+01,
C     ~   0.37955245686153382678D+02,
C     ~  -0.55685624735124136464D+02,
C     ~  -0.86509761148266395292D+00,
C     ~  -0.37778516452277250970D+01,
C     ~  -0.74686459122406816391D+01,
C     ~   0.37523091820409682384D+02,
C     ~  -0.67750312112703170442D+02,
C     ~   0.37807149835021629514D+01,
C     ~   0.20453142129327371634D+02,
C     ~   0.15788458773079129926D+02,
C     ~   0.88840198407937549518D+01,
C     ~   0.59253145460556318369D+02,
C     ~  -0.29649080045840339181D+01,
C     ~  -0.17103608884133656431D+02,
C     ~  -0.15987670755539895140D+02,
C     ~  -0.31371844952042430066D+02,
C     ~  -0.71049048051729997155D+02/
 
 
  	Data (abcd (5, j), j = 1, 25) /        
     ~  -0.99911155916475476429D-04,
     ~  -0.16497770784682154899D-02,
     ~  -0.43294948559606633021D+00,
     ~   0.44153048101428264971D+01,
     ~  -0.71275189851081943715D+01,
     ~   0.79586008001811162083D-01,
     ~   0.41014729592350235299D+00,
     ~   0.16633384478053871991D+01,
     ~  -0.22369799721501067324D+02,
     ~   0.41291100537732880582D+02,
     ~  -0.19782143403119247793D+01,
     ~  -0.11556993496696394530D+02,
     ~  -0.81831519476911562805D+01,
     ~   0.11663048191732711345D+02,
     ~  -0.74936198497502246596D+02,
     ~   0.84827948435163005314D+01,
     ~   0.55050326212560698025D+02,
     ~   0.47358311675180331690D+02,
     ~   0.14227065859524753932D+03,
     ~   0.27681529848104631242D+02,
     ~  -0.65980337848883907625D+01,
     ~  -0.44385714456835309250D+02,
     ~  -0.50397503517239158555D+02,
     ~  -0.12514882115042821908D+03,
     ~  -0.12233187994479561667D+03/

     
   

     


C     L=32
C	Data (abcd (5, j), j = 1, 25)  / 
C     ~       0.38370706169389450224D+02,
C     ~      -0.32479962500340833209D+05,
C     ~      -0.94660043127081774372D+06,
C     ~      -0.56864021527375161468D+06,
C     ~      -0.27814993022411260703D+04,
C     ~       0.20600370116089345984D+00,
C     ~      -0.17025940328754780139D+03,
C     ~      -0.49395843669987780089D+04,
C     ~      -0.30000800342043403290D+04,
C     ~      -0.97079695288155960497D+01,
C     ~      -0.14579821620379142643D+01,
C     ~       0.12204100420071475330D+04,
C     ~       0.35488135357577208495D+05,
C     ~       0.21388387530527799285D+05,
C     ~       0.76407552757476420168D+02,
C     ~      -0.37682332024326319341D+02,
C     ~       0.31912372589853306159D+05,
C     ~       0.93015110639191078689D+06,
C     ~       0.55873046694301766734D+06,
C     ~       0.28055544715537447331D+04,
C     ~       0.56360430778218297831D+00,
C     ~      -0.48256080155038478807D+03,
C     ~      -0.14099258801109653127D+05,
C     ~      -0.84807776710395010866D+04,
C     ~      -0.13800109969955343203D+03/




C     L=64
	Data (abcd (6, j), j = 1, 25)  / 
     ~   -0.88135632848555722773D+01,
     ~    0.11514098078690224724D+05,
     ~    0.21560988124946391941D+06,
     ~    0.45832673059045454522D+05,
     ~    0.21572565347949170799D+03,
     ~    0.54760001008211061090D+00,
     ~   -0.70240562027308603987D+03,
     ~   -0.13144093980819422373D+05,
     ~   -0.29206319879254337479D+04,
     ~    0.16213317181975170910D+02,
     ~   -0.19639794571317899851D+03,
     ~    0.25945417987723824460D+06,
     ~    0.48620565446258918740D+07,
     ~    0.10126003837642234550D+07,
     ~    0.84424842477639039373D+04,
     ~    0.18558502462981778080D+03,
     ~   -0.24493795838444620827D+06,
     ~   -0.45896954969516352562D+07,
     ~   -0.95728945603260022778D+06,
     ~   -0.75585692609229706562D+04,
     ~    0.19078884327836625800D+02,
     ~   -0.25327914385527452445D+05,
     ~   -0.47482700183867816434D+06,
     ~   -0.98233311870422994616D+05,
     ~   -0.12493968682442764261D+04/


C	Peak at 6
C	Data (pq (1, j), j = 1, 10) /
C     ~       0.14671022415161132812D+02,
C     ~       0.17176361083984375000D+02,
C     ~       0.22174074172973634588D+02,
C     ~       0.38154026031494140625D+02,
C     ~       0.73457305908203123223D+02,
C     ~       0.11755837202072143554D+01,
C     ~       0.19926652908325195312D+01,
C     ~       0.38097207546234130859D+01,
C     ~       0.89784660339355468750D+01,
C     ~       0.32120651245117186611D+02/

C     MODE1 PEAK at 4 (valid before Aug 2012)
C	Data (pq (1, j), j = 1, 10) /
C     ~    0.67008705139160156250D+01,
C     ~    0.10257417678833007368D+02,
C     ~    0.17875274658203124111D+02,
C     ~    0.36698616027832029473D+02,
C     ~    0.73508399963378909802D+02,
C     ~    0.14384461641311645507D+01,
C     ~    0.27737910747528076171D+01,
C     ~    0.63922395706176757812D+01,
C     ~    0.35446662902832031250D+02,
C     ~    0.73347290039062498223D+02/

C     MODE1 PEAK at 8 strart at 5
C     	Data (pq (1, j), j = 1, 10) /
C     ~   0.39999859333038330078D+01,
C     ~   0.79999923706054687500D+01,
C     ~   0.16000000000000000000D+02,
C     ~   0.32000000000000000000D+02,
C     ~   0.64000000000000000000D+02,
C     ~   0.39984464645385742188D+01,
C     ~   0.80000057220458984375D+01,
C     ~   0.15999995231628417969D+02,
C     ~   0.32000000000000000000D+02,
C     ~   0.64000000000000000000D+02/

    	Data (pq (1, j), j = 1, 10) /
     ~   0.23783802986145019531D+01,
     ~   0.47568554878234863281D+01,
     ~   0.95136852264404296875D+01,
     ~   0.19027336120605468750D+02,
     ~   0.38054645538330078125D+02,
     ~   0.23773388862609863281D+01,
     ~   0.47568039894104003906D+01,
     ~   0.95136604309082031250D+01,
     ~   0.19027313232421875000D+02,
     ~   0.38054626464843750000D+02/

        

C     ################################
C     Aug 2012
C     MODE1 start at 2; PEAK at 6
C     ################################
C	Data (pq (1, j), j = 1, 10) /
C     ~   0.40000061988830566406D+01,
C     ~   0.80000143051147460938D+01,
C     ~   0.16000020980834960938D+02,
C     ~   0.32000015258789062500D+02,
C     ~   0.64000007629394531250D+02,
C     ~   0.39996228218078613281D+01,
C     ~   0.79999980926513671875D+01,
C     ~   0.16000001907348632812D+02,
C     ~   0.32000000000000000000D+02,
C     ~   0.64000000000000000000D+02/
     
C     ################################
C     Aug 2012
C     MODE1 start at 1.2; PEAK at 5
C     ################################
C	Data (pq (1, j), j = 1, 10) /
C     ~   0.40000061988830566406D+01,
C     ~   0.80000143051147460938D+01,
C     ~   0.16000020980834960938D+02,
C     ~   0.32000015258789062500D+02,
C     ~   0.64000007629394531250D+02,
C     ~   0.39996228218078613281D+01,
C     ~   0.79999980926513671875D+01,
C     ~   0.16000001907348632812D+02,
C     ~   0.32000000000000000000D+02,
C     ~   0.64000000000000000000D+02/     

C	MODE 2
C	Data (pq (2, j), j = 1, 10) / 
C     ~      0.62806210517883300781D+01,
C     ~      0.19055501937866210937D+02,
C     ~      0.16761466979980468750D+02,
C     ~      0.28816335678100584161D+02,
C     ~      0.77478836059570310723D+02,
C     ~      0.12525197267532348632D+01,
C     ~      0.21577339172363281250D+01,
C     ~      0.42458128929138183593D+01,
C     ~      0.10150924682617188388D+02,
C     ~      0.32478603363037108486D+02/
     
C	MODE 2
	Data (pq (2, j), j = 1, 10) /     
     ~   0.19998002052307128906D+01,
     ~   0.40004315376281738281D+01,
     ~   0.80004920959472656250D+01,
     ~   0.16000474929809570312D+02,
     ~   0.32000373840332031250D+02,
     ~   0.19961978197097778320D+01,
     ~   0.39991734027862548828D+01,
     ~   0.80000810623168945313D+01,
     ~   0.16000000000000000000D+02,
     ~   0.32000003814697265625D+02/
     
    
    
    
    
    
    
    
    
    
             
     
C	MODE 2
	Data (pq (3, j), j = 1, 10) /         
     ~   0.19997159242630004883D+01,
     ~   0.40001826286315917969D+01,
     ~   0.80002622604370117188D+01,
     ~   0.16000364303588867188D+02,
     ~   0.32000404357910156250D+02,
     ~   0.19948312044143676758D+01,
     ~   0.39987049102783203125D+01,
     ~   0.80001277923583984375D+01,
     ~   0.16000005722045898438D+02,
     ~   0.32000026702880859375D+02/

    
    
         
      
C	MODE 3
C	Data (pq (3, j), j = 1, 10) / 
C     ~     -0.20120578765869141513D+02,
C     ~      0.20586240768432615411D+02,
C     ~      0.23876194000244139736D+02,
C     ~      0.38100570678710937500D+02,
C     ~      0.71865959167480468750D+02,
C     ~      0.11369531154632568359D+01,
C     ~      0.18996990919113159179D+01,
C     ~      0.35767178535461425781D+01,
C     ~      0.82503862380981445312D+01,
C     ~      0.30094644546508790838D+02/

C	MODE 3
	Data (pq (4, j), j = 1, 10) /
     ~  0.24642193317413330078D+01,
     ~  0.49311742782592773438D+01,
     ~  0.98530044555664062500D+01,
     ~  0.19699954986572265625D+02,
     ~  0.39397197723388671875D+02,
     ~  0.23889598846435546875D+01,
     ~  0.49206199645996093750D+01,
     ~  0.98495092391967773438D+01,
     ~  0.19698287963867187500D+02,
     ~  0.39396602630615234375D+02/
    
    
   
    
       
C	MODE 4
C	Data (pq (4, j), j = 1, 10) / 
C     ~      0.38784885406494140625D+02,
C     ~      0.26744506835937498223D+02,
C     ~      0.27683023452758788174D+02,
C     ~      0.38958595275878904473D+02,
C     ~      0.72836189270019531250D+02,
C     ~      0.11952087879180908203D+01,
C     ~      0.20218055248260498046D+01,
C     ~      0.38115730285644531250D+01,
C     ~      0.88158550262451171875D+01,
C     ~      0.33938060760498047763D+02/

C	MODE 4
C	Data (pq (5, j), j = 1, 10) /
C     ~   0.40072717666625976562D+01,
C     ~   0.80101032257080078125D+01,
C     ~   0.16006877899169921875D+02,
C     ~   0.32003021240234375000D+02,
C     ~   0.64001159667968750000D+02,
C     ~   0.39908657073974609375D+01,
C     ~   0.80005998611450195313D+01,
C     ~   0.16000095367431640625D+02,
C     ~   0.32000057220458984375D+02,
C     ~   0.63999942779541015625D+02/
     
C      	Data (pq (5, j), j = 1, 10) /
C     ~   0.40006184577941894531D+01,
C     ~   0.80003767013549804688D+01,
C     ~   0.16000297546386718750D+02,
C     ~   0.32000118255615234375D+02,
C     ~   0.64000045776367187500D+02,
C     ~   0.39639279842376708984D+01,
C     ~   0.80005216598510742188D+01,
C     ~   0.15999864578247070312D+02,
C     ~   0.32000000000000000000D+02,
C     ~   0.63999992370605468750D+02/


         Data (pq (5, j), j = 1, 10) /                           
     ~   0.40074467658996582031D+01,
     ~   0.80055150985717773438D+01,
     ~   0.16002401351928710938D+02,
     ~   0.32000724792480468750D+02,
     ~   0.64000205993652343750D+02,
     ~   0.39702200889587402344D+01,
     ~   0.80015344619750976563D+01,
     ~   0.15999733924865722656D+02,
     ~   0.32000061035156250000D+02,
     ~   0.63999969482421875000D+02/


C     L=32
C	Data (pq (5, j), j = 1, 10) /
C     ~      0.38836956024169921875D+02,
C     ~      0.24311889648437499111D+02,
C     ~      0.29589153289794921875D+02,
C     ~      0.39441062927246095526D+02,
C     ~      0.74735847473144527697D+02,
C     ~      0.13389585018157958984D+01,
C     ~      0.23527803421020507812D+01,
C     ~      0.45696606636047363281D+01,
C     ~      0.10911170959472655894D+02,
C     ~      0.44013973236083980822D+02/

C     L=64
	

C     L=64
	Data (pq (6, j), j = 1, 10) /
     ~     0.36148849487304688388D+02,
     ~     0.26422071456909179687D+02,
     ~     0.52445343017578123223D+02,
     ~     0.50159130096435546875D+02,
     ~     0.72792282104492187500D+02,
     ~     0.14043723344802856445D+01,
     ~     0.26591603755950927734D+01,
     ~     0.55725831985473632812D+01,
     ~     0.15859733581542969638D+02,
     ~     0.70801498413085939276D+02/



c	Separate out the a. b, c, d and p, q.
	KA=1
	Do 10 I=1,M
	Do 10 J=1,M
	a(i,j)=ABCD(L,KA)
	KA=KA+1
   10 continue 
c	Now the PQs
	Do 14 I=1,M
	J=I+M
	p(i)=pq(L,I)/r0or
	p(j)=pq(L,J)/r0or
   14 continue

c	Now compute the magnetic field from differentiating the potential U.

	call gradU(Bxm,Bym,Bzm,x,y,z,a,p,M,M1,M2)	
		
c
    3 continue
      Return
	End
C*************************************************************************
	Subroutine gradU(Bxm1,Bym1,Bzm1,x1,y1,z1,a,p,M,M1,M2)
	IMPLICIT REAL*8(A-H,O-Z)
c
c	Calculates the gradient of the potential.
c
c
	Dimension a(M,M),p(M2)
	M1=M1
	Bxm1=0.
	Bym1=0.
	Bzm1=0.
	Do 111 i=1,M
	Do 222 k=M+1,M2
	J=K-M
c	Write(*,*)M,M1,M2,p(i),p(k)
	Term1=dsqrt(1.D0/p(i)**2+1.D0/p(k)**2)
	Term2=dexp(Term1*x1)
	Bxm1=Bxm1+Term1*a(i,j)*term2*dcos(y1/p(i))*dsin(z1/p(k))
	Bym1=Bym1- a(i,j)/p(i)*term2*dsin(y1/p(i))*dsin(z1/p(k))
	Bzm1=Bzm1+ a(i,j)/p(k)*term2*dcos(y1/p(i))*dcos(z1/p(k))
  222 continue
  111 continue
	Bxm1=-Bxm1
	Bym1=-Bym1
	Bzm1=-Bzm1
	Return
	End
C*********************************************************************
C******************************************************************
      SUBROUTINE CAR2SPH_MAG(BX,BY,BZ,BR,BTH,BPHI,TH,PHI)
      IMPLICIT REAL*8(A-H,O-Z)
! ARFKEN
      BR  = BX*dsin(TH)*dcos(PHI)+BY*dsin(TH)*dsin(PHI)+BZ*dcos(TH)
      BTH = BX*dcos(TH)*dcos(PHI)+BY*dcos(TH)*dsin(PHI)-BZ*dsin(TH)
      BPHI=-BX*dsin(PHI)+BY*dcos(PHI) 
      RETURN
      END
C******************************************************************
      SUBROUTINE SPH2CAR_MAG(BR,BTH,BPHI, BX,BY,BZ,TH,PHI)
      IMPLICIT REAL*8(A-H,O-Z)
! LIUTAURAS RUSAITIS 02.07.2018
! ARFKEN
      BX  = BR*dsin(TH)*dcos(PHI)+BTH*dcos(TH)*dcos(PHI)-BPHI*dsin(PHI)
      BY  = BR*dsin(TH)*dsin(PHI)+BTH*dcos(TH)*dsin(PHI)+BPHI*dcos(PHI)
      BZ  = BR*dcos(TH)-BTH*dsin(TH) 
      RETURN
      END
C******************************************************************
	SUBROUTINE getIMF_penetration(BY_IMF,BZ_IMF,BY_p,BZ_p)
	IMPLICIT REAL*8(A-H,O-Z)
	DATA e1/0.06816/
	DATA e2/0.55417/

	pi=3.14159265358979


	IF ((BY_IMF .eq. 0.0) .AND. (BZ_IMF .eq. 0.0)) THEN
		THETA=0.0
	ELSE
		THETA=DATAN2(BY_IMF,BZ_IMF)
		IF (THETA .lt. 0.0 ) THEN
			THETA = THETA + 2*pi
		ENDIF
	ENDIF

	BY_p=(e1+e2*(dcos(THETA/2.0))**2)*BY_IMF
	BZ_p=(e1+e2*(dcos(THETA/2.0))**2)*BZ_IMF



	RETURN
	END

C******************************************************************************
	SUBROUTINE checkIfInsideMpSATURN(RR,x,NOUT,Dp,a1,a2,dK)
	IMPLICIT REAL*8(A-H,O-Z)



	r0=a1*Dp**(-a2)

	cosTheta=x/RR
	Rmp=r0*(2/(1+cosTheta))**dK

	IF (Rmp .lt. RR) THEN
		NOUT = 0
	ELSE
		NOUT = 1
	ENDIF


	RETURN
	END
c*********************************************************************

C*********************************************************************************
C
      SUBROUTINE KRONIAN(NM,R,T,F,BR,BT,BF)
C
C     BASED ON THE SUBROUTINE IGRF WRITTEN BY N.A. TSYGANENKO (1979)
C     MODIFIED BY KRISHAN KHURANA, JULY, 1996. AND   NOV. 2004.
C
C     CALCULATES COMPONENTS OF MAIN KRONIAN FIELD IN RIGHT HANDED SPHERICAL
C     COORD SYSTEM. BASED ON THE  SPHERICAL HARMONIC COEFFICIENTS GIVEN BY 
C     ACUNA ET AL. [1983] (Z3 MODEL,  V1+V2 DATA) 
C     MAXIMUM ORDER OF HARMONICS TAKEN INTO ACCOUNT (NOT MORE THAN 0RDER 3)
C
C
C     IT IS ASSUMED THAT THE TRAJECTORY IS IS IN RIGHT HANDED S III COORDINATES.
C     THE OUTPUT IS ALSO IN RTP (RH) COORDINATES.
C
C            INPUT:  NM (INTEGER)- MAXIMUM ORDER OF HARMONICS TAKEN
C                                  INTO ACCOUNT (NM.LE.12)
C
C                    R,T,F (REAL)- POSITION OF DESIRED FIELD VALUE IN
C                                  SPHERICAL JOVIGRAPHIC COORDINATE SYSTEM
C                                  (R IN PLANET RADII, COLATITUDE T AND
C                                   RH LONGITUDE F IN RADIANS)
C
C           OUTPUT: BR,BT,BF (REAL)- COMPONENTS OF THE INTERNAL PORTION
C                                    OF THE MAIN MAGNETIC FIELD IN
C                                    SPHERICAL S III COORD SYSTEM
C                                    (VALUES GIVEN IN GAMMA)
C
C***************************************************************************
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BK,BM
      DIMENSION A(13),B(13),G(91),H(91),REC(91)



	DATA G/0.,21144.,24.,1449.,-179.,-35.,2012.,
     ~       53.,92.,-40.,
     ~       81*0.0/

	DATA H/0.,0.,-0.9110,0.,43.5680,-44.7762,0.,
     ~       140.9273,-92.8876,68.2759,
     ~       81*0.0/



      DATA FIRSTI/0.0/
C     WRITE(1,'(5F15.3)')(G(I),I=1,10)
C     WRITE(1,'(5F14.3)')(H(I),I=1,10)
      IF(FIRSTI.EQ.0.0) GO TO 1
      GO TO 6
 1    FIRSTI=1.0
      G(1)=0.
      H(1)=0.
      KNM=15
      DO 2 N=1,13
      N2=2*N-1
      N2=N2*(N2-2)
      DO 2 M=1,N
      MN=N*(N-1)/2+M
  2   REC(MN)=DBLE((N-M)*(N+M-2))/DBLE(N2)
      S=1.
      DO 5 N=2,13
      MN=N*(N-1)/2+1
      S=S*DBLE(2*N-3)/DBLE(N-1)
      G(MN)=G(MN)*S
      H(MN)=H(MN)*S
      P=S
      DO 5 M=2,N
      AA=1.
      IF (M.EQ.2) AA=2.
      P=P*DSQRT(AA*DBLE(N-M+1)/DBLE(N+M-2))
      MNN=MN+M-1
      G(MNN)=G(MNN)*P
  5   H(MNN)=H(MNN)*P
  6   IF(KNM.EQ.NM) GOTO 61
      KNM=NM
      K=KNM+1
 61   PP=1./R
      P=PP
      DO 7 N=1,K
      P=P*PP
      A(N)=P
  7   B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=DCOS(F)
      SF=DSIN(F)
      C=DCOS(U)
      S=DSIN(U)
      BK=(S.LT.1.D-5)
      DO 12 M=1,K
      BM=(M.EQ.1)
      IF(BM) GOTO 8
      MM=M-1
      W=X
      X=W*CF+Y*SF
      Y=Y*CF-W*SF
      GOTO 9
  8   X=0.
      Y=1.
  9   Q=P
      Z=D
      BI=0.
      P2=0.
      D2=0.
      DO 11 N=M,K
      AN=A(N)
      MN=N*(N-1)/2+M
      E=G(MN)
      HH=H(MN)
      W=E*Y+HH*X
      IF (DABS(P2).LT.1.D-38) P2=0.0
      IF (DABS(Q).LT.1.D-38) Q=0.0
      BBR=BBR+B(N)*W*Q
      BBT=BBT-AN*W*Z
      IF(BM) GOTO 10
      QQ=Q
      IF(BK) QQ=Z
      BI=BI+AN*(E*X-HH*Y)*QQ
  10  XK=REC(MN)
      DP=C*Z-S*Q-XK*D2
      PM=C*Q-XK*P2
      D2=Z
      P2=Q
      Z=DP
      Q=PM
   11 CONTINUE
      D=S*D+C*P
      P=S*P
      IF(BM) GOTO 12
      BI=BI*MM
      BBF=BBF+BI
  12  CONTINUE
      BR=BBR
      BT=BBT
      IF(BK) GOTO 13
      BF=BBF/S
      GOTO 14
  13  IF(C.LT.0.) BBF=-BBF
      BF=BBF
  14  CONTINUE
      RETURN
      END
C
C
  
C*************************************************************************


