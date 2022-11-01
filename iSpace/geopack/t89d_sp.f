c
cc
cc  The small main program below is an example of how to compute field
cc   components with T89D_SP.
cc    See GEOPACK-2008.DOC for a sample field line tracing program.
c cc
c       dimension parmod(10)
c C       read*, x,y,z,ps,iopt      
      
c       parmod(1)=3.
c       parmod(2)=-20.
c       parmod(3)=3.
c       parmod(4)=-5.

c       iopt = 3.0
 
c       x = -6.0
c       y = 1.0
c       z = -1.0
c       ps = 1.0
      
 
c       IYEAR=2000
c       IDAY=350
c       IHOUR=21
c       MIN=0
c       ISEC=0
c c      CALL RECALC (IYEAR,IDAY,IHOUR,MIN,ISEC)
 
c       call T89D_SP (iopt,parmod,ps,x,y,z,bx,by,bz)
c       print *, bx,by,bz
c       end
c

C=======================================================================
C
      SUBROUTINE T89D_SP (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
C          (single precision version)
C
C   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
C   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
C   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
C   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
C   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
C
c   NOTE OF NOV 12, 2014:
C
C   THIS IS THIRD UPGRADE OF THE OLD T89 MODEL, IN WHICH ALL COEFFICIENTS
C    (HENCE, THE MODEL OUTPUT) REMAINED IDENTICAL TO THOSE OF T89C.
C   ALL MODIFICATIONS WERE MADE WITH RESPECT TO THE CODE STRUCTURE, WHICH
C    MADE IT POSSIBLE TO SIGNIFICANTLY BOOST THE SPEED.
C--------------------------------------------------------------------------
C -----------   OLD COMMENTS OF 1992:  ------------------------------------
C
C   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
C     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
C
C   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
C   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
C          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
C
C  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
C       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
C         PP.5-20, 1989.
C
C----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
C
C   IOPT= 1       2        3        4        5        6      7
C                  CORRESPOND TO:
C    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
C
C    PS - GEODIPOLE TILT ANGLE IN RADIANS
C    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
C
C----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
C                        FIELD IN NANOTESLAS
c
c   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
C        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
C           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
C
C   THIS RELEASE OF T89C/D IS DATED  FEB 12, 1996;   UPDATED APR 29, 2013
C--------------------------------------------------------------------------
C
C
C              AUTHOR:     NIKOLAI A. TSYGANENKO
C                          HSTX CORP./NASA GSFC (1992-2007)
C                          SPB STATE UNIVERSITY (2007-present)
C
      DIMENSION PARAM(30,7),A(30),PARMOD(10)
      DATA A02,XLW2,YN,RPI,RT/25.,170.,30.,0.31830989E0,30./
      DATA XD,XLD2/0.,40./
C
C   The last 2 quantities define variation of tail sheet thickness along X
C
      DATA SXC,XLWC2/4.,50./
C
C   The two quantities belong to the function WC which confines tail closure
c    current in X- and Y- direction
C
      DATA DXL/20./

      DATA PARAM/-116.53,-10719.,42.375,59.753,-11363.,
     * 1.7844,30.268,-.35372E-01,-0.66832E-01,0.16456E-01,-1.3024,
     * 0.16529E-02,0.20293E-02,20.289,-0.25203E-01,224.91,-9234.8,
     * 22.788,7.8813,1.8362,-0.27228,8.8184,2.8714,14.468,
     * 32.177,0.01,0.0,7.0459,4.0,20.0,
     * -55.553,-13198.,
     * 60.647,61.072,-16064.,2.2534,34.407,-0.38887E-01,
     * -0.94571E-01,0.27154E-01,-1.3901,0.13460E-02,0.13238E-02,
     * 23.005,-0.30565E-01,55.047,-3875.7,20.178,7.9693,
     * 1.4575,0.89471,9.4039,3.5215,14.474,36.555,0.01,
     * 0.0,7.0787,4.0,20.,
     *-101.34,-13480.,111.35,12.386,
     * -24699.,2.6459,38.948,-0.34080E-01,-0.12404,0.29702E-01,
     * -1.4052,0.12103E-02,0.16381E-02,24.49,-0.37705E-01,-298.32,
     * 4400.9,18.692,7.9064,1.3047,2.4541,9.7012,7.1624,
     * 14.288,33.822,0.01,0.0,6.7442,4.0,20.0,
     *-181.69,
     * -12320.,173.79,-96.664,-39051.,3.2633,44.968,
     * -0.46377E-01,-0.16686,0.048298,-1.5473,0.10277E-02,
     * 0.31632E-02,27.341,-0.50655E-01,-514.10,12482.,16.257,
     * 8.5834,1.0194,3.6148,8.6042,5.5057,13.778,32.373,
     * 0.01,0.0,7.3195,4.0,20.0,
     *-436.54,-9001.0,323.66,
     * -410.08,-50340.,3.9932,58.524,-0.38519E-01,-0.26822,
     * 0.74528E-01,-1.4268,-0.10985E-02,0.96613E-02,27.557,
     * -0.56522E-01,-867.03,20652.,14.101,8.3501,0.72996,
     * 3.8149,9.2908,6.4674,13.729,28.353,0.01,0.0,
     * 7.4237,4.0,20.0,-707.77,-4471.9,432.81,-435.51,
     * -60400.,4.6229,68.178,-0.88245E-01,-0.21002,0.11846,
     * -2.6711,0.22305E-02,.10910E-01,27.547,-0.54080E-01,-424.23,
     * 1100.2,13.954,7.5337,0.89714,3.7813,8.2945,5.174,
     * 14.213,25.237,0.01,0.0,7.0037,4.0,20.0,-1190.4,
     * 2749.9,742.56,-1110.3,-77193.,7.6727,102.05,
     *-0.96015E-01,-0.74507,0.11214,-1.3614,0.15157E-02,0.22283E-01,
     *23.164,-0.74146E-01,-2219.1,48253.,12.714,7.6777,.57138,
     * 2.9633,9.3909,9.7263,11.123,21.558,0.01,0.0,
     * 4.4518,4.0,20.0/

       DATA IOP/10/

       SAVE
C
       IF (IOP.NE.IOPT) THEN
C
       IOP=IOPT
       DO 1 I=1,30
   1   A(I)=PARAM(I,IOPT)
C
       DYC=A(30)
       DYC2=DYC**2
       DX=A(18)
       HA02=0.5*A02
       RDX2M=-1./DX**2
       RDX2=-RDX2M
       RDYC2=1./DYC2
       HLWC2M=-0.5*XLWC2
       DRDYC2=-2.*RDYC2
       DRDYC3=2.*RDYC2*SQRT(RDYC2)
       HXLW2M=-0.5*XLW2
       ADR=A(19)
       D0=A(20)
       DD=A(21)
       RC=A(22)
       G=A(23)
       AT=A(24)
       DT=D0
       DEL=A(26)
       P=A(25)
       Q=A(27)
       SX=A(28)
       GAM=A(29)
       HXLD2M=-0.5*XLD2
       ADSL=0.
       XGHS=0.
       H=0.
       HS=0.
       GAMH=0.
       W1=-0.5/DX
       DBLDEL=2.*DEL
       W2=W1*2.
       W4=-1./3.
       W3=W4/DX
       W5=-0.5
       W6=-3.
       AK1=A(1)
       AK2=A(2)
       AK3=A(3)
       AK4=A(4)
       AK5=A(5)
       AK6=A(6)
       AK7=A(7)
       AK8=A(8)
       AK9=A(9)
       AK10=A(10)
       AK11=A(11)
       AK12=A(12)
       AK13=A(13)
       AK14=A(14)
       AK15=A(15)
       AK16=A(16)
       AK17=A(17)
       SXA=0.
       SYA=0.
       SZA=0.
       AK610=AK6*W1+AK10*W5
       AK711=AK7*W2-AK11
       AK812=AK8*W2+AK12*W6
       AK913=AK9*W3+AK13*W4
       RDXL=1./DXL
       HRDXL=0.5*RDXL
       A6H=AK6*0.5
       A9T=AK9/3.
       YNP=RPI/YN*0.5
       YND=2.*YN
C
       ENDIF
C
       SPS = SIN(PS)
       CPS = COS(PS)
C
       X2=X*X
       Y2=Y*Y
       Z2=Z*Z
       TPS=SPS/CPS
       HTP=TPS*0.5
       GSP=G*SPS
       XSM=X*CPS-Z*SPS
       ZSM=X*SPS+Z*CPS
C
C   CALCULATE THE FUNCTION ZS DEFINING THE SHAPE OF THE TAIL CURRENT SHEET
C    AND ITS SPATIAL DERIVATIVES:
C
       XRC=XSM+RC
       XRC16=XRC**2+16.
       SXRC=SQRT(XRC16)
       Y4=Y2*Y2
       Y410=Y4+1.E4
       SY4=SPS/Y410
       GSY4=G*SY4
       ZS1=HTP*(XRC-SXRC)
       DZSX=-ZS1/SXRC
       ZS=ZS1-GSY4*Y4
       D2ZSGY=-SY4/Y410*4.E4*Y2*Y
       DZSY=G*D2ZSGY
C
C   CALCULATE THE COMPONENTS OF THE RING CURRENT CONTRIBUTION:
C
       XSM2=XSM**2
       DSQT=SQRT(XSM2+A02)
       FA0=0.5*(1.+XSM/DSQT)
       DDR=D0+DD*FA0
       DFA0=HA02/DSQT**3
       ZR=ZSM-ZS
       TR=SQRT(ZR**2+DDR**2)
       RTR=1./TR
       RO2=XSM2+Y2
       ADRT=ADR+TR
       ADRT2=ADRT**2
       FK=1./(ADRT2+RO2)
       DSFC=SQRT(FK)
       FC=FK**2*DSFC
       FACXY=3.0*ADRT*FC*RTR
       XZR=XSM*ZR
       YZR=Y*ZR
       DBXDP=FACXY*XZR
       DER25=FACXY*YZR
       XZYZ=XSM*DZSX+Y*DZSY
       FAQ=ZR*XZYZ-DDR*DD*DFA0*XSM
       DBZDP=FC*(2.*ADRT2-RO2)+FACXY*FAQ
       DER15=DBXDP*CPS+DBZDP*SPS
       DER35=DBZDP*CPS-DBXDP*SPS
C
C  CALCULATE THE TAIL CURRENT SHEET CONTRIBUTION:
C
       DELY2=DEL*Y2
       D=DT+DELY2
       IF (ABS(GAM).LT.1.E-6) GOTO 8
       XXD=XSM-XD
       RQD=1./(XXD**2+XLD2)
       RQDS=SQRT(RQD)
       H=0.5*(1.+XXD*RQDS)
       HS=-HXLD2M*RQD*RQDS
       GAMH=GAM*H
       D=D+GAMH
       XGHS=XSM*GAM*HS
       ADSL=-D*XGHS
   8   D2=D**2
       T=SQRT(ZR**2+D2)
       XSMX=XSM-SX
       RDSQ2=1./(XSMX**2+XLW2)
       RDSQ=SQRT(RDSQ2)
       V=0.5*(1.-XSMX*RDSQ)
       DVX=HXLW2M*RDSQ*RDSQ2
       OM=SQRT(SQRT(XSM2+16.)-XSM)
       OMS=-OM/(OM*OM+XSM)*0.5
       RDY=1./(P+Q*OM)
       OMSV=OMS*V
       RDY2=RDY**2
       FY=1./(1.+Y2*RDY2)
       W=V*FY
       YFY1=2.*FY*Y2*RDY2
       FYPR=YFY1*RDY
       FYDY=FYPR*FY
       DWX=DVX*FY+FYDY*Q*OMSV
       YDWY=-V*YFY1*FY
       DDY=DBLDEL*Y
       ATT=AT+T
       S1=SQRT(ATT**2+RO2)
       F5=1./S1
       F7=1./(S1+ATT)
       F1=F5*F7
       F3=F5**3
       F9=ATT*F3
       FS=ZR*XZYZ-D*Y*DDY+ADSL
       XDWX=XSM*DWX+YDWY
       RTT=1./T
       WT=W*RTT
       BRRZ1=WT*F1
       BRRZ2=WT*F3
       DBXC1=BRRZ1*XZR
       DBXC2=BRRZ2*XZR

       TLT2=PS**2

       DER21=BRRZ1*YZR
       DER22=BRRZ2*YZR
       DER216=DER21*TLT2
       DER217=DER22*TLT2
       WTFS=WT*FS
       DBZC1=W*F5+XDWX*F7+WTFS*F1
       DBZC2=W*F9+XDWX*F1+WTFS*F3
       DER11=DBXC1*CPS+DBZC1*SPS
       DER12=DBXC2*CPS+DBZC2*SPS
       DER31=DBZC1*CPS-DBXC1*SPS
       DER32=DBZC2*CPS-DBXC2*SPS
       DER116=DER11*TLT2
       DER117=DER12*TLT2
       DER316=DER31*TLT2
       DER317=DER32*TLT2
C
C  CALCULATE CONTRIBUTION FROM THE CLOSURE CURRENTS
C
       ZPL=Z+RT
       ZMN=Z-RT
       ROGSM2=X2+Y2
       SPL=SQRT(ZPL**2+ROGSM2)
       SMN=SQRT(ZMN**2+ROGSM2)
       XSXC=X-SXC
       RQC2=1./(XSXC**2+XLWC2)
       RQC=SQRT(RQC2)
       FYC=1./(1.+Y2*RDYC2)
       WC=0.5*(1.-XSXC*RQC)*FYC
       DWCX=HLWC2M*RQC2*RQC*FYC
       DWCY=DRDYC2*WC*FYC*Y
       SZRP=1./(SPL+ZPL)
       SZRM=1./(SMN-ZMN)
       XYWC=X*DWCX+Y*DWCY
       WCSP=WC/SPL
       WCSM=WC/SMN
       FXYP=WCSP*SZRP
       FXYM=WCSM*SZRM
       FXPL=X*FXYP
       FXMN=-X*FXYM
       FYPL=Y*FXYP
       FYMN=-Y*FXYM
       FZPL=WCSP+XYWC*SZRP
       FZMN=WCSM+XYWC*SZRM
       DER13=FXPL+FXMN
       DER14=(FXPL-FXMN)*SPS
       DER23=FYPL+FYMN
       DER24=(FYPL-FYMN)*SPS
       DER33=FZPL+FZMN
       DER34=(FZPL-FZMN)*SPS
C
C   NOW CALCULATE CONTRIBUTION FROM CHAPMAN-FERRARO SOURCES + ALL OTHER
C
       EX=EXP(X/DX)
       EC=EX*CPS
       ES=EX*SPS
       ECZ=EC*Z
       ESZ=ES*Z
       ESZY2=ESZ*Y2
       ESZZ2=ESZ*Z2
       ECZ2=ECZ*Z
       ESY=ES*Y
C
C  FINALLY, CALCULATE NET EXTERNAL MAGNETIC FIELD COMPONENTS,
C    BUT FIRST OF ALL THOSE FOR C.-F. FIELD:
C
       SX1=AK6*ECZ+AK7*ES+AK8*ESY*Y+AK9*ESZ*Z
       SY1=AK10*ECZ*Y+AK11*ESY+AK12*ESY*Y2+AK13*ESY*Z2
       SZ1=AK14*EC+AK15*EC*Y2+AK610*ECZ2+AK711*ESZ+AK812
     * *ESZY2+AK913*ESZZ2
       BXCL=AK3*DER13+AK4*DER14
       BYCL=AK3*DER23+AK4*DER24
       BZCL=AK3*DER33+AK4*DER34
       BXT=AK1*DER11+AK2*DER12+BXCL +AK16*DER116+AK17*DER117
       BYT=AK1*DER21+AK2*DER22+BYCL +AK16*DER216+AK17*DER217
       BZT=AK1*DER31+AK2*DER32+BZCL +AK16*DER316+AK17*DER317
       BX=BXT+AK5*DER15+SX1+SXA
       BY=BYT+AK5*DER25+SY1+SYA
       BZ=BZT+AK5*DER35+SZ1+SZA   
C
       RETURN
       END
c
c=======================================================================

