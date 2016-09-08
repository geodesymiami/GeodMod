      SUBROUTINE  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,                           04610000
     *              AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,                  04620002
     *              UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)  04630002
      IMPLICIT REAL*8 (A-H,O-Z)                                         04640000
      REAL*8   ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, 04650000
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             04660000
C                                                                       04670000
C********************************************************************   04680000
C*****                                                          *****   04690000
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   04700000
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   04710000
C*****                         CODED BY  Y.OKADA ... SEP 1991   *****   04720002
C*****                         REVISED   Y.OKADA ... NOV 1991   *****   04730002
C*****                                                          *****   04740000
C********************************************************************   04750000
C                                                                       04760000
C***** INPUT                                                            04770000
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           04780000
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          04790000
C*****   DEPTH : SOURCE DEPTH                                           04800000
C*****   DIP   : DIP-ANGLE (DEGREE)                                     04810000
C*****   AL1,AL2   : FAULT LENGTH (-STRIKE,+STRIKE)                     04820000
C*****   AW1,AW2   : FAULT WIDTH  ( DOWNDIP, UPDIP)                     04830000
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              04840000
C                                                                       04850000
C***** OUTPUT                                                           04860000
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               04870000
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             04880000
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )04890000
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     04900000
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  04910002
C                                                                       04920000
      COMMON /C0/DUMMY(5),SD,CD ,Dummy2(5)             
      COMMON /C2/XI2,ET2,Q2,R ,Dummy3(20) 
      DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)                   04950000
      DATA  F0/0.D0/                                                    04960000
C-----                                                                  04970000
      IF(Z.GT.0.) WRITE(6,'('' ** POSITIVE Z WAS GIVEN IN SUB-DC3D'')') 04980000
      DO 111 I=1,12                                                     04990000
        U  (I)=F0                                                       05000000
        DUA(I)=F0                                                       05010000
        DUB(I)=F0                                                       05020000
        DUC(I)=F0                                                       05030000
  111 CONTINUE                                                          05040000
      AALPHA=ALPHA                                                      05050000
      DDIP=DIP                                                          05060000
      CALL DCCON0(AALPHA,DDIP)                                          05070000
C======================================                                 05080000
C=====  REAL-SOURCE CONTRIBUTION  =====                                 05090000
C======================================                                 05100000
      D=DEPTH+Z                                                         05110000
      P=Y*CD+D*SD                                                       05120000
      Q=Y*SD-D*CD                                                       05130000
      JXI=0                                                             05140000
      JET=0                                                             05150000
      IF((X+AL1)*(X-AL2).LE.0.) JXI=1                                   05160000
      IF((P+AW1)*(P-AW2).LE.0.) JET=1                                   05170000
      DD1=DISL1                                                         05180000
      DD2=DISL2                                                         05190000
      DD3=DISL3                                                         05200000
C-----                                                                  05210000
      DO 223 K=1,2                                                      05220000
      IF(K.EQ.1) ET=P+AW1                                               05230000
      IF(K.EQ.2) ET=P-AW2                                               05240000
      DO 222 J=1,2                                                      05250000
        IF(J.EQ.1) XI=X+AL1                                             05260000
        IF(J.EQ.2) XI=X-AL2                                             05270000
        CALL DCCON2(XI,ET,Q,SD,CD)                                      05280002
        IF(JXI.EQ.1 .AND. Q.EQ.F0 .AND. ET.EQ.F0) GO TO 99              05290002
        IF(JET.EQ.1 .AND. Q.EQ.F0 .AND. XI.EQ.F0) GO TO 99              05300002
        CALL UA(XI,ET,Q,DD1,DD2,DD3,DUA)                                05310000
C-----                                                                  05320000
        DO 220 I=1,10,3                                                 05330000
          DU(I)  =-DUA(I)                                               05340000
          DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD                              05350000
          DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD                              05360000
          IF(I.LT.10) GO TO 220                                         05370000
          DU(I)  =-DU(I)                                                05380000
          DU(I+1)=-DU(I+1)                                              05390000
          DU(I+2)=-DU(I+2)                                              05400000
  220   CONTINUE                                                        05410000
        DO 221 I=1,12                                                   05420000
          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  05430000
          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  05440000
  221   CONTINUE                                                        05450000
C-----                                                                  05460000
  222 CONTINUE                                                          05470000
  223 CONTINUE                                                          05480000
C=======================================                                05490000
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                05500000
C=======================================                                05510000
      ZZ=Z                                                              05520000
      D=DEPTH-Z                                                         05530000
      P=Y*CD+D*SD                                                       05540000
      Q=Y*SD-D*CD                                                       05550000
      JET=0                                                             05560000
      IF((P+AW1)*(P-AW2).LE.0.) JET=1                                   05570000
C-----                                                                  05580000
      DO 334 K=1,2                                                      05590000
      IF(K.EQ.1) ET=P+AW1                                               05600000
      IF(K.EQ.2) ET=P-AW2                                               05610000
      DO 333 J=1,2                                                      05620000
        IF(J.EQ.1) XI=X+AL1                                             05630000
        IF(J.EQ.2) XI=X-AL2                                             05640000
        CALL DCCON2(XI,ET,Q,SD,CD)                                      05650002
        CALL UA(XI,ET,Q,DD1,DD2,DD3,DUA)                                05660000
        CALL UB(XI,ET,Q,DD1,DD2,DD3,DUB)                                05670000
        CALL UC(XI,ET,Q,ZZ,DD1,DD2,DD3,DUC)                             05680000
C-----                                                                  05690000
        DO 330 I=1,10,3                                                 05700000
          DU(I)=DUA(I)+DUB(I)+Z*DUC(I)                                  05710000
          DU(I+1)=(DUA(I+1)+DUB(I+1)+Z*DUC(I+1))*CD                     05720000
     *           -(DUA(I+2)+DUB(I+2)+Z*DUC(I+2))*SD                     05730000
          DU(I+2)=(DUA(I+1)+DUB(I+1)-Z*DUC(I+1))*SD                     05740000
     *           +(DUA(I+2)+DUB(I+2)-Z*DUC(I+2))*CD                     05750000
          IF(I.LT.10) GO TO 330                                         05760000
          DU(10)=DU(10)+DUC(1)                                          05770000
          DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD                             05780000
          DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD                             05790000
  330   CONTINUE                                                        05800000
        DO 331 I=1,12                                                   05810000
          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  05820000
          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  05830000
  331   CONTINUE                                                        05840000
C-----                                                                  05850000
  333 CONTINUE                                                          05860000
  334 CONTINUE                                                          05870000
C=====                                                                  05880000
      UX=U(1)                                                           05890000
      UY=U(2)                                                           05900000
      UZ=U(3)                                                           05910000
      UXX=U(4)                                                          05920000
      UYX=U(5)                                                          05930000
      UZX=U(6)                                                          05940000
      UXY=U(7)                                                          05950000
      UYY=U(8)                                                          05960000
      UZY=U(9)                                                          05970000
      UXZ=U(10)                                                         05980000
      UYZ=U(11)                                                         05990000
      UZZ=U(12)                                                         06000000
      IRET=0                                                            06010002
      RETURN                                                            06020000
C=======================================                                06030000
C=====  IN CASE OF SINGULAR (R=0)  =====                                06040000
C=======================================                                06050000
   99 UX=F0                                                             06060000
      UY=F0                                                             06070000
      UZ=F0                                                             06080000
      UXX=F0                                                            06090000
      UYX=F0                                                            06100000
      UZX=F0                                                            06110000
      UXY=F0                                                            06120000
      UYY=F0                                                            06130000
      UZY=F0                                                            06140000
      UXZ=F0                                                            06150000
      UYZ=F0                                                            06160000
      UZZ=F0                                                            06170000
      IRET=1                                                            06180002
      RETURN                                                            06190000
      END                                                               06200000



      SUBROUTINE  DCCON0(ALPHA,DIP)                                     09300000
      IMPLICIT REAL*8 (A-H,O-Z)                                         09310000
C                                                                       09320000
C*******************************************************************    09330000
C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    09340000
C*******************************************************************    09350000
C                                                                       09360000
C***** INPUT                                                            09370000
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           09380000
C*****   DIP   : DIP-ANGLE (DEGREE)                                     09390000
C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   09400000
C                                                                       09410000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  09420000
c everything in CO is used in this routine

      DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             09430000
      DATA EPS/1.D-6/                                                   09440000
C-----                                                                  09450000
      ALP1=(F1-ALPHA)/F2                                                09460000
      ALP2= ALPHA/F2                                                    09470000
      ALP3=(F1-ALPHA)/ALPHA                                             09480000
      ALP4= F1-ALPHA                                                    09490000
      ALP5= ALPHA                                                       09500000
C-----                                                                  09510000
      P18=PI2/360.D0                                                    09520000
      SD=DSIN(DIP*P18)                                                  09530000
      CD=DCOS(DIP*P18)                                                  09540000
      IF(DABS(CD).LT.EPS) THEN                                          09550000
        CD=F0                                                           09560000
        IF(SD.GT.F0) SD= F1                                             09570000
        IF(SD.LT.F0) SD=-F1                                             09580000
      ENDIF                                                             09590000
      SDSD=SD*SD                                                        09600000
      CDCD=CD*CD                                                        09610000
      SDCD=SD*CD                                                        09620000
      S2D=F2*SDCD                                                       09630000
      C2D=CDCD-SDSD                                                     09640000
      RETURN                                                            09650000
      END                                                               09660000


      SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U)                       06210000
      IMPLICIT REAL*8 (A-H,O-Z)                                         06220000
      DIMENSION U(12),DU(12)                                            06230000
C                                                                       06240000
C********************************************************************   06250000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   06260000
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   06270000
C********************************************************************   06280000
C                                                                       06290000
C***** INPUT                                                            06300000
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  06310000
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              06320000
C***** OUTPUT                                                           06330000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     06340000
C                                                                       06350000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  06360000
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  06370000
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                06380000
      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/                     06390000
C-----                                                                  06400000
      DO 111  I=1,12                                                    06410000
  111 U(I)=F0                                                           06420000
      XY=XI*Y11                                                         06430000
      QX=Q *X11                                                         06440000
      QY=Q *Y11                                                         06450000
C======================================                                 06460000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 06470000
C======================================                                 06480000
      IF(DISL1.NE.F0) THEN                                              06490000
        DU( 1)=    TT/F2 +ALP2*XI*QY                                    06500000
        DU( 2)=           ALP2*Q/R                                      06510000
        DU( 3)= ALP1*ALE -ALP2*Q*QY                                     06520000
        DU( 4)=-ALP1*QY  -ALP2*XI2*Q*Y32                                06530000
        DU( 5)=          -ALP2*XI*Q/R3                                  06540000
        DU( 6)= ALP1*XY  +ALP2*XI*Q2*Y32                                06550000
        DU( 7)= ALP1*XY*SD        +ALP2*XI*FY+D/F2*X11                  06560000
        DU( 8)=                    ALP2*EY                              06570000
        DU( 9)= ALP1*(CD/R+QY*SD) -ALP2*Q*FY                            06580000
        DU(10)= ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11                  06590000
        DU(11)=                    ALP2*EZ                              06600000
        DU(12)=-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ                            06610000
        DO 222 I=1,12                                                   06620000
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       06630000
      ENDIF                                                             06640000
C======================================                                 06650000
C=====    DIP-SLIP CONTRIBUTION   =====                                 06660000
C======================================                                 06670000
      IF(DISL2.NE.F0) THEN                                              06680000
        DU( 1)=           ALP2*Q/R                                      06690000
        DU( 2)=    TT/F2 +ALP2*ET*QX                                    06700000
        DU( 3)= ALP1*ALX -ALP2*Q*QX                                     06710000
        DU( 4)=        -ALP2*XI*Q/R3                                    06720000
        DU( 5)= -QY/F2 -ALP2*ET*Q/R3                                    06730000
        DU( 6)= ALP1/R +ALP2*Q2/R3                                      06740000
        DU( 7)=                      ALP2*EY                            06750000
        DU( 8)= ALP1*D*X11+XY/F2*SD +ALP2*ET*GY                         06760000
        DU( 9)= ALP1*Y*X11          -ALP2*Q*GY                          06770000
        DU(10)=                      ALP2*EZ                            06780000
        DU(11)= ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ                         06790000
        DU(12)=-ALP1*D*X11          -ALP2*Q*GZ                          06800000
        DO 333 I=1,12                                                   06810000
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       06820000
      ENDIF                                                             06830000
C========================================                               06840000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               06850000
C========================================                               06860000
      IF(DISL3.NE.F0) THEN                                              06870000
        DU( 1)=-ALP1*ALE -ALP2*Q*QY                                     06880000
        DU( 2)=-ALP1*ALX -ALP2*Q*QX                                     06890000
        DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY)                            06900000
        DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32                                06910000
        DU( 5)=-ALP1/R   +ALP2*Q2/R3                                    06920000
        DU( 6)=-ALP1*QY  -ALP2*Q*Q2*Y32                                 06930000
        DU( 7)=-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY                           06940000
        DU( 8)=-ALP1*Y*X11         -ALP2*Q*GY                           06950000
        DU( 9)= ALP1*(D*X11+XY*SD) +ALP2*Q*HY                           06960000
        DU(10)= ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ                           06970000
        DU(11)= ALP1*D*X11         -ALP2*Q*GZ                           06980000
        DU(12)= ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ                           06990000
        DO 444 I=1,12                                                   07000000
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       07010000
      ENDIF                                                             07020000
      RETURN                                                            07030000
      END                                                               07040000


















      SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U)                       07050000
      IMPLICIT REAL*8 (A-H,O-Z)                                         07060000
      DIMENSION U(12),DU(12)                                            07070000
C                                                                       07080000
C********************************************************************   07090000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   07100000
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   07110000
C********************************************************************   07120000
C                                                                       07130000
C***** INPUT                                                            07140000
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  07150000
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              07160000
C***** OUTPUT                                                           07170000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     07180000
C                                                                       07190000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  07200000
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  07210000
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                07220000
      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            07230000
C-----                                                                  07240000
      RD=R+D                                                            07250000
      D11=F1/(R*RD)                                                     07260000
      AJ2=XI*Y/RD*D11                                                   07270000
      AJ5=-(D+Y*Y/RD)*D11                                               07280000
      IF(CD.NE.F0) THEN                                                 07290000
        IF(XI.EQ.F0) THEN                                               07300000
          AI4=F0                                                        07310000
        ELSE                                                            07320000
          X=DSQRT(XI2+Q2)                                               07330000
          AI4=F1/CDCD*( XI/RD*SDCD                                      07340000
     *       +F2*DATAN((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )        07350000
        ENDIF                                                           07360000
        AI3=(Y*CD/RD-ALE+SD*DLOG(RD))/CDCD                              07370000
        AK1=XI*(D11-Y11*SD)/CD                                          07380000
        AK3=(Q*Y11-Y*D11)/CD                                            07390000
        AJ3=(AK1-AJ2*SD)/CD                                             07400000
        AJ6=(AK3-AJ5*SD)/CD                                             07410000
      ELSE                                                              07420000
        RD2=RD*RD                                                       07430000
        AI3=(ET/RD+Y*Q/RD2-ALE)/F2                                      07440000
        AI4=XI*Y/RD2/F2                                                 07450000
        AK1=XI*Q/RD*D11                                                 07460000
        AK3=SD/RD*(XI2*D11-F1)                                          07470000
        AJ3=-XI/RD2*(Q2*D11-F1/F2)                                      07480000
        AJ6=-Y/RD2*(XI2*D11-F1/F2)                                      07490000
      ENDIF                                                             07500000
C-----                                                                  07510000
      XY=XI*Y11                                                         07520000
      AI1=-XI/RD*CD-AI4*SD                                              07530000
      AI2= DLOG(RD)+AI3*SD                                              07540000
      AK2= F1/R+AK3*SD                                                  07550000
      AK4= XY*CD-AK1*SD                                                 07560000
      AJ1= AJ5*CD-AJ6*SD                                                07570000
      AJ4=-XY-AJ2*CD+AJ3*SD                                             07580000
C=====                                                                  07590000
      DO 111  I=1,12                                                    07600000
  111 U(I)=F0                                                           07610000
      QX=Q*X11                                                          07620000
      QY=Q*Y11                                                          07630000
C======================================                                 07640000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 07650000
C======================================                                 07660000
      IF(DISL1.NE.F0) THEN                                              07670000
        DU( 1)=-XI*QY-TT -ALP3*AI1*SD                                   07680000
        DU( 2)=-Q/R      +ALP3*Y/RD*SD                                  07690000
        DU( 3)= Q*QY     -ALP3*AI2*SD                                   07700000
        DU( 4)= XI2*Q*Y32 -ALP3*AJ1*SD                                  07710000
        DU( 5)= XI*Q/R3   -ALP3*AJ2*SD                                  07720000
        DU( 6)=-XI*Q2*Y32 -ALP3*AJ3*SD                                  07730000
        DU( 7)=-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD                           07740000
        DU( 8)=-EY          +ALP3*(F1/R+AJ5)*SD                         07750000
        DU( 9)= Q*FY        -ALP3*(QY-AJ6)*SD                           07760000
        DU(10)=-XI*FZ-Y*X11 +ALP3*AK1*SD                                07770000
        DU(11)=-EZ          +ALP3*Y*D11*SD                              07780000
        DU(12)= Q*FZ        +ALP3*AK2*SD                                07790000
        DO 222 I=1,12                                                   07800000
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       07810000
      ENDIF                                                             07820000
C======================================                                 07830000
C=====    DIP-SLIP CONTRIBUTION   =====                                 07840000
C======================================                                 07850000
      IF(DISL2.NE.F0) THEN                                              07860000
        DU( 1)=-Q/R      +ALP3*AI3*SDCD                                 07870000
        DU( 2)=-ET*QX-TT -ALP3*XI/RD*SDCD                               07880000
        DU( 3)= Q*QX     +ALP3*AI4*SDCD                                 07890000
        DU( 4)= XI*Q/R3     +ALP3*AJ4*SDCD                              07900000
        DU( 5)= ET*Q/R3+QY  +ALP3*AJ5*SDCD                              07910000
        DU( 6)=-Q2/R3       +ALP3*AJ6*SDCD                              07920000
        DU( 7)=-EY          +ALP3*AJ1*SDCD                              07930000
        DU( 8)=-ET*GY-XY*SD +ALP3*AJ2*SDCD                              07940000
        DU( 9)= Q*GY        +ALP3*AJ3*SDCD                              07950000
        DU(10)=-EZ          -ALP3*AK3*SDCD                              07960000
        DU(11)=-ET*GZ-XY*CD -ALP3*XI*D11*SDCD                           07970000
        DU(12)= Q*GZ        -ALP3*AK4*SDCD                              07980000
        DO 333 I=1,12                                                   07990000
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       08000000
      ENDIF                                                             08010000
C========================================                               08020000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               08030000
C========================================                               08040000
      IF(DISL3.NE.F0) THEN                                              08050000
        DU( 1)= Q*QY           -ALP3*AI3*SDSD                           08060000
        DU( 2)= Q*QX           +ALP3*XI/RD*SDSD                         08070000
        DU( 3)= ET*QX+XI*QY-TT -ALP3*AI4*SDSD                           08080000
        DU( 4)=-XI*Q2*Y32 -ALP3*AJ4*SDSD                                08090000
        DU( 5)=-Q2/R3     -ALP3*AJ5*SDSD                                08100000
        DU( 6)= Q*Q2*Y32  -ALP3*AJ6*SDSD                                08110000
        DU( 7)= Q*FY -ALP3*AJ1*SDSD                                     08120000
        DU( 8)= Q*GY -ALP3*AJ2*SDSD                                     08130000
        DU( 9)=-Q*HY -ALP3*AJ3*SDSD                                     08140000
        DU(10)= Q*FZ +ALP3*AK3*SDSD                                     08150000
        DU(11)= Q*GZ +ALP3*XI*D11*SDSD                                  08160000
        DU(12)=-Q*HZ +ALP3*AK4*SDSD                                     08170000
        DO 444 I=1,12                                                   08180000
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       08190000
      ENDIF                                                             08200000
      RETURN                                                            08210000
      END                                                               08220000





      SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U)                     08230000
      IMPLICIT REAL*8 (A-H,O-Z)                                         08240000
      DIMENSION U(12),DU(12)                                            08250000
C                                                                       08260000
C********************************************************************   08270000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   08280000
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   08290000
C********************************************************************   08300000
C                                                                       08310000
C***** INPUT                                                            08320000
C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              08330000
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              08340000
C***** OUTPUT                                                           08350000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     08360000
C                                                                       08370000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  08380000
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  08390000
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                08400000
      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/     08410000
C-----                                                                  08420000
      C=D+Z                                                             08430000
      X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2                     08440000
      Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2                     08450000
      H=Q*CD-Z                                                          08460000
      Z32=SD/R3-H*Y32                                                   08470000
      Z53=F3*SD/R5-H*Y53                                                08480000
      Y0=Y11-XI2*Y32                                                    08490000
      Z0=Z32-XI2*Z53                                                    08500000
      PPY=CD/R3+Q*Y32*SD                                                08510000
      PPZ=SD/R3-Q*Y32*CD                                                08520000
      QQ=Z*Y32+Z32+Z0                                                   08530000
      QQY=F3*C*D/R5-QQ*SD                                               08540000
      QQZ=F3*C*Y/R5-QQ*CD+Q*Y32                                         08550000
      XY=XI*Y11                                                         08560000
      QX=Q*X11                                                          08570000
      QY=Q*Y11                                                          08580000
      QR=F3*Q/R5                                                        08590000
      CQX=C*Q*X53                                                       08600000
      CDR=(C+D)/R3                                                      08610000
      YY0=Y/R3-Y0*CD                                                    08620000
C=====                                                                  08630000
      DO 111  I=1,12                                                    08640000
  111 U(I)=F0                                                           08650000
C======================================                                 08660000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 08670000
C======================================                                 08680000
      IF(DISL1.NE.F0) THEN                                              08690000
        DU( 1)= ALP4*XY*CD           -ALP5*XI*Q*Z32                     08700000
        DU( 2)= ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3                       08710000
        DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)      08720000
        DU( 4)= ALP4*Y0*CD                  -ALP5*Q*Z0                  08730000
        DU( 5)=-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR               08740000
        DU( 6)=-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)    08750000
        DU( 7)=-ALP4*XI*PPY*CD    -ALP5*XI*QQY                          08760000
        DU( 8)= ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD                         08770000
     *                            -ALP5*(CDR*SD-ET/R3-C*Y*QR)           08780000
        DU( 9)=-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD) 08790000
        DU(10)= ALP4*XI*PPZ*CD    -ALP5*XI*QQZ                          08800000
        DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR)   08810000
        DU(12)=         YY0*CD    -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD) 08820000
        DO 222 I=1,12                                                   08830000
  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       08840000
      ENDIF                                                             08850000
C======================================                                 08860000
C=====    DIP-SLIP CONTRIBUTION   =====                                 08870000
C======================================                                 08880000
      IF(DISL2.NE.F0) THEN                                              08890000
        DU( 1)= ALP4*CD/R -QY*SD -ALP5*C*Q/R3                           08900000
        DU( 2)= ALP4*Y*X11       -ALP5*C*ET*Q*X32                       08910000
        DU( 3)=     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32)                   08920000
        DU( 4)=-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD                08930000
        DU( 5)=-ALP4*Y/R3     +ALP5*C*ET*QR                             08940000
        DU( 6)=    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2)                  08950000
        DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR)                08960000
        DU( 8)= ALP4*(X11-Y*Y*X32) -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53) 08970000
        DU( 9)=  XI*PPY*SD+Y*D*X32 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)   08980000
        DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR)                08990000
        DU(11)= ALP4*Y*D*X32       -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53) 09000000
        DU(12)=-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53) 09010000
        DO 333 I=1,12                                                   09020000
  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       09030000
      ENDIF                                                             09040000
C========================================                               09050000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               09060000
C========================================                               09070000
      IF(DISL3.NE.F0) THEN                                              09080000
        DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32)                09090000
        DU( 2)= ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32)                09100000
        DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32)           09110000
        DU( 4)= ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)09120000
        DU( 5)= ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2)             09130000
        DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0)                 09140000
        DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)    09150000
        DU( 8)=-ALP4*F2*XI*PPY*SD-Y*D*X32                               09160000
     *                    +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)            09170000
        DU( 9)=-ALP4*(XI*PPY*CD-X11+Y*Y*X32)                            09180000
     *                    +ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY) 09190000
        DU(10)=  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)  09200000
        DU(11)= ALP4*F2*XI*PPZ*SD-X11+D*D*X32                           09210000
     *                    -ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)            09220000
        DU(12)= ALP4*(XI*PPZ*CD+Y*D*X32)                                09230000
     *                    +ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ) 09240000
        DO 444 I=1,12                                                   09250000
  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       09260000
      ENDIF                                                             09270000
      RETURN                                                            09280000
      END                                                               09290000







      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD)                                 10170002
      IMPLICIT REAL*8 (A-H,O-Z)                                         10180000
C                                                                       10190000
C********************************************************************** 10200000
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 10210000
C********************************************************************** 10220000
C                                                                       10230000
C***** INPUT                                                            10240000
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  10250000
C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                10260000
C                                                                       10270000
C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER010280000
C                                                                       10290000
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  10300000
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                10310000
      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/                          10320000
C-----                                                                  10330000
      IF(DABS(XI).LT.EPS) XI=F0                                         10340000
      IF(DABS(ET).LT.EPS) ET=F0                                         10350000
      IF(DABS( Q).LT.EPS)  Q=F0                                         10360000
      XI2=XI*XI                                                         10370000
      ET2=ET*ET                                                         10380000
      Q2=Q*Q                                                            10390000
      R2=XI2+ET2+Q2                                                     10400000
      R =DSQRT(R2)                                                      10410000
      IF(R.EQ.F0) RETURN                                                10420000
      R3=R *R2                                                          10430000
      R5=R3*R2                                                          10440000
      Y =ET*CD+Q*SD                                                     10450000
      D =ET*SD-Q*CD                                                     10460000
C-----                                                                  10470000
      IF(Q.EQ.F0) THEN                                                  10480000
        TT=F0                                                           10490000
      ELSE                                                              10500000
        TT=DATAN(XI*ET/(Q*R))                                           10510000
      ENDIF                                                             10520000
C-----                                                                  10530000
      IF(XI.LT.F0 .AND. Q.EQ.F0 .AND. ET.EQ.F0) THEN                    10540002
        ALX=-DLOG(R-XI)                                                 10550000
        X11=F0                                                          10560000
        X32=F0                                                          10570000
      ELSE                                                              10580000
        RXI=R+XI                                                        10590002
        ALX=DLOG(RXI)                                                   10600000
        X11=F1/(R*RXI)                                                  10610000
        X32=(R+RXI)*X11*X11/R                                           10620002
      ENDIF                                                             10630000
C-----                                                                  10640000
      IF(ET.LT.F0 .AND. Q.EQ.F0 .AND. XI.EQ.F0) THEN                    10650002
        ALE=-DLOG(R-ET)                                                 10660000
        Y11=F0                                                          10670000
        Y32=F0                                                          10680000
      ELSE                                                              10690000
        RET=R+ET                                                        10700002
        ALE=DLOG(RET)                                                   10710000
        Y11=F1/(R*RET)                                                  10720000
        Y32=(R+RET)*Y11*Y11/R                                           10730002
      ENDIF                                                             10740000
C-----                                                                  10750000
      EY=SD/R-Y*Q/R3                                                    10760000
      EZ=CD/R+D*Q/R3                                                    10770000
      FY=D/R3+XI2*Y32*SD                                                10780000
      FZ=Y/R3+XI2*Y32*CD                                                10790000
      GY=F2*X11*SD-Y*Q*X32                                              10800000
      GZ=F2*X11*CD+D*Q*X32                                              10810000
      HY=D*Q*X32+XI*Q*Y32*SD                                            10820000
      HZ=Y*Q*X32+XI*Q*Y32*CD                                            10830000
      RETURN                                                            10840000
      END                                                               10850000


