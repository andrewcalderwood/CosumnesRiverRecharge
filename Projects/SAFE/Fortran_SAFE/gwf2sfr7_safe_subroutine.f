C-------SUBROUTINE GWF2SFR7SAFE
      SUBROUTINE GWF2SFR7SAFE(Istsg, Wetperm, rchlen, 
     +                        avthk, avhc, IRCH, JRCH, Cstr)
C     ******************************************************************
C     COMPUTE STREAM CONDUCTANCE GIVEN WETTED PERIMTER AND DEPTH USING  
C     8-POINT CROSS SECTION, AQUIFER HYDRAULIC CONDUCTIVITIES, AND
C     PENETRATED AQUIFER PLUS CLOGGING LAYER THICKNESSES
!--------CREATED FOR MODFLOW-2005, MAY 26, 2020
!--------BASED ON SAFE METHOD BY MOREL-SEYTOUX (ET AL.) 2009, 2014, 
!--------2016, 2019 IMPLEMENTED BY CALDERWOOD
C     ******************************************************************
      USE GWFSFRMODULE, ONLY: XSEC
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, 
     +                        HNEW, BOTM, DELR, DELC
      USE GWFLPFMODULE, ONLY: LAYVKA, VKA, HK
      USE GWFHUFMODULE, ONLY: HGUVANI, NHUF, HKHUF=>HK, VKAH=>VKA
      IMPLICIT NONE
      INTRINSIC SQRT, EXP, LOG, DATAN
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER Istsg, IRCH, JRCH
      DOUBLE PRECISION Wetperm, avhc, avthk, Cstr, strlen
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      DOUBLE PRECISION WpN, dpN, Aniso, rho, Daq, maxdepth, maxwidth
      +                stdfar, deltastd, delta, kappa, PI, gflat, 
      +                gflatrect, grectisostd, ganisostd, 
      +                ganisostddelta, ganisostddeltarcl
      REAL aone, atwo
C     ------------------------------------------------------------------
C
CS1-----CALCULATE PARAMETERS USED IN GAMMA (DIMENSIONLESS CONDUCTANCE):  
C       NORMALIZED WETTED PERIMETER AND MAX DEGREE OF PENETRATION, 
C       ANISOTROPY RATIO AND ITS SQUARE ROOT, STANDARD FAR DISTANCE
C
C10-----because BOTM(J,I,0) contains the top elevation of layer 1.
C     BOTM(NCOL, NROW, 0:NBOTM)
      Daq = BOTM(JRCH, IRCH, 0 ) - BOTM(JRCH, IRCH, 1)
C-------CALCULATE THE MAXIMUM DEPTH
      IF( XSEC(9, Istsg) .LE. XSEC(16, Istsg)) THEN
        maxdepth = XSEC(9, Istsg)
      ELSE IF ( XSEC(9, Istsg) .GT. XSEC(16, Istsg) ) THEN
        maxdepth = XSEC(16, Istsg)
      END IF
      maxwidth = XSEC(8, Istsg)
      WpN = Wetperm/Daq
      dpN = maxdepth/Daq
      Aniso = VKA/HK
      rho = SQRT(Aniso)
      stdfar = 2*Daq/rho
      deltastd = stdfar-2*Daq
C     MODFLOW GRID CELL SIZE MAY CHANGE IF DELR!=DELC
      delta = ((DELR/4)-(maxwidth/2))-stdfar
      xi = (1-SQRT(dpN))*(1-rho)
      Rfxi = 1-(0.333*xi)-(0.294*xi*xi)
C
CS2-----CURVE FIT TABLE INFORMATION
      IF (WpN.LT.1 .AND. dpN.LT.0.2) THEN
        aone = 0.89
        atwo = -2.43
      ELSE IF (WpN.LT.1 .AND. dpN.LT.0.5) THEN
        aone = 0.538
        atwo = -0.387
      ELSE IF (WpN.LT.3 .AND. dpN.LT.0.2) THEN
        aone = 0.819
        atwo = -1.34
      ELSE IF (WpN.LT.3 .AND. dpN.LT.0.5) THEN
        aone = 0.672
        atwo = -0.542
      ELSE IF (WpN.LT.3 .AND. dpN.LT.0.9) THEN
        aone = 0.567
        atwo = -0.33
CS3-----DIMENSIONLESS CONDUCTANCE CALCULATIONS
      PI=4.D0*DATAN(1.D0)
C     DIMENSIONLESS CONDUCTANCE ASSUMING FLAT RECHARGE AREA WITH NO
C     PENETRATION
      gflat =1/(2*(1+1/PI*LOG(2/(1-EXP(-1*PI*WpN/2)))))
      IF( WpN.GT.3) THEN
        gflatrect = gflat
C     DIMENSIONLESS CONDUCTANCE ACCOUNTING FOR DEGREE OF PENETRATION
      ELSE gflatrect = gflat*(1+aone*dpN+atwo*dpN*dpN) 
      END IF
C     DIMENSIONLESS CONDUCTANCE WITH CORRECTION FOR STANDARD FAR
C     DISTANCE FROM THE RIVER BANK
      grectisostd = gflatrect/(1+gflatrect*deltastd/Daq)
C     DIMENSIONLESS CONDUCTANCE WITH CORRECTION FOR ANISOTROPY 
      ganisostd = grectisostd*Rfxi
C     DIMENSIONLESS CONDUCTANCE WITH CORRECTION FOR A GRID DISTANCE
C     EXCEEDING THE STANDARD FAR DISTANCE 
      ganisostddelta = ganisostd/(1+ganisostd*delta/Daq)
C     DIMENSIONLESS CONDUCTANCE WITH CORRECTION FOR A CLOGGIN LAYER
      ganisostddeltarcl = ganisostddelta/(1+ganisostddelta*(HK/avhc)
     +                   *(avthk/Wetperm)
C     CALCULATE THE ACTUAL STREAMBED CONDUCTANCE
C     Two times the one side conductance
      Cstr  = 2*HK*rchlen*ganisostddeltarcl
C
CS-----RETURN.
      RETURN
      END SUBROUTINE GWF2SFR7SAFE