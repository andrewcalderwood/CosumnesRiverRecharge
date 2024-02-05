C$LARGE
C       CONNEC3D.FOR            VER.11-NOV-2002
C
C       PROGRAM FOR CONNECTIVITY ANALYSIS OF A 3D INDICATOR MAP
C
C       GIVEN AN INDICATOR FIELDS (VALUES 0 AND 1 ONLY), THE RANDOM SET
C       WITH VALUES 1 IS ANALYSED FOR CONNECTIVITY.
C
C       INPUT PARAMETER FILE WITH
C
C       IPHA : 0 OR 1 FOR CONNECTIVITY ANALYSIS OF PHASE 0 OR 1.
C       ICON : 6, 18 OR 26-CONNECTIVITY.
C       CINP : INPUT FILE WITH INDICATOR VARIABLE (VALUS 0/1 ONLY)
C       NX NY NZ : NUMBER OF POINTS IN X AND Y
C       DX DY DZ : GRID DIMENSIONS IN X AND Y
C       N2    : NUMBER OF LAGS FOR CONNECTIVITY FUNCTION CALCULATION
C       COUT  : OUTPUT FILE WITH STATISTICS
C       COU2  : OUTPUT FILE WITH CONNECTED COMPONENTS
C       COU3  : OUTPUT FILE WITH CONNECTIVITY FUNCITON
C
C       DX,DY,DZ ARE ONLY USED FOR CONVERTING GRID UNITS TO REAL UNITS.
C       NX : NUMBER OF POINTS ALONG THE X DIRECTION
C       NY : NUMBER OF POINTS ALONG THE Y DIRECTION
C       NZ : NUMBER OF POINTS ALONG THE Z DIRECTION
C
C       IND(I,J,L) : 0/1 INDICATOR VALUE
C       AT THE TERMINATION OF THE PROGRAM THE MATRIX IND(I,J) CONTAINS
C       THE CONNECTED COMPONENTS WITH VALUES 1,2,3,4,5,... FOR FIRST,
C       SECOND, THIRD, ... ETC CONNECTED COMMPONENTS.
C       THE VALUE 0 REMAINS 0.
C       COUT : OUTPUT FILE WITH STATISTICS AND CONNECTIVITY FUNCTION.
C       COU2 : OUTPUT FILE WITH THE CONNECTED COMPONENTS.
C       COU3 : OUTPUT FILE WITH CONNECTIVITY FUNCTION
C
C       FOR THE CONNECTIVITY FUNCTION:
C       NPX(K): NUMBER OF PAIRS OF VALUES SEPARATED A DISTANCE K ALONG
C               THE X DIRECTION AND THAT BELONG TO FACIES 1.
C       NPCX(K): NUMBER OF THE PREVIOUS VALUES THAT ARE CONNECTED.
C       NPY(K): NUMBER OF PAIRS OF VALUES SEPARATED A DISTANCE K ALONG
C               THE Y DIRECTION AND THAT BELONG TO FACIES 1.
C       NPCY(K): NUMBER OF THE PREVIOUS VALUES THAT ARE CONNECTED.
C       NPZ(K): NUMBER OF PAIRS OF VALUES SEPARATED A DISTANCE K ALONG
C               THE Z DIRECTION AND THAT BELONG TO FACIES 1.
C       NPCZ(K): NUMBER OF THE PREVIOUS VALUES THAT ARE CONNECTED.
C
C
C       PPHA      : PROPORTION OF FACIES 1.
C       NCC       : NUMBER OF CONNECTED COMPONENTS.
C       RISM      : MEAN CONNECTED COMPONENT SIZE IN PIXELS.
C       RSME      : MRAN SIZE IN REAL UNITS.
C       RTOT      : MEAN SIZE RELATIVE TO SIZE OF FACIES 1.
C       RXME      : MEAN LENGTH ALONG THE X DIRECTION.
C       RYME      : MEAN LENGTH ALONG THE Y DIRECTION.
C       RZME      : MEAN LENGTH ALONG THE Z DIRECTION.
C       ICOM      : NUMBER OF THE LARGEST COMPONENT.
C       ISMA      : SIZE IN PIXELS (OR LARGEST COMPONENT).
C       RSIZ      : SIZE RELATIVE TO SIZE OF FACIES 1.
C       IXMA      : MAXIMUM LENGTH ALONG X.
C       IYMA      : MAXIMUM LENGTH ALONG Y.
C       IZMA      : MAXIMUM LENGTH ALONG Z.
C       ISMI      : SIZE OF SMALLEST COMPONENT.
C       IXMI      : MINIMUM LENGTH ALONG X.
C       IYMI      : MINIMUM LENGTH ALONG Y.
C       IZMI      : MINIMUM LENGTH ALONG Z.
C       IPX       : NUMBER OF PERCOLATING COMPONENTS ALONG X.
C       IPY       : NUMBER OF PERCOLATING COMPONENTS ALONG Y.
C       IPZ       : NUMBER OF PERCOLATING COMPONENTS ALONG Z.
C       CCFU(.,.) : CONNECTIVITY FUNCTION.
C
C        CCFU(.,1): CONNECTIVITY FUNCTION ALONG X (E-W).
C        CCFU(.,2): CONNECTIVITY FUNCTION ALONG Y (N-S).
C        CCFU(.,3): CONNECTIVITY FUNCTION ALONG Z (VERTICAL).
C        CCFU(.,4): CONNECTIVITY FUNCTION ALONG 3D FIRST DIAGONAL.
C        CCFU(.,5): CONNECTIVITY FUNCTION ALONG 3D SECOND DIAGONAL.
C        CCFU(.,6): CONNECTIVITY FUNCTION ALONG 3D THIRD DIAGONAL.
C        CCFU(.,7): CONNECTIVITY FUNCTION ALONG 3D FOURTH DIAGONAL.
C        CCFU(.,8): MEAN CONNECTIVITY FUNCTION ALONG X AND Y.
C        CCFU(.,9): MEAN CONNECTIVITY FUNCTION ALONG X,Y,Z.
C        CCFU(.,10): MEAN CONNECTIVITY FUNCTION ALONG THE 3D DIAGONALS.
C        CCFU(.,11): CONNECTIVITY FUNCTION ALONG FIRST DIAGONAL X-Y PLANE.
C        CCFU(.,12): CONNECTIVITY FUNCTION ALONG SECOND DIAGONAL X-Y PLANE.
C        CCFU(.,13): CONNECTIVITY FUNCTION ALONG FIRST DIAGONAL X-Z PLANE.
C        CCFU(.,14): CONNECTIVITY FUNCTION ALONG SECOND DIAGONAL X-Z PLANE.
C        CCFU(.,15): CONNECTIVITY FUNCTION ALONG FIRST DIAGONAL Y-Z PLANE.
C        CCFU(.,16): CONNECTIVITY FUNCTION ALONG SECOND DIAGONAL Y-Z PLANE.
C
C
C       ----------------------------------------------------------------
C
        DIMENSION IND(500,500,500),NPX(256),NPY(256),NPCX(256),NPCY(256)
        DIMENSION NPZ(256),NPCZ(256),NP1(256),NP2(256),NP3(256)
        DIMENSION NP4(256),NPC1(256),NPC2(256),NPC3(256),NPC4(256)
        DIMENSION CCFU(256,16)
        DIMENSION NXY1(256),NXY2(256),NXZ1(256),NXZ2(256)
        DIMENSION NYZ1(256),NYZ2(256),NCXY1(256),NCXY2(256)
        DIMENSION NCXZ1(256),NCXZ2(256),NCYZ1(256),NCYZ2(256)
C
        CHARACTER*12 CPAR,CINP,COUT,COU2,COU3
C
C       MATRIX DIMENSION LIMITATION
C       Andrew Calderwood increased to 500, compiled with gfortran
C
        NXM=500
        NYM=500
        NZM=500
C
        WRITE (6,400)
        READ (5,100) CPAR
        OPEN (1,FILE=CPAR)
        READ (1,*) IPHA
        IF (IPHA.NE.0.AND.IPHA.NE.1) THEN
          WRITE (6,*)'ERROR IN PARAMETER FILE !!!'
          WRITE (6,*)'THE FIRST LINE IN THE PARAMETER FILE'
          WRITE (6,*)'MUST BE A 0 OR A 1'
          WRITE (6,*)'0 FOR CONNECTIVITY ANALYSIS OF PHASE 0'
          WRITE (6,*)'1 FOR CONNECTIVITY ANALYSIS OF PHASE 1'
          WRITE (6,*)'THE ACTUAL VALUE IN THE PARAMETER FILE IS: ',IPHA
          STOP
        END IF
        READ (1,*) ICON
        IF (ICON.NE.6.AND.ICON.NE.18.AND.ICON.NE.26) THEN
          WRITE (6,*)'ERROR IN PARAMETER FILE !!!'
          WRITE (6,*)'THE SECOND LINE IN THE PARAMETER FILE'
          WRITE (6,*)'MUST BE 6, 18 OR 26'
          WRITE (6,*)'6 MEANS 6-CONNECTIVITY (FACE CONNECTIVITY)'
          WRITE (6,*)'18 MEANS 18-CONNECTIVITY (FACE+EDGE CONNECTIVITY)'
          WRITE (6,*)'26 MEANS 26-CONNECTIVITY (FACE+EDGE+VERTEX CONN.)'
          WRITE (6,*)'THE ACTUAL VALUE IN THE PARAMETER FILE IS: ',ICON
          STOP
        END IF
        READ (1,100) CINP
        READ (1,*) NX,NY,NZ
        READ (1,*) DX,DY,DZ
        READ (1,*) N2
        READ (1,100) COUT
        READ (1,100) COU2
        READ (1,100) COU3
        CLOSE (1)
        IF (NX.GT.NXM) THEN
          WRITE (6,*)'ERROR IN PARAMETER FILE !!!'
          WRITE (6,*)'MAXIMUM ALLOWED NX IS ',NXM
          WRITE (6,*)'ACTUAL VALUE IS ',NX
          STOP
        END IF
        IF (NY.GT.NYM) THEN
          WRITE (6,*)'ERROR IN PARAMETER FILE !!!'
          WRITE (6,*)'MAXIMUM ALLOWED NY IS ',NYM
          WRITE (6,*)'ACTUAL VALUE IS ',NY
          STOP
        END IF
        IF (NZ.GT.NZM) THEN
          WRITE (6,*)'ERROR IN PARAMETER FILE !!!'
          WRITE (6,*)'MAXIMUM ALLOWED NZ IS ',NZM
          WRITE (6,*)'ACTUAL VALUE IS ',NZ
          STOP
        END IF
C
C       READING EXPERIMENTAL INDICATOR DATA 3D FIELD
C
        OPEN (1,FILE=CINP)
        JPHA=0
        DO 1 I=1,NX
          DO 2 J=1,NY
            DO 3 L=1,NZ
            READ (1,*) VAL
            IND(I,J,L)=INT(VAL)
            IF (IPHA.EQ.0) THEN
              IF (INT(VAL).EQ.0) THEN
                IND(I,J,L)=1
              ELSE
                IND(I,J,L)=0
              END IF
            END IF
            IF (IND(I,J,L).NE.1.AND.IND(I,J,L).NE.0) THEN
              WRITE (6,*)'EXPERIMENTAL DATA MUST BE INDICATOR DATA 0/1'
              WRITE (6,*)'ACTUAL VALUE : ',IND(I,J,L)
              STOP
            END IF
            IF (IPHA.EQ.0) THEN
              IF (INT(VAL).EQ.0) THEN
                IND(I,J,L)=1
              ELSE
                IND(I,J,L)=0
              END IF
            END IF
            IF (IND(I,J,L).EQ.1) JPHA=JPHA+1
3           CONTINUE
2         CONTINUE
1       CONTINUE
        CLOSE (1)
        IARE=NX*NY*NZ
        PPHA=FLOAT(JPHA)/IARE
        WRITE (6,*)
        IF (ICON.EQ.6) WRITE (6,*)'6-CONNECTIVITY ANALYSIS'
        IF (ICON.EQ.18) WRITE (6,*)'18-CONNECTIVITY ANALYSIS'
        IF (ICON.EQ.26) WRITE (6,*)'26-CONNECTIVITY ANALYSIS'
        WRITE (6,*)
        WRITE (6,*)'PROPORTION OF PHASE 1 IS: ',PPHA
C
C
C       LOOKING FOR CONNECTED COMPONENTS
C
        NCC=0
C
5       DO 6 I=1,NX
          DO 7 J=1,NY
            DO 67 L=1,NZ
              IF (IND(I,J,L).EQ.1) THEN
                IAI=I
                IAJ=J
                IAL=L
                GOTO 8
              END IF
67          CONTINUE
7         CONTINUE
6       CONTINUE
C       NO MORE CONNECTED COMPONENTS SENDING CONTROL TO LABEL 20
        GOTO 20
C
C       A NEW COMPONEND HAS BEEN FOUND, INCREASE NUMBER OF
C       CONNECTED COMPONENTS BY 1
C
C
8       NCC=NCC+1
        WRITE (6,*)'NCC : ',NCC
C
C       LOOKING FOR COMPONENT NCC WITH LABEL NNC+1
C
        IND(IAI,IAJ,IAL)=NCC+1
        DO 9 I=IAI,NX
          DO 10 J=IAJ,NY
            DO 68 L=IAL,NZ
C
C                 6-CONNECTIVITY
C
              IF ((I-1).GT.0) THEN
                IF (IND(I-1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L)=NCC+1
                END IF
              END IF
              IF (I+1.LE.NX) THEN
                IF (IND(I+1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L)=NCC+1
                END IF
              END IF
              IF (J-1.GT.0) THEN
                IF (IND(I,J-1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L)=NCC+1
                END IF
              END IF
              IF (J+1.LE.NY) THEN
                IF (IND(I,J+1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L)=NCC+1
                END IF
              END IF
              IF (L-1.GT.0) THEN
                IF (IND(I,J,L-1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J,L-1)=NCC+1
                END IF
              END IF
              IF (L+1.LE.NZ) THEN
                IF (IND(I,J,L+1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J,L+1)=NCC+1
                END IF
              END IF
C
C                   18-CONNECTIVITY
C
              IF (ICON.GT.6) THEN
                IF ((J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY) THEN
                  IF (IND(I+1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0) THEN
                  IF (IND(I+1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY) THEN
                  IF (IND(I-1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0) THEN
                  IF (IND(I-1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L)=NCC+1
                  END IF
                END IF
              END IF
C
C                    26-CONNECTIVITY
C
              IF (ICON.GT.18) THEN
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L+1)=NCC+1
                  END IF
                END IF
              END IF
68          CONTINUE
10        CONTINUE
9       CONTINUE
        DO 40 I=IAI,1,-1
          DO 41 J=IAJ,1,-1
            DO 69 L=IAL,1,-1
C
C                 6-CONNECTIVITY
C
              IF ((I-1).GT.0) THEN
                IF (IND(I-1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L)=NCC+1
                END IF
              END IF
              IF (I+1.LE.NX) THEN
                IF (IND(I+1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L)=NCC+1
                END IF
              END IF
              IF (J-1.GT.0) THEN
                IF (IND(I,J-1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L)=NCC+1
                END IF
              END IF
              IF (J+1.LE.NY) THEN
                IF (IND(I,J+1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L)=NCC+1
                END IF
              END IF
              IF (L-1.GT.0) THEN
                IF (IND(I,J,L-1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J,L-1)=NCC+1
                END IF
              END IF
              IF (L+1.LE.NZ) THEN
                IF (IND(I,J,L+1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) IND(I,J,L+1)=NCC+1
                END IF
              END IF
C
C                   18-CONNECTIVITY
C
              IF (ICON.GT.6) THEN
                IF ((J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I,J+1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY) THEN
                  IF (IND(I+1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0) THEN
                  IF (IND(I+1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY) THEN
                  IF (IND(I-1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0) THEN
                  IF (IND(I-1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L)=NCC+1
                  END IF
                END IF
              END IF
C
C                    26-CONNECTIVITY
C
              IF (ICON.GT.18) THEN
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L-1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J-1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I-1,J+1,L+1)=NCC+1
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) IND(I+1,J+1,L+1)=NCC+1
                  END IF
                END IF
              END IF
69          CONTINUE
41        CONTINUE
40      CONTINUE
C
62      NNN=0
        DO 60 I=1,NX
          DO 61 J=1,NY
            DO 70 L=1,NZ
C
C                 6-CONNECTIVITY
C
              IF ((I-1).GT.0) THEN
                IF (IND(I-1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I-1,J,L)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
              IF (I+1.LE.NX) THEN
                IF (IND(I+1,J,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I+1,J,L)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
              IF (J-1.GT.0) THEN
                IF (IND(I,J-1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I,J-1,L)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
              IF (J+1.LE.NY) THEN
                IF (IND(I,J+1,L).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I,J+1,L)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
              IF (L-1.GT.0) THEN
                IF (IND(I,J,L-1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I,J,L-1)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
              IF (L+1.LE.NZ) THEN
                IF (IND(I,J,L+1).EQ.1) THEN
                  IF (IND(I,J,L).EQ.NCC+1) THEN
                    IND(I,J,L+1)=NCC+1
                    NNN=NNN+1
                  END IF
                END IF
              END IF
C
C                   18-CONNECTIVITY
C
              IF (ICON.GT.6) THEN
                IF ((J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I,J-1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I,J+1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I,J-1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I,J+1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY) THEN
                  IF (IND(I+1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J+1,L)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0) THEN
                  IF (IND(I+1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J-1,L)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY) THEN
                  IF (IND(I-1,J+1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J+1,L)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0) THEN
                  IF (IND(I-1,J-1,L).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J-1,L)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
              END IF
C
C                    26-CONNECTIVITY
C
              IF (ICON.GT.18) THEN
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J-1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J-1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J-1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I-1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J+1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L-1).GT.0) THEN
                  IF (IND(I+1,J+1,L-1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J+1,L-1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J-1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J-1).GT.0.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J-1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J-1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I-1).GT.0.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I-1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I-1,J+1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
                IF ((I+1).LE.NX.AND.(J+1).LE.NY.AND.(L+1).LE.NZ) THEN
                  IF (IND(I+1,J+1,L+1).EQ.1) THEN
                    IF (IND(I,J,L).EQ.NCC+1) THEN
                      IND(I+1,J+1,L+1)=NCC+1
                      NNN=NNN+1
                    END IF
                  END IF
                END IF
              END IF
70          CONTINUE
61        CONTINUE
60      CONTINUE
        IF (NNN.GT.0) GOTO 62
        GOTO 5
C
C       THE CONNECTED COMPONENT IS THE COMPONENT WITH VALUES
C       IA(I,J,L)=NCC, WE SUBSTRACT 1 FROM THE LEBEL
C
20      DO 50 I=1,NX
          DO 51 J=1,NY
            DO 71 L=1,NZ
              IF (IND(I,J,L).NE.0) IND(I,J,L)=IND(I,J,L)-1
71          CONTINUE
51        CONTINUE
50      CONTINUE
C
C       STATISTICS ON CONNECTED COMPONENTS
C
        ISME=0
        IXME=0
        IYME=0
        IZME=0
        ISMI=NX*NY*NZ
        ISMA=1
        IXMI=NX
        IXMA=1
        IYMI=NY
        IYMA=1
        IZMI=NZ
        IZMA=1
        IPX=0
        IPY=0
        IPZ=0
        WRITE (6,*)
        WRITE (6,*)'NUMBER OF CONNECTED COMPONENTS: ',NCC
        DO 14 III=1,NCC
          IS=0
          IX1=NX
          IY1=NY
          IZ1=NZ
          IX2=1
          IY2=1
          IZ2=1
          DO 11 I=1,NX
            DO 12 J=1,NY
              DO 72 L=1,NZ
                IF (IND(I,J,L).EQ.III) THEN
                  IS=IS+1
                  IF (I.LT.IX1) IX1=I
                  IF (I.GT.IX2) IX2=I
                  IF (J.LT.IY1) IY1=J
                  IF (J.GT.IY2) IY2=J
                  IF (L.LT.IZ1) IZ1=L
                  IF (L.GT.IZ2) IZ2=L
                END IF
72            CONTINUE
12          CONTINUE
11        CONTINUE
          IXSI=IX2-IX1+1
          IYSI=IY2-IY1+1
          IZSI=IZ2-IZ1+1
          ISME=ISME+IS
          IXME=IXME+IXSI
          IYME=IYME+IYSI
          IZME=IZME+IZSI
          IF (IS.LT.ISMI) ISMI=IS
          IF (IS.GT.ISMA) THEN
            ISMA=IS
            ICOM=III
          END IF
          IF (IXSI.LT.IXMI) IXMI=IXSI
          IF (IXSI.GT.IXMA) IXMA=IXSI
          IF (IYSI.LT.IYMI) IYMI=IYSI
          IF (IYSI.GT.IYMA) IYMA=IYSI
          IF (IZSI.LT.IZMI) IZMI=IZSI
          IF (IZSI.GT.IZMA) IZMA=IZSI
          IF (IXSI.EQ.NX) IPX=IPX+1
          IF (IYSI.EQ.NY) IPY=IPY+1
          IF (IZSI.EQ.NZ) IPZ=IPZ+1
14      CONTINUE
C
C       OUTPUT
C
        OPEN (1,FILE=COUT)
        IF (ICON.EQ.6) WRITE (1,430) CPAR
        IF (ICON.EQ.18) WRITE (1,440) CPAR
        IF (ICON.EQ.26) WRITE (1,450) CPAR
        WRITE (1,110) CINP
        WRITE (1,120) NX
        WRITE (1,130) NY
        WRITE (1,131) NZ
        WRITE (1,140) DX
        WRITE (1,150) DY
        WRITE (1,151) DZ
        WRITE (1,160) IPHA,PPHA
        WRITE (1,170) NCC
        RSME=DX*DY*DZ*ISME/FLOAT(NCC)
        RISM=FLOAT(ISME)/NCC
        RTOT=RISM/JPHA
        WRITE (1,175)
        WRITE (1,180) RISM
        WRITE (1,190) RSME
        WRITE (1,200) IPHA,RTOT
        RXME=FLOAT(IXME)/NCC
        RYME=FLOAT(IYME)/NCC
        RZME=FLOAT(IZME)/NCC
        WRITE (1,210) RXME
        WRITE (1,240) RYME
        WRITE (1,241) RZME
        WRITE (1,177)
        WRITE (1,235) ICOM
        WRITE (1,236) ISMA
        RSIZ=FLOAT(ISMA)/JPHA
        WRITE (1,237) IPHA,RSIZ
        WRITE (1,230) IXMA
        WRITE (1,260) IYMA
        WRITE (1,261) IZMA
        WRITE (1,176)
        WRITE (1,238) ISMI
        WRITE (1,220) IXMI
        WRITE (1,250) IYMI
        WRITE (1,251) IZMI
        WRITE (1,178)
        WRITE (1,270) IPX
        WRITE (1,280) IPY
        WRITE (1,281) IPZ
        WRITE (1,290) COU2
        WRITE (1,295) COU3
C
C       CONNECTIVITY FUNCTION
C       CALCULATED FOR N2 LAGS
C
        DO 595 I=1,N2
          NPX(I)=0
          NPY(I)=0
          NPZ(I)=0
          NP1(I)=0
          NP2(I)=0
          NP3(I)=0
          NP4(I)=0
          NPCX(I)=0
          NPCY(I)=0
          NPCZ(I)=0
          NPC1(I)=0
          NPC2(I)=0
          NPC3(I)=0
          NPC4(I)=0
595     CONTINUE
        DO 30 I=1,NX
          DO 31 J=1,NY
            DO 73 L=1,NZ
              I1=IND(I,J,L)
              IF (I1.EQ.0) GOTO 73
              DO 32 K=1,N2
C
C                ALONG X
C
                IF (I+K.LE.NX) THEN
                  IF (IND(I+K,J,L).NE.0) THEN
                    NPX(K)=NPX(K)+1
                    IF (I1.EQ.IND(I+K,J,L)) NPCX(K)=NPCX(K)+1
                  END IF
                END IF
C
C                ALONG Y
C
                IF (J+K.LE.NY) THEN
                  IF (IND(I,J+K,L).NE.0) THEN
                    NPY(K)=NPY(K)+1
                    IF (I1.EQ.IND(I,J+K,L)) NPCY(K)=NPCY(K)+1
                  END IF
                END IF
C
C                ALONG Z
C
                IF (L+K.LE.NZ) THEN
                  IF (IND(I,J,L+K).NE.0) THEN
                    NPZ(K)=NPZ(K)+1
                    IF (I1.EQ.IND(I,J,L+K)) NPCZ(K)=NPCZ(K)+1
                  END IF
                END IF
C
C                ALONG 3D DIAGONAL 1 (I+, J+, K+)
C
                IF ((J+K.LE.NY).AND.(I+K.LE.NX).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I+K,J+K,L+K).NE.0) THEN
                    NP1(K)=NP1(K)+1
                    IF (I1.EQ.IND(I+K,J+K,L+K)) NPC1(K)=NPC1(K)+1
                  END IF
                END IF
C
C                ALONG 3D DIAGONAL 2 (I+, J-, K+)
C
                IF ((J-K.GT.0).AND.(I+K.LE.NX).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I+K,J-K,L+K).NE.0) THEN
                    NP2(K)=NP2(K)+1
                    IF (I1.EQ.IND(I+K,J-K,L+K)) NPC2(K)=NPC2(K)+1
                  END IF
                END IF
C
C                ALONG 3D DIAGONAL 3 (I-, J-, K+)
C
                IF ((J-K.GT.0).AND.(I-K.GT.0).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I-K,J-K,L+K).NE.0) THEN
                    NP3(K)=NP3(K)+1
                    IF (I1.EQ.IND(I-K,J-K,L+K)) NPC3(K)=NPC3(K)+1
                  END IF
                END IF
C
C                ALONG 3D DIAGONAL 4 (I-, J+, K+)
C
                IF ((J+K.LE.NY).AND.(I-K.GT.0).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I-K,J+K,L+K).NE.0) THEN
                    NP4(K)=NP4(K)+1
                    IF (I1.EQ.IND(I-K,J+K,L+K)) NPC4(K)=NPC4(K)+1
                  END IF
                END IF
C
C               FIRST DIAGONAL IN THE PLANE X-Y
C
                IF ((I+K.LE.NX).AND.(J+K.LE.NY)) THEN
                  IF (IND(I+K,J+K,L).NE.0) THEN
                    NXY1(K)=NXY1(K)+1
                    IF (I1.EQ.IND(I+K,J+K,L)) NCXY1(K)=NCXY1(K)+1
                  END IF
                END IF
C
C               SECOND DIAGONAL IN THE PLANE X-Y
C
                IF ((I+K.LE.NX).AND.(J-K.GT.0)) THEN
                  IF (IND(I+K,J-K,L).NE.0) THEN
                    NXY2(K)=NXY2(K)+1
                    IF (I1.EQ.IND(I+K,J-K,L)) NCXY2(K)=NCXY2(K)+1
                  END IF
                END IF
C
C               FIRST DIAGONAL IN THE PLANE X-Z
C
                IF ((I+K.LE.NX).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I+K,J,L+K).NE.0) THEN
                    NXZ1(K)=NXZ1(K)+1
                    IF (I1.EQ.IND(I+K,J,L+K)) NCXZ1(K)=NCXZ1(K)+1
                  END IF
                END IF
C
C               SECOND DIAGONAL IN THE PLANE X-Z
C
                IF ((I+K.LE.NX).AND.(L-K.GT.0)) THEN
                  IF (IND(I+K,J,L-K).NE.0) THEN
                    NXZ2(K)=NXZ2(K)+1
                    IF (I1.EQ.IND(I+K,J,L-K)) NCXZ2(K)=NCXZ2(K)+1
                  END IF
                END IF
C
C               FIRST DIAGONAL IN THE PLANE Y-Z
C
                IF ((J+K.LE.NY).AND.(L+K.LE.NZ)) THEN
                  IF (IND(I,J+K,L+K).NE.0) THEN
                    NYZ1(K)=NYZ1(K)+1
                    IF (I1.EQ.IND(I,J+K,L+K)) NCYZ1(K)=NCYZ1(K)+1
                  END IF
                END IF
C
C               SECOND DIAGONAL IN THE PLANE Y-Z
C
                IF ((J+K.LE.NY).AND.(L-K.GT.0)) THEN
                  IF (IND(I,J+K,L-K).NE.0) THEN
                    NYZ2(K)=NYZ2(K)+1
                    IF (I1.EQ.IND(I,J+K,L-K)) NCYZ2(K)=NCYZ2(K)+1
                  END IF
                END IF
32            CONTINUE
73          CONTINUE
31        CONTINUE
30      CONTINUE
        WRITE (1,300)
        DO 33 I=1,N2
          CCFU(I,1)=FLOAT(NPCX(I))/AMAX0(1,NPX(I))
          WRITE (1,310) I*DX,CCFU(I,1),NPCX(I),NPX(I)
33      CONTINUE
        WRITE (1,320)
        DO 34 I=1,N2
          CCFU(I,2)=FLOAT(NPCY(I))/AMAX0(1,NPY(I))
          WRITE (1,310) I*DY,CCFU(I,2),NPCY(I),NPY(I)
34      CONTINUE
        WRITE (1,321)
        DO 36 I=1,N2
          CCFU(I,3)=FLOAT(NPCZ(I))/AMAX0(1,NPZ(I))
          WRITE (1,310) I*DZ,CCFU(I,3),NPCZ(I),NPZ(I)
36      CONTINUE
        DXYZ=SQRT(DX*DX+DY*DY+DZ*DZ)
        DXY=SQRT(DX*DX+DY*DY)
        DXZ=SQRT(DX*DX+DZ*DZ)
        DYZ=SQRT(DY*DY+DZ*DZ)
        WRITE (1,322)
        DO 37 I=1,N2
          CCFU(I,4)=FLOAT(NPC1(I))/AMAX0(1,NP1(I))
          WRITE (1,310) I*DXYZ,CCFU(I,4),NPC1(I),NP1(I)
37      CONTINUE
        WRITE (1,323)
        DO 81 I=1,N2
          CCFU(I,5)=FLOAT(NPC2(I))/AMAX0(1,NP2(I))
          WRITE (1,310) I*DXYZ,CCFU(I,5),NPC2(I),NP2(I)
81      CONTINUE
        WRITE (1,324)
        DO 82 I=1,N2
          CCFU(I,6)=FLOAT(NPC3(I))/AMAX0(1,NP3(I))
          WRITE (1,310) I*DXYZ,CCFU(I,6),NPC3(I),NP3(I)
82      CONTINUE
        WRITE (1,325)
        DO 83 I=1,N2
          CCFU(I,7)=FLOAT(NPC4(I))/AMAX0(1,NP4(I))
          WRITE (1,310) I*DXYZ,CCFU(I,7),NPC4(I),NP4(I)
83      CONTINUE
        WRITE (1,330)
        DO 35 I=1,N2
          CCFU(I,8)=(NPCX(I)+NPCY(I))/AMAX0(1,(NPX(I)+NPY(I)))
          WRITE (1,310) I*0.5*(DX+DY),CCFU(I,8),NPCX(I)+NPCY(I)
     *                 ,NPX(I)+NPY(I)
35      CONTINUE
        WRITE (1,328)
        DO 84 I=1,N2
          CCFU(I,9)=(NPCX(I)+NPCY(I)+NPCZ(I))/AMAX0(1,(NPX(I)+NPY(I)
     *     +NPZ(I)))
          WRITE (1,310) I*(DX+DY+DZ)/3.0,CCFU(I,9),NPCX(I)+NPCY(I)
     *                 +NPCZ(I),NPX(I)+NPY(I)+NPZ(I)
84      CONTINUE
        WRITE (1,327)
        DO 38 I=1,N2
        CCFU(I,10)=(NPC1(I)+NPC2(I)+NPC3(I)+NPC4(I))
     *             /AMAX0(1,(NP1(I)+NP2(I)+NP3(I)+NP4(I)))
          WRITE(1,310) I*DXYZ,CCFU(I,10),
     *     NPC1(I)+NPC2(I)+NPC3(I)+NPC4(I),NP1(I)+NP2(I)+NP3(I)+NP4(I)
38      CONTINUE
        WRITE (1,329)
        DO 91 I=1,N2
          CCFU(I,11)=FLOAT(NCXY1(I))/AMAX0(1,NXY1(I))
          WRITE (1,310) I*DXY,CCFU(I,11),NCXY1(I),NXY1(I)
91      CONTINUE
        WRITE (1,331)
        DO 92 I=1,N2
          CCFU(I,12)=FLOAT(NCXY2(I))/AMAX0(1,NXY2(I))
          WRITE (1,310) I*DXY,CCFU(I,12),NCXY2(I),NXY2(I)
92      CONTINUE
        WRITE (1,332)
        DO 93 I=1,N2
          CCFU(I,13)=FLOAT(NCXZ1(I))/AMAX0(1,NXZ1(I))
          WRITE (1,310) I*DXZ,CCFU(I,13),NCXZ1(I),NXZ1(I)
93      CONTINUE
        WRITE (1,333)
        DO 94 I=1,N2
          CCFU(I,14)=FLOAT(NCXZ2(I))/AMAX0(1,NXZ2(I))
          WRITE (1,310) I*DXZ,CCFU(I,14),NCXZ2(I),NXZ2(I)
94      CONTINUE
        WRITE (1,334)
        DO 95 I=1,N2
          CCFU(I,15)=FLOAT(NCYZ1(I))/AMAX0(1,NYZ1(I))
          WRITE (1,310) I*DYZ,CCFU(I,15),NCYZ1(I),NYZ1(I)
95      CONTINUE
        WRITE (1,335)
        DO 96 I=1,N2
          CCFU(I,16)=FLOAT(NCYZ2(I))/AMAX0(1,NYZ2(I))
          WRITE (1,310) I*DYZ,CCFU(I,16),NCYZ2(I),NYZ2(I)
96      CONTINUE
        CLOSE (1)
C
C       OUTPUT FILE WITH CONNECTED COMPONENTS
C
        OPEN (1,FILE=COU2)
C
        DO 15 I=1,NX
          DO 16 J=1,NY
            DO 85 L=1,NZ
              WRITE (1,*) IND(I,J,L)
85          CONTINUE
16        CONTINUE
15      CONTINUE
C
        CLOSE (1)
        CLOSE (2)
        OPEN (1,FILE=COU3)
        DO 17 I=1,N2
          WRITE (1,*) I*DX,CCFU(I,1),NPCX(I),NPX(I)
17      CONTINUE
        WRITE (1,*)
        DO 18 I=1,N2
          WRITE (1,*) I*DY,CCFU(I,2),NPCY(I),NPY(I)
18      CONTINUE
        WRITE (1,*)
        DO 86 I=1,N2
          WRITE (1,*) I*DZ,CCFU(I,3),NPCZ(I),NPZ(I)
86      CONTINUE
        WRITE (1,*)
        DO 19 I=1,N2
          WRITE (1,*) I*DXYZ,CCFU(I,4),NPC1(I),NP1(I)
19      CONTINUE
        WRITE (1,*)
        DO 23 I=1,N2
          WRITE (1,*) I*DXYZ,CCFU(I,5),NPC2(I),NP2(I)
23      CONTINUE
        WRITE (1,*)
        DO 87 I=1,N2
          WRITE (1,*) I*DXYZ,CCFU(I,6),NPC3(I),NP3(I)
87      CONTINUE
        WRITE (1,*)
        DO 88 I=1,N2
          WRITE (1,*) I*DXYZ,CCFU(I,7),NPC4(I),NP4(I)
88      CONTINUE
        WRITE (1,*)
        DO 21 I=1,N2
          WRITE (1,*) I*0.5*(DX+DY),CCFU(I,8),NPCX(I)+NPCY(I)
     *               ,NPX(I)+NPY(I)
21      CONTINUE
        WRITE (1,*)
        DO 89 I=1,N2
          WRITE (1,*) I*(DX+DY+DZ)/3.0,CCFU(I,9),NPCX(I)+NPCY(I)
     *               +NPCZ(I),NPX(I)+NPY(I)+NPZ(I)
89      CONTINUE
        WRITE (1,*)
        DO 22 I=1,N2
          WRITE (1,*) I*DXYZ,CCFU(I,10),NPC1(I)+NPC2(I)+NPC3(I)+NPC4(I)
     *               ,NP1(I)+NP2(I)+NP3(I)+NP4(I)
22      CONTINUE
        WRITE (1,*)
        DO 111 I=1,N2
          CCFU(I,11)=FLOAT(NCXY1(I))/AMAX0(1,NXY1(I))
          WRITE (1,*) I*DXY,CCFU(I,11),NCXY1(I),NXY1(I)
111     CONTINUE
        WRITE (1,*)
        DO 112 I=1,N2
          CCFU(I,12)=FLOAT(NCXY2(I))/AMAX0(1,NXY2(I))
          WRITE (1,*) I*DXY,CCFU(I,12),NCXY2(I),NXY2(I)
112     CONTINUE
        WRITE (1,*)
        DO 113 I=1,N2
          CCFU(I,13)=FLOAT(NCXZ1(I))/AMAX0(1,NXZ1(I))
          WRITE (1,*) I*DXZ,CCFU(I,13),NCXZ1(I),NXZ1(I)
113     CONTINUE
        WRITE (1,*)
        DO 114 I=1,N2
          CCFU(I,14)=FLOAT(NCXZ2(I))/AMAX0(1,NXZ2(I))
          WRITE (1,*) I*DXZ,CCFU(I,14),NCXZ2(I),NXZ2(I)
114     CONTINUE
        WRITE (1,*)
        DO 115 I=1,N2
          CCFU(I,15)=FLOAT(NCYZ1(I))/AMAX0(1,NYZ1(I))
          WRITE (1,*) I*DYZ,CCFU(I,15),NCYZ1(I),NYZ1(I)
115     CONTINUE
        WRITE (1,*)
        DO 116 I=1,N2
          CCFU(I,16)=FLOAT(NCYZ2(I))/AMAX0(1,NYZ2(I))
          WRITE (1,*) I*DYZ,CCFU(I,16),NCYZ2(I),NYZ2(I)
116     CONTINUE
        CLOSE (1)
        WRITE (6,410) COUT
        WRITE (6,420) COU2
        WRITE (6,425) COU3
        STOP
C
100     FORMAT (A)
110     FORMAT (1X,'INPUT INDICATOR FILE : ',A12)
120     FORMAT (1X,'NX : ',I6)
130     FORMAT (1X,'NY : ',I6)
131     FORMAT (1X,'NZ : ',I6)
140     FORMAT (1X,'DX : ',F10.3)
150     FORMAT (1X,'DY : ',F10.3)
151     FORMAT (1X,'DZ : ',F10.3)
160     FORMAT (1X,'PROPORTION OF FACIES ',I1,' : ',F10.4)
170     FORMAT (/1X,'CONNECTED COMPONENTS STATISTICS'//
     *          /1X,'NUMBER OF CONNECTED COMPONENTS : ',I6)
175     FORMAT (/1X,'AVERAGES'/1X,8('='))
176     FORMAT (/1X,'MINIMA'/1X,6('='))
177     FORMAT (/1X,'MAXIMA'/1X,6('='))
178     FORMAT (/1X,'PERCOLATION'/1X,11('='))
180     FORMAT (1X,'MEAN SIZE IN PIXELS : ',F10.4)
190     FORMAT (1X,'MEAN SIZE REAL UNITS : ',F12.4)
200     FORMAT (1X,'MEAN SIZE RELATIVE TO TOTAL AREA OF FACIES ',I1,
     *             ' : ',F10.4)
210     FORMAT (1X,'MEAN LENGTH ALONG X (IN PIXELS) : ',F10.4)
220     FORMAT (1X,'MINIMUM LENGTH ALONG X (IN PIXELS) : ',I6)
230     FORMAT (1X,'MAXIMUM LENGTH ALONG X (IN PIXELS) : ',I6)
235     FORMAT (1X,'THE LARGEST COMPONENT IS NUMBER : ',I6)
236     FORMAT (1X,'WITH MAXIMUM SIZE IN PIXELS : ',I6)
237     FORMAT (1X,'AND RELATIVE TO TOTAL AREA OF FACIES '
     *  ,I1,' : ',F10.4)
238     FORMAT (1X,'SIZE IN PIXELS OF SMALLEST COMPONENT : ',I6)
240     FORMAT (1X,'MEAN LENGTH ALONG Y (IN PIXELS) : ',F10.4)
241     FORMAT (1X,'MEAN LENGTH ALONG Z (IN PIXELS) : ',F10.4)
250     FORMAT (1X,'MIMIMUN LENGTH ALONG Y (IN PIXELS) : ',I6)
260     FORMAT (1X,'MAXIMUM LENGTH ALONG Y (IN PIXELS) : ',I6)
251     FORMAT (1X,'MIMIMUN LENGTH ALONG Z (IN PIXELS) : ',I6)
261     FORMAT (1X,'MAXIMUM LENGTH ALONG Z (IN PIXELS) : ',I6)

270     FORMAT (1X,'NUMBER OF PERCOLATING COMPONENTS IN X : ',I6)
280     FORMAT (1X,'NUMBER OF PERCOLATING COMPONENTS IN Y : ',I6)
281     FORMAT (1X,'NUMBER OF PERCOLATING COMPONENTS IN Z : ',I6)
290     FORMAT (/1X,'OUTPUT FILE WITH CONNECTED COMPONENTS : ',A12)
295     FORMAT (/1X,'OUTPUT FILE WITH CONNECTIVIY FUNCTION : ',A12)
300     FORMAT (/1X,'CONNECTIVITY FUNCTION'/1X,21('=')
     *          //1X,'ALONG THE X DIRECTION (1,0,0)')
310     FORMAT (F10.3,F12.6,2I12)
320     FORMAT (/1X,'ALONG THE Y DIRECTION (0,1,0)')
321     FORMAT (/1X,'ALONG THE Z DIRECTION (0,0,1)')
322     FORMAT (/1X,'ALONG THE FIRST 3D DIAGONAL (1,1,1)')
323     FORMAT (/1X,'ALONG THE SECOND 3D DIAGONAL (1,-1,1)')
324     FORMAT (/1X,'ALONG THE THIRD 3D DIAGONAL (-1,-1,1)')
325     FORMAT (/1X,'ALONG THE FOURTH 3D DIAGONAL (-1,1,1)')
327     FORMAT (/1X,'AVERAGE ALONG THE 3D DIAGONALS')
330     FORMAT (/1X,'AVERAGE ALONG X AND Y')
328     FORMAT (/1X,'AVERAGE ALONG X, Y AND Z')
329     FORMAT (/1X,'FIRST DIAGONAL ON THE X-Y PLANE (1,1,0)')
331     FORMAT (/1X,'SECOND DIAGONAL ON THE X-Y PLANE (1,-1,0')
332     FORMAT (/1X,'FIRST DIAGONAL ON THE X-Z PLANE (1,0,1)')
333     FORMAT (/1X,'SECOND DIAGONAL ON THE X-Z PLANE (1,0,-1)')
334     FORMAT (/1X,'FIRST DIAGONAL ON THE Y-Z PLANE (0,1,1)')
335     FORMAT (/1X,'SECOND DIAGONAL ON THE Y-Z PLANE (0,1,-1)')
400     FORMAT (25(/),1X,'CONNEC3D PROGRAM   VER. 1.0'/1X,27('=')///
     *          1X,'INPUT PARAMETER FILE ---> ',$)
410     FORMAT (/1X,'OUTPUT FILE WITH STATISTICS AND CONNECTIVITY '
     *          'FUNCTION: ',A12)
420     FORMAT (/1X,'OUTPUT FILE WITH CONNECTED COMPONENTS: ',A12)
425     FORMAT (/1X,'OUTPUT FILE WITH CONNECTIVITY FUNCTION: ',A12)
430     FORMAT (1X,'OUTPUT RESULTS OF CONNEC3D'//1X'6-CONNECTIVITY'
     *          ' ANALYSIS'/1X,'PARAMETER FILE WAS : ',A12//)
440     FORMAT (1X,'OUTPUT RESULTS OF CONNEC3D'//1X'18-CONNECTIVITY'
     *          ' ANALYSIS'/1X,'PARAMETER FILE WAS : ',A12//)
450     FORMAT (1X,'OUTPUT RESULTS OF CONNEC3D'//1X'26-CONNECTIVITY'
     *          ' ANALYSIS'/1X,'PARAMETER FILE WAS : ',A12//)
        END
