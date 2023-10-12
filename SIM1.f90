!===Main Programme
!====First layer subroutines used: DATSIM, DECAL, PSHOCK, UPDAT, EXO
            
MODULE DataModule
    implicit none
    real (kind=8), dimension(300) :: A, V, BA, EE
    real (kind=8), dimension(600) :: F, Z
    real (kind=8), dimension(300, 300) ::  X, E, ER,ERX, XX
    real (kind=8), dimension(600, 300) :: XB
    real (kind=8), dimension(1500, 300) :: DUMMY
    real (kind=8), dimension(10, 1000) :: ERR, ZERR
    real (kind=8), dimension(2,4) :: TERM_COEF
    real (kind=8) :: B, BR, BSTEP, P, PR, PX, T, TOL, TOLR, VX, ZERO, P03, P25, P30, P50, P75, ONE, TWO
    real (kind=8) :: CALPH, CBETA, CDELTA, CMIU1, CMIU2, COMGY1, COMGY2, COMGK1, COMGK2, COMGC1, COMGC2, CPSI1, CPSI2
    real (kind=8) :: CYVC, CKVC, CLABOR1, CLABOR2, CPI1, CPI2, CEQR, CEQT, CTHETA, CRHO1, CRHO2, CPHIB1, CPHIB2, CRHO3, CETA1, CETA2
    integer, dimension(300) :: JX
    integer, dimension(10) :: KK, JJ
    integer :: I, J, NG, NIT, NRUNS, IT, NV, IV, IY, IC, LPRINT, KVAR,OBR, KSTART, KEND, KEX, IOPTDA, IOPTB,IOPTC,IOPTAC,IOPTT
    integer :: IB, IG, KKI, KKT, II, NEXH, NEXV, K, KG1, IGP, IVB, KKK, ISSZ, KSL, MG1, NAG,LPER,MAXE,MAXX,MAXP
    integer :: KAG, MAG, NDOG, LEXOG, IPAR, IRAND, LAG, NSTART, ITMX, NPER, NSK, IJX, IPER, IDATE, NEGP, IOPT, IDYN
    integer :: ITR, NARG, MAXITR, NXGP, NCG, IDYN1, IRANDX, IRHO, NFIX2, LDOFIX, NFIXMX, NFIXEND, IBLOC, NCONV, NBLOC
    integer, dimension(41) :: KGROUP, JGROUP
    character(len=8) :: NCKEX
    character(len=80) :: line
    character(len=8), dimension(300) :: NAME
    character(len=8), dimension(600) :: IXOG   
    logical :: LFIX2
END MODULE DataModule

program SIM
  USE DataModule

  IMPLICIT NONE


!=====Step 1: Read in the bench data
OPEN(UNIT=14, FILE='nerror', STATUS='old')
  DO j = 1, 147
     READ(14, *) (ERR(i,j), i=1, 10)
  END DO
  CLOSE(14)

  OPEN(UNIT=17, FILE='term_coef', STATUS='old')
  DO i = 1, 2
     READ(17, *) term_coef(i,:)
  END DO
  CLOSE(17)
      NAG=0
      MAXE=300 
      LPER=200 
      MAXX=300                                                          
      MAXP=NAG+LPER 
      READ(5,15) NEGP,NXGP,NPER,LAG,NSTART,ITMX,IB,NCG,IRAND,NSK,IJX,IOPT,IDYN,IDYN1,NBLOC,IRANDX,IOPTDA,NIT,IOPTAC
      NG=NEGP+NXGP                                                      
      READ(5,15)(KK(I+1),I=1,NG)                                        
      READ(5,15)(JJ(I+1),I=1,NCG), IRHO
      JJ(1)=0                                                           
      KK(1)=0                                                           
DO I = 1, NEGP
    KK(I+1) = KK(I+1) + KK(I)
END DO

NDOG = KK(NEGP+1)
KK(NEGP+1) = 0

DO I = 1, NXGP
    KK(NEGP+I+1) = KK(NEGP+I+1) + KK(NEGP+I)
END DO

LEXOG = KK(NEGP+NXGP+1)

DO I = 1, NCG
    JJ(I+1) = JJ(I+1) + JJ(I)
END DO

IPAR = JJ(NCG+1)
                                                        
      READ(5,21) PX,P,TOL,B,BR,PR,T                                     
      IF(IDYN.GT.0) READ(5,23) NARG,MAXITR,BSTEP,TOLR,IOPTC,IOPTB, IOPTT, LDOFIX
   23 FORMAT(I8,I8,F4.0,F8.0,4I4)                                                   
      IF(NSK.NE.0) CALL PSHOCK                                          
      MAG=LAG+1                                                         
      KAG=LAG+NPER 
!=====Step 2: Write the information to the console(output.out)  
      WRITE(6,30)                                                       
   30 FORMAT(1H1,130(1H*)/,59X,'PROGRAM SIMM',/,1X,130(1H*)//)                                                     
      WRITE(6,20)                                                       
      WRITE(6,19) NDOG,LEXOG,NPER,LAG,NSTART,ITMX,IB,IPAR,IRAND,NSK,IJX,IOPT,IDYN,IDYN1,NBLOC,NARG,IRANDX,IOPTDA,NIT,IOPTAC       
      WRITE(6,22) PX,P,TOL,B,BR,PR,T,BSTEP,TOLR,MAXITR,IOPTC,IOPTB, IOPTT, LDOFIX
      IF((NDOG.GT.MAXE).OR.(LEXOG.GT.MAXX)) WRITE(6,10)                 
      IF(LAG.GT.NAG) WRITE(6,11)                                        
      IF(NPER.GT.LPER) WRITE(6,12)                                      
!=====Step 3: Read initial values of endogenous and exogenous variables and 
!=====residuals, and coefficient values by call to DATSIM.          
DO I = 1, MAXE
    DO J = 1, MAXP
        ER(I, J) = 0.0D0
    END DO
END DO
                                       
      KEX=0                                                                  
      CALL DATSIM(IOPTDA,IOPTAC,KSTART,KEND)                            
      IF(IOPTDA.GT.0) CALL EXO(IOPTAC,KEX)                              
      IF(IOPTDA.EQ.2) NSTART = NSTART + KSTART -1                       
!=====Step 4: Solve system of equations by call to DEBCAL.     
      CALL DEBCAL
!=====Step 5: Rolling forecast                                                          
IF (NIT /= 0 .AND. IOPTDA >= 2) THEN
    NRUNS = KEND - KSTART
ELSE
    NRUNS = NIT
END IF

IF (NRUNS > 0 .AND. IDYN > 0) KEX = 1

DO IT = 1, NRUNS
    CALL UPDAT(IOPTDA, IOPTAC, KSTART)
    CALL EXO(IOPTAC, KEX)
    CALL DEBCAL
END DO
                                                         
!=====Step 6: Save the residuals,shocks, endogenous/exogenous data        
OPEN(UNIT=50, FILE='allout.out')

DO i = 4, 147
    WRITE(50, 1111) e(1, i), e(2, i), e(3, i), e(4, i), e(5, i), e(6, i), &
                     e(7, i), e(8, i), e(9, i), e(10, i), e(11, i), e(12, i), &
                     e(13, i), e(14, i), e(15, i), e(16, i), x(17, i), x(18, i), &
                     x(19, i), x(20, i), x(21, i), x(22, i), x(23, i), x(24, i), x(25, i)
END DO

CLOSE(50)

 1111 format(25f13.8)
 1112 format(7f12.8)
      DO 1113 I = 1,NDOG                                               
      WRITE(13,9002)  ( E(I,J),J=1,147 )                               
 1113 CONTINUE                                                          
      DO 1114 I = 1,LEXOG                                              
      WRITE(13,9002)  ( X(I,J),J=1,147 )                               
 1114 CONTINUE                                                          
 9002 FORMAT(4F18.12) 
!=====Step 7: What is this?                                                     
DO IC = 1, 2
    NV = KK(IC+1) - KK(IC)
    
    DO IV = 1, NV
        LPRINT = 0
        KVAR = KK(IC) + IV
        
        DO IY = 1, KAG
            V(IY) = ER(KVAR, IY)
            
            IF (ABS(V(IY)) > 0.000001) THEN
                LPRINT = 1
            END IF
        END DO
        
        IF (LPRINT == 1) THEN
            WRITE(9, 7000) IC, IV
            WRITE(9, 7500) (ER(KVAR, IY), IY = 1, 147)
        END IF
    END DO
END DO



 7000 FORMAT(2I2)                                                       
 7500 FORMAT(4F18.12)                                                    
 8000 CONTINUE                                                          
 1100 CONTINUE                                                          
      STOP                                                              
   10 FORMAT(46H VARIABLE DIMENSIONS EXCEED MAXIMUM DIMENSIONS)         
   11 FORMAT(50H SPECIFIED NUMBER OF LAGS EXCEED ALLOWABLE MAXIMUM)     
   12 FORMAT(52H SPECIFIED FORECAST PERIOD EXCEEDS ALLOWABLE MAXIMUM)   
   15 FORMAT(20I4)                                                      
   21 FORMAT(20F4.0)                                                    
   16 FORMAT(16H END OF FORECAST)                                       
   19 FORMAT(10X,'NDOG   =',I4/10X,'LEXOG  =',I4/10X,'NPER   =',I4/10X,'LAG    =',I4/10X,'NSTART =',I4/ &     
     10X,'ITMX   =',I4/10X,'IB     =',I4/10X,'IPAR   =',I4/10X,'IRAND  =',I4/10X,'NSK    =',I4/ &           
     10X,'IJX    =',I4/10X,'IOPT   =',I4/10X,'IDYN   =',I4/10X,'IDYN1  =',I4/10X,'IBLOC  =',I4/ &                     
     10X,'NARG   =',I5/10X,'IRANDX =',I4/10X,'IOPTDA =',I4/10X,'NIT    =',I4/10X,'IOPTAC =',I4)                        
                                   
   22 FORMAT(10X,'PX     =',F8.3/10X,'P      =',F8.3/10X,'TOL    =',F8.3/10X,'B      =',F8.3/10X,'BR     =',F8.3/ &
     10X,'PR     =',F8.3/10X,'T      =',F8.3/10X,'BSTEP  =',F8.3/10X,'TOLR   =',F8.3/10X,'MAXITR =',I4/ &           
     10X,'IOPTC  =',I4/10X,'IOPTB  =',I4/10X,'IOPTT  =',I4/10X,'LDOFIX =',I4)                               
   20 FORMAT(8X,'INPUT PARAMETERS')                                     

end program
