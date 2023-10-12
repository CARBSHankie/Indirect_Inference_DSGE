C===Main Programme
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                                                      
      LOGICAL       LFIX2                                               
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/BASDAT/ XB(600,300)                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      COMMON/SIMDIM/NAG,LPER,MAXE,MAXX,MAXP                             
      COMMON/FIX2/ LFIX2, NFIX2, LDOFIX, NFIXMX, NFIXEND                
      COMMON/SHK/err(10,1000)
      COMMON/TERMCF/term_coef(2,4)
C=====Step 1: Read in the bench data
      open(unit=14,file='qerror',status='old')
      do 2001 j=1,147 
         READ(14,*)(ERR(I,j),I=1,10)
 2001 continue
      open(17,file='term_coef',status='old')
      do 7001 i=1,2
         read (17,*)term_coef(i,:)	 
 7001 continue
      NAG=0
      MAXE=300 
      LPER=200 
      MAXX=300                                                          
      MAXP=NAG+LPER 
C=====1st row of the bench_data files
      READ(5,15) NEGP,NXGP,NPER,LAG,NSTART,ITMX,IB,NCG,IRAND            
     1,NSK,IJX,IOPT,IDYN,IDYN1,NBLOC,IRANDX,IOPTDA,NIT,IOPTAC
      NG=NEGP+NXGP
C=====2nd row of the bench_data files                                                         
      READ(5,15)(KK(I+1),I=1,NG)                                        
      READ(5,15)(JJ(I+1),I=1,NCG), IRHO
      JJ(1)=0                                                           
      KK(1)=0                                                           
      DO 1 I=1,NEGP                                                     
    1 KK(I+1)=KK(I+1)+KK(I)                                             
      NDOG=KK(NEGP+1)                                                   
      KK(NEGP+1)=0                                                      
      DO 2 I=1,NXGP                                                     
    2 KK(NEGP+I+1)=KK(NEGP+I+1)+KK(NEGP+I)                              
      LEXOG=KK(NEGP+NXGP+1)                                             
      DO 3 I=1,NCG                                                      
    3 JJ(I+1)=JJ(I+1)+JJ(I)                                             
      IPAR=JJ(NCG+1) 
C=====4th row of the bench_data files                                                        
      READ(5,21) PX,P,TOL,B,BR,PR,T    
C=====5th row of the bench_data files                                  
      IF(IDYN.GT.0) READ(5,23) NARG,MAXITR,BSTEP,TOLR,IOPTC,IOPTB,      
     >                         IOPTT, LDOFIX
   23 FORMAT(I8,I8,F4.0,F8.0,4I4)                                                   
      IF(NSK.NE.0) CALL PSHOCK                                          
      MAG=LAG+1                                                         
      KAG=LAG+NPER 
C=====Step 2: Write the information to the console(output.out)  
      WRITE(6,30)                                                       
   30 FORMAT(1H1,130(1H*)/,59X,'PROGRAM SIMM',/,1X,130(1H*)//)                                                     
      WRITE(6,20)                                                       
      WRITE(6,19) NDOG,LEXOG,NPER,LAG,NSTART,ITMX,IB,IPAR,IRAND,        
     1NSK,IJX,IOPT,IDYN,IDYN1,NBLOC,NARG,IRANDX,IOPTDA,NIT,IOPTAC       
      WRITE(6,22) PX,P,TOL,B,BR,PR,T,BSTEP,TOLR,MAXITR,IOPTC,IOPTB,     
     >            IOPTT, LDOFIX
      IF((NDOG.GT.MAXE).OR.(LEXOG.GT.MAXX)) WRITE(6,10)                 
      IF(LAG.GT.NAG) WRITE(6,11)                                        
      IF(NPER.GT.LPER) WRITE(6,12)                                      
C=====Step 3: Read initial values of endogenous and exogenous variables and residuals, and coefficient values by call to DATSIM.          
      DO 4 I=1,MAXE                                                     
      DO 4 J=1,MAXP                                                     
    4 ER(I,J)=0.0D0                                                     
      KEX=0                                                                  
      CALL DATSIM(IOPTDA,IOPTAC,KSTART,KEND)                            
      IF(IOPTDA.GT.0) CALL EXO(IOPTAC,KEX)                              
      IF(IOPTDA.EQ.2) NSTART = NSTART + KSTART -1                       
C=====Step 4: Solve system of equations by call DEBCAL.     
      CALL DEBCAL
C=====Step 5: Rolling forecast                                                          
      IF(NIT.EQ.0.AND.IOPTDA.LT.2) GOTO 100                             
      NRUNS=NIT                                                         
      IF(IOPTDA.EQ.2) NRUNS=KEND-KSTART                                 
      IF (NRUNS.LE.0) GOTO 100                                          
      IF(IDYN.GT.0) KEX=1                                               
      DO 200 IT=1,NRUNS                                                      
      CALL UPDAT(IOPTDA,IOPTAC,KSTART)                                       
      CALL EXO(IOPTAC,KEX)                                                   
      CALL DEBCAL                                                       
  200 CONTINUE                                                          
  100 CONTINUE                                                          
C=====Step 6: Save the residuals,shocks, endogenous/exogenous data        
      open (unit=50, file='allout.out')
      do 102, i=4,147
         write(50,1111) e(1,i),e(2,i),e(3,i),e(4,i),e(5,i),e(6,i),
     >e(7,i),e(8,i),e(9,i),e(10,i),e(11,i),e(12,i),e(13,i),e(14,i),
     >e(15,i),e(16,i),x(17,i),x(18,i),x(19,i),x(20,i),x(21,i),
     >x(22,i),x(23,i),x(24,i),x(25,i)
 102  continue
      close(50)
 1111 format(25f13.8)
 1112 format(7f12.8)
      DO 1113 I = 1,NDOG                                               
      WRITE(13,9002)  ( E(I,J),J=1,147 )                               
 1113 CONTINUE                                                          
      DO 1114 I = 1,LEXOG                                              
      WRITE(13,9002)  ( X(I,J),J=1,147 )                               
 1114 CONTINUE                                                          
 9002 FORMAT(4F18.12) 
C=====Step 7: What is this?                                                     
      DO 1100 IC = 1,2                                                
      NV  =  KK(IC+1) - KK(IC)                                          
      DO 8000 IV = 1,NV                                                 
      LPRINT = 0                                                        
      KVAR = KK(IC) + IV                                                
      DO 6000 IY = 1,KAG                                                
      V  =  ER(KVAR,IY)                                                 
      IF ( ABS(V) .GT. 0.000001 ) THEN                                  
        LPRINT  =  1                                                    
      END IF                                                            
 6000 CONTINUE
      IF ( LPRINT .EQ. 1 ) THEN                                         
      WRITE(9,7000) IC,IV                                               
      WRITE(9,7500) ( ER(KVAR,IY),IY=1,147 )                            
      END IF                                                            
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
   19 FORMAT(10X,'NDOG   =',I4/10X,'LEXOG  =',I4/10X,'NPER   =',I4/     
     110X,'LAG    =',I4/10X,'NSTART =',I4/10X,'ITMX   =',I4/            
     210X,'IB     =',I4/10X,'IPAR   =',I4/10X,'IRAND  =',I4/            
     310X,'NSK    =',I4/10X,'IJX    =',I4/10X,'IOPT   =',I4/            
     410X,'IDYN   =',I4/10X,'IDYN1  =',I4/10X,'IBLOC  =',I4/            
     510X,'NARG   =',I5/10X,'IRANDX =',I4/10X,'IOPTDA =',I4/            
     610X,'NIT    =',I4/10X,'IOPTAC =',I4)                              
   20 FORMAT(8X,'INPUT PARAMETERS')                                     
   22 FORMAT(10X,'PX     =',F8.3/10X,'P      =',F8.3/10X,'TOL    =',F8.3
     1/10X,'B      =',F8.3/10X,'BR     =',F8.3/10X,'PR     =',F8.3/     
     210X,'T      =',F8.3/10X,'BSTEP  =',F8.3/10X,'TOLR   =',F8.3       
     3/10X,'MAXITR =',I4/10X,'IOPTC  =',I4/10X,'IOPTB  =',I4/10X,       
     4     'IOPTT  =',I4/10X,'LDOFIX =',I4)
      END
C===SUBROUTINE EQN                                                                
      SUBROUTINE EQN                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/V(300),Z(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/K01,K02,K03,K04,K05,K06,K07,K08,K09,K10,K11,K12,K13,
     1K14,K15,K16,K17,K18,K19,K20,K21,K22,K23,K24,K25,K26,K27,K28,K29,  
     2K30,K31,K32,K33,K34,K35,K36,K37,K38,K39,K40,K41                   
      COMMON/JGROUP/J01,J02,J03,J04,J05,J06,J07,J08,J09,J10,J11,J12,J13,
     1J14,J15,J16,J17,J18,J19,J20,J21,J22,J23,J24,J25,J26,J27,J28,J29,  
     2J30,J31,J32,J33,J34,J35,J36,J37,J38,J39,J40,J41                   
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      CALL SWUS
      RETURN                                                            
      END
C===SUBROUTINE SUBPXE(II,J,K,ERX,NN)                                                               
      SUBROUTINE SUBPXE(II,J,K,ERX,NN)                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/V(300),Z(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/JGROUP/JJ(41)                                              
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      L=KK(II)+J                                                        
      IF(NN.EQ.2) GOTO 2                                                
    1 X(L,IPER+K)=ERX                                                   
      GOTO 3                                                            
    2 ER(L,IPER+K)=ERX                                                  
    3 RETURN                                                            
      END
C===SUBROUTINE SUBCHX(II,J,K,ERX,NN)                                                                 
      SUBROUTINE SUBCHX(II,J,K,ERX,NN)                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/V(300),Z(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/JGROUP/JJ(41)                                              
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      L=KK(II)+J                                                        
      IF(NN.EQ.2)GOTO 2                                                 
      XR=ERX                                                            
      X(L,IPER+K)=X(L,IPER+K)+B*XR                                      
      GOTO 3                                                            
    2 XR=ERX                                                            
      ER(L,IPER+K)=ER(L,IPER+K)+B*XR                                    
    3 RETURN                                                            
      END
C===SUBROUTINE LIMHVV(II,J,K,RLL,RUL )                                                               
      SUBROUTINE LIMHVV(II,J,K,RLL,RUL )                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/JGROUP/JJ(41)                                              
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      L=KK(II)                                                          
      IF (K.NE.0.OR.II.GT.NEGP.OR.F(L+J).EQ.0.0D0) GOTO 20              
      QVV=  (F(L+J)-E(L+J,IPER))+E(L+J,IPER)                            
      IF(YLOG(QVV/E(L+J,IPER-1)).LT.RUL) GOTO 10                        
      F(L+J)=E(L+J,IPER-1)*XPN(RUL)                                     
   10 IF(YLOG(QVV/E(L+J,IPER-1)).GT.RLL) GOTO 20                        
      F(L+J)=E(L+J,IPER-1)*XPN(RLL)                                     
   20 CONTINUE         
      RETURN                                                            
      END
C===SUBROUTINE LIMVV(II,J,K,RLL,RUL )                                                                 
      SUBROUTINE LIMVV(II,J,K,RLL,RUL )                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/JGROUP/JJ(41)                                              
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      L=KK(II)                                                          
      IF( K.NE.0.OR.II.GT.NEGP.OR.F(L+J).EQ.0.0D0) GOTO 20              
      QVV=B*(F(L+J)-E(L+J,IPER))+E(L+J,IPER)                            
      IF(YLOG(QVV/E(L+J,IPER-1)).LT.RUL) GOTO 10                        
      E(L+J,IPER)=E(L+J,IPER-1)*XPN(RUL)                                
      F(L+J)=E(L+J,IPER)                                                
   10 IF(YLOG(QVV/E(L+J,IPER-1)).GT.RLL) GOTO 20                        
      E(L+J,IPER)=E(L+J,IPER-1)*XPN(RLL)                                
      F(L+J)=E(L+J,IPER)                                                
   20 CONTINUE                                                          
      RETURN                                                            
      END  
C===SUBROUTINE LIMLVV(II,J,K,RLL,RUL )                                                               
      SUBROUTINE LIMLVV(II,J,K,RLL,RUL )                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/JGROUP/JJ(41)                                              
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      L=KK(II)                                                          
      IF ( K.NE.0.OR.II.GT.NEGP.OR.F(L+J).EQ.0.0D0 ) GOTO 20            
      QVV=1.0D0*(F(L+J)-E(L+J,IPER))+E(L+J,IPER)                        
      IF(QVV.LT.RUL) GOTO 10                                            
      E(L+J,IPER)=  RUL                                                 
      F(L+J)=E(L+J,IPER)                                                
   10 IF(QVV.GT.RLL) GOTO 20                                            
      E(L+J,IPER)= RLL                                                  
      F(L+J)=E(L+J,IPER)                                                
   20 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE EXOG(KEX)                                                               
      SUBROUTINE EXOG(KEX)                                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/V(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I   ,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      DO 10 IC = 1 , 9                                                  
   10 CALL EXCO (IC)                                                    
      CALL EXCOM                                                        
      RETURN                                                            
      END
C===SUBROUTINE EXCO ( ICE )                                                                
      SUBROUTINE EXCO ( ICE )                                           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I   ,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      DV(I,J,K)= VV(I,J,K) - VV(I,J,K-1)                                
      DL(I,J,K)= VV(I,J,K) / VV(I,J,K-1)                                
      N = 1+ I-MAG                                                      
      ICX= 20 + ICE                                                     
      EX1 = VV(ICX,8,-1)                                                
      CALL SUBPXE(ICX,1,0,EX1,1)                                        
      EX2 = VV(ICX,2,-1)*(DL(ICX,2,-1)+DL(ICX,2,-2))/2.0D0              
      CALL SUBPXE(ICX,2,0,EX2,1)                                        
      EX4 = VV(ICX,4,-1)                                                
      CALL SUBPXE(ICX,4,0,EX4,1)                                        
      EX5 = 0.0D0                                                       
      CALL SUBPXE(ICX,5,0,EX5,1)                                        
      EX6 = 0.0D0                                                       
      CALL SUBPXE(ICX,6,0,EX6,1)                                        
      EX7 = VV(ICX,8,-1)                                                
      CALL SUBPXE(ICX,7,0,EX7,1)                                        
      EX8 = VV(ICX,8,-1)                                                
      CALL SUBPXE(ICX,8,0,EX8,1)                                        
      EX9 = 0.0D0                                                       
      CALL SUBPXE(ICX,9,0,EX9,1)                                        
      EX10 = 0.0D0                                                      
      CALL SUBPXE(ICX,10,0,EX10,1)                                      
      EX11 = VV(ICE,15,-N)                                              
      CALL SUBPXE(ICX,11,0,EX11,1)                                      
      EX12 = VV(ICE,15,-N)                                              
      CALL SUBPXE(ICX,12,0,EX12,1)                                      
      EX3  = VV(ICX,2,0)*(DV(ICX,4,0) +(VV(ICE,6,-N)/VV(ICX,2,-N))*     
     1 DV(ICX,8,0) + VV(ICX,3,-1)/VV(ICX,2,-1) )                        
      CALL SUBPXE(ICX, 3,0,EX3 ,1)                                      
      RETURN                                                            
      END
C===SUBROUTINE EXCOM                                                                
      SUBROUTINE EXCOM                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I   ,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      DV(I,J,K)= VV(I,J,K) - VV(I,J,K-1)                                
      DL(I,J,K)= VV(I,J,K) / VV(I,J,K-1)                                
      N = I-MAG                                                         
      CALL SUBPXE(33,1,0,0.0D0,1)                                       
      CALL SUBPXE(33,2,0,0.0D0,1)                                       
      CALL SUBPXE(33,3,0,0.0D0,1)                                       
      CALL SUBPXE(33,4,0,0.0D0,1)                                       
      CALL SUBPXE(33,5,0,0.0D0,1)                                       
      CALL SUBPXE(33,6,0,0.0D0,1)                                       
      CALL SUBPXE(33,7,0,0.0D0,1)                                       
      RETURN                                                            
      END
C===SUBROUTINE DATSIM(IOPTDA,IOPTAC,KSTART,KEND)                                                              
      SUBROUTINE DATSIM(IOPTDA,IOPTAC,KSTART,KEND)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      COMMON/ARFUNC/ ARCOEF(11,50,6)                                    
      IF(IDYN.EQ.0.OR.IOPTC.EQ.1) GOTO 165                              
      K=0                                                               
      DO 160 IG=1,NEGP                                                  
      KKI=KK(IG)+1                                                      
      KKT=KK(IG+1)                                                      
      IF(IG.EQ.NEGP) KKT=NDOG                                           
      IF((KKT-KKI).LT.0) GOTO 160                                       
      DO 155  J=KKI,KKT                                                 
      READ(5,102) NAME(J),NCKEX                                         
      IF(NCKEX.LE.0) GOTO 155                                           
      DO 157 II=1,NCKEX                                                 
      READ(5,103) NEXH,NEXV                                             
      K=K+1                                                             
      NEN(K)=J                                                          
      NEX(K)=KK(IG+NEGP)+NEXV                                           
      NENH(K)=NEXH                                                      
  157 CONTINUE                                                          
  155 CONTINUE                                                          
  160 CONTINUE                                                          
      NEXP=K                                                            
      MXEH=0                                                            
      DO 159 J=1,NEXP                                                   
  159 IF(MXEH.LT.NENH(J)) MXEH=NENH(J)                                  
      NGX=NEXP*NPER                                                     
      IF(IOPTB.EQ.0) NGX=NEXP*(NPER+1)                                  
      NARG=NGX                                                          
      GOTO 175                                                          
  165 CONTINUE                                                          
      DO 50 I=1,NDOG                                                    
   50 READ(5,101)  NAME(I)                                              
  175 CONTINUE                                                          
      DO 55 I=1,LEXOG                                                   
   55 READ(5,101)  IXOG(I)                                              
      IF(IOPTC.EQ.1) GOTO 180                                           
      WRITE(6,1005) NARG,NEXP,MXEH,IOPTB                                
      DO 170 K=1,NEXP                                                   
  170 WRITE(6,1006) NEN(K),NAME(NEN(K)),NENH(K),NEX(K),IXOG(NEX(K))     
 1005 FORMAT(/10X,'NARG    =',I4/10X,'NEXP    =',I4/10X,'MXEH    =',    
     1 I4//,'  ENDOG. VAR.   EXP. HORIZON   EXP. VAR. ("EXO. VAR.")'    
     2 ,I4/)                                                            
 1006 FORMAT(/2X,I4,2X,A8,6X,I4,6X,I4,2X,A8)                            
  180 CONTINUE                                                          
      LAG=MAG-1                                                         
      KG1=KAG                                                           
      IF(IOPTDA.EQ.2) GOTO 200                                          
      IF(IOPTDA.EQ.1) KG1=MAG                                           
      IF(IB.EQ.0) GO TO 20                                              
      DO 1 J=1,MAG                                                      
    1 READ(5,*) (E(I,J),I=1,NDOG)                                       
      DO 10 J=1,KG1                                                     
   10 READ(5,*) (X(I,J),I=1,LEXOG)                                      
      GO TO 60                                                          
   20 CONTINUE                                                          
      DO 22 I=1,NDOG                                                    
   22 READ(5,*) (E(I,J),J=1,KAG)                                        
   60 IF(IRAND.EQ.0) GOTO 26                                            
      DO 25 K=1,IRAND                                                   
      READ(5,2) IGP,IVB                                                 
      KKK=KK(IGP)+IVB                                                   
   25 READ(5,*) (ER(KKK,J),J=1,KAG)                                     
   26 CONTINUE                                                          
      IF(IB.NE.0) GOTO 65                                               
      DO 24 I=1,LEXOG                                                   
   24 READ(5,*) (X(I,J),J=1,KG1)                                        
   65 IF(IRANDX.EQ.0) GOTO 199                                          
      DO 27 K=1,IRANDX                                                  
      READ(5,2) IGP,IVB                                                 
      KKK=KK(IGP)+IVB                                                   
   27 READ(5,*) (ERX(KKK,J),J=1,KAG)                                    
    2 FORMAT(2I2)                                                       
      GOTO 199                                                          
  200 CONTINUE                                                          
      READ(5,3) ISSZ,KSTART,KEND                                        
      DO 201 I=1,NDOG                                                   
  201 READ(5,*)(EE(I,J) ,J=1,ISSZ)                                      
      DO 202 I=1,LEXOG                                                  
  202 READ(5,*)(XX(I,J) ,J=1,ISSZ)                                      
      KSL=KSTART-LAG                                                    
      IF(KSL.GT.0) GOTO 203                                             
    3 FORMAT(3I4)                                                       
      KSTART=1-KSL                                                      
      NSTART=NSTART-KSL+1                                               
      KSL=KSTART-LAG                                                    
  203 CONTINUE                                                          
      MG1=MAG-1                                                         
      DO 204 J=1,NDOG                                                   
      DO 204 I=1,MG1                                                    
  204 E(J,I)=EE(J,KSL-1+I)                                              
      MG1=MAG-1                                                         
      IF(IOPTAC.EQ.1) MG1=MG1+1                                         
      DO 205 J=1,LEXOG                                                  
      DO 205 I=1,MG1                                                    
  205 X(J,I)=XX(J,KSL-1+I)                                              
      DO 206 J=1,NDOG                                                   
  206 E(J,MAG)=E(J,MAG-1)                                               
  199 CONTINUE                                                          
      WRITE(6,709)                                                      
  709 FORMAT(1H1,55X,'COEFFICIENT VALUES',/)                            
      OPEN(UNIT=123,FILE='coef.data',STATUS='old')
      READ(123,*) (A(I),I=1,IPAR)                                         
      CLOSE(123)
      OPEN(UNIT=123,FILE='ar_coeffs.data',STATUS='old')
      READ(123,*) (A(I),I=IPAR+1,IPAR+IRHO)                                         
      CLOSE(123)
      WRITE(6,710) (A(M),M=1,IPAR+IRHO)                                    
  710 FORMAT(10X,10F12.4)                                               
      WRITE(6,304)                                                      
  304 FORMAT(1H1)                                                       
  101 FORMAT(A8)                                                        
  102 FORMAT(A8,I4)                                                     
  103 FORMAT(2I4)                                                       
      RETURN                                                            
      END
C===SUBROUTINE ALISON(NDUM,FX)                                                                
      SUBROUTINE ALISON(NDUM,FX)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      KX(I,J)=600*(J-1)+I                                               
      KVV(IG,J,K)=KX(KK(IG)+J,I+K)                                      
      IF(IDYN.EQ.2) CALL CONT(IDYN1)                                    
      I=KAG                                                             
      MXEH1=MXEH+1                                                      
      GOTO (10,20,30) (IOPTT+1)                                         
   10 CONTINUE                                                          
      DO 101 IE=1,NEXP                                                  
      IF(IE.EQ.1) GOTO 100                                              
      IF(NEN(IE-1).EQ.NEN(IE)) GOTO 101                                 
  100 CONTINUE                                                          
      DO 1 K =1,MXEH1                                                   
    1 E(NEN(IE),I+K)=X(NEX(IE),I)                                       
  101 CONTINUE                                                          
      GOTO 40                                                           
   20 CONTINUE                                                          
      DO 202 IE=1,NEXP                                                  
      IF(IE.EQ.1) GOTO 200                                              
      IF(NEN(IE-1).EQ.NEN(IE)) GOTO 202                                 
  200 CONTINUE                                                          
      DO 2 K =1,MXEH1                                                   
    2 E(NEN(IE),I+K)=E(NEN(IE),I)                                       
  202 CONTINUE                                                          
      GOTO 40                                                           
   30 CONTINUE                                                          
      DO 303 IE=1,NEXP                                                  
      IF(IE.EQ.1) GOTO 300                                              
      IF(NEN(IE-1).EQ.NEN(IE)) GOTO 303                                 
  300 CONTINUE                                                          
      IF(E(NEN(IE),I-1).EQ.0.0D0) GOTO 31                               
      DKG=E(NEN(IE),I)/E(NEN(IE),I-1)                                   
      GOTO 32                                                           
   31 DKG=1.0D0                                                         
   32 CONTINUE                                                          
      DO 3 K =1,MXEH1                                                   
    3 E(NEN(IE),I+K)=E(NEN(IE),I+K-1)*DKG                               
  303 CONTINUE                                                          
   40 CONTINUE                                                          
      DO 50 IE=1,NEXP                                                   
   50 X(NEX(IE),I+1)=E(NEN(IE),I+NENH(IE))                              
      J =0                                                              
      DO 4 IE=1,NEXP                                                    
      DO 4 I =MAG,KAG                                                   
      J=J+1                                                             
      L(J)=KX(NEX(IE),I)                                                
      R(J)=E(NEN(IE),I+NENH(IE))-X(NEX(IE),I)                           
    4 CONTINUE                                                          
      I=MAG-1                                                           
      DO 5 IE=1,NEXP                                                    
      J=J+1                                                             
      L(J)=KX(NEX(IE),I)                                                
      R(J)=E(NEN(IE),I+NENH(IE))-X(NEX(IE),I)                           
    5 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE UPDAT(IOPTDA,IOPTAC,KSTART)                                                               
      SUBROUTINE UPDAT(IOPTDA,IOPTAC,KSTART)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      NSTART=NSTART+1                                                   
      MG2=MAG-2                                                         
      IF(MG2.LE.0) GOTO 50                                              
      DO 10 J=1,NDOG                                                    
      DO 10 I=1,MG2                                                     
   10 E(J,I)=E(J,I+1)                                                   
      IF(IOPTAC.EQ.1) MG2=MG2+1                                         
      DO 20 J=1,LEXOG                                                   
      DO 20 I=1,MG2                                                     
   20 X(J,I)=X(J,I+1)                                                   
   50 CONTINUE                                                          
      MG1=MAG-1                                                         
      IF(IOPTDA.EQ.2) GOTO 40                                           
      READ(5,*) (E(J,MG1),J=1,NDOG)                                     
      IF(IOPTAC.EQ.1) MG1=MG1+1                                         
      READ(5,*) (X(J,MG1),J=1,LEXOG)                                    
      GOTO 45                                                           
   40 CONTINUE                                                          
      DO 23 J=1,NDOG                                                    
   23 E(J,MG1)=EE(J,KSTART)                                             
      KSST=KSTART                                                       
      IF(IOPTAC.EQ.1) MG1=MG1+1                                         
      IF(IOPTAC.EQ.1) KSST=KSTART+1                                     
      DO 24 J=1,LEXOG                                                   
   24 X(J,MG1)=XX(J,KSST)                                               
      KSTART=KSTART+1                                                   
   45 CONTINUE                                                          
      DO 25 J=1,NDOG                                                    
   25 E(J,MAG)=E(J,MAG+1)                                               
      RETURN                                                            
      END 
C===SUBROUTINE EXO(IOPTAC,KEX)                                                                
      SUBROUTINE EXO(IOPTAC,KEX)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/f(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      MG=MAG                                                            
      IF(IOPTAC.EQ.1) MG=MG+1                                           
      IF (IOPTC.GT.0.OR.KEX.EQ.0) GOTO 10                               
      DO 30 J=1,NEXP                                                    
   30 X(NEX(J),MG-1)=X(NEX(J),MG)                                       
   10 CONTINUE                                                          
      DO 20 I=MG,KAG                                                    
      IPER=I                                                            
      CALL EXOG(KEX)                                                    
      IF (IOPTC.GT.0.OR.KEX.EQ.0) GOTO 20                               
      DO 40 J=1,NEXP                                                    
   40 X(NEX(J),I)=X(NEX(J),I+1)                                         
   20 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE DEBCAL              
      SUBROUTINE DEBCAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG, AJX(10)                                 
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300)  ,X(600,300),ER(300,300),
     &          ERX(600,300),EE(300,300),XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      COMMON/FIX2/ LFIX2, NFIX2, LDOFIX, NFIXMX, NFIXEND                
      DIMENSION NJX(10)                                                 
      DIMENSION XEQUIV(120000)                                           
      EQUIVALENCE ( XEQUIV(1),X(1,1) )                                  
      EXTERNAL CALCFX,ALISON                                            
      COMMON/SHK/err(10,1000)
      character*18 errname, gen_errname, xname, gen_xname
  199 IF(NSK.NE.0) CALL SHOCK                                           
      WRITE(6,226)                                                      
  226 FORMAT(1H1,60X,'INPUT DATA')                                      
      DO 225 K=1,NEGP                                                   
      KKK=KK(K+1)                                                       
      IF(K.EQ.NEGP) KKK=NDOG                                            
      KKK=KKK-KK(K)                                                     
      IF(KKK.EQ.0) GOTO 225                                             
      ISTART=NSTART-MAG                                                 
      WRITE(6,409) K                                                    
      DO 222 I=1,MAG,10                                                 
      DO 223 J=1,10                                                     
  223 JX(J)=ISTART+J                                                    
      IE=MIN0(10,MAG-I+1)                                               
      WRITE(6,408)(JX(J),J=1,IE)                                        
      DO 224 J=1,KKK                                                    
      KJ=KK(K)+J                                                        
  224 WRITE(6,405) NAME(KJ),(E(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+10                                                  
  222 CONTINUE                                                          
  225 CONTINUE                                                          
  405 FORMAT(2X,A8,10F12.4)                                             
  408 FORMAT(/14X,I4,9(8X,I4))                                          
  409 FORMAT(//,2X,'GROUP',I3,' ENDOGENOUS VARIABLES')                  
      DO 525 K=1,NXGP                                                   
      KKK=KK(NEGP+K+1)-KK(NEGP+K)                                       
      IF(KKK.EQ.0) GOTO 525                                             
      K2=K+NEGP                                                         
      WRITE(6,509) K2                                                   
  509 FORMAT(//,2X,'GROUP',I3,' EXOGENOUS VARIABLES')                   
      ISTART=NSTART-MAG                                                 
      DO 522 I=1,KAG,10                                                 
      DO 523 J=1,10                                                     
  523 JX(J)=ISTART+J                                                    
      IE=MIN0(10,KAG-I+1)                                               
      WRITE(6,408)(JX(J),J=1,IE)                                        
      DO 524 J=1,KKK                                                    
      KJ=KK(NEGP+K)+J                                                   
  524 WRITE(6,405) IXOG(KJ),(X(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+10                                                  
  522 CONTINUE                                                          
  525 CONTINUE                                                          
      IF(IRAND.EQ.0) GOTO 625                                           
      DO 626 K=1,NEGP                                                   
      KKK=KK(K+1)                                                       
      IF(K.EQ.NEGP) KKK=NDOG                                            
      KKK=KKK-KK(K)                                                     
      IF(KKK.EQ.0) GOTO 626                                             
      WRITE(6,609) K                                                    
  609 FORMAT(//,2X,'GROUP',I3,' RESIDUAL VALUES')                       
      ISTART=NSTART-MAG                                                 
      DO 622 I=1,KAG,10                                                 
      DO 623 J=1,10                                                     
  623 JX(J)=ISTART+J                                                    
      IE=MIN0(10,KAG-I+1)                                               
      WRITE(6,408)(JX(J),J=1,IE)                                        
      DO 624 J=1,KKK                                                    
      KJ=KK(K)+J                                                        
  624 WRITE(6,405) NAME(KJ),(ER(KJ,I+M-1),M=1,IE)                       
      ISTART=ISTART+10                                                  
  622 CONTINUE                                                          
  626 CONTINUE                                                          
  625 CONTINUE                                                          
      IF(IRANDX.EQ.0) GOTO 725                                          
      DO 726 K=1,NEGP                                                   
      KKK=KK(NEGP+K+1)                                                  
      IF(K.EQ.NXGP) KKK=LEXOG                                           
      KKK=KKK-KK(NEGP+K)                                                
      IF(KKK.EQ.0) GOTO 726                                             
      WRITE(6,709) (K+NEGP)                                             
  709 FORMAT(//,2X,'GROUP',I3,' RESIDUAL VALUES')                       
      ISTART=NSTART-MAG                                                 
      DO 722 I=1,KAG,10                                                 
      DO 723 J=1,10                                                     
  723 JX(J)=ISTART+J                                                    
      IE=MIN0(10,KAG-I+1)                                               
      WRITE(6,408)(JX(J),J=1,IE)                                        
      DO 724 J=1,KKK                                                    
      KJ=KK(NEGP+K)+J                                                   
  724 WRITE(6,405) IXOG(KJ),(ERX(KJ,I+M-1),M=1,IE)                      
      ISTART=ISTART+10                                                  
  722 CONTINUE                                                          
  726 CONTINUE                                                          
  725 CONTINUE                                                          
      WRITE(6,204)                                                      
      open(unit=35,file='nlag')
      read(35,*)nlag
      close(35)
      NBGN = nlag+1
      NEND = nlag+1
      open(unit=35,file='lags_endog')
      do 3013 i=1,ndog
 3013    read(35,*)(e(i,j),j=1,nlag)
      close(35)
      DO 1111 ILOOP = NBGN,NEND
         MAG = ILOOP
         LAG = MAG - 1
         KAG = MAG + 49
      open(unit=35,file='lags_exog')
      do 3015 i=1,lexog
 3015    read(35,*)(x(i,j),j=1,nlag)
      close(35)
      ICTRL=0                                                           
      ITR=-1                                                            
      IF(IDYN.GT.0) ITR=0                                               
      IF(IDYN.EQ.2) GOTO 50                                             
  100 CALL CONT(IDYN1)                                                  
      IF(IDYN.EQ.0) CALL PRINT                                             
      IF(IOPTC.GT.0) GOTO 111                                           
      CALL ALISON(NARG,FX)                                              
      GOTO 112                                                          
  111 CALL CALCFX(FX)                                              
  112 CONTINUE                                                          
      DO 602 I=1,50	  
         TOLR = 0.5		 
         IF(XEQUIV(L(I)).EQ.0.0D0) GOTO 982		 
         Z=ABS(R(I)/XEQUIV(L(I)))     		 
         IF(Z.GT.TOLR.AND.ABS(R(I)).GT.TOLR) GOTO 65		 
         GOTO 602
 982  IF(ABS(R(I)).GT.(0.1D0*TOLR)) GOTO 65     
 602  CONTINUE
      DO 603 I=51,100	  
         TOLR = 0.5		 
         IF(XEQUIV(L(I)).EQ.0.0D0) GOTO 983		 
         Z=ABS(R(I)/XEQUIV(L(I)))     		 
         IF(Z.GT.TOLR.AND.ABS(R(I)).GT.TOLR) GOTO 65		 
         GOTO 603
 983  IF(ABS(R(I)).GT.(0.1D0*TOLR)) GOTO 65     
 603  CONTINUE 
      WRITE(6,600) ITR                                                  
  600 FORMAT(' RATIONAL EXPECTATIONS VARIABLES CONVERGED IN ',I4,       
     1' ITERATIONS')                                                    
   66 CALL PRINT                                                        
        write(*,*)'before qerror'
        write(*,*)'model shocks ',(err(i,mag),i=1,10)
      errname = gen_errname(ILOOP)
      open(unit=99,file=errname)
      xname = gen_xname(ILOOP)
      open(unit=98,file=xname)
      IC=1
      DO 8010 IV = 1,16
      WRITE(99,7001) IC,IV                                              
      WRITE(99,7501) ( ER(IV,IY),IY=1,199 )
 8010 continue
      DO 8020 I=1,LEXOG
 8020 WRITE(98,7501) (X(I,J), J=1,199)
 7001 FORMAT(2I2)                                                    
 7501 FORMAT(4F18.12)                                                
      open(unit=35,file='lags_endog')
      do 3016 i=1,ndog
 3016    write(35,*)(e(i,j),j=1,mag)
      close(35)
      open(unit=35,file='lags_exog')
      do 3017 i=1,lexog
 3017    write(35,*)(x(i,j),j=1,mag)
      close(35)
      new_lag=nlag+1
      open(unit=35,file='nlag')
      write(35,*)new_lag
      close(35)
 1111 continue
      return
   65 IF(ITR.LT.MAXITR) GOTO 62                                         
      WRITE(6,601) MAXITR                                               
  601 FORMAT(' RATIONAL EXPECTATIONS VARIABLES FAILED TO CONVERGE IN ', 
     1I4,' ITERATIONS')                                                 
      CALL PRINT                                                       
   62 DO 63 I=1,NARG                                                    
      J=L(I)                                                            
   63 XEQUIV(J)=XEQUIV(J)+BSTEP*R(I)                                    
      IF(ABS(P).NE.1.0D0) GOTO 100                                      
      WRITE(6,73)                                                       
      DO 71 I=1,NARG,4                                                  
      DO 72 J=1,4                                                       
      IT=(L(I+J-1)-1)/600                                               
      IV=L(I+J-1)-IT*600                                                
      IT=IT+1                                                           
      AJX(J)=IXOG(IV)                                                   
   72 NJX(J+4)=NSTART-MAG+IT                                            
      IE=MIN0(4,NARG-I+1)                                               
      WRITE(6,74) (AJX(J),NJX(J+4),R(I+J-1),XEQUIV(L(I+J-1)),J=1,IE)    
   71 CONTINUE                                                          
   73 FORMAT(4(4X,'VARIABLE',3X,'RESIDUAL NEW VALUE'))                  
   74 FORMAT(4(1X,A8,'(',I4,')',2F9.4))                                 
      GOTO 100                                                          
   50 DO 51 I=1,NARG                                                    
   51 DIR(I)=1.0D0                                                      
      ICTRL=0                                                           
      IF(IOPTC.GT.0) GOTO 322                                           
      CALL POW
      GOTO 323                                                          
  322 CALL POW
  323 CONTINUE                                                          
      CALL PRINT                                                          
  101 FORMAT(10A8)                                                      
  201 FORMAT(3(2X,A8,F15.4,10X)/)                                       
  204 FORMAT(1H1)                                                       
      END
C===SUBROUTINE FUNCTO   
      SUBROUTINE FUNCTO                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      ITR=ITR+1                                                         
      IF(ITR.GT.0.AND.ABS(P).EQ.1.0D0) WRITE(6,320) ITR                 
      IDATE  =  NSTART                                                  
      T      =  FLOAT(IDATE-1900)                                       
      DO 3 I=MAG,KAG                                                    
      K=1                                                               
  100 DO 30 J=1,NDOG                                                    
   30 F(J)=0.0D0                                                        
      CALL EQN                                                          
      IF(IOPT.EQ.0) GOTO 55                                             
      DO 51 J=1,NDOG                                                    
      EOLD=E(J,I)                                                       
      F(J)=F(J)-EOLD                                                    
      E(J,I)=E(J,I)+B*(F(J))                                            
      SRES=F(J)/EOLD                                                    
      IF(P.LT.0.0D0) WRITE(6,300) EOLD,E(J,I),SRES,J                    
   51 CONTINUE                                                          
      GO TO 50                                                          
   55 CONTINUE                                                          
      DO 60 J=1,NDOG                                                    
      EOLD=E(J,I)-B*F(J)                                                
      RES=0.0D0                                                         
      IF(EOLD.EQ.0.0D0) GO TO 99                                        
      RES=F(J)/EOLD                                                     
   99 IF (P.LT.0.0D0) WRITE(6,300) EOLD,E(J,I),RES,J                    
   60 CONTINUE                                                          
   50 CONTINUE                                                          
      DO 61 J=1,NDOG                                                    
      IF(E(J,I).EQ.0.0D0) GO TO 98                                      
      Z=ABS(F(J)/E(J,I))                                                
      IF(Z.GE.TOL .AND. ABS(F(J)).GE.0.01D0) GO TO 62                   
      GO TO 61                                                          
   98 IF(ABS(F(J)).GE.0.001D0) GO TO 62                                 
   61 CONTINUE                                                          
      GO TO 63                                                          
   62 K=K+1                                                             
      IF(K.LE.ITMX) GOTO 100                                            
   70 WRITE(6,301) IDATE,K                                              
   63 IF(ABS(P).EQ.1.0D0) WRITE(6,310) IDATE,K                          
      IF(I.EQ.KAG) GO TO 3                                              
      IDATE  =  IDATE + 1                                               
      T      =  FLOAT(IDATE-1900)                                       
      IF(ITR.GT.1) GOTO 3                                               
    3 CONTINUE                                                          
      DO 33 J=1,NDOG                                                    
   33 F(J)=0.0D0                                                        
      RETURN                                                            
  300 FORMAT(1H ,3F15.4,I4)                                             
  310 FORMAT(10X,'PERIOD ',I4,' MODEL CONVERGED IN ',I4,' ITERATIONS')  
  301 FORMAT(10X,'PERIOD ',I4,' MODEL FAILED TO CONVERGE IN ',I4,' ITERA
     1TIONS')                                                           
  320 FORMAT(' RATIONAL EXPECTATIONS ITERATION ',I4)                    
      END 
C===SUBROUTINE CONT(IDY)                                                               
      SUBROUTINE CONT(IDY)                                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      IF(IDY.NE.0) GOTO (1,2,3,4),IDY                                   
      CALL FUNCTO                                                       
      RETURN                                                            
    1 CALL FUNCT                                                        
      RETURN                                                            
    2 CALL FUNCTP                                                       
      RETURN                                                            
    3 CALL FUNCTH                                                       
      RETURN                                                            
    4 CALL FNCTLM                                                       
      RETURN                                                            
      END
C===SUBROUTINE FUNCT                                                               
      SUBROUTINE FUNCT                                                  
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)                              
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      PARAMETER ( MXWITR = 20,                                          
     +            NUNITY =  1,                                          
     +            NZERO  =  0,                                          
     +            P01    = 0.01,                                        
     +            P001   = 0.001,                                       
     +            P0001  = 0.0001,                                      
     +            P00001 = 0.00001,                                     
     +            SMALL  = 0.000001,                                    
     +            UNITY  = 1.0      )                                   
      ITR    =  ITR + NUNITY                                            
      IDATE  =  NSTART                                                  
      T      =  FLOAT(IDATE-1900)                                       
      DO 2001 IPER = MAG,KAG                                           
      I      =  IPER                                                    
      IKW    =  NZERO                                                   
      NCONV  =  NZERO                                                   
 1000 CONTINUE                                                          
      IKW    =  IKW + NUNITY                                            
      NBLOC1 =  NBLOC - NUNITY                                          
 2000 DO 1001 IC = NUNITY,NBLOC1                                       
      K      =  NZERO                                                   
      IBLOC  =  IC                                                      
      KKIB   =  KK(IC) + NUNITY                                         
      KKIB1  =  KK(IC+NUNITY)                                           
 3000 CALL EQN                                                          
      DO 6000 J = KKIB,KKIB1                                            
      YABS  =  ABS( E(J,IPER) )                                         
      FABS  =  ABS( F(J) )                                              
      RATIO =  UNITY                                                    
      IF  ( YABS .GT. P0001 ) RATIO = FABS/YABS                         
      IF ((RATIO.LT.P0001) .OR. (FABS.LT.P0001.AND.YABS.LT.P0001)) THEN 
          GO TO 6000                                                    
      ELSE                                                              
          IF ( K .LT. ITMX ) THEN                                       
          K  =  K + NUNITY                                              
          GO TO 3000                                                    
          ELSE                                                          
          WRITE(6,4000) IDATE, IC, K                                    
 4000 FORMAT(1X,'**** ',I4,' : GROUP ',I2,' FAILED TO CONVERGE IN ',I4, 
     +          ' ITERATIONS ****')                                     
          GO TO 1001                                                   
          END IF                                                        
      END IF                                                            
 6000 CONTINUE                                                          
 1001 CONTINUE                                                          
 1002 IBLOC  =  NBLOC                                                   
      CALL EQN                                                          
      IF  ( NCONV .EQ. NZERO .AND. IKW .LT. MXWITR ) GO TO 1000         
      IF  ( NCONV .EQ. NZERO .AND. IKW .GE. MXWITR ) THEN               
      WRITE(6,1003) IDATE, IKW                                         
 1003 FORMAT(1X,'**** ',I4,' : CONSISTENCY NOT REACHED IN ',I3,         
     +          ' ITERATIONS ****')                                     
      END IF                                                            
      IF  ( IPER .EQ. KAG )  GO TO 2001                                
      IDATE  =  IDATE + NUNITY                                          
      T      =  FLOAT(IDATE-1900)                                       
 1400 CONTINUE                                                          
 2001 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE FUNCTP                                                               
      SUBROUTINE FUNCTP                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      EXTERNAL FCNP                                                     
      ITR=ITR+1                                                         
      IF (ITR.GT.0 .AND. ABS(P).EQ.1.0D0) WRITE(6,320) ITR              
      IDATE  =  NSTART                                                  
      T      =  FLOAT(IDATE-1900)                                       
      MXWITR=20                                                         
      IC=1                                                              
      IBLOC=IC                                                          
      KKIB=1                                                            
      KKIB1=NDOG                                                        
      ND=NDOG                                                           
      DO 3 I=MAG,KAG                                                    
      NCONV=0                                                           
      IKW=0                                                             
  200 CONTINUE                                                          
      IKW=IKW+1                                                         
      K=1                                                               
      DO 250 IC=1,NBLOC                                                 
      IF(NBLOC.EQ.1) GOTO 100                                           
      IBLOC=IC                                                          
      KKIB=KK(IC)+1                                                     
      KKIB1=KK(IC+1)                                                    
      ND=KK(IC+1)-KK(IC)                                                
  100 CALL EQN                                                          
      IF(NBLOC.GT.1.AND.IC.EQ.NBLOC) GOTO 53                            
      DO 51 J=1,ND                                                      
      DIR(J)=B*F(KK(IC)+J)                                              
      F(KK(IC)+J)=0.0D0                                                 
      R(J)=E(KK(IC)+J,I)                                                
   51 CONTINUE                                                          
      CALL POW
      DO 52 J=1,ND                                                      
   52 E(KK(IC)+J,I)=R(J)                                                
   53 IF (P.GE.0.0D0) GOTO 250                                          
      DO 60 J=KKIB,KKIB1                                                
      WRITE(6,300) E(J,I),J,IBLOC                                       
   60 CONTINUE                                                          
  250 CONTINUE                                                          
      IF(NCONV.EQ.1) GOTO 251                                           
      IF(IKW.LE.MXWITR) GOTO 200                                        
      WRITE(6,302)IDATE,IKW                                             
  251 IF(ABS(P).EQ.1.0D0) WRITE(6,303)IDATE,IKW                         
      IF(I.EQ.KAG) GO TO 3                                              
      IDATE  =  IDATE + 1                                               
      T      =  FLOAT(IDATE-1900)                                       
      IF(ITR.GT.1) GOTO 3                                               
      DO 4 J=1,NDOG                                                     
    4 E(J,I+1)=E(J,I)                                                   
    3 CONTINUE                                                          
      DO 33 J=1,NDOG                                                    
   33 F(J)=0.0D0                                                        
      RETURN                                                            
  300 FORMAT(1H ,F15.4,2I4)                                             
  302 FORMAT(10X,'PERIOD ',I4,' CONSITENCY NOT MET IN',I4,' ITERATIONS')
  303 FORMAT(10X,'PERIOD ',I4,' CONSITENCY     MET IN',I4,' ITERATIONS')
  320 FORMAT(' RATIONAL EXPECTATIONS ITERATION ',I4)                    
      END                                                               
      SUBROUTINE FUNCTH                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/HYBLM/XH(300),FVEC(300),WA(3130)                           
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      EXTERNAL FCNH                                                     
      NWRITE=6                                                          
      LWA=3130                                                          
      ITR=ITR+1                                                         
      IF(ITR.GT.0.AND.ABS(P).EQ.1.0D0) WRITE(6,320) ITR                 
      IDATE  =  NSTART                                                  
      T      =  FLOAT(IDATE-1900)                                       
      MXWITR=0                                                      
      IC=1                                                              
      IBLOC=IC                                                          
      KKIB=1                                                            
      KKIB1=NDOG                                                        
      ND=NDOG                                                           
      DO 3 I=MAG,KAG                                                    
      IF ( NBLOC .EQ. 1 ) THEN
      NCONV = 1
      ELSE 
      NCONV=0                                                           
      END IF
      IKW=0                                                             
  200 CONTINUE                                                          
      IKW=IKW+1                                                         
      K=1                                                               
      DO 250 IC=1,NBLOC                                                 
      IF(NBLOC.EQ.1) GOTO 100                                           
      IBLOC=IC                                                          
      KKIB=KK(IC)+1                                                     
      KKIB1=KK(IC+1)                                                    
      ND=KK(IC+1)-KK(IC)                                                
  100 NH=ND                                                             
      IF(NBLOC.GT.1.AND.IC.EQ.NBLOC) GOTO 53                            
      DO 30 J=KKIB,KKIB1                                                
   30 F(J)=0.0D0                                                        
      DO 51 JH=1,ND                                                     
   51 XH(JH)=E(KK(IC)+JH,I)                                             
      TOLH=0.01D0*TOL                                                   
      CALL HYBRD1(FCNH,NH,XH,FVEC,TOLH,INFO,WA,LWA)                     
      FNORM=ENORM(NH,FVEC)                                              
      IF(ABS(P).EQ.1.0D00) WRITE(6,301) IDATE,FNORM,INFO,IBLOC          
      DO 52 JH=1,ND                                                     
   52 E(KK(IC)+JH,I)=XH(JH)                                             
      GOTO 54                                                           
   53 CALL EQN                                                          
   54 IF (P.GE.0.0D0) GOTO 250                                          
      DO 60 J=KKIB,KKIB1                                                
      WRITE(6,300) E(J,I),J,IBLOC                                       
   60 CONTINUE
  250 CONTINUE                                                          
      IF(NCONV.EQ.1) GOTO 251                                           
      IF(IKW.LE.MXWITR) GOTO 200                                        
      WRITE(6,302)IDATE,IKW                                             
  251 IF(ABS(P).EQ.1.0D0) WRITE(6,303)IDATE,IKW                         
      IF(I.EQ.KAG) GO TO 3                                              
      IDATE  =  IDATE + 1                                               
      T      =  FLOAT(IDATE-1900)                                       
    3 CONTINUE                                                          
      DO 33 J=1,NDOG                                                    
   33 F(J)=0.0D0                                                        
      RETURN                                                            
  300 FORMAT(1H ,F15.4,2I4)                                             
  301 FORMAT(10X,'PERIOD ',I4,'FINAL L2 NORM OF RESIDUALS ',E15.7,      
     1           ' EXIT PARAMETER ',I4,' IBLOC ',I4)                    
  302 FORMAT(10X,'PERIOD ',I4,' CONSITENCY NOT MET IN',I4,' ITERATIONS')
  303 FORMAT(10X,'PERIOD ',I4,' CONSITENCY     MET IN',I4,' ITERATIONS')
  320 FORMAT(' RATIONAL EXPECTATIONS ITERATION ',I4)                    
      END
C===SUBROUTINE FNCTLM                                                               
      SUBROUTINE FNCTLM                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/HYBLM/XH(300),FVEC(300),WA(3130)                           
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      EXTERNAL FCNLM                                                    
      DIMENSION IWA(300)                                                
      NWRITE=6                                                          
      LWA=3130                                                          
      ITR=ITR+1                                                         
      IF(ITR.GT.0.AND.ABS(P).EQ.1.0D0) WRITE(6,320) ITR                 
      IDATE=NSTART                                                      
      MXWITR=20                                                         
      IC=1                                                              
      IBLOC=IC                                                          
      KKIB=1                                                            
      KKIB1=NDOG                                                        
      ND=NDOG                                                           
      DO 3 I=MAG,KAG                                                    
      NCONV=0                                                           
      IKW=0                                                             
  200 CONTINUE                                                          
      IKW=IKW+1                                                         
      K=1                                                               
      DO 250 IC=1,NBLOC                                                 
      IF(NBLOC.EQ.1) GOTO 100                                           
      IBLOC=IC                                                          
      KKIB=KK(IC)+1                                                     
      KKIB1=KK(IC+1)                                                    
      ND=KK(IC+1)-KK(IC)                                                
  100 NH=ND                                                             
      MH=NH                                                             
      IF(NBLOC.GT.1.AND.IC.EQ.NBLOC) GOTO 53                            
      DO 30 J=KKIB,KKIB1                                                
   30 F(J)=0.0D0                                                        
      DO 51 JH=1,ND                                                     
   51 XH(JH)=E(KK(IC)+JH,I)                                             
      TOLH=0.01D0*TOL                                                   
      CALL LMDIF1(FCNLM,MH,NH,XH,FVEC,TOLH,INFO,IWA,WA,LWA)             
      FNORM=ENORM(NH,FVEC)                                              
      IF(ABS(P).EQ.1.0D0) WRITE(6,301) IDATE,FNORM,INFO,IBLOC           
      DO 52 JH=1,ND                                                     
   52 E(KK(IC)+JH,I)=XH(JH)                                             
      GOTO 54                                                           
   53 CALL EQN                                                          
   54 IF (P.GE.0.0D0) GOTO 250                                          
      DO 60 J=KKIB,KKIB1                                                
      WRITE(6,300) E(J,I),J,IBLOC                                       
   60 CONTINUE                                                          
  250 CONTINUE                                                          
      IF(NCONV.EQ.1) GOTO 251                                           
      IF(IKW.LE.MXWITR) GOTO 200                                        
      WRITE(6,302)IDATE,IKW                                             
  251 IF(ABS(P).EQ.1.0D0) WRITE(6,303)IDATE,IKW                         
      IF(I.EQ.KAG) GO TO 3                                              
      IDATE=IDATE+1                                                     
      IF(ITR.GT.1) GOTO 3                                               
      DO 4 J=1,NDOG                                                     
    4 E(J,I+1)=E(J,I)                                                   
      T=T+1                                                             
    3 CONTINUE                                                          
      DO 33 J=1,NDOG                                                    
   33 F(J)=0.0D0                                                        
      RETURN                                                            
  300 FORMAT(1H ,F15.4,2I4)                                             
  301 FORMAT(10X,'PERIOD ',I4,'FINAL L2 NORM OF RESIDUALS ',E15.7,      
     1           ' EXIT PARAMETER ',I4,' IBLOC ',I4)                    
  302 FORMAT(10X,'PERIOD ',I4,' CONSITENCY NOT MET IN',I4,' ITERATIONS')
  303 FORMAT(10X,'PERIOD ',I4,' CONSITENCY     MET IN',I4,' ITERATIONS')
  320 FORMAT(' RATIONAL EXPECTATIONS ITERATION ',I4)                    
      END                                                               
      SUBROUTINE FCNH(NH,XH,FVEC,IFLAG)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      DIMENSION XH(NH),FVEC(NH)                                         
      IC=IBLOC                                                          
      DO 10 JH=1,NH                                                     
   10 E(KK(IC)+JH,I)=XH(JH)                                             
      CALL EQN                                                          
      IF(IOPT.EQ.0) GOTO 15                                             
      DO 11 JH=1,NH                                                     
   11 F(KK(IC)+JH)=F(KK(IC)+JH)-E(KK(IC)+JH,I)                          
   15 CONTINUE                                                          
      DO 20 JH=1,NH                                                     
   20 FVEC(JH)=F(KK(IC)+JH)                                             
      RETURN                                                            
      END                                                               
      SUBROUTINE FCNLM(MH,NH,XH,FVEC,IFLAG)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      INTEGER MH,NH,IFLAG                                               
      REAL XH(NH),FVEC(MH)                                              
      IC=IBLOC                                                          
      DO 10 JH=1,NH                                                     
   10 E(KK(IC)+JH,I)=XH(JH)                                             
      CALL EQN                                                          
      IF(IOPT.EQ.0) GOTO 15                                             
      DO 11 JH=1,MH                                                     
   11 F(KK(IC)+JH)=F(KK(IC)+JH)-E(KK(IC)+JH,I)                          
   15 CONTINUE                                                          
      DO 20 JH=1,MH                                                     
   20 FVEC(JH)=F(KK(IC)+JH)                                             
      RETURN                                                            
      END
C===SUBROUTINE FCNP(NH,FX)                                                               
      SUBROUTINE FCNP(NH,FX)                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      IC=IBLOC                                                          
      FX=0.0D0                                                          
      DO 10 JH=1,NH                                                     
   10 E(KK(IC)+JH,I)=R(JH)                                              
      CALL EQN                                                          
      IF(IOPT.EQ.0) GOTO 15                                             
      DO 11 JH=1,NH                                                     
   11 F(KK(IC)+JH)=F(KK(IC)+JH)-E(KK(IC)+JH,I)                          
   15 CONTINUE                                                          
      DO 20 JH=1,NH                                                     
   20 FX=FX+F(KK(IC)+JH)*F(KK(IC)+JH)                                   
      RETURN                                                            
      END
C===SUBROUTINE WORJAC(ILC)                                                               
      SUBROUTINE WORJAC(ILC)                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/V(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      I=IPER                                                            
      KNE=KK(ILC+1)-KK(ILC)                                             
      KIC=KK(ILC)                                                       
      DO 100 IE=1,KNE                                                   
      EOLD=E(KIC+IE,I)                                                  
      V(KIC+IE)=V(KIC+IE)-E(KIC+IE,I)                                   
      E(KIC+IE,I)=E(KIC+IE,I)+B*V(KIC+IE)                               
      ZR=0.0D0                                                          
      IF(EOLD.EQ.0.0D0) GOTO 90                                         
      ZR=V(KIC+IE)/EOLD                                                 
   90 IF(ABS(P).EQ.1.0D0)WRITE(6,300) EOLD,E(KIC+IE,I),ZR,(KIC+IE)      
  100 CONTINUE                                                          
      DO 105 IE=1,KNE  
      IF(E(KIC+IE,I).EQ.0.0D0) GOTO 200                                 
      ZR=ABS(V(KIC+IE)/E(KIC+IE,I))                                     
      IF(ZR.GT.TOL.AND.ABS(V(KIC+IE)).GT.TOL) GOTO 250                  
      GOTO 105                                                          
  200 IF(ABS(V(KIC+IE)).GT.(0.1D0*TOL)) GOTO 250                        
  105 CONTINUE                                                          
      NCONV=1                                                           
  250 CONTINUE                                                          
  300 FORMAT(1H ,3F15.4,I4)                                             
      DO 110 IE=1,KNE                                                   
  110 V(KIC+IE)=E(KIC+IE,I)                                             
      RETURN                                                            
      END
C===SUBROUTINE PRINT               
      SUBROUTINE PRINT                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG, AJX(10)                                 
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      CHARACTER*10  TITLE1, TITLE2                                      
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      DIMENSION XEQUIV(120000)                                           
      DIMENSION  NJX(10)                                                
      EQUIVALENCE ( XEQUIV(1),X(1,1) )                                  
      DO 1005 I = 1,NDOG                                               
      WRITE(12,9000)  ( E(I,J),J=1,147 )                               
 1005 CONTINUE                                                          
      DO 1201 I = 1,LEXOG                                              
      WRITE(12,9000)  ( X(I,J),J=1,147 )                               
 1201 CONTINUE                                                          
 9000 FORMAT(4F18.12)                                                    
      TITLE1='     SIMUL'                                               
      TITLE2='ATION     '                                               
      IF(IDYN.EQ.0) GOTO 501                                            
      WRITE(6,73)                                                       
      DO 71 I=1,NARG,4                                                  
      DO 72 J=1,4                                                       
      IT=(L(I+J-1)-1)/600                                               
      IV=L(I+J-1)-IT*600                                                
      IT=IT+1                                                           
      AJX(J)=IXOG(IV)                                                   
   72 NJX(J+4)=NSTART-MAG+IT                                            
      IE=MIN0(4,NARG-I+1)                                               
      WRITE(6,74) (AJX(J),NJX(J+4),R(I+J-1),XEQUIV(L(I+J-1)),J=1,IE)    
   71 CONTINUE                                                          
   73 FORMAT(4(4X,'VARIABLE',3X,'RESIDUAL     VALUE'))                  
   74 FORMAT(4(1X,A8,'(',I4,')',2F9.4))                                 
      CALL NEWRES                                                       
  501 NSTA=NSTART                                                       
      MG1=MAG                                                           
      IF(MAG.EQ.1) GOTO 500                                             
      MG1=MAG-1                                                         
      NSTA=NSTART-1                                                     
  500 DO 1 K=1,NEGP                                                     
      KK1=KK(K+1)                                                       
      IF(K.EQ.NEGP) KK1=NDOG                                            
      KK1=KK1-KK(K)                                                     
      KK2=0                                                             
      IF(K.LE.NXGP)KK2=KK(NEGP+K+1)-KK(NEGP+K)                          
      IF(KK1.EQ.0.AND.KK2.EQ.0) GOTO 1                                  
      WRITE(6,400) TITLE1,TITLE2,MAG,kag                 
  400 FORMAT(1H1,40X,2A10,'FROM  ',I4,'  TO  ',I4)                      
      IF(KK1.EQ.0) GOTO 2                                               
      ISTART=NSTA  -1                                                   
      ISTART=MAG-2
      WRITE(6,150) K                                                    
  150 FORMAT(//,2X,'GROUP',I3,' ENDOGENOUS VARIABLES')                  
      mgm1=mag-2
      DO 100 I=MG1,KAG,10                                               
      DO 101 J=1,10                                                     
  101 JX(J)=ISTART+J                                                    
      IE=MIN0(10,KAG-I+1)                                               
      WRITE(6,200) (JX(J),J=1,IE)                                       
  200 FORMAT(/14X,I4,9(8X,I4))                                          
      DO 102 J=1,KK1                                                    
      KJ=KK(K)+J                                                        
  102 WRITE(6,300) NAME(KJ),(E(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+10                                                  
  100 CONTINUE                                                          
    2 IF(KK2.EQ.0.OR.PX.GE.0.0) GOTO 1                                  
      K2=K+NEGP                                                         
      ISTART=NSTA  -1                                                   
      WRITE(6,151) K2                                                   
  151 FORMAT(//,2X,'GROUP',I3,' EXOGENOUS VARIABLES')                   
      DO 120 I=MG1,KAG,10                                               
      DO 121 J=1,10                                                     
  121 JX(J)=ISTART+J                                                    
      IE=MIN0(10,KAG-I+1)                                               
      WRITE(6,200) (JX(J),J=1,IE)                                       
      DO 122 J=1,KK2                                                    
      KJ=KK(NEGP+K)+J                                                   
  122 WRITE(6,300) IXOG(KJ),(X(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+10                                                  
  120 CONTINUE                                                          
  300 FORMAT(2X,A8,10F12.4)                                             
    1 CONTINUE                                                          
      IF(PR.NE.1.0D0) GOTO 3                                            
      DO 4 I=MAG,KAG                                                    
    4 WRITE(7,1000) (E(J,I),J=1 ,NDOG)                                  
 1000 FORMAT(5F16.5)                                                    
    3 IF(BR.NE.1.0D0) RETURN                                            
      CALL BASE                                                         
      XX = 0.0D0                                                        
      BR = 0.0D0                                                        
      PR = 0.0D0                                                        
      TITLE1='        DI'                                               
      TITLE2='FFERENCES '                                               
      GOTO 500                                                          
      END 
C===SUBROUTINE PRINTX                                                               
      SUBROUTINE PRINTX                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG, AJX(10)                                 
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      CHARACTER*10  TITLE1, TITLE2                                      
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      DIMENSION XEQUIV(120000)                                           
      DIMENSION  NJX(10)                                                
      EQUIVALENCE ( XEQUIV(1),X(1,1) )                                  
      DO 1006 I = 1,NDOG                                               
      WRITE(12,9000)  ( E(I,J),J=1,KAG )                               
 1006 CONTINUE                                                          
      DO 1207 I = 1,LEXOG                                              
      WRITE(12,9000)  ( X(I,J),J=1,KAG )                               
 1207 CONTINUE                                                          
 9000 FORMAT(4F18.12)                                                    
      TITLE1='     SIMUL'                                               
      TITLE2='ATION     '                                               
      IF(IDYN.EQ.0) GOTO 501                                            
      WRITE(6,73)                                                       
      DO 71 I=1,NARG,4                                                  
      DO 72 J=1,4                                                       
      IT=(L(I+J-1)-1)/600                                               
      IV=L(I+J-1)-IT*600                                                
      IT=IT+1                                                           
      AJX(J)=IXOG(IV)                                                   
   72 NJX(J+4)=NSTART-MAG+IT                                            
      IE=MIN0(4,NARG-I+1)                                               
      WRITE(6,74) (AJX(J),NJX(J+4),R(I+J-1),XEQUIV(L(I+J-1)),J=1,IE)    
   71 CONTINUE                                                          
   73 FORMAT(4(4X,'VARIABLE',3X,'RESIDUAL     VALUE'))                  
   74 FORMAT(4(1X,A8,'(',I4,')',2F9.4))                                 
      CALL NEWRES                                                       
  501 NSTA=NSTART                                                       
      MG1=MAG                                                           
      IF(MAG.EQ.1) GOTO 500                                             
      MG1=MAG-1                                                         
      NSTA=NSTART-1                                                     
  500 DO 1 K=1,NEGP                                                     
      KK1=KK(K+1)                                                       
      IF(K.EQ.NEGP) KK1=NDOG                                            
      KK1=KK1-KK(K)                                                     
      KK2=0                                                             
      IF(K.LE.NXGP)KK2=KK(NEGP+K+1)-KK(NEGP+K)                          
      IF(KK1.EQ.0.AND.KK2.EQ.0) GOTO 1                                  
      WRITE(6,400) TITLE1,TITLE2,NSTART,(NSTART+NPER-1)                 
  400 FORMAT(1H1,40X,2A10,'FROM  ',I4,'  TO  ',I4)                      
      IF(KK1.EQ.0) GOTO 2                                               
      ISTART=NSTA  -1                                                   
      WRITE(6,150) K                                                    
  150 FORMAT(//,2X,'GROUP',I3,' ENDOGENOUS VARIABLES')                  
      DO 100 I=MG1,KAG,5                                               
      DO 101 J=1,5
  101 JX(J)=ISTART+J                                                    
      IE=MIN0(5,KAG-I+1)                                               
      WRITE(6,200) (JX(J),J=1,IE)                                       
  200 FORMAT(/14X,I4,9(8X,I4))                                          
      DO 102 J=1,KK1                                                    
      KJ=KK(K)+J                                                        
  102 WRITE(6,300) NAME(KJ),(E(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+5
  100 CONTINUE                                                          
    2 IF(KK2.EQ.0.OR.PX.GE.0.0) GOTO 1                                  
      K2=K+NEGP                                                         
      ISTART=NSTA  -1                                                   
      WRITE(6,151) K2                                                   
  151 FORMAT(//,2X,'GROUP',I3,' EXOGENOUS VARIABLES')                   
      DO 120 I=MG1,KAG,5
      DO 121 J=1,5
  121 JX(J)=ISTART+J                                                    
      IE=MIN0(5,KAG-I+1)                                               
      WRITE(6,200) (JX(J),J=1,IE)                                       
      DO 122 J=1,KK2                                                    
      KJ=KK(NEGP+K)+J                                                   
  122 WRITE(6,300) IXOG(KJ),(X(KJ,I+M-1),M=1,IE)                        
      ISTART=ISTART+5
  120 CONTINUE                                                          
  300 FORMAT(2X,A8,5F12.4)                                             
    1 CONTINUE                                                          
      IF(PR.NE.1.0D0) GOTO 3                                            
      DO 4 I=MAG,KAG                                                    
    4 WRITE(7,1000) (E(J,I),J=1 ,NDOG)                                  
 1000 FORMAT(5F16.5)                                                    
    3 IF(BR.NE.1.0D0) RETURN                                            
      CALL BASE                                                         
      XX = 0.0D0                                                        
      BR = 0.0D0                                                        
      PR = 0.0D0                                                        
      TITLE1='        DI'                                               
      TITLE2='FFERENCES '                                               
      GOTO 500                                                          
      END
C===SUBROUTINE BASE                                                               
      SUBROUTINE BASE                                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      NGP  =  13                                                        
      DO 1000 I = MAG,KAG                                               
      READ(4,1500)  ( EE(IV,I), IV = 1,NDOG )                           
 1000 CONTINUE                                                          
 1500 FORMAT(5F16.5)                                                    
      DO 9000 I = MAG,KAG                                               
      DO 9000 IV = 1,NDOG                                               
      E(IV,I)  =  E(IV,I) - EE(IV,I)                                    
 9000 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE SHOCK                                                               
      SUBROUTINE SHOCK                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)                                              
      KJ=KAG                                                            
      IF(NSK.EQ.1) KJ=MAG                                               
      DO 25 J=MAG,KJ                                                    
      DO 24 K=1,IJX                                                     
      IGP=JX(K)/100                                                     
      IVB=JX(K)-IGP*100                                                 
   24 X(KK(IGP)+IVB,J)=X(KK(IGP)+IVB,J)+VX(K)                           
   25 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE PSHOCK                                                               
      SUBROUTINE PSHOCK                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      READ(5,20) (JX(I),I=1,IJX)                                        
      READ(5,*) (VX(I),I=1,IJX)                                         
      IF(NSK.EQ.1) GO TO 5                                              
      WRITE(6,24)                                                       
      GO TO 7                                                           
    5 WRITE(6,25)                                                       
    7 WRITE(6,26) (JX(I),VX(I),I=1,IJX)                                 
      RETURN                                                            
   20 FORMAT(20I4)                                                      
   24 FORMAT(20X,'EXOGENOUS VARIABLES SHOCK:CONTINUOUS')                
   25 FORMAT(20X,'EXOGENOUS VARIABLES SHOCK:IMPACT')                    
   26 FORMAT(20X,'EXOGENOUS VARIABLE',I4,20X,'VALUES=',F10.4)           
      END
C===DOUBLE PRECISION FUNCTION VV(I,J,K)                                                                
      DOUBLE PRECISION FUNCTION VV(I,J,K)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      L=KK(I)                                                           
      IF(I.LE.NEGP) VV=E(L+J,IPER+K)                                    
      IF(I.GT.NEGP) VV=X(L+J,IPER+K)                                    
      RETURN                                                            
      END
C===DOUBLE PRECISION FUNCTION RR(I,J,K)                                                                
      DOUBLE PRECISION FUNCTION RR(I,J,K)                               
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)                              
      CHARACTER*8   NAME, IXOG                                          
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      L=KK(I)+J                                                         
      RR=ER(L,IPER+K)                                                   
      RETURN                                                            
      END
C===SUBROUTINE POW                                                                
      SUBROUTINE POW 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      dummy = 0.0
      END 
C===SUBROUTINE HYBRD                                                              
      SUBROUTINE HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG,      
     *                 MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,R,LR,   
     *                 QTF,WA1,WA2,WA3,WA4)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),WA1(N),
     *     WA2(N),WA3(N),WA4(N)                                         
      EXTERNAL FCN                                                      
      DIMENSION IWA(1)                                                  
      LOGICAL JEVAL,SING                                                
      DATA ONE,P1,P5,P001,P0001,ZERO                                    
     *     /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/                    
      EPSMCH = SPMPAR(1)                                                
      INFO = 0                                                          
      IFLAG = 0                                                         
      NFEV = 0                                                          
      IF (N .LE. 0 .OR. XTOL .LT. ZERO .OR. MAXFEV .LE. 0               
     *    .OR. ML .LT. 0 .OR. MU .LT. 0 .OR. FACTOR .LE. ZERO           
     *    .OR. LDFJAC .LT. N .OR. LR .LT. (N*(N + 1))/2) GO TO 300      
      IF (MODE .NE. 2) GO TO 20                                         
      DO 10 J = 1, N                                                    
         IF (DIAG(J) .LE. ZERO) GO TO 300                               
   10    CONTINUE                                                       
   20 CONTINUE                                                          
      IFLAG = 1                                                         
      CALL FCN(N,X,FVEC,IFLAG)                                          
      NFEV = 1                                                          
      IF (IFLAG .LT. 0) GO TO 300                                       
      FNORM = ENORM(N,FVEC)                                             
      MSUM = MIN0(ML+MU+1,N)                                            
      ITER = 1                                                          
      NCSUC = 0                                                         
      NCFAIL = 0                                                        
      NSLOW1 = 0                                                        
      NSLOW2 = 0                                                        
   30 CONTINUE                                                          
         JEVAL = .TRUE.                                                 
         IFLAG = 2                                                      
         CALL FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,   
     *               WA2)                                               
         NFEV = NFEV + MSUM                                             
         IF (IFLAG .LT. 0) GO TO 300                                    
         CALL QRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)          
         IF (ITER .NE. 1) GO TO 70                                      
         IF (MODE .EQ. 2) GO TO 50                                      
         DO 40 J = 1, N                                                 
            DIAG(J) = WA2(J)                                            
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE                         
   40       CONTINUE                                                    
   50    CONTINUE                                                       
         DO 60 J = 1, N                                                 
            WA3(J) = DIAG(J)*X(J)                                       
   60       CONTINUE                                                    
         XNORM = ENORM(N,WA3)                                           
         DELTA = FACTOR*XNORM                                           
         IF (DELTA .EQ. ZERO) DELTA = FACTOR                            
   70    CONTINUE                                                       
         DO 80 I = 1, N                                                 
            QTF(I) = FVEC(I)                                            
   80       CONTINUE                                                    
         DO 120 J = 1, N                                                
            IF (FJAC(J,J) .EQ. ZERO) GO TO 110                          
            SUM = ZERO                                                  
            DO 90 I = J, N                                              
               SUM = SUM + FJAC(I,J)*QTF(I)                             
   90          CONTINUE                                                 
            TEMP = -SUM/FJAC(J,J)                                       
            DO 100 I = J, N                                             
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP                         
  100          CONTINUE                                                 
  110       CONTINUE                                                    
  120       CONTINUE                                                    
         SING = .FALSE.                                                 
         DO 150 J = 1, N                                                
            L = J                                                       
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 140                                   
            DO 130 I = 1, JM1                                           
               R(L) = FJAC(I,J)                                         
               L = L + N - I                                            
  130          CONTINUE                                                 
  140       CONTINUE                                                    
            R(L) = WA1(J)                                               
            IF (WA1(J) .EQ. ZERO) SING = .TRUE.                         
  150       CONTINUE                                                    
         CALL QFORM(N,N,FJAC,LDFJAC,WA1)                                
         IF (MODE .EQ. 2) GO TO 170                                     
         DO 160 J = 1, N                                                
            DIAG(J) = MAX(DIAG(J),WA2(J))                               
  160       CONTINUE                                                    
  170    CONTINUE                                                       
  180    CONTINUE                                                       
            IF (NPRINT .LE. 0) GO TO 190                                
            IFLAG = 0                                                   
            IF (MOD(ITER-1,NPRINT) .EQ. 0) CALL FCN(N,X,FVEC,IFLAG)     
            IF (IFLAG .LT. 0) GO TO 300                                 
  190       CONTINUE                                                    
            CALL DOGLEG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)              
            DO 200 J = 1, N                                             
               WA1(J) = -WA1(J)                                         
               WA2(J) = X(J) + WA1(J)                                   
               WA3(J) = DIAG(J)*WA1(J)                                  
  200          CONTINUE                                                 
            PNORM = ENORM(N,WA3)                                        
            IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)                   
            IFLAG = 1                                                   
            CALL FCN(N,WA2,WA4,IFLAG)                                   
            NFEV = NFEV + 1                                             
            IF (IFLAG .LT. 0) GO TO 300                                 
            FNORM1 = ENORM(N,WA4)                                       
            ACTRED = -ONE                                               
            IF (FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2     
            L = 1                                                       
            DO 220 I = 1, N                                             
               SUM = ZERO                                               
               DO 210 J = I, N                                          
                  SUM = SUM + R(L)*WA1(J)                               
                  L = L + 1                                             
  210             CONTINUE                                              
               WA3(I) = QTF(I) + SUM                                    
  220          CONTINUE                                                 
            TEMP = ENORM(N,WA3)                                         
            PRERED = ZERO                                               
            IF (TEMP .LT. FNORM) PRERED = ONE - (TEMP/FNORM)**2         
            RATIO = ZERO                                                
            IF (PRERED .GT. ZERO) RATIO = ACTRED/PRERED                 
            IF (RATIO .GE. P1) GO TO 230                                
               NCSUC = 0                                                
               NCFAIL = NCFAIL + 1                                      
               DELTA = P5*DELTA                                         
               GO TO 240                                                
  230       CONTINUE                                                    
               NCFAIL = 0                                               
               NCSUC = NCSUC + 1                                        
               IF (RATIO .GE. P5 .OR. NCSUC .GT. 1)                     
     *            DELTA = MAX(DELTA,PNORM/P5)                           
               IF (ABS(RATIO-ONE) .LE. P1) DELTA = PNORM/P5             
  240       CONTINUE                                                    
            IF (RATIO .LT. P0001) GO TO 260                             
            DO 250 J = 1, N                                             
               X(J) = WA2(J)                                            
               WA2(J) = DIAG(J)*X(J)                                    
               FVEC(J) = WA4(J)                                         
  250          CONTINUE                                                 
            XNORM = ENORM(N,WA2)                                        
            FNORM = FNORM1                                              
            ITER = ITER + 1                                             
  260       CONTINUE                                                    
            NSLOW1 = NSLOW1 + 1                                         
            IF (ACTRED .GE. P001) NSLOW1 = 0                            
            IF (JEVAL) NSLOW2 = NSLOW2 + 1                              
            IF (ACTRED .GE. P1) NSLOW2 = 0                              
            IF (DELTA .LE. XTOL*XNORM .OR. FNORM .EQ. ZERO) INFO = 1    
            IF (INFO .NE. 0) GO TO 300                                  
            IF (NFEV .GE. MAXFEV) INFO = 2                              
            IF (P1*MAX(P1*DELTA,PNORM) .LE. EPSMCH*XNORM) INFO = 3      
            IF (NSLOW2 .EQ. 5) INFO = 4                                 
            IF (NSLOW1 .EQ. 10) INFO = 5                                
            IF (INFO .NE. 0) GO TO 300                                  
            IF (NCFAIL .EQ. 2) GO TO 290                                
            DO 280 J = 1, N                                             
               SUM = ZERO                                               
               DO 270 I = 1, N                                          
                  SUM = SUM + FJAC(I,J)*WA4(I)                          
  270             CONTINUE                                              
               WA2(J) = (SUM - WA3(J))/PNORM                            
               WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)                
               IF (RATIO .GE. P0001) QTF(J) = SUM                       
  280          CONTINUE                                                 
            CALL R1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)                      
            CALL R1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)                        
            CALL R1MPYQ(1,N,QTF,1,WA2,WA3)                              
            JEVAL = .FALSE.                                             
            GO TO 180                                                   
  290    CONTINUE                                                       
         GO TO 30                                                       
  300 CONTINUE                                                          
      IF (IFLAG .LT. 0) INFO = IFLAG                                    
      IFLAG = 0                                                         
      IF (NPRINT .GT. 0) CALL FCN(N,X,FVEC,IFLAG)                       
      RETURN                                                            
      END
C===SUBROUTINE HYBRD1                                                               
      SUBROUTINE HYBRD1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),FVEC(N),WA(LWA)                                    
      EXTERNAL FCN                                                      
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/                          
      INFO = 0                                                          
      IF (N .LE. 0 .OR. TOL .LT. ZERO .OR. LWA .LT. (N*(3*N + 13))/2)   
     *   GO TO 20                                                       
      MAXFEV = 200*(N + 1)                                              
      XTOL = TOL                                                        
      ML = N - 1                                                        
      MU = N - 1                                                        
      EPSFCN = ZERO                                                     
      MODE = 2                                                          
      DO 10 J = 1, N                                                    
         WA(J) = ONE                                                    
   10    CONTINUE                                                       
      NPRINT = 0                                                        
      LR = (N*(N + 1))/2                                                
      INDEX = 6*N + LR                                                  
      CALL HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,WA(1),MODE,      
     *           FACTOR,NPRINT,INFO,NFEV,WA(INDEX+1),N,WA(6*N+1),LR,    
     *           WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))       
      IF (INFO .EQ. 5) INFO = 4                                         
   20 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE LMDIF                                                               
      SUBROUTINE LMDIF(FCN,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,     
     *                 DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,   
     *                 IPVT,QTF,WA1,WA2,WA3,WA4)                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION IPVT(N)                                                 
      DIMENSION X(N),FVEC(M),DIAG(N),FJAC(LDFJAC,N),QTF(N),WA1(N),      
     *     WA2(N),WA3(N),WA4(M)                                         
      EXTERNAL FCN                                                      
      DATA ONE,P1,P5,P25,P75,P0001,ZERO                                 
     *     /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/             
      EPSMCH = SPMPAR(1)                                                
      INFO = 0                                                          
      IFLAG = 0                                                         
      NFEV = 0                                                          
      IF (N .LE. 0 .OR. M .LT. N .OR. LDFJAC .LT. M                     
     *    .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO   
     *    .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300           
      IF (MODE .NE. 2) GO TO 20                                         
      DO 10 J = 1, N                                                    
         IF (DIAG(J) .LE. ZERO) GO TO 300                               
   10    CONTINUE                                                       
   20 CONTINUE                                                          
      IFLAG = 1                                                         
      CALL FCN(M,N,X,FVEC,IFLAG)                                        
      NFEV = 1                                                          
      IF (IFLAG .LT. 0) GO TO 300                                       
      FNORM = ENORM(M,FVEC)                                             
      PAR = ZERO                                                        
      ITER = 1                                                          
   30 CONTINUE                                                          
         IFLAG = 2                                                      
         CALL FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)       
         NFEV = NFEV + N                                                
         IF (IFLAG .LT. 0) GO TO 300                                    
         IF (NPRINT .LE. 0) GO TO 40                                    
         IFLAG = 0                                                      
         IF (MOD(ITER-1,NPRINT) .EQ. 0) CALL FCN(M,N,X,FVEC,IFLAG)      
         IF (IFLAG .LT. 0) GO TO 300                                    
   40    CONTINUE                                                       
         CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)          
         IF (ITER .NE. 1) GO TO 80                                      
         IF (MODE .EQ. 2) GO TO 60                                      
         DO 50 J = 1, N                                                 
            DIAG(J) = WA2(J)                                            
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE                         
   50       CONTINUE                                                    
   60    CONTINUE                                                       
         DO 70 J = 1, N                                                 
            WA3(J) = DIAG(J)*X(J)                                       
   70       CONTINUE                                                    
         XNORM = ENORM(N,WA3)                                           
         DELTA = FACTOR*XNORM                                           
         IF (DELTA .EQ. ZERO) DELTA = FACTOR                            
   80    CONTINUE                                                       
         DO 90 I = 1, M                                                 
            WA4(I) = FVEC(I)                                            
   90       CONTINUE                                                    
         DO 130 J = 1, N                                                
            IF (FJAC(J,J) .EQ. ZERO) GO TO 120                          
            SUM = ZERO                                                  
            DO 100 I = J, M                                             
               SUM = SUM + FJAC(I,J)*WA4(I)                             
  100          CONTINUE                                                 
            TEMP = -SUM/FJAC(J,J)                                       
            DO 110 I = J, M                                             
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP                         
  110          CONTINUE                                                 
  120       CONTINUE                                                    
            FJAC(J,J) = WA1(J)                                          
            QTF(J) = WA4(J)                                             
  130       CONTINUE                                                    
         GNORM = ZERO                                                   
         IF (FNORM .EQ. ZERO) GO TO 170                                 
         DO 160 J = 1, N                                                
            L = IPVT(J)                                                 
            IF (WA2(L) .EQ. ZERO) GO TO 150                             
            SUM = ZERO                                                  
            DO 140 I = 1, J                                             
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)                     
  140          CONTINUE                                                 
            GNORM = MAX(GNORM,ABS(SUM/WA2(L)))                          
  150       CONTINUE                                                    
  160       CONTINUE                                                    
  170    CONTINUE                                                       
         IF (GNORM .LE. GTOL) INFO = 4                                  
         IF (INFO .NE. 0) GO TO 300                                     
         IF (MODE .EQ. 2) GO TO 190                                     
         DO 180 J = 1, N                                                
            DIAG(J) = MAX(DIAG(J),WA2(J))                               
  180       CONTINUE                                                    
  190    CONTINUE                                                       
  200    CONTINUE                                                       
            CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2,   
     *                 WA3,WA4)                                         
            DO 210 J = 1, N                                             
               WA1(J) = -WA1(J)                                         
               WA2(J) = X(J) + WA1(J)                                   
               WA3(J) = DIAG(J)*WA1(J)                                  
  210          CONTINUE                                                 
            PNORM = ENORM(N,WA3)                                        
            IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)                   
            IFLAG = 1                                                   
            CALL FCN(M,N,WA2,WA4,IFLAG)                                 
            NFEV = NFEV + 1                                             
            IF (IFLAG .LT. 0) GO TO 300                                 
            FNORM1 = ENORM(M,WA4)                                       
            ACTRED = -ONE                                               
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2  
            DO 230 J = 1, N                                             
               WA3(J) = ZERO                                            
               L = IPVT(J)                                              
               TEMP = WA1(L)                                            
               DO 220 I = 1, J                                          
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP                      
  220             CONTINUE                                              
  230          CONTINUE                                                 
            TEMP1 = ENORM(N,WA3)/FNORM                                  
            TEMP2 = (SQRT(PAR)*PNORM)/FNORM                             
            PRERED = TEMP1**2 + TEMP2**2/P5                             
            DIRDER = -(TEMP1**2 + TEMP2**2)                             
            RATIO = ZERO                                                
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED                 
            IF (RATIO .GT. P25) GO TO 240                               
               IF (ACTRED .GE. ZERO) TEMP = P5                          
               IF (ACTRED .LT. ZERO)                                    
     *            TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)                 
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1    
               DELTA = TEMP*MIN(DELTA,PNORM/P1)                         
               PAR = PAR/TEMP                                           
               GO TO 260                                                
  240       CONTINUE                                                    
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250        
               DELTA = PNORM/P5                                         
               PAR = P5*PAR                                             
  250          CONTINUE                                                 
  260       CONTINUE                                                    
            IF (RATIO .LT. P0001) GO TO 290                             
            DO 270 J = 1, N                                             
               X(J) = WA2(J)                                            
               WA2(J) = DIAG(J)*X(J)                                    
  270          CONTINUE                                                 
            DO 280 I = 1, M                                             
               FVEC(I) = WA4(I)                                         
  280          CONTINUE                                                 
            XNORM = ENORM(N,WA2)                                        
            FNORM = FNORM1                                              
            ITER = ITER + 1                                             
  290       CONTINUE                                                    
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL            
     *          .AND. P5*RATIO .LE. ONE) INFO = 1                       
            IF (DELTA .LE. XTOL*XNORM) INFO = 2                         
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL            
     *          .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3     
            IF (INFO .NE. 0) GO TO 300                                  
            IF (NFEV .GE. MAXFEV) INFO = 5                              
            IF (ABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH        
     *          .AND. P5*RATIO .LE. ONE) INFO = 6                       
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7                       
            IF (GNORM .LE. EPSMCH) INFO = 8                             
            IF (INFO .NE. 0) GO TO 300                                  
            IF (RATIO .LT. P0001) GO TO 200                             
         GO TO 30                                                       
  300 CONTINUE                                                          
      IF (IFLAG .LT. 0) INFO = IFLAG                                    
      IFLAG = 0                                                         
      IF (NPRINT .GT. 0) CALL FCN(M,N,X,FVEC,IFLAG)                     
      RETURN                                                            
      END
C===SUBROUTINE LMDIF1                                                               
      SUBROUTINE LMDIF1(FCN,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION IWA(N)                                                  
      DIMENSION X(N),FVEC(M),WA(LWA)                                    
      EXTERNAL FCN                                                      
      INTEGER MAXFEV,MODE,MP5N,NFEV,NPRINT                              
      REAL EPSFCN,FACTOR,FTOL,GTOL,XTOL,ZERO                            
      DATA FACTOR,ZERO /1.0D2,0.0D0/                                    
      INFO = 0                                                          
      IF (N .LE. 0 .OR. M .LT. N .OR. TOL .LT. ZERO                     
     *    .OR. LWA .LT. M*N + 5*N + M) GO TO 10                         
      MAXFEV = 200*(N + 1)                                              
      FTOL = TOL                                                        
      XTOL = TOL                                                        
      GTOL = ZERO                                                       
      EPSFCN = ZERO                                                     
      MODE = 1                                                          
      NPRINT = 0                                                        
      MP5N = M + 5*N                                                    
      IF (INFO .EQ. 8) INFO = 4                                         
   10 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)                                                                
      SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION R(LR),DIAG(N),QTB(N),X(N),WA1(N),WA2(N)                 
      DATA ONE,ZERO /1.0D0,0.0D0/                                       
      EPSMCH = SPMPAR(1)                                                
      JJ = (N*(N + 1))/2 + 1                                            
      DO 50 K = 1, N                                                    
         J = N - K + 1                                                  
         JP1 = J + 1                                                    
         JJ = JJ - K                                                    
         L = JJ + 1                                                     
         SUM = ZERO                                                     
         IF (N .LT. JP1) GO TO 20                                       
         DO 10 I = JP1, N                                               
            SUM = SUM + R(L)*X(I)                                       
            L = L + 1                                                   
   10       CONTINUE                                                    
   20    CONTINUE                                                       
         TEMP = R(JJ)                                                   
         IF (TEMP .NE. ZERO) GO TO 40                                   
         L = J                                                          
         DO 30 I = 1, J                                                 
            TEMP = MAX(TEMP,ABS(R(L)))                                  
            L = L + N - I                                               
   30       CONTINUE                                                    
         TEMP = EPSMCH*TEMP                                             
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH                              
   40    CONTINUE                                                       
         X(J) = (QTB(J) - SUM)/TEMP                                     
   50    CONTINUE                                                       
      DO 60 J = 1, N                                                    
         WA1(J) = ZERO                                                  
         WA2(J) = DIAG(J)*X(J)                                          
   60    CONTINUE                                                       
      QNORM = ENORM(N,WA2)                                              
      IF (QNORM .LE. DELTA) GO TO 140                                   
      L = 1                                                             
      DO 80 J = 1, N                                                    
         TEMP = QTB(J)                                                  
         DO 70 I = J, N                                                 
            WA1(I) = WA1(I) + R(L)*TEMP                                 
            L = L + 1                                                   
   70       CONTINUE                                                    
         WA1(J) = WA1(J)/DIAG(J)                                        
   80    CONTINUE                                                       
      GNORM = ENORM(N,WA1)                                              
      SGNORM = ZERO                                                     
      ALPHA = DELTA/QNORM                                               
      IF (GNORM .EQ. ZERO) GO TO 120                                    
      DO 90 J = 1, N                                                    
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)                                
   90    CONTINUE                                                       
      L = 1                                                             
      DO 110 J = 1, N                                                   
         SUM = ZERO                                                     
         DO 100 I = J, N                                                
            SUM = SUM + R(L)*WA1(I)                                     
            L = L + 1                                                   
  100       CONTINUE                                                    
         WA2(J) = SUM                                                   
  110    CONTINUE                                                       
      TEMP = ENORM(N,WA2)                                               
      SGNORM = (GNORM/TEMP)/TEMP                                        
      ALPHA = ZERO                                                      
      IF (SGNORM .GE. DELTA) GO TO 120                                  
      BNORM = ENORM(N,QTB)                                              
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)                 
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2                     
     *       + SQRT((TEMP-(DELTA/QNORM))**2                             
     *              +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))    
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP            
  120 CONTINUE                                                          
      TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)                            
      DO 130 J = 1, N                                                   
         X(J) = TEMP*WA1(J) + ALPHA*X(J)                                
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END 
C===DOUBLE PRECISION FUNCTION SPMPAR(I)                                                               
      DOUBLE PRECISION FUNCTION SPMPAR(I)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      INTEGER MCHEPS(4)                                                 
      INTEGER MINMAG(4)                                                 
      INTEGER MAXMAG(4)                                                 
      DIMENSION RMACH(3)                                                
      EQUIVALENCE (RMACH(1),MCHEPS(1))                                  
      EQUIVALENCE (RMACH(2),MINMAG(1))                                  
      EQUIVALENCE (RMACH(3),MAXMAG(1))                                  
      IF ( I .EQ. 1 ) THEN                                              
      SPMPAR = 1.1102230246252D-16                                      
      ELSE IF ( I .EQ. 2 ) THEN                                         
      SPMPAR = 2.22507385850721D-308                                     
      ELSE IF ( I .EQ. 3 ) THEN                                         
      SPMPAR = 1.7976931348623D+308                                     
      ELSE                                                              
      WRITE(*,*) 'WRONG ARGUMENT FOR SPMPAR'                            
      END IF                                                            
      RETURN                                                            
      END
C===DOUBLE PRECISION FUNCTION ENORM(N,X)                                                                
      DOUBLE PRECISION FUNCTION ENORM(N,X)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N)                                                    
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/      
      S1 = ZERO                                                         
      S2 = ZERO                                                         
      S3 = ZERO                                                         
      X1MAX = ZERO                                                      
      X3MAX = ZERO                                                      
      FLOATN = N                                                        
      AGIANT = RGIANT/FLOATN                                            
      DO 90 I = 1, N                                                    
         XABS = ABS(X(I))                                               
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70          
            IF (XABS .LE. RDWARF) GO TO 30                              
               IF (XABS .LE. X1MAX) GO TO 10                            
                  S1 = ONE + S1*(X1MAX/XABS)**2                         
                  X1MAX = XABS                                          
                  GO TO 20                                              
   10          CONTINUE                                                 
                  S1 = S1 + (XABS/X1MAX)**2                             
   20          CONTINUE                                                 
               GO TO 60                                                 
   30       CONTINUE                                                    
               IF (XABS .LE. X3MAX) GO TO 40                            
                  S3 = ONE + S3*(X3MAX/XABS)**2                         
                  X3MAX = XABS                                          
                  GO TO 50                                              
   40          CONTINUE                                                 
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2         
   50          CONTINUE                                                 
   60       CONTINUE                                                    
            GO TO 80                                                    
   70    CONTINUE                                                       
            S2 = S2 + XABS**2                                           
   80    CONTINUE                                                       
   90    CONTINUE                                                       
      IF (S1 .EQ. ZERO) GO TO 100                                       
         ENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)                        
         GO TO 130                                                      
  100 CONTINUE                                                          
         IF (S2 .EQ. ZERO) GO TO 110                                    
            IF (S2 .GE. X3MAX)                                          
     *         ENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))             
            IF (S2 .LT. X3MAX)                                          
     *         ENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))              
            GO TO 120                                                   
  110    CONTINUE                                                       
            ENORM = X3MAX*SQRT(S3)                                      
  120    CONTINUE                                                       
  130 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE FDJAC1                                                               
      SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,    
     *                  WA1,WA2)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),FVEC(N),FJAC(LDFJAC,N),WA1(N),WA2(N)               
      DATA ZERO /0.0D0/                                                 
      EPSMCH = SPMPAR(1)                                                
      EPS = SQRT(MAX(EPSFCN,EPSMCH))                                    
      MSUM = ML + MU + 1                                                
      IF (MSUM .LT. N) GO TO 40                                         
         DO 20 J = 1, N                                                 
            TEMP = X(J)                                                 
            H = EPS*ABS(TEMP)                                           
            IF (H .EQ. ZERO) H = EPS                                    
            X(J) = TEMP + H                                             
            CALL FCN(N,X,WA1,IFLAG)                                     
            IF (IFLAG .LT. 0) GO TO 30                                  
            X(J) = TEMP                                                 
            DO 10 I = 1, N                                              
               FJAC(I,J) = (WA1(I) - FVEC(I))/H                         
   10          CONTINUE                                                 
   20       CONTINUE                                                    
   30    CONTINUE                                                       
         GO TO 110                                                      
   40 CONTINUE                                                          
         DO 90 K = 1, MSUM                                              
            DO 60 J = K, N, MSUM                                        
               WA2(J) = X(J)                                            
               H = EPS*ABS(WA2(J))                                      
               IF (H .EQ. ZERO) H = EPS                                 
               X(J) = WA2(J) + H                                        
   60          CONTINUE                                                 
            CALL FCN(N,X,WA1,IFLAG)                                     
            IF (IFLAG .LT. 0) GO TO 100                                 
            DO 80 J = K, N, MSUM                                        
               X(J) = WA2(J)                                            
               H = EPS*ABS(WA2(J))                                      
               IF (H .EQ. ZERO) H = EPS                                 
               DO 70 I = 1, N                                           
                  FJAC(I,J) = ZERO                                      
                  IF (I .GE. J - MU .AND. I .LE. J + ML)                
     *               FJAC(I,J) = (WA1(I) - FVEC(I))/H                   
   70             CONTINUE                                              
   80          CONTINUE                                                 
   90       CONTINUE                                                    
  100    CONTINUE                                                       
  110 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE QFORM(M,N,Q,LDQ,WA)                                                                
      SUBROUTINE QFORM(M,N,Q,LDQ,WA)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION Q(LDQ,M),WA(M)                                          
      DATA ONE,ZERO /1.0D0,0.0D0/                                       
      MINMN = MIN0(M,N)                                                 
      IF (MINMN .LT. 2) GO TO 30                                        
      DO 20 J = 2, MINMN                                                
         JM1 = J - 1                                                    
         DO 10 I = 1, JM1                                               
            Q(I,J) = ZERO                                               
   10       CONTINUE                                                    
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      NP1 = N + 1                                                       
      IF (M .LT. NP1) GO TO 60                                          
      DO 50 J = NP1, M                                                  
         DO 40 I = 1, M                                                 
            Q(I,J) = ZERO                                               
   40       CONTINUE                                                    
         Q(J,J) = ONE                                                   
   50    CONTINUE                                                       
   60 CONTINUE                                                          
      DO 120 L = 1, MINMN                                               
         K = MINMN - L + 1                                              
         DO 70 I = K, M                                                 
            WA(I) = Q(I,K)                                              
            Q(I,K) = ZERO                                               
   70       CONTINUE                                                    
         Q(K,K) = ONE                                                   
         IF (WA(K) .EQ. ZERO) GO TO 110                                 
         DO 100 J = K, M                                                
            SUM = ZERO                                                  
            DO 80 I = K, M                                              
               SUM = SUM + Q(I,J)*WA(I)                                 
   80          CONTINUE                                                 
            TEMP = SUM/WA(K)                                            
            DO 90 I = K, M                                              
               Q(I,J) = Q(I,J) - TEMP*WA(I)                             
   90          CONTINUE                                                 
  100       CONTINUE                                                    
  110    CONTINUE                                                       
  120    CONTINUE                                                       
      RETURN                                                            
      END
C===SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)                                                                
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION IPVT(LIPVT)                                             
      LOGICAL PIVOT                                                     
      DIMENSION A(LDA,N),RDIAG(N),ACNORM(N),WA(N)                       
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/                            
      EPSMCH = SPMPAR(1)                                                
      DO 10 J = 1, N                                                    
         ACNORM(J) = ENORM(M,A(1,J))                                    
         RDIAG(J) = ACNORM(J)                                           
         WA(J) = RDIAG(J)                                               
         IF (PIVOT) IPVT(J) = J                                         
   10    CONTINUE                                                       
      MINMN = MIN0(M,N)                                                 
      DO 110 J = 1, MINMN                                               
         IF (.NOT.PIVOT) GO TO 40                                       
         KMAX = J                                                       
         DO 20 K = J, N                                                 
            IF (RDIAG(K) .GT. RDIAG(KMAX)) KMAX = K                     
   20       CONTINUE                                                    
         IF (KMAX .EQ. J) GO TO 40                                      
         DO 30 I = 1, M                                                 
            TEMP = A(I,J)                                               
            A(I,J) = A(I,KMAX)                                          
            A(I,KMAX) = TEMP                                            
   30       CONTINUE                                                    
         RDIAG(KMAX) = RDIAG(J)                                         
         WA(KMAX) = WA(J)                                               
         K = IPVT(J)                                                    
         IPVT(J) = IPVT(KMAX)                                           
         IPVT(KMAX) = K                                                 
   40    CONTINUE                                                       
         AJNORM = ENORM(M-J+1,A(J,J))                                   
         IF (AJNORM .EQ. ZERO) GO TO 100                                
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM                         
         DO 50 I = J, M                                                 
            A(I,J) = A(I,J)/AJNORM                                      
   50       CONTINUE                                                    
         A(J,J) = A(J,J) + ONE                                          
         JP1 = J + 1                                                    
         IF (N .LT. JP1) GO TO 100                                      
         DO 90 K = JP1, N                                               
            SUM = ZERO                                                  
            DO 60 I = J, M                                              
               SUM = SUM + A(I,J)*A(I,K)                                
   60          CONTINUE                                                 
            TEMP = SUM/A(J,J)                                           
            DO 70 I = J, M                                              
               A(I,K) = A(I,K) - TEMP*A(I,J)                            
   70          CONTINUE                                                 
            IF (.NOT.PIVOT .OR. RDIAG(K) .EQ. ZERO) GO TO 80            
            TEMP = A(J,K)/RDIAG(K)                                      
            RDIAG(K) = RDIAG(K)*SQRT(MAX(ZERO,ONE-TEMP**2))             
            IF (P05*(RDIAG(K)/WA(K))**2 .GT. EPSMCH) GO TO 80           
            RDIAG(K) = ENORM(M-J,A(JP1,K))                              
            WA(K) = RDIAG(K)                                            
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
         RDIAG(J) = -AJNORM                                             
  110    CONTINUE                                                       
      RETURN                                                            
      END
C===SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)                                                               
      SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION A(LDA,N),V(N),W(N)                                      
      DATA ONE /1.0D0/                                                  
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 50                                          
      DO 20 NMJ = 1, NM1                                                
         J = N - NMJ                                                    
         IF (ABS(V(J)) .GT. ONE) COS = ONE/V(J)                         
         IF (ABS(V(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)                 
         IF (ABS(V(J)) .LE. ONE) SIN = V(J)                             
         IF (ABS(V(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)                 
         DO 10 I = 1, M                                                 
            TEMP = COS*A(I,J) - SIN*A(I,N)                              
            A(I,N) = SIN*A(I,J) + COS*A(I,N)                            
            A(I,J) = TEMP                                               
   10       CONTINUE                                                    
   20    CONTINUE                                                       
      DO 40 J = 1, NM1                                                  
         IF (ABS(W(J)) .GT. ONE) COS = ONE/W(J)                         
         IF (ABS(W(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)                 
         IF (ABS(W(J)) .LE. ONE) SIN = W(J)                             
         IF (ABS(W(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)                 
         DO 30 I = 1, M                                                 
            TEMP = COS*A(I,J) + SIN*A(I,N)                              
            A(I,N) = -SIN*A(I,J) + COS*A(I,N)                           
            A(I,J) = TEMP                                               
   30       CONTINUE                                                    
   40    CONTINUE                                                       
   50 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)                                                               
      SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL SING                                                      
      DIMENSION S(LS),U(M),V(N),W(M)                                    
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/                  
      GIANT = SPMPAR(3)                                                 
      JJ = (N*(2*M - N + 1))/2 - (M - N)                                
      L = JJ                                                            
      DO 10 I = N, M                                                    
         W(I) = S(L)                                                    
         L = L + 1                                                      
   10    CONTINUE                                                       
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 70                                          
      DO 60 NMJ = 1, NM1                                                
         J = N - NMJ                                                    
         JJ = JJ - (M - J + 1)                                          
         W(J) = ZERO                                                    
         IF (V(J) .EQ. ZERO) GO TO 50                                   
         IF (ABS(V(N)) .GE. ABS(V(J))) GO TO 20                         
            COTAN = V(N)/V(J)                                           
            SIN = P5/SQRT(P25+P25*COTAN**2)                             
            COS = SIN*COTAN                                             
            TAU = ONE                                                   
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS                  
            GO TO 30                                                    
   20    CONTINUE                                                       
            TAN = V(J)/V(N)                                             
            COS = P5/SQRT(P25+P25*TAN**2)                               
            SIN = COS*TAN                                               
            TAU = SIN                                                   
   30    CONTINUE                                                       
         V(N) = SIN*V(J) + COS*V(N)                                     
         V(J) = TAU                                                     
         L = JJ                                                         
         DO 40 I = J, M                                                 
            TEMP = COS*S(L) - SIN*W(I)                                  
            W(I) = SIN*S(L) + COS*W(I)                                  
            S(L) = TEMP                                                 
            L = L + 1                                                   
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
      DO 80 I = 1, M                                                    
         W(I) = W(I) + V(N)*U(I)                                        
   80    CONTINUE                                                       
      SING = .FALSE.                                                    
      IF (NM1 .LT. 1) GO TO 140                                         
      DO 130 J = 1, NM1                                                 
         IF (W(J) .EQ. ZERO) GO TO 120                                  
         IF (ABS(S(JJ)) .GE. ABS(W(J))) GO TO 90                        
            COTAN = S(JJ)/W(J)                                          
            SIN = P5/SQRT(P25+P25*COTAN**2)                             
            COS = SIN*COTAN                                             
            TAU = ONE                                                   
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS                  
            GO TO 100                                                   
   90    CONTINUE                                                       
            TAN = W(J)/S(JJ)                                            
            COS = P5/SQRT(P25+P25*TAN**2)                               
            SIN = COS*TAN                                               
            TAU = SIN                                                   
  100    CONTINUE                                                       
         L = JJ                                                         
         DO 110 I = J, M                                                
            TEMP = COS*S(L) + SIN*W(I)                                  
            W(I) = -SIN*S(L) + COS*W(I)                                 
            S(L) = TEMP                                                 
            L = L + 1                                                   
  110       CONTINUE                                                    
         W(J) = TAU                                                     
  120    CONTINUE                                                       
         IF (S(JJ) .EQ. ZERO) SING = .TRUE.                             
         JJ = JJ + (M - J + 1)                                          
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      L = JJ                                                            
      DO 150 I = N, M                                                   
         S(L) = W(I)                                                    
         L = L + 1                                                      
  150    CONTINUE                                                       
      IF (S(JJ) .EQ. ZERO) SING = .TRUE.                                
      RETURN                                                            
      END
C===SUBROUTINE FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)                                                               
      SUBROUTINE FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),FVEC(M),FJAC(LDFJAC,N),WA(M)                       
      DATA ZERO /0.0D0/                                                 
      EPSMCH = SPMPAR(1)                                                
      EPS = SQRT(MAX(EPSFCN,EPSMCH))                                    
      DO 20 J = 1, N                                                    
         TEMP = X(J)                                                    
         H = EPS*ABS(TEMP)                                              
         IF (H .EQ. ZERO) H = EPS                                       
         X(J) = TEMP + H                                                
         CALL FCN(M,N,X,WA,IFLAG)                                       
         IF (IFLAG .LT. 0) GO TO 30                                     
         X(J) = TEMP                                                    
         DO 10 I = 1, M                                                 
            FJAC(I,J) = (WA(I) - FVEC(I))/H                             
   10       CONTINUE                                                    
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      RETURN                                                            
      END
C===SUBROUTINE LMPAR                                                               
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,WA1,     
     *                 WA2)                                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      INTEGER IPVT(N)                                                   
      DIMENSION R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA1(N),WA2(N)     
      DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/                           
      DWARF = SPMPAR(2)                                                 
      NSING = N                                                         
      DO 10 J = 1, N                                                    
         WA1(J) = QTB(J)                                                
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1         
         IF (NSING .LT. N) WA1(J) = ZERO                                
   10    CONTINUE                                                       
      IF (NSING .LT. 1) GO TO 50                                        
      DO 40 K = 1, NSING                                                
         J = NSING - K + 1                                              
         WA1(J) = WA1(J)/R(J,J)                                         
         TEMP = WA1(J)                                                  
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 30                                       
         DO 20 I = 1, JM1                                               
            WA1(I) = WA1(I) - R(I,J)*TEMP                               
   20       CONTINUE                                                    
   30    CONTINUE                                                       
   40    CONTINUE                                                       
   50 CONTINUE                                                          
      DO 60 J = 1, N                                                    
         L = IPVT(J)                                                    
         X(L) = WA1(J)                                                  
   60    CONTINUE                                                       
      ITER = 0                                                          
      DO 70 J = 1, N                                                    
         WA2(J) = DIAG(J)*X(J)                                          
   70    CONTINUE                                                       
      DXNORM = ENORM(N,WA2)                                             
      FP = DXNORM - DELTA                                               
      IF (FP .LE. P1*DELTA) GO TO 220                                   
      PARL = ZERO                                                       
      IF (NSING .LT. N) GO TO 120                                       
      DO 80 J = 1, N                                                    
         L = IPVT(J)                                                    
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)                               
   80    CONTINUE                                                       
      DO 110 J = 1, N                                                   
         SUM = ZERO                                                     
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 100                                      
         DO 90 I = 1, JM1                                               
            SUM = SUM + R(I,J)*WA1(I)                                   
   90       CONTINUE                                                    
  100    CONTINUE                                                       
         WA1(J) = (WA1(J) - SUM)/R(J,J)                                 
  110    CONTINUE                                                       
      TEMP = ENORM(N,WA1)                                               
      PARL = ((FP/DELTA)/TEMP)/TEMP                                     
  120 CONTINUE                                                          
      DO 140 J = 1, N                                                   
         SUM = ZERO                                                     
         DO 130 I = 1, J                                                
            SUM = SUM + R(I,J)*QTB(I)                                   
  130       CONTINUE                                                    
         L = IPVT(J)                                                    
         WA1(J) = SUM/DIAG(L)                                           
  140    CONTINUE                                                       
      GNORM = ENORM(N,WA1)                                              
      PARU = GNORM/DELTA                                                
      IF (PARU .EQ. ZERO) PARU = DWARF/MIN(DELTA,P1)                    
      PAR = MAX(PAR,PARL)                                               
      PAR = MIN(PAR,PARU)                                               
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM                             
  150 CONTINUE                                                          
         ITER = ITER + 1                                                
         IF (PAR .EQ. ZERO) PAR = MAX(DWARF,P001*PARU)                  
         TEMP = SQRT(PAR)                                               
         DO 160 J = 1, N                                                
            WA1(J) = TEMP*DIAG(J)                                       
  160       CONTINUE                                                    
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SDIAG,WA2)                  
         DO 170 J = 1, N                                                
            WA2(J) = DIAG(J)*X(J)                                       
  170       CONTINUE                                                    
         DXNORM = ENORM(N,WA2)                                          
         TEMP = FP                                                      
         FP = DXNORM - DELTA                                            
         IF (ABS(FP) .LE. P1*DELTA                                      
     *       .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP                     
     *            .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220     
         DO 180 J = 1, N                                                
            L = IPVT(J)                                                 
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)                            
  180       CONTINUE                                                    
         DO 210 J = 1, N                                                
            WA1(J) = WA1(J)/SDIAG(J)                                    
            TEMP = WA1(J)                                               
            JP1 = J + 1                                                 
            IF (N .LT. JP1) GO TO 200                                   
            DO 190 I = JP1, N                                           
               WA1(I) = WA1(I) - R(I,J)*TEMP                            
  190          CONTINUE                                                 
  200       CONTINUE                                                    
  210       CONTINUE                                                    
         TEMP = ENORM(N,WA1)                                            
         PARC = ((FP/DELTA)/TEMP)/TEMP                                  
         IF (FP .GT. ZERO) PARL = MAX(PARL,PAR)                         
         IF (FP .LT. ZERO) PARU = MIN(PARU,PAR)                         
         PAR = MAX(PARL,PAR+PARC)                                       
         GO TO 150                                                      
  220 CONTINUE                                                          
      IF (ITER .EQ. 0) PAR = ZERO                                       
      RETURN                                                            
      END 
C===SUBROUTINE QRSOLV                                                              
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      INTEGER IPVT(N)                                                   
      DIMENSION R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA(N)             
      DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/                            
      DO 20 J = 1, N                                                    
         DO 10 I = J, N                                                 
            R(I,J) = R(J,I)                                             
   10       CONTINUE                                                    
         X(J) = R(J,J)                                                  
         WA(J) = QTB(J)                                                 
   20    CONTINUE                                                       
      DO 100 J = 1, N                                                   
         L = IPVT(J)                                                    
         IF (DIAG(L) .EQ. ZERO) GO TO 90                                
         DO 30 K = J, N                                                 
            SDIAG(K) = ZERO                                             
   30       CONTINUE                                                    
         SDIAG(J) = DIAG(L)                                             
         QTBPJ = ZERO                                                   
         DO 80 K = J, N                                                 
            IF (SDIAG(K) .EQ. ZERO) GO TO 70                            
            IF (ABS(R(K,K)) .GE. ABS(SDIAG(K))) GO TO 40                
               COTAN = R(K,K)/SDIAG(K)                                  
               SIN = P5/SQRT(P25+P25*COTAN**2)                          
               COS = SIN*COTAN                                          
               GO TO 50                                                 
   40       CONTINUE                                                    
               TAN = SDIAG(K)/R(K,K)                                    
               COS = P5/SQRT(P25+P25*TAN**2)                            
               SIN = COS*TAN                                            
   50       CONTINUE                                                    
            R(K,K) = COS*R(K,K) + SIN*SDIAG(K)                          
            TEMP = COS*WA(K) + SIN*QTBPJ                                
            QTBPJ = -SIN*WA(K) + COS*QTBPJ                              
            WA(K) = TEMP                                                
            KP1 = K + 1                                                 
            IF (N .LT. KP1) GO TO 70                                    
            DO 60 I = KP1, N                                            
               TEMP = COS*R(I,K) + SIN*SDIAG(I)                         
               SDIAG(I) = -SIN*R(I,K) + COS*SDIAG(I)                    
               R(I,K) = TEMP                                            
   60          CONTINUE                                                 
   70       CONTINUE                                                    
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         SDIAG(J) = R(J,J)                                              
         R(J,J) = X(J)                                                  
  100    CONTINUE                                                       
      NSING = N                                                         
      DO 110 J = 1, N                                                   
         IF (SDIAG(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1       
         IF (NSING .LT. N) WA(J) = ZERO                                 
  110    CONTINUE                                                       
      IF (NSING .LT. 1) GO TO 150                                       
      DO 140 K = 1, NSING                                               
         J = NSING - K + 1                                              
         SUM = ZERO                                                     
         JP1 = J + 1                                                    
         IF (NSING .LT. JP1) GO TO 130                                  
         DO 120 I = JP1, NSING                                          
            SUM = SUM + R(I,J)*WA(I)                                    
  120       CONTINUE                                                    
  130    CONTINUE                                                       
         WA(J) = (WA(J) - SUM)/SDIAG(J)                                 
  140    CONTINUE                                                       
  150 CONTINUE                                                          
      DO 160 J = 1, N                                                   
         L = IPVT(J)                                                    
         X(L) = WA(J)                                                   
  160    CONTINUE                                                       
      RETURN                                                            
      END
C===DOUBLE PRECISION FUNCTION XPN(X)                                                               
      DOUBLE PRECISION FUNCTION XPN(X)                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      XMIN  = -75.0D0                                                   
      XMAX  =  75.0D0                                                   
      IF  ( X .LT. XMIN )  X = XMIN                                     
      IF  ( X .GT. XMAX )  X = XMAX                                     
      XPN   =  EXP(X)                                                   
      RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION YLOG(X)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      IF ( X .LT. 0.000001D0 ) X = 0.000001D0                           
      YLOG = LOG(X)                                                     
      RETURN                                                            
      END
C===SUBROUTINE NEWRES                                                               
      SUBROUTINE NEWRES                                                 
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z )                             
      CHARACTER*8   NAME, IXOG                                          
      LOGICAL       LFIX2                                               
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),A(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,IPER,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/BASDAT/ XB(600,300)                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      COMMON/SIMDIM/NAG,LPER,MAXE,MAXX,MAXP                             
      COMMON/FIX2/ LFIX2, NFIX2, LDOFIX, NFIXMX, NFIXEND                
      DO 9000 IC = 1,2                                                
      NV  =  KK(IC+1) - KK(IC)                                          
      DO 8000 IV = 1,NV                                                 
      LPRINT = 0                                                        
      KVAR = KK(IC) + IV                                                
      DO 6000 IY = 1,KAG                                                
      V  =  ER(KVAR,IY)                                                 
      IF ( ABS(V) .GT. 0.000001 ) THEN                                  
        LPRINT  =  1                                                    
      END IF                                                            
 6000 CONTINUE
      IF ( LPRINT .EQ. 1 ) THEN                                         
      WRITE(8,7000) IC,IV                                               
      WRITE(8,7500) ( ER(KVAR,IY),IY=1,147 )                            
      END IF                                                            
 7000 FORMAT(2I2)                                                       
 7500 FORMAT(4F18.12)                                                    
 8000 CONTINUE                                                          
 9000 CONTINUE                                                          
      RETURN                                                            
      END 
C===SUBROUTINE SWUS   
      SUBROUTINE SWUS                                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      CHARACTER*8   NAME, IXOG                                          
      LOGICAL       LFIX2                                               
      COMMON/ALPHA/ NAME(300), IXOG(600)                                
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          DUMMY(1500,300)
      COMMON/HELEN/V(300),Z(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,XX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR   
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/K01,K02,K03,K04,K05,K06,K07,K08,K09,K10,K11,K12,K13,
     1K14,K15,K16,K17,K18,K19,K20,K21,K22,K23,K24,K25,K26,K27,K28,K29,  
     2K30,K31,K32,K33,K34,K35,K36,K37,K38,K39,K40,K41                   
      COMMON/JGROUP/J01,J02,J03,J04,J05,J06,J07,J08,J09,J10,J11,J12,J13,
     1J14,J15,J16,J17,J18,J19,J20,J21,J22,J23,J24,J25,J26,J27,J28,J29,  
     2J30,J31,J32,J33,J34,J35,J36,J37,J38,J39,J40,J41                   
      COMMON/FIX2/ LFIX2, NFIX2, LDOFIX, NFIXMX, NFIXEND   
      COMMON/SHK/err(10,1000)
      double precision zerr(10,1000)
      PARAMETER( ZERO  = 0.00D0,
     &           P03   = 0.03D0,
     &           P25   = 0.25D0,
     &           P30   = 0.30D0,
     &           P50   = 0.50D0,
     &           P75   = 0.75D0,
     &           ONE   = 1.00D0,
     &           TWO   = 2.00D0 )
      DV(I,J,K)=VV(I,J,K)-VV(I,J,K-1)                                   
      DL(I,J,K)=VV(I,J,K)/VV(I,J,K-1)      
      CALPH  = Z(1)
      CBETA  = Z(2)
      CDELTA = Z(3)
      CMIU1  = Z(4)
      CMIU2  = ONE-CMIU1
      COMGY1 = Z(5)
      COMGY2 = Z(6)
      COMGK1 = Z(7)
      COMGK2 = Z(8)
      COMGC1 = Z(9)
      COMGC2 = Z(10)
      CPSI1  = Z(11)
      CPSI2  = Z(12)
      CYVC   = Z(13)
      CKVC   = Z(14)
      CLABOR1= Z(15)
      CLABOR2= Z(16)
      CPI1   = Z(17)
      CPI2   = Z(18)
      CEQR   = Z(19)
      CEQT   = Z(20)  
      CTHETA = Z(21)
      CRHO1  = Z(22)
      CRHO2  = Z(23)
      CPHIB1 = Z(24)
      CPHIB2 = Z(25)
      CRHO3  = Z(26)
      CETA1  = (ONE-CEQT)*CALPH*(CYVC/CKVC)/(CDELTA+CEQR)	  
      CETA2  = ONE/(CDELTA+CEQR)-CETA1/CPSI1
      if(i.eq.mag) then
      do 3110, k=1,10
      do 3111, j=1,147
         zerr(k,j)=err(k,j)
 3111 continue
 3110 continue
      else
      do 2110, k=1,10
      do 2111, j=1,147
         zerr(k,j)=0.0d0
 2111    continue
 2110 continue
      endif
      X(17,I)=VV(21,17,0)
      X(18,I)=VV(21,18,-1)+Z(31)*(VV(21,18,-1)
     >       -VV(21,18,-2))+zerr(1,i)
      X(19,I)=Z(32)*VV(21,19,-1)+zerr(2,i)
      X(20,I)=Z(33)*VV(21,20,-1)+zerr(3,i)
      X(21,I)=Z(34)*VV(21,21,-1)+zerr(4,i)
      X(22,I)=Z(35)*VV(21,22,-1)+zerr(5,i)	  
      X(23,I)=Z(36)*VV(21,23,-1)+zerr(6,i)
      X(24,I)=Z(37)*VV(21,24,-1)+zerr(7,i)
      X(25,I)=Z(38)*VV(21,25,-1)+zerr(8,i)
      X(26,I)=Z(39)*VV(21,26,-1)+zerr(9,i)
      X(27,I)=Z(40)*VV(21,27,-1)+zerr(10,i)
      X(28,I)= ZERO
      X(29,I)= ZERO
      X(30,I)= ZERO
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(1,I) =VV(21,1,0)-CPSI1*(VV(21,32,0)-VV(1,10,0))
      ENDIF  
      V(1) =CPSI1*(VV(21,32,0)-VV(1,10,0))
     > +RR(1,1,0) 
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(2,I) =VV(21,2,0)-(COMGY1*VV(1,5,0)+COMGY2*VV(1,6,0))
     > -VV(21,20,0)	 
      ENDIF
      V(2) =COMGY1*VV(1,5,0)+COMGY2*VV(1,6,0)
     > +VV(21,20,0)	 
     > +RR(1,2,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(3,I) =VV(21,3,0)-(COMGK1*VV(1,7,0)+COMGK2*VV(1,8,0))
     > -VV(21,21,0)	 
      ENDIF
      V(3) =COMGK1*VV(1,7,0)+COMGK2*VV(1,8,0)
     > +VV(21,21,0)	 
     > +RR(1,3,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(4,I)= VV(21,4,0)-(CYVC*((ONE-CEQT)*VV(1,2,0)
     > -VV(21,17,0))-CKVC*(VV(1,3,0)-(ONE-CDELTA)*VV(1,3,-1)))
     > -VV(21,22,0)	 
      ENDIF
      V(4)= CYVC*((ONE-CEQT)*VV(1,2,0)-VV(21,17,0))
     > -CKVC*(VV(1,3,0)-(ONE-CDELTA)*VV(1,3,-1))
     > +VV(21,22,0)	 
     > +RR(1,4,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(5,I)= VV(21,5,0)-CALPH*VV(1,7,-1)-(ONE-CALPH)*VV(1,11,0)
     > -VV(1,15,-1)	  
      ENDIF
      V(5)= CALPH*VV(1,7,-1)+(ONE-CALPH)*VV(1,11,0)+VV(1,15,-1)
     > +RR(1,5,0) 
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(6,I)= VV(21,6,0)-CALPH*VV(1,8,-1)-(ONE-CALPH)*VV(1,12,0)
     > -VV(1,16,-1)	  
      ENDIF
      V(6)= CALPH*VV(1,8,-1)+(ONE-CALPH)*VV(1,12,0)+VV(1,16,-1)
     > +RR(1,6,0)	 
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(7,I)= VV(21,7,0)-(VV(1,5,0)+VV(1,1,0)/CPSI1-CKVC
     > /CYVC*(VV(1,1,0)+CDELTA)/CALPH/(ONE-VV(21,17,0)))
      ENDIF	  
      V(7)= VV(1,5,0)+VV(1,1,0)/CPSI1-CKVC/CYVC*(VV(1,1,0)
     > +CDELTA)/CALPH/(ONE-VV(21,17,0))
     > +RR(1,7,0)	 
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(8,I)= VV(21,8,0)-(VV(1,6,0)+VV(1,1,0)/CPSI1-CKVC
     > /CYVC*(VV(1,1,0)+CDELTA)/CALPH/(ONE-VV(21,17,0)))
      ENDIF	  
      V(8)= VV(1,6,0)+VV(1,1,0)/CPSI1-CKVC/CYVC*(VV(1,1,0)
     > +CDELTA)/CALPH/(ONE-VV(21,17,0))
     > +RR(1,8,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(9,I)= VV(21,9,0)-(VV(21,31,0)-VV(1,1,0)/CPSI1)
     > -VV(21,23,0)
      ENDIF
      V(9)= VV(21,31,0)-VV(1,1,0)/CPSI1
     > +VV(21,23,0)
     > +RR(1,9,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(10,I)= VV(21,10,0)-(VV(1,4,0)-COMGC1*VV(1,9,0))/COMGC2
     > -VV(21,24,0)
      ENDIF
      V(10)= (VV(1,4,0)-COMGC1*VV(1,9,0))/COMGC2
     > +VV(21,24,0)
     > +RR(1,10,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(11,I)=VV(21,11,0)-(VV(1,5,0)-CPSI1*VV(1,9,0)-(TWO
     > *CPSI2/CTHETA)*(VV(1,15,0)-VV(1,15,-1)))/(ONE+CPSI2)
     > -VV(21,25,0)
      ENDIF    
      V(11)=(VV(1,5,0)-CPSI1*VV(1,9,0)-(TWO*CPSI2/CTHETA)
     > *(VV(1,15,0)-VV(1,15,-1)))/(ONE+CPSI2)
     > +VV(21,25,0)
     > +RR(1,11,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(12,I)= VV(21,12,0)-(VV(1,6,0)-CPSI1*VV(1,10,0)-(TWO
     > *CPSI2/CTHETA)*(VV(1,16,0)-VV(1,16,-1)))/(ONE+CPSI2)
     > -VV(21,26,0)
      ENDIF    
      V(12)=(VV(1,6,0)-CPSI1*VV(1,10,0)-(TWO*CPSI2/CTHETA)
     > *(VV(1,16,0)-VV(1,16,-1)))/(ONE+CPSI2)	 
     > +VV(21,26,0)
     > +RR(1,12,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(13,I)= VV(21,13,0)-(CRHO1*VV(1,13,-1)-CRHO2*CMIU1
     > /COMGY1*(exp(VV(1,7,-2)-VV(1,3,-2))**TWO)) 
     > -VV(21,19,0)	 
      ENDIF 
      V(13)= CRHO1*VV(1,13,-1)-CRHO2*CMIU1
     > /COMGY1*(exp(VV(1,7,-2)-VV(1,3,-2))**TWO) 
     > +VV(21,19,0)	 
     > +RR(1,13,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(14,I)= VV(21,14,0)-(CRHO1*VV(1,14,-1)-CRHO2*CMIU2
     > /COMGY2*(exp(VV(1,8,-2)-VV(1,3,-2))**TWO)) 
     > +CRHO3*VV(21,27,0)-VV(21,19,0)	 
      ENDIF 
      V(14)= CRHO1*VV(1,14,-1)-CRHO2*CMIU2
     > /COMGY2*(exp(VV(1,8,-2)-VV(1,3,-2))**TWO) 
     > -CRHO3*VV(21,27,0)+VV(21,19,0)
     > +RR(1,14,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(15,I)= VV(21,15,0)-VV(1,15,-1)+CPHIB1*(VV(1,13,0)
     > +VV(21,17,0)/(ONE-CEQT))
      ENDIF
      V(15)= VV(1,15,-1)-CPHIB1*(VV(1,13,0)+VV(21,17,0)
     > /(ONE-CEQT))+zerr(1,i)	 
     > +RR(1,15,0)
      IF ( LDOFIX .EQ. 1 ) THEN
      ER(16,I)= VV(21,16,0)-VV(1,16,-1)+CPHIB2*(VV(1,14,0)
     > +VV(21,17,0)/(ONE-CEQT))
      ENDIF 
      V(16)= VV(1,16,-1)-CPHIB2*(VV(1,14,0)+VV(21,17,0)
     > /(ONE-CEQT))+zerr(1,i)
     > +RR(1,16,0)
      RETURN                                                            
      END
C===SUBROUTINE CALCFX( FX )              
      SUBROUTINE CALCFX( FX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      COMMON/L2/E(300,300),   X(600,300),  ER(300,300),
     1          ERX(600,300), EE(300,300), XX(600,300)
      COMMON/HELEN/F(300),Z(600),BA(300),JX(10),VX(10),                 
     1 KAG,MAG,NDOG,LEXOG,IPAR,IRAND,B,NSTART,ITMX,PX,P,TOL,TOLR,BR     
     2,PR,NPER,NSK,IJX,I   ,IDATE,NEGP,T,IOPT,IDYN,ITR,NARG,BSTEP,MAXITR
     3,NXGP,NCG,IDYN1,IRANDX,IB,IRHO                                         
      COMMON/KGROUP/KK(41)/JGROUP/JJ(41)                                
      COMMON/WOREQN/IBLOC,NCONV,NBLOC                                   
      COMMON/NEWCAL/NEN(100),NEX(100),NENH(100),NEXP,MXEH,IOPTT,IOPTB,  
     1IOPTC                                                             
      COMMON/SIMDIM/NAG,LPER,MAXE,MAXX,MAXP                             
      COMMON/POWELL/R(200000),DIR(200000),L(200000),ICTRL                     
      COMMON/TERMCF/term_coef(2,4)   
      KX(I,J)=600*(J-1)+I                                               
      KVV(IG,J,K)=KX(KK(IG)+J,I+K)                                      
      J = 0
      ZERO   = 0.00D0
      P25    = 0.25D0
      P33    = 1.00D0/3.00D0
      P50    = 0.50D0
      ONE    = 1.00D0
      TWO    = 2.00D0	  	  
      RC1TERM=term_coef(1,1)*E(15,KAG-1)+term_coef(1,2)
     > *E(16,KAG-1)+term_coef(1,3)*X(27,KAG)
      RC2TERM=term_coef(2,1)*E(15,KAG-1)+term_coef(2,2)
     > *E(16,KAG-1)+term_coef(2,3)*X(27,KAG)
      MAGM1  = MAG - 1
      MAGM2  = MAG - 2
      MAGM3  = MAG - 3
      MAGM4  = MAG - 4
      KAGM1  = KAG - 1
      KAGM2  = KAG - 2
      KAGM3  = KAG - 3
      KAGM4  = KAG - 4
      KAGM5  = KAG - 5
      KAGM6  = KAG - 6
      DO 9815 I = MAG,KAG
      J = J + 1
      IF ( I .EQ. KAG ) THEN
      L(J) = KVV(21,31,0)
      R(J) = RC1TERM - VV(21,31,0)
      ELSE
      L(J) = KVV(21,31,0)
      R(J) = VV(1,9,1) - VV(21,31,0)
      END IF
 9815 CONTINUE
      DO 9816 I = MAG,KAG
      J = J + 1
      IF ( I .EQ. KAG ) THEN
      L(J) = KVV(21,32,0)
      R(J) = RC2TERM - VV(21,32,0)
      ELSE
      L(J) = KVV(21,32,0)
      R(J) = VV(1,10,1) - VV(21,32,0)
      END IF
 9816 CONTINUE
      I = MAG - 1
      J = J + 1
      L(J) = KVV(21,31,0)
      R(J) = VV(1,9,1) - VV(21,31,0)
      J = J + 1
      L(J) = KVV(21,32,0)
      R(J) = VV(1,10,1) - VV(21,32,0)
      RETURN
      END
      function gen_errname(i)
      implicit none
      integer i,i1,j,k
      character*10 digit
      character*18 gen_errname
      digit='0123456789'
      gen_errname(1:9)='tmp/errnc'
      gen_errname(14:18)='.data'
      i1=i
      do 10, j=13,10,-1
         k = i1 - int(i1/10)*10 + 1
         i1 = int(i1/10)
         gen_errname(j:j) = digit(k:k)
 10   continue
      return
      end
      function gen_xname(i)
      implicit none
      integer i,i1,j,k
      character*10 digit
      character*18 gen_xname
      digit='0123456789'
      gen_xname(1:9)='tmp/exonc'
      gen_xname(14:18)='.data'
      i1=i
      do 10, j=13,10,-1
         k = i1 - int(i1/10)*10 + 1
         i1 = int(i1/10)
         gen_xname(j:j) = digit(k:k)
 10   continue
      return
      end