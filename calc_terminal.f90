!!=========================================================================
!! This new version has one difference that the credit now has time trend, 
!! because the credit data now are replaced by "M4 stocdk / GDP"
!! So "credit" now will also affect the terminal condition
!!=========================================================================

 PROGRAM MAIN
 IMPLICIT NONE
 INTEGER i,j
!***************************** COEFFICIENTS ****************************
 DOUBLE PRECISION ALPHA,BETA,DELTA,MIU1,MIU2,OMGY1,OMGY2,OMGK1,OMGK2
 DOUBLE PRECISION OMGC1,OMGC2,PSI1,PSI2,YVC,KVC,LABOR1,LABOR2,PI1,PI2
 DOUBLE PRECISION EQR,EQT,THETA,RHO1,RHO2,PHIB1,PHIB2,PHI3,ETA1,ETA2
 DOUBLE PRECISION YVW1,YVW2,PPP1,PPP2,YVC1,KVC1 
!***********************************************************************
! Here I set Pi' in EQ is constant
 DOUBLE PRECISION,DIMENSION(30)::    coef
 DOUBLE PRECISION,DIMENSION(13,13):: coef_A, inv_A
! DOUBLE PRECISION,DIMENSION(13,2)::  coef_B, coef_TM1 
 DOUBLE PRECISION,DIMENSION(13,3)::  coef_B, coef_TM1
 DOUBLE PRECISION,DIMENSION(13,1)::  coef_D, coef_TM2
! DOUBLE PRECISION,DIMENSION(2,3)::   term_coef
 DOUBLE PRECISION,DIMENSION(2,4)::   term_coef
 
 open(10,file='coef.data',status='old')
     read(10,*)coef
 close(10)
!========================= Terminal conditions =========================
! Yang-- As we have too many equations, we read updated parameter values
! Yang-- and calculate terminal conditions below each time by A*X=B*Y+D
      ALPHA = coef(1)
      BETA  = coef(2)
      DELTA = coef(3)
      MIU1  = coef(4)
      MIU2  = 1.0D0-MIU1
      OMGY1 = coef(5)
      OMGY2 = coef(6)
      OMGK1 = coef(7)
      OMGK2 = coef(8)
      OMGC1 = coef(9)
      OMGC2 = coef(10)
      PSI1  = coef(11)
      PSI2  = coef(12)
      YVC   = coef(13)
      KVC   = coef(14)
      LABOR1= coef(15)
      LABOR2= coef(16)
      PI1   = coef(17)
      PI2   = coef(18)
      EQR   = coef(19)
      EQT   = coef(20)  
      THETA = coef(21)
      RHO1  = coef(22)
      RHO2  = coef(23)
      PHIB1 = coef(24)
      PHIB2 = coef(25)
      PHI3  = coef(26)

      ETA1  = (1.0D0-EQT)*ALPHA*(YVC/KVC)/(DELTA+EQR)	  
      ETA2  = 1.0D0/(DELTA+EQR)-ETA1/PSI1
  
      YVC1  = 1.710D0
      KVC1  = 5.431D0
  
!      YVW1  = LABOR1/(1.0D0-ALPHA)	  
!      PHIB1 =(1.0D0-EQT)*BETA*THETA*YVW1/(1.0D0-BETA)/((1.0D0-EQT+PI1)**2.0D0)
!      YVW2  = LABOR2/(1.0D0-ALPHA)	  
!      PHIB2 =(1.0D0-EQT)*BETA*THETA*YVW2/(1.0D0-BETA)/((1.0D0-EQT+PI2)**2.0D0)
	  
! Yang-- x: Endogenous VECTOR; y: Non-stationary Exogenous Vector; D: Constant Vector
! Yang--    1  2  3  4   5   6   7   8   9   10  11  12  13 
! Yang-- x=[Y; K; C; Y1; Y2; K1; K2; C1; C2; N1; N2; A1; A2]
! Yang-- y=[LA1;LA2;Crd] where Crd means Credit
 do i=1,13
    do j=1,13
       coef_A(i,j)=0.0D0
    end do
    do j=1,2
       coef_B(i,j)=0.0D0
    end do 
       coef_D(i,1)=0.0D0
 end do

 do i=1,13
    coef_A(i,i)=1.0D0
 end do 

!1  Y  = OMGY1*Y1 + OMGY2*Y2
!2  K  = OMGK1*K1 + OMGK2*K2
!3  C  = YVC*[(1-EQT)Y-EQT] - KVC*DELTA*K
!  Here I minus 0 instead of EQT because the exo stationary tax =EQT in steady state
!4  Y1 = ALPHA*K1 + (1-ALPHA)*N1 + LA1
!5  Y2 = ALPHA*K2 + (1-ALPHA)*N2 + LA2
!6  K1 = Y1 + 0 - KVC*[(0+DELTA)/ALPHA/(1-EQT)]
!7  K2 = Y2 + 0 - KVC*[(0+DELTA)/ALPHA/(1-EQT)]
!8  C1 = YVC*[(1-EQT)*Y1-EQT] - KVC*DELTA*K1
!9  C2 = (C-OMGC1*C1)/OMGC2
!10 N1 = [Y1-PSI1*C1-2*PSI2/THETA*(A1-LA1)]/(1+PSI2)
!11 N2 = [Y2-PSI1*C2-2*PSI2/THETA*(A2-LA2)]/(1+PSI2)
!12 A1 = LA1 - PHIB1*[PI1+EQT/(1-EQT)]
!New
!13 A2 = LA2 - PHIB2*[(PI2-PHI3*Crd)+EQT/(1-EQT)]

 coef_A(1,4) = -OMGY1
 coef_A(1,5) = -OMGY2
 coef_A(2,6) = -OMGK1
 coef_A(2,7) = -OMGK2
 coef_A(3,1) = -(1.0D0-EQT)*YVC
 coef_A(3,2) =  KVC*DELTA 
 coef_A(4,6) = -ALPHA
 coef_A(4,10)=  ALPHA-1.0D0
 coef_A(5,7) = -ALPHA
 coef_A(5,11)=  ALPHA-1.0D0 
 coef_A(6,4) = -1.0D0
 coef_A(7,5) = -1.0D0
 coef_A(8,4) = -YVC*(1.0D0-EQT)
 coef_A(8,6) =  KVC*DELTA
 coef_A(9,3) = -1.0D0/OMGC2
 coef_A(9,8) =  OMGC1/OMGC2
 coef_A(10,4)= -1.0D0/(1.0D0+PSI2)
 coef_A(10,8)=  PSI1/(1.0D0+PSI2)
 coef_A(10,12)= 2.0*PSI2/THETA/(1.0D0+PSI2)
 coef_A(11,5)= -1.0D0/(1.0D0+PSI2)
 coef_A(11,9)=  PSI1/(1.0D0+PSI2)
 coef_A(11,13)= 2.0*PSI2/THETA/(1.0D0+PSI2)
 
 coef_B(4,1) = 1.0D0 
 coef_B(5,2) = 1.0D0
 coef_B(10,1)= 2.0*PSI2/THETA/(1.0D0+PSI2)
 coef_B(11,2)= 2.0*PSI2/THETA/(1.0D0+PSI2)
 coef_B(12,1)= 1.0D0
 coef_B(13,2)= 1.0D0
!New
 coef_B(13,3)= PHIB2*PHI3

 coef_D(3,1) = -YVC*EQT
 coef_D(6,1) = -KVC*((0.0+DELTA)/ALPHA/(1.0-EQT))
 coef_D(7,1) = -KVC*((0.0+DELTA)/ALPHA/(1.0-EQT))
 coef_D(8,1) = -YVC*EQT
 coef_D(12,1)= -PHIB1*(PI1+EQT)
 coef_D(13,1)= -PHIB2*(PI2+EQT)
 
! Yang------------------ coef_TERM=inv(A)*B*y+inv(A)*D ------------------
 call matinv2(coef_A,13)
 inv_A = coef_A
 coef_TM1 = matmul(inv_A,coef_B)
 coef_TM2 = matmul(inv_A,coef_D)  
!New
! RC1TERM = coef_TM1(8,1)*E(15,KAG-1)+coef_TM1(8,2)*E(16,KAG-1)+coef_TM1(8,3)*E(27,KAG)+coef_TM2(8,1) 
! RC2TERM = coef_TM1(9,1)*E(15,KAG-1)+coef_TM1(9,2)*E(16,KAG-1)+coef_TM1(9,3)*E(27,KAG)+coef_TM2(9,1)	  

 do i=1,2
    do j=1,3
       term_coef(i,j) = coef_TM1(8+i-1,j)
    end do
    term_coef(i,4) = coef_TM2(9+i-1,1)
 end do
    
 open(15,file='term_coef')
 do i=1,2
    write(15,99)(term_coef(i,j),j=1,4)
 end do	
 close(15)
 99 format(F12.7) 

 END


! ============================================================================	  
! Yang-- I add the subroutine here to calculate terminal conditions

SUBROUTINE matinv2(a,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  double precision, DIMENSION(n,n), INTENT(INOUT)  :: a
  double precision, ALLOCATABLE :: y(:,:)
  double precision :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=1.
  ENDDO
  !     LU decompose the matrix just once
  CALL  lu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL lu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE matinv2

!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE lu_decompose(a,n,indx,d)
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  double precision :: sum , tiny, aamax, dum, d
  double precision, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  double precision, ALLOCATABLE :: vv(:)

  tiny=1.0e-20
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=0.
     DO j=1,n
        IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == 0.) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=0.
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (dum >= aamax) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == 0.)  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE lu_decompose

!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion


SUBROUTINE lu_linear_equation(a,n,indx,b)
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  double precision :: sum 
  double precision, DIMENSION(n,n) :: a
  double precision, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE lu_linear_equation
