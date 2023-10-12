
 program main
 implicit none
 call calc_growth
 end
!============================================================================= 
 subroutine calc_growth
 implicit none
 integer i,j
 integer,parameter:: nper=147,nact=25,next=380,nlag=3,nout=244
!! nper: periods of actual data
!! nact: number of variables in the file "act_data.data"
!! next: number of periods to extend to
!! nlag: number of lags when stimulate data
!! nout: number of periods of growth to print out
 real(8),dimension(nper,nact):: act_data
 real(8),dimension(next,1):: Y_extent    !! next: Number of periods to extend
 real(8),dimension(2,1):: beta_trend
 real(8),dimension(nper,2):: X_trend
 real(8),dimension(nper,1):: Y_trend, trend_resid
 real(8),dimension(next-1,2):: g_act
 open(11,file='act_data.data',status='old')
 do i=1,nper
    read(11,*)act_data(i,:)
 end do
 close(11)
!! Order of Y in "act_sim_data" is 2  
 do i=1,nper
    Y_extent(i,1) = act_data(i,2)
 end do	
!! Extend data to 380 periods from 147 periods
 do i=1,nper
    X_trend(i,1)=1
    X_trend(i,2)=i
 end do
 Y_trend(:,1)=act_data(:,2)
 call OLS(Y_trend,X_trend,nper,2,beta_trend,trend_resid)
 do j=nper+1,next	
    Y_extent(j,1)=beta_trend(1,1)+beta_trend(2,1)*j
 end do
!! Calculate the net growth rate of the extended Y by two ways
!! column 1 is by difference in lnY; column 2 is by Y'/Y-1
 do i=1,next-1
    g_act(i,1)=Y_extent(i+1,1)-Y_extent(i,1)
    g_act(i,2)=exp(Y_extent(i+1,1)-Y_extent(i,1))-1.0	
 end do
 open(15,file='act_growth.txt')
 do i=nlag+1,nout+nlag
    write(15,'(2(f12.9,5X))') g_act(i,:)
 end do
 close(15)
 end subroutine calc_growth
!=========================================================================================  
 SUBROUTINE OLS(Y,X,nper,nexog,beta,residual)
 IMPLICIT NONE
! Estimate OLS regression Y on X
! INPUTS: 
! Y=endogenous variables
! X=exogenous variables
! nper=number of periods
! nexog=number of variables in X
 integer nper, nexog, i
 real(8), dimension(nper,nexog) :: X
 real(8), dimension(nper,1) :: Y, Xbeta, residual
 real(8), dimension (nexog,nexog) ::  XTX, XTXI, XTY
 real(8), dimension (nexog,nper) :: XT
 real(8), dimension (nexog,1) :: beta
 integer, dimension(nexog) :: INDX
  
 XT=transpose(X)
 XTX=matmul(XT,X)
! call inverse(XTX,XTXI,nexog)
! call matinv(nexog,XTX,XTXI)
! call MIGS(XTX,nexog,XTXI,INDX)
 call matinv2(XTX,nexog)
 XTXI=XTX

 XTY=matmul(XT,Y)
 beta=matmul(XTXI,XTY)
   
 Xbeta=matmul(X,beta)
 
 do i=1,nper
   residual(i,1)=Y(i,1)-Xbeta(i,1)
 end do

 END SUBROUTINE OLS
 
!=========================================================================================
 SUBROUTINE VARE_NOCONST(series,nper,nvar,beta)
 IMPLICIT NONE
! Estimate AR parameter for data:
! INPUTS: 
! series=data
! nper=number of periods
! nvar=number of variables
 integer nper, nvar, i, j
 real(8),dimension(nper,nvar) :: series
 real(8),dimension(nvar,nvar) ::  beta, XTX, XTXI, XTY
 real(8),dimension(nper-1,nvar) :: Y, X
 real(8),dimension(nvar,nper-1) :: XT
 integer,dimension(nvar) :: INDX
 
 do i=1,nper-1
   do j=1,nvar
     Y(i,j)=series(i+1,j)
     X(i,j)=series(i,j)
   end do
 end do
 XT=transpose(X)

 XTX=matmul(XT,X)

! call inverse(XTX,XTXI,nvar)
! call matinv(nvar,XTX,XTXI)
! call MIGS(XTX,nvar,XTXI,INDX)
 call matinv2(XTX,nvar)
 XTXI=XTX

 XTY=matmul(XT,Y)
 beta=matmul(XTXI,XTY)

 END SUBROUTINE VARE_NOCONST
 
!=========================================================================================
 SUBROUTINE VARE(series,nper,nvar,beta)
 IMPLICIT NONE
! Estimate AR parameter for data:
! INPUTS: 
! series=data
! nper=number of periods
! nvar=number of variables
 integer nper, nvar, i, j
 real(8), dimension (nper,nvar) :: series
 real(8), dimension (nvar,nvar) ::  beta, XTX, XTXI, XTY
 real(8), dimension (nper-1,nvar) :: Y, X
 real(8), dimension (nvar,nper-1) :: XT
 integer, dimension(nvar) :: INDX
 
 do i=1,nper-1
   do j=1,nvar
     Y(i,j)=series(i+1,j)
     X(i,j)=series(i,j)
   end do
   X(i,nvar+1)=1.0
 end do
 XT=transpose(X)

 XTX=matmul(XT,X)

! call inverse(XTX,XTXI,nvar)
! call matinv(nvar,XTX,XTXI)
! call MIGS(XTX,nvar,XTXI,INDX)
 call matinv2(XTX,nvar)
 XTXI=XTX

 XTY=matmul(XT,Y)
 beta=matmul(XTXI,XTY)

 END SUBROUTINE VARE
 
!=========================================================================================
  subroutine inverse(aa,cc,nn)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! aa(nn,nn) - array of coefficients for matrix A
! nn      - dimension
! output ...
! cc(nn,nn) - inverse matrix of A
! comments ...
! the original matrix aa(nn,nn) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer nn
real(8) aa(nn,nn), cc(nn,nn), aatmp(nn,nn)
real(8) LL(nn,nn), UU(nn,nn), bb(nn), dd(nn), xx(nn)
real(8) coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
LL=0.0
UU=0.0
bb=0.0

! duplicate a matrix so it won't be destroyed
aatmp=aa

! step 1: forward elimination
do k=1, nn-1
   do i=k+1,nn
      coeff=aatmp(i,k)/aatmp(k,k)
      LL(i,k) = coeff
      do j=k+1,nn
         aatmp(i,j) = aatmp(i,j)-coeff*aatmp(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,nn
  LL(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,nn
  do i=1,j
    UU(i,j) = aatmp(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix cc
do k=1,nn
  bb(k)=1.0
  dd(1) = bb(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,nn
    dd(i)=bb(i)
    do j=1,i-1
      dd(i) = dd(i) - LL(i,j)*dd(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  xx(nn)=dd(nn)/UU(nn,nn)
  do i = nn-1,1,-1
    xx(i) = dd(i)
    do j=nn,i+1,-1
      xx(i)=xx(i)-UU(i,j)*xx(j)
    end do
    xx(i) = xx(i)/UU(i,i)
  end do
! Step 3c: fill the solutions x(nn) into column k of cc
  do i=1,nn
    cc(i,k) = xx(i)
  end do
  bb(k)=0.0
end do
end subroutine inverse

! Matrix inversion: B = Inv(A) by Gauss-Jordan method
! A and B are n by n matrices
Subroutine Matinv(n,A,B)  
  ! Labels: 10, 20, 30
  parameter(MMAX=25,NMAX=10)
  integer n
  real(8)  A(MMAX,MMAX), B(MMAX,2*MMAX)  
  integer i,j,k
  real(8) bb

  do i = 1, n
    do j = 1, n
      B(i,j + n) = 0.d0
      B(i,j) = A(i,j)
    end do
    B(i,i + n) = 1.d0
  end do 

  do k = 1, n
    if (k.eq.n) goto 10
    m = k
    do i = k+1, n
      if (abs(B(i,k)) > abs(B(m,k)))  m = i
    end do
    if (m == k) goto 10
    do j = k, 2*n
      bb = B(k,j)
      B(k,j) = B(m,j)
      B(m,j) = bb
    end do
10  do j = k+1, 2*n 
      B(k,j) = B(k,j) / B(k,k)
    end do
    if (k.eq.1) goto 20
    do i = 1, k-1
      do j = k+1, 2*n
        B(i,j) = B(i,j) - B(i,k) * B(k,j)
      end do
    end do
    if (k.eq.n) goto 30
20  do i = k+1, n
      do j = k+1, 2*n
        B(i,j) = B(i,j) - B(i,k) * B(k,j)
      end do
    end do
  end do    ! k loop

30  do i = 1, n
      do j = 1, n   
        B(i,j) = B(i,j + n)
      end do
    end do
  return
end ! Matinv()


 SUBROUTINE SORT(ARR,n)
!SORT ARRAY INTO NUMERICAL ORDER
 IMPLICIT NONE
 INTEGER i,j,k,n
 real(8), dimension(n,1), intent(INOUT) :: arr 
! real(8), dimension(:), intent(INOUT) :: arr 
 real(8) a
 
! n=size(arr,1)
 do j=2,n
   a=arr(j,1)
   do i=j-1,1,-1
     if (arr(i,1) <= a ) exit
     arr(i+1,1)=arr(i,1)
   end do
   arr(i+1,1)=a
 end do

 END SUBROUTINE SORT


SUBROUTINE matinv2(a,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  real(8), DIMENSION(n,n), INTENT(INOUT)  :: a
  real(8), ALLOCATABLE :: y(:,:)
  real(8) :: d
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
  real(8) :: sum , tiny, aamax, dum, d
  real(8), DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  real(8), ALLOCATABLE :: vv(:)

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
  real(8) :: sum 
  real(8), DIMENSION(n,n) :: a
  real(8), DIMENSION(n) :: b
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
 
 
 
