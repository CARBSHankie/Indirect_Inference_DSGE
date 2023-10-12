!!=========================================================================
!! This new version has one difference that the credit now has time trend, 
!! because the credit data now are replaced by "M4 stocdk / GDP"
!!=========================================================================

 program main
 implicit none
 
 call CALC_RESIDS
 call CALC_SHOCKS
  
 end

 subroutine CALC_RESIDS
 implicit none
 integer i,j,k
!***********************************************************************
 double precision ALPHA,BETA,DELTA,MIU1,MIU2,OMGY1,OMGY2
 double precision OMGK1,OMGK2,OMGC1,OMGC2,PSI1,PSI2
 double precision YVC,KVC,LABOR1,LABOR2,PI1,PI2,EQR,EQT
 double precision THETA,RHO1,RHO2,PHIB1,PHIB2,ETA1,ETA2
 double precision YVW1,YVW2
!***********************************************************************
 integer, parameter:: ncoeffs=30,nper=147,nact=25,nsh=10 
!!Yang-- nact: Number of variables in the "act_data.data" 
!!Yang-- nsh:  Number of shocks I will back out
 real(8) coef(ncoeffs)
 real(8),dimension(nper,nact):: act_data,act_sim_data
 real(8),dimension(nper,nsh):: f,lr
!!Yang-- real(8),dimension(nper-1,nsh):: resids !!Yang-- Original
 real(8),dimension(nper-2,nsh):: resids 
 
 open(10,file='coef.data',status='old')
 read(10,*)coef
 close(10)
 
 ALPHA = coef(1)     !!Yang-- Share of capital in production α
 BETA  = coef(2)     !!Yang-- Utility Discount Factor β
 DELTA = coef(3)     !!Yang-- Capital Discount Factor δ  
 MIU1  = coef(4)     !!Yang-- Group1's Population Share
 MIU2  = 1.0D0-MIU1  !!Yang-- Group2's Population Share 
 OMGY1 = coef(5)     !!Yang-- Top 10% Income Share ωy,1
 OMGY2 = coef(6)     !!Yang-- Top 90% Income Share ωy,2 
 OMGK1 = coef(7)     !!Yang-- Top 10% Capital Share ωk,1
 OMGK2 = coef(8)     !!Yang-- Top 90% Capital Share ωk,2 
 OMGC1 = coef(9)     !!Yang-- Top 10% Consumption Share ωc,1 
 OMGC2 = coef(10)    !!Yang-- Top 90% Consumption Share ωc,2
 PSI1  = coef(11)    !!Yang-- The power of consumption in the utility
 PSI2  = coef(12)    !!Yang-- The power of leisure in the utility
 YVC   = coef(13)    !!Yang-- Steady State Y/C Ratio
 KVC   = coef(14)    !!Yang-- Steady State K/C Ratio
 LABOR1= coef(15)    !!Yang-- Steady State Labour of the rich
 LABOR2= coef(16)    !!Yang-- Steady State Labour of the poor
 PI1   = coef(17)    !!Yang-- Steady State Penalty Cost of the rich, π'1
 PI2   = coef(18)    !!Yang-- Steady State Penalty Cost of the poor, π'2
 EQR   = coef(19)    !!Yang-- Steady State Real Interest Rate 
 EQT   = coef(20)    !!Yang-- Steady State marginal income tax rate  
 THETA = coef(21)    !!Yang-- θ2 in the equation of Z
 RHO1  = coef(22)    !!Yang-- ρ1 in the equation: π'i(t)=ρ1*π'i(t-1)-ρ2*[μi/wyi*(ki/K)^2]
 RHO2  = coef(23)    !!Yang-- ρ2 in the equation: π'i(t)=ρ1*π'i(t-1)-ρ2*[μi/wyi*(ki/K)^2]
 PHIB1 = coef(24)    !!Yang-- φb1 in the equation: dlnA1(t+1)=...-φb1*π'1(t)+...
 PHIB2 = coef(25)    !!Yang-- φb2 in the equation: dlnA2(t+1)=...-φb2*π'2(t)+... 

 ETA1  = (1.0D0-EQT)*ALPHA*(YVC/KVC)/(DELTA+EQR) !!Yang-- Coefficient of lnYt in the equation of lnkt	  
 ETA2  = 1.0D0/(DELTA+EQR)-ETA1/PSI1             !!Yang-- Minus Coefficient of Rt in the equation of lnkt
	  
! YVW1  = LABOR1/(1.0D0-ALPHA)	  
! PHIB1 =(1.0D0-EQT)*BETA*THETA*YVW1/(1.0D0-BETA)/((1.0D0-EQT+PI1)**2.0D0)
! YVW2  = LABOR2/(1.0D0-ALPHA)	  
! PHIB2 =(1.0D0-EQT)*BETA*THETA*YVW2/(1.0D0-BETA)/((1.0D0-EQT+PI2)**2.0D0)
 	  
!=====================================================================================
 
 open(11,file='act_data.data',status='old')
 do i=1,nper
    read(11,*)act_data(i,:)
 end do
 close(11)

!!Yang--"act_sim_data" is introduced in order to create expectations  
 do i=1,nper
    do j=1,nact
       act_sim_data(i,j) = act_data(i,j)
    end do
 end do

!======================================================================================   
!!Yang-- The following aims to calculate the value of RHS (or some transformed RHS)
!======================================================================================   
! Calculate RHS
!--------- act_sim_data order:
!  1  2  3  4  5   6   7   8   9   10  11  12  13  14  15  16  
!  R, Y, K, C, Y1, Y2, K1, K2, C1, C2, N1, N2, P1, P2, A1, A2,
!  17   18  19  20      21    22   23   24   25
!  TAU, A,  P', Credit, None, EK1, EK2, EC1, EC2 
!
!--------- shocks I need to back out:
!  1   2    3     4     5      6      7      8      9    10
!  A,  P', SH_Y, SH_K, SH_C, SH_C1, SH_C2, SH_N1, SH_N2, Credit
 do j=1,nper
    if (j==1) then
       do k=1,nsh
          f(j,k)=0.0
       end do
    else
!!Yang--Aggregate Productivity lnAt
   f(j,1)=0.0D0
!!Yang--p'
   f(j,2)=0.0D0
!!Yang-- Y=Weighted Sum of yi
   f(j,3)=OMGY1*act_sim_data(j,5)+OMGY2*act_sim_data(j,6)
!!Yang-- K=Weighted Sum of ki
   f(j,4)=OMGK1*act_sim_data(j,7)+OMGK2*act_sim_data(j,8)
!!Yang-- C (Market Clearing)
!   f(j,5)=(1.0D0-act_sim_data(j,17))*YVC*act_sim_data(j,2)-KVC*(act_sim_data(j,3)-(1.0D0-DELTA)*act_sim_data(j-1,3))
!   f(j,5)=YVC*((1.0D0-EQT)*act_sim_data(j,2)-act_sim_data(j,17))-KVC*(act_sim_data(j,3)-(1.0D0-DELTA)*act_sim_data(j-1,3))
   f(j,5)=YVC*((1.0D0-EQT)*act_sim_data(j,2)-act_sim_data(j,17))-KVC*(act_sim_data(j,3)-(1.0D0-DELTA)*act_sim_data(j-1,3))
!!Yang-- C1 Individual 1's Euler
   f(j,6)=act_sim_data(j,24)-act_sim_data(j,1)/PSI1
!!Yang-- Consumption Aggregation Equation
   f(j,7)=(act_sim_data(j,4)-OMGC1*act_sim_data(j,9))/OMGC2
!!Yang-- N1
   f(j,8)=(act_sim_data(j,5)-PSI1*act_sim_data(j,9)- 2.0*PSI2/THETA*(act_sim_data(j,15)-act_sim_data(j-1,15)))/(1.0D0+PSI2)
!!Yang-- N2
   f(j,9)=(act_sim_data(j,6)-PSI1*act_sim_data(j,10)-2.0*PSI2/THETA*(act_sim_data(j,16)-act_sim_data(j-1,16)))/(1.0D0+PSI2) 
!!Yang-- Credit
   f(j,10)=0.0D0 
    end if
 end do
!------------------------------------------------------------------------------------------------
! Calculate LHS-RHS
 do j=1,nper
    if (j==1) then
       do k=1,nsh
          lr(j,k)=0.0
       end do
    else
!!Yang-- lnA	 
       lr(j,1) = act_sim_data(j,18)-f(j,1)
!!Yang-- P'	 
       lr(j,2) = act_sim_data(j,19)-f(j,2)
!!Yang-- Y
       lr(j,3) = act_sim_data(j,2) -f(j,3)
!!Yang-- K	 
       lr(j,4) = act_sim_data(j,3) -f(j,4)
!!Yang-- C	 
       lr(j,5) = act_sim_data(j,4) -f(j,5)
!!Yang-- C1	 
       lr(j,6) = act_sim_data(j,9) -f(j,6)
!!Yang-- C2	 
       lr(j,7) = act_sim_data(j,10)-f(j,7)
!!Yang-- N1		 
       lr(j,8) = act_sim_data(j,11)-f(j,8)
!!Yang-- N2	 
       lr(j,9) = act_sim_data(j,12)-f(j,9)
!!Yang-- Credit	 
       lr(j,10)= act_sim_data(j,20)-f(j,10)	 
    end if
 end do

 do i=1,nper-1
    resids(i,1) = lr(i+1,1)
    resids(i,2) = lr(i+1,2)
    resids(i,3) = lr(i+1,3)
    resids(i,4) = lr(i+1,4)   
    resids(i,5) = lr(i+1,5)
    resids(i,6) = lr(i+1,6)
    resids(i,7) = lr(i+1,7)
    resids(i,8) = lr(i+1,8)
    resids(i,9) = lr(i+1,9)
    resids(i,10)= lr(i+1,10)
 end do

!!Yang-- I discard the 1st and last period of resids
 open(unit=12,file='resids.data')
!!Yang-- do i=1,nper-1 !!Yang-- Original
 do i=2,nper-1 
    write(12,222) resids(i,:)
 end do
 close(12)
 
!!Yang!!
! open(unit=13,file='resids_147.data')
! do i=1,nper
!   write(13,222) lr(i,:)
! end do
! close(13) 
  
 222 format(10f20.8)

 END SUBROUTINE CALC_RESIDS

!==============================================================================================
 SUBROUTINE CALC_SHOCKS
 implicit none
 integer i,j,k
!!Yang!!---Note: Here "nper" below = "nper" in the subroutine "CALC_RESIDS"-1
!!Yang!!---& "nvar" is the number of shocks 
! integer,parameter :: nper=146, nvar=10
 integer,parameter :: nper=145, nvar=10
 real(8),dimension(nper,2):: X_trend
 real(8),dimension(nper,1):: Y_trend, trend_resid
 real(8),dimension(nper-1,1):: Y_ar, X_ar, ar_resid, d_prod,Y_ar_crd
 real(8),dimension(nper-1,2):: X_ar_crd
 real(8),dimension(nper-1,2):: XPP_ar              !!Yang-- This is only for Pi'
 real(8),dimension(nper-2,2):: d_prod_X, err1_X
 real(8),dimension(nper-2,1):: d_prod_Y, prod_shock,credit_shock, err1_Y, shock1
 real(8),dimension(2,1):: beta_trend, beta_ar_PP,beta_ar_crd  !!Yang-- beta_ar_PP is only for Pi'
 real(8),dimension(1,1):: beta_ar, bgp_const
 real(8),dimension(nvar,1):: all_beta_ar 
 real(8),dimension(nper,nvar):: resids, resids_dt
 real(8),dimension(nper-1,nvar):: shocks
 
! resids from period 2 to nper
 open(10,file='resids.data',status='old')
 do i=1,nper
    read(10,*)resids(i,:)
 end do
 close(10)
 101 format(10f20.8)

 do i=1,nper
    X_trend(i,1)=1
    X_trend(i,2)=i
 end do
 
 do i=1,nvar
! New 
  if(i==3 .or. i==4 .or. i==7) then
    shocks(:,i)=0.0
    all_beta_ar(i,1)=0.0
  else	
 
! detrend residuals
    Y_trend(:,1)=resids(:,i)
    call OLS(Y_trend,X_trend,nper,2,beta_trend,trend_resid)
    resids_dt(:,i)=trend_resid(:,1)
!!Yang-- As the new added exo-variable "credit" has no trend at all, here I do not detrend it
    resids_dt(:,10)=resids(:,10)   

!!Yang-- The 1st residual is the lnAt & take difference firstly and then regress on a constant & deA t-1      
    if(i==1) then
       d_prod(:,1)=resids(2:nper,i)-resids(1:nper-1,i)
       do j=1,nper-1
          d_prod_X(j,1)=1
          d_prod_X(j,2)=d_prod(j,1)
          d_prod_Y(j,1)=d_prod(j+1,1)
       end do
       call OLS(d_prod_Y,d_prod_X,nper-2,2,beta_trend,prod_shock)
       bgp_const(1,1)=beta_trend(1,1)
       open(11,file='bgp_const.data')
           write(11,91)bgp_const
 91    format(f15.9)	 
       close(11)
!!Yang-- "beta_trend" for each variable is a 2 by 1 vector. But the 1st value is a constant 
!!Yang-- & only the 2nd is the coefficient of a time trend	 
       all_beta_ar(i,1)=beta_trend(2,1)
       shocks(1,i)=0.0
       shocks(2:nper-1,i)=prod_shock(:,1)
	 
    else if(i==2) then                 !!Yang-- Do not detrend Pi'! with constant
       Y_ar(:,1)=resids(2:nper,i)
	   do j=1,nper-1
          XPP_ar(j,1)=1.0D0
       end do  
       XPP_ar(:,2)=resids(1:nper-1,i)	 
       call OLS(Y_ar,XPP_ar,nper-1,2,beta_ar_PP,ar_resid)
       shocks(:,i)=ar_resid(:,1)        !!Yang-- calculate shocks
       if(beta_ar_PP(2,1)>=1.0) then
          beta_ar_PP(2,1)=0.99999
          do j=2,nper-1
            shocks(j,i)=Y_ar(j,1)-beta_ar_PP(2,1)*XPP_ar(j,2)
          end do
       end if
       all_beta_ar(i,1)=beta_ar_PP(2,1)	 
!    else if(i==2) then                 !!Yang-- Do not detrend Pi'!  without constant
!      Y_ar(:,1)=resids(2:nper,i)
!      X_ar(:,1)=resids(1:nper-1,i)
!      call OLS(Y_ar,X_ar,nper-1,1,beta_ar,ar_resid)
!      shocks(:,i)=ar_resid(:,1)        !!Yang-- calculate shocks
!      if(beta_ar(1,1)>=1.0) then
!         beta_ar(1,1)=0.99999
!         do j=2,nper-1
!            shocks(j,i)=Y_ar(j,1)-beta_ar(1,1)*X_ar(j,1)
!         end do
!      end if
!      all_beta_ar(i,1)=beta_ar(1,1)   

!    else if(i==10) then
!!Yang-- Suppose "credit" following I(2) as A, use "else if(i==10) ..."
!       d_prod(:,1)=resids(2:nper,i)-resids(1:nper-1,i)
!       do j=1,nper-1
!          d_prod_X(j,1)=1
!          d_prod_X(j,2)=d_prod(j,1)
!          d_prod_Y(j,1)=d_prod(j+1,1)
!       end do
!       call OLS(d_prod_Y,d_prod_X,nper-2,2,beta_trend,credit_shock)
!       all_beta_ar(i,1)=beta_trend(2,1)
!       shocks(1,i)=0.0
!       shocks(2:nper-1,i)=credit_shock(:,1)
!!      ----------------------------------------------------------	   
    else
!!Yang-- calculate AR coefficient   
       Y_ar(:,1)=resids_dt(2:nper,i)
       X_ar(:,1)=resids_dt(1:nper-1,i)
       call OLS(Y_ar,X_ar,nper-1,1,beta_ar,ar_resid)
!!Yang-- calculate shocks
       shocks(:,i)=ar_resid(:,1)
       if(beta_ar(1,1)>=1.0) then
          beta_ar(1,1)=0.99999
          do j=2,nper-1
            shocks(j,i)=Y_ar(j,1)-beta_ar(1,1)*X_ar(j,1)
          end do
       end if
       all_beta_ar(i,1)=beta_ar(1,1)
    end if
! New	
  end if
  
 end do
! New
 all_beta_ar(7,1)=all_beta_ar(6,1) 
 
 open(11,file='ar_coeffs.data')
 do i=1,nvar
    write(11,102)all_beta_ar(i,1)
 end do
 close(11)
 102 format(1f15.9)

!!Yang-- I discard the 1st shock with high magnitude.
 open(12,file='shocks.data')
 do i=2,nper-1 
    write(12,103)shocks(i,:)
 end do
 close(12)
 
!!Yang----------------------------- 
! open(13,file='resids_detrend.data')
! do i=2,nper-1 
!   write(13,103)resids_dt(i,:)
! end do
! close(13) 
!!Yang----------------------------- 

 103 format(10f16.10)

 END SUBROUTINE CALC_SHOCKS
 
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
