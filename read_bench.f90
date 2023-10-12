program main
 use toolbox90
 implicit none
 integer i,j,k,m
 integer,parameter:: nper=147,nact=25,nend=16,nexo=32,next=380,ncol=4, nrow=next/ncol 
 real(8),dimension(nper,nact):: act_data
 real(8),dimension(nend,nrow,ncol):: endo     
 real(8),dimension(nexo,nrow,ncol):: exog
 real(8),dimension(next,nact):: forecast         !!next: Number of total periods in the model data file
 real(8),dimension(2,1):: beta_trend
 real(8),dimension(nper,2):: X_trend
 real(8),dimension(nper,1):: Y_trend, trend_resid
 real(8),dimension(2,2):: var_beta
 integer,dimension(15) :: non_stationary_indices
 real :: random_fluctuation, stddev
 !STEP 1: Load the act_data.data 
  open(11,file='act_data.data',status='old')
  do i=1,nper
    read(11,*)act_data(i,:)
  end do
  close(11) 
 !STEP 2: Extend the data to forecasting periods 
  do i=1,nper
    forecast(i,:) = act_data(i,:)           !!First, load the first nper periods actual data
  end do	
 stddev = 0.05 
 do i=nper+1,next                          !!Second, extend data to next periods from the (nper+1) period
     call random_seed()
     call random_number(random_fluctuation)
    random_fluctuation = (random_fluctuation - 0.5) * 2.0 * stddev 
   forecast(i,:) = act_data(nper,:)+ random_fluctuation         !!stationary data equals the value at the last actual periods + a random component
  end do	
  do i=1,nper
    X_trend(i,1)=1        !1s
    X_trend(i,2)=i        !Time trend 
  end do
  non_stationary_indices = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 18, 20, 24, 25]
  do i=1, size(non_stationary_indices)
    Y_trend(:, 1) = act_data(:, non_stationary_indices(i))
    call OLS(Y_trend, X_trend, nper, 2, beta_trend, trend_resid) !!non-stationary data equals the value from a time trend regression
    do j=nper+1, next
      forecast(j, non_stationary_indices(i)) = beta_trend(1, 1) + beta_trend(2, 1) * j
    end do
  end do
 !STEP 3: Write in the bench data files 
  do j=1,nrow
    do k=1,ncol
       do i=1,nend
          endo(i,j,k)=0.0D0
       end do
       do m=1,nexo		  
          exog(m,j,k)=0.0D0
       end do			 
    end do
  end do
  do j=1,nrow
    do k=1,ncol
       do i=1,nend                                 
          endo(i,j,k) = forecast((j-1)*ncol+k,i)
       end do                                      
       do i=1,nend+1     
          exog(i,j,k) = forecast((j-1)*ncol+k,i)   !the 17th series is the TAU(also 17th in forecast())
       end do
          exog(27,j,k) = forecast((j-1)*ncol+k,20) !the 27th series is the Credit(20th in forecast())
       do i=31,32
          exog(i,j,k) = forecast((j-1)*ncol+k,i-7) !the 31th/32th is the EC1/EC2(24th 25th forecast())
       end do                                      !All others are keeping as zeros 
    end do                                         !Data was written like this because of the rigid set-up of the old Fortran Programme
  end do 
  open(15,file='bench_end.data')
  do i=1,nend
    do j=1,nrow 
       write(15,'(X,4f14.9,X)')(endo(i,j,k),k=1,ncol)  !data was written by extended series:Total 1520 rows=380 periods *16 endos / 4 columns
	end do   
  end do
  close(15)
  open(16,file='bench_exo.data')                      ! bench_exo.data include the bench_end.data in the first 1520 rows
  do i=1,nexo
    do j=1,nrow 
       write(16,'(X,4f14.9,X)')(exog(i,j,k),k=1,ncol) !data was written by extended series:Total 3040 rows=380 periods * 32 exogs / 4 columns
	end do                                         
  end do
  close(16) 
 end program
