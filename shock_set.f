      PROGRAM MAIN
      INTEGER N,M
      PARAMETER(NPER=400,NVAR=10)
C.  NPER: Number of shock periods
C.  NVAR: Number of shocks 	  
      REAL Shock(NPER,NVAR), Radm(NPER,5)
      open(unit=11,file='random.data',status='old')
	  do 10 i=1,NPER
          read(11,*) Radm(i,:)
 10   continue		  
      close(11) 	  
	  do 30 i=1,NPER
         do 20 j=1,NVAR
            Shock(i,j) = 0.0
 20      continue
 30   continue
      do 50 i=1,NPER	
C     Shocks to the individual labor*********************           
         Shock(i,8) = Radm(i,3)
         Shock(i,9) = Radm(i,4)				
 50   continue
      open(unit=13,file='bootshocks.txt')
      do 70 i=1,NPER 
         write(13,99) Shock(i,:)
 70   continue
 99   format(10(F13.9,2X))
      close(13)
	  
      end
   
