#Step 1: Generate new data files for identical groups===========
module load matlab/R2019a
matlab < load_data.m
#Step 2: Format the data files for Fortran======================
gfortran toolbox90.f90 read_bench.f90 -o read_bench.exe 
./read_bench.exe
cp  bench_data1a.data bench_data1.data
cp  bench_data2a.data bench_data2.data
cat bench_end.data >> bench_data1.data
cat bench_exo.data >> bench_data1.data 
cat bench_datab.data >> bench_data1.data  
cat bench_end.data >> bench_data2.data
cat bench_exo.data >> bench_data2.data 
cat bench_datab.data >> bench_data2.data
#Step 3: Compile all relevant programmes=========================
gfortran calc_terminal.f90 -o calc_terminal.exe
gfortran calc_shocks.f90 -o calc_shocks.exe
gfortran bench_mod1.f -o bench_mod1.exe
#gfortran toolbox77.f SIM1.f90 -o bench_mod1.exe
gfortran bench_mod2.f -o bench_mod2.exe
gfortran shock_set.f -o shock_set.exe
gfortran calc_act_gy.f90 -o calc_act_gy.exe
c++ normal.cpp -o normal.exe
#Step 4: Run the main programme==================================
./calc_shocks.exe
./calc_terminal.exe
cp nerror qerror
cp start_lags_endog lags_endog
cp start_lags_exog lags_exog
echo 3 > nlag
for i in {4..147}        #Controls the periods of simulation
   do                    #number of files in tmp
     rm allout.out	
     ./bench_mod1.exe < bench_data1.data > output.out
     grep ITER output.out		
   done
cp allout.out base
./normal.exe
./shock_set.exe
cp bootshocks.txt qerror

echo 3 > nlag
for i in {4..147}
    do
      rm -f allout.out
      ./bench_mod2.exe < bench_data2.data > output.out
	  grep ITER output.out
    done
cp allout.out simulation
#Step 5: Plot the tendency graphs=================================
./calc_act_gy.exe
matlab < tendency_plot.m
#Step 6: Clear the intermediate files for next simulation========
rm -rf calc_terminal.exe calc_shocks.exe bench_mod1.exe bench_mod2.exe shock_set.exe calc_act_gy.exe normal.exe
rm act_growth.txt resids.data bgp_const.data ar_coeffs.data shocks.data term_coef  base random.data bootshocks.txt simulation
rm -f fort.* nlag lags_endog lags_exog start_lags_endog start_lags_exog allout.out all_model gather_sim.data qerror
rm read_bench.exe bench_end.data bench_exo.data
rm bench_data1.data bench_data2.data act_data.data
rm $PWD/tmp/* output.out
