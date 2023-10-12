// This code aims to randomly draw a number from a uniform distribution then transfromed to normal distribution
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <cmath>
using namespace std;
double u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,x[5];
double mean=0.01;             ////the value of the mean will significantly inpact the results of simulations
int ramseed;
// declare the variable "top" which in an integer
int main (int argc, char**argv ) {
ifstream infile("randseed.data",ios::in);	
ofstream outfile("random.data",ios::out);
infile >> ramseed;
srand(ramseed);
// initialise the random
    for(int m=0; m<400; m++) {	
//  random draw m times with the for loop
   u1 = rand()/(double)(RAND_MAX);
   u2 = rand()/(double)(RAND_MAX);
   u3 = rand()/(double)(RAND_MAX);   
   u4 = rand()/(double)(RAND_MAX);
   u5 = rand()/(double)(RAND_MAX);
   u6 = rand()/(double)(RAND_MAX);
   u7 = rand()/(double)(RAND_MAX);
   u8 = rand()/(double)(RAND_MAX);
   u9 = rand()/(double)(RAND_MAX);
   u10= rand()/(double)(RAND_MAX);
//Box-Muller transformation, an old method to transform the uniformly distributed numbers to normally distributed ones    
   x[1] = mean*sqrt(-2*log(u1))*cos(2*3.1415926*u3)-mean;
   x[2] = mean*sqrt(-2*log(u2))*cos(2*3.1415926*u4)+mean;   
   x[3] = mean*sqrt(-2*log(u5))*cos(2*3.1415926*u7)+mean;
   x[4] = mean*sqrt(-2*log(u6))*cos(2*3.1415926*u8)-mean;   
   x[5] = mean*sqrt(-2*log(u9))*cos(2*3.1415926*u10);  
   outfile<<setprecision(2)<<x[1]<<" "<<setprecision(2)<<x[2]<<" "<<setprecision(2)<<x[3]<<" "<<setprecision(2)<<x[4]<<" "<<setprecision(2)<<x[5];  
       outfile<<endl;	
       }	
infile.close();   
outfile.close();
return 0;
}
