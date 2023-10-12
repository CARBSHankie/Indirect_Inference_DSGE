%STEP 1: Load all the actual data from Excel================================================================
clear;
[DATA,txt,~]=xlsread('UK1870','identical','A1:AD148');  %size: 148*30 
YEAR = DATA(:,1);
ACT_DATA = DATA(:,2:26);
VARIABLE = txt(1,2:26);
%STEP 2: Use VAR(1) to generate the inial data of expectation variables Expected Values=============
% Here assume both groups are initially identical, estimate EC1 only by VAR(1) and set EC2=EC1.
    CMP=ACT_DATA(:,9);                  % Individual Consumption Reading
    EDG =CMP;                           % EDG: Endogenous variables in the VAR, set to CMP
    FST = 1;                            % FST: Number of the forecasting period by the estimated VAR, set to 1
    N = size(EDG, 1);
    Y = EDG(2:N, :);
    X = EDG(1:N-1, :);
    EXPT = zeros(N+FST-1, size(EDG, 2));
    B = (inv(X' * X)) * X' * Y;
    for i = 1:N
        EXPT(i, :) = EDG(i, :) * B;
    end
    if FST > 1
       for i = N+1:N+FST-1
           EXPT(i, :) = EXPT(i-1, :) * B;
       end
    end
    EC1 = EXPT;
EC(:,1) = EC1;
EC(:,2) = EC1;
ACT_NEW = ACT_DATA;
ACT_NEW(:,24:25) = EC;              % Update the EC values of the act_data
%STEP 3: Write the actual data with expectation data to "act_data.data"=====================================
fid1=fopen('act_data.data','wt'); % The path of writing-to
[m,n]=size(ACT_NEW);
for i=1:m
    for j=1:n
        if j==n
           fprintf(fid1,'%14.9f\n',ACT_NEW(i,j));
        else
           fprintf(fid1,'%14.9f\t',ACT_NEW(i,j));
        end
    end
end
fclose(fid1);
%STEP 4: Save the the starting endo and exo data into plain data files=======================================
START_LAG = 3;
N_EDG = 16;
N_EXG = 32;
START_EDG = zeros(N_EDG,START_LAG);
START_EXG = zeros(N_EXG,START_LAG);
for j=1:START_LAG
    for i=1:16        
        START_EDG(i,j) = ACT_NEW(j,i);
        START_EXG(i,j) = ACT_NEW(j,i);        
    end
    for i=17:19
        START_EXG(i,j) = ACT_NEW(j,i);
    end
    for i=27:32
        START_EXG(i,j) = ACT_NEW(j,i-7);
    end
end
fid2=fopen('start_lags_endog','wt'); % The path of writing-to
fid3=fopen('start_lags_exog','wt');
for i=1:N_EDG
    for j=1:START_LAG
        if j==START_LAG
           fprintf(fid2,'%14.9f\n',START_EDG(i,j));
        else
           fprintf(fid2,'%14.9f\t',START_EDG(i,j));
        end
    end      
end
for i=1:N_EXG
    for j=1:START_LAG
        if j==START_LAG
           fprintf(fid3,'%14.9f\n',START_EXG(i,j));
        else
           fprintf(fid3,'%14.9f\t',START_EXG(i,j));
        end
    end        
end
fclose(fid2);
fclose(fid3);





