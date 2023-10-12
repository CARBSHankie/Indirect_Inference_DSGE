clear;
NSIM = 144;
NPER = 143;
NVAR = 25;
NOUT = 20;
MIU = 0.5;
SIM_IN = zeros(NSIM, NVAR);
SIM_OUT = zeros(NPER, NOUT);
%STEP 1: Read in data from simulation/act_growth data files  
fid = fopen('simulation', 'r');
for i = 1:NSIM
    SIM_IN(i, :) = fscanf(fid, '%f', NVAR);
end
fclose(fid);
g_act = importdata('act_growth.txt');
%STEP 2: generate capital/income share of G1 and logarithmic/ arithmetic growth rates 
for i = 1:NPER
    SIM_OUT(i, 1:16) = SIM_IN(i, 1:16);                          %16 endogenous variables
    SIM_OUT(i, 17) = MIU * exp(SIM_IN(i, 7) - SIM_IN(i, 3));     %Capital share of G1
    SIM_OUT(i, 18) = MIU * exp(SIM_IN(i, 5) - SIM_IN(i, 2));     %Income share of G1
    SIM_OUT(i, 19) = SIM_IN(i + 1, 2) - SIM_IN(i, 2);            %Logarithmic growth rate
    SIM_OUT(i, 20) = exp(SIM_IN(i + 1, 2) - SIM_IN(i, 2)) - 1.0; %Arithmetic growth rate
end
sim   = SIM_OUT;
nper  = size(sim,1);
nvar  = size(sim,2);
%STEP 3: Calculate relevant series to be plotted 
g_compare = zeros(nper,2);
changes = zeros(nper-1,3);
for i=1:nper
  g_compare(i,1) = sim(i,19)-g_act(i,1); %Difference of logarithmic growth and actual growth   
  g_compare(i,2) = sim(i,20)-g_act(i,2); %Difference of arithmetic growth and actual growth
end
for i= 1:nper-1
 changes(i,1) = sim(i+1,17) - sim(i,17);  % change of the capital share of G1
 changes(i,2) = sim(i+1,18) - sim(i,18);  % change of the income share of G1
 changes(i,3) = sim(i+1,19) - sim(i,19);  % change of the logarithmic growth rate
 changes(i,4) = sim(i+1,20) - sim(i,20);  % change of the arithmetic growth rate 
end
%%STEP 4: Plot the relevant figures 
%===========================================================================================
figure(1);                                  % Plot simulated IneqK & Loagrithmic growth rate 
subplot(2,1,1);
    plot(1:nper,sim(1:nper,17),'k','LineWidth',1.5);
    axis_value=axis;
    title('Capital share of G1');
    axis([1 nper axis_value(3) axis_value(4)]);   
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
subplot(2,1,2);
    plot(1:nper,sim(1:nper,19),'k','LineWidth',1.5);
    axis_value=axis;
    title('Economic growth rate');
    axis([1 nper axis_value(3) axis_value(4)]);
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0.25 0.25 10 10]);
print(gcf,'-dpng','-r1080','tendency1.jpg');
%===========================================================================================
figure(2);                                   % Plot simulated IneqK & Arithmetic growth rate
subplot(2,1,1);
    plot(1:nper,sim(1:nper,17),'k','LineWidth',1.5);
    axis_value=axis;
    title('Capital share of G1');
    axis([1 nper axis_value(3) axis_value(4)]);   
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
subplot(2,1,2);
    plot(1:nper,sim(1:nper,20),'k','LineWidth',1.5);
    axis_value=axis;
    title('Economic growth rate');
    axis([1 nper axis_value(3) axis_value(4)]);
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0.25 0.25 10 10]);
print(gcf,'-dpng','-r1080','tendency2.jpg');
%============================================================================================
figure(3);   % Plot simulated IneqK & Logarithmic growth rate deviated from the actual growth
subplot(2,1,1);
    plot(1:nper,sim(1:nper,17),'k','LineWidth',1.5);
    axis_value=axis;
    title('Capital share of G1');
    axis([1 nper axis_value(3) axis_value(4)]); 
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
subplot(2,1,2);
    plot(1:nper,g_compare(1:nper,1),'k','LineWidth',1.5);
    axis_value=axis;
    title('Difference of simulated & actual growth rate');
    axis([1 nper axis_value(3) axis_value(4)]); 
    set(gca,'yticklabel',get(gca,'ytick'));   
    set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0.25 0.25 8 8]);
print(gcf,'-dpng','-r1080','tendency3.jpg');
%=============================================================================================
figure(4);     % Plot simulated IneqK & Arithmetic growth rate deviated from the actual growth 
subplot(2,1,1);
    plot(1:nper,sim(1:nper,17),'k','LineWidth',1.5);
    axis_value=axis;
    title('Capital share of G1');
    axis([1 nper axis_value(3) axis_value(4)]);   
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
subplot(2,1,2);
    plot(1:nper,g_compare(1:nper,2),'k','LineWidth',1.5);
    axis_value=axis;
    title('Difference of simulated & actual growth rate');
    axis([1 nper axis_value(3) axis_value(4)]);
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0.25 0.25 10 10]);
print(gcf,'-dpng','-r1080','tendency4.jpg');
%=============================================================================================
figure(5);            % Plot change in capital share of G1 & change of logarithmic growth rate 
subplot(2,1,1);
    plot(1:nper-1,changes(1:nper-1,1),'k','LineWidth',1.5);
    axis_value=axis;
    title('Change in capital share of G1');
    axis([1 nper axis_value(3) axis_value(4)]);   
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
subplot(2,1,2);
    plot(1:nper-1,changes(1:nper-1,3),'k','LineWidth',1.5);
    axis_value=axis;
    title('Change in aggregate growth rate');
    axis([1 nper axis_value(3) axis_value(4)]);
    set(gca,'yticklabel',get(gca,'ytick'));     
    set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0.25 0.25 10 10]);
print(gcf,'-dpng','-r1080','tendency5.jpg');
