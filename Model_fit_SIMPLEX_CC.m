clear all
close all

coculture_data = load('co_culture_data_fit.txt')
coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
stdev_PhTri = coculture_data_stdev.Stdev(6:10)


coculture_data(:,1)

vTime=coculture_data(1:end,1);
vTimeNew=vTime(1):2.5:vTime(end);
PhTAC125 = coculture_data(:,2);
PhTri = coculture_data(:,3);

vTimeNew=vTime(1):3:vTime(end);
%vPhTAC125inter = spline(vTime,PhTAC125,vTimeNew);
vPhTAC125inter = interp1(vTime,PhTAC125,vTimeNew, 'cubic');
% vPhTAC125inter(1) = PhTAC125(1);
% vPhTAC125inter(3) = PhTAC125(2);
% vPhTAC125inter(5) = PhTAC125(3);
% vPhTAC125inter(7) = PhTAC125(4);
% vPhTAC125inter(9) = PhTAC125(5);

vPhTriinter = interp1(vTime,PhTri,vTimeNew, 'cubic');
Data=[vTimeNew' vPhTAC125inter' vPhTriinter'];

%% plotting real data
figure(1)
plot(Data(:,1), Data(:,2), 'Marker', 'o', 'Color', '[1, 0, 0]', 'LineWidth', 1.5)
hold on
plot(Data(:,1), Data(:,3), 'Marker', 'o', 'Color', '[0.4940, 0.1840, 0.5560]','LineWidth', 1.5)
xlabel('Time (days)', 'FontSize', 18);
ylabel('Log of cell counts', 'FontSize', 18)
coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
stdev_PhTri = coculture_data_stdev.Stdev(6:10)
legend( 'PhTAC125', 'PhTri')

figure(2)
errorbar(coculture_data(:,1) , coculture_data(:,2), stdev_PhTAC125, '-o', 'Color', '[1, 0, 0]')
hold on
errorbar(coculture_data(:,1) , coculture_data(:,3), stdev_PhTri,'-o', 'Color', '[0.4940, 0.1840, 0.5560]')
title('Co-culture original data', 'FontSize', 22)
xlabel('Time (days)', 'FontSize', 18);
ylabel('Log of cell counts', 'FontSize', 18)
coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
stdev_PhTri = coculture_data_stdev.Stdev(6:10)
legend( 'PhTAC125', 'PhTri')




%% plotting biomass ratio

% figure(4)
% plot(Data(:,1),Data(:,3)./Data(:,2), 'Marker', 'o')
% title('Co-culture biomass ratio from interpolated data', 'FontSize', 22)
% xlabel('Time (days)', 'FontSize', 18);
% ylabel('Log of cell counts ratio', 'FontSize', 18)


%% output interpolated data to file

fileID = fopen('interpolated_data.txt','w');
nbytes = fprintf(fileID,'%5d %5d %5d\n' ,Data')
fclose(fileID);
type('interpolated_data.txt')


v_mu_D = .19;
CC_D = .18;
v_delta_D = .02;
v_mu_B = .32;
CC_B = .01;
K_doc1 = .24;
K_doc2 = .24;
v_delta_B = .22;
lambda = .1;
delta_doc1 = .000001;
delta_doc2 = .000001;
rAB = .01;


ModelParameters = [v_mu_D, CC_D, v_delta_D, v_mu_B, CC_B, K_doc1,v_delta_B, lambda, K_doc2, delta_doc1, delta_doc2, rAB];


%ModelParameters = load('ParBoni.dat');
Data_PC=Data;
parCal = importdata('SetPar.dat')
ModelParameters = parCal;
ModelParameters(9) = K_doc2;
ModelParameters(10) = 0.0000001;
ModelParameters(11) = 0.0000001;
ModelParameters(12) = .001;

[nPoints,nData]=size(Data_PC);

options=foptions;
options(1)=1;
options(2)=.1;
options(3)=.1;
options(14)=10000;


if(0)
    [parCal,~] = SIMPLEXL('objectiveFunction',ModelParameters,options,[],Data_PC);
else
    [parCal,FVAL,EXITFLAG,OUTPUT] = SIMPSA('objectiveFunction',ModelParameters',...
    zeros(size(ModelParameters)),20*ones(size(ModelParameters)),[],Data_PC);
end


y0 = [PhTAC125(1) PhTri(1) 1 1];
[time,sol] = ode45(@(t,y) odeSystem(t,y,parCal), vTimeNew, y0);



% figure(5);
% subplot(1,2,1)
% plot(time,sol(:,2), 'Color', '[0.4940, 0.1840, 0.5560]')
% hold on
% e=errorbar(coculture_data(:,1) , coculture_data(:,3), stdev_PhTri, 'o', 'Color', '[0.4940, 0.1840, 0.5560]')
% e.Marker = 'o';
% e.MarkerSize = 5;
% e.CapSize = 5;
% xlabel('Time (days)', 'FontSize', 18);
% ylabel('Log of cell counts', 'FontSize', 18)
% coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
% stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
% stdev_PhTri = coculture_data_stdev.Stdev(6:10)
% legend( 'Model', 'Experimental', 'Location', 'southeast')
% subplot(1,2,2)
% plot(time,sol(:,1),'Color', '[1, 0, 0]')
% hold on
% e=errorbar(coculture_data(:,1) , coculture_data(:,2), stdev_PhTri, 'o', 'Color', '[1, 0, 0]')
% e.Marker = 'o';
% e.MarkerSize = 5;
% e.CapSize = 5;
% xlabel('Time (days)', 'FontSize', 18);
% ylabel('Log of cell counts', 'FontSize', 18)
% coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
% stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
% stdev_PhTri = coculture_data_stdev.Stdev(6:10)
% legend( 'Model', 'Experimental', 'Location', 'southeast')

stdev_PhTAC125_inter = [stdev_PhTAC125(1), 0, stdev_PhTAC125(2), 0, stdev_PhTAC125(3), 0, stdev_PhTAC125(4), 0, stdev_PhTAC125(5), 0]
stdev_PhTri_inter= [stdev_PhTri(1), 0, stdev_PhTri(2), 0, stdev_PhTri(3), 0, stdev_PhTri(4), 0, stdev_PhTri(5), 0]
stdev_cell_count_ratio= [0,0,0,0,0,0,0,0,0,0]
% vPhTAC125inter(1) = PhTAC125(1);
% vPhTAC125inter(3) = PhTAC125(2);
% vPhTAC125inter(5) = PhTAC125(3);
% vPhTAC125inter(7) = PhTAC125(4);
% vPhTAC125inter(9) = PhTAC125(5);



figure
subplot(3,1,1)
plot(time,sol(:,1),'Color', '[1, 0, 0]', 'LineWidth', 1.5)
hold on
e=errorbar(vTimeNew , vPhTAC125inter, stdev_PhTAC125_inter, 'o', 'Color', '[1, 0, 0]', 'LineWidth', 1.5)
e.Marker = 'o';
e.MarkerSize = 5;
e.CapSize = 5;
ylabel('Log of cell counts', 'FontSize', 12)
legend( 'Model', 'Experimental (interpolated)', 'Location', 'southeast')
title('A')

subplot(3,1,2)
plot(time,sol(:,2), 'Color', '[0.4660, 0.6740, 0.1880]', 'LineWidth', 1.5)
hold on
e=errorbar(vTimeNew , vPhTriinter, stdev_PhTri_inter, 'o', 'Color', '[0.4660, 0.6740, 0.1880]', 'LineWidth', 1.5)
e.Marker = 'o';
e.MarkerSize = 5;
e.CapSize = 5;
ylabel('Log of cell counts', 'FontSize', 12)
coculture_data_stdev = readtable('co_culture_aggregate_log.txt')
stdev_PhTAC125 = coculture_data_stdev.Stdev(1:5)
stdev_PhTri = coculture_data_stdev.Stdev(6:10)
legend( 'Model', 'Experimental (interpolated)', 'Location', 'southeast')
title('B')

subplot(3,1,3)
plot(time,sol(:,2)./sol(:,1), 'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 1.5)
hold on
errorbar(Data(:,1),Data(:,3)./Data(:,2), stdev_cell_count_ratio, 'o', 'Color', '[0.9290 0.6940 0.1250]', 'LineWidth', 1.5)
xlabel('Time (days)', 'FontSize', 18);
ylabel('Cell counts ratio', 'FontSize', 12)
legend( 'Model', 'Experimental (interpolated)', 'Location', 'southeast')
title('C')

save('parDef2.dat', 'parCal');

ModelParametersNames = ["v_mu_D", "CC_D", "v_delta_D","v_mu_B", "CC_B" ,"K_doc1","v_delta_B" ,"lambda", "K_doc2","delta_doc1", "delta_doc2","rAB"]



table(ModelParameters' , ModelParametersNames')

