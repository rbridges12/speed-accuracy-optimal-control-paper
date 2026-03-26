clc; close all, clear;

%For k_u = 0.05, k_t = 1
EE_init_x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
EE_init_y = [.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3];
EE_target_x = [0,0,0,0,0,0,0,0,.1,.1,.1,.1,.1,.1,.1,.1];
EE_target_y = [.35,.36,.38,.40,.45,.5,.55,.6,.35,.36,.38,.40,.45,.5, .55, .6];
target_radius = [.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03,.03];
target_radius2 = [.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04,.04];
time03 = [0.141458976785484, 0.154322467081960, 0.176855420661658, 0.196491924214380,.238186716519532, 0.274711512733866, 0.309957382135256, 0.349918702136348...
    0.196550630652105,0.193061077938777, 0.191265661483701, 0.195628614006720, 0.224730785551278, 0.260138481070536, 0.297582789245952, 0.342523782667534];
time04 = [0.141420480572573, 0.154315429856226, 0.176816241902346, 0.196500173515422, 0.238272643202502 ,0.274521902632999, 0.310203947683087, 0.349815032394272...
    0.196558330865217, 0.193070159749635, 0.191252175476443, 0.195665413251651,0.224773492051688, 0.260152322793365, 0.297651554092767, 0.342421478211354];





x_diff = EE_target_x - EE_init_x;
y_diff = EE_target_y - EE_init_y;
D = sqrt(x_diff.^2 + y_diff.^2);
W1 = target_radius * 2;
W2 = target_radius2 * 2;
ID1 = log2(2.*D./W1);
ID2 = log2(2.*D./W2);

figure()
plot(ID1, time03, 'o', LineWidth=1.5)
hold on
%plot(ID1, time04, 'o')


a_paper = 0.0234;
a_lower = 0.0047;
a_upper = 0.5239;
b_paper = 0.0550;
b_lower = 0.0393;
b_upper = 0.1987;


ID_range = linspace(0,4,20);
fitts_time_paper = a_paper + b_paper.*ID_range;


% plot(ID_range, fitts_time_paper)

const = polyfit(ID1, time03, 1);
a = const(2)
b = const(1)

fitts_time = a + b.*ID_range;
plot(ID_range, fitts_time, LineWidth=1.5);

xlabel('Index of Difficulty')
ylabel('Movement Duration (s)')

legend('Model', "Fitt's Law", Location="northwest")
title("")
%% 


mdl = fitlm(ID1, time03);
mdl.Rsquared.Ordinary

