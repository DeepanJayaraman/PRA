% Data for failure and non failure event
clear;
load('all_realizations_data.mat')

F = Failure_incidents;
Nf = Normal_incidents;

% Parent distributoin identification
% 
% % Using C-moment
% [M,S,skew,kurt] = moment_pears(Nf);
% [M_f,S_f,skew_f,kurt_f] = moment_pears(F);
% 
% P_Nf = pearsrnd(X,M,S,skew,kurt);


% Using L-moments
[P,Dist,L_sample,D,D_min] = parameter_identify(Nf(1:125),1);
[P1,Dist1,L_sample1,D1,D_min1] = parameter_identify(F(1:8),1);


X_Nf = Random_l(char(Dist),P.P,1e4,1);
X_f = Random_l(char(Dist1),P1.P,1e4,1);

% Generate random numbers

