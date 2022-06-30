clear;
load('Design_surrogate_RBF.mat','Sur_F','Sur_G2');
load('data_399.mat')
% Sample size
NN = [25,50,75,100];

% Design Variables
X1 = data_399(:,1:7);

lb = [121, 21,  21, 21, 31, 151, 21];        % lower bound specification
ub = [147.5, 77.5, 78, 78, 88, 205, 48];     % upper bound specification

% Noise Parameters (Material)
X2 = data_399(:,8:17);

% Bounds
LB = [lb min(X2)];
UB = [ub max(X2)];

% Constraints
C = data_399(:,20:23);

% Failure criterian: Maximum total deformation = 35
FC = C(:,2);

% Response - Von mises stress
Von = data_399(:,19);

std_P = 0.05;
% Prior probability
Alpha = sum(FC>35)/399; % 399 - Total samples, 101 - failure samples

% Failure data
TD = [X1(147,:) X2(147,:)];

[Pos_cmom_act, Pos_lmom_act] = actprob(TD,LB,UB,Sur_F,Sur_G2,FC);

for i = 1:length(NN)
    N = NN(i);
    for j = 1:100
        
        T = rand_T(LB,UB,N);
        [V_Nf,V_f,V_f_ext] = sep_data(T,LB,UB,Sur_F,Sur_G2);
        
        % Moments estimation for pearson system
        [M,S,skew,kurt] = moment_pears(V_Nf);
        [M_f,S_f,skew_f,kurt_f] = moment_pears(V_f);
        [M_f_ext,S_f_ext,skew_f_ext,kurt_f_ext] = moment_pears(V_f_ext);
        
        % Parameter estimation for Lmoment
        [P,Dist,L_sample,D,D_min] = parameter_identify(V_Nf,1);
        [Pf,Distf,L_sample1,D1,D_min1] = parameter_identify(V_f,1);
        [Pf_ext,Distf_ext,L_sample1_ext,D1_ext,D_min1_ext] = parameter_identify(V_f_ext,1);
        
        TDD = rand_n(TD,0.025,1e4);
        W = Vonm(TDD,Sur_F);
        Alpha = normrnd(sum(FC>35)/399,(sum(FC>35)/399)*0.05,1e4,1);
        % Posterior probability
        Pos_cmom = basyesian_cmom(W,Alpha,M,S,skew,kurt,M_f,S_f,skew_f,kurt_f);
        Pos_lmom = basyesian_lmom(W,Alpha,P,Dist,Pf,Distf);
        Pos_cmom_ext = basyesian_cmom(W,Alpha,M,S,skew,kurt,M_f_ext,S_f_ext,skew_f_ext,kurt_f_ext);
        Pos_lmom_ext = basyesian_lmom(W,Alpha,P,Dist,Pf_ext,Distf_ext);
        Pos_cmom(isnan(Pos_cmom))=[];Pos_lmom(isnan(Pos_lmom))=[];
        Pos_cmom_ext(isnan(Pos_cmom_ext))=[];Pos_lmom_ext(isnan(Pos_lmom_ext))=[];
        
        switch N
            case 25
                M_pos_cmom_25(j) = mean(Pos_cmom);
                M_pos_lmom_25(j) = mean(Pos_lmom);
                M_pos_cmom_ext_25(j) = mean(Pos_cmom_ext);
                M_pos_lmom_ext_25(j) = mean(Pos_lmom_ext);
                Pos_cmom_25 = Pos_cmom;
                Pos_lmom_25 = Pos_lmom;
                Pos_cmom_ext_25 = Pos_cmom_ext;
                Pos_lmom_ext_25 = Pos_lmom_ext;
        
            case 50
                M_pos_cmom_50(j) = mean(Pos_cmom);
                M_pos_lmom_50(j) = mean(Pos_lmom);
                M_pos_cmom_ext_50(j) = mean(Pos_cmom_ext);
                M_pos_lmom_ext_50(j) = mean(Pos_lmom_ext);
                Pos_cmom_50 = Pos_cmom;
                Pos_lmom_50 = Pos_lmom;
                Pos_cmom_ext_50 = Pos_cmom_ext;
                Pos_lmom_ext_50 = Pos_lmom_ext;
            case 75
                M_pos_cmom_75(j) = mean(Pos_cmom);
                M_pos_lmom_75(j) = mean(Pos_lmom);
                M_pos_cmom_ext_75(j) = mean(Pos_cmom_ext);
                M_pos_lmom_ext_75(j) = mean(Pos_lmom_ext);
                Pos_cmom_75 = Pos_cmom;
                Pos_lmom_75 = Pos_lmom;
                Pos_cmom_ext_75 = Pos_cmom_ext;
                Pos_lmom_ext_75 = Pos_lmom_ext;
            case 100
                M_pos_cmom_100(j) = mean(Pos_cmom);
                M_pos_lmom_100(j) = mean(Pos_lmom);
                M_pos_cmom_ext_100(j) = mean(Pos_cmom_ext);
                M_pos_lmom_ext_100(j) = mean(Pos_lmom_ext);
                Pos_cmom_100 = Pos_cmom;
                Pos_lmom_100 = Pos_lmom;
                Pos_cmom_ext_100 = Pos_cmom_ext;
                Pos_lmom_ext_100 = Pos_lmom_ext;
        end
        
    end
    
end


save('Rotordisk_Bayesian_inference_17_sam147.mat')



% L-moment approach
function Posterior = basyesian_lmom(X,Alpha,P,Dist,P_f,Dist_f)

T_f = f_fail(X,P_f,Dist_f).*Alpha;
T = f(X,P,Dist).*(1-Alpha);
D = T_f+T;
P = T_f./D;

for i = 1:length(X)
    if and(T_f(i)==0,D(i)==0)==1
        P(i)= 0;
        Posterior = P;
    else
        Posterior = P;
    end
end
    function PDF = f(X,P,Dist)
        PDF = PDF_l(X,char(Dist),P.P);
    end
    function PDF1 = f_fail(X,P_f,Dist_f)
        PDF1 = PDF_l(X,char(Dist_f),P_f.P);
    end
end

% C-moment approach
function Posterior = basyesian_cmom(X,Alpha,M,S,skew,kurt,M_f,S_f,skew_f,kurt_f)

T_f = f_fail(X,M_f,S_f,skew_f,kurt_f).*Alpha;
T = f(X,M,S,skew,kurt).*(1-Alpha);
D  = T_f+T;
P = T_f./D;

for i = 1:length(X)
    if and(T_f(i)==0,D(i)==0)==1
        P(i)= 0;
        Posterior = P;
    else
        Posterior = P;
    end
end
    function PDF = f(X,M,S,skew,kurt)
        PDF = pearspdf(X,M,S,skew,kurt);
    end
    function PDF = f_fail(X,M_f,S_f,skew_f,kurt_f)
        PDF = pearspdf(X,M_f,S_f,skew_f,kurt_f);
    end
end

% Vonmises stress
function Von = Vonm(t,Sur_F)
Von = srgtsRBFEvaluate(t,Sur_F);
end
% Data separation
function [V_Nf,V_f,V_f_ext] = sep_data(T,LB,UB,Sur_F,Sur_G2)

Dec = srgtsRBFEvaluate(T,Sur_G2); % Decision constraint
[V_Nf,V_f] = sep_data1(Sur_F,Dec,35);
V_f_ext = [V_f;max(Vonm(rand_T(LB,UB,1e6),Sur_F))];

    function [V_Nf,V_f] = sep_data1(Sur_F,Dec,threshold)
        L = threshold;
        idx_nf = Dec<L;
        W = Vonm(T,Sur_F);
        V_Nf = W(idx_nf);
        V_f = W(~idx_nf);
    end
end

function random_T = rand_T(lb,ub,N)
random_T = uni(lb,ub,N);
    function RR = uni(A,B,N)
        LA = length(A);
        for ii = 1:LA
            RR(:,ii) = unifrnd(A(ii),B(ii),N,1);
        end
    end
end

function RRR = rand_n(T,std_P,N)
for ij = 1: length(T)
    RRR(:,ij) = normrnd(T(ij),std_P*T(ij),N,1);
end
end


function [Pos_cmom_act, Pos_lmom_act] = actprob(TD,LB,UB,Sur_F,Sur_G2,FC)
        T = rand_T(LB,UB,1e6);
        [V_Nf,V_f,~] = sep_data(T,LB,UB,Sur_F,Sur_G2);
        
        % Moments estimation for pearson system
        [M,S,skew,kurt] = moment_pears(V_Nf);
        [M_f,S_f,skew_f,kurt_f] = moment_pears(V_f);
        
        % Parameter estimation for Lmoment
        [P,Dist,L_sample,D,D_min] = parameter_identify(V_Nf,1);
        [Pf,Distf,L_sample1,D1,D_min1] = parameter_identify(V_f,1);
        
        TDD = TD; %rand_n(TD,0.05,1e4);
        W = Vonm(TDD,Sur_F);
        Alpha = sum(FC>35)/399; %normrnd(sum(FC>35)/399,(sum(FC>35)/399)*0.05,1e4,1);
        % Posterior probability
        Pos_cmom_act = basyesian_cmom(W,Alpha,M,S,skew,kurt,M_f,S_f,skew_f,kurt_f);
        Pos_lmom_act = basyesian_lmom(W,Alpha,P,Dist,Pf,Distf);
        Pos_cmom_act(isnan(Pos_cmom_act))=[];Pos_lmom_act(isnan(Pos_lmom_act))=[];
end