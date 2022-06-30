% -----------------------------------------------------------------------------
% Reference
% Metamodel-based lightweight design of B-pillar with TWB structure via support vector regression
% Feng Pan, Ping Zhu*, Yu Zhang The
% -------------------------------------------------------------------------------
clear;
load('HICmodel.mat','HIC_rbf');
load('Data.mat');
% Sample size
NN = [25,50,75,100];

% Details of RVs
lb = [1 1 1]; ub = [2.5 2.5 2.5];
std_P = 0.1;
SB = makedist('Normal',1.75,std_P*1.75);
TB = truncate(SB,1,2.5);

std_P = 0.05;
% Prior probability
Alpha = 7/15; % 15 - Total samples, 7 - failure samples
SO1 = makedist('Normal',1.615,std_P*1.615);
SO2 = makedist('Normal',1.822,std_P*1.822);
SO3 = makedist('Normal',1.276,std_P*1.276);

TO1 = truncate(SO1,1,2.5);TO2 = truncate(SO2,1,2.5);TO3 = truncate(SO3,1,2.5);

DD = Data(:,1:3); % Samples 2, 6, 8, 9, 10, 14, 15

Fail =15; %[ 9, 10, 14, 15]; %2, 6, 8,
for p = 1:length(Fail)
    
    Sam = Fail(p);
    % Deterministic optima
    TD =  DD(Sam,:);          % Deterministic - [1.615,1.822,1.276];
    % TD = [2.00310000000000,2.06470000000000,1.00260000000000];
    % TD = [1.93040000000000,1.33040000000000,1.21770000000000];
    
    
    [Pos_cmom_act, Pos_lmom_act] = actprob(TB,TD,HIC_rbf);
    
    
    for i = 1:length(NN)
        N = NN(i);
        for j = 1:100
            
            T = random(TB,N,3);
            [W_Nf,W_f,W_f_ext] = sep_data(T,HIC_rbf,TB);
            
            % Moments estimation for pearson system
            [M,S,skew,kurt] = moment_pears(W_Nf);
            [M_f,S_f,skew_f,kurt_f] = moment_pears(W_f);
            [M_f_ext,S_f_ext,skew_f_ext,kurt_f_ext] = moment_pears(W_f_ext);
            
            % Parameter estimation for Lmoment
            [P,Dist,L_sample,D,D_min] = parameter_identify(W_Nf,1);
            [Pf,Distf,L_sample1,D1,D_min1] = parameter_identify(W_f,1);
            [Pf_ext,Distf_ext,L_sample1_ext,D1_ext,D_min1_ext] = parameter_identify(W_f_ext,1);
            
            %TO = [random(TO1,1e4,1),random(TO2,1e4,1),random(TO3,1e4,1)];
            W = weight(TD)*ones(1e4,1);
            Alpha = normrnd(7/15,0.4667*0.025,1e4,1);
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
    
    
    save(['Side_impact_Bayesian_inference_Deter',num2str(Sam),'.mat'])
    
end

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

% Weight of B-piller
function W = weight(t)
W = 0.518.*t(:,1)+0.8733.*t(:,2)+1.0413.*t(:,3);
end

function [W_Nf,W_f,W_f_ext] = sep_data(T,HIC_rbf,TB)
HIC = srgtsRBFEvaluate(T,HIC_rbf);
[W_Nf,W_f] = sep_data1(HIC,330);
W_f_ext = [W_f;max(weight(random(TB,1e6,3)))];

    function [W_Nf,W_f] = sep_data1(HIC,threshold)
        L = threshold;
        idx_nf = HIC<L;
        W = weight(T);
        W_Nf = W(idx_nf);
        W_f = W(~idx_nf);
    end
end


function [Pos_cmom_act, Pos_lmom_act] = actprob(TB,TD,HIC_rbf)

T = random(TB,1e6,3);
[W_Nf,W_f,~] = sep_data(T,HIC_rbf,TB);

% Moments estimation for pearson system
[M,S,skew,kurt] = moment_pears(W_Nf);
[M_f,S_f,skew_f,kurt_f] = moment_pears(W_f);

% Parameter estimation for Lmoment
[P,Dist] = parameter_identify(W_Nf,1);
[Pf,Distf] = parameter_identify(W_f,1);

%TO = [random(TO1,1e4,1),random(TO2,1e4,1),random(TO3,1e4,1)];
W = weight(TD);
Alpha = 7/15; % normrnd(7/15,0.4667*0.05,1e4,1);
% Posterior probability
Pos_cmom_act = basyesian_cmom(W,Alpha,M,S,skew,kurt,M_f,S_f,skew_f,kurt_f);
Pos_lmom_act = basyesian_lmom(W,Alpha,P,Dist,Pf,Distf);
Pos_cmom_act(isnan(Pos_cmom_act))=[];Pos_lmom_act(isnan(Pos_lmom_act))=[];
end