clear;

delx = 0.001;
x = 25:delx:84;
% PDF
f1 = wblpdf3(x,2.14,23.25,84); % Failure
f0 = logwblpdf3(x,10.56,3.11,90); % non Failure

%% Actual population generation
NP = 1e6; % Population size

% CDF
cdf = cumsum(f0).*delx;
cdf1  = cumsum(f1).*delx;

% Random number - Non failure events
C = cdf;
[CC, index] = unique(C);
X_act = sort(interp1(CC, x(index), rand(NP,1)));
X_act(isnan(X_act))=[];


% Random number - Failure events
C1 = cdf1;
[CC1, index1] = unique(C1);
X1_act = sort(interp1(CC1, x(index1), rand(NP,1)));
X1_act(isnan(X1_act))=[];

% Posterior prob. from ref
        T = 57;
        E = normrnd(T,T*0.01,1e4,1);
        
        Nof  = round(normrnd(9,0.9,1e4,1));
        Alpha = Nof./138;%*ones(1e4,1);

Posterior_act = PosteriorProb(f0,f1,Alpha);
NofFailure_act = Posterior_act*6;


%% Sample
NN = [5,7,9]; % Sample size
I = 100; % No.of iterations

%Preallocation
% L = length(E);
% NofFailure_cmom = ones(L,I); NofFailure_lmom = ones(L,I);
% NofFailure_cmom_ext = ones(L,I); NofFailure_lmom_ext = ones(L,I);

for i = 1:length(NN)
    N = NN(i);
    for j = 1:I
        % Sample generation
        % Random number - Non failure events
        X_nonfail = sort(interp1(CC, x(index), rand(N*5,1)));
        X_nonfail(isnan(X_nonfail))=[];
        
        % Random number - Failure events
        X_fail = sort(interp1(CC1, x(index1), rand(N,1)));
        X_fail(isnan(X_fail))=[];
        Ext_f = min(interp1(CC1, x(index1), rand(1e5,1)));
        RM = randperm(N,1);
        Xf = X_fail;
        Xf(RM) = Ext_f;
        X_fail_ext = Xf;
        
        % Moments for pearson
        X_f = X_fail; X_f(isnan(X_f))=[];
        X_nf = X_nonfail;
        X_nf(isnan(X_nf))=[];
        
        [M,S,skew,kurt] = moment_pears(X_nf);
        [M_f,S_f,skew_f,kurt_f] = moment_pears(X_f);
        
        % With extreme
        X_f = X_fail_ext; X_f(isnan(X_f))=[];
        
        [M_f_ext,S_f_ext,skew_f_ext,kurt_f_ext] = moment_pears(X_f);
        
        % Dist identification using Lmom
        [P,Dist,L_sample,D,D_min] = parameter_identify_noln(X_nonfail,1);
        [P1,Dist1,L_sample1,D1,D_min1] = parameter_identify_noln(X_fail,1);
        [P1_ext,Dist1_ext,L_sample1_ext,D1_ext,D_min1_ext] = parameter_identify_noln(X_fail_ext,1);
        
        %% Posterioer probability
        
        Posterior_lmom = basyesian_lmom(E,Alpha,P,char(Dist),...
            P1,Dist1);
        NofFailure_lmom(:,j) = mode(round(Posterior_lmom*6));
        
        Posterior_lmom_ext = basyesian_lmom(E,Alpha,P,char(Dist),...
            P1_ext,Dist1_ext);
        NofFailure_lmom_ext(:,j) = mode(round(Posterior_lmom_ext*6));
        
        Posterior_cmom = basyesian_cmom(E,Alpha,M,S,skew,...
            kurt,M_f,S_f,skew_f,kurt_f);
        NofFailure_cmom(:,j) = mode(round(Posterior_cmom*6));
        
        Posterior_cmom_ext = basyesian_cmom(E,Alpha,M,S,skew,...
            kurt,M_f_ext,S_f_ext,skew_f_ext,kurt_f_ext);
        NofFailure_cmom_ext(:,j) = mode(round(Posterior_cmom_ext*6));
        
        
    end
    save(['31temp_failures_',num2str(N),'_samples'],'Posterior_act','Posterior_cmom','Posterior_lmom','NofFailure_cmom',...
        'NofFailure_lmom','NofFailure_act',...
        'Posterior_cmom_ext','Posterior_lmom_ext','NofFailure_cmom_ext',...
        'NofFailure_lmom_ext','E','Alpha');
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
        %         PDF = pearspdfforall(X,M_f,S_f,skew_f,kurt_f,1);
    end
end
function [P,Distribution_type,L_sample,D,D_min,Parameter1] = parameter_identify_noln(X,K)
X(isnan(X))=[];

[Distribution_type,L_sample,D,D_min] = Identify_dist_noln(X,K);
L1= L_sample(1);L2 = L_sample(2);T3 = L_sample(3);T4 = L_sample(4);
for k = 1:K
    Parameter1 = Parameter_estimation(X,Distribution_type(k),L1,L2,T3,T4);
    C = max(size(Parameter1));
    switch C
        case 3
        P(1,k)=struct('P',Parameter1);
        case 2
        P(1,k)=struct('P',[Parameter1,0]);
        case 1
        P(1,k)=struct('P',[Parameter1,0,0]);
    end
end
end