% clear;
%Input data

load ('data_399')
X = data_399(1:390,1:17);   % Design variable
Y = data_399(1:390,19);     % Objective functions
Z = data_399(1:390,20:23);  % Constraints

IdxFolds = srgtsGetKfolds(Y, 5, 'ClusterBal');
opt = {'BH','MQ' ,'IMQ' ,'TPS','G'};
% opt = {'MQ'};
% Von mises stress

for i = 1
    for j = 1:10 
        options = srgtsRBFSetOptions(X,Y(:,1), @rbf_build, [], opt{i}, j*0.1, 0);
        rbf (i,j) = srgtsRBFFit(options);
        PRESSRMS_rbf(i,j) = srgtsCrossValidation(options);
    end
end

[best_opt,best_ord] = find(PRESSRMS_rbf==min(PRESSRMS_rbf(:)),1); % [2,10]
Sur_F = rbf (best_opt,best_ord);
PRESSRMS_RBF_F = PRESSRMS_rbf(best_opt,best_ord);

%% constraint 1

Y = Z(:,1);

for i = 1:5
    for j = 1:10
        options = srgtsRBFSetOptions(X,Y(:,1), @rbf_build, 1, opt{i}, j*0.1, 1);
        rbf(i,j) = srgtsRBFFit(options);
        PRESSRMS_rbf(i,j) = srgtsCrossValidation(options,10,IdxFolds);
    end
end
[best_opt,best_ord] = find(PRESSRMS_rbf==min(PRESSRMS_rbf(:)),1); %[ 1,1]
Sur_G1 = rbf (best_opt,best_ord);
PRESSRMS_RBF_G1 = PRESSRMS_rbf(best_opt,best_ord);
%% constraint 2

Y = Z(:,2);

for i = 1:5
    for j = 1:10
        options = srgtsRBFSetOptions(X,Y(:,1), @rbf_build, 1, opt{i}, j*0.1, 1);
        rbf(i,j) = srgtsRBFFit(options);
        PRESSRMS_rbf(i,j) = srgtsCrossValidation(options,10,IdxFolds);
    end
end

[best_opt,best_ord] = find(PRESSRMS_rbf==min(PRESSRMS_rbf(:)),1); %[2,10]
Sur_G2 = rbf (best_opt,best_ord); 
PRESSRMS_RBF_G2 = PRESSRMS_rbf(best_opt,best_ord);

%% constraint 3

Y = Z(:,3);

for i = 1:5
    for j = 1:10
        options = srgtsRBFSetOptions(X,Y(:,1), @rbf_build, 1, opt{i}, j*0.1, 1);
        rbf(i,j) = srgtsRBFFit(options);
        PRESSRMS_rbf(i,j) = srgtsCrossValidation(options,10,IdxFolds);
    end
end

[best_opt,best_ord] = find(PRESSRMS_rbf==min(PRESSRMS_rbf(:)),1); %[2,10]

Sur_G3 = rbf (best_opt,best_ord);
PRESSRMS_RBF_G3 = PRESSRMS_rbf(best_opt,best_ord);

%% constraint 4

Y = Z(:,4);

for i = 1:5
    for j = 1:10
        options = srgtsRBFSetOptions(X,Y(:,1), @rbf_build, 1, opt{i}, j*0.1, 1);
        rbf(i,j) = srgtsRBFFit(options);
        PRESSRMS_rbf(i,j) = srgtsCrossValidation(options,10,IdxFolds);
    end
end

[best_opt,best_ord] = find(PRESSRMS_rbf==min(PRESSRMS_rbf(:)),1); %[2,10]
Sur_G4 = rbf (best_opt,best_ord);
PRESSRMS_RBF_G4 = PRESSRMS_rbf(best_opt,best_ord);
% ssTot = sum((Z(:,3)-mean(Z(:,3))).^2);
% ssRes = min(PRESSRMS_RBF_G3(:)).^2;
% r2Pred_F= 1-ssRes./ssTot
save('Design_surrogate_RBF_1','X','Y','Z','Sur_F','Sur_G1','Sur_G2','Sur_G3','Sur_G4',...
    'PRESSRMS_RBF_F','PRESSRMS_RBF_G1','PRESSRMS_RBF_G2','PRESSRMS_RBF_G3','PRESSRMS_RBF_G4')