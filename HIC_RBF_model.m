clear;

load('Data.mat');
% Total 15 data points are available
x = Data(:,1:3); % Input data t1,t2,t3
y = Data(:,4); % Head injury criterion

% KRG model fo HIC

[HIC_rbf,PRESSRMS_HIC] = rbf_forall(x,y);

SST = sum((y-mean(y)).^2);
SSR = (PRESSRMS_HIC^2)*15;
pred_R2 = 1-(SSR/SST);

M = mean(y);
save('HICmodel.mat','HIC_rbf','PRESSRMS_HIC','pred_R2','M')



% options = srgtsSVRSetOptions(x, y,'LinearSpline',[],inf,'einsensitive',0);
% [HIC_rbf,PRESSRMS_HIC] = srgtsSVRFit(options);
% SST = sum((y-mean(y)).^2);
% SSR = (PRESSRMS_HIC^2)*15;
% pred_R2 = 1-(SSR/SST)
% P=PRESSRMS_HIC/mean(y)
% 
% Mdl = fitrsvm(x,y,'KernelFunction','gaussian','KernelScale','auto',...
%     'Standardize',true);
% yfit = predict(Mdl,x);
% SST = sum((y-mean(y)).^2);
% SSR = sum((yfit-y).^2);
% pred_R2 = 1-(SSR/SST)

