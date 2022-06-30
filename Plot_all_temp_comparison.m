clear;
addpath('D:\Risk analysis\Challenger example\Different sample size\')
NN = [5,10,15,20];
%T = 42;
for i = 1:length(NN)
    N = NN(i);
    load(['All_temperature_failures_',num2str(N),'_samples'],'NofFailure_cmom',...
        'NofFailure_lmom','NofFailure_act',...
        'NofFailure_cmom_ext',...
        'NofFailure_lmom_ext','E','Alpha');
    delx = 0.001;
    x = 25:delx:84;
    
    NofFailure_cmom(NofFailure_cmom> 6) = 6;
    NofFailure_cmom(NofFailure_cmom< 0) = 0;
    NofFailure_lmom(NofFailure_lmom> 6) = 6;
    NofFailure_lmom(NofFailure_lmom< 0) = 0;
    NofFailure_cmom_ext(NofFailure_cmom_ext> 6) = 6;
    NofFailure_cmom_ext(NofFailure_cmom_ext< 0) = 0;
    NofFailure_lmom_ext(NofFailure_lmom_ext> 6) = 6;
    NofFailure_lmom_ext(NofFailure_lmom_ext< 0) = 0;
    
    % NofFailure_cmom=round(NofFailure_cmom);
    % NofFailure_cmom=round(NofFailure_cmom);
    % NofFailure_lmom=round(NofFailure_lmom);
    % NofFailure_lmom=round(NofFailure_lmom);
    %
    % NofFailure_cmom_ext=round(NofFailure_cmom_ext);
    % NofFailure_cmom_ext=round(NofFailure_cmom_ext);
    % NofFailure_lmom_ext=round(NofFailure_lmom_ext);
    % NofFailure_lmom_ext=round(NofFailure_lmom_ext);
    
    %% [Min_cmom, Max_cmom, Med_cmom] = minmaxmedian(NofFailure_cmom');
    % [Min_lmom, Max_lmom, Med_lmom] = minmaxmedian(NofFailure_lmom');
    % [Min_cmom_ext, Max_cmom_ext, Med_cmom_ext] = minmaxmedian(NofFailure_cmom_ext');
    % [Min_lmom_ext, Max_lmom_ext, Med_lmom_ext] = minmaxmedian(NofFailure_lmom_ext');
    %
    % Legend = {'Population','Minimum','Median','Maximum'};
    %
    % Area_shade_plot(x,NofFailure_act,E,Min_cmom',Med_cmom',Max_cmom','C-moment',Legend);
    % saveas(gcf,['Area_cmom_minmax_',num2str(N),'_samples.fig'])
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
    % print('-dpng', ['Area_cmom_minmax_',num2str(N),'_samples.png'], '-r100')
    %
    % Area_shade_plot(x,NofFailure_act,E,Min_lmom',Med_lmom',Max_lmom','L-moment',Legend);
    % saveas(gcf,['Area_lmom_minmax_',num2str(N),'_samples.fig'])
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
    % print('-dpng', ['Area_lmom_minmax_',num2str(N),'_samples.png'], '-r100')
    %
    % Area_shade_plot(x,NofFailure_act,E,Min_cmom_ext',Med_cmom_ext',Max_cmom_ext','C-moment',Legend);
    % saveas(gcf,['Area_cmom_minmax_',num2str(N),'_samples_ext.fig'])
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
    % print('-dpng', ['Area_cmom_minmax_',num2str(N),'_samples_ext.png'], '-r100')
    %
    % Area_shade_plot(x,NofFailure_act,E,Min_lmom_ext',Med_lmom_ext',Max_lmom_ext','L-moment',Legend);
    % saveas(gcf,['Area_lmom_minmax_',num2str(N),'_samples_ext.fig'])
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
    % print('-dpng',['Area_lmom_minmax_',num2str(N),'_samples_ext.png'], '-r100')
    
    
    %%
    PP = [2575];
    for j = 1:length(PP)
        P = PP(j);
        switch P
            case 1090
                
                P = [10,50,90];
                Cmom_ptile = prctile(NofFailure_cmom',P)';
                Lmom_ptile = prctile(NofFailure_lmom',P)';
                Cmom_ptile_ext = prctile(NofFailure_cmom_ext',P)';
                Lmom_ptile_ext = prctile(NofFailure_lmom_ext',P)';
                
                %% Plot
                
                Cmom_10 = Cmom_ptile(:,1); Lmom_10 = Lmom_ptile(:,1);
                Cmom_10_ext = Cmom_ptile_ext(:,1); Lmom_10_ext = Lmom_ptile_ext(:,1);
                
                Cmom_50 = Cmom_ptile(:,2); Lmom_50 = Lmom_ptile(:,2);
                Cmom_50_ext = Cmom_ptile_ext(:,2); Lmom_50_ext = Lmom_ptile_ext(:,2);
                
                Cmom_90 = Cmom_ptile(:,3); Lmom_90 = Lmom_ptile(:,3);
                Cmom_90_ext = Cmom_ptile_ext(:,3); Lmom_90_ext = Lmom_ptile_ext(:,3);
                
                Legend = {'Population','10th percentile','Median','90th percentile'};
                
                Area_shade_plot(x,NofFailure_act,E,Cmom_10,Cmom_50,Cmom_90,'C-moment',Legend);
                saveas(gcf,'Area_Cmom_10_50_90_5sample.fig')
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', 'Area_Cmom_10_50_90_5sample.png', '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Lmom_10,Lmom_50,Lmom_90,'L-moment',Legend);
                saveas(gcf,'Area_Lmom_10_50_90_5sample.fig')
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', 'Area_Lmom_10_50_90_5sample.png', '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Cmom_10_ext,Cmom_50_ext,Cmom_90_ext,'C-moment',Legend);
                
                saveas(gcf,'Area_Cmom_10_50_90_5sample_ext.fig')
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', 'Area_Cmom_10_50_90_5sample_ext.png', '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Lmom_10_ext,Lmom_50_ext,Lmom_90_ext,'L-moment',Legend);
                
                saveas(gcf,'Area_Lmom_10_50_90_5sample_ext.fig')
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', 'Area_Lmom_10_50_90_5sample_ext.png', '-r100')
                
            case 2575
                P = [25,50,75];
                Cmom_ptile = prctile(NofFailure_cmom',P)';
                Lmom_ptile = prctile(NofFailure_lmom',P)';
                Cmom_ptile_ext = prctile(NofFailure_cmom_ext',P)';
                Lmom_ptile_ext = prctile(NofFailure_lmom_ext',P)';
                
                %% Plot
                
                Cmom_25 = Cmom_ptile(:,1); Lmom_25 = Lmom_ptile(:,1);
                Cmom_25_ext = Cmom_ptile_ext(:,1); Lmom_25_ext = Lmom_ptile_ext(:,1);
                
                Cmom_50 = Cmom_ptile(:,2); Lmom_50 = Lmom_ptile(:,2);
                Cmom_50_ext = Cmom_ptile_ext(:,2); Lmom_50_ext = Lmom_ptile_ext(:,2);
                
                Cmom_75 = Cmom_ptile(:,3); Lmom_75 = Lmom_ptile(:,3);
                Cmom_75_ext = Cmom_ptile_ext(:,3); Lmom_75_ext = Lmom_ptile_ext(:,3);
                
                Legend = {'Population','25th percentile','Median','75th percentile'};
                
                Area_shade_plot(x,NofFailure_act,E,Cmom_25,Cmom_50,Cmom_75,'C-moment',Legend);
                saveas(gcf,['Area_Cmom_25_50_75_',num2str(N),'_sample.fig'])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', ['Area_Cmom_25_50_75_',num2str(N),'_sample.png'], '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Lmom_25,Lmom_50,Lmom_75,'L-moment',Legend);
                saveas(gcf,['Area_lmom_25_50_75_',num2str(N),'_sample.fig'])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', ['Area_lmom_25_50_75_',num2str(N),'_sample.png'], '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Cmom_25_ext,Cmom_50_ext,Cmom_75_ext,'C-moment',Legend);
                
                saveas(gcf,['Area_Cmom_25_50_75_',num2str(N),'_sample_ext.fig'])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', ['Area_Cmom_25_50_75_',num2str(N),'_sample_ext.png'], '-r100')
                
                
                Area_shade_plot(x,NofFailure_act,E,Lmom_25_ext,Lmom_50_ext,Lmom_75_ext,'L-moment',Legend);
                
                saveas(gcf,['Area_lmom_25_50_75_',num2str(N),'_sample_ext.fig'])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0  4.73 4.42])
                print('-dpng', ['Area_lmom_25_50_75_',num2str(N),'_sample_ext.png'], '-r100')
        end
    end
end



