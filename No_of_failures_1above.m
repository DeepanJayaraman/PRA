clear;
NN = [5,7,9];
Nof = [0.5,1.5,2.5,3.5,4.5,5.5];
for j = 1:6
    Nfail = Nof(j);
    for i = 1:length(NN)
        N = NN(i);
        load(['31temp_failures_',num2str(N),'_samples'],...
            'NofFailure_cmom',...
            'NofFailure_lmom',...
            'NofFailure_cmom_ext',...
            'NofFailure_lmom_ext','E','Alpha');
        
        %F = find(E==31);
        N_C = sum(NofFailure_cmom >= Nfail);
        N_L = sum(NofFailure_lmom >= Nfail);
                N_C_ext = sum(NofFailure_cmom_ext >= Nfail);
        N_L_ext = sum(NofFailure_lmom_ext >= Nfail);
        
        switch N
            case 5
                N5 = [N_C N_L N_C_ext N_L_ext];
            case 7
                N7 = [N_C N_L N_C_ext N_L_ext];
            case 9
                N9 = [N_C N_L N_C_ext N_L_ext];
        end
    end
    
    switch Nfail
        case 0.5
            N_1_above = [N5;N7;N9];
        case 1.5
            N_2_above = [N5;N7;N9];
        case 2.5
            N_3_above = [N5;N7;N9];
        case 3.5
            N_4_above = [N5;N7;N9];
        case 4.5
            N_5_above = [N5;N7;N9];
        case 5.5
            N_6_above = [N5;N7;N9];
    end
    
    
    
end

save('Numberoffailures_forall_the_above.mat','N_1_above','N_2_above',...
    'N_3_above','N_4_above','N_5_above','N_6_above')