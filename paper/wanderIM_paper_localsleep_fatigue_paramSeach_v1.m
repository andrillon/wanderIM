%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))

%% load data
prticle_Thr_all=[80 90 95 99];
AmpCriterion_all=[4 9]; AmpCriterion_labels={'P2P','NP'};
thrElec_all={'byCh','allCh'};
for nThr=1:length(prticle_Thr_all)
    for nCrit=1:length(AmpCriterion_all)
        for nCrit2=1:length(thrElec_all)
            
            filename = [root_path filesep 'hddm' filesep 'HDDM_WIM_localsleep_amp_pup_Dec21_perprobe_thr' num2str(prticle_Thr_all(nThr)) '_' thrElec_all{nCrit2} '_' AmpCriterion_labels{nCrit} '_v5.txt'];
            delimiter = ',';
            startRow = 2;
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            res_table = table;
            res_table.SubID = dataArray{:, 1};
            res_table.ContProbeN = dataArray{:, 2};
            res_table.BlockN = dataArray{:, 3};
            res_table.ProbeN = dataArray{:, 4};
            res_table.Task = dataArray{:, 5};
            res_table.W_Fp1 = dataArray{:, 6};
            res_table.W_Fz = dataArray{:, 7};
            res_table.W_F3 = dataArray{:, 8};
            res_table.W_F7 = dataArray{:, 9};
            res_table.W_FT9 = dataArray{:, 10};
            res_table.W_FC5 = dataArray{:, 11};
            res_table.W_FC1 = dataArray{:, 12};
            res_table.W_C3 = dataArray{:, 13};
            res_table.W_T7 = dataArray{:, 14};
            res_table.W_TP9 = dataArray{:, 15};
            res_table.W_CP5 = dataArray{:, 16};
            res_table.W_CP1 = dataArray{:, 17};
            res_table.W_Pz = dataArray{:, 18};
            res_table.W_P3 = dataArray{:, 19};
            res_table.W_P7 = dataArray{:, 20};
            res_table.W_O1 = dataArray{:, 21};
            res_table.W_Oz = dataArray{:, 22};
            res_table.W_O2 = dataArray{:, 23};
            res_table.W_P4 = dataArray{:, 24};
            res_table.W_P8 = dataArray{:, 25};
            res_table.W_TP10 = dataArray{:, 26};
            res_table.W_CP6 = dataArray{:, 27};
            res_table.W_CP2 = dataArray{:, 28};
            res_table.W_Cz = dataArray{:, 29};
            res_table.W_C4 = dataArray{:, 30};
            res_table.W_T8 = dataArray{:, 31};
            res_table.W_FT10 = dataArray{:, 32};
            res_table.W_FC6 = dataArray{:, 33};
            res_table.W_FC2 = dataArray{:, 34};
            res_table.W_F4 = dataArray{:, 35};
            res_table.W_F8 = dataArray{:, 36};
            res_table.W_Fp2 = dataArray{:, 37};
            res_table.W_AF7 = dataArray{:, 38};
            res_table.W_AF3 = dataArray{:, 39};
            res_table.W_F1 = dataArray{:, 40};
            res_table.W_F5 = dataArray{:, 41};
            res_table.W_FT7 = dataArray{:, 42};
            res_table.W_FC3 = dataArray{:, 43};
            res_table.W_C1 = dataArray{:, 44};
            res_table.W_C5 = dataArray{:, 45};
            res_table.W_TP7 = dataArray{:, 46};
            res_table.W_CP3 = dataArray{:, 47};
            res_table.W_P1 = dataArray{:, 48};
            res_table.W_P5 = dataArray{:, 49};
            res_table.W_PO7 = dataArray{:, 50};
            res_table.W_PO3 = dataArray{:, 51};
            res_table.W_POz = dataArray{:, 52};
            res_table.W_PO4 = dataArray{:, 53};
            res_table.W_PO8 = dataArray{:, 54};
            res_table.W_P6 = dataArray{:, 55};
            res_table.W_P2 = dataArray{:, 56};
            res_table.W_CPz = dataArray{:, 57};
            res_table.W_CP4 = dataArray{:, 58};
            res_table.W_TP8 = dataArray{:, 59};
            res_table.W_C6 = dataArray{:, 60};
            res_table.W_C2 = dataArray{:, 61};
            res_table.W_FC4 = dataArray{:, 62};
            res_table.W_FT8 = dataArray{:, 63};
            res_table.W_F6 = dataArray{:, 64};
            res_table.W_AF8 = dataArray{:, 65};
            res_table.W_AF4 = dataArray{:, 66};
            res_table.W_F2 = dataArray{:, 67};
            res_table.W_FCz = dataArray{:, 68};
            res_table.A_Fp1 = dataArray{:, 69};
            res_table.A_Fz = dataArray{:, 70};
            res_table.A_F3 = dataArray{:, 71};
            res_table.A_F7 = dataArray{:, 72};
            res_table.A_FT9 = dataArray{:, 73};
            res_table.A_FC5 = dataArray{:, 74};
            res_table.A_FC1 = dataArray{:, 75};
            res_table.A_C3 = dataArray{:, 76};
            res_table.A_T7 = dataArray{:, 77};
            res_table.A_TP9 = dataArray{:, 78};
            res_table.A_CP5 = dataArray{:, 79};
            res_table.A_CP1 = dataArray{:, 80};
            res_table.A_Pz = dataArray{:, 81};
            res_table.A_P3 = dataArray{:, 82};
            res_table.A_P7 = dataArray{:, 83};
            res_table.A_O1 = dataArray{:, 84};
            res_table.A_Oz = dataArray{:, 85};
            res_table.A_O2 = dataArray{:, 86};
            res_table.A_P4 = dataArray{:, 87};
            res_table.A_P8 = dataArray{:, 88};
            res_table.A_TP10 = dataArray{:, 89};
            res_table.A_CP6 = dataArray{:, 90};
            res_table.A_CP2 = dataArray{:, 91};
            res_table.A_Cz = dataArray{:, 92};
            res_table.A_C4 = dataArray{:, 93};
            res_table.A_T8 = dataArray{:, 94};
            res_table.A_FT10 = dataArray{:, 95};
            res_table.A_FC6 = dataArray{:, 96};
            res_table.A_FC2 = dataArray{:, 97};
            res_table.A_F4 = dataArray{:, 98};
            res_table.A_F8 = dataArray{:, 99};
            res_table.A_Fp2 = dataArray{:, 100};
            res_table.A_AF7 = dataArray{:, 101};
            res_table.A_AF3 = dataArray{:, 102};
            res_table.A_F1 = dataArray{:, 103};
            res_table.A_F5 = dataArray{:, 104};
            res_table.A_FT7 = dataArray{:, 105};
            res_table.A_FC3 = dataArray{:, 106};
            res_table.A_C1 = dataArray{:, 107};
            res_table.A_C5 = dataArray{:, 108};
            res_table.A_TP7 = dataArray{:, 109};
            res_table.A_CP3 = dataArray{:, 110};
            res_table.A_P1 = dataArray{:, 111};
            res_table.A_P5 = dataArray{:, 112};
            res_table.A_PO7 = dataArray{:, 113};
            res_table.A_PO3 = dataArray{:, 114};
            res_table.A_POz = dataArray{:, 115};
            res_table.A_PO4 = dataArray{:, 116};
            res_table.A_PO8 = dataArray{:, 117};
            res_table.A_P6 = dataArray{:, 118};
            res_table.A_P2 = dataArray{:, 119};
            res_table.A_CPz = dataArray{:, 120};
            res_table.A_CP4 = dataArray{:, 121};
            res_table.A_TP8 = dataArray{:, 122};
            res_table.A_C6 = dataArray{:, 123};
            res_table.A_C2 = dataArray{:, 124};
            res_table.A_FC4 = dataArray{:, 125};
            res_table.A_FT8 = dataArray{:, 126};
            res_table.A_F6 = dataArray{:, 127};
            res_table.A_AF8 = dataArray{:, 128};
            res_table.A_AF4 = dataArray{:, 129};
            res_table.A_F2 = dataArray{:, 130};
            res_table.A_FCz = dataArray{:, 131};
            res_table.MA_Fp1 = dataArray{:, 132};
            res_table.MA_Fz = dataArray{:, 133};
            res_table.MA_F3 = dataArray{:, 134};
            res_table.MA_F7 = dataArray{:, 135};
            res_table.MA_FT9 = dataArray{:, 136};
            res_table.MA_FC5 = dataArray{:, 137};
            res_table.MA_FC1 = dataArray{:, 138};
            res_table.MA_C3 = dataArray{:, 139};
            res_table.MA_T7 = dataArray{:, 140};
            res_table.MA_TP9 = dataArray{:, 141};
            res_table.MA_CP5 = dataArray{:, 142};
            res_table.MA_CP1 = dataArray{:, 143};
            res_table.MA_Pz = dataArray{:, 144};
            res_table.MA_P3 = dataArray{:, 145};
            res_table.MA_P7 = dataArray{:, 146};
            res_table.MA_O1 = dataArray{:, 147};
            res_table.MA_Oz = dataArray{:, 148};
            res_table.MA_O2 = dataArray{:, 149};
            res_table.MA_P4 = dataArray{:, 150};
            res_table.MA_P8 = dataArray{:, 151};
            res_table.MA_TP10 = dataArray{:, 152};
            res_table.MA_CP6 = dataArray{:, 153};
            res_table.MA_CP2 = dataArray{:, 154};
            res_table.MA_Cz = dataArray{:, 155};
            res_table.MA_C4 = dataArray{:, 156};
            res_table.MA_T8 = dataArray{:, 157};
            res_table.MA_FT10 = dataArray{:, 158};
            res_table.MA_FC6 = dataArray{:, 159};
            res_table.MA_FC2 = dataArray{:, 160};
            res_table.MA_F4 = dataArray{:, 161};
            res_table.MA_F8 = dataArray{:, 162};
            res_table.MA_Fp2 = dataArray{:, 163};
            res_table.MA_AF7 = dataArray{:, 164};
            res_table.MA_AF3 = dataArray{:, 165};
            res_table.MA_F1 = dataArray{:, 166};
            res_table.MA_F5 = dataArray{:, 167};
            res_table.MA_FT7 = dataArray{:, 168};
            res_table.MA_FC3 = dataArray{:, 169};
            res_table.MA_C1 = dataArray{:, 170};
            res_table.MA_C5 = dataArray{:, 171};
            res_table.MA_TP7 = dataArray{:, 172};
            res_table.MA_CP3 = dataArray{:, 173};
            res_table.MA_P1 = dataArray{:, 174};
            res_table.MA_P5 = dataArray{:, 175};
            res_table.MA_PO7 = dataArray{:, 176};
            res_table.MA_PO3 = dataArray{:, 177};
            res_table.MA_POz = dataArray{:, 178};
            res_table.MA_PO4 = dataArray{:, 179};
            res_table.MA_PO8 = dataArray{:, 180};
            res_table.MA_P6 = dataArray{:, 181};
            res_table.MA_P2 = dataArray{:, 182};
            res_table.MA_CPz = dataArray{:, 183};
            res_table.MA_CP4 = dataArray{:, 184};
            res_table.MA_TP8 = dataArray{:, 185};
            res_table.MA_C6 = dataArray{:, 186};
            res_table.MA_C2 = dataArray{:, 187};
            res_table.MA_FC4 = dataArray{:, 188};
            res_table.MA_FT8 = dataArray{:, 189};
            res_table.MA_F6 = dataArray{:, 190};
            res_table.MA_AF8 = dataArray{:, 191};
            res_table.MA_AF4 = dataArray{:, 192};
            res_table.MA_F2 = dataArray{:, 193};
            res_table.MA_FCz = dataArray{:, 194};
            res_table.State = dataArray{:, 195};
            res_table.Vig = dataArray{:, 196};
            res_table.Pup = dataArray{:, 197};
            res_table.pcPup = dataArray{:, 198};
            
            clearvars filename delimiter startRow formatSpec fileID dataArray ans;
            
            res_table.W_all=nanmean(table2array(res_table(:,0*63+[6:68])),2);
            res_table.A_all=nanmean(table2array(res_table(:,1*63+[6:68])),2);
            res_table.MA_all=nanmean(table2array(res_table(:,2*63+[6:68])),2);
            
            %%
            % clean 334
%             warning('removing 334')
            res_table(res_table.SubID==334,:)=[];
            
            res_table.State=categorical(res_table.State);
            res_table.SubID=categorical(res_table.SubID);
            res_table.Task=categorical(res_table.Task);
            
            %% LME models
%             mdl_0 = fitlme(res_table,sprintf('W_all~1+(1|SubID)'));
%             mdl_1 = fitlme(res_table,sprintf('W_all~1+Task+(1|SubID)'));
%             mdl_2a= fitlme(res_table,sprintf('W_all~1+Task+pcPup+(1|SubID)'));
            mdl_2b= fitlme(res_table,sprintf('W_all~1+Task+Vig  +(1|SubID)'));
            mdl_2c= fitlme(res_table,sprintf('W_all~1+Task+State  +(1|SubID)'));
            
            fprintf('... ... %g thr on %s for %s: %g waves across electrodes per probes\n- effect of vigilance:\n\tEstimate = %g\n\tt-Value = %g\n\tp-Value = %g\n',...
                prticle_Thr_all(nThr),thrElec_all{nCrit2},AmpCriterion_labels{nCrit},nanmean(sum(table2array(res_table(:,6:68)),2)),...
                double(mdl_2b.Coefficients(3,2)),double(mdl_2b.Coefficients(3,4)),double(mdl_2b.Coefficients(3,6)));
  
            fprintf('- effect of MW:\n\tEstimate = %g\n\tt-Value = %g\n\tp-Value = %g\n- effect of MB:\n\tEstimate = %g\n\tt-Value = %g\n\tp-Value = %g\n\n',...
                double(mdl_2c.Coefficients(3,2)),double(mdl_2c.Coefficients(3,4)),double(mdl_2c.Coefficients(3,6)),...
                double(mdl_2c.Coefficients(4,2)),double(mdl_2c.Coefficients(4,4)),double(mdl_2c.Coefficients(4,6)));
          
            mdl_results(nThr,nCrit,nCrit2)=double(mdl_2b.Coefficients(3,4));
            mdl_results2(nThr,nCrit,nCrit2,1)=double(mdl_2c.Coefficients(3,4));
            mdl_results2(nThr,nCrit,nCrit2,2)=double(mdl_2c.Coefficients(4,4));
        end
    end
end

%%
figure; 
for nCrit2=1:2
subplot(1,2,nCrit2); format_fig;
imagesc(squeeze(mdl_results(:,:,nCrit2))); colorbar;
title(thrElec_all{nCrit2})
xlabel('Amp Crit')
ylabel('Thr')
set(gca,'ytick',1:4,'YTickLabel',prticle_Thr_all,'xtick',1:2,'xTickLabel',AmpCriterion_labels)
caxis([-1 1]*7)
end

%%
figure; 
for nCrit2=1:2
subplot(1,2,nCrit2); format_fig;
imagesc(squeeze(mdl_results2(:,:,nCrit2,1))); colorbar;
title(thrElec_all{nCrit2})
xlabel('Amp Crit')
ylabel('Thr')
set(gca,'ytick',1:4,'YTickLabel',prticle_Thr_all,'xtick',1:2,'xTickLabel',AmpCriterion_labels)
caxis([-1 1]*5)
end

%%
figure; 
for nCrit2=1:2
subplot(1,2,nCrit2); format_fig;
imagesc(squeeze(mdl_results2(:,:,nCrit2,2))); colorbar;
title(thrElec_all{nCrit2})
xlabel('Amp Crit')
ylabel('Thr')
set(gca,'ytick',1:4,'YTickLabel',prticle_Thr_all,'xtick',1:2,'xTickLabel',AmpCriterion_labels)
caxis([-1 1]*5)
end