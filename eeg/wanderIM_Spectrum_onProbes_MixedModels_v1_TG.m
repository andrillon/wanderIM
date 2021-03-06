%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Teigane's general notes:
%D.conditions - probe conditions
%D.chantype - shows where each channel is coming from - general)
%D.channels - all 66 channels with information about them
%D.fsample - sampling rate
%D.time - time within the epoch - we have from -32 seconds to 32 seconds. - this is the variable we will need to modify to get the right time window for our analysis

%D.indsample - (stands for index of the sample) is a function that will cut
%the data to the time window I need. To be modified when looking at
%baseline vs around probe (baseline cut into 3 10s blocks)

%% Init
clear all;
close all;
try
    run localdef_wanderIM
catch
    run ../localdef_wanderIM
end
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_probes_mat=[];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    % I probably won't use the left and right flickers individually for now
    % but note: the left and right are showing inverse to logic (needs to
    % be checked)
    left_freq=SubjectInfo.FlickerL;
    right_freq=SubjectInfo.FlickerR;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1; %20s window
%    these_times=D.indsample(-10):D.indsample(0)-1; %10s window
    temp_data=D(:,these_times,:); % D contains the data with channels * time * trials
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);

    %%% aggregate trials between probes
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
%             if npr==1
%                 last_pr_tridx=0;
%             else
%                 last_pr_tridx=these_probes(npr-1,6);
%             end  
            last_pr_tridx=this_pr_tridx-20;
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look            
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
%             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            temp_corr_go=nanmean(temp_testres(:,12));
            temp_corr_nogo=nanmean(temp_testres(:,11));
            temp_rt_go=nanmean(temp_testres(:,10)-temp_testres(:,8));
            % Colum order:
            
            temp_behav=[str2num(SubID) nbl these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go left_freq right_freq];
            
            % select probe
            this_trial_label=sprintf('B%g_P%g_',nbl,npr);
            idx_probe_eeg=find_trials(D.conditions,this_trial_label);
            if length(idx_probe_eeg)~=1
                warning('Problem')
                continue
            end
            
            FOI=[6 7.5 12 15 13.5];
            FOIlab={'f1','f2','2f1','2f2','IM'};
            temp_EEGres=[];
            temp_EEGgroup=[];
            temp_SNR=[];
            temp_PWR=[];
            for nf=1:length(FOI)
                [~,fidx]=findclosest(faxis,FOI(nf));
                for nE=0:63
                    if nE>0
                    temp_SNR=squeeze(logSNR(nE,fidx,idx_probe_eeg));
                    temp_PWR=squeeze(logpow(nE,fidx,idx_probe_eeg));
                    temp_EEGres=[temp_EEGres ; temp_SNR];
                    temp_EEGgroup=[temp_EEGgroup ; [nf nE]];
                    else
                        temp_EEGres=[temp_EEGres ; 0];
                    temp_EEGgroup=[temp_EEGgroup ; [nf nE]];
                    end
                end
            end
%             temp_EEGres=[temp_SNR(1,:) temp_SNR(2,:) temp_PWR(1,:) temp_PWR(2,:)];
            
            all_probes_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq','SNR','Freq','Chan'};
            
            all_probes_mat=[all_probes_mat ; [repmat(temp_behav,length(temp_EEGres),1) temp_EEGres temp_EEGgroup]];
        end
    end
end



%% transform into tables and export
% all_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'};

tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.L_freq=categorical(tbl_probe.L_freq);
tbl_probe.R_freq=categorical(tbl_probe.R_freq);
tbl_probe.Freq=categorical(tbl_probe.Freq);
tbl_probe.Chan=categorical(tbl_probe.Chan);

writetable(tbl_probe,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh.txt']);

%%
tbl_probe(tbl_probe.State=="4",:)=[];
tbl_probe.State=removecats(tbl_probe.State);
tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});

% model comparison
% mdl_Full= fitlme(tbl_probe,'SNR~Task * Chan * State * Vig * Freq + (1| SubID)');

%%
pF=tbl_probe(tbl_probe.Freq == "1" | tbl_probe.Freq == "2",:);
p2F=tbl_probe(tbl_probe.Freq == "3" | tbl_probe.Freq == "4",:);
pIM=tbl_probe(tbl_probe.Freq == "5",:);

pF_F=tbl_probe(tbl_probe.Task == "1" & (tbl_probe.Freq == "1" | tbl_probe.Freq == "2"),:);
p2F_F=tbl_probe(tbl_probe.Task == "1" & (tbl_probe.Freq == "3" | tbl_probe.Freq == "4"),:);
pIM_F=tbl_probe(tbl_probe.Task == "1" & tbl_probe.Freq == "5",:);

pF_D=tbl_probe(tbl_probe.Task == "2" & (tbl_probe.Freq == "1" | tbl_probe.Freq == "2"),:);
p2F_D=tbl_probe(tbl_probe.Task == "2" & (tbl_probe.Freq == "3" | tbl_probe.Freq == "4"),:);
pIM_D=tbl_probe(tbl_probe.Task == "2" & tbl_probe.Freq == "5",:);


pF.Freq=removecats(pF.Freq);
p2F.Freq=removecats(p2F.Freq);
pIM.Freq=removecats(pIM.Freq);

pF_F.Freq=removecats(pF_F.Freq);
p2F_F.Freq=removecats(p2F_F.Freq);
pIM_F.Freq=removecats(pIM_F.Freq);

pF_D.Freq=removecats(pF_D.Freq);
p2F_D.Freq=removecats(p2F_D.Freq);
pIM_D.Freq=removecats(pIM_D.Freq);

pF2=pF;
pIM2=pIM;
p2F2=p2F;
pF2.State=reordercats(pF2.State,{'2','1','3'});
pIM2.State=reordercats(pIM2.State,{'2','1','3'});
p2F2.State=reordercats(p2F2.State,{'2','1','3'});

pF2_F=pF_F;
pIM2_F=pIM_F;
p2F2_F=p2F_F;
pF2_F.State=reordercats(pF2_F.State,{'2','1','3'});
pIM2_F.State=reordercats(pIM2_F.State,{'2','1','3'});
p2F2_F.State=reordercats(p2F2_F.State,{'2','1','3'});

pF2_D=pF_D;
pIM2_D=pIM_D;
p2F2_D=p2F_D;
pF2_D.State=reordercats(pF2_D.State,{'2','1','3'});
pIM2_D.State=reordercats(pIM2_D.State,{'2','1','3'});
p2F2_D.State=reordercats(p2F2_D.State,{'2','1','3'});

% Overall models
mdl_IM= fitlme(pIM,'SNR~Chan * State + (1| SubID)');
mdl_F= fitlme(pF,'SNR~Chan * State + (1| SubID)');
mdl_2F= fitlme(p2F,'SNR~Chan * State + (1| SubID)');

mdl2_IM= fitlme(pIM2,'SNR~Chan * State + (1| SubID)');
mdl2_F= fitlme(pF2,'SNR~Chan * State + (1| SubID)');
mdl2_2F= fitlme(p2F2,'SNR~Chan * State + (1| SubID)');

% Models in Face Task
mdl_IM_F= fitlme(pIM_F,'SNR~Chan * State + (1| SubID)');
mdl_F_F= fitlme(pF_F,'SNR~Chan * State + (1| SubID)');
mdl_2F_F= fitlme(p2F_F,'SNR~Chan * State + (1| SubID)');

mdl2_IM_F= fitlme(pIM2_F,'SNR~Chan * State + (1| SubID)');
mdl2_F_F= fitlme(pF2_F,'SNR~Chan * State + (1| SubID)');
mdl2_2F_F= fitlme(p2F2_F,'SNR~Chan * State + (1| SubID)');

% Models in Digit Task
mdl_IM_D= fitlme(pIM_D,'SNR~Chan * State + (1| SubID)');
mdl_F_D= fitlme(pF_D,'SNR~Chan * State + (1| SubID)');
mdl_2F_D= fitlme(p2F_D,'SNR~Chan * State + (1| SubID)');

mdl2_IM_D= fitlme(pIM2_D,'SNR~Chan * State + (1| SubID)');
mdl2_F_D= fitlme(pF2_D,'SNR~Chan * State + (1| SubID)');
mdl2_2F_D= fitlme(p2F2_D,'SNR~Chan * State + (1| SubID)');



%% IM
pVal=double(mdl_IM.Coefficients(:,6));
beta=double(mdl_IM.Coefficients(:,2));
MWvsONidx=find_trials(mdl_IM.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_IM.CoefficientNames,'State_3:');

pVal2=double(mdl2_IM.Coefficients(:,6));
beta2=double(mdl2_IM.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_IM.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - IM')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - IM')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - IM')
caxis([-1 1]*0.5)
colorbar;
%% F
pVal=double(mdl_F.Coefficients(:,6));
beta=double(mdl_F.Coefficients(:,2));
MWvsONidx=find_trials(mdl_F.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_F.CoefficientNames,'State_3:');

pVal2=double(mdl2_F.Coefficients(:,6));
beta2=double(mdl2_F.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_F.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - F')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - F')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - F')
caxis([-1 1]*0.5)
colorbar;


%% 2F
pVal=double(mdl_2F.Coefficients(:,6));
beta=double(mdl_2F.Coefficients(:,2));
MWvsONidx=find_trials(mdl_2F.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_2F.CoefficientNames,'State_3:');

pVal2=double(mdl2_2F.Coefficients(:,6));
beta2=double(mdl2_2F.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_2F.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - 2F')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - 2F')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - 2F')
caxis([-1 1]*0.5)
colorbar;

%% IM - Faces
pVal=double(mdl_IM_F.Coefficients(:,6));
beta=double(mdl_IM_F.Coefficients(:,2));
MWvsONidx=find_trials(mdl_IM_F.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_IM_F.CoefficientNames,'State_3:');

pVal2=double(mdl2_IM_F.Coefficients(:,6));
beta2=double(mdl2_IM_F.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_IM_F.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig;     
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - IM - Faces')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - IM - Faces')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - IM - Faces')
caxis([-1 1]*0.5)
colorbar;
%% F - Faces
pVal=double(mdl_F_F.Coefficients(:,6));
beta=double(mdl_F_F.Coefficients(:,2));
MWvsONidx=find_trials(mdl_F_F.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_F_F.CoefficientNames,'State_3:');

pVal2=double(mdl2_F_F.Coefficients(:,6));
beta2=double(mdl2_F_F.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_F_F.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - F - Faces')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - F - Faces')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - F - Faces')
caxis([-1 1]*0.5)
colorbar;


%% 2F - Faces
pVal=double(mdl_2F_F.Coefficients(:,6));
beta=double(mdl_2F_F.Coefficients(:,2));
MWvsONidx=find_trials(mdl_2F_F.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_2F_F.CoefficientNames,'State_3:');

pVal2=double(mdl2_2F_F.Coefficients(:,6));
beta2=double(mdl2_2F_F.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_2F_F.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - 2F - Faces')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - 2F - Faces')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - 2F - Faces')
caxis([-1 1]*0.5)
colorbar;

%% IM - Digits
pVal=double(mdl_IM_D.Coefficients(:,6));
beta=double(mdl_IM_D.Coefficients(:,2));
MWvsONidx=find_trials(mdl_IM_D.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_IM_D.CoefficientNames,'State_3:');

pVal2=double(mdl2_IM_D.Coefficients(:,6));
beta2=double(mdl2_IM_D.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_IM_D.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - IM - Digits')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - IM - Digits')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - IM - Digits')
caxis([-1 1]*0.5)
colorbar;
%% F - Digits
pVal=double(mdl_F_D.Coefficients(:,6));
beta=double(mdl_F_D.Coefficients(:,2));
MWvsONidx=find_trials(mdl_F_D.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_F_D.CoefficientNames,'State_3:');

pVal2=double(mdl2_F_D.Coefficients(:,6));
beta2=double(mdl2_F_D.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_F_D.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - F - Digits')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - F - Digits')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - F - Digits')
caxis([-1 1]*0.5)
colorbar;


%% 2F - Digits
pVal=double(mdl_2F_D.Coefficients(:,6));
beta=double(mdl_2F_D.Coefficients(:,2));
MWvsONidx=find_trials(mdl_2F_D.CoefficientNames,'State_2:');
MBvsONidx=find_trials(mdl_2F_D.CoefficientNames,'State_3:');

pVal2=double(mdl2_2F_D.Coefficients(:,6));
beta2=double(mdl2_2F_D.Coefficients(:,2));
MWvsMBidx2=find_trials(mdl2_2F_D.CoefficientNames,'State_3:');

load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
subplot(1,3,1); format_fig; 
temp_topo=beta(MWvsONidx);
temp_pV=pVal(MWvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs ON - 2F - Digits')
caxis([-1 1]*0.5)

subplot(1,3,2); format_fig; 
temp_topo=beta(MBvsONidx);
temp_pV=pVal(MBvsONidx);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MB vs ON - 2F - Digits')
caxis([-1 1]*0.5)

subplot(1,3,3); format_fig; 
temp_topo=beta2(MWvsMBidx2);
temp_pV=pVal2(MWvsMBidx2);
temp_topo(temp_pV>0.05)=0;

simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1])
title('MW vs MB - 2F - Digits')
caxis([-1 1]*0.5)
colorbar;