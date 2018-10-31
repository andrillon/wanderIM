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
run ../localdef_wanderIM

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
IMclusterE=[19 20 50 51 54];
F2clusterE=[44 45 46 47];
all_cluster_SNR=[];
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
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);

    %%% aggregate trials between probes
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end  
%             last_pr_tridx=this_pr_tridx-20;
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look            
            
            temp_testresgo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testresgo(end-19:end,12))';%/corr_go(n,these_probes(npr,5));
            temp_testresnogo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testresnogo(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            
%             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            temp_corr_go=nanmean(tcorr_go);
            temp_corr_nogo=nanmean(tcorr_nogo);
            temp_rt_go=nanmean(temp_testresgo(:,10)-temp_testresgo(:,8));
            temp_corr_go2=tcorr_go;
            temp_corr_nogo2=tcorr_nogo;
            [dprime, crit]=calc_dprime(temp_corr_go2,temp_corr_nogo2==0);
            % Colum order
            
            temp_behav=[str2num(SubID) nbl npr these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go left_freq right_freq dprime crit];
            
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
            
            
            % add alpha
            for nE=0:63
                if nE>0
                    temp_PWR=polyarea([faxis(faxis==8) faxis((faxis>=8 & faxis<10.3) | (faxis>10.7 & faxis<=11.5)) faxis(faxis==11.5)],...
                        [0 squeeze(logpow(nE,(faxis>=8 & faxis<10.3) | (faxis>10.7 & faxis<=11.5) ,idx_probe_eeg)) 0]);
                    temp_PWR0=polyarea([faxis(faxis==8) faxis(faxis==8) faxis(faxis==11.5) faxis(faxis==11.5)],...
                        [0 squeeze(logpow(nE,(faxis==8) ,idx_probe_eeg)) squeeze(logpow(nE,(faxis==11.5) ,idx_probe_eeg)) 0]);
                    temp_EEGres=[temp_EEGres ; temp_PWR-temp_PWR0];
                    temp_EEGgroup=[temp_EEGgroup ; [length(FOI)+1 nE]];
                else
                    temp_EEGres=[temp_EEGres ; 0];
                    temp_EEGgroup=[temp_EEGgroup ; [length(FOI)+1 nE]];
                end
            end
            
            all_probes_headers={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq','dp','crit','SNR','Freq','Chan'};
            
            all_probes_mat=[all_probes_mat ; [repmat(temp_behav,length(temp_EEGres),1) temp_EEGres temp_EEGgroup]];
            
            [~,idxIM]=findclosest(faxis,FOI(nf));
            cluster_IM=squeeze(mean(logSNR(IMclusterE,idxIM,idx_probe_eeg),1));
            [~,idx2F1]=findclosest(faxis,FOI(4));
            [~,idx2F2]=findclosest(faxis,FOI(3));
            cluster_2f=squeeze(mean(mean(logSNR(F2clusterE,[idx2F1 idx2F2],idx_probe_eeg),2),1));
            all_cluster_SNR=[all_cluster_SNR ; [temp_behav cluster_IM cluster_2f]];
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


tbl_probe(tbl_probe.State=="4",:)=[];
tbl_probe.State=removecats(tbl_probe.State);
tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});

writetable(tbl_probe,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh3.txt']);

%% acorss all Electrodes
mat_mdl_comp=[];
fprintf('E%2.0f\n')
for nE=1:63
fprintf('\b\b\b\bE%2.0f\n',nE)
    temp_tbl=tbl_probe(tbl_probe.Chan==num2str(nE) & (tbl_probe.Freq == "1" | tbl_probe.Freq == "2"),:);
    temp_tbl.Freq=removecats(temp_tbl.Freq);
    mdl_IM1= fitlme(temp_tbl,'SNR~nBlock + Task * Freq + (1| SubID)');
    mdl_IM2= fitlme(temp_tbl,'SNR~nBlock + State * Task * Freq + (1| SubID)');
    mdl_comp=compare(mdl_IM1,mdl_IM2);
    mat_F_comp(nE,:)=double(mdl_comp(2,6:8));
    
    temp_tbl=tbl_probe(tbl_probe.Chan==num2str(nE) & (tbl_probe.Freq == "3" | tbl_probe.Freq == "4"),:);
    temp_tbl.Freq=removecats(temp_tbl.Freq);
    mdl_IM1= fitlme(temp_tbl,'SNR~nBlock + Task * Freq + (1| SubID)');
    mdl_IM2= fitlme(temp_tbl,'SNR~nBlock + State * Task * Freq + (1| SubID)');
    mdl_comp=compare(mdl_IM1,mdl_IM2);
    mat_2F_comp(nE,:)=double(mdl_comp(2,6:8));
    
    temp_tbl=tbl_probe(tbl_probe.Chan==num2str(nE) & (tbl_probe.Freq == "5"),:);
    temp_tbl.Freq=removecats(temp_tbl.Freq);
    mdl_IM1= fitlme(temp_tbl,'SNR~nBlock + Task + (1| SubID)');
    mdl_IM2= fitlme(temp_tbl,'SNR~nBlock + State * Task + (1| SubID)');
    mdl_comp=compare(mdl_IM1,mdl_IM2);
    mat_IM_comp(nE,:)=double(mdl_comp(2,6:8));
    
%     %%%% or
%      temp_tbl=tbl_probe(tbl_probe.Chan==num2str(nE) & (tbl_probe.Freq == "1" | tbl_probe.Freq == "2"),:);
%     temp_tbl.Freq=removecats(temp_tbl.Freq);
%     mdl_IM1= fitlme(temp_tbl,'SNR~nBlock + Task + (1| SubID)');
%     mdl_IM2= fitlme(temp_tbl,'SNR~nBlock + State * Task + (1| SubID)');
%     mdl_comp=compare(mdl_IM1,mdl_IM2);
%     mat_F_comp2(nE,:)=double(mdl_comp(2,6:8));
%     
%     temp_tbl=tbl_probe(tbl_probe.Chan==num2str(nE) & (tbl_probe.Freq == "3" | tbl_probe.Freq == "4"),:);
%     temp_tbl.Freq=removecats(temp_tbl.Freq);
%     mdl_IM1= fitlme(temp_tbl,'SNR~nBlock + Task + (1| SubID)');
%     mdl_IM2= fitlme(temp_tbl,'SNR~nBlock + State * Task + (1| SubID)');
%     mdl_comp=compare(mdl_IM1,mdl_IM2);
%     mat_2F_comp2(nE,:)=double(mdl_comp(2,6:8));
end

%% plot
filename = '/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh3_results.txt';
data=load_res_csv_comparemodels(filename);
% stat_thr=fdr([data(:,3)],0.1);
scaleVal=40;
load(['EasyCap64_layout'])

figure;
titlesplots={'F','2F','IM'};
for nplot=1:3
    subplot(1,3,nplot); format_fig;
    tempplot=data(data(:,5)==nplot,2);
    tempplotpV=data(data(:,5)==nplot,3);
    stat_thr=fdr(tempplotpV,0.05);
    addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','w',24,2},'whitebk','on');
    rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    caxis([-1 1]*scaleVal)
    title(titlesplots{nplot})
    colorbar;
    
    EOI=find(tempplotpV<stat_thr)
    for nstate=1:3
        for ntask=1:2
            for nFreq=1:5
                if ~isempty(EOI)
            tempVal=tbl_probe.SNR(ismember(double(tbl_probe.Chan),EOI) & tbl_probe.Freq==num2str(nFreq) & tbl_probe.Task==num2str(ntask) & tbl_probe.State==num2str(nstate));
            meanTag(ntask,nstate,nFreq)=mean(tempVal);
                else
            meanTag(ntask,nstate,nFreq)=NaN;
                end
            end
        end
    end
end

%%
pF=tbl_probe(tbl_probe.Freq == "1" | tbl_probe.Freq == "2",:);
p2F=tbl_probe(tbl_probe.Freq == "3" | tbl_probe.Freq == "4",:);
pIM=tbl_probe(tbl_probe.Freq == "5",:);
pA=tbl_probe(tbl_probe.Freq == "6",:);

pF.Freq=removecats(pF.Freq);
p2F.Freq=removecats(p2F.Freq);
pIM.Freq=removecats(pIM.Freq);
pA.Freq=removecats(pA.Freq);

pF2=pF;
pIM2=pIM;
p2F2=p2F;
pA2=pA;
pF2.State=reordercats(pF2.State,{'2','1','3'});
pIM2.State=reordercats(pIM2.State,{'2','1','3'});
p2F2.State=reordercats(p2F2.State,{'2','1','3'});
pA2.State=reordercats(pA2.State,{'2','1','3'});

%%
mdl_IM= fitlme(pIM,'SNR~Chan * State + nBlock + (1| SubID)');
mdl_F= fitlme(pF,'SNR~Chan * State + nBlock + (1| SubID)');
mdl_2F= fitlme(p2F,'SNR~Chan * State + nBlock  +  (1| SubID)');
mdl_A= fitlme(pA,'SNR~Chan * State + nBlock  +  (1| SubID)');

mdl2_IM= fitlme(pIM2,'SNR~Chan * State + nBlock + (1| SubID)');
mdl2_F= fitlme(pF2,'SNR~Chan * State + nBlock + (1| SubID)');
mdl2_2F= fitlme(p2F2,'SNR~Chan * State + nBlock + (1| SubID)');
mdl2_A= fitlme(pA2,'SNR~Chan * State + nBlock + (1| SubID)');

fullmdl_IM= fitlme(pIM,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl_F= fitlme(pF,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl_2F= fitlme(p2F,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl_A= fitlme(pA,'SNR~Task * Chan * State + nBlock + (1| SubID)');

fullmdl2_IM= fitlme(pIM2,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl2_F= fitlme(pF2,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl2_2F= fitlme(p2F2,'SNR~Task * Chan * State + nBlock + (1| SubID)');
fullmdl2_A= fitlme(pA2,'SNR~Task * Chan * State + nBlock + (1| SubID)');

fullfullmdl_F= fitlme(pF,'SNR~L_freq * Task * Chan + nBlock + (1| SubID)');
fullfullmdl_IM= fitlme(pIM,'SNR~L_freq *Task * Chan + nBlock + (1| SubID)');

load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified

%% Check Topos
stat_thr=0.01;


this_model=fullfullmdl_F;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
ChanIdx=find_trials(this_model.CoefficientNames,'^Chan');

figure;
subplot(2,2,1); format_fig;
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('F')
colorbar;

this_model=fullfullmdl_F;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
ChanIdx=find_trials(this_model.CoefficientNames,'^L_freq_15:Chan');

subplot(2,2,2); format_fig; 
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('F*Fside')
colorbar;

this_model=fullfullmdl_IM;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
ChanIdx=find_trials(this_model.CoefficientNames,'^Chan');

subplot(2,2,3); format_fig; 
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('IM')
colorbar;

this_model=fullfullmdl_IM;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
ChanIdx=find_trials(this_model.CoefficientNames,'^L_freq_15:Chan');

subplot(2,2,4); format_fig; 
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('IM*Fside')
colorbar;

%% For Both Task
scaleVal=1.5;
stat_thr=0.05;
%%% F
this_model=fullmdl_F;
this_model2=fullmdl2_F;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
MWvsONidx=find_trials(this_model.CoefficientNames,'^State_2:');
MBvsONidx=find_trials(this_model.CoefficientNames,'^State_3:');

pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));
MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^State_3:');

figure;
subplot(3,2,1); format_fig; 
tempplot=beta(MWvsONidx);
tempplotpV=pVal(MWvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs ON - F')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,3); format_fig; 
tempplot=beta(MBvsONidx);
tempplotpV=pVal(MBvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MB vs ON - F')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,5); format_fig; 
tempplot=beta(MWvsMBidx2);
tempplotpV=pVal(MWvsMBidx2);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs MB - F')
caxis([-1 1]*scaleVal)
colorbar;

% %%% 2F
% this_model=fullmdl_2F;
% this_model2=fullmdl2_2F;
% pVal=double(this_model.Coefficients(:,6));
% beta=double(this_model.Coefficients(:,2));
% MWvsONidx=find_trials(this_model.CoefficientNames,'^State_2:');
% MBvsONidx=find_trials(this_model.CoefficientNames,'^State_3:');
% 
% pVal2=double(this_model2.Coefficients(:,6));
% beta2=double(this_model2.Coefficients(:,2));
% MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^State_3:');
% 
% subplot(3,3,2); format_fig; 
% tempplot=beta(MWvsONidx);
% tempplotpV=pVal(MWvsONidx);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MW vs ON - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;
% 
% subplot(3,3,5); format_fig; 
% tempplot=beta(MBvsONidx);
% tempplotpV=pVal(MBvsONidx);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MB vs ON - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;
% 
% subplot(3,3,8); format_fig; 
% tempplot=beta(MWvsMBidx2);
% tempplotpV=pVal(MWvsMBidx2);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MW vs MB - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;

%%% A
this_model=fullmdl_IM;
this_model2=fullmdl2_IM;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
MWvsONidx=find_trials(this_model.CoefficientNames,'^State_2:');
MBvsONidx=find_trials(this_model.CoefficientNames,'^State_3:');

pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));
MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^State_3:');

subplot(3,2,2); format_fig; 
tempplot=beta(MWvsONidx);
tempplotpV=pVal(MWvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs ON - IM')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,4); format_fig; 
tempplot=beta(MBvsONidx);
tempplotpV=pVal(MBvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MB vs ON - IM')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,6); format_fig; 
tempplot=beta(MWvsMBidx2);
tempplotpV=pVal(MWvsMBidx2);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs MB - IM')
caxis([-1 1]*scaleVal)
colorbar;

%% For Both Task - Intreactions
scaleVal=1;
stat_thr=0.05;
%%% F
this_model=fullmdl_F;
this_model2=fullmdl2_F;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
MWvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_2:');
MBvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_3:');

pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));
MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^Task_2:State_3:');

figure;
subplot(3,2,1); format_fig; 
tempplot=beta(MWvsONidx);
tempplotpV=pVal(MWvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs ON - F')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,3); format_fig; 
tempplot=beta(MBvsONidx);
tempplotpV=pVal(MBvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MB vs ON - F')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,5); format_fig; 
tempplot=beta(MWvsMBidx2);
tempplotpV=pVal(MWvsMBidx2);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs MB - F')
caxis([-1 1]*scaleVal)
colorbar;

% %%% 2F
% this_model=fullmdl_2F;
% this_model2=fullmdl2_2F;
% pVal=double(this_model.Coefficients(:,6));
% beta=double(this_model.Coefficients(:,2));
% MWvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_2:');
% MBvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_3:');
% 
% pVal2=double(this_model2.Coefficients(:,6));
% beta2=double(this_model2.Coefficients(:,2));
% MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^Task_2:State_3:');
% 
% subplot(3,3,2); format_fig; 
% tempplot=beta(MWvsONidx);
% tempplotpV=pVal(MWvsONidx);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MW vs ON - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;
% 
% subplot(3,3,5); format_fig; 
% tempplot=beta(MBvsONidx);
% tempplotpV=pVal(MBvsONidx);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MB vs ON - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;
% 
% subplot(3,3,8); format_fig; 
% tempplot=beta(MWvsMBidx2);
% tempplotpV=pVal(MWvsMBidx2);
% addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
% rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
% title('MW vs MB - 2F')
% caxis([-1 1]*scaleVal)
% colorbar;

%%% A
this_model=fullmdl_IM;
this_model2=fullmdl2_IM;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
MWvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_2:');
MBvsONidx=find_trials(this_model.CoefficientNames,'^Task_2:State_3:');

pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));
MWvsMBidx2=find_trials(this_model2.CoefficientNames,'^Task_2:State_3:');

subplot(3,2,2); format_fig; 
tempplot=beta(MWvsONidx);
tempplotpV=pVal(MWvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs ON - IM')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,4); format_fig; 
tempplot=beta(MBvsONidx);
tempplotpV=pVal(MBvsONidx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MB vs ON - IM')
caxis([-1 1]*scaleVal)
colorbar;

subplot(3,2,6); format_fig; 
tempplot=beta(MWvsMBidx2);
tempplotpV=pVal(MWvsMBidx2);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('MW vs MB - IM')
caxis([-1 1]*scaleVal)
colorbar;
