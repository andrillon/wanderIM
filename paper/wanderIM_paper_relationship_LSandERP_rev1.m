%% Initilaise - clear all variables and figures
clear all
close all

%% run local file to set up paths
run ../localdef_wanderIM
addpath(genpath(lscpTools_path)) % Thomas' general toolkit
addpath(genpath(spm12_path)) % SMP12 toolbox (EEG)

data_path=[root_path filesep 'behav/']; % path of behavioural data
eeg_path=[root_path filesep 'preproc_eeg'];
pupil_path=[root_path filesep 'eyetracker'];
eeg_path2=[root_path filesep 'preproc_ica'];

files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

load([pwd filesep '..' filesep '..' filesep 'WanderIM' filesep 'paper' filesep 'paper_SubID'])

%% Parameter slow-wave detection
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
thisChannel=17;
art_ampl=200;
max_posampl=75;
max_Freq=7;
fixThr=[];
myElecs=1:63;
Fs_EEG=500;

%%
all_onsets=[];
all_onsets_perm=[];
waves_perProbes_perElec=[];
nSc=0;
totperm=1;
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    if exist([pupil_path filesep 'wanderIM_eyelink_S' SubID '_clean.mat'])==0 || exist([eeg_path2 filesep 'triggers_S' SubID '.mat'])==0
        continue;
    end
    fprintf('... %s\n',SubID)
    if ~ismember(SubID,GoodSudID)
        continue;
    end
    %     CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
    % SART
    %  1: num block
    %  2: block cond (1: Faces / 2: Squares)
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    
    
    %%% Gathering behavioural data
    %     RTs=[RTs; test_res(:,10)-test_res(:,8)];
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
    end
    temp_perf=min(test_res(:,11:12),[],2);
    temp_cat=(test_res(:,5)==test_res(:,6));
    code_resp=nan(length(temp_perf),1);
    code_resp(temp_perf==1 & temp_cat==0)=1;
    code_resp(temp_perf==1 & temp_cat==1)=0;
    code_resp(temp_perf==0 & temp_cat==0)=0;
    code_resp(temp_perf==0 & temp_cat==1)=1;
    temp_RT=(test_res(:,10)-test_res(:,8));
    
    % Local sleep
    load([eeg_path2 filesep 'wanderIM_twa5_noica_bigwnd_' SubID])
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs_EEG);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,4)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>max_posampl| abs(all_Waves(:,15))>art_ampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,4)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>max_posampl| abs(all_Waves(:,15))>art_ampl,:)=[];
    
    
    slow_Waves=[];
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        if ~isempty(fixThr)
            thr_Wave=fixThr;
        else
            thr_Wave=prctile(thisE_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        
        slow_Waves=[slow_Waves ; (thisE_Waves(temp_p2p>thr_Wave,:))];
        
    end
    
    load([eeg_path2 filesep 'triggers_S' SubID])
    
    %%% Waves per probe
    nSc=nSc+1;
    for nProbe=1:60
        these_Waves=slow_Waves(slow_Waves(:,2)==nProbe,:);
        nout=hist(these_Waves(:,3),1:63);
        waves_perProbes_perElec=[waves_perProbes_perElec ; [repmat([n nProbe probe_res(nProbe,5) probe_res(nProbe,4) nProbe-10*(probe_res(nProbe,4)-1)],length(nout),1) (1:length(nout))' nout']];
    end
    
    temp_onsets=[];
    for nW=1:size(slow_Waves,1)
        this_onset=slow_Waves(nW,5)+start_probe(slow_Waves(nW,2));
        [closestvalue,index]=findclosest(clean_start_trial,this_onset);
        this_diff=(closestvalue-this_onset)/500;
        temp_onsets=[temp_onsets ; [n 0 slow_Waves(nW,2) slow_Waves(nW,3) slow_Waves(nW,4) this_diff]];
    end
    temp_onsets=temp_onsets(abs(temp_onsets(:,6))<0.3,:);
    
    for nProbe=1:60
        mean_perElec=[];
        for nEl=1:63
            mean_perElec(nEl)=mean(temp_onsets(temp_onsets(:,3)==nProbe & temp_onsets(:,4)==nEl,6));
        end
        all_onsets=[all_onsets ; [repmat([n 0 nProbe probe_res(nProbe,5) probe_res(nProbe,4) nProbe-10*(probe_res(nProbe,4)-1)],length(mean_perElec),1) (1:length(mean_perElec))' mean_perElec']];
    end
    
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for nperm=1:totperm
        fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',nperm,totperm)
        slow_Waves_perm=slow_Waves;
        slow_Waves_perm(:,2)=slow_Waves_perm(randperm(size(slow_Waves_perm,1)),2);
        
        
        temp_onsets=[];
        for nW=1:size(slow_Waves_perm,1)
            this_onset=slow_Waves_perm(nW,5)+start_probe(slow_Waves_perm(nW,2));
            [closestvalue,index]=findclosest(clean_start_trial,this_onset);
            this_diff=(closestvalue-this_onset)/500;
            temp_onsets=[temp_onsets ; [n 0 slow_Waves_perm(nW,2) slow_Waves_perm(nW,3) slow_Waves_perm(nW,4) this_diff]];
        end
        temp_onsets=temp_onsets(abs(temp_onsets(:,6))<0.3,:);
        
        for nProbe=1:60
            mean_perElec=[];
            for nEl=1:63
                mean_perElec(nEl)=mean(temp_onsets(temp_onsets(:,3)==nProbe & temp_onsets(:,4)==nEl,6));
            end
            all_onsets_perm=[all_onsets_perm ; [repmat([n nperm nProbe probe_res(nProbe,5) probe_res(nProbe,4) nProbe-10*(probe_res(nProbe,4)-1)],length(mean_perElec),1) (1:length(mean_perElec))' mean_perElec']];
        end
        
    end
end

%%
bins=-0.3:0.025:0.3;
uniqueS=unique(all_onsets(:,1));
for nS=1:length(uniqueS)
    nout=histc(all_onsets_perm(all_onsets_perm(:,1)==uniqueS(nS),end),bins);
    distrib_perm(nS,:)=nout/sum(nout)*100;
    
    nout=histc(all_onsets(all_onsets(:,1)==uniqueS(nS),end),bins); %& all_onsets(:,7)==match_str(layout.label,'CP4')
    distrib(nS,:)=nout/sum(nout)*100;
end

for nEl=1:63
    temp_perm=all_onsets_perm(all_onsets_perm(:,7)==nEl,end);
    temp=all_onsets(all_onsets(:,7)==nEl,end);
    
    [~,pV,~,stats]=ttest2(temp,temp_perm);
    diff_lag(nEl,:)=[pV stats.tstat];
end
%%
figure; hold on;
simpleTplot(bins,distrib_perm,0,[1 1 1]*0.7,0,'-',0.5,1,0,1,1);
simpleTplot(bins,distrib,0,[1 0 0]*0.7,0,'-',0.5,1,0,1,1);
line([0 0],ylim,'Color',[1 1 1]*0.7,'LineStyle','--');
format_fig;
ylim([0 9]);
xlim([-0.3+0.025 0.3-0.025])

export_fig([path_fig filesep 'LocalSleep_compERPs_LagWaves.fig'])
export_fig([path_fig filesep 'LocalSleep_compERPs_LagWaves.eps'],'-r 300')

%%
rmpath(genpath('/Users/tand0009/Work/local/spm12/'))
figure;
cmap=cbrewer('div','RdBu',256); cmap=flipud(cmap);
addpath((path_fieldtrip)); ft_defaults;
temp_topo=diff_lag(:,2);
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,0);
colormap(cmap);
caxis([-1 1]*7)
hb=colorbar;
format_fig;

export_fig([path_fig filesep 'LocalSleep_compERPs_LagWaves_Topo.fig'])
export_fig([path_fig filesep 'LocalSleep_compERPs_LagWaves_Topo.eps'],'-r 300')
