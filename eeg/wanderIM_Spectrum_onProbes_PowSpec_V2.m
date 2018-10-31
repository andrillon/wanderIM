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
all_probes_mat=nan(75600,31);
all_probes_mat2=nan(75600,821);
lastline=0;lastline2=0;
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
            fprintf('... Block %g Probe %g\n',nbl,npr)
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
            
            paramSpec=[];
            paramSpec.alpha_band=[8 11.5];
            paramSpec.beta_band=[15 30];
            paramSpec.StopFreqs=[1.5:1.5:30];
            alpha_theta=nan(63,1);
            vig_index=nan(63,1);
            W_index=nan(63,1);
            NREM_index=nan(63,1);
            REM_index=nan(63,1);
            a_pow=nan(63,1); t_pow=nan(63,1); d_pow=nan(63,1); s_pow=nan(63,1); b_pow=nan(63,1); 
            for nE=1:63
                [alpha_theta(nE),vig_index(nE),W_index(nE),NREM_index(nE),REM_index(nE),powbyband]=get_sleep_vig_indexes(squeeze(temp_data(nE,:,idx_probe_eeg)),D.fsample,paramSpec);
                a_pow(nE)=powbyband.alpha;
                t_pow(nE)=powbyband.theta;
                d_pow(nE)=powbyband.delta;
                s_pow(nE)=powbyband.spindle;
                b_pow(nE)=powbyband.beta;
            end
            
            all_probes_headers={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq','dp','crit','Chan',...
                'alpha_theta','vig_index','W_index','NREM_index','REM_index',...
                'alpha','theta','delta','sigma','beta'};
            all_probes_mat(lastline2+1:lastline2+size(alpha_theta,1),:)=[repmat(temp_behav,length(NREM_index),1) (1:63)' alpha_theta vig_index W_index NREM_index REM_index a_pow t_pow d_pow s_pow b_pow];
            lastline2=lastline2+size(alpha_theta,1);
            all_probes_mat2(lastline+1:lastline+size(logpow,1),:)=[repmat(temp_behav,size(logpow,1),1) (1:size(logpow,1))' squeeze(logpow(:,faxis<40,idx_probe_eeg))];
            lastline=lastline+size(logpow,1);
        end
    end
end

%%
figure; set(gcf,'Position',[680   161   570   817])
subplot(2,1,1)
for nstate=1:3
    temp=mean(all_probes_mat2(all_probes_mat2(:,7)==nstate & all_probes_mat2(:,21)==match_str(D.chanlabels,'F4'),22:end));
    simpleTplot(faxis(faxis<40),temp,0,Colors(nstate,:),0,'-',1,0.5,5,[],3);
end
xlabel('Freq (Hz)')
ylabel('Power (log)')
legend({'ON','MW','MB'})
title('Fz'); format_fig;
subplot(2,1,2)
for nstate=1:3
    temp=mean(all_probes_mat2(all_probes_mat2(:,7)==nstate & all_probes_mat2(:,21)==match_str(D.chanlabels,'Oz'),22:end));
    simpleTplot(faxis(faxis<40),temp,0,Colors(nstate,:),0,'-',1,0.5,5,[],3);
end
xlabel('Freq (Hz)')
ylabel('Power (log)')
title('Oz'); format_fig;
for nstate=1:3
    for nE=1:63
        temp_toplot{nstate}(nE,:,:)=squeeze(all_probes_mat2(all_probes_mat2(:,7)==nstate & all_probes_mat2(:,21)==nE,22:end));
    end
end

%%
freqW1=[15 30];
freqW2=[8 11.5];
figure;
load(['EasyCap64_layout'])
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
faxis2=faxis(faxis<40);
for nstate=1:3
    for nstate2=nstate:3
        subplot(3,3,3*(nstate-1)+nstate2)
        if nstate==nstate2
            tempplot1=mean(temp_toplot{nstate}(:,:,faxis2>=freqW1(1) & faxis2<=freqW1(2)),3);
            topoplot((nanmean(tempplot1,2)), layout.chaninfo);
            caxis([5 10])
else
            tempplot1=mean(temp_toplot{nstate}(:,:,faxis2>=freqW1(1) & faxis2<=freqW1(2)),3);
            tempplot2=mean(temp_toplot{nstate2}(:,:,faxis2>=freqW1(1) & faxis2<=freqW1(2)),3);
            topoplot(nanmean(tempplot1,2)-nanmean(tempplot2,2), layout.chaninfo);
                  caxis([-1 1]*1)
  end
    end
end


figure;
load(['EasyCap64_layout'])
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
faxis2=faxis(faxis<40);
for nstate=1:3
    for nstate2=nstate:3
        subplot(3,3,3*(nstate-1)+nstate2)
        if nstate==nstate2
            tempplot1=mean(temp_toplot{nstate}(:,:,faxis2>=freqW2(1) & faxis2<=freqW2(2)),3);
            topoplot((nanmean(tempplot1,2)), layout.chaninfo);
            caxis([5 10])
        else
            tempplot1=mean(temp_toplot{nstate}(:,:,faxis2>=freqW2(1) & faxis2<=freqW2(2)),3);
            tempplot2=mean(temp_toplot{nstate2}(:,:,faxis2>=freqW2(1) & faxis2<=freqW2(2)),3);
            topoplot(nanmean(tempplot1,2)-nanmean(tempplot2,2), layout.chaninfo);
            caxis([-1 1]*1)
        end
    end
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));

%% transform into tables and export
% all_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'};

tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
tbl_probe(tbl_probe.State==4,:)=[];
tbl_probe.ONOFF=nan(size(tbl_probe,1),1);
tbl_probe.MWMB=nan(size(tbl_probe,1),1);
tbl_probe.ONOFF(tbl_probe.State==1)=1;
tbl_probe.ONOFF(tbl_probe.State~=1)=0;
tbl_probe.MWMB(tbl_probe.State==2)=1;
tbl_probe.MWMB(tbl_probe.State==3)=0;
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.L_freq=categorical(tbl_probe.L_freq);
tbl_probe.R_freq=categorical(tbl_probe.R_freq);
% tbl_probe.Chan=categorical(tbl_probe.Chan);


% % writetable(tbl_probe,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh2.txt']);
%%
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
mdl_Alpha= fitlme(tbl_probe,'alpha~nBlock+nTrial+State+Task+(1| SubID)');
clear Alpha_* Beta_*
for nE=1:63
subtbl_probe=tbl_probe(tbl_probe.Chan==nE,:);
    mdl_Beta_perE= fitlme(subtbl_probe,'beta~nBlock+Task*State+(1| SubID)');
    Beta_B_tval(:,nE)=double(mdl_Beta_perE.Coefficients(:,4));
    Beta_B_pval(:,nE)=double(mdl_Beta_perE.Coefficients(:,6));
    EffectNames=mdl_Beta_perE.CoefficientNames;
    
    subtbl_probe=tbl_probe(tbl_probe.Chan==nE,:);
    mdl_Alpha_perE= fitlme(subtbl_probe,'alpha~nBlock+Task*State+(1| SubID)');
    Alpha_B_tval(:,nE)=double(mdl_Alpha_perE.Coefficients(:,4));
    Alpha_B_pval(:,nE)=double(mdl_Alpha_perE.Coefficients(:,6));
    EffectNames2=mdl_Alpha_perE.CoefficientNames;
end

%%

figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
varNames={'Task_2','State_2','State_3'};
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=Alpha_B_tval(match_str(EffectNames2,varNames{nplot}),:);
    tempplotpV=Alpha_B_pval(match_str(EffectNames2,varNames{nplot}),:);
    % tempplot(tempplotpV>0.005)=0;
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<fdr(tempplotpV,0.05)),'.','k',10,2},'whitebk','on');
%     topoplot(find(tempplotpV<0.05), layout.chaninfo,'style','blank','electrodes','off','emarker',{'.','k',10,2});
%     topoplot(find(tempplotpV<fdr(tempplotpV,0.05)), layout.chaninfo,'electrodes','on','style','blank','emarker',{'o','k',20,2});
%     
    caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));



figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
varNames={'Task_2','State_2','State_3'};
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=Beta_B_tval(match_str(EffectNames,varNames{nplot}),:);
    tempplotpV=Beta_B_pval(match_str(EffectNames,varNames{nplot}),:);
    % tempplot(tempplotpV>0.005)=0;
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<fdr(tempplotpV,0.05)),'.','k',10,2},'whitebk','on');
%     topoplot(find(tempplotpV<0.05), layout.chaninfo,'style','blank','electrodes','off','emarker',{'.','k',10,2});
%     topoplot(find(tempplotpV<fdr(tempplotpV,0.05)), layout.chaninfo,'electrodes','on','style','blank','emarker',{'o','k',20,2});
%     
    caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));

%%
for nE=1:63
    subtbl_probe=tbl_probe(tbl_probe.Chan==nE,:);
    mdl_ONOFF= fitlme(subtbl_probe,'ONOFF~(alpha_theta+vig_index+W_index+NREM_index+REM_index+alpha+theta+delta+sigma+beta)+ (1| SubID)');
    ONOFF_tval(:,nE)=mdl_ONOFF.Coefficients(:,4);
    ONOFF_pval(:,nE)=mdl_ONOFF.Coefficients(:,6);
    
    mdl_MWMB= fitlme(subtbl_probe,'MWMB~(alpha_theta+vig_index+W_index+NREM_index+REM_index+alpha+theta+delta+sigma+beta)+ (1| SubID)');
    MWMB_tval(:,nE)=mdl_MWMB.Coefficients(:,4);
    MWMB_pval(:,nE)=mdl_MWMB.Coefficients(:,6);
end

%%
load(['EasyCap64_layout'])
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
figure;
for nVar=2:11
    subplot(2,10,nVar-1)
    temp_toplot=double(ONOFF_pval(nVar,:))<0.01;
    topoplot(temp_toplot, layout.chaninfo);
    colormap;
    title(mdl_ONOFF.CoefficientNames{nVar})
    
    subplot(2,10,10+nVar-1)
    temp_toplot=double(ONOFF_tval(nVar,:));
    topoplot(temp_toplot, layout.chaninfo);
    caxis([-4 4])
end

figure;
for nVar=2:11
    subplot(2,10,nVar-1)
    temp_toplot=double(MWMB_pval(nVar,:))<0.01;
    topoplot(temp_toplot, layout.chaninfo);
    colormap;
    title(mdl_MWMB.CoefficientNames{nVar})
    
    subplot(2,10,10+nVar-1)
    temp_toplot=double(MWMB_tval(nVar,:));
    topoplot(temp_toplot, layout.chaninfo);
    caxis([-4 4])
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
