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

% % writetable(tbl_probe,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh2.txt']);
tbl_probe(tbl_probe.State=="4",:)=[];
tbl_probe.State=removecats(tbl_probe.State);
tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});
tbl_probe2=tbl_probe;
tbl_probe2.State=reordercats(tbl_probe2.State,{'2','1','3'});


%% model on performance - by electrode
clear pval_* beta_*
fprintf('%2.0f/63\n',0)
for n=1:63
    fprintf('\b\b\b\b\b\b%2.0f/63\n',n)
    temp_p=tbl_probe(tbl_probe.Chan == num2str(n) & ismember(tbl_probe.Freq,{'1','2'}),:);
    temp_mdl= fitlme(temp_p,'SNR~nBlock+Task*State + (1| SubID)');
    pval_F(n,:)=double(temp_mdl.Coefficients(:,6));
    beta_F(n,:)=double(temp_mdl.Coefficients(:,2));
    var_F=temp_mdl.CoefficientNames;
    
    temp_p=tbl_probe(tbl_probe.Chan == num2str(n) & ismember(tbl_probe.Freq,{'3','4'}),:);
    temp_mdl= fitlme(temp_p,'SNR~nBlock+Task*State + (1| SubID)');
    pval_2F(n,:)=double(temp_mdl.Coefficients(:,6));
    beta_2F(n,:)=double(temp_mdl.Coefficients(:,2));
    var_2F=temp_mdl.CoefficientNames;
    
    temp_p=tbl_probe(tbl_probe.Chan == num2str(n) & ismember(tbl_probe.Freq,{'5'}),:);
    temp_mdl= fitlme(temp_p,'SNR~nBlock+Task*State + (1| SubID)');
    pval_IM(n,:)=double(temp_mdl.Coefficients(:,6));
    beta_IM(n,:)=double(temp_mdl.Coefficients(:,2));
    var_IM=temp_mdl.CoefficientNames;
end

%%
load(['EasyCap64_layout'])
figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
varNames={'Task_2','Task_2:State_2','Task_2:State_3'};
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=beta_IM(:,match_str(var_F,varNames{nplot}));
    tempplotpV=pval_IM(:,match_str(var_F,varNames{nplot}));
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<0.05),'.','k',10,2},'whitebk','on');
%     caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));

