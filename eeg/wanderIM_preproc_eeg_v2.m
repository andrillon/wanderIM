%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(genpath(sleepTools_path));
addpath(genpath(spm12_path));
preproc_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
files=dir([preproc_path filesep 'EEG_S*.mat']);

%% General parameters
trial_window=[-0.2 1.3];
hb_window=[-0.1 0.4];
probe_window=[-15 15];

%% Loop on files
redo.hp=0;
redo.notch=1;

for n=1:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    
    filebehav=dir([behav_path filesep 'wanderIM_behavres_s' subID '*.mat']);
    Behav=load([behav_path filesep filebehav.name]);
    D=spm_eeg_load([preproc_path filesep filename]);
    
   
    %%% High-pass filter
    if exist([preproc_path filesep 'f' D.fname])==0 || redo.hp==1
        % update channel type for ECG
        D = chantype(D, match_str(D.chanlabels,'ECG'), 'ECG');
        D.save;
        
        
        type = 'butterworth';
        order = 5;
        dirfilt = 'twopass';
        S = [];
        S.D = D;
        S.band = 'high';
        S.type = type;
        S.order = order;
        S.dir = dirfilt;
        S.freq = 0.1;
        S.save=1;
        D = spm_eeg_filter(S);
    else
        D=spm_eeg_load([preproc_path filesep 'f' D.fname]);
    end
    
    %%% NOTCH FILTERS
    if exist([preproc_path filesep 'n' D.fname])==0 || redo.notch==1
        type = 'butterworth';
        order = 4;
        dirfilt = 'twopass';
        S = [];
        S.prefix='n';
        S.D = D;
        S.band = 'stop';
        S.type = type;
        S.order = order;
        S.dir = dirfilt;
        S.freq = [45 55];
        S.save=1;
        D = spm_eeg_filter(S);
    else
        D=spm_eeg_load([preproc_path filesep 'n' D.fname]);
    end
    
    %%% Find Triggers
    clear din_*
    din_idx=match_str(D.chanlabels,'D1');
    din_chan=D(din_idx,:,1);
    % take away the first 10% to normalise
    din_chan=(din_chan-min(din_chan(round(0.1*length(din_chan)):end)))/max(din_chan(round(0.1*length(din_chan)):end)-min(din_chan(round(0.1*length(din_chan)):end)));
    % find the din
    din_thr=0.2;
    din_start=find(diff(din_chan>din_thr)==1)-2;
    din_end=find(diff(din_chan>din_thr)==-1)+2;
    % eliminate incomplete din (no end or no start)
    while din_end(1)<din_start(1)
        din_end(1)=[];
    end
    if length(din_start)>length(din_end)
        din_start=din_start(1:length(din_end));
    end
    if length(din_end)>length(din_end)
        din_end=din_end(1:length(din_start));
    end
    din_dur=(din_end-din_start)/D.fsample;
    din_dist=din_start(1:end-1)-din_start(2:end);
    % eliminate din too long or too short
    discard_din=din_dur>10 & din_dur<0.005;
    din_start(discard_din)=[];
    din_end(discard_din)=[];
    din_dur(discard_din)=[];
    din_dist(discard_din)=[];
    fprintf('... ... %g triggers found\n',length(din_start));
    
    %%% Correct trial timing with triggers
    % {'B'} Block {'C'} Probe End {'E'} End {'K'} End Block {'P'} Probe {'Q'} Question
    % {'S'} Start {'T'} Trial
    myevents=D.events;
    fprintf('... ... %g blocks found (end %g)\n',length(match_str({myevents.type},'B')),length(match_str({myevents.type},'K')));
    fprintf('... ... %g probes found (end %g)\n',length(match_str({myevents.type},'P')),length(match_str({myevents.type},'C')));
    fprintf('... ... %g trials found \n',length(match_str({myevents.type},'T')));
    
    if length(match_str({myevents.type},'B'))~=14 % 8 baseline blocks + 6 test blocks
        warning('Incorrect number of blocks!')
        continue;
    end
    if length(match_str({myevents.type},'P'))~=60 % 8 baseline blocks + 6 test blocks
        warning('Incorrect number of probes!')
        continue;
    end
    
    % Identify start/end blocks
    idBlock_start=match_str({myevents.type},'B');
    idBlock_end=match_str({myevents.type},'K');
    for nB=1:8
        thist=D.indsample(myevents(idBlock_start(nB)).time);
        [thisd, thisidx]=findclosest(din_start,thist);
        lag_block_baseline(nB)=thist-thisd;
        start_block_baseline(nB)=thisd;
    end
    for nB=9:14
        thist=D.indsample(myevents(idBlock_start(nB)).time);
        [thisd, thisidx]=findclosest(din_start,thist);
        lag_block_test(nB-8)=thist-thisd;
        start_block_test(nB-8)=thisd;
        
        thist=D.indsample(myevents(idBlock_end(nB)).time);
        [thisd, thisidx]=findclosest(din_start,thist);
        lag_block_test2(nB-8)=thist-thisd;
        end_block_test(nB-8)=thisd;
    end
    fprintf('... ... max lag for blocks: start %g / %g (in samples)\n',max(abs(lag_block_test)),max(abs(lag_block_test2)))
    
    % Identify start trials
    idTrial_start=match_str({myevents.type},'T');
    Trial_start=[myevents(idTrial_start).time];
    lag_trial=nan(1,length(idTrial_start));
    start_trial=nan(1,length(idTrial_start));
    for nT=1:length(idTrial_start)
        thist=D.indsample(myevents(idTrial_start(nT)).time);
        [thisd, thisidx]=findclosest(din_start,thist);
        lag_trial(nT)=thist-thisd;
        start_trial(nT)=thisd;
        durdin_trial(nT)=din_dur(thisidx);
    end
    training_trials=start_trial<start_block_test(1);
    % eliminate training trials
    clean_lag_trial=lag_trial(~training_trials);
    clean_durdin_trial=durdin_trial(~training_trials);
    clean_start_trial=start_trial(~training_trials);
    
    % eliminate end block trials
    clean_lag_trial(clean_durdin_trial>2)=[];
    clean_start_trial(clean_durdin_trial>2)=[];
    
    if length(clean_start_trial)~=size(Behav.test_res,1)
        warning('Problem with din and number of trials!')
        continue;
    end
    fprintf('... ... max lag for trials: start %g (in samples)\n',max(abs(clean_lag_trial)))
    
    % Identify start probes
    idProbes_start=match_str({myevents.type},'P');
    lag_probe=nan(1,length(idProbes_start));
    start_probe=nan(1,length(idProbes_start));
    for nT=1:length(idProbes_start)
        thist=D.indsample(myevents(idProbes_start(nT)).time);
        [thisd, thisidx]=findclosest(din_start,thist);
        lag_probe(nT)=thist-thisd;
        start_probe(nT)=thisd;
    end
    fprintf('... ... max lag for probes: start %g (in samples)\n',max(abs(lag_probe)))
    
    
    %%%%% Epoch baseline
    if exist([preproc_path filesep 'basel_' D.fname])==0 || redo==1
        pretrig  =  -5 * D.fsample;
        posttrig =  30 * D.fsample;
        trl=[];
        trllabels=[];
        name_baseline_blocks={'Faces_IM','Faces_NoIM','Squares_IM','Squares_NoIM','Faces_IM','Faces_NoIM','Squares_IM','Squares_NoIM'};
        for ntr=1:length(start_block_baseline)
            trlbegin = start_block_baseline(ntr) + pretrig;
            trlend   = start_block_baseline(ntr) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend offset];
            trl      = [trl; newtrl];
            trllabels{ntr}= name_baseline_blocks(ntr);
        end
        S=[];
        S.prefix = 'basel_';
        S.D=D;
        S.bc=1;
        S.trl=trl;
        S.conditionlabels=trllabels;
        S.save=1;
        D_baseline=spm_eeg_epochs(S);
    end
    
    %%%%% Epoch probes
    if exist([preproc_path filesep 'probe_' D.fname])==0 || redo==1
        pretrig  =  probe_window(1) * D.fsample;
        posttrig =  probe_window(2) * D.fsample;
        trl=[];
        trllabels=[];
        task_type={'FA','DT'};
        mindstate_type={'ON','MW','MB','XX'};
        for ntr=1:length(start_probe)
            trlbegin = start_probe(ntr) + pretrig;
            trlend   = start_probe(ntr) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend offset];
            trl      = [trl; newtrl];
            
            this_probe=Behav.probe_res(ntr,:);
            trllabels{ntr}= sprintf('B%g_P%g_%s_%s',this_probe(4),this_probe(1),task_type{this_probe(5)},mindstate_type{this_probe(32)});
        end
        S=[];
        S.prefix = 'probe_';
        S.D=D;
        S.bc=1;
        S.trl=trl;
        S.conditionlabels=trllabels;
        S.save=1;
        D_probes=spm_eeg_epochs(S);
    end
    
    %%%%% Epoch trials
    if exist([preproc_path filesep 'trial_' D.fname])==0 || redo==1
        pretrig  =  trial_window(1) * D.fsample;
        posttrig =  trial_window(2) * D.fsample;
        trllabels=[];
        trl=[];
        task_type={'FA','DT'};
        fprintf('... ... %3.0f\n',0)
        %     erp_try=nan(5,1501,length(clean_start_trial));
        for ntr=1:length(clean_start_trial)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b... ... %3.0f\n',round(ntr/length(clean_start_trial)*100))
            trlbegin = clean_start_trial(ntr) + pretrig;
            trlend   = clean_start_trial(ntr) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend offset];
            trl      = [trl; newtrl];
            
            this_trial=Behav.test_res(ntr,:);
            if ~isnan(this_trial(11))
                trllabels{ntr}= sprintf('B%g_T%g_%s_NG%g',this_trial(1),this_trial(4),task_type{this_trial(2)},this_trial(11));
            elseif ~isnan(this_trial(12))
                trllabels{ntr}= sprintf('B%g_T%g_%s_GO%g',this_trial(1),this_trial(4),task_type{this_trial(2)},this_trial(12));
            end
            
            %         erp_try(:,:,ntr)=D(match_str(D.chanlabels,{'Oz','Pz','Cz','TP9','TP8'}),trlbegin:trlend,1);
        end
        %     go1=find(~cellfun(@isempty,regexpi(trllabels,'GO1')));
        %     nogo1=find(~cellfun(@isempty,regexpi(trllabels,'NG1')));
        %     go0=find(~cellfun(@isempty,regexpi(trllabels,'GO0')));
        %     nogo0=find(~cellfun(@isempty,regexpi(trllabels,'NG0')));
        %     erp_time=(pretrig:posttrig)/D.fsample;
        %     erp_try=erp_try-repmat(mean(erp_try(:,erp_time>-0.1 & erp_time<0,:),2),[1 size(erp_try,2) 1]);
        S=[];
        S.prefix = 'trial_';
        S.D=D;
        S.bc=1;
        S.trl=trl;
        S.conditionlabels=trllabels;
        S.save=1;
        D_trials=spm_eeg_epochs(S);
    end
    
    %%%%% Epoch HeartBeat
    if exist([preproc_path filesep 'hb_' D.fname])==0 || redo==1
        hb_times=detect_heartbeat(D,2,0);
        
        pretrig  =  hb_window(1) * D.fsample;
        posttrig =  hb_window(2) * D.fsample;
        trllabels=[];
        trl=[];
        fprintf('... ... %3.0f\n',0)
        for ntr=1:length(hb_times)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b... ... %3.0f\n',round(ntr/length(hb_times)*100))
            trlbegin = hb_times(ntr) + pretrig;
            trlend   = hb_times(ntr) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend offset];
            trl      = [trl; newtrl];
            trllabels{ntr}='HB';
        end
        S=[];
        S.prefix = 'hb_';
        S.D=D;
        S.bc=1;
        S.trl=trl;
        S.conditionlabels=trllabels;
        S.save=1;
        D_hb=spm_eeg_epochs(S);
    end
end

%%
%
% figure
% plot(erp_time,mean(erp_try(2,:,nogo1),3),'r')
% hold on
% plot(erp_time,mean(erp_try(2,:,go1),3),'b')
% plot(erp_time,mean(erp_try(2,:,go0),3),'b--')
% plot(erp_time,mean(erp_try(2,:,nogo0),3),'r--')
