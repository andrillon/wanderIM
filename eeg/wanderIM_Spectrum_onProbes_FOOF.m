%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(genpath(fooof_path))

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
myElecs=1:63;
foof_bg_probes=nan(length(bsl_files)*60*length(myElecs),6);
info_probes=nan(length(bsl_files)*60*length(myElecs),15);
psdall_probes=[];
vigindexes_probes=nan(length(bsl_files)*60*length(myElecs),13);
psd_probes=nan(length(bsl_files),length(myElecs),513,60);
ncp=0;
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
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    mastref_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials
    
    %%%%%%
    f_range = [2, 40];
    settings = struct();  % Use defaults
    fprintf('%2.0f-%2.0f\n',0,0)
    for nPr=1:size(mastref_data,3)
        for nE=myElecs
            ncp=ncp+1;
            fprintf('\b\b\b\b\b\b%2.0f-%2.0f\n',nPr,nE)
            % Calculate a power spectrum with Welch's method
            [psd_probes(n,nE,:,nPr), freqs] = pwelch(squeeze(mastref_data(nE,:,nPr)), 1000, [], [], D.fsample);
            
            % Run FOOOF
            fooof_results = fooof(freqs', squeeze(psd_probes(n,nE,:,nPr))', f_range, settings);
            foof_bg_probes(ncp,:)=[n nPr nE fooof_results.background_params fooof_results.error];
            info_probes(ncp,:)=[n nPr nE probe_res(nPr,[1 4:6 31:38])];
            %             psdall_probes=[psdall_probes ; squeeze(psd_probes(n,nE,:,nPr))'];
            
            [alpha_theta,vig_index,W_index,NREM_index,REM_index, powbyband]=get_sleep_vig_indexes_ps(squeeze(psd_probes(n,nE,:,nPr))',freqs',[]);
            vigindexes_probes(ncp,:)=[n nPr nE alpha_theta vig_index,W_index,NREM_index,REM_index powbyband.delta powbyband.theta powbyband.alpha powbyband.spindle powbyband.beta];
        end
    end
end

%%
tbl=array2table([foof_bg_probes vigindexes_probes(:,4:end) info_probes(:,[6:15])],'VariableNames',...
    {'SubID','nPr','nE','offset','slope','error',...
    'a_t','vigidx','Widx','NREMidx','REidx','dt','th','al','sp','be',...
    'task','ntrial','look','state','ori','awa','int','eng','perf','alert'});
tbl(tbl.state==4,:)=[];
tbl.state2=tbl.state;
tbl.nE=categorical(tbl.nE);
tbl.task=categorical(tbl.task);
tbl.state=categorical(tbl.state);
tbl.SubID=categorical(tbl.SubID);

tbl2=tbl(tbl.nE=="2" | tbl.nE=="17",:);
tbl2.nE=removecats(tbl2.nE);
%% plot variable of interests
voi={'offset','slope','a_t'}; %,'vigidx','Widx','NREMidx','REidx','dt','th','al','sp','be'};
myElecs2=[30];
for nplot=1:length(voi)
    figure;
    for nE=1:length(myElecs2)
        subplot(1,length(myElecs2),nE);
        format_fig;
        eval(sprintf('temp=tbl.%s;',voi{nplot}))
        for nstate=1:3
            tempplot=temp(tbl.state2==nstate & tbl.nE==num2str(myElecs(nE)));
            simpleBarPlot(nstate,tempplot,Colors(nstate,:),0.9,'k');
        end
        set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'});
        title([D.chanlabels{myElecs2(nE)}]);
        ylabel(voi{nplot})
        xlim([0.2 3.8])
    end
end

figure;
stateNames={'ON','MW','MB'};
for nstate=1:3
    for nstate2=nstate:3
        if nstate==nstate2
            temp_topo=[]; tempplotpV=[];
            for nE=1:63
                temp_topo(nE)=mean(tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate));
                [~,tempplotpV(nE)]=ttest(tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate),1);
            end
            
        else
            temp_topo=[]; tempplotpV=[];
            for nE=1:63
                temp_topo(nE)=mean(tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate2))-mean(tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate));
                [~,tempplotpV(nE)]=ttest2(tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate2),tbl.a_t(tbl.nE==num2str(nE) & tbl.state2==nstate));
            end
        end
        subplot(3,3,3*(nstate-1)+nstate2)
        addpath(genpath(path_eeglab));
        stat_thr=fdr(tempplotpV,0.05);
        topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{[],'.','w',24,2});
        rmpath(genpath(path_eeglab));
        if nstate==nstate2
            caxis([0 1.8])
            title(stateNames{nstate})
        else
            caxis([-1 1]*0.5)
            title({stateNames{nstate2},'vs',stateNames{nstate}})

        end
        colormap('parula'); colorbar;
        format_fig
    end
end
