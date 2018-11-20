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
recompute=0;
all_probes_mat=nan(384*60*length(bsl_files),23);
all_probes_mat2=nan(5*60*length(bsl_files),22);
all_probes_mat3=nan(5*60*length(bsl_files),10021);
linecount1=0; linecount2=0; linecount3=0;
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
    left_freq(n)=SubjectInfo.FlickerL;
    right_freq(n)=SubjectInfo.FlickerR;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
    
    %%%% compute RESS
    myFreqs=[6 7.5 12 15 13.5];
    myFreqsN={'6','7','12','15','13'};
    onprobe_logSNR_RESS=[];
    for nT=1:2
        if nT==1
            condName='FA';
            idx_trials=find_trials(D.conditions,'FA');
        elseif nT==2
            condName='DT';
            idx_trials=find_trials(D.conditions,'DT');
        end
        temp_data2=D(1:63,these_times,idx_trials); % D contains the data with channels * time * trials
        for nf=1:length(myFreqs)
            savename=sprintf('wanderIM_RESS_%s_%s_F%s',SubID,condName,myFreqsN{nf});
            if exist([eeg_path filesep savename '.mat'])==0 || recompute==1
                fprintf('... ... compute RESS for %s - %s - %g Hz',SubID,condName,myFreqs(nf))
                tic;
                paramRESS=[];
                paramRESS.peakfreq1=myFreqs(nf);
                paramRESS.peakwidt=0.2;
                paramRESS.neighfreq=1;
                paramRESS.neighwidt=1;
                paramRESS.fft_res=0.1; %0.1 default
                paramRESS.mindist=0.3; %0.5Hz default
                paramRESS.lowerfreq=1; %2Hz default
                [snrR, snrE, faxis2, maps, components]=get_logSNR_RESS(temp_data2,D.fsample,paramRESS);
                save([eeg_path filesep savename],'snrR','snrE','faxis2','maps','components');
                tictoc=toc;
                fprintf('... took %g sec.\n',tictoc)
            else
                fprintf('... ... load RESS for %s - %s - %g Hz\n',SubID,condName,myFreqs(nf))
                load([eeg_path filesep savename])
            end
            param=[];
            param.method='fft'; % fast fourier transform
            param.mindist=1;
            RESS_comp(1,:,:)=components;
            [logSNR_comp, faxis_comp, logpow_comp]=get_logSNR(RESS_comp,D.fsample,param);
            onprobe_logSNR_RESS(n,nf,:,idx_trials)=logSNR_comp;
            onprobe_maps_RESS(n,nf,nT,:)=maps;
            onprobe_components_RESS(n,nf,:,idx_trials)=components;
            onprobe_snrR_RESS(n,nf,nT,:)=snrR;
            onprobe_snrE_RESS(n,nf,nT,:,:)=snrE;
        end
    end
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
            
            temp_behav=[str2num(SubID) nbl npr these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go left_freq(n) right_freq(n) dprime crit];
            
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
            
            FOI=[6 7.5 12 15 13.5];
            FOIlab={'f1','f2','2f1','2f2','IM'};
            temp_EEGres2=[];
            temp_EEGgroup2=[];
            temp_EEGres3=[];
            temp_SNR=[];
            for nf=1:length(FOI)
                [~,fidx]=findclosest(faxis,FOI(nf));
                temp_SNR=onprobe_logSNR_RESS(n,nf,fidx,idx_probe_eeg);
                temp_EEGres2=[temp_EEGres2 ; temp_SNR];
                temp_EEGgroup2=[temp_EEGgroup2 ; [nf]];
                
                temp_Comp=squeeze(onprobe_components_RESS(n,nf,:,idx_probe_eeg));
                temp_EEGres3=[temp_EEGres3 ; temp_Comp'];
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
            all_probes_mat(linecount1+(1:length(temp_EEGres)),:)=[repmat(temp_behav,length(temp_EEGres),1) temp_EEGres temp_EEGgroup];
            linecount1=linecount1+length(temp_EEGres);
            
            all_probes_headers2={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq','dp','crit','SNR','Comp'};
            all_probes_mat2(linecount2+(1:length(temp_EEGres2)),:)=[repmat(temp_behav,length(temp_EEGres2),1) temp_EEGres2 temp_EEGgroup2];
            linecount2=linecount2+length(temp_EEGres2);
            
            all_probes_headers3={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq','dp','crit','Comp'};
            all_probes_mat3(linecount3+(1:length(temp_EEGres2)),:)=[repmat(temp_behav,length(temp_EEGres2),1) temp_EEGgroup2 temp_EEGres3];
            linecount3=linecount3+length(temp_EEGres2);
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

%%
tbl_ress=array2table(all_probes_mat2,'VariableNames',all_probes_headers2);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_ress.SubID=categorical(tbl_ress.SubID);
tbl_ress.Task=categorical(tbl_ress.Task);
tbl_ress.Look=categorical(tbl_ress.Look);
tbl_ress.State=categorical(tbl_ress.State);
tbl_ress.Orig=categorical(tbl_ress.Orig);
tbl_ress.L_freq=categorical(tbl_ress.L_freq);
tbl_ress.R_freq=categorical(tbl_ress.R_freq);
tbl_ress.Comp=categorical(tbl_ress.Comp);

% % writetable(tbl_ress,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh2.txt']);
tbl_ress(tbl_ress.State=="4",:)=[];
tbl_ress.State=removecats(tbl_ress.State);
tbl_ress.State=reordercats(tbl_ress.State,{'1','2','3'});
tbl_ress2=tbl_ress;
tbl_ress2.State=reordercats(tbl_ress2.State,{'2','1','3'});

%%
figure;
plotTitles={'6','7.5','12','15','IM (13.5)'};
for nComp=1:5
    subplot(1,5,nComp); format_fig;
    for ntask=1:2
        for nstate=1:3
            if ntask==1
                tempC=Colors(nstate,:);
            else
                tempC=[1 1 1;Colors(nstate,:)];
            end
            simpleBarPlot(nstate+0.2*(2*ntask-3),tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & tbl_ress.Comp==num2str(nComp)),tempC,0.35,'k');
        end
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title(plotTitles{nComp})
    ylim([2 5])
    xlim([0.2 3.8])
end

%%

figure;
plotTitles={'F','IM (13.5)'};
for nComp=1:2
    subplot(1,2,nComp); format_fig;
    for ntask=1:2
        for nstate=1:3
            if ntask==1
                tempC=Colors(nstate,:);
            else
                tempC=[1 1 1;Colors(nstate,:)];
            end
            if nComp==1
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & (tbl_ress.Comp=='1' |  tbl_ress.Comp=='2'));
                %             elseif nComp==2
                %             toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & (tbl_ress.Comp=='3' |  tbl_ress.Comp=='4'))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(1) & (tbl_ress.Comp=='3' |  tbl_ress.Comp=='4')));
            elseif nComp==2
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & tbl_ress.Comp==num2str(5));
            end
            simpleBarPlot(nstate+0.2*(2*ntask-3),toplot,tempC,0.35,'k',[],2);
        end
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title(plotTitles{nComp})
    ylim([150 450])
    xlim([0.2 3.8])
    %     line(xlim,[100 100],'Color','k','linestyle','--');
    ylabel('RESS SNR')
end

figure;
plotTitles={'F','IM (13.5)'};
for nComp=1:2
    subplot(1,2,nComp); format_fig;
    for ntask=1:2
        for nstate=2:3
            if ntask==1
                tempC=Colors(nstate,:);
            else
                tempC=[1 1 1;Colors(nstate,:)];
            end
            if nComp==1
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & (tbl_ress.Comp=='1' |  tbl_ress.Comp=='2'))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(1)  & (tbl_ress.Comp=='1' |  tbl_ress.Comp=='2')));
                %             elseif nComp==2
                %             toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & (tbl_ress.Comp=='3' |  tbl_ress.Comp=='4'))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(1) & (tbl_ress.Comp=='3' |  tbl_ress.Comp=='4')));
            elseif nComp==2
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & tbl_ress.Comp==num2str(5))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(1) & tbl_ress.Comp==num2str(5)));
            end
            simpleBarPlot(nstate+0.2*(2*ntask-3),toplot,tempC,0.35,'k',{1 100 0.05/4},2);
        end
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title(plotTitles{nComp})
    ylim([70 110])
    xlim([1.2 3.8])
    line(xlim,[100 100],'Color','k','linestyle','--');
    ylabel('RESS SNR (% of ON)')
end

figure;
plotTitles={'F','IM (13.5)'};
for nComp=1:2
    subplot(1,2,nComp); format_fig;
    for nstate=1:2
        for ntask=1:2
            if ntask==1
                tempC=Colors(3,:);
            else
                tempC=[1 1 1;Colors(3,:)];
            end
            if nComp==1
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(3) & (tbl_ress.Comp=='1' |  tbl_ress.Comp=='2'))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate)  & (tbl_ress.Comp=='1' |  tbl_ress.Comp=='2')));
            elseif nComp==2
                toplot=100*tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(3) & tbl_ress.Comp==num2str(5))/mean(tbl_ress.SNR(tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate) & tbl_ress.Comp==num2str(5)));
            end
            simpleBarPlot(nstate+0.2*(2*ntask-3),toplot,tempC,0.35,'k',{1 100 0.05/4},2);
        end
    end
    set(gca,'XTick',1:2,'XTickLabel',{'MB vs ON''MB vs MW'})
    title(plotTitles{nComp})
    ylim([70 110])
    xlim([0.2 2.8])
    line(xlim,[100 100],'Color','k','linestyle','--');
    ylabel('RESS SNR (% of ON)')
end

%%
%% RESS component
% F1, F2, IM for left
% F1, F2, IM for right
% myFreqs=[6 7.5 12 15 13.5];
L_F_FMaps=[];
R_F_FMaps=[];
IM_FMaps=[];
L_F_DMaps=[];
R_F_DMaps=[];
IM_DMaps=[];
L_2F_FMaps=[];
R_2F_FMaps=[];
L_2F_DMaps=[];
R_2F_DMaps=[];
for n=1:9 %length(bsl_files)
    if left_freq(n)==15
        L_F_FMaps=[L_F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,2,1,:),2))))'];
        R_F_FMaps=[R_F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,1,1,:),2))))'];
        
        L_2F_FMaps=[L_2F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,4,1,:),2))))'];
        R_2F_FMaps=[R_2F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,3,1,:),2))))'];
    else
        L_F_FMaps=[L_F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,1,1,:),2))))'];
        R_F_FMaps=[R_F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,2,1,:),2))))'];
        
        L_2F_FMaps=[L_2F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,3,1,:),2))))'];
        R_2F_FMaps=[R_2F_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,4,1,:),2))))'];
    end
    IM_FMaps=[IM_FMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,5,1,:),2))))'];
    
    if left_freq(n)==15
        L_F_DMaps=[L_F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,2,2,:),2))))'];
        R_F_DMaps=[R_F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,1,2,:),2))))'];
        
        L_2F_DMaps=[L_2F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,4,2,:),2))))'];
        R_2F_DMaps=[R_2F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,3,2,:),2))))'];
    else
        L_F_DMaps=[L_F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,1,2,:),2))))'];
        R_F_DMaps=[R_F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,2,2,:),2))))'];
        
        L_2F_DMaps=[L_2F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,3,2,:),2))))'];
        R_2F_DMaps=[R_2F_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,4,2,:),2))))'];
    end
    IM_DMaps=[IM_DMaps ; zscore(real(squeeze(mean(onprobe_maps_RESS(n,5,2,:),2))))'];
end

figure;
subplot(2,5,1); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Right Fund')

subplot(2,5,2); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Left Fund')

subplot(2,5,3); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_2F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Right 2F')

subplot(2,5,4); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_2F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Left 2F')


subplot(2,5,5); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(IM_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - IM: F2-F1')

subplot(2,5,5+1); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Right Fund')

subplot(2,5,5+2); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Left Fund')

subplot(2,5,5+3); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_2F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Right 2F')

subplot(2,5,5+4); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_2F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Left 2F')

subplot(2,5,5+5); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(IM_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - IM: F2-F1')

%%
tbl_ress=array2table(all_probes_mat2,'VariableNames',all_probes_headers2);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_ress.SubID=categorical(tbl_ress.SubID);
tbl_ress.Task=categorical(tbl_ress.Task);
tbl_ress.Look=categorical(tbl_ress.Look);
tbl_ress.State=categorical(tbl_ress.State);
tbl_ress.Orig=categorical(tbl_ress.Orig);
tbl_ress.L_freq=categorical(tbl_ress.L_freq);
tbl_ress.R_freq=categorical(tbl_ress.R_freq);
tbl_ress.Comp=categorical(tbl_ress.Comp);

% % writetable(tbl_ress,[behav_path filesep 'WanderIM_ProbeResults_EEG_allCh2.txt']);
tbl_ress(tbl_ress.State=="4",:)=[];
tbl_ress.State=removecats(tbl_ress.State);
tbl_ress.State=reordercats(tbl_ress.State,{'1','2','3'});

figure; format_fig;
for nComp=1:5
        subplot(1,5,nComp); format_fig;
    for ntask=1:2
        for nstate=1:3
%             if nstate~=2
                tempplot=(tbl_ress.SNR(tbl_ress.Comp==num2str(nComp) & tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate)));
%             else
%                 tempplot=(tbl_ress.SNR(tbl_ress.Comp==num2str(nComp) & tbl_ress.Task==num2str(ntask) & tbl_ress.State==num2str(nstate)  & tbl_ress.Orig==num2str(nstate)));
%             end
            if ntask==1
                simpleBarPlot(nstate+0.2*(ntask*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
            elseif ntask==2
                simpleBarPlot(nstate+0.2*(ntask*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
            end
        end
    end
    xlim([0.2 3.8])
    
end
