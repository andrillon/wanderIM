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
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=9; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)

%% loop across trials for baseline blocks
all_Waves_byProbes=[];
all_Waves_byProbes2=[];
all_len=[];
all_Waves_byProbes3=[];
nc=0;
for n=1:length(bsl_files)
    %     if n==5
    %         continue;
    %     end
    nc=nc+1;
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    fprintf('... processing subject %s\n',SubID)
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    load([eeg_path filesep 'triggers_S' SubID]);
    if size(test_res,1)~=length(clean_start_trial)
        warning('mismatch with trigers')
    end
    load([eeg_path filesep 'wanderIM_cont_twa2_' SubID '_selected'])
    if AmpCriterionIdx==9
        slow_Waves(:,AmpCriterionIdx)=-slow_Waves(:,AmpCriterionIdx);
    end
    slow_Waves(slow_Waves(:,4)<75/2,:)=[];
    for nE=1:63
        thisE_Waves=slow_Waves(slow_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        
        num_Waves(nc,nE)=size(thisE_Waves,1);
        upSlope_Waves(nc,nE)=mean(thisE_Waves(:,12));
        downSlope_Waves(nc,nE)=mean(thisE_Waves(:,11));
        
        probe_res=[probe_res (probe_res(:,4)-1)*10+probe_res(:,1)];
        for npr=1:60
            % Gather behaviour
            these_probes=probe_res(probe_res(:,end)==npr,:);
            these_trials=test_res(test_res(:,1)==these_probes(4),:);
            this_pr_tridx=these_probes(6);
            this_pr_tridx_all=find(test_res(:,1)==these_probes(1) & test_res(:,4)==this_pr_tridx);
            if these_probes(1)==1
                last_pr_tridx=0;
            else
                last_pr_tridx=probe_res(probe_res(:,end)==npr-1,6);
            end
            
            thesegotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)~=3);
            thesenogotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)==3);
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            
            tcorr_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),12));%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),11));%/corr_nogo(n,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(~isnan(temp_testres(:,12)),12)));
            num_nogo=sum(~isnan(temp_testres(~isnan(temp_testres(:,11)),11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(~isnan(temp_testres(:,12)),12)),10)-temp_testres(~isnan(temp_testres(~isnan(temp_testres(:,12)),12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(~isnan(temp_testres(:,11)),11)),10)-temp_testres(~isnan(temp_testres(~isnan(temp_testres(:,11)),11)),8));
            
            tcorr_goall=(temp_testres(~isnan(temp_testres(:,12)),12));
            tcorr_nogoall=(temp_testres(~isnan(temp_testres(:,11)),11));
            [dprime, crit]=calc_dprime2(tcorr_goall,tcorr_nogoall==0);
            
            last_pr_timedx_all=find(test_res(:,1)==these_probes(1) & test_res(:,4)==min(temp_testres(:,4)));
            last_pr_timedx=clean_start_trial(last_pr_timedx_all);
            this_pr_timeidx=clean_start_trial(this_pr_tridx_all);
            thisE_Waves2=thisE_Waves(thisE_Waves(:,5)>last_pr_timedx & thisE_Waves(:,5)<this_pr_timeidx,:);
            if ~isempty(thisE_Waves2)
                temp_len=abs((thisE_Waves2(:,5)-thisE_Waves2(:,7)))./500;
                temp_freq=1./temp_len;
                temp_abs=1./temp_len;
                temp_p2p=thisE_Waves2(:,AmpCriterionIdx);
                
                tp_num_Start=thisE_Waves2(:,5);
                tp_num_Waves=size(thisE_Waves2,1);
                tp_Amp_Waves=nanmean(temp_p2p);
                tp_P2P_Waves=nanmean(thisE_Waves2(:,4));
                tp_Glob_Waves=nanmean(thisE_Waves2(:,end));
                
                tp_upSlope_Waves=nanmean(thisE_Waves2(:,12));
                tp_downSlope_Waves=nanmean(thisE_Waves2(:,11));
                
                
                all_Waves_byProbes=[all_Waves_byProbes ; [n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves tp_P2P_Waves tp_Amp_Waves rt_go tcorr_go tcorr_nogo]];
                all_Waves_byProbes2=[all_Waves_byProbes2 ; [n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves tp_P2P_Waves tp_Amp_Waves]];
                all_Waves_byProbes3=[all_Waves_byProbes3 ; [repmat([n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32)],length(tp_num_Start),1) tp_num_Start]];
            else
                all_Waves_byProbes=[all_Waves_byProbes ; [n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) 0 NaN NaN NaN NaN rt_go tcorr_go tcorr_nogo]];
                all_Waves_byProbes2=[all_Waves_byProbes2 ; [n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32) 0 NaN NaN NaN NaN]];
%                 all_Waves_byProbes3=[all_Waves_byProbes3 ; [repmat([n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32)],length(tp_num_Start),1) tp_num_Start]];
            end
        end
        
%         if nE==24 %Cz=24
%             
%             thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
%             temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
%             temp_freq=1./temp_len;
%             temp_abs=1./temp_len;
%             temp_p2p=thisE_Waves(:,4);
%             onset_Waves=thisE_Waves(temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),5);
%             probes_Waves=thisE_Waves(temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),2);
%             
%             temp_ERP_Waves=[];
%             for nW=1:length(onset_Waves)
%                 if onset_Waves(nW)-1*D.fsample<1 || onset_Waves(nW)+2*D.fsample>size(temp_data,2)-1
%                     temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
%                     continue;
%                 end
%                 temp=temp_data(nE,onset_Waves(nW)+(-1*D.fsample:2*D.fsample),probes_Waves(nW));
%                 if max(abs(temp))>200
%                     temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
%                     continue;
%                 end
%                 temp_ERP_Waves(nW,:)=temp;
%                 temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
%             end
%             numERP_Waves(nc)=size(temp_ERP_Waves,1);
%             if numERP_Waves(nc)<50
%                 ERP_Waves(nc,:)=nan;
%             else
%                 ERP_Waves(nc,:)=nanmean(temp_ERP_Waves);
%             end
%         end
    end
end

%%
Fs=500;
figure;
format_fig;
temp_topo=nanmean(num_Waves(setdiff(1:size(num_Waves,1),5),:))/20/60;
addpath(genpath(path_eeglab));
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
caxis([0 1])
rmpath(genpath(path_eeglab));
cb=colorbar; set(cb,'FontSize',24)

% figure; format_fig
% simpleTplot(-1:1/Fs:2,squeeze((ERP_Waves(:,:))),0,'k',[0],'-',0.5,1,2,[],2);
% xlim([-0.3 0.6]);
% xlabel('Time from wave onset (s)')
% ylabel('Amplitude (\muV)')
% title('20% higher amplitude')

%%
%[n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves]
tbl=array2table(all_Waves_byProbes,'VariableNames',{'SubID','nProbe','nBlock','Task','nTrial','State','Vig','Chan','State2','nWave','UpSlope','DownSlope','P2PAmp','NegAmp','RT','GO','NOGO'});
tbl(tbl.State==4,:)=[];
% tbl.nBlock=categorical(tbl.nBlock);
tbl.Task=categorical(tbl.Task);
tbl.State=categorical(tbl.State);
% tbl.Chan=categorical(tbl.Chan);
tbl.pWave=tbl.nWave;
tbl.pWave(tbl.nWave>0)=1;
tbl.pWave(tbl.nWave<=0)=0;
myS=unique(tbl.SubID);
for nE=1:63
    for nS=1:length(myS)
        this_tbl=tbl(tbl.Chan==nE & tbl.SubID==myS(nS),:);
        [r_NOGO(nS,nE) pV_NOGO(nS,nE)]=corr(this_tbl.nWave,1-this_tbl.NOGO,'type','pearson','rows','pairwise');
        [r_MISS(nS,nE) pV_MISS(nE)]=corr(this_tbl.nWave,1-this_tbl.GO,'type','pearson','rows','pairwise');
        
%         [b_NOGO(nS,nE,:)]=glmfit(this_tbl.nWave,1-this_tbl.NOGO);
%         [b_MISS(nS,nE,:)]=glmfit(this_tbl.nWave,1-this_tbl.GO);
        
        [r_RT(nS,nE) pV_RT(nS,nE)]=corr(this_tbl.nWave(this_tbl.RT>0.2),this_tbl.RT(this_tbl.RT>0.2),'type','pearson','rows','pairwise');
        [r_Vig(nS,nE) pV_Vig(nS,nE)]=corr(this_tbl.nWave,this_tbl.Vig,'type','pearson','rows','pairwise');

    end
    
    this_tbl=tbl(tbl.Chan==nE,:);
    %     mdl_3= fitlme(this_tbl,'Vig~Task');
    %     mdl_0= fitlme(this_tbl,'Vig~1+(1|SubID)');
    %
    %
    %     temp=tbl.RT(tbl.Chan==nE);
    %     temp(temp<0.2)=nan;
end
% filename = '/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results3.txt';
% res=import_localSleepRes_fromR(filename);

%%
figure;
subplot(1,3,1)
addpath(genpath(path_eeglab));
temp_topo=nanmean(r_MISS,1);
[~,tempplotpV]=ttest(r_MISS,0);
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); %caxis([-1 1]*0.1)
format_fig;
title('Misses')

subplot(1,3,2)
addpath(genpath(path_eeglab));
temp_topo=nanmean(r_NOGO,1);
[~,tempplotpV]=ttest(r_NOGO,0);
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); caxis([-1 1]*0.1)
format_fig;
title('False Alarms')

subplot(1,3,3)
addpath(genpath(path_eeglab));
temp_topo=nanmean(r_RT,1);
[~,tempplotpV]=ttest(r_RT,0);
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); caxis([-0.2 0.2])
format_fig;title('RT')

thisChannel=17;
figure; format_fig;
simpleCorPlotsetbin(tbl.nWave(tbl.Chan==thisChannel),tbl.RT(tbl.Chan==thisChannel),0:2:8,{'o','k',[1 1 1]*0.7,72,3});
xlabel('Number of Waves')
ylabel('RT')
set(gca,'XTick',0:2:8,'XTickLabel',{'0','2','4','6','8+'})
xlim([-0.5 8.5])

%%
figure;
hb=[];
hb(1)=simpleBarPlot(1-0.2,tbl.RT(tbl.Vig<=2 & tbl.Chan==17),'b',0.35,'k',[],4);
simpleBarPlot(2-0.2,tbl.RT(tbl.nWave==0 & tbl.Chan==17),[1 1 1;0 0 1],0.35,'k',[],4);

hb(2)=simpleBarPlot(1+0.2,tbl.RT(tbl.Vig>2 & tbl.Chan==17),'r',0.35,'k',[],4);
simpleBarPlot(2+0.2,tbl.RT(tbl.nWave>0 & tbl.Chan==17),[1 1 1;1 0 0],0.35,'k',[],4);
xlim([0.2 2.8])
ylim([0.5 0.6])
format_fig
set(gca,'XTick',1:2,'XTickLabel',{'subj','obj'})
xlabel('Tiredness')
ylabel('RT (s)')
legend(hb,{'alert','tired'})

%%
filename = '/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results3.txt';
res=import_localSleepRes_fromR(filename);

figure;
% subplot(1,3,1)
addpath(genpath(path_eeglab));
temp_topo=res.Chi;
tempplotpV=res.pVal;
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));


% figure;
% % subplot(1,3,1)
% addpath(genpath(path_eeglab));
% temp_topo=res.MBvsON; temp_topo(temp_topo>100)=NaN;
% tempplotpV=res.pMBvsON;
% stat_thr=0.05; %fdr(tempplotpV,0.05);
% topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
% rmpath(genpath(path_eeglab));
% caxis([-1 1]*3)
%%
EOI=[2 35 62 63]; %find(tempplotpV<stat_thr);
EOI2=[13 46 47 48]; %[16 17 18 45:49]';

figure;
for ntask=1:2
    subplot(2,2,ntask);
    format_fig;
    for nstate=1:3
        temp=[];
        for nE=EOI'
            temp=[temp tbl.P2P(tbl.Task==num2str(ntask) & tbl.State==num2str(nstate) & ismember(tbl.Chan,nE))];
        end
        simpleBarPlot(nstate,mean(temp,2),Colors(nstate,:),0.9,'k',[],2);
    end
    %     ylim([20 35])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    xlim([0.2 3.8])
    ylabel('ampl. theta-waves')
end
for ntask=1:2
    subplot(2,2,2+ntask);
    format_fig;
    for nstate=1:3
        temp=[];
        for nE=EOI2'
            temp=[temp tbl.P2P(tbl.Task==num2str(ntask) & tbl.State==num2str(nstate) & ismember(tbl.Chan,nE))];
        end
        simpleBarPlot(nstate,mean(temp,2),Colors(nstate,:),0.9,'k',[],2);
    end
    %     ylim([20 35])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    xlim([0.2 3.8])
    ylabel('ampl. theta-waves')
end

figure;
for ntask=1:2
    %     subplot(1,2,ntask);
    format_fig;
    for nstate=1:3
        temp=[];
        for nE=EOI'
            temp=[temp tbl.nWave(tbl.Task==num2str(ntask) & tbl.State==num2str(nstate) & ismember(tbl.Chan,nE))];
        end
        if ntask==1
            thisC=Colors(nstate,:);
        else
            thisC=[1 1 1;Colors(nstate,:)];
        end
        simpleBarPlot(nstate+0.2*(2*ntask-3),mean(temp,2),thisC,0.35,'k',[],4);
    end
    ylim([4 8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    xlim([0.2 3.8])
    ylabel('Numb. theta-waves')
end

%%
for ntask=1:2
    subplot(2,2,2+ntask);
    format_fig;
    for nstate=1:3
        temp=[];
        for nE=EOI2'
            temp=[temp tbl.nWave(tbl.Task==num2str(ntask) & tbl.State==num2str(nstate) & ismember(tbl.Chan,nE))];
        end
        simpleBarPlot(nstate,mean(temp,2),Colors(nstate,:),0.9,'k',[],2);
    end
    %     ylim([20 35])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    xlim([0.2 3.8])
    ylabel('Numb. theta-waves')
end
subplot(1,3,2)
addpath(genpath(path_eeglab));
temp_topo=res.MWvsON;
tempplotpV=res.pMWvsON;
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));


subplot(1,3,3)
addpath(genpath(path_eeglab));
temp_topo=res.MBvsON;
tempplotpV=res.pMBvsON;
stat_thr=fdr(tempplotpV,0.05);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));
%%
% 8: resp Q4 Awareness (1 (fully) to 4 (not at all))
% 9: resp Q5 Intention (1 (fully) to 4 (not at all))
% 10: resp Q6 Engagement (1 (not) to 4 (very))
% 11: resp Q7 Performance (1 (bad) to 4 (good))
% 12: resp Q8 Alterness (1 (alert) to 4 (sleepy))
corrScale=[-1 -1 1 1 1];
Nameplot={'Awareness','Intention','Engagement','Performance','Sleepiness'};
figure;
for np=1:5
    subplot(2,5,np)
    r=[]; pV=[];
    for nE=1:63
        if corrScale(np)==1
            [r(nE) pV(nE)]=corr(double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,15)),double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,7+np)),'rows','pairwise','type','spearman');
        else
            [r(nE) pV(nE)]=corr(double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,15)),5-double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,7+np)),'rows','pairwise','type','spearman');
        end
    end
    addpath(genpath(path_eeglab));
    temp_topo=r;
    tempplotpV=pV;
    stat_thr=fdr(tempplotpV,0.05);
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','k',24,2});
    rmpath(genpath(path_eeglab));
    caxis([-1 1]*0.2)
    title(Nameplot{np})
    format_fig; colorbar;
    
    subplot(2,5,5+np)
    r=[]; pV=[];
    for nE=1:63
        if corrScale(np)==1
            [r(nE) pV(nE)]=corr(double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,18)),double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,7+np)),'rows','pairwise','type','spearman');
        else
            [r(nE) pV(nE)]=corr(double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,18)),5-double(all_Waves_byProbes2(all_Waves_byProbes2(:,[13])==nE,7+np)),'rows','pairwise','type','spearman');
        end
    end
    addpath(genpath(path_eeglab));
    temp_topo=r;
    tempplotpV=pV;
    stat_thr=fdr(tempplotpV,0.05);
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','k',24,2});
    rmpath(genpath(path_eeglab));
    caxis([-1 1]*0.7)
    format_fig; colorbar;
end

%%
prop_wave=[];
bintime=0:500:10000;
for nt=1:2
    for nstate=1:3
        temp=all_Waves_byProbes3(all_Waves_byProbes3(:,4)==nt & all_Waves_byProbes3(:,6)==nstate,end);
        ntemp=histc(temp,bintime);
        prop_wave(nt,nstate,:)=ntemp(1:end-1)./sum(ntemp)*100;
    end
end

figure; format_fig; Styles={'-','--'}; hold on;
for nt=1:2
    for nstate=1:3
        plot(bintime(1:end-1)/500,squeeze(prop_wave(nt,nstate,:)),'Color',Colors(nstate,:),'LineStyle',Styles{nt},'LineWidth',2)
    end
end