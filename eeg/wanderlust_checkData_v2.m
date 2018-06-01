%%
clear all
% close all

addpath(genpath('/Users/Thomas/Work/local/toolbox/spm12/'));
eeg_path='/Users/Thomas/temp_data/WanderIM/';

files={'MWI303_30_05_2018.eeg'};

save_names={{'MWI303',''}};
%% Loop on files

n=1;
fprintf('... importing %s - %s\n',save_names{n}{1},save_names{n}{2})
%%% Get headers and events
S = [];
S.dataset = [eeg_path filesep files{n}];
S.outfile = [eeg_path filesep 'EEG_' save_names{n}{1} '_' save_names{n}{2}];
S.channels = 'all';
S.timewindow = [];
S.blocksize = 3276800;
S.checkboundary = 1;
S.usetrials = 1;
S.datatype = 'float32-le';
S.eventpadding = 0;
S.saveorigheader = 0;
S.conditionlabel = {'Undefined'};
S.inputformat = [];
S.continuous = true;
S.autoloc = false;
S.units='uV'; % it will be lost at montage anyway...
D = spm_eeg_convert(S);

%% find triggers
myev=D.events;
din_idx=match_str(D.chanlabels,'D1');
din_chan=D(din_idx,:,1);
din_chan=(din_chan-min(din_chan))/max(din_chan-min(din_chan));
din_thr=0.2;
%     end

din_start=find(diff(din_chan>din_thr)==1)-2; din_start(end)=[];
din_end=find(diff(din_chan>din_thr)==-1)+2; din_end(1)=[];
din_dur=(din_end-din_start)/D.fsample;
din_dist=din_start(1:end)-[0 din_start(1:end-1)];

fprintf('... ... %g triggers found\n',length(din_start));
 figure;
 plot(din_chan);
 hold on
 scatter(din_start,din_chan(din_start))
 
%%
trialonset=[myev(match_str({myev.type},'T')).time];
trialidx=trialonset*D.fsample+1;
erp_diode1=nan(length(trialonset),length((-0.2*D.fsample:10*D.fsample)));
for nTr=1:length(trialonset)
    [closestdin, closestidx]=findclosest(din_start,trialidx(nTr));
    diffdin(nTr)=closestdin-trialidx(nTr);
    temp=squeeze(D(17,closestdin+(-0.2*D.fsample:10*D.fsample),1))-mean(squeeze(D(17,closestdin+(-0.2*D.fsample:0),1)));
    erp_diode1(nTr,:)=bandpass(temp,D.fsample,0.1,30,4);
end

%%
blockonset=[myev(match_str({myev.type},'B')).time];
blockoffset=[myev(match_str({myev.type},'K')).time];
param=[];
param.method='fft';
param.mindist=1.5;
baseline_cond1=[];
baseline_cond2=[];
baseline_cond3=[];
baseline_cond4=[];
baseline_all=[];
for nBase=1:8
    pow_tp=[];
    snr_tp=[];
    kc=0;
%     for nCh=2
        for k=1:3
            kc=kc+1;
            temp=D(17,D.indsample(blockonset(nBase))+(k-1)*10*D.fsample+(1:10*D.fsample),1);
            [faxis, pow_tp(kc,:)]=get_PowerSpec(temp,D.fsample,0,0);
            [snr_tp(kc,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
        end
%     end
    if nBase==1 || nBase==5
        baseline_cond1=[baseline_cond1 ; pow_tp];
    elseif nBase==2 || nBase==6
        baseline_cond2=[baseline_cond2 ; pow_tp];
    elseif nBase==3 || nBase==7
        baseline_cond3=[baseline_cond3 ; pow_tp];
    elseif nBase==4 || nBase==8
        baseline_cond4=[baseline_cond4 ; pow_tp];
    end
    baseline_all=[baseline_all ; pow_tp];
end

%%
data=D(1:64,:,:);
data=data-repmat(mean(data,1),[64 1 1]);
baseline_allE=[];
baseline_allC=[];
baseline_SNR_all=[];
kc=0;
for nBase=1:8
    fprintf('... Baseline block %g\n',nBase)
    pow_tp=[];
    snr_tp=[];
    for nCh=1:64
        for k=1:3
            temp=data(nCh,D.indsample(blockonset(nBase))+(k-1)*10*D.fsample+(1:10*D.fsample),1);
            [faxis, pow_tp(nCh,k,:)]=get_PowerSpec(temp,D.fsample,0,0);
            [snr_tp(nCh,k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
        end
    end
    if nBase==1 || nBase==5
        baseline_allC=[baseline_allC ; 1];
    elseif nBase==2 || nBase==6
        baseline_allC=[baseline_allC ; 2];
    elseif nBase==3 || nBase==7
        baseline_allC=[baseline_allC ; 3];
    elseif nBase==4 || nBase==8
        baseline_allC=[baseline_allC ; 4];
    end
    kc=kc+1;
    baseline_allE(kc,:,:)=squeeze(mean(pow_tp,2));
    baseline_SNR_all(kc,:,:)=squeeze(mean(snr_tp,2));
end

%%
tempall=baseline_SNR_all;
load('/Users/Thomas/Work/PostDoc/Monash/Wanderlust/Analyses/BrainVision/CMA-64_REF.mat')
load('/Users/Thomas/Work/PostDoc/Monash/Wanderlust/Analyses/BrainVision/myLayout_BV64.mat')
lay=lay_BV64;
figure;
% topo 6Hz
subplot(2,2,1)
[thisf,fidx]=findclosest(faxis,6);
temptopo=squeeze(mean(tempall(:,:,fidx),1));%your 64 values;
% lay=lay_BV64;
labels=elec.label(3:65);
pos=[-elec.chanpos(3:65,2),elec.chanpos(3:65,1)];
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

% topo 6Hz
subplot(2,2,2)
[thisf,fidx]=findclosest(faxis,7.5);
temptopo=squeeze(mean(tempall(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

% topo 6Hz
subplot(2,2,3)
[thisf,fidx]=findclosest(faxis,12);
temptopo=squeeze(mean(tempall(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

% topo 6Hz
subplot(2,2,4)
[thisf,fidx]=findclosest(faxis,15);
temptopo=squeeze(mean(tempall(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

figure;
[thisf,fidx]=findclosest(faxis,13.5);
temptopo=squeeze(mean(tempall(:,:,fidx),1));%your 64 values;
% lay=lay_BV64;
labels=elec.label(3:65);
pos=[-elec.chanpos(3:65,2),elec.chanpos(3:65,1)];
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

%%
figure;
% topo 6Hz
subplot(2,2,1)
[thisf,fidx]=findclosest(faxis,13.5);
temptopo=squeeze(mean(tempall(baseline_allC==1,:,fidx),1));%your 64 values;
% lay=lay_BV64;
labels=elec.label(3:65);
pos=[-elec.chanpos(3:65,2),elec.chanpos(3:65,1)];
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])


% topo 6Hz
subplot(2,2,2)
temptopo=squeeze(mean(tempall(baseline_allC==2,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])


% topo 6Hz
subplot(2,2,3)
temptopo=squeeze(mean(tempall(baseline_allC==3,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])

% topo 6Hz
subplot(2,2,4)
temptopo=squeeze(mean(tempall(baseline_allC==4,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
%%
figure;
plot(faxis,mean(baseline_cond1),'r'); 
hold on
plot(faxis,mean(baseline_cond2),'r--'); 
plot(faxis,mean(baseline_cond3),'b'); 
plot(faxis,mean(baseline_cond4),'b--'); 
xlim([5 20])

figure; format_fig
plot(faxis,mean(baseline_all)); 
xlim([3 20])

%%

figure; 
subplot(1,3,1);
[thisf,fidx]=findclosest(faxis,6);
format_fig
simpleBarPlot(1-0.2,baseline_cond1(:,fidx),'r',0.35,'k'); 
simpleBarPlot(1+0.2,baseline_cond2(:,fidx),[1 1 1 ; 1 0 0],0.35,'k'); 

simpleBarPlot(2-0.2,baseline_cond3(:,fidx),'b',0.35,'k'); 
simpleBarPlot(2+0.2,baseline_cond4(:,fidx),[1 1 1 ; 0 0 1],0.35,'k'); 
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Face','Digit'})
ylabel('SNR F1 (6Hz)')

subplot(1,3,2);
[thisf,fidx]=findclosest(faxis,7.5);
format_fig
simpleBarPlot(1-0.2,baseline_cond1(:,fidx),'r',0.35,'k'); 
simpleBarPlot(1+0.2,baseline_cond2(:,fidx),[1 1 1 ; 1 0 0],0.35,'k'); 

simpleBarPlot(2-0.2,baseline_cond3(:,fidx),'b',0.35,'k'); 
simpleBarPlot(2+0.2,baseline_cond4(:,fidx),[1 1 1 ; 0 0 1],0.35,'k'); 
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Face','Digit'})
ylabel('SNR F2 (7.5Hz)')

subplot(1,3,3);
[thisf,fidx]=findclosest(faxis,13.5);
format_fig
simpleBarPlot(1-0.2,baseline_cond1(:,fidx),'r',0.35,'k'); 
simpleBarPlot(1+0.2,baseline_cond2(:,fidx),[1 1 1 ; 1 0 0],0.35,'k'); 

simpleBarPlot(2-0.2,baseline_cond3(:,fidx),'b',0.35,'k'); 
simpleBarPlot(2+0.2,baseline_cond4(:,fidx),[1 1 1 ; 0 0 1],0.35,'k'); 
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Face','Digit'})
ylabel('SNR IM (13.5Hz)')
%%
param=[];
param.method='fft';
param.mindist=1;

probeonset=[myev(match_str({myev.type},'P')).time];
probeoffset=[myev(match_str({myev.type},'C')).time];
IPI=probeonset(2:end)-probeoffset(1:end-1);
all_IPIpow_01=[];
all_IPIpow_02=[];
for nP=1:length(IPI)-1
    pow_tp=[];
    snr_tp=[];
    k=1;
    beg=D.indsample(probeoffset(nP+1)); beg0=beg;
    endd=D.indsample(probeoffset(nP+1))+10*D.fsample;
    temp=D(18,beg:endd,1);
    [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
    %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
    while endd+10*D.fsample<beg0+IPI(nP)*D.fsample
        beg=beg+10*D.fsample;
        endd=endd+10*D.fsample;
        temp=D(18,beg:endd,1);
        k=k+1;
        [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
        %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
    end
    all_IPIpow_02=[all_IPIpow_02 ; snr_tp];
    
    pow_tp=[];
    snr_tp=[];
    
    k=1;
    beg=D.indsample(probeoffset(nP+1)); beg0=beg;
    endd=D.indsample(probeoffset(nP+1))+10*D.fsample;
    temp=D(16,beg:endd,1);
    %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
    [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
    while endd+10*D.fsample<beg0+IPI(nP)*D.fsample
        beg=beg+10*D.fsample;
        endd=endd+10*D.fsample;
        temp=D(16,beg:endd,1);
        k=k+1;
        [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
        %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
    end
    all_IPIpow_01=[all_IPIpow_01 ; snr_tp];
end
figure;
format_fig; hold on;
plot(faxis2,mean(all_IPIpow_01),'r'); xlim([3 20])
plot(faxis2,mean(all_IPIpow_02),'b'); xlim([3 20])
xlabel('Freq (Hz)')
ylabel('SNR')

%%
probeonset=[myev(match_str({myev.type},'P')).time];
probeoffset=[myev(match_str({myev.type},'C')).time];
IPI=probeonset(2:end)-probeoffset(1:end-1);
all_IPIpow=nan(length(IPI)-1,64,2501);
fprintf('%2.0f/%2.0f\n',0,0);
for nP=1:length(IPI)-1
    for nCh=1:64
        fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nP,nCh);

        pow_tp=[];
        snr_tp=[];
        k=1;
        beg=D.indsample(probeoffset(nP+1)); beg0=beg;
        endd=D.indsample(probeoffset(nP+1))+10*D.fsample;
        temp=D(nCh,beg:endd,1);
        [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
        %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
        while endd+10*D.fsample<beg0+IPI(nP)*D.fsample
            beg=beg+10*D.fsample;
            endd=endd+10*D.fsample;
            temp=D(nCh,beg:endd,1);
            k=k+1;
            [snr_tp(k,:), faxis2, ~]=get_logSNR(temp,D.fsample,param);
            %     [faxis, pow_tp(k,:)]=get_PowerSpec(temp,D.fsample,0,0);
        end
    all_IPIpow(nP,nCh,:)=mean(snr_tp);
    end
end

%%
figure;
% topo 6Hz
subplot(2,2,1)
[thisf,fidx]=findclosest(faxis,6);
temptopo=squeeze(nanmean(all_IPIpow(:,:,fidx),1));%your 64 values;
% lay=lay_BV64;
labels=elec.label(3:65);
pos=[-elec.chanpos(3:65,2),elec.chanpos(3:65,1)];
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-1 1])

% topo 6Hz
subplot(2,2,2)
[thisf,fidx]=findclosest(faxis,7.5);
temptopo=squeeze(nanmean(all_IPIpow(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-1 1])

% topo 6Hz
subplot(2,2,3)
[thisf,fidx]=findclosest(faxis,12);
temptopo=squeeze(nanmean(all_IPIpow(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-1 1])

% topo 6Hz
subplot(2,2,4)
[thisf,fidx]=findclosest(faxis,15);
temptopo=squeeze(nanmean(all_IPIpow(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-1 1])

figure;
[thisf,fidx]=findclosest(faxis,13.5);
temptopo=squeeze(nanmean(all_IPIpow(:,:,fidx),1));%your 64 values;
simpleTopoPlot2(temptopo(1:63), pos', labels,0,[],0,lay,[]);
caxis([-1 1])

%%
trialonset=[myev(match_str({myev.type},'P')).time];
erp_diode2=[];
erp_pow2=nan(length(trialonset),2501);
erp_pow3=nan(length(trialonset),2501);

for nTr=1:length(trialonset)
%     erp_diode2=[erp_diode2 ; squeeze(D(66,D.indsample(trialonset(nTr))+(-0.2*D.fsample:1.5*D.fsample),1))-mean(squeeze(D(66,D.indsample(trialonset(nTr))+(-0.1*D.fsample:0),1)))];

    [faxis, pow]=get_PowerSpec(data(17,D.indsample(trialonset(1))+(-10*D.fsample:-1))-mean(data(17,D.indsample(trialonset(1))+(-10*D.fsample:-1))),D.fsample,0,0);
    erp_pow2(nTr,:)=(pow);
    [faxis, pow]=get_PowerSpec(data(17,D.indsample(trialonset(1))+(1:10*D.fsample))-mean(data(17,D.indsample(trialonset(1))+(1:10*D.fsample))),D.fsample,0,0);
    erp_pow3(nTr,:)=log(pow);
end

