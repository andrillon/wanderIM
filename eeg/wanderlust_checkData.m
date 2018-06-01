%%
clear all
% close all

addpath(genpath('/Users/Thomas/Work/local/toolbox/spm12/'));
eeg_path='/Users/Thomas/temp_data/WanderIM/';

files={'MWI301_DB.eeg'};

save_names={{'MWI001','1'}};
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

%%
[faxis, pow]=get_PowerSpec(D(17,:,1),D.fsample,1,1); xlim([0.5 30])

%%
trialonset=[myev(match_str({myev.type},'T')).time];
erp_diode1=nan(length(trialonset),length((-0.2*D.fsample:1.5*D.fsample)));
erp_pow=nan(length(trialonset),2501);
for nTr=1:length(trialonset)
    erp_diode1(nTr,:)=squeeze(D(17,D.indsample(trialonset(nTr))+(-0.2*D.fsample:1.5*D.fsample),1))-mean(squeeze(D(17,D.indsample(trialonset(nTr))+(-0.1*D.fsample:0),1)));
    

end

%%
data=D(1:63,:,1);
data=data-repmat(mean(data,1),[size(data,1) 1 1]);
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
    erp_pow3(nTr,:)=(pow);
end

%% find triggers
clear din_*
%     if strcmp(subID,'S301')
din_idx=match_str(D.chanlabels,'Diode 1');
din_chan=D(din_idx,:,1);
din_chan=(din_chan-min(din_chan))/max(din_chan-min(din_chan));
din_thr=0.2;
%     end

din_start=find(diff(din_chan>din_thr)==1)-2; din_start(end)=[];
din_end=find(diff(din_chan>din_thr)==-1)+2; din_end(1)=[];
din_dur=(din_end-din_start)/D.fsample;
din_dist=din_start(1:end)-[0 din_start(1:end-1)];

din_trans=din_start(din_dist>200);

fprintf('... ... %g triggers found\n',length(din_start));
 figure;
 plot(din_chan);
 hold on
 scatter(din_trans,din_chan(din_trans))
 
 %%
 SR=D.fsample;
 param.method='fft';
 param.mindist=1;
 param.taper=1;
 figure;
 cch=0;
 logSNR=[];
%  for nCh=[29 30 31 61 62 63]
     cch=cch+1;
 for n=1:8
     templogSNR=[];
     for nw=1:3
         data=mean(D([29 30 31],din_trans(n)+((nw-1)*10*D.fsample:nw*10*D.fsample),1),1);
         [templogSNR(nw,:), faxis, logpow]=get_logSNR(data,SR,param);
     end
     if n<5
         logSNR(n,:,1,cch)=mean(templogSNR);
     else
         logSNR(n-4,:,2,cch)=mean(templogSNR);
     end
     %      [faxis, pow]=get_PowerSpec(D(30,din_trans(n)+(0:30*D.fsample),1),D.fsample,1,0);
 end
%  end
 for nsub=1:4
     subplot(1,4,nsub)
     plot(faxis,mean(mean(logSNR(nsub,:,:,:),4),3));
     xlim([0.5 30])
 ylim([-2 3])    
 end
 
%%
din_tr=din_start(find(din_dist>150)+1);
epoched_data=[];
for n=1:length(din_tr)-1
    for nCh=1:66
        epoched_data(nCh,:,n)=D(nCh,din_tr(n)-D.fsample:din_tr(n)+30*D.fsample,1);
        epoched_data(nCh,:,n)=epoched_data(nCh,:,n)-mean(epoched_data(nCh,1:D.fsample,n));
    end
end

%%
condIM=[2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
stimIM=[1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2];
clear snr pow
for n=1:size(epoched_data,3)
    tp=squeeze(epoched_data(30,D.fsample:31*D.fsample,n));
    for k=1:1
        tp2=tp((k-1)*15000+1:k*15000);
        kernel=[-.25 -.25 0 0 1 0 0 -.25 -.25];
        [faxis, tp_powspec]=get_PowerSpec_new(tp2,1/D.fsample,length(tp2)/D.fsample,0,0);
        snr(n,:,k)= conv(log(tp_powspec), kernel, 'same');
        pow(n,:,k)=tp_powspec;
    end
end
figure;
plot(faxis,squeeze(mean(mean(snr,3),1)));
xlim([1 15])

figure;
subplot(1,2,1)
plot(faxis,squeeze(mean(mean(snr(stimIM==1,:,:),3),1)),'b'); hold on;
plot(faxis,squeeze(mean(mean(snr(stimIM==2,:,:),3),1)),'r');
xlim([1 15])

subplot(1,2,2)
plot(faxis,squeeze(mean(mean(snr(condIM==1,:,:),3),1)),'b'); hold on;
plot(faxis,squeeze(mean(mean(snr(condIM==2,:,:),3),1)),'r');
xlim([1 15])

%%
files={'Thomas Round 3 20.04.18.eeg'};

save_names={{'test','IM2'}};
% Loop on files
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

%%
data=D(1:64,:,1);
% data=data-repmat(mean(data,1),[64 1]);
[myfaxis, mypow]=get_PowerSpec(data(30,:),D.fsample,0,0); xlim([0.5 30])
% kernel=[-.25 -.25 0 0 0 1 0 0 0 0-.25 -.25];
% mysnr= conv(log(pow), kernel, 'same');
figure
plot(myfaxis,(mypow))
xlim([2 25])
% 
% figure
% plot(myfaxis,(mysnr))
% xlim([2 25])
