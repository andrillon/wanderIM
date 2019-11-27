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

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
myStates={'ON','MW','MB'};
pbmeDisp={'325','326','327','328','332'};
onprobe_logSNR_chan_byState=cell(2,3);
onprobe_logPow_chan_byState=cell(2,3);
onprobe_cat_chan_byState=cell(2,3);
nc=0;
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if ismember(SubID,pbmeDisp)
        fprintf('... ... PBME DISPLAY %s\n',SubID)
        continue;
    end
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
%     these_times=D.indsample(0)+1:D.indsample(20);
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    [logSNR, faxis, logpow]=get_logSNR2(temp_data,D.fsample,param);
    
%     locked_logSNR_IM=[];
%     locked_logSNR_F=[];
%     mindist=1; df=1/D.fsample;
%     length_kernel=length(-mindist+df:df:mindist-df);
%     kernel=ones(1,length_kernel);
%     kernel((-mindist+df:df:mindist-df)>-df & (-mindist+df:df:mindist-df)<df)=0;
%     kernel=kernel/sum(kernel);
%     for nE=1:63
%         for nP=1:60
%             nic=0;
%             nic2=0;
%             for ni=1:20
%                 if rem(1.5*ni,6)~=0 && rem(1.5*ni,7.5)~=0
%                     nic=nic+1;
%                     templocked=logpow(nE,faxis>=(1.5*ni-1) & faxis<=(1.5*ni+1),nP);
% %                     convlogpow= conv(templocked, kernel, 'same');
%                     locked_logSNR_IM(nic,nE,:,nP)=templocked-convlogpow;
%                 else
%                     nic2=nic2+1;
%                     templocked=logpow(nE,faxis>=(1.5*ni-1) & faxis<=(1.5*ni+1),nP);
% %                     convlogpow= conv(templocked, kernel, 'same');
%                     locked_logSNR_F(nic2,nE,:,nP)=templocked-convlogpow;
%                 end
%             end
%         end
%     end
%     onprobe_logSNR_lockedIM(n,:,:)=squeeze(mean(mean(locked_logSNR_IM,1),3));
%     onprobe_logSNR_lockedF(n,:,:)=squeeze(mean(mean(locked_logSNR_F,1),3));

    nc=nc+1;
    onprobe_logSNR(nc,:,:)=mean(logSNR,3);
    onprobe_logPow(nc,:,:)=mean(logpow,3);
    
    % Split between Face and Digit
    idx_trials=find_trials(D.conditions,'FA');
    onprobe_logSNR_FA(nc,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_FA(nc,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT');
    onprobe_logSNR_DG(nc,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_DG(nc,:,:)=mean(logpow(:,:,idx_trials),3);
    
   
        
    for nT=1:2
        if nT==1
            idx_trials=find_trials(D.conditions,'FA');
        elseif nT==2
            idx_trials=find_trials(D.conditions,'DT');
        end
        
        onprobe_logSNR_chan(nc,nT,:,:,:)=logSNR(:,:,idx_trials);
        
        for nState=1:3
            if nT==1
                idx_trials=find_trials(D.conditions,['FA_' myStates{nState}]);
            elseif nT==2
                idx_trials=find_trials(D.conditions,['DT_' myStates{nState}]);
            end
            onprobe_logSNR_chan_byState{nT,nState}=cat(3,onprobe_logSNR_chan_byState{nT,nState},logSNR(:,:,idx_trials));
            onprobe_logPow_chan_byState{nT,nState}=cat(3,onprobe_logPow_chan_byState{nT,nState},logpow(:,:,idx_trials));
            onprobe_cat_chan_byState{nT,nState}=[onprobe_cat_chan_byState{nT,nState} ; repmat(n,length(idx_trials),1)];
        end
    end
end

%%
figure;
subplot(1,2,1); format_fig;
plot(faxis,squeeze(mean(onprobe_logSNR(:,match_str(D.chanlabels,'Oz'),:),1)),'b')
xlim([1 30]); format_fig;


subplot(1,2,2); format_fig;
plot(faxis,squeeze(nanmean(onprobe_logPow(:,match_str(D.chanlabels,'Oz'),:),1)),'r')
xlim([1 30]); format_fig;


subplot(1,2,1); format_fig;
plot(faxis,squeeze(mean(onprobe_logSNR_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logSNR_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'}); format_fig;
xlabel('Freq (Hz)')
ylabel('SNR')

subplot(1,2,2); format_fig;
plot(faxis,squeeze(mean(onprobe_logPow_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logPow_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'}); format_fig;
xlabel('Freq (Hz)')
ylabel('Power')

%% by-state
figure; set(gcf,'position',[4    12   542   762])
for nT=1:2
    subplot(2,1,nT); format_fig; hold on;
    for nState=1:3
        temp=squeeze(onprobe_logPow_chan_byState{nT,nState}(match_str(D.chanlabels,'Oz'),:,:));
        temp=temp-repmat(mean(temp(faxis>15 & faxis<28,:),1),[size(temp,1) 1]);
        plot(faxis,fastsmooth(mean(temp,2),3,3,0),'Color',Colors(nState,:),'LineWidth',3)
    end
    xlim([1 30]); format_fig;
xlim([1 28])
% ylim([7 11])
format_fig;
xlabel('Freq (Hz)')
ylabel('Power')
end

figure; set(gcf,'position',[4    12   542   762])
for nT=1:2
    subplot(2,1,nT); format_fig; hold on;
    for nState=1:3
        temp=squeeze(onprobe_logSNR_chan_byState{nT,nState}(match_str(D.chanlabels,'Oz'),:,:));
%         temp=temp-repmat(mean(temp(faxis>15 & faxis<28,:),1),[size(temp,1) 1]);
        plot(faxis,mean(temp,2),'Color',Colors(nState,:),'LineWidth',3)
    end
    xlim([1 30]); format_fig;
xlim([1 28])
% ylim([7 11])
format_fig;
xlabel('Freq (Hz)')
ylabel('SNR')
end

%%
myElecs=[16:18 45:49];
f1=6;
f2=7.5;
figure; set(gcf,'position',[4    12   443   762])
subplot(2,1,1); format_fig; hold on;
forlme2=[];
for nT=1:2
    mean_temp=[];
    sem_temp=[];
    for nstate=1:3
        [~,idx1]=findclosest(faxis,f1);
        [~,idx2]=findclosest(faxis,f2);
        [~,idx3]=findclosest(faxis,f1*2);
        [~,idx4]=findclosest(faxis,f2*2);
        temp=squeeze(mean(mean(onprobe_logSNR_chan_byState{nT,nstate}(myElecs,[idx1 idx2 idx3 idx4],:),1),2));
        tempS=squeeze(onprobe_cat_chan_byState{nT,nstate});
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(temp(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
         mean_temp(nstate)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nstate)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
       
                forlme2=[forlme2 ; [repmat([nT nstate],length(temp),1) [temp tempS]]];

        %         temp=temp-repmat(mean(temp(faxis>15 & faxis<28,:),1),[size(temp,1) 1]);
        %         plot(faxis,mean(temp,2),'Color',Colors(nState,:),'LineWidth',3)
        %         simple_violinplot(temp,[],nState,'ViolinColor',Colors(nState,:),'Width',0.5,'BoxWidth',0.1);
%         simpleBarPlot(nstate+(2*nT-3)*0.2,temp,tempC,0.33,'k',[],4);
     line(nstate*[1 1]+(2*nT-3)*0.1,[-1 1]*sem_temp(nstate)+mean_temp(nstate),'Color',Colors(nstate,:),'LineWidth',3)
    end
    hold on;
    plot((1:3)+(2*nT-3)*0.1,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nstate=1:3
        if nT==1
        scatter(nstate+(2*nT-3)*0.1,mean_temp(nstate),'MarkerFaceColor',Colors(nstate,:),'MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        else
        scatter(nstate+(2*nT-3)*0.1,mean_temp(nstate),'MarkerFaceColor','w','MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        end
    end
    format_fig;
    ylim([0.5 2])
    format_fig;
    % xlabel('Freq (Hz)')
    ylabel('SNR')
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])
end
title('SNR (SSVEP)')

subplot(2,1,2); format_fig; hold on;
forlme=[];
for nT=1:2
     mean_temp=[];
    sem_temp=[];
    for nstate=1:3
    myS=unique(onprobe_cat_chan_byState{nT,nstate});
    if nT==1
            tempC=Colors(nstate,:);
        else
            tempC=[ [1 1 1]; Colors(nstate,:)];
        end
        [~,idx1]=findclosest(faxis,9.5);
        [~,idx1a]=findclosest(faxis,10.2);
        [~,idx1b]=findclosest(faxis,10.8);
        [~,idx2]=findclosest(faxis,11.5);
        [~,idx3]=findclosest(faxis,15);
        [~,idx4]=findclosest(faxis,28);
        temp=squeeze(mean(mean(onprobe_logPow_chan_byState{nT,nstate}(myElecs,[idx1:idx1a idx1b:idx2],:),1),2))-...
            squeeze(mean(mean(onprobe_logPow_chan_byState{nT,nstate}(myElecs,idx3:idx4,:),1),2));
        tempS=squeeze(onprobe_cat_chan_byState{nT,nstate});
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(temp(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
         mean_temp(nstate)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nstate)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
       
        forlme=[forlme ; [repmat([nT nstate],length(temp),1) [temp tempS]]];
        %         temp=temp-repmat(mean(temp(faxis>15 & faxis<28,:),1),[size(temp,1) 1]);
        %         plot(faxis,mean(temp,2),'Color',Colors(nState,:),'LineWidth',3)
        %         simple_violinplot(temp,[],nState,'ViolinColor',Colors(nState,:),'Width',0.5,'BoxWidth',0.1);
%         simpleBarPlot(nstate+(2*nT-3)*0.2,temp,tempC,0.33,'k',[],4);
     line(nstate*[1 1]+(2*nT-3)*0.1,[-1 1]*sem_temp(nstate)+mean_temp(nstate),'Color',Colors(nstate,:),'LineWidth',3)
    end
    hold on;
    plot((1:3)+(2*nT-3)*0.1,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nstate=1:3
        if nT==1
        scatter(nstate+(2*nT-3)*0.1,mean_temp(nstate),'MarkerFaceColor',Colors(nstate,:),'MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        else
        scatter(nstate+(2*nT-3)*0.1,mean_temp(nstate),'MarkerFaceColor','w','MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        end
    end
    
    format_fig;
    ylim([1 1.8])
    format_fig;
    % xlabel('Freq (Hz)')
    ylabel('Power')
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])
end
title('POWER (ALPHA)')

tbl_alpha=array2table(forlme,'VariableNames',{'Task','MS','ALPHA','SubID'});
tbl_alpha.MS=categorical(tbl_alpha.MS);
tbl_alpha.Task=categorical(tbl_alpha.Task);
tbl_alpha.MS=reordercats(tbl_alpha.MS,[2 1 3]);
mdl_0= fitlme(tbl_alpha,'ALPHA~1+(1|SubID)');
mdl_1= fitlme(tbl_alpha,'ALPHA~Task+(1|SubID)');
mdl_2= fitlme(tbl_alpha,'ALPHA~Task+MS+(1|SubID)');
mdl_3= fitlme(tbl_alpha,'ALPHA~Task*MS+(1|SubID)');


tbl_f=array2table(forlme2,'VariableNames',{'Task','MS','F','SubID'});
tbl_f.MS=categorical(tbl_f.MS);
tbl_f.Task=categorical(tbl_f.Task);
tbl_f.MS=reordercats(tbl_f.MS,[2 1 3]);
mdl2_0= fitlme(tbl_f,'F~1+(1|SubID)');
mdl2_1= fitlme(tbl_f,'F~Task+(1|SubID)');
mdl2_2= fitlme(tbl_f,'F~Task+MS+(1|SubID)');
mdl2_3= fitlme(tbl_f,'F~Task*MS+(1|SubID)');