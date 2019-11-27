%% initial commands
clear all
clc

%% Where is the data to be processed (path and groups)

%         score_path =  '/DATA/WhiteNoiseLearning/PSN/sleepScore2/';
%         data_path =  '/DATA/WhiteNoiseLearning/PSN/EEGdata_NF/dataFiles_Test_ERP';
%         addpath('/home/andrillon/Work/toolbox/spm8/')
%         addpath(genpath('/home/andrillon/EEGgit/LSCPtools/'))
%         addpath(genpath('/home/andrillon/Work/toolbox/NoiseTools/'))
%         addpath(genpath('/home/andrillon/Work/toolbox/ExportFig/'))
%         addpath(genpath('/home/andrillon/EEGgit/eeghub/'))
%         score_path =  '/DATA/WhiteNoiseLearning/PSN/sleepScore2/';

data_path =  '/Volumes/Dell Portable Hard Drive/BackUp_DATA_tamarix/WhiteNoiseLearning/PSN/EEGdata_NF/dataFiles_Test_ERP';
addpath(genpath('/Users/tand0009/Work/local/spm12/'))
addpath(genpath('/Users/tand0009/WorkGit/LSCPtools/'))
addpath(genpath('/Users/tand0009/Work/local/NoiseTools/'))
addpath(genpath('/Users/tand0009/WorkGit/eeghub/'))
score_path =  '/Volumes/Dell Portable Hard Drive/BackUp_DATA_tamarix/WhiteNoiseLearning/PSN/sleepScore2/';

Colors2={[0.8 0.2 0],[0.5 0.3 0.8], [0.8 0.5 0.8], [0.5 0.5 0.5]};
Colors ={[1 0.5 0.1],[0.4 0.65 0.8],[0.6 0.6 0.6], [0.45 0.45 0.7]};

data_prefix='BadbeffeffMWNL2*';

% data_prefix='BadbeffeffMWNL2*';

ListSubj=dir([data_path filesep data_prefix '.mat']);

discardSWs=1;
discardSPs=1;
discardREMs=1;

%%
AEP_byCond=cell(1,4);
AEP_Full=cell(1,4);
AEP_byCond_allE=cell(1,4);
AEP_FullGA=cell(1,4);
NumT=[];
for nS=1:length(ListSubj)
    
    % load EEG
    thisSub=ListSubj(nS).name;
    fprintf('... %s ...\n',thisSub)
    D=spm_eeg_load([data_path filesep thisSub]);
    EEGchannels=D.chanlabels(match_str(D.chantype,'EEG'));
    myChan=match_str(EEGchannels,{'Cz','C3','C4'});
    
    data=D(match_str(D.chantype,'EEG'),:,:);
    
    databs=repmat(mean(data(:,D.indsample(-0.5):D.indsample(0),:),2),[1 size(data,2) 1]);
    data=data-databs;
    
    
    % load Behav
    SubName=thisSub; SubName=SubName(findstr(SubName,'WNL'):findstr(SubName,'.mat')-1);
    load([score_path 'CodeTrial_' SubName])
    
    conditions=D.conditions;
    %     if length(conditions)~=length(codeTrial.score)
    %         error('Discr')
    %     end
    %     if discardSPs==1
    %         conditions(find(max(codeTrial.Spflag(:,match_str(D.chanlabels(D.meegchannels),{'Cz','Fz'}))')==1 & codeTrial.list>=2))={'disc'};
    %     end
    %     if discardSWs==1
    %         conditions(find(max(codeTrial.SWflag(:,match_str(D.chanlabels(D.meegchannels),{'Cz','Fz'}))')==1 & codeTrial.list>=2))={'disc'};
    %     end
    %     if discardREMs==1
    %         conditions(find(codeTrial.phasicREM & codeTrial.list==3))={'disc'};
    %     elseif discardREMs==-1
    %         conditions(find(~codeTrial.phasicREM & codeTrial.list==3))={'disc'};
    %     end
    
    %     for nB=1:5
    for nLc=1:4
        if nLc==1
            myStage='0'; nL=1; maxVal=80;
            discardSWs=0;
            discardSPs=0;
            discardREMs=0;
        elseif nLc==2
            myStage='[2]';  nL=2; maxVal=150;
            discardSWs=0;
            discardSPs=0;
            discardREMs=0;
        elseif nLc==3
            myStage='[3]';  nL=2; maxVal=150;
            discardSWs=0;
            discardSPs=0;
            discardREMs=0;
        elseif nLc==4
            myStage='5'; nL=3; maxVal=80;
            discardSWs=0;
            discardSPs=0;
            discardREMs=0;
        end
        
        myT=find(cellfun(@any, regexp(conditions,sprintf('Ref[0-1]_R[0-1]_L%g_B[1-5]_S%s_MA0_E[1 3]_R[0 1 3]_C[0 1]',nL,myStage)))); %'Ref?_R?_L?_B?_S?_E?_R?_C?'
        
        if discardSWs==1 && nL==2
            myT(sum(codeTrial.evSWflag(myT,:),2)>0)=[];
        elseif discardSWs==-1 && nL==2
            myT(sum(codeTrial.evSWflag(myT,:),2)==0)=[];
        end
        if discardSPs==1 && nL==2
            myT(sum(codeTrial.evSpflag(myT,:)+codeTrial.evSpflag_F(myT,:),2)>0)=[];
        elseif discardSPs==-1 && nL==2
            myT(sum(codeTrial.evSpflag(myT,:)+codeTrial.evSpflag_F(myT,:),2)==0)=[];
        end
        if discardREMs==1 && nL==3
            myT(codeTrial.REMflag(myT)>0)=[];
        elseif discardREMs==-1 && nL==3
            myT(codeTrial.REMflag(myT)==0)=[];
        end
        
        temp=squeeze(mean(data(myChan,:,myT),1)); badT=find(max(abs(temp),[],1)>maxVal); myT(badT)=[];
        fprintf('... List %g ... numT: %g\n',nLc,length(myT))
        
        thisDIN=D.indsample(-0.1):D.indsample(0.8);
        thisDINb=D.indsample(-0.1):D.indsample(0);
        
        tempdata=[];
        tempdata(1,:,:)=squeeze(mean(data(myChan,:,myT),1)-repmat(mean(mean(data(myChan,thisDINb,myT),2),1),[1 size(data,2) 1]));
        %         filt=nan(size(tempdata));
        %         for n=1:size(tempdata,3)
        %             filt(:,:,n)=ft_preproc_highpassfilter(tempdata(:,:,n),200,1,4,'but','twopass');
        %         end
        %         tempdata=filt;
        %
        if length(myT)>20
            AEP_byCond{nLc}=[AEP_byCond{nLc} ; mean(tempdata(:,thisDIN,:),3)];
            AEP_Full{nLc}=[AEP_Full{nLc} ; mean(tempdata(:,:,:),3)];
        end
        
        tempdata=squeeze(data(:,thisDIN,myT)-repmat(mean(data(:,thisDINb,myT),2),[1 length(thisDIN) 1]));
        if length(myT)>20
            AEP_byCond_allE{nLc}=cat(3,AEP_byCond_allE{nLc}, mean(tempdata,3));
        end
        AEP_FullGA{nLc}=[AEP_FullGA{nLc} ; squeeze(mean(tempdata(myChan,:,:),1))'];
        %                 ERP_byCond_GA{nLc}{nSt,nDin}=[ERP_byCond_GA{nLc}{nSt,nDin} ; (squeeze(mean(data(myChan,thisDIN,myT),1))-repmat(mean(squeeze(mean(data(myChan,thisDINb,myT),1))),length(thisDIN),1))'];
    end
    %             ERP_whole_byCond{nLc,nSt}=[ERP_whole_byCond{nLc,nSt} ; squeeze(mean(mean(data(myChan,:,myT),1),3))];
    %             ERP_whole_byCond_GA{nLc,nSt}=[ERP_whole_byCond_GA{nLc,nSt} ; squeeze(mean(data(myChan,:,myT),1))'];
    
    NumT=[NumT ; length(myT)];
    
end

%%
finalCondNames={'WAKE','N2','N3','REM'};
myTime=-0.1:1/D.fsample:0.8;
nLc2=0;
ColorsF=[1 2 2 3];
figure; set(gcf,'Position',[ 566         601        1129         364])
for nL2=[1 2 3 4]
    nLc2=nLc2+1;
    subplot(1,4,nLc2);     format_fig;
    simpleTplot(myTime,AEP_byCond{nL2},0,Colors2{ColorsF(nL2)},0,'-',0.5,1,5);
    xlim([-0.1 0.5])
    %     ylim([-2.4 4])
    if nL2==1
        ylabel('MEPs (\muV)')
    end
    title(finalCondNames{nL2})
    xlabel('Time (s)')
    
end
% export_fig([figPath filesep 'sWNL_Test_AEP.png'])


myTimes=[[0.12 0.16];[0.19 0.25];[0.19 0.25];[0.19 0.25];[0.19 0.25]];
figure; set(gcf,'Position',[ 231         150        1464         364])
load layout_ft_1020.mat
addpath(('/Users/tand0009/Work/local/fieldtrip/'))
nLc2=0;
% nLc2=nLc2+1;
% subplot(1,5,nLc2);     format_fig;
% simpleTopoPlot_ft(squeeze(mean(mean(AEP_byCond_allE{1}(:,myTime>=myTimes(nLc2,1) & myTime<=myTimes(nLc2,2),:),2),3)),layout,'off','jet',0);

for nL2=1:4
    nLc2=nLc2+1;
    subplot(1,4,nLc2);     format_fig;
    simpleTopoPlot_ft(squeeze(mean(mean(AEP_byCond_allE{nL2}(:,myTime>=myTimes(nLc2,1) & myTime<=myTimes(nLc2,2),:),2),3)),layout,'off','jet',0);
    
end
h=colorbar('Position',[.9    0.8    0.05    0.15],'FontSize',14);

% export_fig([figPath filesep 'sWNL_Test_AEP_Topo.png'])

%
%%
finalCondNames={'WAKE','N2','N3','REM'};
myTime=-0.1:1/D.fsample:0.8;
nLc2=0;
ColorsF=[1 2 2 3];
figure; set(gcf,'Position',[ 566         601        1129         364])
for nL2=[1 2 3 4]
    nLc2=nLc2+1;
    subplot(1,4,nLc2);     format_fig;
    simpleTplot(myTime,mean(AEP_FullGA{nL2}(max(AEP_FullGA{nL2}(:,myTime>0 & myTime<0.400),[],2)<40,:)),0,Colors2{ColorsF(nL2)},0,'-',0.5,1,10);
    xlim([-0.1 0.5])
    %     ylim([-2.4 4])
    if nL2==1
        ylabel('MEPs (\muV)')
    end
    title(finalCondNames{nL2})
    xlabel('Time (s)')
    
end
% % export_fig([figPath filesep 'sWNL_Test_AEP.png'])

