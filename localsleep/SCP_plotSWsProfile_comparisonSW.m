%% Init
clear all;
% close all;


% localdef;
user='tand009';
Expe = 'SCP';
if strcmp(user,'GuillaumeConscience')
    addpath('/usr/local/src/NoiseTools/');
    addpath(genpath('/home/guillaume/LSCPtools'));
    addpath(genpath('/home/guillaume/Projects/toolbox/NoiseTools'),'-end')
    addpath(genpath('/home/guillaume/Projects/toolbox/spm8'),'-end')
    addpath(genpath('/home/guillaume/Projects/toolbox/fieldtrip'),'-end')
    addpath(genpath('/home/guillaume/Projects/toolbox/EEGLAB'),'-end')
    addpath(genpath('/home/guillaume/Projects/eeghub'),'-end')
    addpath(genpath('/home/guillaume/Projects/sleeptools'),'-end')
    
    % Generic Path :
    genericPath = '/data/SleepAttentionAllocation/Expe/';
    
elseif strcmp(user,'Guillaume')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\NoiseTools'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\LSCPtools\'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\spm8'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\projects\eeghub\'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\projects\sleepattention\'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\projects\sleepattentionadd\'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\projects\sleeptools\'),'-end')
    addpath(genpath('C:\Users\Guillaume\Documents\MATLAB\fieldtrip\'),'-end')
    
    % Generic Path :
    genericPath = 'D:\LSCPData\';
    
elseif strcmp(user,'tamarix')
    addpath(genpath('/home/andrillon/EEGgit/LSCPtools'));
    
    addpath(genpath('/home/andrillon/EEGgit/eeghub'),'-end')
    addpath(genpath('/home/andrillon/EEGgit/sleeptools'),'-end')
    
    addpath(genpath('/home/andrillon/Work/toolbox/NoiseTools'),'-end')
    addpath(genpath('/home/andrillon/Work/toolbox/spm8'),'-end')
    %     addpath(genpath('/home/andrillon/Work/toolbox/fieldtrip'),'-end')
    %     addpath(genpath('/home/andrillon/Work/toolbox/EEGLAB'),'-end')
    % Generic Path :
    genericPath = '/DATA/SleepAttentionAllocation/Expe/';
    
elseif strcmp(user,'tand009')
    
    addpath(genpath('/Users/tand0009/WorkGit/LSCPtools/'));
    
    addpath(genpath('/Users/tand0009/WorkGit/eeghub/'),'-end')
    
    addpath(genpath('/Users/tand0009/Work/local/NoiseTools/'),'-end')
    addpath(genpath('/Users/tand0009/Work/local/spm12/'),'-end')
    %     addpath(genpath('/home/andrillon/Work/toolbox/fieldtrip'),'-end')
    %     addpath(genpath('/home/andrillon/Work/toolbox/EEGLAB'),'-end')
    % Generic Path :
    genericPath = '/Volumes/Dell Portable Hard Drive/BackUp_DATA_tamarix/SleepAttentionAllocation/Expe/';
    
end

dataPath=[genericPath, Expe,filesep,'EEG'];
rawPath=[genericPath, Expe,filesep,'EEG',filesep,'Segmented'];
savePath=[genericPath, Expe,filesep,'DataDecoding'];
savesPath=[genericPath, Expe,filesep,'EEG',filesep,'Saves'];
modelPath=[genericPath, Expe,filesep,'EEG',filesep,'Models'];
stimPath=[genericPath, Expe,filesep,'StimsSCP'];
scoringPath=[genericPath, Expe,filesep,'Scoring',filesep,'Thomas'];
scoringSavePath=[genericPath, Expe,filesep 'Scoring' filesep 'RhythmDetection'];

questPath = [genericPath, Expe,filesep,'Behav'];


ListSubj=dir([dataPath filesep '*.raw']);
subList = [1:10,12:13,16:19,21:28];


%% Pre-definition of variables

rho_att_SleepS = cell(length(subList),4);
rho_ign_SleepS = cell(length(subList),4);

windowPeriod = cell(length(subList),4);

% Parameters

filterType = 'Ind';       % Filter to use ('Oracle' for grand average filter or 'Ind' for individual filter)

smallEpoch = 1000;           % Size of the window for the computation of pearson's R between reconstruction and real streams (correlation window)

trialPoints = [1,6001];      % Minimum and maximum timings for the beginning of correlation windows
binSize = 3000;              % Bin size for the division of the trial

% Division of the trial in sections
trialSections = trialPoints(1):binSize:trialPoints(2);

allKC=[];
allSW=[];

%% Loop on subjects
nSc=0; nSc2=0;
for nS=1:numel(subList)
    %%% Import in SPM
    subName=ListSubj(subList(nS)).name;
    SubID=subName(1:end-4);
    SubjectBehavData = load([dataPath filesep '..' filesep 'Behav' filesep SubID filesep 'Result_' SubID '.mat']);
    load([dataPath filesep '..' filesep 'Behav' filesep SubID filesep 'Excepts_' SubID '.mat']);
    fprintf('... Subject: %s\n',SubID);
    saveName=sprintf('byTr2_sleepScoring_TA_%s',SubID);
    
    load([scoringPath filesep saveName])%,'resStage','modStage')
    
    
    %%% Sleep Data
    load([scoringPath filesep saveName]);
    load([scoringSavePath filesep sprintf('byTr_SPandSW_%s',SubID)]);
    %%% Retrive behavioral data
    load([dataPath filesep '..' filesep 'Behav' filesep SubID filesep 'Result_' SubID '.mat'])
    fprintf('... ... behavioral data retrived\n')
    
    %%% Epoching and saving
    hdr=ft_read_header([dataPath filesep SubID '.raw']);
    
    
    nEc=0;
    iteChange = ones(1,4);
    tempSW=[];
    tempKC=[]; nSc2=nSc2+1;
    for nE=13:length(SubjectBehavData.TrialsCaracs)
        if exepts(nE)==1
            nEc=nEc+1;
            %       fprintf('.. %g/40 ..\n',nE)
            saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);
            load([savePath filesep saveName])
            
            %%%%% LOAD BROADBAND
            saveName2=sprintf('%s_T%02.0f_filtDataBroad',SubID,nE);
            Broad=load([savePath filesep saveName2]);
                Broad.filtData.eeg=Broad.filtData.eeg-repmat(mean(Broad.filtData.eeg(:,[62 63]),2),1,65);

            %%%
            % Keep only the slow-waves in the current trial
            trial_SO = new_slowWaves(find(new_slowWaves(:,1)==nE & ismember(new_slowWaves(:,24),[2 3])),:);
            if isempty(trial_SO)
                continue;
            end
            %%% Divide trial in epochs of 1000 ms
            hypno = fullStage{nE}(1:5:end);
            for nSW=1:size(trial_SO,1)
                onset=trial_SO(nSW,2);
                timepoints = [(-100:100)+onset];
                if min(timepoints)<1 || max(timepoints)>size(Broad.filtData.eeg,1)
                    %%% Compute delta power
                    continue;
                else
                    if max(hypno(timepoints))==3 && min(hypno(timepoints))>=2
                        tempSW(end+1,:,:)=Broad.filtData.eeg(timepoints,:)-repmat(mean(Broad.filtData.eeg(onset+(-100:-80),:)),length(timepoints),1);
                    elseif max(hypno(timepoints))==2 && min(hypno(timepoints))>=2
                        tempKC(end+1,:,:)=Broad.filtData.eeg(timepoints,:)-repmat(mean(Broad.filtData.eeg(onset+(-100:-80),:)),length(timepoints),1);
                    end
                end
            end
        else
            trial_SO=[];
        end
        fprintf('%d/%d ... %g\n',nE,length(SubjectBehavData.TrialsCaracs),size(trial_SO,1));
        SOelecs=new_slowWaves(:,1);SOelecs(SOelecs==0)=[];
           countSO2(nSc2,:)=hist(SOelecs,1:65);
 end
    if ~isempty(tempSW) && ~isempty(tempKC)
        nSc=nSc+1;
        allKC(nSc,:,:)=squeeze(nanmean(tempKC,1));
        allSW(nSc,:,:)=squeeze(nanmean(tempSW,1));
        countKC(nSc)=size(tempKC,1);
        countSW(nSc)=size(tempSW,1);
        
         allSO(nSc,:,:)=squeeze(nanmean(cat(1,tempSW,tempKC),1));
        countSO(nSc)=size(cat(1,tempSW,tempKC),1);
    end
end

%%
timeP=-1.0:1/100:1.0;
figure; format_fig;
% simpleTplot(timeP,squeeze(allKC(countKC>16,:,65)),0,[1 1 1]*0.5,0,'-',0.5,1,4);
simpleTplot(timeP,squeeze(allSO(countSO>16,:,65)),0,'k',0,'-',0.5,1,4);
% export_fig('/Users/tand0009/Desktop/SW_ERP_fromSCP.eps')
save('/Users/tand0009/Data/WanderIM/SWsleepComp/ERP_SW_Sleep','allSO','timeP','countSO')

%%
addpath(genpath('/Users/tand0009/Work/local/fieldtrip/'))
rmpath(genpath('/Users/tand0009/Work/local/spm12/'))
rmpath('/Users/tand0009/Work/local/fieldtrip/external/dmlt/external/utils/');
figure;
load('/Users/tand0009/WorkGit/projects/done/whitenoise/Pilot/ProjectData/EGI64_ft_layout.mat')
temp_topo=squeeze(nanmedian(countSO2(countSO>16,:),1))';
simpleTopoPlot_ft(temp_topo, layout,'off','parula',0,0);
% export_fig('/Users/tand0009/Desktop/SW_topo_fromSCP.eps')
