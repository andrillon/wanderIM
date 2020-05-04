%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

% state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
% cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

load([pwd filesep 'paper_SubID'])

load([root_path filesep 'surveys' filesep 'WanderIM_Surveys_April2020.mat'])


%%
Age=nan(1,length(GoodSudID));
GenderF=nan(1,length(GoodSudID));
ASRS=nan(1,length(GoodSudID));
ESS=nan(1,length(GoodSudID));
MWQ=nan(1,length(GoodSudID));
ASRS_Flag=nan(1,length(GoodSudID));
ASRS_code=[3 3 3 4 4 4];
for nS=1:length(GoodSudID)
    findSubj=find(ismember(surveys.ExternalReference,str2num(GoodSudID{nS})));
    findSubj2=find(~cellfun(@isempty,regexp({files.name},GoodSudID{nS})));
    load([files(findSubj2).folder filesep files(findSubj2).name]);
    
    Age(nS)=str2num(SubjectInfo.Age);
    GenderF(nS)=SubjectInfo.Gender=='F';
    if ~isempty(findSubj)
        ESS(nS)=sum(table2array(surveys(findSubj,find(~cellfun(@isempty,regexp(surveys.Properties.VariableNames,'^ESS_'))))),2);
        ASRS(nS)=mean(table2array(surveys(findSubj,find(~cellfun(@isempty,regexp(surveys.Properties.VariableNames,'^ASRS_'))))),2);
        MWQ(nS)=mean(table2array(surveys(findSubj,find(~cellfun(@isempty,regexp(surveys.Properties.VariableNames,'^MWQ_'))))),2);
        
        ASRS_Q=find(~cellfun(@isempty,regexp(surveys.Properties.VariableNames,'^ASRS_')));
        ARSR_PartA=ASRS_Q(1:6);
        temp=[];
        for k=1:length(ARSR_PartA)
            temp(k)=table2array(surveys(findSubj,ARSR_PartA(k)))>=ASRS_code(k);
        end
        ASRS_Flag(nS)=sum(temp)>=4;
    else
        warning(sprintf('missing subject %s!\n',GoodSudID{nS}));
    end
    
    Duration(nS)=(all_GrandEnd-all_GrandStart)/60;
end