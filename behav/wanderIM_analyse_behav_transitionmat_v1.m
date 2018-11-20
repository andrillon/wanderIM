%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

%load([data_path filesep 'CARS_quest'])
%%
% transMat=zeros(20,3,3);
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    
    % probe
    %  1: num probe within block
    %  2: probe onset (expected)
    %  3: probe onset (actual)
    %  4: num block
    %  5: block type (1: FACE / 2: SQUARE)
    %  6: num trial
    %  7: resp Q1 Looking (1: Yes / 2: No)
    %  8: resp Q2 Mind-State (1: On / 2: MW / 3: MB / 4: ?)
    %  9: resp Q3 Origine (1: room / 2: personal / 3: task)
    % 10: resp Q4 Awareness (1 (fully) to 4 (not at all))
    % 11: resp Q5 Intention (1 (fully) to 4 (not at all))
    % 12: resp Q6 Engagement (1 (not) to 4 (very))
    % 13: resp Q7 Performance (1 (bad) to 4 (good))
    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))
    % 15-22: onset R1-8
    % 23-30: onset Q1-8
    % 31-38: resp number Q1-8
    for nbl=1:6
        %         for npr=1:9
        %             for nstate=1:3
        %                 for nstate2=1:3
        x=probe_res(probe_res(:,4)==nbl,32);
        p=get_transitionmat(x,1:4,0);
        transMat(n,nbl,:,:)=p;
        for nperm=1:1000
            x=probe_res(probe_res(:,4)==nbl,32);
            x=x(randperm(length(x)));
            p=get_transitionmat(x,1:4,0);
            transMatperm(n,nbl,:,:,nperm)=p;
        end
        %                 end
        %             end
        %         end
    end
%     p=get_transitionmat(probe_res(:,32),1:4);
%     transMat2(n,:,:)=p;
end

%%
figure
heatmap(squeeze(nanmean(nanmean(transMat(:,:,1:3,1:3),1),2))-squeeze(nanmean(nanmean(nanmean(transMatperm(:,:,1:3,1:3,:),5),1),2))); format_fig

for i=1:3
    for j=1:3
        myval=nanmean(nanmean(transMat(:,:,i,j),1),2);
        permval=squeeze(nanmean(nanmean(transMatperm(:,:,i,j,:),1),2));
        allval=sort([myval permval']);
        idx=find(myval==allval);
        pmc(i,j)=idx(1)/length(allval);
    end
end
