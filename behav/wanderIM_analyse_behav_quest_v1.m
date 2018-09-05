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
all_behav_mat=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
%    CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
    % SART
    %  1: num block
    %  2: block cond (1: Faces / 2: Squares)
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    
    %%% cleaning too fast RTs
    RTs=[test_res(:,10)-test_res(:,8)];
    %     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %     test_res(findFalseStarts,[10 11 12])=NaN;
    %     warning('correcting for false starts')
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
        
        rt_gos(n,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,12))));
        rt_nogos(n,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,11))));
    end
    temp_behav1=[SubID corr_go(n,1) corr_go(n,2) corr_nogo(n,1) corr_nogo(n,2) rt_gos(n,1) rt_gos(n,2)];
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
    for nbt=1:2
        temp=probe_res(probe_res(:,5)==nbt,[5 31:38]);
        clear pption pption2
        for nstate=1:4
            pption(nstate)=mean(temp(:,3)==nstate);
            pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
        end
        pption_MS(n,nbt,:)=pption;
        pption_Ori(n,nbt,:)=pption2;
    end
    temp_behav2=[squeeze(pption_MS(n,1,:)) squeeze(pption_MS(n,2,:))];
    all_behav_mat=[all_behav_mat ; [temp_behav1 temp_behav2]];
end

