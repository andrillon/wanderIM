%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);
dataADHD_path=[root_path_adhd filesep 'behav/'];
filesADHD=dir([dataADHD_path filesep 'MWADHD_behavres_s4*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

%load([data_path filesep 'CARS_quest'])
%%
SART_perf_controls=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    
        age=str2num(SubjectInfo.Age);
gender=strcmp(SubjectInfo.Gender,'F');

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
    %%% Only digit blocks
    test_res(test_res(:,2)==1,:)=[];
    
    %%% cleaning too fast RTs
    RTs=test_res(:,10)-test_res(:,8);
    RTs(RTs<0.2)=NaN;
    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test, crit_test]=calc_dprime2((tp_gos==1),(tp_nogos==0));
    corr_go=nanmean(tp_gos);
    corr_nogo=nanmean(tp_nogos);
    
    rt_gos=nanmean(RTs(~isnan(test_res(:,12))));
    rt_nogos=nanmean(RTs(~isnan(test_res(:,11))));
    
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
    probe_res(probe_res(:,5)==1,:)=[];
    temp_probe=probe_res(:,[5 31:38]);
    mystates={1,2,[3 4]};
    clear pption pption2
    for nstate=1:3
        pption(nstate)=nanmean(ismember(temp_probe(:,3),mystates{nstate}));
        pption2(nstate)=nanmean(temp_probe(temp_probe(:,3)==2,4)==nstate);
    end
    pption_MS=pption;
    pption_Ori=pption2.*pption(:,2);
    vigilance=nanmean(temp_probe(:,end));
    
    SART_perf_controls=[SART_perf_controls ; [str2num(SubID) age gender corr_go corr_nogo rt_gos rt_nogos pption_MS pption_Ori vigilance]];
end
SART_perf_controls_table=array2table(SART_perf_controls,'VariableNames',{'SubID','Age','Gender','PerfGO','PerfNOGO','rtGO','rtNOGO','ppON','ppOFF','ppMB','ppDIST','ppMW','ppINT','VIG'});

%%
SART_perf_ADHD=[];
for n=1:length(filesADHD)
    % load
    load([dataADHD_path filesep filesADHD(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    
    age=str2num(SubjectInfo.Age);
gender=strcmp(SubjectInfo.Gender,'F');

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
    %%% Only digit blocks
    test_res(test_res(:,2)==1,:)=[];
    
    %%% cleaning too fast RTs
    RTs=test_res(:,10)-test_res(:,8);
    RTs(RTs<0.2)=NaN;
    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test, crit_test]=calc_dprime2((tp_gos==1),(tp_nogos==0));
    corr_go=nanmean(tp_gos);
    corr_nogo=nanmean(tp_nogos);
    
    rt_gos=nanmean(RTs(~isnan(test_res(:,12))));
    rt_nogos=nanmean(RTs(~isnan(test_res(:,11))));
    
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
    probe_res(probe_res(:,5)==1,:)=[];
    temp_probe=probe_res(:,[5 31:38]);
    mystates={1,2,[3 4]};
    clear pption pption2
    for nstate=1:3
        pption(nstate)=nanmean(ismember(temp_probe(:,3),mystates{nstate}));
        pption2(nstate)=nanmean(temp_probe(temp_probe(:,3)==2,4)==nstate);
    end
    pption_MS=pption;
    pption_Ori=pption2.*pption(:,2);
    vigilance=nanmean(temp_probe(:,end));
    
    SART_perf_ADHD=[SART_perf_ADHD ; [str2num(SubID) age gender corr_go corr_nogo rt_gos rt_nogos pption_MS pption_Ori vigilance]];
end
SART_perf_ADHD_table=array2table(SART_perf_ADHD,'VariableNames',{'SubID','Age','Gender','PerfGO','PerfNOGO','rtGO','rtNOGO','ppON','ppOFF','ppMB','ppDIST','ppMW','ppINT','VIG'});