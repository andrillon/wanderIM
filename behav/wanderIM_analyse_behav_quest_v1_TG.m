%%
clear all
close all

run localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

%load([data_path filesep 'CARS_quest'])
%%
all_behav_mat=[];
all_behav_mat2=[];
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
    temp_behav1=[];
    temp_behav1(1,:)=[str2num(SubID) 1 corr_go(n,1) corr_nogo(n,1) rt_gos(n,1) dprime_test(n,1) crit_test(n,1)];
    temp_behav1(2,:)=[str2num(SubID) 2 corr_go(n,2) corr_nogo(n,2) rt_gos(n,2) dprime_test(n,2) crit_test(n,2)];

    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test2, crit_test2]=calc_dprime((tp_gos==1),(tp_nogos==0));
    corr_go2=nanmean(tp_gos);
    corr_nogo2=nanmean(tp_nogos);
    
    rt_gos2=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,12))));
    temp_behav1b=[str2num(SubID) 1 corr_go2 corr_nogo2 rt_gos2 dprime_test2 crit_test2];
    
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
    temp_behav2=[];
    temp_behav2(1,:)=[squeeze(pption_MS(n,1,1:3))'];
    temp_behav2(2,:)=[squeeze(pption_MS(n,2,1:3))'];
    all_behav_mat=[all_behav_mat ; [temp_behav1 temp_behav2]];
    
    temp_behav2b=squeeze(mean(pption_MS(n,:,1:3),2))';
    all_behav_mat2=[all_behav_mat2 ; [temp_behav1b temp_behav2b]];
end
%%
M = readtable([pwd '/wanderIM/survey/data/output/MWI_Computed_Scales.csv']);
TBL=array2table(all_behav_mat,'VariableNames',{'Sub','Task','CorrGo','CorrNoGo','RTGo','dp','crit','ON','MW','MB'});
mysubs=TBL.Sub;

for ncol=1:size(M,2)
    if isnumeric(eval(sprintf('M.%s(1)',M.Properties.VariableNames{ncol})))
        eval(sprintf('TBL.%s=nan(length(mysubs),1);',M.Properties.VariableNames{ncol}))
    else
        eval(sprintf('TBL.%s=cell(length(mysubs),1);',M.Properties.VariableNames{ncol}))
        
    end
end

for ns=1:length(mysubs)
    theserows=find(M.Participant_No==mysubs(ns));
    for ncol=1:size(M,2)
        eval(sprintf('TBL.%s(TBL.Sub==mysubs(ns))=M.%s(theserows);',M.Properties.VariableNames{ncol},M.Properties.VariableNames{ncol}))
    end
end

lme_1= fitlme(TBL,'CorrGo~1+(1|Sub)');
lme_2= fitlme(TBL,'CorrGo~1+Task+(1+Task|Sub)');
% lme_3= fitlme(TBL,'CorrGo~1+Task+D2___Age+(1+Task+D2___Age|Sub)');
% lme_4= fitlme(TBL,'CorrGo~1+Task+D1___Gender+(1+Task+D1___Gender|Sub)');
TBL.D1___Gender=categorical(TBL.D1___Gender);
TBL.D3___mindfulness=categorical(TBL.D3___mindfulness);
TBL.D3___mindfulness=reordercats(TBL.D3___mindfulness,{'I have never practiced mindfulness','A few times per year or less','A few times per month','A few times per week'});
%%
fullformula='CorrNoGo~1+Task';
allvarnames=TBL.Properties.VariableNames;
for ncol=[10 11 13 14 19 21 29]
    fullformula=[fullformula '+' allvarnames{ncol}];
end
fullformula=[fullformula '+(1|Sub)'];
lme_full= fitlme(TBL,fullformula);

%% Correlation approach
TBL2=array2table(all_behav_mat2,'VariableNames',{'Sub','Task','CorrGo','CorrNoGo','RTGo','dp','crit','ON','MW','MB'});

var_tests={'CorrGo','CorrNoGo','RTGo','dp','crit','ON','MW','MB'};
for nrow=1:length(var_tests)
    eval(sprintf('varia=TBL.%s;',var_tests{nrow}))
    for ncol=1:size(M,2)
        if isnumeric(eval(sprintf('M.%s(1)',M.Properties.VariableNames{ncol})))
            eval(sprintf('predictor=TBL.%s;',M.Properties.VariableNames{ncol}))
            [r, pV]=corr(varia,predictor,'type','spearman');
            Mat_Coeff(nrow,ncol)=r;
            Mat_PVal(nrow,ncol)=pV;
        elseif strcmp(M.Properties.VariableNames{ncol},'D3___mindfulness')
            eval(sprintf('predictor=TBL.%s;',M.Properties.VariableNames{ncol}))
            predictor2(match_str(predictor,'I have never practiced mindfulness'))=1;
            predictor2(match_str(predictor,'A few times per year or less'))=2;
            predictor2(match_str(predictor,'A few times per month'))=3;
            predictor2(match_str(predictor,'A few times per week'))=4;
            [r, pV]=corr(varia,predictor2','type','spearman');
            Mat_Coeff(nrow,ncol)=r;
            Mat_PVal(nrow,ncol)=pV;
        else
            eval(sprintf('predictor=TBL.%s;',M.Properties.VariableNames{ncol}))
            [p,anovatab,stats] = kruskalwallis(varia,predictor,'off');
            Mat_Coeff(nrow,ncol)=anovatab{2,5};
            Mat_PVal(nrow,ncol)=anovatab{2,6};
        end
    end
end