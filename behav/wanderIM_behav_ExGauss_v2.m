%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))
addpath(genpath(exgauss_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

% load([data_path filesep 'CARS_quest'])
%%
exgauss_param=[];
    these_RTs=cell(2,3);

for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    %     CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
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
%     for nbl=1:6
%         temp_RT=test_res(test_res(:,1)==nbl,10)-test_res(test_res(:,1)==nbl,8);
%         [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
%         exgauss_blocks=[exgauss_blocks ; [n nbl unique(test_res(test_res(:,1)==nbl,2)) X]];
%     end
    
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
    test_res(:,13:14)=nan;
    for nTask=1:2
        these_tridx=find(test_res(:,2)==nTask & test_res(:,5)~=3);
        these_trials=test_res(test_res(:,2)==nTask & test_res(:,5)~=3,:);
        
        temp_RT=these_trials(:,10)-these_trials(:,8);
        temp_RT(temp_RT<0.3)=NaN;
        [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
        
        gen_RT=X(1)+X(2)*randn(1,100000);
        flag_lapse=zeros(length(temp_RT),1);
%         flag_lapse(temp_RT>prctile(gen_RT,99.99))=1;
%         flag_lapse(temp_RT<prctile(gen_RT,0.01))=-1;
        flag_lapse(temp_RT>prctile(temp_RT,90))=1;
        flag_lapse(temp_RT<prctile(temp_RT,10))=-1;
        
        test_res(these_tridx,end)=flag_lapse;
        test_res(these_tridx,end-1)=temp_RT;
        
        exgauss_param=[exgauss_param ; [n nTask X prctile(temp_RT,10) prctile(temp_RT,90)]];
    end
    countpr=0;
    for nblock=1:6
        for npr=1:10
            these_probes=probe_res(probe_res(:,4)==nblock & probe_res(:,1)==npr,:);
            if these_probes(32)==4
                these_probes(32)=3;
            end
            these_trials=test_res(test_res(:,1)==nblock,:);
            countpr=countpr+1;
            this_pr_tridx=these_probes(6);

            temp_testres=these_trials(these_trials(:,1)==nblock & these_trials(:,4)>this_pr_tridx-21 & these_trials(:,4)<this_pr_tridx,:);
            temp_RT=temp_testres(:,10)-temp_testres(:,8);
            temp_Lapse=temp_testres(:,end);
            these_RTs{these_probes(5),these_probes(32)}=[these_RTs{these_probes(5),these_probes(32)} ; [nanmean(temp_RT) nanmean(temp_Lapse==1) nanmean(temp_Lapse==-1)]];
            
%             [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
            
        end
    end
%     for nstate=1:3
%         if length(these_RTs{nstate})>100
%             tempRT=these_RTs{nstate};
%             [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
%             exgauss_probes=[exgauss_probes ; [n nstate X length(these_RTs{nstate})]];
%         end
%     end
end

%%