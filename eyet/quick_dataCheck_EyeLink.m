%% import
clear all;

addpath(('/Users/Thomas/Work/local/toolbox/edf-converter-master/'))
filename='/media/tLab_BackUp1/Monash/WanderIM/eyetracker/wanderIM_eyelink_s303_30May2018-1017.edf';
myedf=Edf2Mat(filename);

%% select event
events=myedf.Events.Messages.info;
allP_onsets=[];
for nP=1:10
    allP_onsets=[allP_onsets ; match_str(events,sprintf('P%g',nP))];
end
allP_onsets=myedf.Events.Messages.time(sort(allP_onsets));

TA_probes=[];
for nP=1:length(allP_onsets)
    [tval,tidx]=findclosest(myedf.Samples.time,allP_onsets(nP));
    TA_probes=[TA_probes ; myedf.Samples.pupilSize((-1*1000:5*1000)+tidx,1)'];
end

%%
figure;
format_fig;
plot(-1:1/1000:5,mean(TA_probes,1),'LineWidth',2);
xlabel('Time from probe (s)')
ylabel('Pupil Size')
title('Probe-locked')
%%
allT_onsets=myedf.Events.Messages.time(find(~cellfun(@isempty,regexp(events,'B[1-6]_T'))));
TA_trials=nan(length(allT_onsets),length((-0.2*1000:5*1000)));
fprintf('%4.0f\n',0)
for nT=1:length(allT_onsets)
   fprintf('\b\b\b\b\b%4.0f\n',nT)
 [tval,tidx]=findclosest(myedf.Samples.time,allT_onsets(nT));
    TA_trials(nT,:)=myedf.Samples.pupilSize((-0.2*1000:5*1000)+tidx,1)';
end
%%
figure;
plot(-0.2:1/1000:5,nanmean(TA_trials,1),'LineWidth',2);
xlabel('Time from face/digit onset (s)')
ylabel('Pupil Size')
title('Stim-locked')

xlim([-0.2 5])
format_fig;
