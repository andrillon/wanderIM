%%
clear all
% close all

addpath(genpath('/Users/tand0009/Work/local/spm12/'));
eeg_path='/Users/tand0009/Data/WanderIM_ADHD/';

files={'Test_ADHD_Feb5.eeg'};

save_names={{'TestFeb5',''}};
%% Loop on files

n=1;
fprintf('... importing %s - %s\n',save_names{n}{1},save_names{n}{2})
%%% Get headers and events
S = [];
S.dataset = [eeg_path filesep files{n}];
S.outfile = [eeg_path filesep 'EEG_' save_names{n}{1} '_' save_names{n}{2}];
S.channels = 'all';
S.timewindow = [];
S.blocksize = 3276800;
S.checkboundary = 1;
S.usetrials = 1;
S.datatype = 'float32-le';
S.eventpadding = 0;
S.saveorigheader = 0;
S.conditionlabel = {'Undefined'};
S.inputformat = [];
S.continuous = true;
S.autoloc = false;
S.units='uV'; % it will be lost at montage anyway...
D = spm_eeg_convert(S);

%% find triggers
myev=D.events;
din_idx=match_str(D.chanlabels,'Diode 1');
din_chan=D(din_idx,:,1);
din_chan=(din_chan-min(din_chan))/max(din_chan-min(din_chan));
din_thr=0.2;
%     end

din_start=find(diff(din_chan>din_thr)==1)-1; %din_start(end)=[];
din_end=find(diff(din_chan>din_thr)==-1)+2; %din_end(1)=[];
din_dur=(din_end-din_start)/D.fsample;
din_dist=din_start(1:end)-[0 din_start(1:end-1)];

din_block=din_start(din_dur>0.28 & din_dur<0.35);
fprintf('... ... %g triggers found\n',length(din_start));
figure;
plot(D.time,din_chan);
hold on
scatter(D.time(din_start),din_chan(din_start),'or')
scatter(D.time(din_end),din_chan(din_end),'ob')
scatter(D.time(din_block),din_chan(din_block),'dk')

%% Check precision freq tagging
trigg_base1=din_start(din_start>din_block(1) & din_start<din_block(2));
trigg_base2=din_start(din_start>din_block(2) & din_start<din_block(3));
trigg_base3=din_start(din_start>din_block(3) & din_start<din_block(4));
trigg_base4=din_start(din_start>din_block(4) & din_start<din_block(5));

iti_base1=diff(trigg_base1);
iti_base2=diff(trigg_base2);
iti_base3=diff(trigg_base3);
iti_base4=diff(trigg_base4);

f_base1=1./(iti_base1/D.fsample);
f_base2=1./(iti_base2/D.fsample);
f_base3=1./(iti_base3/D.fsample);
f_base4=1./(iti_base4/D.fsample);

fprintf('... ... mean trigger for baseline block 1: %2.4f Hz (range period: %g-%g samples)\n',mean(f_base1),min(iti_base1),max(iti_base1))
fprintf('... ... mean trigger for baseline block 2: %2.4f Hz (range period: %g-%g samples)\n',mean(f_base2),min(iti_base2),max(iti_base2))
fprintf('... ... mean trigger for baseline block 3: %2.4f Hz (range period: %g-%g samples)\n',mean(f_base3),min(iti_base3),max(iti_base3))
fprintf('... ... mean trigger for baseline block 4: %2.4f Hz (range period: %g-%g samples)\n',mean(f_base4),min(iti_base4),max(iti_base4))

