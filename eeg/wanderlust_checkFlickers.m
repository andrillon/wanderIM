%%
clear all
% close all

addpath(genpath('/Users/Thomas/Work/local/toolbox/spm12/'));
eeg_path='/Users/Thomas/temp_data/WanderIM/';

files={'TestFlickers1.eeg'};

save_names={{'TestFlickers',''}};
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
din_idx=match_str(D.chanlabels,'D1');
din_chan=D(din_idx,:,1);
din_chan=(din_chan-min(din_chan))/max(din_chan-min(din_chan));
din_thr=0.2;
%     end

din_start=find(diff(din_chan>din_thr)==1)-2; din_start(end)=[];
din_end=find(diff(din_chan>din_thr)==-1)+2; din_end(1)=[];
din_dur=(din_end-din_start)/D.fsample;
din_dist=din_start(1:end)-[0 din_start(1:end-1)];

fprintf('... ... %g triggers found\n',length(din_start));
figure;
plot(din_chan);
hold on
scatter(din_start,din_chan(din_start))

%%
match_str({myev.type},'B')
begBlocks=match_str({myev.type},'B');
endBlocks=match_str({myev.type},'K');
myTimings=[myev.time];
BlockTimings=myTimings(begBlocks(1:8));
BlockEndings=myTimings(endBlocks(1:8));

for nB=1:8
   tempDINs=find(din_start>D.indsample(BlockTimings(nB)) & din_start<D.indsample(BlockEndings(nB)));
   durDINs=din_end(tempDINs)-din_start(tempDINs);
   keep=durDINs>30 & durDINs<100;
    DIN_perblock{nB}={din_start(tempDINs(keep)) din_end(tempDINs(keep))};
end

%%
Base1_beg=34220;
Base1_end=49167;

Base2_beg=52816;
Base2_end=67785;

Base3_beg=80778;
Base3_end=85616;

Base4_beg=91754;
Base4_end=104765;

Base5_beg=114415;
Base5_end=128227;

Base6_beg=133949;
Base6_end=148895;

Base7_beg=158237;
Base7_end=173083;

Base8_beg=34220;
Base8_end=49167;



