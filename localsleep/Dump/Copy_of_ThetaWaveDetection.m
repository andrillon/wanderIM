%% Step 1: 2.5 Hz high-passed filtering (%% only on the good segments, use EEGfilt for filtering)
%% Init
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_amp_Waves=[];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
%     left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
%     right_freq(n)=SubjectInfo.FlickerL;
        
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials
fs=128; %sampling rate changes for decimated signal
Wp=[1.0 10.0]/(fs/2); % Filtering parameters
Ws=[0.1 15]/(fs/2); % Filtering parameters
Rp=3;
Rs=25;
[n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bbp,abp]=cheby2(n,Rs,Wn); % Loses no more than 3 dB in pass band and has at least 10 dB attenuation in stop band
clear pass* stop* Rp Rs W* n;    
    


%% Step 3: Theta wave detection (output: wave structure)
for seg=1:60
    % Detect Theta Waves and Save Wave Structure at single channel level
    for ch=1:63
        datainput = squeeze(temp_data(:,:,seg));
        allavg=nanmean(datainput,2);
        dataref=datainput-repmat(allavg,[1,size(datainput,2)]);
        LoadedEEG.data=dataref; clear dataref;  
        datax = squeeze(LoadedEEG.data(ch,1:size(LoadedEEG.data,2),1));
        orig_fs = 500;
        fs=128;
        signal = resample(double(datax),fs,orig_fs);
        EEG=filtfilt(bbp, abp, signal);
%         EEG=dataFfl{seg,1}(ch,:);
        if ~isnan(nanmean(EEG,2))
            
            epochl=4;
            seg=1;
            eval(['waves=WaveDetect(EEG',',',num2str(fs),',',num2str(epochl),',',num2str(seg),');']);
            eval(['channels(',num2str(ch),',',num2str(1),').segment(',num2str(seg),',',num2str(1),').waves=waves;']);
        else
            eval(['channels(',num2str(ch),',',num2str(1),').segment(',num2str(seg),',',num2str(1),').waves=nan;']);
        end;
    end;
end;
eval(['save ',outputpath,'wavestructure_',sessname,'_',ordername,'.mat channels -mat;']);
end

%% Step 4: extract the wvinfo (amp, freq) from the wavestructure, and get the amplitude threshold 
%  Step 4.1: find waves of 6-9 Hz (using the 1/period)
%  Step 4.2: extract the maxnegpkamp of all 6-9 Hz waves in each channel and each session. 
%  Step 4.3: concatenate the maxnegpkamp values of all 6-9 waves in each channel of all sessions, Then determine the amplituade threshold of each channels (Top 20% of the amplitude distribution).
%  Step 4.4: apply the amplitude threshold to extract the big (top 20%) 6-9Hz waves

%% Then you can do whatever analysis you like (wave number, wave amplitude,....)






