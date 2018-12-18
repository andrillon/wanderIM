%% Step 1: 2.5 Hz high-passed filtering (%% only on the good segments, use EEGfilt for filtering)
filepath=[path,subject{i,1},'_',exp{e,1},'\',condition{c,1},'\raw\'];
filename=files(j,1).name;
EEG=load([filepath filename]);
name=fieldnames(EEG);
data=eval(['EEG.',name{1,1},'(',num2str(1),',:)']); % data = 4s segments(185x2000 double) x Segment Number)
dataFf=cell(size(data,2),1);  % start of chebyshev filter on the raw data (data=4s segments:185x2000 double)
for seg=1:size(data,2)
    eeg=data{1,seg};
    Ff=nan(185,2000);
    for ch=1:185
        if ismember(ch, badchidx)==1
            Ff(ch,:)=nan(1,2000);
        else
            fs=500; % EEG sampling rate
            Ff(ch,:)=eegfilt(eeg(ch,:),fs,2.5,0);
        end;
    end; clear ch;
    dataFf{seg,1}=Ff;
end; clear seg;

%% Step 2: Linked mastoid referencing (A1: ch94, A2: ch190)
for seg=1:size(dataFf,1)
    load insidenew185;
    load insidenew;
    chidx=insidenew(inside185);
    tempdata=dataFf{seg,1}-repmat(nanmean(dataFf{seg,1}([find(chidx==94),find(chidx==190)],:),1),185,1);
    dataFfl{seg,1}=tempdata; clear tempdata;    
end; clear seg;

%% Step 3: Theta wave detection (output: wave structure)
for seg=1:size(dataFfl,1)
    % Detect Theta Waves and Save Wave Structure at single channel level
    for ch=1:185
        EEG=dataFfl{seg,1}(ch,:);
        if ~isnan(nanmean(EEG,2))
            fs=500;
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

%% Step 4: extract the wvinfo (amp, freq) from the wavestructure, and get the amplitude threshold 
%  Step 4.1: find waves of 6-9 Hz (using the 1/period)
%  Step 4.2: extract the maxnegpkamp of all 6-9 Hz waves in each channel and each session. 
%  Step 4.3: concatenate the maxnegpkamp values of all 6-9 waves in each channel of all sessions, Then determine the amplituade threshold of each channels (Top 20% of the amplitude distribution).
%  Step 4.4: apply the amplitude threshold to extract the big (top 20%) 6-9Hz waves

%% Then you can do whatever analysis you like (wave number, wave amplitude,....)






