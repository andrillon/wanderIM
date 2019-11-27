%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);
eyet_path=[root_path filesep 'eyetracker'];

files=dir([eyet_path filesep 'wanderIM_eyelink_s3*.edf']);

%% Loop on files
redo=0;
for n=1:length(files)
    subID=files(n).name;
    bound=findstr(subID,'wanderIM_eyelink_s');
    subID=subID(length('wanderIM_eyelink_s')+(1:3));
    savename=['wanderIM_eyelink_S' subID];
    
    if exist([eyet_path filesep savename '.mat'])==0 || redo==1
        fprintf('... converting %s from EDF to .mat file\n',subID)
        % import into matlab
        filename=files(n).name;
        myedf=Edf2Mat([eyet_path filesep filename]);
        
        % clean from useless info
        %     data=myedf;
        %     list_fields_toremove={'AUTHOR','AUTHOREMAIL','COPYRIGHT','VERSION','VERSIONDATE','CHANGELOG',...
        %         'RECORDING_STATES',};
        %     data=rmfiled(data,'AUTHOR');
        EL_headers=[];
        EL_headers=myedf.Header;
        EL_headers.Fs=unique(diff(myedf.timeline))*1000'; % in Hertz
        
        EL_data=[];
        EL_data.time=myedf.Samples.time;
        EL_data.pupilSize=myedf.Samples.pupilSize;
        EL_data.posX=myedf.Samples.posX;
        EL_data.posY=myedf.Samples.posY;
        
        EL_events=[];
        EL_events.Events.time=myedf.Events.Messages.time;
        EL_events.Events.type=myedf.Events.Messages.info;
        EL_events.StartRec=myedf.Events.Start.time;
        EL_events.EndRec=myedf.Events.End.time;
        
        EL_events.Fix.start=myedf.Events.Efix.start;
        EL_events.Fix.end=myedf.Events.Efix.end;
        EL_events.Fix.duration=myedf.Events.Efix.duration;
        EL_events.Fix.posX=myedf.Events.Efix.posX;
        EL_events.Fix.posY=myedf.Events.Efix.posY;
        EL_events.Fix.pupilSize=myedf.Events.Efix.pupilSize;
        
        EL_events.Blinks.start=myedf.Events.Eblink.start;
        EL_events.Blinks.end=myedf.Events.Eblink.end;
        EL_events.Blinks.duration=myedf.Events.Eblink.duration;
        
        
        EL_events.Sacc.start=myedf.Events.Esacc.start;
        EL_events.Sacc.end=myedf.Events.Esacc.end;
        EL_events.Sacc.duration=myedf.Events.Esacc.duration;
        EL_events.Sacc.posX_start=myedf.Events.Esacc.posX;
        EL_events.Sacc.posY_start=myedf.Events.Esacc.posY;
        EL_events.Sacc.posX_end=myedf.Events.Esacc.posXend;
        EL_events.Sacc.posY_end=myedf.Events.Esacc.posYend;
        EL_events.Sacc.velo=myedf.Events.Esacc.pvel;
        
        save([eyet_path filesep savename],'EL_headers','EL_data','EL_events');
    else
        fprintf('... load %s from .mat file\n',subID)
        load([eyet_path filesep savename])
    end
    % correct doubling
    if length(EL_data.time(1:2:end))==length(EL_data.time(2:2:end)) && min(EL_data.time(1:2:end)-EL_data.time(2:2:end))==0 && max(EL_data.time(1:2:end)-EL_data.time(2:2:end))==0
        warning('time doubled! correcting')
        tim=EL_data.time(1:2:end);
        pupSize=nanmean([EL_data.pupilSize(1:2:end) EL_data.pupilSize(2:2:end)],2);
        posX=nanmean([EL_data.posX(1:2:end) EL_data.posX(2:2:end)],2);
        posY=nanmean([EL_data.posY(1:2:end) EL_data.posY(2:2:end)],2);
        
        EL_data.time=tim;
        EL_data.pupilSize=pupSize;
        EL_data.posX=posX;
        EL_data.posY=posY;
    end
    
    
    % Clean ET data
    blinks=EL_events.Blinks;
    fprintf('... ... cleaning blinks (NaN replacement)\n')
    %     fprintf('%3.0f%%\n',0)
    %     data_time=EL_data.time;
    %     data_pupil=EL_data.pupilSize;
    %     data_posX=EL_data.posX;
    %     data_posY=EL_data.posY;
    %     %         beg_ERP=nan(length(blinks.start),2000);
    %     %         end_ERP=nan(length(blinks.start),2000);
    %     for nbl=1:length(blinks.start)
    %         fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
    %         startB=max(find(data_time==blinks.start(nbl))-100,1);
    %         endB=min(find(data_time==blinks.end(nbl))+100,length(data_time));
    %
    %         %%%%%
    %         %             beg_ERP(nbl,:)=data_pupil(startB(1)+(-1000:999))-nanmean(data_pupil(startB(1)+(-1000:-500)));
    %         %             end_ERP(nbl,:)=data_pupil(endB(1)+(-999:1000))-nanmean(data_pupil(endB(1)+(500:1000)));
    %
    %         data_pupil(startB(1):endB(1))=nan;
    %         data_posX(startB(1):endB(1))=nan;
    %         data_posY(startB(1):endB(1))=nan;
    %     end
    %     fprintf('... ... done\n')
    data_pupil=get_EyeLink_cleanpupil(EL_data.pupilSize,EL_headers.Fs,EL_data.time,EL_events);
    EL_data.clean_pupilSize=data_pupil;
    EL_data.filt_pupilSize=lowpass(EL_data.clean_pupilSize, EL_headers.Fs, 6, 4);
        %     EL_data.clean_posX=data_posX;
    %     EL_data.clean_posY=data_posY;
    
    save([eyet_path filesep savename '_clean'],'EL_headers','EL_data','EL_events');
    %     else
    %         fprintf('... %s already imported\n',subID)
    % end
end

