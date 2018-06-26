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
        fprintf('... %s already imported\n',subID)
    end
end

