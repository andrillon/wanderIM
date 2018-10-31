% EEGLAB history file generated on the 18-Oct-2018
% ------------------------------------------------

EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_fileio('/Users/tand0009/Data/WanderIM/eeg/MWI322_27_06_2018.eeg');
EEG.setname='tr';
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','/Users/tand0009/Work/local/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );
