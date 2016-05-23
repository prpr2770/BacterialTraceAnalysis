% Compute -changeScore for each dataTrace signal.


%%
% close all; clear all;
dAnz = 'H:\KraljLab\2016-02-18-PROPS_CALC_Ecoli_vs_Salmonella\Noise_Results_22-Mar-2016';
% dAnz = saveDir; % while executing processVideo.m
% clear saveDir;
cd(dAnz)

saveDir = [dAnz filesep 'ChangeScores_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end

flist = dir('*.mat');
nfiles = length(flist);

% File names to be saved.
saveN = [saveDir filesep 'Extracted_data.mat'];


%% Determine sampleID from fileNames
condList = {'Ecoli';'Styph';'Mix'};
for fID = 1:nfiles
    for j = 1:length(condList)
        if regexp(flist(fID).name,condList{j})
            sampleID(fID) = j;
            break;
        end
    end
end


%% Iterate over each file to do the following:
% 1. Read the dataTraces
% 2. compute ChangeScores
% 3. archive ChangeScores in appropriate folder, with proper filenames.

for fID = 1:nfiles
    
    fname = flist(fID).name;
    fdata = load(fname);
    
    fields = fieldnames(fdata);
    if numel(fields) == 1
        fields{1}
        data = fdata(1).(fields{1});
        
    else
        warning('File has more fields than 1.')
    end
    
    ncells = length(data);
    
    for idx = 1:ncells;
        sig = data(idx).timeTrace;
        centroid = data(idx).centroid;
        
        % remove dc-component
        sig = sig - mean(sig);
        
        % extract denoised signal
        sig_dn_wav = getWaveletDenoisedTrace(sig);
        
        % obtain change-score
        sig_chgScr = getChangeDetectionScore(sig);
        sig_dn_wav_chgScr = getChangeDetectionScore(sig_dn_wav);
        
        
        % archive
        volt_sigs(idx).centroid = centroid;
        volt_sigs(idx).sig = sig ;
        volt_sigs(idx).sig_chgScr = sig_chgScr  ;
        volt_sigs(idx).sig_dn_wav = sig_dn_wav  ;
        volt_sigs(idx).sig_dn_wav_chgScr = sig_dn_wav_chgScr  ;
        
    end
    
    % save the dataFile
    fName = [saveDir filesep fname];
    save(fName,'volt_sigs');
    
    % clear the fileData for next file.
    clear volt_sigs
    
end



%%

