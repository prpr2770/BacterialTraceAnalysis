% Script to read the SPECTG files and generate NORM_SPECTG data.

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016';
cd(dAnz)

% Accessing Data Files
flist = dir('*.mat');
nfiles = length(flist);

% Determining SpectG directory
allFiles = dir();
for fid = 1:length(allFiles)
    if length(regexp(allFiles(fid).name,'spectG'))
        saveDir = allFiles(fid).name;
    end
end
        
saveDir     % Directory containing SpectG data for SetA of files.


%% iterate over each file

for fID = 1:nfiles
    fID
    fname = flist(fID).name
    
    % --------------------------------------------------------------------
    % File names to be saved.
    fDIR_windows = [saveDir filesep fname(1:end-4)]
    
    % create file to store AggregateInfo from all Cells in File.
    aggSpectg_fileName = [fDIR_windows filesep 'aggSPECTG.mat'];
    aggSpectg_mFile = matfile(aggSpectg_fileName,'Writable',false);
   
     % --------------------------------------------------------------------
   
    % extract the aggINFO
    MEAN_SPECTG = aggSpectg_mFile.MEAN_SPECTG;
    STD_SPECTG = aggSpectg_mFile.STD_SPECTG;
    
    % --------------------------------------------------------------------
    
    % Determine number of Cells
    fdata = load(fname);
    % extract data
    fields = fieldnames(fdata);
    if numel(fields) == 1
        fields{1};
        data = fdata(1).(fields{1});
    else
        warning('File has more fields than 1.')
    end
    
    % read the denoised signals
    ncells = length(data);
    clear data fdata fields
    % --------------------------------------------------------------------
   % Read SPECTG information for each cell.
   for idx = 1:ncells

       % access SPECTG info
        trcName = sprintf('SPECTG_%d.mat',idx);
        spectg_fileName = [fDIR_windows filesep trcName];
        spectg_mFile = matfile(spectg_fileName,'Writable',false);
        
        SPECTG = spectg_mFile.SPECTG;
        numFrames = size(SPECTG,2);

                % Compute Norm_SPECTG

        MEAN_SPECTG_mat = repmat(MEAN_SPECTG,1,numFrames);
        STD_SPECTG_mat = repmat(STD_SPECTG,1,numFrames);
        
        NORM_SPECTG = (SPECTG - MEAN_SPECTG_mat)./STD_SPECTG_mat;
    
        % Archive NORM_SPECTG
        trcName = sprintf('NORM_SPECTG_%d.mat',idx);
        nrm_spectg_fileName = [fDIR_windows filesep trcName];
        nrm_spectg_mFile = matfile(nrm_spectg_fileName,'Writable',true);

        nrm_spectg_mFile.NORM_SPECTG = NORM_SPECTG;
        
        % close the matfile.
        clear nrm_spectg_mFile
        
   end
    
   

   
end