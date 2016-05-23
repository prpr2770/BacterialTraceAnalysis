% Script to read the FEATURE files and generate NORM_FEATURE data.

dAnz = 'H:\KraljLab\2016-02-18-PROPS_CALC_Ecoli_vs_Salmonella\ChangeScores';
cd(dAnz)

% Accessing Data Files
flist = dir('*.mat');
nfiles = length(flist);

% Determining FEATRE directory
allFiles = dir();
for fid = 1:length(allFiles)
    % ASSUMING ONLY ONE FEATURE FOLDER EXISTS
    if length(regexp(allFiles(fid).name,featureType));
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
    aggFEATURE_fileName = [fDIR_windows filesep 'aggFEATURE.mat'];
    aggFEATURE_mFile = matfile(aggFEATURE_fileName,'Writable',false);
   
     % --------------------------------------------------------------------
   
    % extract the aggINFO
    MEAN_FEATURE = aggFEATURE_mFile.MEAN_FEATURE;
    STD_FEATURE = aggFEATURE_mFile.STD_FEATURE;
    
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
   % Read FEATURE information for each cell.
   for idx = 1:ncells

       % access FEATURE info
        trcName = sprintf('FEATURE_%d.mat',idx);
        feature_fileName = [fDIR_windows filesep trcName];
        feature_mFile = matfile(feature_fileName,'Writable',false);
        
        FEATURE = feature_mFile.FEATURE;
        numFrames = size(FEATURE,2);

                % Compute Norm_FEATURE

        MEAN_FEATURE_mat = repmat(MEAN_FEATURE,1,numFrames);
        STD_FEATURE_mat = repmat(STD_FEATURE,1,numFrames);
        
        NORM_FEATURE = (FEATURE - MEAN_FEATURE_mat)./STD_FEATURE_mat;
    
        % Archive NORM_FEATURE
        trcName = sprintf('NORM_FEATURE_%d.mat',idx);
        nrm_feature_fileName = [fDIR_windows filesep trcName];
        nrm_feature_mFile = matfile(nrm_feature_fileName,'Writable',true);

        nrm_feature_mFile.NORM_FEATURE = NORM_FEATURE;
        
        % close the matfile.
        clear nrm_feature_mFile
        
   end
    
   

   
end