% Extract traces from the Salmonella and Ecoli dataset.

clear
close all

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella';
cd(dAnz)

saveDir = [dAnz filesep 'DenoisedResults_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end
saveN = [saveDir filesep 'Extracted_data.mat'];

flist = dir('*.tif');
nfiles = length(flist);


xKeep = 91:360;
yKeep = 91:420;

%%
condList = {'Ecoli';'Styph';'Mix'};
for f = 1:nfiles
    for j = 1:length(condList)
        if regexp(flist(f).name,condList{j})
            cond(f) = j;
            break;
        end
    end
end

threshold = 0.001;

for f = 3:nfiles
    flist(f).name
    
    dat = double(bigread2(flist(f).name,1));
    % denoise the signal.
    threshold = 0.001;
    dat = wavelet_denoise3D(dat,threshold);
    
    %%  Look for intensity dependence
    illumImg = zeros(length(xKeep),length(yKeep));
    tmp = dat(xKeep,yKeep,1:25);
    illumImg = mean(tmp,3);
    
    H = fspecial('gaussian',10);
    illumImg2 = imfilter(illumImg,H);
    illumImg3 = imopen(illumImg2,strel('disk',20));
    
    %%
    datT = dat(xKeep,yKeep,:);
    clear dat
    [ysize,xsize,nframes] = size(datT);
    
    datT2 = datT./repmat(illumImg3,[1 1 nframes]);
    
        covImg = covar_image(datT2);
    avgImg = mean(datT2,3);
    colorimg = zeros(ysize,xsize,3);
    colorimg(:,:,1) = 1.5*mat2gray(avgImg);
    colorimg(:,:,3) = 1.5*mat2gray(covImg);
    
    figure
    imshow(avgImg,[])
    
    %%  Extract
    [L,num] = segment_gcamp_tester(covImg);
    numE = 250;
    evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    
    %%  Catch Errors and initialize zero vectors.
    
    try
    catch
        Traces = zeros(1,nframes);
        Df = zeros(1,nframes);
        ExImage = [];
    end
    
%     TracesAll{f} = Traces;
%     ExImageAll{f} = ExImage;
%     colorimgAll{f} = colorimg;
    
    %% -------------------------------------------------
    % Save the dataTraces into respective variables and storing/archiving.
    ncells = length(ExImage);
    sampleType = cond(f);
    switch sampleType
        case 1 % ECOLI
            for indx = 1:ncells
                evalc(sprintf('ecoli_%d_TraceData(indx).timeTrace = Traces(indx,:);',f));
                evalc(sprintf('ecoli_%d_TraceData(indx).centroid = ExImage(indx).centroid;',f));
            end
            dTitle = flist(f).name(1:end-4);
            fTitle = sprintf('ecoli_%d_Traces.mat',f);
            fName = [saveDir filesep dTitle '_' fTitle];
            save(fName,sprintf('ecoli_%d_TraceData',f));
        case 2 % STYPH
            for indx = 1:ncells
                evalc(sprintf('styph_%d_TraceData(indx).timeTrace = Traces(indx,:);',f));
                evalc(sprintf('styph_%d_TraceData(indx).centroid = ExImage(indx).centroid;',f));
            end
            dTitle = flist(f).name(1:end-4);
            fTitle = sprintf('styph_%d_Traces.mat',f);
            fName = [saveDir filesep dTitle '_' fTitle];
            save(fName,sprintf('styph_%d_TraceData',f));
        case 3 % MIXTURE
            for indx = 1:ncells
                evalc(sprintf('mix_%d_TraceData(indx).timeTrace = Traces(indx,:);',f));
                evalc(sprintf('mix_%d_TraceData(indx).centroid = ExImage(indx).centroid;',f));
            end
            dTitle = flist(f).name(1:end-4);
            fTitle = sprintf('mix_%d_Traces.mat',f);
            fName = [saveDir filesep dTitle '_' fTitle];
            save(fName,sprintf('mix_%d_TraceData',f));
            
        otherwise
            warning('Error: Writing Data Files')
    end
    
    %% Clear variables and empty space
    
    clearvars -except *All flist nfiles saveDir dAnz f cond saveN threshold *Keep
    
    save(saveN)
    close all; 
end

















