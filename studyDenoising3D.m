% 3D-wavelet denoising : For extracting noise-free data.

close all; clear all;
loadData=0;% comment below line after  1 run.

if loadData == 0
    try
        inDIR = 'H:\KraljLab\denoising3D\';
        fname = strcat(inDIR,'PROPS_118_0-2_561nm_200ms_05.tif');
        data = double(bigread2(fname,1));
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

outDIR = 'H:\KraljLab\denoising3D\';
% =========================================================================

threshold = 0.001;
filtered_data = wavelet_denoise3D(data,threshold);
saveImages = [outDIR filesep 'all_frames_data.mat'];
save(saveImages,'data','filtered_data');


% ---------------------------------------------------
% write the filtered data as a tiff file. 
fname = strcat(inDIR,'PROPS_118_0-2_561nm_200ms_05_3DwavDenoised_double.tif');

% % % using TIFF library
% % fileTIF = Tiff(fname,'w');
% % fileTIF.write(filtered_data);
% % fileTIF.close();
% % % imagesc(imread(fname));

% convert data into uint16 format
filtered_data_uint16 = uint16(filtered_data);

% export image_data into TIFF format. 
image_data = filtered_data;
nframes = size(image_data,3);
imwrite(image_data(:,:,1),fname);
for j = 2:nframes
  imwrite(image_data(:,:,j),fname,'WriteMode','append')
end


% =======================================================================
saveDir = [inDIR filesep 'Results_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end
saveN = [saveDir filesep 'Extracted_data.mat'];
saveTraces = [saveDir filesep 'cell_traces.mat'];

% =======================================================================
xKeep = 91:240;
yKeep = 91:320;

%% Look for intensity dependence
illumImg = zeros(length(xKeep),length(yKeep));
illumImg = mean(filtered_data(xKeep,yKeep,1:25),3);

H = fspecial('gaussian',10);
illumImg2 = imfilter(illumImg,H);

illumImg3 = imopen(illumImg2,strel('disk',20));

%%

for f=1:2
    if f==1
        datT = data(xKeep,yKeep,:);
    elseif f==2
        datT = filtered_data(xKeep,yKeep,:);
    end
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
    
    % ---------------------------------------
    % which one to use?
    
    [L,num] = segment_gcamp_tester(covImg);
    numE = 250;
    evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    
    % ---------------------------------------
    
    Traces = extract_traces(L,datT2,num);
    
    
    try
    catch
        Traces = zeros(1,nframes);
        Df = zeros(1,nframes);
        ExImage = [];
    end
    
    TracesAll{f} = Traces;
    ExImageAll{f} = ExImage;
    colorimgAll{f} = colorimg;
    
    
    save(saveN)
end




%%
for f = 1:2
    ncells = length(ExImageAll{f});
    switch f
        case 1
            for indx = 1:ncells
                denoisedTraces(indx).timeTrace = TracesAll{f}(indx,:);
                denoisedTraces(indx).centroid = ExImageAll{f}(indx).centroid;
            end
        case 2
            
            for indx = 1:ncells
                noisyTraces(indx).timeTrace = TracesAll{f}(indx,:);
                noisyTraces(indx).centroid = ExImageAll{f}(indx).centroid;
            end
        otherwise
            warning('invailid file-numbers')
    end
end

save(saveTraces,'denoisedTraces', 'noisyTraces');