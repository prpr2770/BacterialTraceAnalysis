clear
close all

dAnz = 'Z:\Data and Analysis\KraljLabCustom\2016-02-29 PROPS vs density 2';
cd(dAnz)

saveDir = [dAnz filesep 'Results_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end
saveN = [saveDir filesep 'Extracted_data.mat'];

flist = dir('*.tif');
nfiles = length(flist);
xKeep = 91:360;
yKeep = 91:420;

%% Look for intensity dependence
illumImg = zeros(length(xKeep),length(yKeep));
for f = 1:nfiles
    tmp = double(bigread2(flist(f).name,1,25));
    tmpT = mean(tmp(xKeep,yKeep,:),3);
    illumImg = illumImg+tmpT;
end
illumImg = illumImg/nfiles;

H = fspecial('gaussian',10);
illumImg2 = imfilter(illumImg,H);

illumImg3 = imopen(illumImg2,strel('disk',20));
%%
condList = {'001';'005';'025';'125'};
for f = 1:nfiles
    for j = 1:length(condList)
        if regexp(flist(f).name,condList{j})
            cond(f) = j;
            break;
        end
    end
end

for f = 1
    flist(f).name
    dat = double(bigread2(flist(f).name,1));
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
    
    [L,num] = segment_gcamp_tester(covImg);
    
    if cond(f) == 1
        numE = 250;
        evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    elseif cond(f) == 2
        numE = 150;
        evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    elseif cond(f) == 3
        [L,num] = segment_gcamp_tester(avgImg);
        Traces = extract_traces(L,datT2,num);
    elseif cond(f) == 4
        [L,num] = segment_gcamp_tester(avgImg);
        Traces = extract_traces(L,datT2,num);
    end
    
    try
    catch
        Traces = zeros(1,nframes);
        Df = zeros(1,nframes);
        ExImage = [];
    end
    
    TracesAll{f} = Traces;
    ExImageAll{f} = ExImage;
    colorimgAll{f} = colorimg;
    
    clearvars -except *All flist nfiles saveDir dAnz f cond saveN *Keep
    
    
    save(saveN)
end

%%
indx = 0;
for f = 1:3
    ncells = length(ExImageAll{f})
    
    for c = 1:ncells
        indx = indx+1;
        if f == 1
            traceData(indx).density = 'High';
        elseif f == 2
            traceData(indx).density = 'Medium';
        elseif f == 3
            traceData(indx).density = 'Low';
        end
        
        traceData(indx).timeTrace = TracesAll{f}(c,:);
        traceData(indx).centroid = ExImageAll{1}(1).centroid;
    end
end

%%
f = 1;
figure
for j = 1:10
    plot(.8*j+mat2gray(Traces(j,:)))
    hold all
end

%%
cd(saveDir)
sflist = dir('*.mat');
load(sflist(1).name);


%%
for f = 1:3
    ncells = 241
    [tOn,tOff,numB] = find_PROPS_blinks(Traces,.13);
    blinksPerCell(f) = numB/ncells;
end


plot_Ca_blinks(Traces,tOn,tOff,[8 2])
   
    
    
    
    
    
    
    
    
    
    
    
