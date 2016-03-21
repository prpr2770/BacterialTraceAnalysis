clear
close all

dAnz = 'Z:\Data and Analysis\KraljLabCustom\2016-02-29 CaPR tests';
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
for f = 2
    flist(f).name
    dat = double(bigread2(flist(f).name,1));
    G = dat(xKeep,yKeep,:);
    R = dat(xKeep+464,yKeep,:);
    
    colorimg = zeros(length(xKeep),length(yKeep),3);
    colorimg(:,:,1) = 1.5*mat2gray(mean(R,3));
    colorimg(:,:,2) = 1.5*mat2gray(mean(G,3));
    
    figure
    imshow(colorimg)
    
    [ysize,xsize,nframes] = size(R);
    
    covImg = covar_image(R);
    avgImg = mean(R,3);
    colorimg = zeros(ysize,xsize,3);
    colorimg(:,:,1) = 1.5*mat2gray(avgImg);
    colorimg(:,:,3) = 1.5*mat2gray(covImg);
    
    figure
    imshow(avgImg,[])
    
    [L,num] = segment_gcamp_tester(avgImg);
    tic 
    [TracesR,TracesG,ExImage] = Extract_temporal_components_two_color(G,R,30,5);
    toc
    
    time = .2*(1:900);
    figure
    for j = 1:24
        plot(time,Tracesf(j,:),time,Tracesm(j,:))
        pause
    end
    
    if cond(f) == 1
        numE = 250;
        evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    elseif cond(f) == 2
        numE = 150;
        evalc('[Traces,Df,ExImage] = Extract_temporal_components(datT2,numE,5);');
    elseif cond(f) == 3
        [L,num] = segment_gcamp_tester(avgImg);
        TracesR = extract_traces(L,R,num);
        TracesG = extract_traces(L,G,num);
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
time = (1/5)*(1:nframes);
figure
for c = 1:num
    plotyy(time,TracesR(c,:),time,TracesG(c,:))
    pause
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
   
    
    
    
    
    
    
    
    
    
    
    
