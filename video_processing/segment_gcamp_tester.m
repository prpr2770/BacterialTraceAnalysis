function [L, num] = segment_gcamp_tester(dat,disp_images,tester)

% FUNCTION [L, NUM] = SEGMENT_GCAMP(DAT)
%
% This function is written to segment bacteria that are expressing a strong
% fluorescent marker.  The current parameters are optimized for a 100x
% resolution (160 nm per pixel).
%
% L   = Labeled Image from watershed command
% NUM = Number of cells found in image
%
% DISP_IMAGES - 'yes' or 'no' to display the output images
%
% Created by Joel Kralj, 2014-12-12

if nargin < 1
    fprint('No image selected for segmentation');
elseif nargin < 2
    disp_images = 'yes';
end

if strcmpi(disp_images,'yes')
    dispB = 1;
else
    dispB = 0;
end


dat = double(dat);
dat2 = mat2gray(dat);

bg = imopen(dat2,strel('disk',20));
I = mat2gray(dat2-bg);
I2 = I;
I2(I < .071) = 0;
I3 = wiener2(I2,[5,5]);
I4 = mat2gray(I3.^0.3);
%%
se1 = strel('square',11);
se2 = strel('square',4);

y1 = 2*I4 - imdilate(I4,se1);
y1(y1<0) = 0;
y1(y1>1) = 1;
y2 = imdilate(I4,se2) - y1;
y3 = mat2gray(y2.^1);
y4 = im2bw(y3,.0483 );
y5 = bwareaopen(y4,49);

L = watershed(y5);
tooSmall = find(cell2mat(struct2cell(regionprops(L,'Area'))) < 25);
tooBig = find(cell2mat(struct2cell(regionprops(L,'Area'))) > 1000);
for k = 1:length(tooSmall)
    badpixels = find(L(:) == tooSmall(k));
    L(badpixels) = 0;
end
for k = 1:length(tooBig)
    badpixels = find(L(:) == tooBig(k));
    L(badpixels) = 0;
end

L(find(L == 1)) = 0;
labeledImage = label2rgb(L, @spring, 'c', 'shuffle');

[L, num] = bwlabel(L);

mask = im2bw(L, 1);
colorimgtmp = 1.1*repmat(I,[1 1 3]);
mask2 = .2*grs2rgb(double(mask),copper,0,1);
colorimg = colorimgtmp+mask2;

if dispB
    figure, imshow([labeledImage 255*colorimg])
    drawnow
end




%% Old code

% bw = im2bw(I, graythresh(I)*.5);
% % bw = im2bw(dat3,graythresh(dat3)-.05);
% bw2 = imfill(bw, 'holes');
% bw3 = imopen(bw2, strel('disk',7));
% bw4 = bwareaopen(bw3, 49);
% 
% bw4_perim = bwperim(bw4);
% 
% maxs = imextendedmax(dat2,  .01);
% maxs2 = imclose(maxs, strel('disk',1));
% maxs3 = imfill(maxs2, 'holes');
% maxs4 = bwareaopen(maxs3, 1);
% 
% Jc = imcomplement(I);
% I_mod = imimposemin(Jc, ~bw4 | maxs4);
% 
