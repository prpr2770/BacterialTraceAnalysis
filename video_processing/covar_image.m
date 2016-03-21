function covimg = covar_image(dat)
%COVIMG = COVAR_IMAGE(DAT)
%
% This function takes in a movie, dat, and calculates the covariance in
% time.  The output image looks at the covariance between neighboring
% pixels in x and y.
%
% Created by Joel Kralj - 2014-11-28

[ysize,xsize,nframes] = size(dat);
avgimg = mean(dat,3);

xcovimg = (zeros(ysize,xsize));
ycovimg = (zeros(ysize,xsize));
for j = 1:nframes
    xcovimg(:,1:xsize-1) = xcovimg(:,1:xsize-1) + dat(:,1:xsize-1,j).*dat(:,2:xsize,j);
    ycovimg(1:ysize-1,:) = ycovimg(1:ysize-1,:) + dat(1:ysize-1,:,j).*dat(2:ysize,:,j);
end;
xcovimg = xcovimg/j;
ycovimg = ycovimg/j;
xcovimg(:,1:xsize-1) = xcovimg(:,1:xsize-1) - avgimg(:,1:xsize-1).*avgimg(:,2:xsize);
ycovimg(1:ysize-1,:) = ycovimg(1:ysize-1,:) - avgimg(1:ysize-1,:).*avgimg(2:ysize,:);
covimg = abs((xcovimg + ycovimg)/2).^.5;

end

