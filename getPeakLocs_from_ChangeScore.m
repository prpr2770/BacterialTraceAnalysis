function [ pks, locs] = getPeakLocs_from_ChangeScore(sig,leftMargin,minPeakHeight, minPeakProminence)
% output: 
% pks: 
% locs: 

T = size(sig,2);
rightMargin = T-leftMargin;
[maxVal, maxIdx] = max(sig(1,leftMargin:rightMargin));
%obtain peak values

try
[pks, locs] = findpeaks(sig(1,leftMargin:rightMargin),'MinPeakHeight',minPeakHeight*maxVal,'MinPeakProminence',minPeakProminence*maxVal);
catch
   warning('Error finding peaks: getPeakLocs_from_ChangeScore() ') 
end

end
