% analysing volt-cal changeScore
% Script uses information about peak-locations in the changeScore for
% calcium and volt signals,
% the Ca-events  caused by volt-events : group A
% the Ca-events  NOT-caused by volt-events : group B

close all; 
clear all;

inDirName = 'H:\KraljLab\voltCal_results\';
fName = strcat(inDirName ,'changeScore_Analysis.mat');
load(fName);
totalsigs = size(data,2)


bandWidth = 16; % determines interval that determines overlap

for idx = 1:totalsigs
    idx
    v1_score_fin = data(idx).voltChangeScore ;
    v1_score_fin = data(idx).calcChangeScore ;
    
    c1_loc = data(idx).calcPeakLocs ;
    v1_loc = data(idx).voltPeakLocs ;
    
    [volt_groupA, calc_groupA, volt_groupB, calc_groupB] = bifurcateSequences(v1_loc, c1_loc, bandWidth);
    
    data(idx).voltProlongedEventLocs = volt_groupA;
    data(idx).calcProlongedEventLocs = calc_groupA;
    
    data(idx).voltImpulseEventLocs = volt_groupB;
    data(idx).calcImpulseEventLocs = calc_groupB;
    
end

save(fName,'data');


