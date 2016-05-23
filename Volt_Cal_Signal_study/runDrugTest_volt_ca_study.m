%{
DRUG TESTS:

Script:
1. Extracts/Reads the Volt-Calcium signals.
2. Cleans up signal:
    remove DC component
    wavelet denoising.
3. For each cell obtain the TransferFunction - model-parameters

Plots:
1. Time-series of Pre-Drug and Post-Drug : Volt-Ca
2. Update into the same Data-Structure - denoising/ transfer-function
params.
3. Model-paramters: before and after denoising
4.

data():
    .drug
    .volt
    .ca
    .volt_dn
    .ca_dn
    .tf_params

%}


% %{
clear all; close all; clc

%% Read the data file

inDIR = 'H:\KraljLab\studyDenoising\drugShock_Play2';
file_nm = 'volt_ca_data.mat';
load([inDIR filesep file_nm]);


% create a file to write into all the data
file_nm = 'volt_ca_processed.mat';
matfile_nm = fullfile(inDIR, file_nm);
matfile_obj = matfile(matfile_nm, 'Writable', true);


[numExperiments, numVideos, numCells] = size(data);

Ts = 1; % 1s sampling time.

featurespace = {};

for expID = 1:numExperiments
    expName = sprintf('expID_%d',expID);
    expDIR = [inDIR filesep expName];
    if ~exist(expDIR)
        mkdir(expDIR)
    end
    
    
    for videoID = 1:numVideos
        vidName = sprintf('videoID_%d',videoID);
        vidDIR = [expDIR filesep vidName];
        if ~exist(vidDIR)
            mkdir(vidDIR)
        end
        
        model_params = [];
        model_fitness = [];
        for cellID = 1:numCells

%             cellName = sprintf('cellID_%d',cellID);
%             cellDIR = [inDIR filesep cellName];
%             if ~exist(cellDIR)
%                 mkdir(cellDIR)
%             end

            
            sprintf('Accessing file: (expID, videoID, cellID) = (%d, %d, %d)',expID, videoID, cellID)
            volt = data(expID, videoID, cellID).volt;
            ca = data(expID, videoID, cellID).ca;
            
            if (length(volt) ~= 0 && length(ca) ~= 0 )
                %% SIGNAL DENOISING WITH DC REJECTION
                
                volt = volt - mean(volt);
                ca = ca - mean(ca);
                
                % determine the denoised signals
                dn_volt = getWaveletDenoisedTrace(volt);
                dn_ca = getWaveletDenoisedTrace(ca);
                
                % remove DC-Component from denoised signal
                dn_volt = dn_volt - mean(dn_volt);
                dn_ca = dn_ca - mean(dn_ca);
                
                
%{
                % GENERATE PLOTS AND STORE IT
                
                fig1 = figure(1)
                subplot(2,2,1)
                plot(volt,'color',[1 0 0]);
                title('voltage signal')
                subplot(2,2,3)
                plot(dn_volt,'color',[1 1 0]);
                

                subplot(2,2,2)
                plot(ca,'color',[0 0 1]);
                title('calcium signal')
                subplot(2,2,4)
                plot(dn_ca,'color',[0 1 1]);
                
                
                plt_nm = sprintf('volt_ca_signals_%d.png',cellID);
                saveas(fig1, fullfile(vidDIR, plt_nm), 'png');
                close all;
%}
                
                %% TRANSFER FUNCTION MODEL PARAMETER EXTRACTION
                
                % extract the transfer-function parameters
                num_poles = 5;
                num_zeros = 3;
                iodelay = NaN; % unknown transport delay
                
                u = dn_volt';    % input: volt
                y = dn_ca';      % output: calcium
                
                model_data = iddata(y,u,Ts);
                
                % assign the data-subsets for tasks.
                estimation = model_data;
                validation = model_data;  % determine validation dataset.
                
                % estimated - model parameters.
                
                model_est = tfest(estimation,num_poles, num_zeros, iodelay,'Ts',model_data.Ts);
                %         model_est = tfest(data,num_poles,'Ts',data.Ts);
                [model_response,fit_goodness,~] = compare(estimation, model_est);
                

                %{
                fig_model = figure();
                mp1 = plot(model_response,'r'); hold on;
                mp2 = plot(y,'g'); hold off
                legend('pred', 'msrmt');
                legend('boxoff')
                plot_title = sprintf('Model fitness: %f',fit_goodness);
                plt_nm = sprintf('model_fitness_cellID_%d_vidID_%d.png',cellID,videoID);
                title(plot_title);
                saveas(fig_model, fullfile(expDIR, plt_nm), 'png');
                
                %}
                
                [pvec, dpvec] =  getpvec(model_est);
                close all; 
                
                
                %% UPDATE THE DATA STRUCTURE
                
                data(expID, videoID, cellID).volt_dn = dn_volt;
                data(expID, videoID, cellID).ca_dn = dn_ca;
                data(expID, videoID, cellID).model_est_pvec = pvec(2:end-1);
                data(expID, videoID, cellID).model_est_dpvec = dpvec(2:end-1);
                data(expID, videoID, cellID).model_est_fit = fit_goodness;
                
                
                model_params = [model_params pvec(2:end-1)];
                model_fitness = [model_fitness fit_goodness];
            else
                warning('Error: Empty time-series found in: (expID, videoID, cellID): (%d,%d,%d) ', expID, videoID, cellID)
            end
            
        end
        
        featurespace(expID,videoID).model_params = model_params;
        featurespace(expID,videoID).model_fitness = model_fitness;
        clear model_params;
        clear model_fitness;
    end
end

% STORE DATASET INTO FILE.
save(matfile_nm,'data','featurespace','-v7.3');

% % }

%% Call Script to obtain tSNE visualization of model-parameters. 
runTSNEVis_drug_volt_ca

