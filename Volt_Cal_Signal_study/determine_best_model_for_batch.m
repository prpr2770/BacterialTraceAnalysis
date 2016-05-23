%{
Script to determine the best family of mathematical models, to represent
the volt-ca dynamics in the cells.
%}


% %{
clear all; close all; clc

%% Read the data file

inDIR = 'H:\KraljLab\studyDenoising\drugShock_Play2';
file_nm = 'volt_ca_processed.mat';

load([inDIR filesep file_nm]);

[numExperiments, numVideos, numCells] = size(data);
sprintf('numExperiments: %d, numVideos: %d, numCells: %d',numExperiments, numVideos, numCells);
pause

Ts = 1; % 1s sampling time.

max_poles = 6;

% data-structure to store the fitness matrix for each cell. 
model_family_fitness = zeros(max_poles, max_poles, numCells, numVideos, numExperiments);% (poles, zeros, cellID) = data


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
            
            cellName = sprintf('cellID_%d',cellID);
            cellDIR = [inDIR filesep cellName];
            if ~exist(cellDIR)
                mkdir(cellDIR)
            end
            
   
            %% Iterate over every cell: 
            
            sprintf('Accessing file: (expID, videoID, cellID) = (%d, %d, %d)',expID, videoID, cellID)
            
            
            volt = data(expID, videoID, cellID).volt_dn';
            ca = data(expID, videoID, cellID).ca_dn';
            
            if length(ca) ~= 0 && length(volt) ~= 0
                        %% TRANSFER FUNCTION MODEL PARAMETER EXTRACTION
                
                for num_poles = 1:max_poles
                    for num_zeros = (num_poles-1):-1:0
                        
                        
                        % extract the transfer-function parameters
                        % num_poles
                        % num_zeros
                        sprintf('Processing: (cellID, num_poles, num_zeros) = (%d, %d, %d)',cellID, num_poles, num_zeros)

                        iodelay = NaN; % unknown transport delay
                        
                        u = volt;    % input: volt
                        y = ca;      % output: calcium
                        
                        model_data = iddata(y,u,Ts);
                        
                        % assign the data-subsets for tasks.
                        estimation = model_data;
                        validation = model_data;  % determine validation dataset.
                        
                        % estimated - model parameters.
                        
                        model_est = tfest(estimation,num_poles, num_zeros, iodelay,'Ts',model_data.Ts);
                        %         model_est = tfest(data,num_poles,'Ts',data.Ts);
                        [model_response,fit_goodness,~] = compare(estimation, model_est);
                        
                        model_family_fitness(num_poles, (num_zeros + 1) ,cellID, videoID, expID) = fit_goodness;
                        
                        
                        %{
                        fig_model = figure();
                        mp1 = plot(model_response,'r'); hold on;
                        mp2 = plot(y,'g'); hold off
                        legend('pred', 'msrmt');
                        legend('boxoff')
                        plot_title = sprintf('Model fitness[cell-ID: %d](poles:%d, zeros:%d): %f',cellID,num_poles, num_zeros, fit_goodness);
                        title(plot_title);
                        
                        plt_nm = sprintf('fitness_vidID_%d_poles_%d_zeros_%d.png',videoID, num_poles, num_zeros);
                        saveas(fig_model, fullfile(cellDIR, plt_nm), 'png');
                        %}
                        
                        [pvec, dpvec] =  getpvec(model_est);
                        close all;
                        
                        
                        
                        
                    end
                end
            end

            %{
            %% Plot the fitness function for the cells with variation in model type
            %% as a grid. 
            XY_data = model_family_fitness(:,:,cellID);
            poles_row_indices = 1:max_poles;
            zeros_col_indices = 0:max_poles-1;
            
            modelfits_plot = surf(zeros_col_indices, poles_row_indices,XY_data );
            ylabel('poles');
            xlabel('zeros');
            plt_nm = sprintf('model_fitness_expID_%d_vidID_%d_cellID_%d.png',expID,videoID,cellID);
            saveas(modelfits_plot, fullfile(cellDIR, plt_nm), 'png');
            %}
            
            
        end
    end
end



% STORE DATASET INTO FILE.

% create a file to write into all the data
file_nm = 'model_order_study.mat';
matfile_nm = fullfile(inDIR, file_nm);
matfile_obj = matfile(matfile_nm, 'Writable', true);
save(matfile_nm,'model_family_fitness','-v7.3');


%% Determine best model-order for the cells. 
A = (model_family_fitness > 50);
count_good_models = sum(sum(sum(A,3),4),5);

poles_row_indices = 1:max_poles;
zeros_col_indices = 0:max_poles-1;

modelfits_plot = surf(zeros_col_indices, poles_row_indices,count_good_models);
ylabel('poles');
xlabel('zeros');
plt_nm = sprintf('model_fitness_allCells.png');
saveas(modelfits_plot, fullfile(inDIR, plt_nm), 'png');
