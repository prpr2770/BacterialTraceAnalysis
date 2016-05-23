%%
% EXTRACT FEATURESPACE MODEL-PARAMETERS AND TSNE-VISUALIZE THEM.


% ABSENCE OF DRUGS:

close all;

numVideos = size(data,2);
if numVideos == 2
    
    
    MODEL_DATA_BEFORE_DRUG = [];
    MODEL_DATA_AFTER_DRUG = [];
    
    for expID = 1:size(data,1)
        expID
        
        % Creating directory for given experiment
        expName = sprintf('expID_%d',expID);
        expDIR = [inDIR filesep expName];
        if ~exist(expDIR)
            mkdir(expDIR)
        end
        
        colorCode = expID;
        
        switch colorCode
            case 0
                colorChoice = [0 0 0];
            case 1
                colorChoice = [0 0 1];
            case 2
                colorChoice = [0 1 0];
            case 3
                colorChoice = [0 1 1];
            case 4
                colorChoice = [1 0 0];
            case 5
                colorChoice = [1 0 1];
            case 6
                colorChoice = [1 1 0];
            case 7
                colorChoice = [1 1 1];
            otherwise
                warning('ColorCode out of range!')
        end
        
        
        
        
        param_data_before_drug = featurespace(expID,1).model_params;
        fitness_data_before_drug = featurespace(expID,1).model_fitness;
        
        param_data_after_drug = featurespace(expID,3).model_params;
        fitness_data_after_drug = featurespace(expID,3).model_fitness;
        
        
        %% Ensuring that datasets exist for atleast 2 videos: Before Drug and After Drug
        if length(param_data_before_drug) ~= 0 && length(param_data_before_drug) ~=0
            
            model_fit_threshold = 50;
            
            % determine the index of cells satsifying threshold.
            good_params_before_drug_indx = find(fitness_data_before_drug > model_fit_threshold);
            num_good_params_before = length(good_params_before_drug_indx);
            
            
            good_params_after_drug_indx = find(fitness_data_after_drug > model_fit_threshold);
            num_good_params_after = length(good_params_after_drug_indx);
            
            % find intersection of the indices common to before and after drug.
            good_indx = intersect(good_params_before_drug_indx, good_params_after_drug_indx);
            
            if length(good_indx) > 0
                % extract the good_cell_models ensuring same_cellID in both before and after datasets.
                color_data_before_drug = colorCode*ones(1,length(good_indx));
                color_data_after_drug = colorCode*ones(1,length(good_indx));
                
                good_param_data_before_drug = param_data_before_drug(:,good_indx);
                good_param_data_after_drug = param_data_after_drug(:,good_indx);
                
                
                % -----------------------------------------------------------------------------
                % CREATE DATA-STRUCTURE TO COMPILE ALL THE DATASETS BASED ON
                % BEFORE-AFTER DRUG USE
                
                MODEL_DATA_BEFORE_DRUG = [MODEL_DATA_BEFORE_DRUG good_param_data_before_drug];
                MODEL_DATA_AFTER_DRUG = [MODEL_DATA_AFTER_DRUG good_param_data_after_drug];
                
                % -----------------------------------------------------------------------------
                % FOR EACH EXPERIMENT, COMPARE VISUALIZATION OF THE BEFORE AND AFTER
                % DRUG USE DATASET.
                % Implement tSNE on the dataset.
                
                
                % implement tSNE
                perplexity = 30;
                out_dims = 3;
                initial_dims = 8;
                
                figure;
                EXP_DATA = [good_param_data_before_drug good_param_data_after_drug]'; % row-vectors
                
                size(EXP_DATA)
                size(param_data_before_drug)
                size(good_param_data_before_drug)
                size(param_data_after_drug)
                size(good_param_data_after_drug)
                
                REDX_DATA = tsne(EXP_DATA, [color_data_before_drug 0.5*color_data_after_drug]', out_dims, initial_dims, perplexity);
                PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params for ExpID: %d', expID);
                
                % scatter plot for all the data
                fig = figure(expID);
                colormap jet
                h1 = scatter3(REDX_DATA(1:length(color_data_before_drug),1),REDX_DATA(1:length(color_data_before_drug),2),REDX_DATA(1:length(color_data_before_drug),3),'MarkerFaceColor',colorChoice);
                hold on
                h2 = scatter3(REDX_DATA(1+length(color_data_before_drug):end,1),REDX_DATA(1+length(color_data_before_drug):end,2),REDX_DATA(1+length(color_data_before_drug):end,3),'MarkerFaceColor', 0.5*colorChoice);
                hold off
                title(PLOT_TITLE)
                legend([h1,h2],'before','after')
                
                figName = [expDIR filesep 'model_est_params.fig'];
                savefig(figName);
                
                plt_nm = sprintf('model_est_params.png');
                saveas(fig, fullfile(expDIR, plt_nm), 'png');
            else
                warning('No common good_model_parameter cells found for Before/After Tests in ExpID: %d.', expID);
            end
        end
        
    end
    
    
    % -----------------------------------------------------------------------------
    % FOR ALL EXPERIMENT, COMPARE VISUALIZATION OF THE BEFORE AND AFTER
    % DRUG USE DATASET.
    
    % implement tSNE
    perplexity = 30;
    out_dims = 3;
    initial_dims = 8;
    
    figure;
    EXP_DATA = [MODEL_DATA_BEFORE_DRUG MODEL_DATA_AFTER_DRUG]'; % row-vectors
    CLUSTER_DATA = [1*ones(1,size(MODEL_DATA_BEFORE_DRUG,2)) 2*ones(1,size(MODEL_DATA_AFTER_DRUG,2))];
    REDX_DATA = tsne(EXP_DATA, CLUSTER_DATA', out_dims, initial_dims, perplexity);
    PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params for all experiments: Before Drug');
    
    len_before_data = size(MODEL_DATA_BEFORE_DRUG,2);
    % scatter plot for all the data
    fig = figure;
    colormap jet
    h1 = scatter3(REDX_DATA(1:len_before_data,1),REDX_DATA(1:len_before_data,2),REDX_DATA(1:len_before_data,3),'MarkerFaceColor',colorChoice);
    hold on
    h2 = scatter3(REDX_DATA(1+len_before_data:end,1),REDX_DATA(1+len_before_data:end,2),REDX_DATA(1+len_before_data:end,3),'MarkerFaceColor', 0.5*colorChoice);
    hold off
    title(PLOT_TITLE)
    legend([h1,h2],'before','after')
    
    figName = [inDIR filesep 'allexp_model_est_params.fig'];
    savefig(figName);
    
    plt_nm = sprintf('allexp_model_est_params.png');
    saveas(fig, fullfile(inDIR, plt_nm), 'png');
    
    
elseif numVideos ==3
    %% Clustering/ Data Visualization with 3 videos:
    
    
    MODEL_DATA_BEFORE_DRUG = [];
    MODEL_DATA_DURING_DRUG = [];
    MODEL_DATA_AFTER_DRUG = [];
    
    for expID = 1:size(data,1)
        expID
        
        % Creating directory for given experiment
        expName = sprintf('expID_%d',expID);
        expDIR = [inDIR filesep expName];
        if ~exist(expDIR)
            mkdir(expDIR)
        end
        
        colorCode = expID;
        
        switch colorCode
            case 0
                colorChoice = [0 0 0];
            case 1
                colorChoice = [0 0 1];
            case 2
                colorChoice = [0 1 0];
            case 3
                colorChoice = [0 1 1];
            case 4
                colorChoice = [1 0 0];
            case 5
                colorChoice = [1 0 1];
            case 6
                colorChoice = [1 1 0];
            case 7
                colorChoice = [1 1 1];
            otherwise
                warning('ColorCode out of range!')
        end
        
        
        
        
        param_data_before_drug = featurespace(expID,1).model_params;
        fitness_data_before_drug = featurespace(expID,1).model_fitness;
        
        param_data_during_drug = featurespace(expID,2).model_params;
        fitness_data_during_drug = featurespace(expID,2).model_fitness;
        
        
        param_data_after_drug = featurespace(expID,3).model_params;
        fitness_data_after_drug = featurespace(expID,3).model_fitness;
        
        
        %% Ensuring that datasets exist for atleast 2 videos: Before Drug and After Drug
        if length(param_data_before_drug) ~= 0 && length(param_data_before_drug) ~=0 && length(param_data_during_drug) ~=0
            
            model_fit_threshold = 45;
            
            % determine the index of cells satsifying threshold.
            good_params_before_drug_indx = find(fitness_data_before_drug > model_fit_threshold);
            num_good_params_before = length(good_params_before_drug_indx);
            
            good_params_during_drug_indx = find(fitness_data_during_drug > model_fit_threshold);
            num_good_params_during = length(good_params_during_drug_indx);
            
            good_params_after_drug_indx = find(fitness_data_after_drug > model_fit_threshold);
            num_good_params_after = length(good_params_after_drug_indx);
            
            % find intersection of the indices common to before and after drug.
            good_indx = intersect(good_params_before_drug_indx, good_params_after_drug_indx);
            good_indx = intersect(good_indx, good_params_during_drug_indx)
            
            if length(good_indx) > 0
                % extract the good_cell_models ensuring same_cellID in both before and after datasets.
                color_data_before_drug = colorCode*ones(1,length(good_indx));
                color_data_during_drug = colorCode*ones(1,length(good_indx));
                color_data_after_drug = colorCode*ones(1,length(good_indx));
                
                good_param_data_before_drug = param_data_before_drug(:,good_indx);
                good_param_data_during_drug = param_data_during_drug(:,good_indx);
                good_param_data_after_drug = param_data_after_drug(:,good_indx);
                
                
                % -----------------------------------------------------------------------------
                % CREATE DATA-STRUCTURE TO COMPILE ALL THE DATASETS BASED ON
                % BEFORE-AFTER DRUG USE
                
                MODEL_DATA_BEFORE_DRUG = [MODEL_DATA_BEFORE_DRUG good_param_data_before_drug];
                MODEL_DATA_DURING_DRUG = [MODEL_DATA_DURING_DRUG good_param_data_during_drug];
                MODEL_DATA_AFTER_DRUG = [MODEL_DATA_AFTER_DRUG good_param_data_after_drug];
                
                % -----------------------------------------------------------------------------
                % FOR EACH EXPERIMENT, COMPARE VISUALIZATION OF THE BEFORE AND AFTER
                % DRUG USE DATASET.
                % Implement tSNE on the dataset.
                
                
                % implement tSNE
                perplexity = 30;
                out_dims = 3;
                initial_dims = 8;
                
                figure;
                EXP_DATA = [good_param_data_before_drug good_param_data_during_drug good_param_data_after_drug]'; % row-vectors
                
                %                 size(EXP_DATA)
                %                 size(param_data_before_drug)
                %                 size(good_param_data_before_drug)
                %                 size(param_data_after_drug)
                %                 size(good_param_data_after_drug)
                
                REDX_DATA = tsne(EXP_DATA, [color_data_before_drug 0.66*color_data_during_drug 0.33*color_data_after_drug]', out_dims, initial_dims, perplexity);
                PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params for ExpID: %d', expID);
                
                % scatter plot for all the data
                fig = figure(expID);
                colormap jet
                h1 = scatter3(REDX_DATA(1:length(color_data_before_drug),1),REDX_DATA(1:length(color_data_before_drug),2),REDX_DATA(1:length(color_data_before_drug),3),'MarkerFaceColor',colorChoice);
                hold on
                h2 = scatter3(REDX_DATA(1+length(color_data_before_drug):length([color_data_before_drug color_data_during_drug]),1),REDX_DATA(1+length(color_data_before_drug):length([color_data_before_drug color_data_during_drug]),2),REDX_DATA(1+length(color_data_before_drug):length([color_data_before_drug color_data_during_drug]),3),'MarkerFaceColor', 0.6*colorChoice);
                hold on
                h3 = scatter3(REDX_DATA(1+length([color_data_before_drug color_data_during_drug]):end,1),REDX_DATA(1+length([color_data_before_drug color_data_during_drug]):end,2),REDX_DATA(1+length([color_data_before_drug color_data_during_drug]):end,3),'MarkerFaceColor', 0.3*colorChoice);
                hold off
                title(PLOT_TITLE)
                legend([h1,h2,h3],'before','during','after');
                
                figName = [expDIR filesep 'model_est_params.fig'];
                savefig(figName);
                
                plt_nm = sprintf('model_est_params.png');
                saveas(fig, fullfile(expDIR, plt_nm), 'png');
                
                % histogram of distances travelled
                before_redx_data = REDX_DATA(1:length(color_data_before_drug), :);
                during_redx_data = REDX_DATA(1+length(color_data_before_drug):length([color_data_before_drug color_data_during_drug]),:);
                after_redx_data = REDX_DATA(1+length([color_data_before_drug color_data_during_drug]):end,:);
                
                before_during_dist = sqrt(sum((before_redx_data - during_redx_data).^2,2));
                before_after_dist = sqrt(sum((before_redx_data - after_redx_data).^2,2));

                hist1 = histogram(before_during_dist,16);
                ylim([0 15]);
                plt_nm = sprintf('hist_before_during_redx_dist_expID_%d.png',expID);
                PLOT_TITLE = sprintf('histogram : before/during drug [ExpID: %d]',expID);
                title(PLOT_TITLE)
                saveas(hist1, fullfile(expDIR, plt_nm), 'png');

                hist2 = histogram(before_after_dist,16);
                ylim([0 15]);
                plt_nm = sprintf('hist_before_after_redx_dist_expID_%d.png',expID);
                PLOT_TITLE = sprintf('histogram : before/after drug [ExpID: %d]',expID);
                title(PLOT_TITLE)
                saveas(hist2, fullfile(expDIR, plt_nm), 'png');
                
                
            else
                warning('No common good_model_parameter cells found for Before/After Tests in ExpID: %d.', expID);
            end
        end
        
    end
    
    
    % -----------------------------------------------------------------------------
    % FOR ALL EXPERIMENT, COMPARE VISUALIZATION OF THE BEFORE AND AFTER
    % DRUG USE DATASET.
    
    % implement tSNE
    perplexity = 30;
    out_dims = 3;
    initial_dims = 8;
    
    figure;
    EXP_DATA = [MODEL_DATA_BEFORE_DRUG MODEL_DATA_DURING_DRUG MODEL_DATA_AFTER_DRUG]'; % row-vectors
    CLUSTER_DATA = [1*ones(1,size(MODEL_DATA_BEFORE_DRUG,2)) 2*ones(1,size(MODEL_DATA_DURING_DRUG,2)) 3*ones(1,size(MODEL_DATA_AFTER_DRUG,2))];
    REDX_DATA = tsne(EXP_DATA, CLUSTER_DATA', out_dims, initial_dims, perplexity);
    PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params for all experiments');
    
    len_before_data = size(MODEL_DATA_BEFORE_DRUG,2);
    len_before_during_data = size([MODEL_DATA_BEFORE_DRUG MODEL_DATA_DURING_DRUG],2);
    % scatter plot for all the data
    fig = figure;
    colormap jet
    h1 = scatter3(REDX_DATA(1:len_before_data,1),REDX_DATA(1:len_before_data,2),REDX_DATA(1:len_before_data,3),'MarkerFaceColor',colorChoice);
    hold on
    h2 = scatter3(REDX_DATA(1+len_before_data:len_before_during_data,1),REDX_DATA(1+len_before_data:len_before_during_data,2),REDX_DATA(1+len_before_data:len_before_during_data,3),'MarkerFaceColor', 0.66*colorChoice);
    hold on
    h3 = scatter3(REDX_DATA(1+len_before_during_data:end,1),REDX_DATA(1+len_before_during_data:end,2),REDX_DATA(1+len_before_during_data:end,3),'MarkerFaceColor', 0.33*colorChoice);
    hold off
    title(PLOT_TITLE)
    legend([h1,h2,h3],'before','during','after')
    
    figName = [inDIR filesep 'allexp_model_est_params.fig'];
    savefig(figName);
    
    plt_nm = sprintf('allexp_model_est_params.png');
    saveas(fig, fullfile(inDIR, plt_nm), 'png');
    
    % histograms
    before_data = REDX_DATA(1:len_before_data,:);
    during_data = REDX_DATA(1+len_before_data:len_before_during_data,:)
    after_data = REDX_DATA(1+len_before_during_data:end,:);
    
    before_during_dist = sqrt(sum((before_data - during_data).^2,2));
    before_after_dist = sqrt(sum((before_data - after_data).^2,2));
    
    hist1 = histogram(before_during_dist,16);
    ylim([0 30]);
    plt_nm = sprintf('hist_before_during_redx_dist_allData.png');
    PLOT_TITLE = sprintf('histogram : before/during drug [all experiments]');
    title(PLOT_TITLE)
    saveas(hist1, fullfile(inDIR, plt_nm), 'png');
    
    hist2 = histogram(before_after_dist,16);
    ylim([0 30]);
    plt_nm = sprintf('hist_before_after_redx_dist_allData.png');
    PLOT_TITLE = sprintf('histogram : before/after drug [all experiments]');
    title(PLOT_TITLE)
    saveas(hist2, fullfile(inDIR, plt_nm), 'png');
    

    %% SVM visualization of the data. 
    EXP_DATA = [MODEL_DATA_BEFORE_DRUG MODEL_DATA_AFTER_DRUG]'; % row-vectors
    CLUSTER_DATA = [-1*ones(1,size(MODEL_DATA_BEFORE_DRUG,2)) 1*ones(1,size(MODEL_DATA_AFTER_DRUG,2))]';
    
    runSVM_Classifier(EXP_DATA, CLUSTER_DATA);
    
    

end