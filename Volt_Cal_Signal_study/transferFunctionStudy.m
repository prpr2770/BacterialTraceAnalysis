% Studying the transfer function relationship between the Calcium and
% Voltage signals.


% Procedure:
% 1. Read each dataset pair of volt-ca signals.
% 2. create the iddata object
% 3. obtain transfer function relationship containing one/two poles.
% (Q: How does the number of poles of the system matter? What about the number of zeros?)
% alternate the input and output relationship: volt-ca, ca-volt
% 4. What do you understand from linear systems analysis?

% -------------------------------------------------------------------------

% read all the signals.
close all;

% %{
% clear all;
% 
% inDIR = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\ecoli';
% file_nm = 'ecoli_volt_ca_sigs.mat';
% load([inDIR filesep file_nm]);





%%

outDIR = [inDIR filesep 'transferFunctionStudy'];
if ~exist(outDIR)
    mkdir(outDIR)
end

caDIR = [outDIR filesep 'inputCA\'];
voltDIR = [outDIR filesep 'inputVOLT\'];

if ~exist(caDIR) || ~exist(voltDIR)
    mkdir(caDIR)
    mkdir(voltDIR)
end

% %}


Fs = 5;
Ts = 1/Fs;



% ensure size of the variable-structures are equal.
if length(volt_sigs_noisy) ~= length(ca_sigs_noisy)
    error('Unequal dataset lengths.')
    
else
    
    model_est_pvecs = [];
    model_est_dvecs = [];
    model_est_fit = [];
    
    for idx = 1:length(volt_sigs_noisy)
        sprintf('Processing >> (%d /%d) from (%s, %s)', idx, length(volt_sigs_noisy), cellFamily, inputType )
        volt_sig = volt_sigs_noisy(idx).sig_dn_wav';
        ca_sig = ca_sigs_noisy(idx).sig_dn_wav';
        
        
        
        switch inputType
            
            case 'ca'
                saveDIR = caDIR;
                u = ca_sig - mean(ca_sig);
                y = volt_sig - mean(volt_sig);
                
                num_poles = 6;
                num_zeros = 6;
                iodelay = NaN; % unknown transport delay
            case 'volt'
                saveDIR = voltDIR;
                u = volt_sig - mean(volt_sig);
                y = ca_sig - mean(ca_sig);
                
                num_poles = 5;
                num_zeros = 3;
                iodelay = NaN; % unknown transport delay
        end
        
        
        %% Find the correlation and the delay, and plot the magnitudes
        
        %{
        [corr_sig, lag] = xcorr(y,u);
        [mVal, mInd] = max(abs(corr_sig));
        tau = lag(mInd);
        delay = abs(tau);
        %         tau = ceil(length(corr_sig)/2 - mInd);

        % extract delayed signal
        if tau < 0 % y leads u
            y_tau = y(1:end-delay);
            u_tau = u(1+delay:end);
        else % y lags u
            u_tau = u(1:end-delay);
            y_tau = y(1+delay:end);
        end
        
        
        
        
        fig = figure;
        % plot the delayed signals
        subplot(2,2,1)
        p1 = plot(u,'g'); hold on;
        p2 = plot(y,'r'); hold off;
        ylabel('msrmt')
        legend([p1 p2],'inp','oup')
        subplot(2,2,2)
        plot(u_tau*norm(y_tau)/norm(u_tau),'g'); hold on; %
        plot(y_tau,'r'); hold off;
        ylabel(sprintf('delay: %d',tau))
        subplot(2,2,3)
        plot(corr_sig)
        ylabel('xcorr')
        subplot(2,2,4)
        plot(y_tau - (u_tau*norm(y_tau)/norm(u_tau)))
        ylabel('diff')
        
        
        plt_nm = sprintf('delayed_io_%d.png',idx);
        saveas(fig, fullfile(saveDIR, plt_nm), 'png');
        
        %}
        
        %%
        
        %%{
        
        % extract data-subsets for evaluation
        data = iddata(y,u,Ts);
        
        %         data_len = length(data);
        %         data_est = data(1:ceil(0.7*data_len));
        %         data_valid = data(ceil(0.5*data_len): end);
        %%
        % assign the data-subsets for tasks.
        estimation = data;
        validation = data;  % determine validation dataset.
        
        %%
        % estimated - model parameters.
        
        model_est = tfest(estimation,num_poles, num_zeros, iodelay,'Ts',data.Ts);
        %         model_est = tfest(data,num_poles,'Ts',data.Ts);
        [~,fit_goodness,~] = compare(estimation, model_est);
        [pvec, dpvec] =  getpvec(model_est);
        model_est_pvecs = [model_est_pvecs pvec];
        model_est_dvecs = [model_est_dvecs dpvec];
        model_est_fit = [model_est_fit fit_goodness];
        
        
        % obtain the bode-plot and nyquist plot
        model_idfrd = idfrd(model_est);
        
        fig_freqresp = figure; 
        subplot(1,2,1); bode(model_idfrd); 
        subplot(1,2,2); nyquist(model_idfrd);
        
        plt_nm = sprintf('freqResponse_model_%d.png',idx);
        saveas(fig_freqresp, fullfile(saveDIR, plt_nm), 'png');
        
        %%
        
        %{
        % Test the model
        I = 1:length(validation);
        opt = compareOptions('Samples',I);
        
        fig1 = figure;
        %         compare(validation,model_est, opt);
        compare(validation,model_est);
        
        pause
        plt_nm = sprintf('compare_io_model_%d.png',idx);
        %         saveas(fig1, fullfile(saveDIR, plt_nm), 'png');
        
        
        %}
        %%
        
        
        
        %{
        %% alternative method

        
        fig2 = figure(2)
        plot(data)
        
        z = [-0.5+1i, -0.5-1i, -1];
        p = [-1.11+2i, -1.11-2i, -3.01, -4.01, -0.02];
        k = 10.1;
        parameters = {z,p,k};
        Ts = 0;
        odefun = @zpkestODE;
        init_sys = idgrey(odefun,parameters,'cd',{},Ts);
        
        fig3 = figure(3);
        compare(data, init_sys)
        
        opt = greyestOptions('InitialState','zero','DisturbanceModel','none','SearchMethod','gna');
        model2 = greyest(data,init_sys,opt);

        fig4 = figure(4)
        compare(data,init_sys,model2, model_est)
        
        fig5 = figure(5)
        compare(data, model_est)
        
        %}
        
        %%
        close all;
    end
    
    fig = figure;
    subplot(2,1,1)
    plot(model_est_pvecs')
    title('model-parameters')
    subplot(2,1,2)
    plot(model_est_dvecs')
    title('parameter-deviation')
    
        plt_nm = sprintf('model_paramters.png');
        saveas(fig, fullfile(saveDIR, plt_nm), 'png');
        
file_nm = [inputType '_input_model_parameters.mat'];
save(fullfile(saveDIR, file_nm),'model_est_dvecs','model_est_pvecs','model_est_fit');
    

%{
%% tSNE - visualization

% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 8;

figure;
DATA = tsne(model_est_pvecs', [], out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params');

% scatter plot for all the data
fig = figure;
colormap jet
h = scatter3(DATA(:,1),DATA(:,2),DATA(:,3),'filled');
title(PLOT_TITLE)

figName = [saveDIR filesep 'model_est_params.fig'];
savefig(figName);

plt_nm = sprintf('model_est_params.png');
saveas(fig, fullfile(saveDIR, plt_nm), 'png');

%}
end