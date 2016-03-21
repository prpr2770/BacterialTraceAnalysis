function find_motif()   % (c) Abdullah Al Mueen and Eamonn Keogh 2008

    % get details about the file. 
    [FileName,PathName] = uigetfile('*.txt','Select the data file');
        FileNameMatlab = strcat(PathName,FileName);
        raw_time_series = load(FileNameMatlab); 
    FileName = strcat('"',PathName,FileName,'"');
    
    
    % get details about the time-series
    prompt = {'m : Length of the Timeseries','s : Minimum Length of the Motif','t : maximum Length of the Motif'};
    dlg_title = 'Input for MK';
    num_lines = 1;
    def = {'10000','128','256'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    c = sprintf('!moen.exe %s %s %s %s > out.txt',FileName,char(answer{1}),char(answer{2}),char(answer{3}));
    result = evalc(c);
    disp(result);
    
    
    
    %%% At this point, the algorithm is finished, the data has been written to disk
    %%% Ask the user if she wants to see the output in a pretty plot.  
 
    L = csvread('out.txt');

    plotMotifs(raw_time_series,L(:,2),L(:,3),L(:,1));
  






    
    
    
    
    
    
    
    
    
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All the code below here is just used for ploting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





