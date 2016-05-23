function runSVM_Classifier(data, class_data)
%Train the SVM Classifier

highd_svm_model = fitcsvm(data,class_data,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[-1,1]);
highd_svm_model.KernelParameters

%{

%% SVM over 8-dim space.  

%Predict scores over the grid

min_data = min(data);
max_data = max(data);
d = min(max_data, max((max_data - min_data)/10)); % Grid Width

columnGridSquares = zeros(1,size(data,2));

% determine GridMarkers for each dimension
for i = 1:size(data,2) % determine dimensions of data
    % construct grid for each dimension
    eval(sprintf('x%d_Coords = min(data(:,i)):d:max(data(:,i));',i));
    eval(sprintf('columnGridSquares(%d)= length(x%d_Coords);',i,i));
end

% determine GridMarkers based GRID for all dimensions
for dim_data = 1:size(data,2)
    repopts = '';
    for j = 1:length(columnGridSquares)
        
        if j~=dim_data
            repopts = [repopts ', ' num2str(columnGridSquares(j))];
        else
            repopts = [repopts ', 1'];
        end
    end
    
    
    % CLARIFY??? confirm ??
    eval( sprintf( 'x%dGrid = repmat( x%d_Coords %s ) ;', dim_data, dim_data, repopts ) );
end

% Combine all the GridMarker based Grids into a single GRID BLOCK of
% dimensions as the original data.
xGrid = [];
for dim = 1:size(data,2)
    % modify this for a generic sentence.
    eval(sprintf('xGrid = [xGrid x%dGrid(:)];'));
end

% obtain the predictions using SVM for the GridPoints.
[labels,scores] = predict(cl,xGrid); % [labels, Score] = predict()


%}


%%
% Now obtain the DIMENSION-REDUCTION based representation for the above
% high-dimensional data, using tSNE.


% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 8;

figure;
EXP_DATA = data; % row-vectors
CLUSTER_DATA = class_data;
REDX_DATA = tsne(EXP_DATA, CLUSTER_DATA, out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of Model-Est Params for all experiments');

%% SVM over 2-dim space

%Train the SVM Classifier
red_svm_model = fitcsvm(REDX_DATA,CLUSTER_DATA ,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[-1,1]);

% Predict scores over the grid
min_data = min(REDX_DATA);
max_data = max(REDX_DATA);
d = min(max_data, max((max_data - min_data)/10)); % Grid Width

% d = 0.02;
x1_Coords = min(REDX_DATA(:,1)):d:max(REDX_DATA(:,1));
x2_Coords = min(REDX_DATA(:,2)):d:max(REDX_DATA(:,2));
x3_Coords = min(REDX_DATA(:,3)):d:max(REDX_DATA(:,3));
[x1Grid,x2Grid, x3Grid ] = meshgrid(x1_Coords,x2_Coords,x3_Coords);

xGrid = [x1Grid(:),x2Grid(:),x3Grid(:)];
[~,scores] = predict(red_svm_model,xGrid);


% Plot the data and the decision boundary
figure;
h(1:2) = gscatter(REDX_DATA(:,1),REDX_DATA(:,2),class_data,'rb','.');
hold on
h(3) = plot(REDX_DATA(red_svm_model.IsSupportVector,1),REDX_DATA(red_svm_model.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'-1','+1','Support Vectors'});
axis equal
hold off







end
