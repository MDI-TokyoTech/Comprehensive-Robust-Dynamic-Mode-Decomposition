function [X, h, v] = load_cylinder_data(rootDir, dataName)
    % Load data specific to the project structure
    % Inputs:
    %   rootDir: Root directory of the project
    %   dataName: Name of the dataset (e.g., 'cylinder')
    
    if nargin < 2
        dataName = 'cylinder'; % Default
    end

    switch lower(dataName)
        case 'cylinder'
            dataFilePath = fullfile(rootDir, 'data', 'cylinder', 'cylinder.mat');
            % Data Dimensions hardcoded for this dataset
            h = 449;
            v = 199;
        
        % Add other datasets here
        % case 'other_data'
        %     dataFilePath = fullfile(rootDir, 'data', 'other', 'data.mat');
        %     h = ...; v = ...;

        otherwise
            error('Unknown dataset name: %s', dataName);
    end

    if ~isfile(dataFilePath)
        error('Data file not found at: %s', dataFilePath);
    end

    fprintf('Loading data (%s) from %s...\n', dataName, dataFilePath);
    loadedData = load(dataFilePath);
    varNames = fieldnames(loadedData);
    X = loadedData.(varNames{1}); % Assume the first variable is the data

    [N, M] = size(X);
    fprintf('Data Matrix Size: %d x %d (State x Time)\n', N, M);
end