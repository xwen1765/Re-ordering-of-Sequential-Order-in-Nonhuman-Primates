%% License
% Dataset License: CC BY-NC-ND 4.0 + Custom Restriction
% This dataset is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0).
% You are free to:
% Share — copy and redistribute the material in any medium or format Under the following terms:
% 
% Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. 
% NonCommercial — You may not use the material for commercial purposes. NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material. 
% Additional Restrictions: Academic use is not permitted without prior written permission. This includes any use for research projects, academic publications, presentations, or classroom instruction. 
% If you wish to use this dataset and code for academic or commercial purposes, please contact: Xuan Wen Vanderbilt University xuan.wen@vanderbilt.edu
% 
% Data for data analysis is avaliable at https://osf.io/2vk97/, with DOI 10.17605/OSF.IO/2VK97.


%% Data Loading
% touch_table = readtable();
% block_table = readtable();
% WM_WWW_Interaction = readtable();


%% SUP
% Figure S1A
% accuracy_matrices = generateTouchAccuracyPlots(touch_table);
% plotLastFourTrialsAccuracy(accuracy_matrices, touch_table)
% Figure S1B
% plotProportionTrialsByState(touch_table);
% Figure S1C
% plotTrialTo80PercentByState(touch_table);
% Figure S1D
% plotNewErrorRates(touch_table);
% Figure S1E
% plotErrorDropRates(touch_table);


% Figure S2A
% calculateLearningPerformanceSwap(touch_table);
% Figure S2B
% plotLearningPerformanceSwapBoxPlot(touch_table);
% Figure S2C
% plotAvgSearchTimeByOrdinalPositionLearning(touch_table);
% Figure S2D
% plotAvgSearchTimeByStateAndSwap(touch_table);

% Figure S3A
% plotFirstTrialPosition80Repeat(touch_table);
% Figure S3B
% plotAvgSearchTimeCorrectTouchesRepeatConditions(touch_table, block_table)

% Figure S4A
% plotWM_Performance(WM_WWW_Interaction)
% Figure S4B
% plotWM_Performance_Individual(WM_WWW_Interaction)


% Figure S5A
% plotErrorDecreasingRateSwap(touch_table);
% Figure S5B
% plotErrorDecreasingRateMemory(touch_table);
% Figure S5C
% plotDistractorErrorFrequency(touch_table);




function plotWM_Performance_Individual(WM_WWW_Interaction)
    % plotWM_Performance_Individual plots the average accuracy for each subject
    % across different delay conditions and similarity conditions.
    %
    %   plotWM_Performance_Individual(WM_WWW_Interaction) takes a table
    %   WM_WWW_Interaction as input. Each row of the table represents a
    %   session's data. The table must contain the following columns:
    %   - 'Subject': A column identifying each subject (e.g., double, string, categorical).
    %   - 'Delay': A cell array where each cell contains a 1xN double array
    %              of delay times for that session (e.g., [0.5, 1.3, 1.8]).
    %              N can be different for each session/row.
    %   - 'Condition': A cell array where each cell contains a 1xN double
    %                  array of condition codes for that session (0 for Low
    %                  Similarity, 2 for High Similarity).
    %                  N can be different for each session/row.
    %   - 'Accuracy': A cell array where each cell contains a 1xN double
    %                 array of accuracy values (0 or 1) for that session.
    %                 N can be different for each session/row.
    %
    %   The function generates a plot with:
    %   - X-axis: Equally spaced numerical points (1, 2, 3) corresponding to
    %             delay conditions (labeled as 0.5, 1.25, 1.75 seconds).
    %   - Y-axis: Average Accuracy for each subject.
    %   - Lines: Separate lines for each subject. Dashed lines for Low
    %            Similarity (Condition 0) and solid lines for High Similarity
    %            (Condition 2).
    %   - Figure size set to 300x400 pixels.

    % 1. Input Validation
    if ~istable(WM_WWW_Interaction)
        error('plotWM_Performance_Individual:InvalidInput', 'Input must be a MATLAB table.');
    end

    requiredCols = {'Subject', 'Delay', 'Condition', 'Accuracy'};
    for i = 1:length(requiredCols)
        if ~ismember(requiredCols{i}, WM_WWW_Interaction.Properties.VariableNames)
            error('plotWM_Performance_Individual:MissingColumn', ...
                  'Table must contain a column named ''%s''.', requiredCols{i});
        end
    end

    % Validate cell array columns
    cellCols = {'Delay', 'Condition', 'Accuracy'};
     for i = 1:length(cellCols)
        if ~iscell(WM_WWW_Interaction.(cellCols{i}))
            error('plotWM_Performance_Individual:ColumnTypeMismatch', ...
                  'Column ''%s'' must be a cell array, expecting 1xN double arrays within cells.', cellCols{i});
        end
    end


    % Get unique subject IDs
    % If 'Subject' is a cell array, unique will return a cell array of unique IDs
    subjectIDs = unique(WM_WWW_Interaction.Subject);
    numSubjects = length(subjectIDs);

    if numSubjects == 0
        warning('plotWM_Performance_Individual:NoSubjects', 'No subjects found in the data.');
        return; % Exit if no subjects
    end

    % Flatten all delays to get unique delay values across all data
    allDelaysFlat = [WM_WWW_Interaction.Delay{:}];
    dataDelays = unique(allDelaysFlat);
    dataDelays = sort(dataDelays); % Ensures processing order is consistent

    if length(dataDelays) ~= 3
         warning('plotWM_Performance_Individual:UnexpectedDelays', ...
                 'Expected 3 unique delay values (0.5, 1.3, 1.8), but found %d.', length(dataDelays));
         % Proceed assuming the found delays are the ones to plot, but labels might be off
    end


    % Define the condition codes
    lowSimCondition = 0;
    highSimCondition = 2;

    % Define x-coordinates for plotting as equally spaced points (1, 2, 3, ...)
    plotX = 1:length(dataDelays);

    % Define x-axis labels as specified (0.5, 1.25, 1.75)
    xTickLabels = {'0.5', '1.25', '1.75'};

    % Create a new figure window and set its position (size)
    figure('Name', 'Individual Subject Working Memory Performance', ...
           'Position', [100, 100, 300, 400]); % [left, bottom, width, height]

    hold on; % Allow multiple plots on the same axes

    % Define a set of distinct colors for subjects
    colors = lines(numSubjects); % Use built-in 'lines' colormap for distinct colors

    % Initialize arrays to store legend handles
    % We need two handles per subject for the two line styles
    legendHandles = gobjects(numSubjects * 2, 1);
    legendLabels = cell(numSubjects * 2, 1);

    lineIndex = 1;

    % Iterate through each subject
    for s = 1:numSubjects
        currentSubjectID = subjectIDs{s}; % Access subject ID from cell array

        % Fix for the 'Operator '==' is not supported for operands of type 'cell'.' error
        % Use cellfun with isequal to compare elements within the 'Subject' cell array
        idx_subject = cellfun(@(x) isequal(x, currentSubjectID), WM_WWW_Interaction.Subject);
        subjectData = WM_WWW_Interaction(idx_subject, :);

        % Flatten data for the current subject
        subjectDelays = [subjectData.Delay{:}];
        subjectConditions = [subjectData.Condition{:}];
        subjectAccuracies = [subjectData.Accuracy{:}];

        % Initialize arrays to store average accuracies for current subject
        avgAccuracy_LowSim_Sub = zeros(size(dataDelays));
        avgAccuracy_HighSim_Sub = zeros(size(dataDelays));

        % Calculate average accuracy for each delay and condition for the current subject
        for i = 1:length(dataDelays)
            currentDelay = dataDelays(i);

            % Filter data for the current delay and subject
            idxForCurrentDelay = (subjectDelays == currentDelay);

            % Calculate average accuracy for Low Similarity (Condition 0)
            idxLowSim = idxForCurrentDelay & (subjectConditions == lowSimCondition);
            if any(idxLowSim)
                avgAccuracy_LowSim_Sub(i) = mean(subjectAccuracies(idxLowSim));
            else
                avgAccuracy_LowSim_Sub(i) = NaN; % Assign NaN if no data for this group
            end

            % Calculate average accuracy for High Similarity (Condition 2)
            idxHighSim = idxForCurrentDelay & (subjectConditions == highSimCondition);
            if any(idxHighSim)
                avgAccuracy_HighSim_Sub(i) = mean(subjectAccuracies(idxHighSim));
            else
                avgAccuracy_HighSim_Sub(i) = NaN; % Assign NaN if no data for this group
            end
        end

        % Determine label for current subject (robust for numeric or string IDs)
        if ischar(currentSubjectID)
            subjectLabel = currentSubjectID;
        else
            subjectLabel = num2str(currentSubjectID);
        end

        % Plot Low Similarity data for the current subject (dashed line)
        h_low = plot(plotX, avgAccuracy_LowSim_Sub, '--', ...
                     'Color', colors(s,:), ...    % Use subject-specific color
                     'LineWidth', 1.5, ...
                     'Marker', 'o', ...         % Use a marker for points
                     'MarkerSize', 6);
        legendHandles(lineIndex) = h_low;
        legendLabels{lineIndex} = sprintf('Subject %s - Low Sim', subjectLabel);
        lineIndex = lineIndex + 1;

        % Plot High Similarity data for the current subject (solid line)
        h_high = plot(plotX, avgAccuracy_HighSim_Sub, '-', ...
                      'Color', colors(s,:), ...    % Use subject-specific color
                      'LineWidth', 1.5, ...
                      'Marker', 's', ...         % Use a different marker
                      'MarkerSize', 6);
        legendHandles(lineIndex) = h_high;
        legendLabels{lineIndex} = sprintf('Subject %s - High Sim', subjectLabel);
        lineIndex = lineIndex + 1;
    end

    % Remove unused legend entries (in case some subjects had missing data or less than 2 conditions)
    validLegendIndices = ~cellfun('isempty', legendLabels);
    legendHandles = legendHandles(validLegendIndices);
    legendLabels = legendLabels(validLegendIndices);

    % Set plot properties
    xlabel('Delay (seconds)'); % Label still refers to the conceptual delay
    ylabel('Average Accuracy'); % Y-axis is average accuracy *for each subject*
    title('Individual Subject Working Memory Performance');

    % Set x-axis ticks to the equally spaced plot points and apply custom labels
    xticks(plotX);         % Set tick locations to the equally spaced points
    xticklabels(xTickLabels); % Apply the custom string labels

    % Set y-axis limits and ticks
    ylim([-0.05, 1.05]); % Slightly extend limits for better visualization of 0 and 1
    yticks(0:0.1:1);     % Set y-axis ticks at 0.1 intervals

    % Add legend
    % Use the collected handles and labels
    if ~isempty(legendHandles)
        legend(legendHandles, legendLabels, 'Location', 'southwest', 'AutoUpdate', 'off');
    end


    % Remove grid (default is no grid, but explicitly turn off if it were on)
    grid off;

    hold off; % Release the hold on the current axes
end


function plotWM_Performance(WM_WWW_Interaction)
    % plotWM_Performance plots the average accuracy for a working memory task
    % across different delay conditions and similarity conditions, with SE error bars.
    %
    %   plotWM_Performance(WM_WWW_Interaction) takes a table
    %   WM_WWW_Interaction as input. Each row of the table represents a
    %   session's data. The table must contain the following columns:
    %   - 'Delay': A cell array where each cell contains a 1xN double array
    %              of delay times for that session (e.g., [0.5, 1.3, 1.8]).
    %              N can be different for each session/row.
    %   - 'Condition': A cell array where each cell contains a 1xN double
    %                  array of condition codes for that session (0 for Low
    %                  Similarity, 2 for High Similarity).
    %                  N can be different for each session/row.
    %   - 'Accuracy': A cell array where each cell contains a 1xN double
    %                 array of accuracy values (0 or 1) for that session.
    %                 N can be different for each session/row.
    %
    %   The function generates a plot with:
    %   - X-axis: Equally spaced numerical points (1, 2, 3) corresponding to
    %             delay conditions (labeled as 0.5, 1.25, 1.75 seconds).
    %   - Y-axis: Average Accuracy with Standard Error (SE) bars.
    %   - Two lines: one for Low Similarity (Condition 0) and one for High
    %     Similarity (Condition 2).
    %   - Figure size set to 300x400 pixels.

    % 1. Input Validation
    if ~istable(WM_WWW_Interaction)
        error('plotWM_Performance:InvalidInput', 'Input must be a MATLAB table.');
    end

    requiredCols = {'Delay', 'Condition', 'Accuracy'};
    for i = 1:length(requiredCols)
        if ~ismember(requiredCols{i}, WM_WWW_Interaction.Properties.VariableNames)
            error('plotWM_Performance:MissingColumn', ...
                  'Table must contain a column named ''%s''.', requiredCols{i});
        end
        % Check if the column is a cell array, which is expected for 1xN double per row
        if ~iscell(WM_WWW_Interaction.(requiredCols{i}))
            error('plotWM_Performance:ColumnTypeMismatch', ...
                  'Column ''%s'' must be a cell array, expecting 1xN double arrays within cells.', requiredCols{i});
        end
    end

    % 2. Data Flattening
    % The '{:}' operator unpacks the cell array contents into a comma-separated list.
    % Enclosing this list in square brackets `[]` performs horizontal concatenation
    % of all the 1xN arrays into a single long 1xM vector. This works even if N varies.
    allDelays = [WM_WWW_Interaction.Delay{:}];
    allConditions = [WM_WWW_Interaction.Condition{:}];
    allAccuracies = [WM_WWW_Interaction.Accuracy{:}];

    % Ensure flattened vectors are of the same length
    if ~isequal(length(allDelays), length(allConditions), length(allAccuracies))
        error('plotWM_Performance:DataMismatch', ...
              'Flattened Delay, Condition, and Accuracy vectors must have the same length.');
    end

    % Get unique delay values from the data and sort them
    dataDelays = unique(allDelays);
    dataDelays = sort(dataDelays); % Ensures processing order is consistent

    % Define the condition codes
    lowSimCondition = 0;
    highSimCondition = 2;

    % Initialize arrays to store average accuracies and standard errors for each condition across delays
    avgAccuracy_LowSim = zeros(size(dataDelays));
    seAccuracy_LowSim = zeros(size(dataDelays));
    avgAccuracy_HighSim = zeros(size(dataDelays));
    seAccuracy_HighSim = zeros(size(dataDelays));

    % 3. Data Aggregation and SE Calculation
    % Calculate average accuracy and SE for each delay and condition
    for i = 1:length(dataDelays)
        currentDelay = dataDelays(i);

        % Filter data for the current delay
        idxForCurrentDelay = (allDelays == currentDelay);

        % Calculate average accuracy and SE for Low Similarity (Condition 0)
        idxLowSim = idxForCurrentDelay & (allConditions == lowSimCondition);
        if any(idxLowSim)
            dataPoints = allAccuracies(idxLowSim);
            avgAccuracy_LowSim(i) = mean(dataPoints);
            seAccuracy_LowSim(i) = std(dataPoints) / sqrt(length(dataPoints));
        else
            avgAccuracy_LowSim(i) = NaN; % Assign NaN if no data for this group
            seAccuracy_LowSim(i) = NaN;  % Assign NaN for SE as well
        end

        % Calculate average accuracy and SE for High Similarity (Condition 2)
        idxHighSim = idxForCurrentDelay & (allConditions == highSimCondition);
        if any(idxHighSim)
            dataPoints = allAccuracies(idxHighSim);
            avgAccuracy_HighSim(i) = mean(dataPoints);
            seAccuracy_HighSim(i) = std(dataPoints) / sqrt(length(dataPoints));
        else
            avgAccuracy_HighSim(i) = NaN; % Assign NaN if no data for this group
            seAccuracy_HighSim(i) = NaN;  % Assign NaN for SE as well
        end
    end

    % Define x-coordinates for plotting as equally spaced points (1, 2, 3, ...)
    plotX = 1:length(dataDelays);

    % Define x-axis labels as specified (0.5, 1.25, 1.75)
    % These labels correspond to the dataDelays [0.5, 1.3, 1.8] plotted at plotX [1, 2, 3]
    xTickLabels = {'0.5', '1.25', '1.75'};

    % 4. Plotting
    % Create a new figure window and set its position (size)
    figure('Name', 'Working Memory Performance Across Delay Conditions', ...
           'Position', [100, 100, 300, 400]); % [left, bottom, width, height]

    hold on; % Allow multiple plots on the same axes

    % Plot Low Similarity data with error bars
    h1 = errorbar(plotX, avgAccuracy_LowSim, seAccuracy_LowSim, 'bo-', ...
                  'LineWidth', 1.5, ...       % Line width
                  'MarkerSize', 8, ...        % Marker size
                  'DisplayName', sprintf('Low Similarity (Condition %d)', lowSimCondition));
    % Customize error bar appearance (optional)
    set(h1, 'CapSize', 10);

    % Plot High Similarity data with error bars
    h2 = errorbar(plotX, avgAccuracy_HighSim, seAccuracy_HighSim, 'rs-', ...
                  'LineWidth', 1.5, ...       % Line width
                  'MarkerSize', 8, ...        % Marker size
                  'DisplayName', sprintf('High Similarity (Condition %d)', highSimCondition));
     % Customize error bar appearance (optional)
    set(h2, 'CapSize', 10);

    % Set plot properties
    xlabel('Delay (seconds)'); % Label still refers to the conceptual delay
    ylabel('Average Accuracy');
    title('Working Memory Performance Across Delay Conditions');

    % Set x-axis ticks to the equally spaced plot points and apply custom labels
    xticks(plotX);         % Set tick locations to the equally spaced points
    xticklabels(xTickLabels); % Apply the custom string labels

    % Set y-axis limits and ticks
    ylim([-0.05, 1.05]); % Slightly extend limits for better visualization of 0 and 1
    yticks(0:0.1:1);     % Set y-axis ticks at 0.1 intervals

    % Add legend
    legend('Location', 'southwest', 'AutoUpdate', 'off'); % 'AutoUpdate', 'off' prevents legend from adding new plots automatically

    % Remove grid (default is no grid, but explicitly turn off if it were on)
    grid off;

    hold off; % Release the hold on the current axes
end




function plotDistractorErrorFrequency(touchTable)
    % Extract relevant data
    distractorErrorRows = touchTable(touchTable.DistractorError == true & touchTable.isRepetition == 0, :); % Filter rows with distractor error
    
    % Initialize an array for adjusted states
    adjustedStates = zeros(height(distractorErrorRows), 1);

    % Loop through the filtered rows and adjust the states
    for i = 1:height(distractorErrorRows)
        % Check if this is the first touch in the trial
        if distractorErrorRows.TouchNumber(i) == 1
            adjustedStates(i) = distractorErrorRows.CurrentState(i) + 1;
        else
            % Check if the previous row (in the same trial) had state 0
            currentTrial = distractorErrorRows.TrialNumberInBlock(i);
            currentSub = distractorErrorRows.Subject(i);
            currentSession = distractorErrorRows.SessionNumber(i);
            currentBlock = distractorErrorRows.BlockNumber(i);
            previousRowIdxs = find( ...
                                    touchTable.SessionNumber == currentSession & ...
                                    touchTable.TrialNumberInBlock == currentTrial & ...
                                    touchTable.BlockNumber == currentBlock & ...
                                    touchTable.TouchNumber == distractorErrorRows.TouchNumber(i) - 1 ...
                                    );
            tempidx = find(strcmp(touchTable(previousRowIdxs,:).Subject, currentSub));
            previousRowIdx = previousRowIdxs(tempidx);
            if touchTable.CurrentState(previousRowIdx) == 0
                adjustedStates(i) = distractorErrorRows.CurrentState(i) + 1;
            else
                adjustedStates(i) = distractorErrorRows.CurrentState(i) + 2;
            end
        end
    end

    % Get unique adjusted states
    uniqueStates = unique(adjustedStates);

    % Split data into before and after swap
    beforeSwapRows = distractorErrorRows.isSwapped == false;
    afterSwapRows = distractorErrorRows.isSwapped == true;

    % Initialize arrays to store frequencies and CI for both groups
    frequenciesBefore = zeros(length(uniqueStates), 1);
    ci95Before = zeros(length(uniqueStates), 1);
    
    frequenciesAfter = zeros(length(uniqueStates), 1);
    ci95After = zeros(length(uniqueStates), 1);

    % Initialize array to store p-values
    p_values = zeros(length(uniqueStates), 1);

    % Calculate frequencies, CI, and perform Z-test for proportions for each state
    for i = 1:length(uniqueStates)
        state = uniqueStates(i);
        
        % For before swap
        stateRowsBefore = adjustedStates == state & beforeSwapRows;
        freqBefore = sum(stateRowsBefore) / sum(beforeSwapRows);
        frequenciesBefore(i) = freqBefore;
        nBefore = sum(beforeSwapRows);
        ci95Before(i) = 1.96 * sqrt((freqBefore * (1 - freqBefore)) / nBefore);
        fprintf('%d: %.2f ± %.2f; \n',i, freqBefore, ci95Before(i));
        % For after swap
        stateRowsAfter = adjustedStates == state & afterSwapRows;
        freqAfter = sum(stateRowsAfter) / sum(afterSwapRows);
        frequenciesAfter(i) = freqAfter;
        nAfter = sum(afterSwapRows);
        ci95After(i) = 1.96 * sqrt((freqAfter * (1 - freqAfter)) / nAfter);
		
		fprintf('%d: %.2f ± %.2f; \n',i, freqAfter, ci95After(i));

        % Perform Z-test for proportions
        p_pool = (sum(stateRowsBefore) + sum(stateRowsAfter)) / (nBefore + nAfter);
        z_stat = (freqBefore - freqAfter) / sqrt(p_pool * (1 - p_pool) * (1 / nBefore + 1 / nAfter));
        p_values(i) = 2 * (1 - normcdf(abs(z_stat))); % Two-tailed test

		fprintf('p = %.4f \n', p_values(i));

    end

     % Define custom colors
    colorBefore = [162, 34, 91] / 255;
    colorAfter = [34,74,137] / 255;

    % Create the dot plot with error bars
    figure;
    % Plot before swap
    errorbar(uniqueStates, frequenciesBefore, ci95Before, '-o', 'DisplayName', 'Before Swap', ...
        'MarkerFaceColor', colorBefore, 'MarkerEdgeColor', 'none', 'Color', colorBefore);
    hold on;
    
    % Plot after swap
    errorbar(uniqueStates, frequenciesAfter, ci95After, '-o', 'DisplayName', 'After Swap', ...
        'MarkerFaceColor', colorAfter, 'MarkerEdgeColor', 'none', 'Color', colorAfter);
    
    % Add significance stars based on p-values
	max_y = 0;
	for i = 1:length(uniqueStates)
		max_yx = max(frequenciesBefore(i) + ci95Before(i), frequenciesAfter(i) + ci95After(i)); % Get max y-value for star placement
		max_y = max(max_y, max_yx);
	end
    for i = 1:length(uniqueStates)
%         max_y = max(frequenciesBefore(i) + ci95Before(i), frequenciesAfter(i) + ci95After(i)); % Get max y-value for star placement
        if p_values(i) <= 0.001
            text(uniqueStates(i), max_y + 0.02, '***', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
        elseif p_values(i) <= 0.01
            text(uniqueStates(i), max_y + 0.02, '**', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
        elseif p_values(i) <= 0.05
            text(uniqueStates(i), max_y + 0.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
        end
    end
    
    hold off;
    
    % Plot settings
    xlabel('State (Adjusted State)', 'FontSize', 12);
    ylabel('Frequency of Distractor Errors', 'FontSize', 12);
    title('Frequency of Distractor Errors Before and After Swap', 'FontSize', 12);
    xlim([0.5, 5.5]);
    % Set plot properties
    set(gca, 'TickDir', 'out', 'FontSize', 12); % Ticks out, font size 12
    box off; % No box around the plot

	legend_obj = legend('Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    hold off;
    
    % Adjust figure size to accommodate legend if it exists
    canvas_size = [300, 250]; % Define canvas size
    if exist('legend_obj', 'var')
        legend_pos = get(legend_obj, 'Position');
        legend_width = legend_pos(3);
        new_fig_width = canvas_size(1) + legend_width * canvas_size(1);
        set(gcf, 'Position', [350, 250, new_fig_width, canvas_size(2)]);
    else
        set(gcf, 'Position', [350, 250, canvas_size(1), canvas_size(2)]); % Canvas size
    end

    % Ensure the renderer is set to 'Painters'
    set(gcf, 'renderer', 'Painters');
end


function plotErrorDecreasingRateMemory(touch_table)
    % Set canvas size
    canvas_size = [100, 100, 450, 250];

    % Merge all subjects for analysis
    subjects = unique(touch_table.Subject);
    num_subjects = length(subjects);

    % Filter for non-repetitive trials
%     touch_table = touch_table(touch_table.isRepetition == 0, :);

    % Initialize storage for individual subject error rates
    subject_error_rates = cell(num_subjects, 31);  % {subject, trial} for each trial per subject

    % Process data for each session
    sessions = unique(touch_table.SessionNumber);
    
    for i = 1:length(sessions)
        session = sessions(i);
        session_rows = touch_table(touch_table.SessionNumber == session, :);
        blocks = unique(session_rows.BlockNumber);
        
        for block = blocks'
            block_rows = session_rows(session_rows.BlockNumber == block, :);

            for trial = 1:15  % For trials 1–15 with isSwapped == 0
                for s = 1:num_subjects
                    subject_rows = block_rows(strcmp(block_rows.Subject, subjects{s}) & block_rows.TrialNumberInBlock == trial & block_rows.isRepetition == 0, :);
                    total_touches = height(subject_rows);

                    if total_touches == 0
                        continue  % Skip if there were no touches for this trial
                    end

                    % Calculate Distractor Errors
                    distractor_touches = sum(subject_rows.TouchObjectCorrectOrdinalPosition == 0);
                    error_rate = distractor_touches / total_touches;

                    % Save error rate as {subject, trial}
                    subject_error_rates{s, trial} = [subject_error_rates{s, trial}; error_rate];
                end
            end

            for trial = 17:31  % For trials 16–30 with isSwapped == 1
                for s = 1:num_subjects
                    subject_rows = block_rows(strcmp(block_rows.Subject, subjects{s}) & block_rows.TrialNumberInBlock == trial - 16 & block_rows.isRepetition == 1 & block_rows.isNewRepetition == 0, :);
                    total_touches = height(subject_rows);

                    if total_touches == 0
                        continue  % Skip if there were no touches for this trial
                    end

                    % Calculate Distractor Errors
                    distractor_touches = sum(subject_rows.TouchObjectCorrectOrdinalPosition == 0);
                    error_rate = distractor_touches / total_touches;

                    % Save error rate as {subject, trial}
                    subject_error_rates{s, trial} = [subject_error_rates{s, trial}; error_rate];
                end
            end
        end
    end

    % Plotting
    fig = figure;
    set(gcf, 'renderer', 'Painters');  % Set renderer for high-quality vector graphics

    hold on;
    colors = [39, 93, 44;
              58, 112, 175;
              233, 166, 64;
              199, 54, 55;
              107, 37, 110;
              247, 211, 76] / 255;  % Color matching for subjects

    % Plot individual subject lines
    for s = 1:num_subjects
        subject_mean = zeros(1, 30);  % Initialize mean error rates for 30 trials

        for trial = 1:31
            if ~isempty(subject_error_rates{s, trial})
                subject_mean(trial) = mean(subject_error_rates{s, trial}, 'omitnan');
            else
                subject_mean(trial) = NaN;  % Handle missing data
            end
        end

        % Plot the subject's line
        plot(1:31, subject_mean, '-', 'Color', [colors(s,:), 0.7], 'LineWidth', 1, ...
            'DisplayName', ['Subject ', subjects{s}]);
    end

    % Calculate mean error rates across all subjects
    mean_rate = zeros(1, 31);
    sem = zeros(1, 31);
    ci = zeros(1, 31);

    for trial = 1:31
        all_subject_errors = [];
        for s = 1:num_subjects
            if ~isempty(subject_error_rates{s, trial})
                all_subject_errors = [all_subject_errors; subject_error_rates{s, trial}];
            end
		end
		if ~isempty(all_subject_errors)
        	mean_rate(trial) = mean(all_subject_errors, 'omitnan');
        	sem(trial) = std(all_subject_errors, 0, 1, 'omitnan') / sqrt(size(all_subject_errors, 1));
        	ci(trial) = 1.96 * sem(trial);
		else
			mean_rate(trial) = NaN;
			sem(trial) = NaN;
			ci(trial) = NaN;
		end
	end

	% Plot triangle for average value of all trials at x = 15.5 and x = 31.5
    avg_value_15 = mean(mean_rate(1:15), 'omitnan');
    avg_value_31 = mean(mean_rate(17:31), 'omitnan');
    err_15 = std(mean_rate(1:15), 'omitnan') / sqrt(15);
    err_31 = std(mean_rate(17:31), 'omitnan') / sqrt(15);

    % Plot triangle for x = 15.5
    plot(15.5, avg_value_15, '<', 'MarkerSize', 8, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(5,:));
    errorbar(15.5, avg_value_15, 1.96*err_15, 'Color', colors(5,:), 'CapSize', 8, 'LineWidth', 1);
	fprintf('%.3f +- %.3f \n',avg_value_15, 1.96*err_15);
    % Plot triangle for x = 31.5
    plot(31.5, avg_value_31, '<', 'MarkerSize', 8, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(5,:));
    errorbar(31.5, avg_value_31, 1.96*err_31, 'Color', colors(5,:), 'CapSize', 8, 'LineWidth', 1);
	fprintf('%.3f +- %.3f \n',avg_value_31, 1.96*err_31);


    % Plot averaged error rate (thicker line)
    plot(1:31, mean_rate, '-', 'Color', colors(5,:), 'LineWidth', 2, 'DisplayName', 'Average');
    fill([1:15, fliplr(1:15)], [mean_rate(1:15) - ci(1:15), fliplr(mean_rate(1:15) + ci(1:15))], colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
	hold on;
	fill([17:31, fliplr(17:31)], [mean_rate(17:31) - ci(17:31), fliplr(mean_rate(17:31) + ci(17:31))], colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % Dashed grey line at x = 15
    plot([15, 15], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');

    % Customize plot
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Proportion of Distractor Errors', 'FontSize', 12);
    title('Proportion of Distractor Errors Across Trials (Repetition)', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

    % Custom x-tick labels for trials 16–30
    xticks([1, 5, 10, 15, 17, 21, 26, 31]);
    xticklabels({'1', '5', '10', '15', '1', '5', '10', '15'});

    xlim([0.5, 31.5]);
	ylim([0, 0.14]);
    % Set canvas size
    fig.Position = canvas_size;

    % Create custom legend with subject lines
    legendHandles = gobjects(1, num_subjects + 1);  % Including the average line
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, '-', 'Color', colors(s,:), 'LineWidth', 1);
    end

    % Add a handle for the averaged line
    legendHandles(end) = plot(NaN, NaN, '-', 'Color', colors(5,:), 'LineWidth', 2);

    % Create the legend with subject names and 'Average'
    legend_obj = legend(legendHandles, [cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false); {'Average'}], ...
                        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Box = 'off';

    % Adjust figure size to accommodate the legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function plotErrorDecreasingRateSwap(touch_table)
    % Set canvas size
    canvas_size = [100, 100, 450, 250];

    % Merge all subjects for analysis
    subjects = unique(touch_table.Subject);
    num_subjects = length(subjects);

    % Filter for non-repetitive trials
    touch_table = touch_table(touch_table.isRepetition == 0, :);

    % Initialize storage for individual subject error rates
    subject_error_rates = cell(num_subjects, 31);  % {subject, trial} for each trial per subject

    % Process data for each session
    sessions = unique(touch_table.SessionNumber);
    
    for i = 1:length(sessions)
        session = sessions(i);
        session_rows = touch_table(touch_table.SessionNumber == session, :);
        blocks = unique(session_rows.BlockNumber);
        
        for block = blocks'
            block_rows = session_rows(session_rows.BlockNumber == block, :);

            for trial = 1:15  % For trials 1–15 with isSwapped == 0
                for s = 1:num_subjects
                    subject_rows = block_rows(strcmp(block_rows.Subject, subjects{s}) & block_rows.TrialNumberInBlock == trial & block_rows.isSwapped == 0, :);
                    total_touches = height(subject_rows);

                    if total_touches == 0
                        continue  % Skip if there were no touches for this trial
                    end

                    % Calculate Distractor Errors
                    distractor_touches = sum(subject_rows.TouchObjectCorrectOrdinalPosition == 0);
                    error_rate = distractor_touches / total_touches;

                    % Save error rate as {subject, trial}
                    subject_error_rates{s, trial} = [subject_error_rates{s, trial}; error_rate];
                end
            end

            for trial = 17:31  % For trials 16–30 with isSwapped == 1
                for s = 1:num_subjects
                    subject_rows = block_rows(strcmp(block_rows.Subject, subjects{s}) & block_rows.TrialNumberInBlock == trial - 16 & block_rows.isSwapped == 1, :);
                    total_touches = height(subject_rows);

                    if total_touches == 0
                        continue  % Skip if there were no touches for this trial
                    end

                    % Calculate Distractor Errors
                    distractor_touches = sum(subject_rows.TouchObjectCorrectOrdinalPosition == 0);
                    error_rate = distractor_touches / total_touches;

                    % Save error rate as {subject, trial}
                    subject_error_rates{s, trial} = [subject_error_rates{s, trial}; error_rate];
                end
            end
        end
    end

    % Plotting
    fig = figure;
    set(gcf, 'renderer', 'Painters');  % Set renderer for high-quality vector graphics

    hold on;
    colors = [39, 93, 44;
              58, 112, 175;
              233, 166, 64;
              199, 54, 55;
              107, 37, 110;
              247, 211, 76] / 255;  % Color matching for subjects

    % Plot individual subject lines
    for s = 1:num_subjects
        subject_mean = zeros(1, 30);  % Initialize mean error rates for 30 trials

        for trial = 1:31
            if ~isempty(subject_error_rates{s, trial})
                subject_mean(trial) = mean(subject_error_rates{s, trial}, 'omitnan');
            else
                subject_mean(trial) = NaN;  % Handle missing data
            end
        end

        % Plot the subject's line
        plot(1:31, subject_mean, '-', 'Color', [colors(s,:), 0.7], 'LineWidth', 1, ...
            'DisplayName', ['Subject ', subjects{s}]);
    end

    % Calculate mean error rates across all subjects
    mean_rate = zeros(1, 31);
    sem = zeros(1, 31);
    ci = zeros(1, 31);

    for trial = 1:31
        all_subject_errors = [];
        for s = 1:num_subjects
            if ~isempty(subject_error_rates{s, trial})
                all_subject_errors = [all_subject_errors; subject_error_rates{s, trial}];
            end
		end
		if ~isempty(all_subject_errors)
        	mean_rate(trial) = mean(all_subject_errors, 'omitnan');
        	sem(trial) = std(all_subject_errors, 0, 1, 'omitnan') / sqrt(size(all_subject_errors, 1));
        	ci(trial) = 1.96 * sem(trial);
		else
			mean_rate(trial) = NaN;
			sem(trial) = NaN;
			ci(trial) = NaN;
		end
	end

    % Plot triangle for average value of all trials at x = 15.5 and x = 31.5
    avg_value_15 = mean(mean_rate(1:15), 'omitnan');
    avg_value_31 = mean(mean_rate(17:31), 'omitnan');
    err_15 = std(mean_rate(1:15), 'omitnan') / sqrt(15);
    err_31 = std(mean_rate(17:31), 'omitnan') / sqrt(15);

    % Plot triangle for x = 15.5
    plot(15.5, avg_value_15, '<', 'MarkerSize', 8, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(5,:));
    errorbar(15.5, avg_value_15, 1.96*err_15, 'Color', colors(5,:), 'CapSize', 8, 'LineWidth', 1);
	fprintf('%.3f +- %.3f \n',avg_value_15, 1.96*err_15);
    % Plot triangle for x = 31.5
    plot(31.5, avg_value_31, '<', 'MarkerSize', 8, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(5,:));
    errorbar(31.5, avg_value_31, 1.96*err_31, 'Color', colors(5,:), 'CapSize', 8, 'LineWidth', 1);
	fprintf('%.3f +- %.3f \n',avg_value_31, 1.96*err_31);

    % Plot averaged error rate (thicker line)
    plot(1:31, mean_rate, '-', 'Color', colors(5,:), 'LineWidth', 2, 'DisplayName', 'Average');
    fill([1:15, fliplr(1:15)], [mean_rate(1:15) - ci(1:15), fliplr(mean_rate(1:15) + ci(1:15))], colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
	hold on;
	fill([17:31, fliplr(17:31)], [mean_rate(17:31) - ci(17:31), fliplr(mean_rate(17:31) + ci(17:31))], colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % Dashed grey line at x = 15
    plot([16, 16], 0, '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');

    % Customize plot
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Proportion of Distractor Errors', 'FontSize', 12);
    title('Proportion of Distractor Errors Across Trials (isSwapped)', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');

    % Custom x-tick labels for trials 16–30
    xticks([1, 5, 10, 15, 17, 21, 26, 31]);
    xticklabels({'1', '5', '10', '15', '1', '5', '10', '15'});

    xlim([0.5, 31.5]);

    % Set canvas size
    fig.Position = canvas_size;

    % Create custom legend with subject lines
    legendHandles = gobjects(1, num_subjects + 1);  % Including the average line
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, '-', 'Color', colors(s,:), 'LineWidth', 1);
    end

    % Add a handle for the averaged line
    legendHandles(end) = plot(NaN, NaN, '-', 'Color', colors(5,:), 'LineWidth', 2);

    % Create the legend with subject names and 'Average'
    legend_obj = legend(legendHandles, [cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false); {'Average'}], ...
                        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Box = 'off';
	ylim([0, 0.14]);

    % Adjust figure size to accommodate the legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end


function plotAvgSearchTimeCorrectTouchesRepeatConditions(touch_table, block_table)
    % Define conditions
    conditions = {'New Early', 'Repeat', 'New Later'};
%     colors = { hex2rgb('#507EBA'),hex2rgb('#D19723'), hex2rgb('#B32927')};
colors = {[20, 52, 94] ./ 255,[213, 197, 131] ./ 255,  [213,174, 200]./255};
    shapes = {'o', 's', 'd', '^'};

    % Get unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);
    numConditions = length(conditions);

    % Initialize storage for session-wise data
    sessionData = cell(numConditions, 1);

    % Process data
    for s = 1:numSubjects
        subject = subjects{s};
        subjectData = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subjectData.SessionNumber);
        
        for sess = sessions'
			sessionData{1} = [sessionData{1}; processCondition(subjectData, block_table, sess, 0, 0, subject)];
            sessionData{2} = [sessionData{2}; processCondition(subjectData, block_table, sess, 1, 0, subject)];
            sessionData{3} = [sessionData{3}; processCondition(subjectData, block_table, sess, 1, 1, subject)];
        end
    end

    % Define canvas size
    canvas_size = [100, 100, 280, 300];

    % Create figure
    fig = figure('Position', canvas_size);
    hold on;

    % Plot data for each condition
    for i = 1:numConditions
        % Add shapes for each subject
        for s = 1:numSubjects
            subjectIndices = strcmp(sessionData{i}(:,2), subjects{s});
            sw = swarmchart(repmat(i, sum(subjectIndices), 1), cell2mat(sessionData{i}(subjectIndices, 1)), 30, shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.3,'MarkerEdgeColor', 'none', 'LineWidth', 1);
			sw.XJitterWidth = 0.5;
		end

        % Add box plot
        boxplot(cell2mat(sessionData{i}(:,1)), 'Positions', i, 'Colors', colors{i}, 'Width', 0.3, 'Symbol', '');
        set(findobj(gca,'type','line'),'linew',2);

        % Calculate mean and 95% CI
        meanVal = mean(cell2mat(sessionData{i}(:,1)), 'omitnan');
        ci = 1.96 * std(cell2mat(sessionData{i}(:,1)), 'omitnan') / sqrt(sum(~isnan(cell2mat(sessionData{i}(:,1)))));
        errorbar(i, meanVal, ci, 'k', 'LineWidth', 2, 'CapSize', 10);
		fprintf('%d, %.2f ± %.2f\n', i, meanVal, ci); 
    end

    % Perform Welch's t-test for each pair
    pairs = [1 2; 2 3; 1 3]; % Define pairs in the desired order
    maxY = max(cellfun(@(x) max(cell2mat(x(:,1))), sessionData));
    minY = min(cellfun(@(x) min(cell2mat(x(:,1))), sessionData));
    baseY = maxY * 1.1;
    stepY = maxY * 0.05;

    % Initialize sigY for each pair
    sigY = baseY + [0; stepY; 2*stepY];

    for i = 1:size(pairs, 1)
        [~, p] = ttest2(cell2mat(sessionData{pairs(i,1)}(:,1)), cell2mat(sessionData{pairs(i,2)}(:,1)), 'Vartype', 'unequal');
        if p < 0.001
            sigString = '***';
        elseif p < 0.01
            sigString = '**';
        elseif p < 0.05
            sigString = '*';
        else
            sigString = 'n.s.';
        end
        plot(pairs(i,:), [1 1]*sigY(i)-stepY/3, '-k');
		fprintf('%d vs %d, p-value: %.10f\n', pairs(i,1), pairs(i,2), p);
        text(mean(pairs(i,:)), sigY(i)+stepY/5, sigString, 'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    % Customize the plot
    xlabel('Condition', 'FontSize', 12);
    ylabel('Average Search Time (s)', 'FontSize', 12);
    title('Average Search Time for Correct Touches in Repeat Conditions', 'FontSize', 12);
    set(gca, 'XTick', 1:numConditions, 'XTickLabel', conditions, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    ylim([minY * 0.8, (max(sigY) + stepY) *1.05]);
	xlim([0.5, 3.5]);
	set(gcf,'renderer','Painters')

    % Create dummy plots for legend
    subjectHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        subjectHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end

    % Create legend for subjects only
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
	function data = processCondition(touch_table, block_table, session, isRepetition, isNewRepetition, subject)
    % Filter data for the current subject, session, and condition
    blockFilter = strcmp(block_table.Subject, subject) & ...
                  block_table.SessionNumber == session & ...
                  block_table.isRepetition == isRepetition & ...
                  block_table.isNewRepetition == isNewRepetition;
    
    validBlocks = block_table(blockFilter, {'Subject', 'SessionNumber', 'BlockNumber'});
    
    % Find correct touches for the valid blocks
    trialIndices = ismember(touch_table(:, {'Subject', 'SessionNumber', 'BlockNumber'}), validBlocks) & ...
                   (strcmp(touch_table.TouchCategory, 'correctSelection') | ...
                    strcmp(touch_table.TouchCategory, 'Correct'));
    
    correctTouches = touch_table(trialIndices, :);
    
    % Calculate average search time for correct touches
    avgSearchTime = mean(correctTouches.SearchTime, 'omitnan');
    
    data = [avgSearchTime, {subject}];
	end

end



function plotFirstTrialPosition80Repeat(touch_table)
    % Define trial range
    maxTrialNumber = 15;

    % Get unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define conditions
    conditions = {'New Early', 'Repeat', 'New Later'};
    numConditions = length(conditions);

    % Define colors and shapes
%     colors = {hex2rgb('#507EBA'),hex2rgb('#D19723'),  hex2rgb('#B32927')};
colors = {[20, 52, 94] ./ 255,[213, 197, 131] ./ 255,  [213,174, 200]./255};
    shapes = {'o', 's', 'd', '^'};

    % Initialize storage for session-wise data
    sessionData = cell(numConditions, 1);

    % Process data
    for s = 1:numSubjects
        subject = subjects{s};
        subjectData = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subjectData.SessionNumber);
        
        for sess = sessions'
            sessionData{2} = [sessionData{2}; processCondition(subjectData, sess, 1, 0, subject)];
            sessionData{1} = [sessionData{1}; processCondition(subjectData, sess, 0, 0, subject)];
            sessionData{3} = [sessionData{3}; processCondition(subjectData, sess, 1, 1, subject)];
        end
    end

    % Define canvas size
    canvas_size = [100, 100, 280, 300];

    % Create figure
    fig = figure('Position', canvas_size);
    hold on;

    % Plot data for each condition
    for i = 1:numConditions
        % Add shapes for each subject
        for s = 1:numSubjects
            subjectIndices = strcmp(sessionData{i}(:,2), subjects{s});
            sm = swarmchart(repmat(i, sum(subjectIndices), 1), cell2mat(sessionData{i}(subjectIndices, 1)), 30, shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceColor', colors{i},'MarkerFaceAlpha',0.3,'MarkerEdgeColor','none', 'LineWidth', 1);
            sm.XJitterWidth = 0.75;
        end

        % Add box plot
        boxplot(cell2mat(sessionData{i}(:,1)), 'Positions', i, 'Colors', colors{i}, 'Width', 0.3, 'Symbol', '');
        set(findobj(gca,'type','line'),'linew',2);
		hold on;
        % Calculate mean and 95% CI
        meanVal = nanmean(cell2mat(sessionData{i}(:,1)));
        ci = 1.96 * nanstd(cell2mat(sessionData{i}(:,1))) / sqrt(size(sessionData{i}, 1));
		fprintf('%d: %.2f, %.2f\n', i, meanVal, ci);
        errorbar(i, meanVal, ci, 'k', 'LineWidth', 2);
    end

    % Perform Welch's t-test for each pair
    pairs = [1 2; 2 3; 1 3];  % Define pairs in the desired order
    maxY = max(cellfun(@(x) max(cell2mat(x(:,1))), sessionData));
    baseY = maxY * 1.15;
    stepY = maxY * 0.1;
    
    % Initialize sigY for each pair
    sigY = baseY + [0; stepY; 2*stepY];

    for i = 1:size(pairs, 1)
        [~, p] = ttest2(cell2mat(sessionData{pairs(i,1)}(:,1)), cell2mat(sessionData{pairs(i,2)}(:,1)), 'Vartype', 'unequal');
        
        if p < 0.001
            sigString = '***';
        elseif p < 0.01
            sigString = '**';
        elseif p < 0.05
            sigString = '*';
        else
            sigString = 'n.s.';
        end
        
        plot(pairs(i,:), [1 1]*sigY(i), '-k');
		fprintf('%d vs %d: p-value: %.10f\n', pairs(i,1), pairs(i,2), p);
        text(mean(pairs(i,:)), sigY(i)+stepY/4, sigString, 'HorizontalAlignment', 'center', 'FontSize', 20);
    end

    % Customize the plot
    xlabel('Condition', 'FontSize', 12);
    ylabel('Trial Position Reaching 80% Completion', 'FontSize', 12);
    title('First Trial Position Reaching 80% Completion Rate', 'FontSize', 12);
    set(gca, 'XTick', 1:numConditions, 'XTickLabel', conditions, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    ylim([0, (max(sigY) + stepY)*1.05]);
	xlim([0.5, 3.5]);
    yticks(1:2:15);
	set(gcf,'renderer','Painters')

    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = scatter(NaN, NaN, 36, shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceColor','k','MarkerEdgeColor', 'k','LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);


    hold off;
end

function data = processCondition(subjectData, session, isRepetition, isNewRepetition, subject)
    sessionData = subjectData(subjectData.SessionNumber == session, :);
    blockData = sessionData(sessionData.isRepetition == isRepetition & sessionData.isNewRepetition == isNewRepetition, :);
    
    completedTrials = zeros(15, 1);
    totalTrials = zeros(15, 1);
    
    for trial = 1:15
        trialData = blockData(blockData.TrialNumberInBlock == trial, :);
        if ~isempty(trialData)
            totalTrials(trial) = height(trialData);
            completedTrials(trial) = sum(trialData.isRewarded);
        end
    end
    
    proportion = completedTrials ./ totalTrials;
    firstPosition80 = find(proportion >= 0.8, 1, 'first');
    if isempty(firstPosition80)
        firstPosition80 = nan;
    end
    
    data = [firstPosition80, {subject}];
end



function plotAvgSearchTimeByStateAndSwap(data)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];
    
    % Filter only correct touches
    correctTouches = data(strcmp(data.TouchCategory, 'Correct') | strcmp(data.TouchCategory, 'correctSelection'), :);
    
    % Define states and subjects
    states = 1:5;
    subjects = unique(correctTouches.Subject);
    
    % Colors for before and after swap
    colorBeforeSwap = '#4A7298';
    colorAfterSwap = '#F3C846';
    
    % Markers for subjects
    markers = {'o', 's', 'd', '^'};
    
    % Initialize cell array to store search times
    searchTimes = cell(2, 5);
    
    % Collect search times for each state and swap condition
    for i = 1:5
        searchTimes{1, i} = correctTouches.SearchTime(correctTouches.CurrentState == i & correctTouches.isSwapped == 0);
        searchTimes{2, i} = correctTouches.SearchTime(correctTouches.CurrentState == i & correctTouches.isSwapped == 1);
    end
    
    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;
    
    % Create box plot
    boxplot_data = [];
    group_labels = [];
    for i = 1:5
        boxplot_data = [boxplot_data; searchTimes{1,i}; searchTimes{2,i}];
        group_labels = [group_labels; repmat({sprintf('Ordinal Position %d Before', i)}, length(searchTimes{1,i}), 1); ...
            repmat({sprintf('Ordinal Position %d After', i)}, length(searchTimes{2,i}), 1)];
    end
    h = boxplot(boxplot_data, group_labels, 'Colors', [hex2rgb(colorBeforeSwap); hex2rgb(colorAfterSwap)], 'Width', 0.7);
    
    % Remove the '+' signs (outliers) from the boxplot
    outliers = findobj(h, 'Tag', 'Outliers');
    delete(outliers);
    
    % Calculate average search times and standard errors
    avg_before = cellfun(@mean, searchTimes(1,:));
    avg_after = cellfun(@mean, searchTimes(2,:));
    se_before = cellfun(@(x) std(x)/sqrt(length(x)), searchTimes(1,:));
    se_after = cellfun(@(x) std(x)/sqrt(length(x)), searchTimes(2,:));
    
    % Plot average lines
    plot(1:2:9, avg_before, 'Color', colorBeforeSwap, 'LineWidth', 2);
    plot(2:2:10, avg_after, 'Color', colorAfterSwap, 'LineWidth', 2);
    
    % Plot average points with error bars
    errorbar(1:2:9, avg_before, 1.96*se_before, 'o', 'Color', colorBeforeSwap, 'MarkerFaceColor', colorBeforeSwap, 'MarkerSize', 3, 'LineWidth', 1.5);
    errorbar(2:2:10, avg_after, 1.96*se_after, 'o', 'Color', colorAfterSwap, 'MarkerFaceColor', colorAfterSwap, 'MarkerSize', 3, 'LineWidth', 1.5);
    
    % Plot individual subject data points
    legendEntries = gobjects(1, length(subjects));
    for s = 1:length(subjects)
        subject = subjects{s};
        for i = 1:5
            % Before swap
            subjectData = correctTouches(strcmp(correctTouches.Subject, subject) & correctTouches.CurrentState == i & correctTouches.isSwapped == 0, :);
            if ~isempty(subjectData)
                h = scatter(2*i-1 + (rand-0.5)*0.3, mean(subjectData.SearchTime), 50, markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', colorBeforeSwap, 'MarkerFaceAlpha', 0.3);
                if i == 1 % Only add to legend for the first state
                    legendEntries(s) = h;
                end
            end
            % After swap
            subjectData = correctTouches(strcmp(correctTouches.Subject, subject) & correctTouches.CurrentState == i & correctTouches.isSwapped == 1, :);
            if ~isempty(subjectData)
                scatter(2*i + (rand-0.5)*0.3, mean(subjectData.SearchTime), 50, markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', colorAfterSwap, 'MarkerFaceAlpha', 0.3);
            end
        end
    end
    
    % Customize plot
    xlabel('Ordinal Position', 'FontSize', 12);
    ylabel('Search Time (s)', 'FontSize', 12);
    title('Search Time by State: Before vs After Swap', 'FontSize', 14);
    set(gca, 'XTick', 1.5:2:9.5, 'XTickLabel', {'1', '2', '3', '4', '5'});
    set(gcf,'renderer','Painters')
   % Add legend with black markers and average lines
	legendEntries = gobjects(1, length(subjects) + 2);
	for s = 1:length(subjects)
    	legendEntries(s) = scatter(NaN, NaN, 50, markers{s}, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
	end
	legendEntries(end-1) = plot(NaN, NaN, 'Color', colorBeforeSwap, 'LineWidth', 2);
	legendEntries(end) = plot(NaN, NaN, 'Color', colorAfterSwap, 'LineWidth', 2);
	legend_labels = cell(1,6);
	temp = [cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
                 	];
	for i = 1:length(subjects)
		legend_labels{i} = temp{i};
	end
	legend_labels{length(subjects) + 1} = 'Before Swap Average'; 
	legend_labels{length(subjects) + 2} = 'After Swap Average' ;
	legend_obj = legend(legendEntries, legend_labels, 'FontSize', 12);
	legend_obj.EdgeColor = 'none';
	legend_obj.Location = 'eastoutside';
    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');
    ylim([0, 3]);
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
    
    hold off;
end

function plotAvgSearchTimeByOrdinalPositionLearning(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 400, 300];

    % Constants for the plot
    ordinalPositions = 1:5; % Assuming 5 as the maximum ordinal position
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;
    subjects = unique(touch_table.Subject);
    sessions = unique(touch_table.SessionNumber);
    accuracy_threshold = 0.8;

    % Define shapes for each subject
    shapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
    if length(subjects) > length(shapes)
        shapes = repmat(shapes, 1, ceil(length(subjects)/length(shapes)));
    end

    % Create the plot
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Prepare data for swarm chart and error bars
    allDataBefore = cell(1, length(ordinalPositions));
    allDataAfter = cell(1, length(ordinalPositions));
    maxSearchTime = 0;

    for s = 1:length(subjects)
        subject_data = touch_table(strcmp(touch_table.Subject, subjects{s}), :);
        
        for sess = sessions'
            session_data = subject_data(subject_data.SessionNumber == sess, :);
            blocks = unique(session_data.BlockNumber);
            session_completed_trials = [];
            session_block_performances = [];
            
            % Calculate learning speed for the session
            for blk = 1:numel(blocks)
                block_data = session_data(session_data.BlockNumber == blocks(blk), :);
                if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                    unique_trials = unique(block_data.TrialNumberInBlock);
                    numTrials = numel(unique_trials);
                    completed_trials = zeros(1, numTrials);
                    for trial = 1:numTrials
                        trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                        completed_trials(trial) = any(trial_data.isRewarded == 1);
                    end
                    padded_completed_trials = nan(1, 15);
                    padded_completed_trials(1:length(completed_trials)) = completed_trials;
                    session_completed_trials = [session_completed_trials; padded_completed_trials];
                    block_performance = sum(completed_trials) / numTrials;
                    session_block_performances = [session_block_performances; block_performance];
                end
            end
            
            % Calculate learning speed
            meanCompletionRate = nanmean(session_completed_trials, 1);
            learning_speed = 16;
            for i = 1:length(meanCompletionRate)
                if mean(meanCompletionRate(i:end)) >= accuracy_threshold && meanCompletionRate(i) >= accuracy_threshold
                    learning_speed = i;
                    break;
                end
            end
            
            % Separate data before and after learning
            correctTouchesBefore = session_data((strcmp(session_data.TouchCategory, 'Correct') | strcmp(session_data.TouchCategory, 'correctSelection')) & ...
                                                session_data.isRepetition == 0 & ...
                                                session_data.isSwapped == 0 & ...
                                                session_data.TrialNumberInBlock < learning_speed, :);
            correctTouchesAfter = session_data((strcmp(session_data.TouchCategory, 'Correct') | strcmp(session_data.TouchCategory, 'correctSelection')) & ...
                                               session_data.isRepetition == 0 & ...
                                               session_data.isSwapped == 0 & ...
                                               session_data.TrialNumberInBlock >= learning_speed, :);
            
            % Process data for each ordinal position
            for position = ordinalPositions
                % Before learning
                beforeData = correctTouchesBefore(correctTouchesBefore.TouchObjectCorrectOrdinalPosition == position, :);
                if ~isempty(beforeData)
                    meanSearchTime = mean(beforeData.SearchTime, 'omitnan');
                    allDataBefore{position} = [allDataBefore{position}; position-0.2, meanSearchTime, s];
                    maxSearchTime = max(maxSearchTime, meanSearchTime);
                end
                
                % After learning
                afterData = correctTouchesAfter(correctTouchesAfter.TouchObjectCorrectOrdinalPosition == position, :);
                if ~isempty(afterData)
                    meanSearchTime = mean(afterData.SearchTime, 'omitnan');
                    allDataAfter{position} = [allDataAfter{position}; position+0.2, meanSearchTime, s];
                    maxSearchTime = max(maxSearchTime, meanSearchTime);
                end
            end
        end
    end

    % Plot data and error bars
    meanValuesBefore = zeros(1, length(ordinalPositions));
    meanValuesAfter = zeros(1, length(ordinalPositions));
    pValues = zeros(1, length(ordinalPositions));

    for position = ordinalPositions
        % Before learning (empty shapes)
        beforeData = cell2mat(allDataBefore(position));
        if ~isempty(beforeData)
            for s = 1:length(subjects)
                subjectData = beforeData(beforeData(:,3) == s, :);
                if ~isempty(subjectData)
                    sw = swarmchart(subjectData(:,1), subjectData(:,2), 36, colors(position,:), shapes{mod(s-1, length(shapes))+1}, ...
                            'MarkerEdgeColor', colors(position,:), 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
                    sw.XJitter = 'randn';
                    sw.XJitterWidth = 0.25;
                end
            end
            meanValue = mean(beforeData(:,2));
            meanValuesBefore(position) = meanValue;
            stdError = std(beforeData(:,2)) / sqrt(size(beforeData, 1));
            ciError = 1.96 * stdError;
            errorbar(position-0.2, meanValue, ciError, 'Color',  hex2rgb('#DB432C'), 'LineWidth', 2, 'CapSize', 10);
			fprintf('%.2f +- %.2f\n', meanValue, ciError);
        end
        
        % After learning (filled shapes)
        afterData = cell2mat(allDataAfter(position));
        if ~isempty(afterData)
            for s = 1:length(subjects)
                subjectData = afterData(afterData(:,3) == s, :);
                if ~isempty(subjectData)
                    sw = swarmchart(subjectData(:,1), subjectData(:,2), 36, colors(position,:), shapes{mod(s-1, length(shapes))+1}, ...
                            'MarkerFaceColor', colors(position,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2);
                    sw.XJitter = 'randn';
                    sw.XJitterWidth = 0.25;
                end
            end
            meanValue = mean(afterData(:,2));
            meanValuesAfter(position) = meanValue;
            stdError = std(afterData(:,2)) / sqrt(size(afterData, 1));
            ciError = 1.96 * stdError;
            errorbar(position+0.2, meanValue, ciError, 'Color',  hex2rgb('#6C96CC'), 'LineWidth', 2, 'CapSize', 10);
			fprintf('%.2f +- %.2f\n', meanValue, ciError);
        end
        
        % Perform t-test
        [~, pValues(position)] = ttest2(beforeData(:,2), afterData(:,2));
    end

    % Connect all before learning points
    plot(ordinalPositions-0.2, meanValuesBefore, '-', 'Color', hex2rgb('#DB432C'), 'LineWidth', 2);

    % Connect all after learning points
    plot(ordinalPositions+0.2, meanValuesAfter, '-', 'Color', hex2rgb('#6C96CC'), 'LineWidth', 2);

    % Add significance annotations
    yLim = get(gca, 'YLim');
    yRange = yLim(2) - yLim(1);
    starY = yLim(2) + 0.05 * yRange;
    
    for position = ordinalPositions
        xPos = position;
        p = pValues(position);
        fprintf('Position %d: %.4f\n',xPos, p);
        if p < 0.001
            starText = '***';
        elseif p < 0.01
            starText = '**';
        elseif p < 0.05
            starText = '*';
        else
            starText = 'ns';
        end
        
        % Draw line
        line([xPos-0.2, xPos+0.2], [starY, starY], 'Color', 'k', 'LineWidth', 1);
        
        % Add star or 'ns'
        text(xPos, starY + 0.02 * yRange, starText, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    end

    % Adjust y-axis limit to accommodate stars
    set(gca, 'YLim', [yLim(1), starY + 0.1 * yRange]);

    % Customize the plot
    title('Search Time Distribution by Ordinal Position', 'FontSize', 14);
    xlabel('Ordinal Position', 'FontSize', 12);
    ylabel('Average Search Time (s)', 'FontSize', 12);
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    xlim([0.5, 5.5]);
    set(gcf, 'renderer', 'Painters');
    
    % Add legend
    legendHandles = gobjects(1, length(subjects) * 2);
    for s = 1:length(subjects)
        legendLabels(s*2) = {['Subject ' subjects{s}(1) ' (After)' ]};
        legendLabels(s*2-1) = {['Subject ' subjects{s}(1) ' (Before)']};

        emptyMarker = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 8);
        filledMarker = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        legendHandles(s*2-1) = emptyMarker;
        legendHandles(s*2) = filledMarker;
    end
    
    legend_obj = legend(legendHandles, legendLabels, 'FontSize', 10, 'Location', 'eastoutside');
    legend_obj.EdgeColor = 'none';

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function plotLearningPerformanceSwapBoxPlot(touch_table)
       % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    sessions = unique(touch_table.SessionNumber);

    % Initialize storage variables
    trials_80percent_not_swapped = cell(numel(subjects), 1);
    trials_80percent_swapped = cell(numel(subjects), 1);

    for s = 1:numel(subjects)
        subject = subjects{s};
        subject_data_not_swapped = [];
        subject_data_swapped = [];
        
        for sess = 1:numel(sessions)
            session = sessions(sess);
            
            % Process not swapped condition
            data_not_swapped = processCondition(touch_table, subject, session, 0, 0, 0);
            subject_data_not_swapped = [subject_data_not_swapped; cell2mat(data_not_swapped(1))];
            
            % Process swapped condition
            data_swapped = processCondition(touch_table, subject, session, 0, 0, 1);
            subject_data_swapped = [subject_data_swapped; cell2mat(data_swapped(1))];
        end
        
        trials_80percent_not_swapped{s} = subject_data_not_swapped;
        trials_80percent_swapped{s} = subject_data_swapped;
    end

    % Print results for each subject
    fprintf('Results by subject:\n');
    for s = 1:numel(subjects)
        subject = subjects{s};
        data_not_swapped = trials_80percent_not_swapped{s};
        data_swapped = trials_80percent_swapped{s};
        
        mean_not_swapped = mean(data_not_swapped, 'omitnan');
        se_not_swapped = std(data_not_swapped, 'omitnan') / sqrt(sum(~isnan(data_not_swapped)));
        
        mean_swapped = mean(data_swapped, 'omitnan');
        se_swapped = std(data_swapped, 'omitnan') / sqrt(sum(~isnan(data_swapped)));
        
        fprintf('Subject %s:\n', subject);
        fprintf('  Initial Blocks: %.2f ± %.2f\n', mean_not_swapped, se_not_swapped);
        fprintf('  Swapped Blocks: %.2f ± %.2f\n', mean_swapped, se_swapped);
    end

    % Calculate overall results
    all_not_swapped = cell2mat(trials_80percent_not_swapped);
    all_swapped = cell2mat(trials_80percent_swapped);

    overall_mean_not_swapped = mean(all_not_swapped, 'omitnan');
    overall_se_not_swapped = std(all_not_swapped, 'omitnan') / sqrt(sum(~isnan(all_not_swapped)));

    overall_mean_swapped = mean(all_swapped, 'omitnan');
    overall_se_swapped = std(all_swapped, 'omitnan') / sqrt(sum(~isnan(all_swapped)));

    % Perform Welch's t-test
    [~, p_value, ~, stats] = ttest2(all_not_swapped, all_swapped, 'Vartype', 'unequal');

    % Print overall results
    fprintf('\nOverall results:\n');
    fprintf('Initial Blocks: %.2f ± %.2f\n', overall_mean_not_swapped, overall_se_not_swapped);
    fprintf('Swapped Blocks: %.2f ± %.2f\n', overall_mean_swapped, overall_se_swapped);
    fprintf('p-value: %.4f\n', p_value);

    % Initialize storage variables
    trials_80percent_not_swapped = [];
    trials_80percent_swapped = [];
    subject_indices = [];
    session_indices = [];

    for s = 1:numel(subjects)
        subject = subjects{s};
        for sess = 1:numel(sessions)
            session = sessions(sess);
            
            % Process not swapped condition
            data_not_swapped = processCondition(touch_table, subject, session, 0, 0, 0);
            trials_80percent_not_swapped = [trials_80percent_not_swapped; cell2mat(data_not_swapped(1))];
            
            % Process swapped condition
            data_swapped = processCondition(touch_table, subject, session, 0, 0, 1);
            trials_80percent_swapped = [trials_80percent_swapped; cell2mat(data_swapped(1))];
            
            subject_indices = [subject_indices; s];
            session_indices = [session_indices; sess];
        end
    end

	  canvas_size = [100, 100, 300, 250];

    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    
    
    % Add individual data points
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x'};
    colors = {[0.2 0.4 0.6], [0.8 0.6 0.2]};
    
    for s = 1:numel(subjects)
        subj_data_not_swapped = trials_80percent_not_swapped(subject_indices == s);
        subj_data_swapped = trials_80percent_swapped(subject_indices == s);
        subj_sessions = session_indices(subject_indices == s);
        
        % Plot points for not swapped condition
        sw1 = swarmchart(ones(size(subj_data_not_swapped)), subj_data_not_swapped, 50, markers{s}, ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{1}, 'DisplayName', ['Subject ', subjects{s}], ...
            'MarkerFaceAlpha', 0.5);
        sw1.XJitterWidth = 0.5;
        
        % Plot points for swapped condition
        sw2 = swarmchart(2*ones(size(subj_data_swapped)), subj_data_swapped, 50, markers{s}, ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{2}, 'DisplayName', ['Subject ', subjects{s}], ...
            'MarkerFaceAlpha', 0.5);
        sw2.XJitterWidth = 0.5;
        
        % Connect dots from same session with gray lines
        for i = 1:length(subj_sessions)
            x1 = sw1.XData(i);
            y1 = sw1.YData(i);
            x2 = sw2.XData(i);
            y2 = sw2.YData(i);
            line([x1, x2], [y1, y2], 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1);
        end
	end

	% Create box plot
    boxplot([trials_80percent_not_swapped, trials_80percent_swapped], 'Labels', {'Initial Blocks', 'Swapped Blocks'}, 'Colors', 'k', 'Symbol','');
    set(findobj(gca,'type','line'),'linew',1.5);

    % Add labels and title
    xlabel('Block Type', 'FontSize', 12);
    ylabel('Trials to Reach 80% Completion Rate', 'FontSize', 12);
    title('Trials to Reach 80% Completion Rate: Initial vs Swapped Blocks', 'FontSize', 14);

    % Perform Welch's t-test
    [~, p_value, ~, stats] = ttest2(trials_80percent_not_swapped, trials_80percent_swapped, 'Vartype', 'unequal');
    
    % Calculate mean and standard error for both conditions
    mean_not_swapped = mean(trials_80percent_not_swapped, 'omitnan');
    se_not_swapped = std(trials_80percent_not_swapped, 'omitnan') / sqrt(sum(~isnan(trials_80percent_not_swapped)));
    
    mean_swapped = mean(trials_80percent_swapped, 'omitnan');
    se_swapped = std(trials_80percent_swapped, 'omitnan') / sqrt(sum(~isnan(trials_80percent_swapped)));
    
    % Print results
    fprintf('Initial Blocks: %.2f ± %.2f\n', mean_not_swapped, se_not_swapped);
    fprintf('Swapped Blocks: %.2f ± %.2f\n', mean_swapped, se_swapped);
    fprintf('p-value: %.10f\n', p_value);
	
    % Add significance line and star if significant
    if p_value < 0.05
        y_max = max([trials_80percent_not_swapped; trials_80percent_swapped]);
        line([1, 2], [y_max*1.1, y_max*1.1], 'Color', 'k', 'LineWidth', 1.5);
        if p_value < 0.001
            text(1.5, y_max*1.105, '***', 'FontSize', 12, 'HorizontalAlignment', 'center');
        elseif p_value < 0.01
            text(1.5, y_max*1.105, '**', 'FontSize', 12, 'HorizontalAlignment', 'center');
        elseif p_value < 0.05
            text(1.5, y_max*1.105, '*', 'FontSize', 12, 'HorizontalAlignment', 'center');
        end
    end

    % Add error bars
    for j = 1:2
        if j == 1
            data = trials_80percent_not_swapped;
        else
            data = trials_80percent_swapped;
        end
        SEM = nanstd(data) / sqrt(length(data(~isnan(data))));
        errorbar(j, nanmean(data), SEM, 'k', 'LineWidth', 2, 'CapSize', 10, 'HandleVisibility', 'off');
    end

    % Add legend for subjects
    legendHandles = gobjects(1, 4);
    for s = 1:4
        legendHandles(s) = scatter(NaN, NaN, 36, markers{mod(s-1, length(markers))+1}, 'MarkerFaceColor','k','MarkerEdgeColor', 'k','LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';
	set(gcf,'renderer','Painters')
    % Adjust figure properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');
    ylim([0, y_max*1.15]);
    xlim([0.5, 2.5]);
    yticks(1:2:15);

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;

	
	function data = processCondition(touch_table, subject, session, isRepetition, isNewRepetition, isSwapped)
    sessionData = touch_table(strcmp(touch_table.Subject, subject) & ...
                              touch_table.SessionNumber == session & ...
                              touch_table.isRepetition == isRepetition & ...
                              touch_table.isNewRepetition == isNewRepetition & ...
                              touch_table.isSwapped == isSwapped, :);
    
    completedTrials = zeros(15, 1);
    totalTrials = zeros(15, 1);
    
    for trial = 1:15
        trialData = sessionData(sessionData.TrialNumberInBlock == trial, :);
        if ~isempty(trialData)
            totalTrials(trial) = height(trialData);
            completedTrials(trial) = sum(trialData.isRewarded);
        end
    end
    
    proportion = completedTrials ./ totalTrials;
    firstPosition80 = find(proportion >= 0.8, 1, 'first');
    if isempty(firstPosition80)
        firstPosition80 = nan;
    end
    
    data = [firstPosition80, {subject}];
end


end


function calculateLearningPerformanceSwap(touch_table)
% Filter unique rows
unique_rows = unique(touch_table(:, {'Subject', 'SessionNumber', 'BlockNumber', 'TrialNumberInBlock', 'isRepetition', 'isSwapped', 'isRewarded'}));

% Extract unique subjects
subjects = unique(unique_rows.Subject);
accuracy_threshold = 0.8;
% Initialize storage variables
all_completed_trials = [];
all_swapped_statuses = [];
all_subjects = []; % Store subject information
swap_session_learning_speeds = [];
nonswap_session_learning_speeds = [];
% Iterate through subjects and sessions to process data
for s = 1:numel(subjects)
	subject = subjects{s};
	subject_rows = unique_rows(strcmp(unique_rows.Subject, subject), :);
	sessions = unique(subject_rows.SessionNumber);

	for i = 1:numel(sessions)
		session = sessions(i);
		session_rows = subject_rows(subject_rows.SessionNumber == session, :);
		blocks = unique(session_rows.BlockNumber);
		swap_session_complete_trials = [];
		nonswap_session_complete_trials = [];
		for block = blocks'
			block_indices = find(session_rows.BlockNumber == block);
			completed_trials = session_rows.isRewarded(block_indices);
			if any(session_rows.isRepetition(block_indices) == 1)
				continue
			end
			padded_completed_trials = nan(1, 15);
			padded_completed_trials(1:length(completed_trials)) = completed_trials';
			all_completed_trials = [all_completed_trials; padded_completed_trials];

			
			% Assuming isSwapped is constant for all trials within a block
			swapped_status = session_rows.isSwapped(block_indices(1));
			all_swapped_statuses = [all_swapped_statuses; swapped_status];
			all_subjects = [all_subjects; {subject}];

			if swapped_status
				swap_session_complete_trials = [swap_session_complete_trials; padded_completed_trials];
			else
				nonswap_session_complete_trials = [nonswap_session_complete_trials; padded_completed_trials];
			end
		end

		            % Calculate learning speed for the session
            if ~isempty(swap_session_complete_trials)
                mean_completion_rate = nanmean(swap_session_complete_trials, 1);
                learning_speed = NaN;
                for i = 1:length(mean_completion_rate)
                    % if mean(mean_completion_rate(i:end)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                    %     learning_speed = i;
                    %     break;
                    % end
					if i < length(mean_completion_rate) - 3
                    	if mean(mean_completion_rate(i:i+3)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                        	learning_speed = i;
                        	break;
						end
					else
						if mean(mean_completion_rate(i:end)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                        	learning_speed = i;
                        	break;
						end
					end
                end
                if ~isnan(learning_speed)
                    swap_session_learning_speeds = [swap_session_learning_speeds, learning_speed];

                end
			end

            if ~isempty(nonswap_session_complete_trials)
                mean_completion_rate = nanmean(nonswap_session_complete_trials, 1);
                learning_speed = NaN;
                for i = 1:length(mean_completion_rate)
                    if i < length(mean_completion_rate) - 3
                    	if mean(mean_completion_rate(i:i+3)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                        	learning_speed = i;
                        	break;
						end
					else
						if mean(mean_completion_rate(i:end)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                        	learning_speed = i;
                        	break;
						end
					end
                end
                if ~isnan(learning_speed)
                    nonswap_session_learning_speeds = [nonswap_session_learning_speeds, learning_speed];

                end
			end

	end
end

% Calculate learning trials
proportions = cumsum(all_completed_trials, 2) ./ max(cumsum(all_completed_trials, 2), [], 2);
[M, N] = size(proportions);
firstColumnIndex = NaN(M, 1);

for i = 1:M
	indices = find(proportions(i,:) >= 0.8, 1, 'first');
	if ~isempty(indices)
		firstColumnIndex(i) = indices;
	end
end
learning_trial = firstColumnIndex;
avg_speed_learning = mean(learning_trial, 'omitnan');
se_speed_learning = std(learning_trial, 'omitnan') / sqrt(length(learning_trial));

% Separate indices based on swap status
not_swapped_indices = all_swapped_statuses == 0;
swapped_indices = all_swapped_statuses == 1;

% Calculate averages and standard errors for swapped and not swapped
avg_not_swapped_blocks = mean(learning_trial(not_swapped_indices), 'omitnan');
se_not_swapped_blocks = std(learning_trial(not_swapped_indices), 'omitnan') / sqrt(sum(not_swapped_indices));

avg_swapped_blocks = mean(learning_trial(swapped_indices), 'omitnan');
se_swapped_blocks = std(learning_trial(swapped_indices), 'omitnan') / sqrt(sum(swapped_indices));

% Mean completion rates
meanCompletionRateNotSwapped = nanmean(all_completed_trials(not_swapped_indices,:), 1);
meanCompletionRateSwapped = nanmean(all_completed_trials(swapped_indices,:), 1);

trials = 1:size(all_completed_trials, 2);


% Sigmoid model fitting
sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', ...
	'independent', 'x', ...
	'coefficients', {'a', 'b', 'c'});
initialGuess = [1, 0.1, trials(round(end / 2))];
[fitResultNotSwapped, ~] = fit(trials', meanCompletionRateNotSwapped', sigmoidModel, 'StartPoint', initialGuess);
[fitResultSwapped, ~] = fit(trials', meanCompletionRateSwapped', sigmoidModel, 'StartPoint', initialGuess);

% Calculate SEM and CI bounds for not swapped and swapped
semNotSwapped = std(all_completed_trials(not_swapped_indices,:), [], 1, 'omitnan') / sqrt(sum(not_swapped_indices));
semSwapped = std(all_completed_trials(swapped_indices,:), [], 1, 'omitnan') / sqrt(sum(swapped_indices));
ciNotSwapped = 1.96 * semNotSwapped;
ciSwapped = 1.96 * semSwapped;
lowerBoundNotSwapped = meanCompletionRateNotSwapped - ciNotSwapped;
upperBoundNotSwapped = meanCompletionRateNotSwapped + ciNotSwapped;
lowerBoundSwapped = meanCompletionRateSwapped - ciSwapped;
upperBoundSwapped = meanCompletionRateSwapped + ciSwapped;

    % Calculate session-wise learning points
    [session_learning_points_not_swapped, session_learning_points_swapped] = calculateSessionLearningPoints(touch_table);

    fig1 = figure;
    canvas_size = [100, 100, 500, 250];
    fig1.Position = canvas_size;
    hold on;
    subject_handles = gobjects(1, 2+length(subjects));

    % Plot shaded areas
    fill([trials, fliplr(trials)], [lowerBoundNotSwapped, fliplr(upperBoundNotSwapped)], hex2rgb('#4A7298'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([trials+16, fliplr(trials+16)], [lowerBoundSwapped, fliplr(upperBoundSwapped)], hex2rgb('#F3C846'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Plot overall data fits
    subject_handles(1) = plot(trials, feval(fitResultNotSwapped, trials), 'Color', '#4A7298', 'LineWidth', 2);
    subject_handles(2) = plot(trials+16, feval(fitResultSwapped, trials), 'Color', '#F3C846', 'LineWidth', 2);
    block_types = {'Before Swapped Blocks', 'After Swapped Blocks'};

    % Plot individual subject dots with different shapes
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x'};
    subject_labels = cell(numel(subjects), 1);
    for s = 1:numel(subjects)
        subject = subjects{s};
        subj_indices = strcmp(all_subjects, subject);
        subject_labels{s} = ['Subject ', subject(1)];

        % Before swap
        subj_trials_not_swapped = all_completed_trials(subj_indices & not_swapped_indices, :);
        mean_subj_not_swapped = nanmean(subj_trials_not_swapped, 1);
        scatter(trials, mean_subj_not_swapped, 50, markers{mod(s-1, numel(markers))+1}, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#4A7298', 'MarkerFaceAlpha',0.5,'HandleVisibility', 'off');

        % After swap
        subj_trials_swapped = all_completed_trials(subj_indices & swapped_indices, :);
        mean_subj_swapped = nanmean(subj_trials_swapped, 1);
        scatter(trials+16, mean_subj_swapped, 50, markers{mod(s-1, numel(markers))+1}, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#F3C846', 'MarkerFaceAlpha',0.5,'HandleVisibility', 'off');
    end

    % Add inverse triangles for session-wise learning points
    errorbar(mean(nonswap_session_learning_speeds), 1.1, 1.96*  std(nonswap_session_learning_speeds)/sqrt(length(nonswap_session_learning_speeds)),1.96*  std(nonswap_session_learning_speeds)/sqrt(length(nonswap_session_learning_speeds)), 'horizontal', 'Color', '#4A7298', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    scatter(mean(nonswap_session_learning_speeds), 1.1, 100, 'v', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#4A7298', 'HandleVisibility', 'off');
    fprintf('before: %.2f  %.2f \n', mean(nonswap_session_learning_speeds),  1.96*  std(nonswap_session_learning_speeds)/sqrt(length(nonswap_session_learning_speeds)));
    errorbar(mean(swap_session_learning_speeds)+16, 1.1, 1.96* std(swap_session_learning_speeds)/sqrt(length(swap_session_learning_speeds)), 1.96* std(swap_session_learning_speeds)/sqrt(length(swap_session_learning_speeds)), 'horizontal', 'Color', '#F3C846', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    scatter(mean(swap_session_learning_speeds)+16, 1.1, 100, 'v', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', '#F3C846', 'HandleVisibility', 'off');
	fprintf('after: %.2f  %.2f \n', mean(swap_session_learning_speeds),  1.96*  std(swap_session_learning_speeds)/sqrt(length(swap_session_learning_speeds)));

    % Add lines for reference
    xline(15, '--', 'Color', [0.8,0.8,0.8], 'HandleVisibility', 'off');
    yline(0.8, '--', 'HandleVisibility', 'off');
    ylim([0.2, 1.2]);
	yticks([0.4, 0.7, 1]);
	xlim([0, 32])
	set(gcf,'renderer','Painters')

    % Add labels and title
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Completion Rate', 'FontSize', 12);
    title('Trial Completion Rates: Before vs After Swapped', 'FontSize', 14);
    set(gca, 'XTick', [[1,5:5:15], [1,5:5:15]+16], 'XTickLabel', [[1,5:5:15], [1,5:5:15]]);

    % Add legend for subjects
    
    for s = 3:2+length(subjects)
        subject_handles(s) = scatter(NaN, NaN, 100, 'k', markers{s-2}, 'filled');
	end
	subject_names = cellfun(@(x) ['Subject ' x(1)], subjects, 'UniformOutput', false);
    legend_obj = legend(subject_handles, [block_types, subject_names'], 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';

    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    hold off;
end


function [session_learning_points_not_swapped, session_learning_points_swapped] = calculateSessionLearningPoints(touch_table)
    subjects = unique(touch_table.Subject);
    accuracy_threshold = 0.8;
    session_learning_points_not_swapped = [];
    session_learning_points_swapped = [];

    for subj = 1:numel(subjects)
        subject_data = touch_table(strcmp(touch_table.Subject, subjects{subj}) & touch_table.isRepetition == 0, :);
        sessions = unique(subject_data.SessionNumber);
        
        for sess = 1:numel(sessions)
            session_data = subject_data(subject_data.SessionNumber == sessions(sess), :);
            blocks = unique(session_data.BlockNumber);
            session_completed_trials_beforeswap = [];
			session_completed_trials_afterswap = [];
            session_swapped = false;

            for blk = 1:numel(blocks)
                block_data = session_data(session_data.BlockNumber == blocks(blk), :);
                
					
					unique_trials = unique(block_data.TrialNumberInBlock);
                    numTrials = numel(unique_trials);
                    completed_trials = zeros(1, numTrials);
                    for trial = 1:numTrials
                        trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                        completed_trials(trial) = any(trial_data.isRewarded == 1);
                    end
                    padded_completed_trials = nan(1, 15);
                    padded_completed_trials(1:length(completed_trials)) = completed_trials;
                   

					if all(block_data.isSwapped == 0)
						session_completed_trials_beforeswap = [session_completed_trials_beforeswap; padded_completed_trials];
						session_swapped = false;
                	elseif all(block_data.isSwapped == 1)
						 session_completed_trials_afterswap = [session_completed_trials_afterswap; padded_completed_trials];
                    	session_swapped = true;
					end
				
            end

            if ~isempty(session_completed_trials_beforeswap)
                mean_completion_rate = nanmean(session_completed_trials_beforeswap, 1);
                learning_speed = NaN;
                for i = 1:length(mean_completion_rate)
                    if all(mean_completion_rate(i:end) >= accuracy_threshold)
                        learning_speed = i;
                        break;
                    end
                end
                if ~isnan(learning_speed) 
					session_learning_points_not_swapped = [session_learning_points_not_swapped, learning_speed];
                end
			end

			if ~isempty(session_completed_trials_afterswap)
                mean_completion_rate = nanmean(session_completed_trials_afterswap, 1);
                learning_speed = NaN;
                for i = 1:length(mean_completion_rate)
                    if all(mean_completion_rate(i:end) >= accuracy_threshold)
                        learning_speed = i;
                        break;
                    end
                end
                if ~isnan(learning_speed)
                    session_learning_points_swapped = [session_learning_points_swapped, learning_speed];
                end
            end
        end
    end
end


function plotErrorDropRates(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 400];

    % Define subjects
    subjects = {'Bard', 'Sindri', 'Jotun', 'Kyrre'};
    
    % Define error types
    error_types = {'ExplorationError', 'RuleBreakingError', 'DistractorError', 'PerseverativeError'};

    % Define colors and shapes
    colors = [hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A')];
    shapes = {'o', 's', 'd', '^'};

    % Initialize storage variables
    drop_rates = cell(length(subjects), 1);

    % Process data
    for s = 1:length(subjects)
        subject = subjects{s};
        subject_rows = touch_table(strcmp(touch_table.Subject, subject) & touch_table.isRepetition == 0 & touch_table.isSwapped == 0, :);
        sessions = unique(subject_rows.SessionNumber);
        drop_rates{s} = zeros(length(sessions), length(error_types));
        
        for sess = 1:length(sessions)
            session_rows = subject_rows(subject_rows.SessionNumber == sessions(sess), :);
            
            for e = 1:length(error_types)
                session_error_rates = zeros(15, 1);
                
                for trial = 1:15
                    trial_rows = session_rows(session_rows.TrialNumberInBlock == trial, :);
                    total_touches = height(trial_rows);
                    
                    if total_touches > 0
                        session_error_rates(trial) = sum(trial_rows.(error_types{e})) / total_touches;
                    end
                end
                
                % Fit exponential decay function for each session
                trials = 1:15;
                ft = fittype('a*exp(-b*x) + c', 'independent', 'x');
                [curve, ~] = fit(trials', session_error_rates, ft, 'StartPoint', [max(session_error_rates), 0.1, min(session_error_rates)]);
                drop_rates{s}(sess, e) = curve.b;
            end
        end
    end

    % Combine all drop rates
    all_drop_rates = cell2mat(drop_rates);

    % Create figure
    fig = figure;
    fig.Position = canvas_size;

    hold on;

    y_max = 0;

    % Print out mean and 95% CI for each error type
    fprintf('\n Error Type \t Mean \t 95 CI \n');
    for e = 1:length(error_types)
        % Remove outliers
        current_rates = all_drop_rates(:, e);
        [~, TF] = rmoutliers(current_rates);
        clean_rates = current_rates(~TF);
        
        % Calculate mean and 95% CI
        mean_rate = mean(clean_rates);
        ci = tinv(0.975, length(clean_rates)-1) * std(clean_rates) / sqrt(length(clean_rates));
        
        fprintf('%s\t%s\t[%s, %s]\n', error_types{e}, formatNumber(mean_rate), formatNumber(mean_rate-ci), formatNumber(mean_rate+ci));
        
        
        % Add box plot for this error type
        y_max = max([y_max,max(clean_rates)]);
        boxplot_data = clean_rates;
        h = boxplot(boxplot_data, 'Positions', e, 'Colors', colors(e,:), 'Width', 0.3, 'Symbol', '');
        set(h, 'LineWidth', 1.5);

        % Add shapes for each subject
        for s = 1:length(subjects)
            subject_drop_rates = drop_rates{s}(:, e);
            subject_drop_rates = subject_drop_rates(~TF(sum(cellfun(@length, drop_rates(1:s-1)))+1:sum(cellfun(@length, drop_rates(1:s)))));
            swarmchart(repmat(e, length(subject_drop_rates), 1), subject_drop_rates, 36,colors(e,:), 'filled', shapes{s},'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none', 'LineWidth', 1);
        end

        % Add error bar (95% CI)
        errorbar(e, mean_rate, ci, 'k', 'LineWidth', 2);
    end
    
    % Perform Welch's t-test for each pair of conditions
    num_conditions = length(error_types);
    p_values = zeros(num_conditions, num_conditions);
    fprintf('\nPairwise Comparisons (p-values):\n');
    for i = 1:num_conditions
        for j = i+1:num_conditions
            condition1 = all_drop_rates(:,i);
            [~, TF1] = rmoutliers(condition1);
            condition1 = condition1(~TF1);
            condition2 = all_drop_rates(:,j);
            [~, TF2] = rmoutliers(condition2);
            condition2 = condition2(~TF2);

            [~, p_values(i,j)] = ttest2(condition1, condition2, 'Vartype', 'unequal');
           fprintf('%s vs %s: %s\n', error_types{i}, error_types{j}, formatNumber(p_values(i,j)));
        end
    end
    
    % Add significance lines and stars for condition comparisons
    y_start = y_max * 1.1;
    y_step = y_max * 0.1;
    for i = 1:num_conditions
        for j = i+1:num_conditions
            p = p_values(i,j);
            line_shorten = 0.1;
            y = y_start + (j-i-1)*y_step;
            if i == 2 && j == 4
                y = y + y_step;
            elseif i == 1 && j == 4
                y = y + y_step;
            end
            plot([i+line_shorten, j-line_shorten], [y, y], 'k', 'LineWidth', 1.5);
            if p < 0.001
                text((i+j)/2, y+y_step/2, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
            elseif p < 0.01
                text((i+j)/2, y+y_step/2, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
            elseif p < 0.05
                text((i+j)/2, y+y_step/2, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
            else
                text((i+j)/2, y+y_step/2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
            end
        end
    end
    
    % Add significance test for comparison to 0
    fprintf('\nComparison to 0 (p-values):\n');
    for e = 1:length(error_types)
        % Remove outliers
        current_rates = all_drop_rates(:, e);
        [~, TF] = rmoutliers(current_rates);
        clean_rates = current_rates(~TF);

        % Perform one-sample t-test against 0
        [~, p] = ttest(clean_rates);
        fprintf('%s: %s\n', error_types{e}, formatNumber(p));
        
        % Add star and line to the left side of each condition
        if p < 0.001
            text(e-0.4, mean(clean_rates), '***', 'FontSize', 12, 'HorizontalAlignment', 'right');
        elseif p < 0.01
            text(e-0.4, mean(clean_rates), '**', 'FontSize', 12, 'HorizontalAlignment', 'right');
        elseif p < 0.05
            text(e-0.4, mean(clean_rates), '*', 'FontSize', 12, 'HorizontalAlignment', 'right');
        else
            text(e-0.4, mean(clean_rates), 'ns', 'FontSize', 12, 'HorizontalAlignment', 'right');
        end
        
        if p < 0.05
            plot([e-0.25, e-0.25], [mean(clean_rates)*1/2, mean(clean_rates)*3/2], 'k', 'LineWidth', 1.5);
        end
    end

    ylabel('Exponential Drop Rate', 'FontSize', 12);
    title('Exponential Drop Rates by Error Type', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    set(gca, 'XTick', 1:length(error_types), 'XTickLabel', error_types);
    set(gcf,'renderer','Painters');
    % Add legend for subjects
    subject_handles = gobjects(1, length(subjects));
    for s = 1:length(subjects)
        subject_handles(s) = scatter(NaN, NaN, 100, 'k', shapes{s}, 'filled');
    end
    legend_obj = legend(subject_handles, cellfun(@(x) ['Subject ' x(1)], subjects, 'UniformOutput', false), 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';
    
    yline(0,'--', 'HandleVisibility', 'off');
    ylim([-0.15, y_max*1.5]);  % Increased upper limit to accommodate significance lines

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;

	    function str = formatNumber(num)
        if num == 0
            str = '0';
        else
            % Find the exponent
            exponent = floor(log10(abs(num)));
            
            % Normalize the number
            normalized = num / 10^exponent;
            
            % Find the number of decimal places needed
            if abs(normalized) >= 1
                decimal_places = 1;
            else
                decimal_places = 2 - floor(log10(abs(normalized)));
            end
            
            % Create the format string
            format_str = sprintf('%%.%df%%+03d', decimal_places);
            
            % Format the number
            str = sprintf(format_str, normalized, exponent);
        end
    end
end


function plotNewErrorRates(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];
    
    % Define subjects
    subjects = {'Bard', 'Sindri', 'Jotun', 'Kyrre'};
    
    % Define error types
    error_types = {'ExplorationError', 'RuleBreakingError', 'DistractorError', 'PerseverativeError'};
    
    % Define colors and shapes
    colors = [hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A')];
    shapes = {'o', 's', 'd', '^'};
    
    % Initialize storage variables
    all_error_counts = cell(length(error_types), 1);
    all_touch_counts = zeros(length(subjects), 15);
    overall_error_counts = zeros(length(subjects), 15);
    drop_rates = zeros(length(subjects), length(error_types));
    overall_decay_rates = zeros(1, length(error_types));
    decreasing_rates = zeros(1, length(error_types));
    
    % Process data
    for s = 1:length(subjects)
        subject = subjects{s};
        subject_rows = touch_table(strcmp(touch_table.Subject, subject) & touch_table.isRepetition == 0 & touch_table.isSwapped==0, :);
        for trial = 1:15
            trial_rows = subject_rows(subject_rows.TrialNumberInBlock == trial, :);
            total_touches = height(trial_rows);
            if total_touches > 0
                all_touch_counts(s, trial) = total_touches;
                % Count overall errors (any error counts as 1)
                overall_error_counts(s, trial) = sum(any(trial_rows{:, error_types} > 0, 2));
                % Count individual error types
                for e = 1:length(error_types)
                    all_error_counts{e}(s, trial) = sum(trial_rows.(error_types{e}));
                end
            end
        end
    end
    
    % Calculate error rates
    all_error_rates = cellfun(@(x) x ./ all_touch_counts, all_error_counts, 'UniformOutput', false);
    overall_error_rates = overall_error_counts ./ all_touch_counts;
    
    % Calculate overall decreasing rate
    mean_overall_rate = mean(overall_error_rates, 1, 'omitnan');
    overall_decreasing_rate = mean_overall_rate(end) / mean_overall_rate(1);
    overall_percent_decrease = (1 - overall_decreasing_rate) * 100;
    
    % Plotting
    fig = figure;
    fig.Position = canvas_size;
    hold on;
    legend_handles = gobjects(1, length(subjects) + length(error_types));
    
    % Plot overall error rate
%     mean_overall_rate = mean(overall_error_rates, 1, 'omitnan');
%     plot(1:15, mean_overall_rate, '-k', 'LineWidth', 2);
%     legend_handles(end) = plot(NaN, NaN, '-k', 'LineWidth', 2);
    
    for e = 1:length(error_types)
        % Calculate mean error rates
        mean_rate = mean(all_error_rates{e}, 1, 'omitnan');
        trials = 1:15;
        
        % Calculate decreasing rate
        first_trial_value = mean_rate(1);
        last_trial_value = mean_rate(end);
        decreasing_rates(e) = last_trial_value / first_trial_value;
        
        % Plot individual subject data points and fit exponential
        for s = 1:length(subjects)
            scatter(trials, all_error_rates{e}(s, :), 30, shapes{s}, 'MarkerEdgeColor', colors(e,:), 'HandleVisibility', 'off');
            if e == 1
                legend_handles(s) = scatter(NaN, NaN, 50, 'k', shapes{s}, 'MarkerEdgeColor', 'k');
            end
            % Fit exponential decay function for each subject
            ft = fittype('a*exp(-b*x) + c', 'independent', 'x');
            [curve, ~] = fit(trials', all_error_rates{e}(s, :)', ft, 'StartPoint', [max(all_error_rates{e}(s, :)), 0.1, min(all_error_rates{e}(s, :))]);
            
			drop_rates(s, e) = curve.b;
        end
        
        % Plot mean data points
        scatter(trials, mean_rate, 30, colors(e,:), 'filled', 'HandleVisibility', 'off');
        
        % Fit exponential decay function for mean (overall data)
        [curve, ~] = fit(trials', mean_rate', ft, 'StartPoint', [max(mean_rate), 0.1, min(mean_rate)]);
        disp(curve)
		overall_decay_rates(e) = curve.b;
        
        % Plot fitted curve
        legend_handles(length(subjects) + e) = plot(trials, curve(trials), '-', 'Color', colors(e,:), 'LineWidth', 2);
    end
    
    % Customize plot
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Error Rate', 'FontSize', 12);
    title('Error Rates Across Trials for All Subjects', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    ylim([0, max(max(overall_error_rates)) * 1.1]);
    xlim([0, 16]);
    xticks(1:2:15);
    set(gcf,'renderer','Painters')
    % Create legend
    subject_labels = cellfun(@(x) ['Subject ' x(1)], subjects, 'UniformOutput', false);
    legend_labels = [subject_labels, error_types, 'Overall Error Rate'];
    legend_obj = legend(legend_handles, legend_labels, 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
    hold off;
    
    % Print out the overall decay rates and decreasing rates for each error condition
    fprintf('Results for each error condition:\n');
    for e = 1:length(error_types)
        fprintf('%s:\n', error_types{e});
        fprintf('  Overall decay rate: %.4f\n', overall_decay_rates(e));
        fprintf('  Decreasing rate (last/first): %.4f\n', decreasing_rates(e));
        fprintf('  Percent decrease: %.2f%%\n\n', (1 - decreasing_rates(e)) * 100);
    end
    
    % Print out the overall error rate drop
    fprintf('Overall Error Rate:\n');
    fprintf('  Decreasing rate (last/first): %.4f\n', overall_decreasing_rate);
    fprintf('  Percent decrease: %.2f%%\n', overall_percent_decrease);
end

function plotTrialTo80PercentByState(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects and ensure we have 5 states
    subjects = unique(touch_table.Subject);
    numStates = 5;
    maxTrialNumber = max(touch_table.TrialNumberInBlock);

    % Initialize cell array to store crossing points for each session
    allCrossingPoints = cell(length(subjects), 1);

    % Process data for each subject and state
    for s = 1:length(subjects)
        subjectData = touch_table(strcmp(touch_table.Subject, subjects{s}), :);
        sessions = unique(subjectData.SessionNumber);
        sessionCrossingPoints = nan(length(sessions), numStates);

        for sess = 1:length(sessions)
            sessionData = subjectData(subjectData.SessionNumber == sessions(sess), :);
            blocks = unique(sessionData.BlockNumber);

            for state = 1:numStates
                stateData = zeros(maxTrialNumber, 1);
                trialCounts = zeros(maxTrialNumber, 1);

                for b = 1:length(blocks)
                    blockData = sessionData(sessionData.BlockNumber == blocks(b), :);
                    if all(blockData.isSwapped == 1) || all(blockData.isRepetition == 1)
                        continue
                    end

                    for trialNumber = 1:maxTrialNumber
                        currentTrialData = blockData(blockData.TrialNumberInBlock == trialNumber, :);
                        if ~isempty(currentTrialData)
                            trialCounts(trialNumber) = trialCounts(trialNumber) + 1;
                            if max(currentTrialData.CurrentState) >= state
                                stateData(trialNumber) = stateData(trialNumber) + 1;
                            end
                        end
                    end
                end

                % Calculate proportion
                proportion = stateData ./ trialCounts;
                crossingPoint = find(proportion >= 0.7, 1);
                if ~isempty(crossingPoint)
                    sessionCrossingPoints(sess, state) = crossingPoint;
                end
            end
        end
        allCrossingPoints{s} = sessionCrossingPoints;
    end

    % Create the plot
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Define colors for each ordinal position
   colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;

    % Define shapes for each subject
    shapes = {'o', 's', 'd', '^'};
    allStateData = cell(1, numStates);
    
    % Plot swarm for each state
    for state = 1:numStates
        allStatePoints = [];
        subjectIndices = [];
        for s = 1:length(subjects)
            statePoints = allCrossingPoints{s}(:, state);
            allStatePoints = [allStatePoints; statePoints];
            subjectIndices = [subjectIndices; repmat(s, length(statePoints), 1)];
        end

        allStateData{state} = allStatePoints;
       
        % Add shapes for each subject
        for s = 1:length(subjects)
            subjectPoints = allStatePoints(subjectIndices == s);
            swarmchart(repmat(state, length(subjectPoints), 1), subjectPoints, 36, colors(state,:), 'filled', shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', colors(state,:), 'MarkerEdgeAlpha', 0.2);
        end

        % Calculate mean and 95% CI
        meanPoint = nanmean(allStatePoints);
        sePoint = nanstd(allStatePoints) / sqrt(length(allStatePoints));
        ciPoint = 1.96 * sePoint;
		
		fprintf('state: %d, mean: %.2f ± %.2f\n',state,meanPoint,ciPoint);
		
        % Plot mean and CI
        errorbar(state, meanPoint, ciPoint, 'k', 'LineWidth', 2);
        plot(state, meanPoint, 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
    end

    % Add box plot
    boxplot(cell2mat(allStateData), 'Colors', colors, 'Symbol', '', 'Width', 0.3, 'Positions', 1:numStates);
	set(findobj(gca,'type','line'),'linew',2);
    % Add legend for subjects
    legendHandles = gobjects(1, length(subjects));
    for s = 1:length(subjects)
        legendHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 8);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';

    % Customize the plot
    xlabel('Ordinal Position', 'FontSize', 12);
    ylabel('Trial Number to Reach 80% Completion', 'FontSize', 12);
    title('Trial Number to Reach 80% Completion Rate for Each State', 'FontSize', 14);
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
	set(gcf,'renderer','Painters')
    % Set axis limits
    xlim([0.5 numStates+0.5]);
    ylim([0 max(cellfun(@(x) max(x(:)), allCrossingPoints))+1]);

    % Set x-axis ticks
    xticks(1:numStates);
    xticklabels({'1', '2', '3', '4', '5'});

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function plotProportionTrialsByState(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects, sessions, and blocks
    subjects = unique(touch_table.Subject);
    numStates = 5; % Assuming there are 5 states
    maxTrialNumber = max(touch_table.TrialNumberInBlock);
    
    % Initialize matrices to store counts
    trialsReachedStateCount = zeros(maxTrialNumber, numStates, length(subjects));
    totalTrialsPerNumber = zeros(maxTrialNumber, length(subjects));
    
    % Step 1: Data Processing
    % Loop through subjects, sessions, and blocks
    for s = 1:length(subjects)
        subjectData = touch_table(strcmp(touch_table.Subject, subjects{s}), :);
        sessions = unique(subjectData.SessionNumber);
        
        for sess = 1:length(sessions)
            sessionData = subjectData(subjectData.SessionNumber == sessions(sess), :);
            blocks = unique(sessionData.BlockNumber);
            
            for b = 1:length(blocks)
                blockData = sessionData(sessionData.BlockNumber == blocks(b), :);
                
                if any(blockData.isSwapped == 1) || any(blockData.isRepetition == 1)
                    continue
                end
                
                for trialNumber = 1:maxTrialNumber
                    currentTrialData = blockData(blockData.TrialNumberInBlock == trialNumber, :);
                    
                    if ~isempty(currentTrialData)
                        highestState = max(currentTrialData.CurrentState);
                        
                        % Update total trials count for this trial number
                        totalTrialsPerNumber(trialNumber, s) = totalTrialsPerNumber(trialNumber, s) + 1;
                        
                        % Update counts for all states up to and including the highest state
                        for state = 1:numStates
                            if state <= highestState
                                trialsReachedStateCount(trialNumber, state, s) = trialsReachedStateCount(trialNumber, state, s) + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Calculate the proportion of trials reaching each state for each subject
    proportionTrialsReachedState = zeros(maxTrialNumber, numStates, length(subjects));
    for s = 1:length(subjects)
        proportionTrialsReachedState(:,:,s) = trialsReachedStateCount(:,:,s) ./ totalTrialsPerNumber(:,s);
    end
    
    % Calculate average proportion across subjects
    avgProportionTrialsReachedState = mean(proportionTrialsReachedState, 3, 'omitnan');
    
    fig = figure;
    fig.Position = canvas_size;
    hold on;
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255; % Get unique colors for each state
    markers = {'o', 's', 'd', '^'}; % Different markers for each subject
    
    % Define the sigmoid model
    sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', ...
                       'independent', 'x', ...
                       'coefficients', {'a', 'b', 'c'});
    
    % Initialize matrix to store 90% maximum points
    ninetyPercentPoints = zeros(length(subjects), numStates);
    
    % Plot proportions
    legendHandles = gobjects(numStates + length(subjects), 1);
    legendLabels = cell(numStates + length(subjects), 1);
    
    % Initialize cell array to store fit results
    fitResults = cell(length(subjects), numStates);
    
    for state = 1:numStates
        yData = avgProportionTrialsReachedState(:, state);
        xData = (1:maxTrialNumber)';
        
        % Plot data points for each subject
        for s = 1:length(subjects)
            h = plot(xData, proportionTrialsReachedState(:, state, s), markers{s}, 'Color', colors(state, :), 'MarkerSize', 4);
            if state == 1
                legendHandles(numStates + s) = h;
                legendLabels{numStates + s} = ['Subject ', subjects{s}(1)];
            end
            
            % Fit the sigmoid model to the data
            [fitResult, ~] = fit(xData, proportionTrialsReachedState(:, state, s), sigmoidModel, 'StartPoint', [1, 0.1, 7.5]);
            % Store the fit result
            fitResults{s, state} = fitResult;
            
            % Calculate the effective maximum (bounded by 1)
            effectiveMax = min(fitResult.a, 1);
            % Function to find the difference from 90% of effective maximum
            diffFrom90Percent = @(x) fitResult(x) - 0.9 * effectiveMax;
            
            % Check if the function crosses 90% of effective maximum within [1, 15]
            if diffFrom90Percent(1) * diffFrom90Percent(15) > 0
                if diffFrom90Percent(15) < 0
                    % If it doesn't reach 90%, set to 15
                    ninetyPercentPoint = 15;
                else
                    % If it's always above 90%, set to 1
                    ninetyPercentPoint = 1;
                end
            else
                % Find the point where the fitted curve reaches 90% of effective maximum
                ninetyPercentPoint = fzero(diffFrom90Percent, [1, 15]);
            end
            
            % Constrain the point between 1 and 15 (shouldn't be necessary, but just in case)
            ninetyPercentPoint = max(1, min(15, ninetyPercentPoint));
            
            ninetyPercentPoints(s, state) = ninetyPercentPoint;
        end
        
        % Plot average data points
        plot(xData, yData, '.', 'Color', colors(state, :), 'HandleVisibility', 'off', 'MarkerSize', 12, 'LineWidth', 2);
        
        % Fit the sigmoid model to the average data
        [avgFitResult, ~] = fit(xData, yData, sigmoidModel, 'StartPoint', [1, 0.1, 7.5]);
        
        % Plot the fitted sigmoid line
        fittedxData = linspace(1, maxTrialNumber, 100);
        fittedyData = avgFitResult(fittedxData);
        h = plot(fittedxData, fittedyData, '-', 'Color', colors(state, :), 'LineWidth', 1);
        legendHandles(state) = h;
        legendLabels{state} = sprintf('Ordinal Position %d', state);
	end

	
	for s = 1:length(subjects)
        legendHandles(numStates + s) = plot(NaN, NaN, markers{mod(s-1, length(markers))+1}, 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 8);
		subjectname = cell2mat(subjects(s));
		legendLabels{numStates + s} = sprintf('Subject %s',subjectname(1));
	end

    
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';
    set(gcf,'renderer','Painters')
    % Customize the plot
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Proportion of Trials Reaching Each Ordinal Position', 'FontSize', 12);
    title('Proportion of Trials Reaching Each Ordinal Position by Trial Number', 'FontSize', 14);
    legend_obj = legend(legendHandles, legendLabels, 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';
    ylim([0.4, 1.15]);
    xticks(1:2:15);
	yticks(0.4:0.2:1);
    ax = gca;  % Get current axis handle
    ax.FontSize = 12;  % Set font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end


function plotLastFourTrialsAccuracyDots(touch_table)
    % Get unique subjects
    subjects = unique(touch_table.Subject);
    num_subjects = length(subjects);
    colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A'); [107, 37, 110]/255};
    markers = {'o'; 's'; 'd'; '^'};
    chance_levels1 = [1/6, 1/5, 1/4, 1/3, 1/2];
    chance_levels2 = zeros(1,7) +1/6;

    % Initialize figure
    fig = figure('Position', [100 100 300 250]);
    hold on;

    % Prepare data for plotting
    all_data = []; % Store all individual accuracies for swarmplot

    for subj_idx = 1:num_subjects
        current_subject = subjects(subj_idx);
        subject_data = touch_table(strcmp(touch_table.Subject, current_subject) & ...
                                   touch_table.isRepetition == 0 & ...
                                   touch_table.isSwapped == 0, :);

        sessions = unique(subject_data.SessionNumber);

        for sess = sessions'
            sess_data = subject_data(subject_data.SessionNumber == sess, :);
            blocks = unique(sess_data.BlockNumber);

            total_correct = zeros(5, 15);
            total_errors = zeros(5, 15);

            for blk = blocks'
                block_data = sess_data(sess_data.BlockNumber == blk, :);

                for trial = 12:15
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);
                    if max(trial_data.CurrentState) < 5; continue; end

                    for state = 1:5
                        start_idx = find(trial_data.CurrentState == state - 1, 1, 'first');
                        end_idx = find(trial_data.CurrentState == state, 1, 'first');

                        if isempty(start_idx)
                            start_idx = 1;
                        end
                        if isempty(end_idx)
                            end_idx = height(trial_data);
                        end

                        trial_touches = trial_data(start_idx:end_idx - 1, :);

                        correct_touches = 1;
                        exploration_errors = sum(trial_touches.ExplorationError == 1);

                        total_correct(state, trial) = total_correct(state, trial) + correct_touches;
                        total_errors(state, trial) = total_errors(state, trial) + exploration_errors;
                    end
                end
            end

            % Calculate accuracy per state (ordinal position)
            accuracy_per_state = zeros(5, 1);
            total_touches = total_correct + total_errors;
            for state = 1:5
                if sum(total_touches(state, :)) > 0
                    accuracy_per_state(state) = sum(total_correct(state, :)) / sum(total_touches(state, :));
                else
                    accuracy_per_state(state) = NaN;
                end
            end

            % Store data for swarmplot
            all_data = [all_data; (1:5)', accuracy_per_state, repmat(subj_idx, 5, 1)];
        end
    end

    % Swarmplot using swarmchart (MATLAB native function)
    for state = 1:5
        state_indices = all_data(:,1) == state;
        swarmchart(repmat(state, sum(state_indices), 1), ...
                   all_data(state_indices,2), ...
                   30, cell2mat(colors(all_data(state_indices, 3),:)), 'filled', 'MarkerFaceAlpha', 0.5, 'XJitterWidth',0.5, 'HandleVisibility','off');
    end

    % Plot mean accuracies with error bars (95% CI)
    mean_accuracies = arrayfun(@(state) nanmean(all_data(all_data(:,1) == state,2)), 1:5);
    ci_accuracies = arrayfun(@(state) nanstd(all_data(all_data(:,1) == state,2)) / sqrt(sum(~isnan(all_data(all_data(:,1) == state,2)))) * tinv(0.975, sum(~isnan(all_data(all_data(:,1) == state,2))) - 1), ...
                             1:5);

    errorbar(1:5, mean_accuracies, ci_accuracies, 'k', 'LineWidth', 2, 'CapSize',10, 'MarkerSize',2, 'DisplayName','Average Accuracy');
	

	% Plot chance level
	plot(1:5, chance_levels1, '--','Color','[0.1,0.1,0.1]', 'DisplayName', 'Chance Performance w/o Choosing Previous Objects');
    plot(0:6, chance_levels2, 'k--', 'DisplayName', 'Chance Performance');
    
	% Add legend for subjects
    legendHandles = gobjects(1, num_subjects+3);
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, 'o', 'MarkerFaceColor', colors{s}, ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
	end
	legendHandles(num_subjects+1) = errorbar(NaN, NaN, 'k', 'LineWidth', 2, 'CapSize',10, 'MarkerSize',2, 'DisplayName','Average Accuracy');
	legendHandles(num_subjects+2) = plot(NaN, NaN, '--','Color','[0.1,0.1,0.1]', 'DisplayName', 'Chance Performance w/o Choosing Previous Objects');
	legendHandles(num_subjects+3) = plot(NaN, NaN,  'k--', 'DisplayName', 'Chance Performance');
	result = [cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false); 'Average Accuracy';'Chance Performance w/o Choosing Previous Objects';'Chance Performance'];
    % result = vertcat(result{:});
	legend_obj = legend(legendHandles, result, ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Customize plot
    xlabel('Ordinal Position');
    ylabel('Accuracy');
    title('Average Accuracy for Trials (12-15)');
    xlim([0.5 5.5]);
    ylim([0 1.1]);
    xticks(1:5);
    legend_obj = legend('Location', 'eastoutside');
	legend_obj.EdgeColor = 'none';
    set(gca, 'TickDir', 'out', 'FontSize', 12);
	   
    hold off;


     % Adjust figure size to accommodate legend if it exists
    canvas_size = [100, 100, 300, 250];
        % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    new_fig_size = [canvas_size(1), canvas_size(2), ...
                    canvas_size(3) + canvas_size(3) * legend_pos(3), canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    % Ensure the renderer is set to 'Painters'
    set(gcf, 'renderer', 'Painters');

	   % Print additional information
    fprintf('\n--- Summary Information ---\n');

    % Average accuracy at each ordinal position for each subject
    for subj_idx = 1:num_subjects
        current_subject = subjects(subj_idx);
        subject_data = all_data(all_data(:,3) == subj_idx, :);
        fprintf('Subject %s:\n', current_subject{1});
        for state = 1:5
            state_data = subject_data(subject_data(:,1) == state, 2);
            if ~isempty(state_data)
                mean_accuracy = nanmean(state_data);
                ci_accuracy = nanstd(state_data) / sqrt(sum(~isnan(state_data))) * ...
                              tinv(0.975, sum(~isnan(state_data)) - 1);
                fprintf('  Ordinal Position %d: %.2f ± %.3f\n', state, mean_accuracy, ci_accuracy);
            else
                fprintf('  Ordinal Position %d: No data\n', state);
            end
        end

        % Number of valid sessions for the subject
        valid_sessions = unique(touch_table.SessionNumber(strcmp(touch_table.Subject, current_subject) & ...
                                                          touch_table.isRepetition == 0 & ...
                                                          touch_table.isSwapped == 0));
        fprintf('  Number of valid sessions: %d\n\n', numel(valid_sessions));
    end

    % Overall average accuracy and 95% CI at each ordinal position
    fprintf('Overall Average Accuracy:\n');
    for state = 1:5
        state_data = all_data(all_data(:,1) == state, 2);
        if ~isempty(state_data)
            mean_accuracy = nanmean(state_data);
            ci_accuracy = nanstd(state_data) / sqrt(sum(~isnan(state_data))) * ...
                          tinv(0.975, sum(~isnan(state_data)) - 1);
            fprintf('  Ordinal Position %d: %.2f ± %.3f\n', state, mean_accuracy, ci_accuracy);
        else
            fprintf('  Ordinal Position %d: No data\n', state);
        end
	end

end





function accuracy_matrices = generateTouchAccuracyPlots(touch_table)
% Get unique subjects
subjects = unique(touch_table.Subject);

% Initialize output matrix
accuracy_matrices = zeros(5, 15, length(subjects));
% Initialize figure
	figure('Position', [100 100 400 240*length(subjects)]);

% Process each subject
for subj_idx = 1:4
	current_subject = subjects(subj_idx);
	subject_data = touch_table(strcmp(touch_table.Subject, current_subject) & touch_table.isRepetition == 0 & touch_table.isSwapped == 0, :);

	% Initialize accuracy matrix for current subject
	accuracy_matrix = zeros(5, 15);

	% Get unique sessions and blocks for current subject
	sessions = unique(subject_data.SessionNumber);

	% Initialize counters for accumulating touches across sessions/blocks
	total_correct = zeros(5, 15);
	total_errors = zeros(5, 15);

	% Loop through sessions
	for sess = sessions'
		sess_data = subject_data(subject_data.SessionNumber == sess, :);
		blocks = unique(sess_data.BlockNumber);

		% Loop through blocks
		for blk = blocks'
			block_data = sess_data(sess_data.BlockNumber == blk, :);

			% Calculate accuracy for each state and trial
			for trial = 1:15
				% Get relevant touches
				trial_data = block_data(block_data.TrialNumberInBlock == trial,:);
				if max(trial_data.CurrentState) < 5; continue; end

				for state = 1:5
					start_idx = find(trial_data.CurrentState == state-1, 1, 'first');
            		end_idx = find(trial_data.CurrentState == state, 1, 'first');
            		
            		% Handle edge cases
            		if isempty(start_idx)
                		start_idx = 1;
            		end
            		if isempty(end_idx)
                		end_idx = height(trial_data);
					end
					trial_touches = trial_data(start_idx:end_idx-1, :);

					% Count correct touches and exploration errors
					correct_touches = 1;
					exploration_errors = sum(trial_touches.ExplorationError == 1);

					% Accumulate counts
					total_correct(state, trial) = total_correct(state, trial) + correct_touches;
					total_errors(state, trial) = total_errors(state, trial) + exploration_errors;
				end
			end
		end
	end

	% Calculate final accuracy matrix
	total_touches = total_correct + total_errors;
	accuracy_matrix = zeros(5, 15);
	for state = 1:5
		for trial = 1:15
			if total_touches(state, trial) > 0
				accuracy_matrix(state, trial) = total_correct(state, trial) / total_touches(state, trial);
			else
				accuracy_matrix(state, trial) = 0;
			end
		end
	end

	% Store accuracy matrix
	accuracy_matrices(:,:,subj_idx) = accuracy_matrix;

	
	% Create subplot for current subject
	subplot(length(subjects), 1, subj_idx);

	% Plot heatmap
	imagesc(accuracy_matrix, [0.3,0.85]);

	colorbar;

	% Set labels and title
	xlabel('Trial Number');
	ylabel('Ordinal Position');
	title(['Subject ' char(current_subject(1)) '']);
	
	set(gcf, 'renderer', 'Painters');
	% Set axis properties
	xticks(1:15);
	yticks(1:5);
	axis xy;
end

end