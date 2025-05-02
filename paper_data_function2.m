%% Data Loading
cd /Users/wenxuan/Documents/MATLAB/WWW_Project/ 
% 
% filename = "/Users/wenxuan/Documents/MATLAB/WWW_Project/combined_touch_data.csv";
% opts = detectImportOptions(filename);
% touch_table = readtable(filename, opts);
% touch_table = fixRuleBreakingErrors(touch_table);
% % [dataWithPerf, avgSearchTimesByBlock] = classifyPerformance(touch_table, 0.8);

% filename = "/Users/wenxuan/Documents/MATLAB/WWW_Project/combined_block_data.csv";
% opts = detectImportOptions(filename);
% block_table = readtable(filename, opts);
% [block_table] = calculateAnticipatedSwapping(touch_table, block_table);
% [block_table] = calculateOrdinalObjectIdentityInference(touch_table, block_table);

% load('/Users/wenxuan/Documents/MATLAB/WWW_Project/WM_MAT01_METRICS/WWW_WM_WWW_resultsWM.mat');

%% Plotting

% height(block_table(strcmp(block_table.Subject, 'Kyrre') & block_table.BlockLevelAccuracy >= 0.8 , :))/ height(block_table(strcmp(block_table.Subject, 'Kyrre'), :))
%% Figure 1

% Figure 1C
% plotLastFourTrialsAccuracyDots(touch_table);

% Figure 1D
% calculateSessionLearningSpeed(touch_table, block_table, 0.8);

% Figure 1E
% calculateSessionBasedLearningSpeed(touch_table, block_table, 0.8);


%% Figure 2

% Figure 2B
% calculateLearningPerformanceSwap(touch_table);

% Figure 2C
% plotLearningPerformanceSwapBoxPlot(touch_table);

% Figure 2D
% plotSessionWiseTransitionProbability(touch_table);

% Figure 2E
% plotNextErrorFrequency(touch_table, 0.8);

% Figure 2F
% plotSessionWiseAnticipatedSwappingTransitionProbability(touch_table);

% Figure 2G
% plotSessionWiseSecondToThirdTransitionProbability(touch_table);

%% Figure 3

% Figure 3ABC-DEF
% plotTransitionProbabilitiesWithAccuracyFilter(touch_table)

% Figure 3G
% plotTransitionProbabilitiesFirstTrial(touch_table)

% Figure 3H
% plotTransitionDifferencesByTrial(touch_table)



%% Figure 4
% Figure 4A
% plotTrialCompletionRatesRepeat(touch_table);
% Figure 4B
% plotAnticipatedSwappingThreeConditions(block_table)
% Figure 4C
% plotErrorInferenceInRepeat(block_table);

% Figure 4D
% plotInitialPerformaceVsAnticipatedSwappingTwoConditions(block_table)
% Figure 4E
% plotBlockCompletionVsAnticipatedSwapping(block_table);
% Figure 4H
% plotAnticipatedSwappingOverSession(block_table);

% Figure 4F
% plotInitialPerformaceVsOrdinalObjectIdentityTwoConditions(block_table);
% Figure 4J
% plotBlockCompletionVsErrorCorrection(block_table);
% Figure 4I
% plotOrdinalPositionInferenceOverSession(block_table);


%% Figure 5
% Figure 5D
% plot_block_accuracy_difference(block_table);
% % Figure 5E
% plot_ordinal_identity_inference_difference(block_table);
% % Figure 5F
% plot_anticipated_swapping_difference(block_table);



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
% plotAvgSearchTimeByOrdinalPositionLearning(touch_table);
% Figure S2B
% plotAvgSearchTimeByStateAndSwap(touch_table);

% Figure S3A
% plotFirstTrialPosition80Repeat(touch_table);
% Figure S3B
% plotAvgSearchTimeCorrectTouchesRepeatConditions(touch_table, block_table)


% Figure S5A
% plotErrorDecreasingRateSwap(touch_table);
% Figure S5B
% plotErrorDecreasingRateMemory(touch_table);
% Figure S5C
% plotDistractorErrorFrequency(touch_table);

%% Sup5
% plotErrorDecreasingRate(touch_table);

% plotAvgSearchTimeByTrialAndPerformance(touch_table, block_table, 0.8);
% plotDistractorErrorRates(touch_table)
% plotDistractorErrorProportionsByState(touch_table);
% plotOrdinalPositionInference(block_table);
 


% plotContextConditionAccuracy(block_table);
% plotBackgroundErrorInference(block_table);

% plotASTransitionProbPerformanceSeperated(touch_table);

% plotCompletionRateVsBlocks(block_table);
% plotAvgSearchTimeByOrdinalPosition(touch_table);



function plotTransitionDifferencesByTrial(touch_table)
    canvas_size = [100, 100, 600, 250];
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);
    colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#808080')};
    numTrials = 15;
    numStates = 6; % All 6 transition states
    
    % Initialize array to store all session data
    all_data = [];
    
    % Loop through each subject
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subject_data.SessionNumber);
        
        % Process each session
        for sess_idx = 1:length(sessions)
            session = sessions(sess_idx);
            session_data = subject_data(subject_data.SessionNumber == session, :);
            blocks = unique(session_data.BlockNumber);
            
            % Initialize transitions for this session
            session_transitions = zeros(numStates, numTrials);
            
            % Process each trial number
            for trial_num = 1:numTrials
                transitions_count = zeros(1, numStates);
                total_transitions = 0;
                
                % Process each block
                for block = blocks'
                    block_data = session_data(session_data.BlockNumber == block, :);
                    trial_data = block_data(block_data.TrialNumberInBlock == trial_num, :);
                    
                    % Check for swap condition
                    if ~isempty(trial_data) && trial_data.isSwapped(1) == 1 && trial_data.isRepetition(1) == 0
                        for touch = 1:height(trial_data) - 1
                            if trial_data.CurrentState(touch) == 1 && ...
                                (strcmpi(trial_data.TouchCategory{touch}, 'correct') || ...
                                strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
                                
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data)
                                    nextPos = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextPos == 4
                                        transitions_count(2) = transitions_count(2) + 1;
                                    elseif nextPos == 3
                                        transitions_count(3) = transitions_count(3) + 1;
                                    elseif nextPos == 2
                                        transitions_count(4) = transitions_count(4) + 1;
                                    elseif nextPos == 1
                                        transitions_count(1) = transitions_count(1) + 1;
                                    elseif nextPos == 5
                                        transitions_count(5) = transitions_count(5) + 1;
                                    elseif nextPos == 0
                                        transitions_count(6) = transitions_count(6) + 1;
                                    end
                                    total_transitions = total_transitions + 1;
                                    break;
                                end
                            end
                        end
                    end
                end
                
                % Calculate probabilities if transitions exist
                if total_transitions > 0
                    session_transitions(:, trial_num) = transitions_count / total_transitions;
                end
            end
            
            % Store session data with subject identifier
            all_data = cat(3, all_data, session_transitions);
        end
    end
    
    % Calculate mean transitions across all sessions
    mean_transitions = mean(all_data, 3);
    
    % Calculate differences
    D_B_diff = mean_transitions(4, :) - mean_transitions(2, :);
    D_C_diff = mean_transitions(4, :) - mean_transitions(3, :);
    C_B_diff = mean_transitions(3, :) - mean_transitions(2, :);
    
    % Calculate standard errors using all sessions
    se_D_B = std(squeeze(all_data(4, :, :) - all_data(2, :, :)), [], 2) / sqrt(size(all_data, 3));
    se_D_C = std(squeeze(all_data(4, :, :) - all_data(3, :, :)), [], 2) / sqrt(size(all_data, 3));
    se_C_B = std(squeeze(all_data(3, :, :) - all_data(2, :, :)), [], 2) / sqrt(size(all_data, 3));

    % Create figure
    figure('Position', canvas_size);
    hold on;
    
    errorbar(1:numTrials, D_B_diff, se_D_B', '-', 'Color', colors{1}, 'LineWidth', 2);
    errorbar(1:numTrials, D_C_diff, se_D_C', '-', 'Color', colors{2}, 'LineWidth', 2);
    errorbar(1:numTrials, C_B_diff, se_D_C', '-', 'Color', colors{3}, 'LineWidth', 2);

    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Difference in Choice Probability', 'FontSize', 12);
    title('D-B and D-C Differences Across Trials', 'FontSize', 14);
    legend('D-B', 'D-C','C-B', 'Location', 'eastoutside');
    xlim([0.5 numTrials+0.5]);
    yline(0, '--k', HandleVisibility='off');
    grid off;
    box off;
    set(gca, 'FontSize', 12, 'TickDir', 'out');
    
    hold off;
end


function plotTransitionProbabilitiesWithAccuracyFilter(touch_table)
    % Define canvas size and basic parameters
    canvas_size = [100, 100, 300, 250];
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);
    
    % Define visual elements
    colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); 
             hex2rgb('#4E937A'); [107, 37, 110]/255};
    markers = {'o'; 's'; 'd'; '^'};
    
    % Initialize data storage
    allSessionProbs_Swap0_LastTrial = [];
    allSessionProbs_Swap1_FirstTrial = [];
    subjectData_Swap0_LastTrial = cell(numSubjects, 1);
    subjectData_Swap1_FirstTrial = cell(numSubjects, 1);
    subjectData_Difference = cell(numSubjects, 1);
    
    % Process each subject
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subject_data.SessionNumber);
        sessionProbs_Swap0_LastTrial = [];
        sessionProbs_Swap1_FirstTrial = [];
        
        % Process each session
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            blocks = unique(session_data.BlockNumber);
            
           
            
    		% Initialize block inclusion map
    		includeBlock = zeros(size(blocks));
			includeBlockBad = ones(size(blocks));
            for i = 1:length(blocks)
        		block = blocks(i);
        		if mod(block, 2) == 1  % Odd Block
            		block_data = session_data(session_data.BlockNumber == block, :);
            		unique_trials = unique(block_data.TrialNumberInBlock);
            		numTrials = numel(unique_trials);
            		
					total_number_errors = 0;
        			for trial = 1:numTrials
            			trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
						total_number_errors = total_number_errors + sum(contains(trial_data.TouchCategory, 'Erro'));
					end
					block_level_error_rate = total_number_errors/15;
            		 
            		if block_level_error_rate > 6.6% or > 6.6
						includeBlock(i) = 1;
                		includeBlock(i + 1) = 1;  % Exclude next block if it exists
					end
        		end
    		end
            
            % Calculate transition probabilities
            transitionProbs_Swap0_LastTrial = zeros(1, numStates);
            transitionProbs_Swap1_FirstTrial = zeros(1, numStates);
            totalTrials_Swap0_LastTrial = 0;
            totalTrials_Swap1_FirstTrial = 0;
            
            % Process blocks
            for i = 1:length(blocks)
				block = blocks(i);
				if ~includeBlock(i)
					continue
				end
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);
                
                % Process last trial of non-swapped blocks
                if ~isempty(trials)
                    trial_data_last = block_data(block_data.TrialNumberInBlock == max(trials), :);
                    if trial_data_last.isSwapped(1) == 0 && trial_data_last.isRepetition(1) == 0
                        [transitionProbs_Swap0_LastTrial, totalTrials_Swap0_LastTrial] = ...
                            processTrialData(trial_data_last, transitionProbs_Swap0_LastTrial, ...
                            totalTrials_Swap0_LastTrial, false);
                    end
                    
                    % Process first trial of swapped blocks
                    trial_data_first = block_data(block_data.TrialNumberInBlock == 1, :);
                    if trial_data_first.isSwapped(1) == 1 && trial_data_first.isRepetition(1) == 0
                        [transitionProbs_Swap1_FirstTrial, totalTrials_Swap1_FirstTrial] = ...
                            processTrialData(trial_data_first, transitionProbs_Swap1_FirstTrial, ...
                            totalTrials_Swap1_FirstTrial, true);
                    end
                end
            end
            
            % Normalize and store probabilities
            if totalTrials_Swap0_LastTrial > 0
                sessionProbs_Swap0_LastTrial(end+1, :) = ...
                    transitionProbs_Swap0_LastTrial / totalTrials_Swap0_LastTrial;
            end
            if totalTrials_Swap1_FirstTrial > 0
                sessionProbs_Swap1_FirstTrial(end+1, :) = ...
                    transitionProbs_Swap1_FirstTrial / totalTrials_Swap1_FirstTrial;
            end
        end
        
        % Store subject data
        subjectData_Swap0_LastTrial{s} = sessionProbs_Swap0_LastTrial;
        subjectData_Swap1_FirstTrial{s} = sessionProbs_Swap1_FirstTrial;
        subjectData_Difference{s} = sessionProbs_Swap1_FirstTrial - sessionProbs_Swap0_LastTrial;
	end

	meanSwap0_LastTrial = zeros(1, numStates);
seSwap0_LastTrial = zeros(1, numStates);
meanSwap1_FirstTrial = zeros(1, numStates);
seSwap1_FirstTrial = zeros(1, numStates);
meanDifference = zeros(1, numStates);
seDifference = zeros(1, numStates);


% Calculate for each state
for state = 1:numStates
    % Collect all data points for this state
    swap0_data = [];
    swap1_data = [];
    diff_data = [];
    
    for s = 1:numSubjects
        if ~isempty(subjectData_Swap0_LastTrial{s})
            swap0_data = [swap0_data; subjectData_Swap0_LastTrial{s}(:,state)];
            swap1_data = [swap1_data; subjectData_Swap1_FirstTrial{s}(:,state)];
            diff_data = [diff_data; subjectData_Difference{s}(:,state)];
        end
    end
    
    % Calculate means
    meanSwap0_LastTrial(state) = mean(swap0_data);
    meanSwap1_FirstTrial(state) = mean(swap1_data);
    meanDifference(state) = mean(diff_data);
    
    % Calculate standard errors
    seSwap0_LastTrial(state) = std(swap0_data) / sqrt(length(swap0_data));
    seSwap1_FirstTrial(state) = std(swap1_data) / sqrt(length(swap1_data));
    seDifference(state) = std(diff_data) / sqrt(length(diff_data));
end
    
    % Create plots using the existing createFigure function
    createFigure(meanSwap0_LastTrial, seSwap0_LastTrial, canvas_size, colors, markers, ...
        transitionStates, 'Swap=0 & Last Trial', ...
        subjectData_Swap0_LastTrial, subjects);
    createFigure(meanSwap1_FirstTrial, seSwap1_FirstTrial, canvas_size, colors, markers, ...
        transitionStates, 'Swap=1 & First Trial', ...
        subjectData_Swap1_FirstTrial, subjects);
    createFigure(meanDifference, seDifference, canvas_size, colors, markers, ...
        transitionStates, 'Difference', subjectData_Difference, subjects);
end

function [transitionProbs, totalTrials] = processTrialData(trial_data, transitionProbs, ...
    totalTrials, isSwapped)
    for touch = 1:height(trial_data) - 1
        if trial_data.CurrentState(touch) == 1 && ...
                (strcmpi(trial_data.TouchCategory{touch}, 'correct') || ...
                strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
            nextTouchIndex = touch + 1;
            if nextTouchIndex <= height(trial_data)
                nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                if nextErrorObjectPosition == 0
                    nextErrorObjectPosition = 6;
                elseif nextErrorObjectPosition == 2 && isSwapped
                    nextErrorObjectPosition = 4;
                elseif nextErrorObjectPosition == 4 && isSwapped
                    nextErrorObjectPosition = 2;
                end
                transitionProbs(nextErrorObjectPosition) = ...
                    transitionProbs(nextErrorObjectPosition) + 1;
                totalTrials = totalTrials + 1;
                break;
            end
        end
    end
end




function plotTransitionProbabilitiesFirstTrial(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Define colors for subjects
    colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A'); [107, 37, 110]/255};
	markers = {'o'; 's'; 'd'; '^'};


    % Initialize arrays to store session-wise data
    allSessionProbs_Swap0_LastTrial = [];
    allSessionProbs_Swap1_FirstTrial = [];
    subjectData_Swap0_LastTrial = cell(numSubjects, 1);
    subjectData_Swap1_FirstTrial = cell(numSubjects, 1);
	subjectData_Difference = cell(numSubjects, 1);

    % Loop through each subject and calculate probabilities
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);

        % Aggregate data across all sessions for this subject
        sessions = unique(subject_data.SessionNumber);
        sessionProbs_Swap0_LastTrial = [];
        sessionProbs_Swap1_FirstTrial = [];

        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            blocks = unique(session_data.BlockNumber);

            % Initialize variables for this session
            transitionProbs_Swap0_LastTrial = zeros(1, numStates);
            transitionProbs_Swap1_FirstTrial = zeros(1, numStates);
            totalTrials_Swap0_LastTrial = 0;
            totalTrials_Swap1_FirstTrial = 0;

            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);

                % Collect data for Swap == 0 and last trial
                if ~isempty(trials) 
                    trial_data_last = block_data(block_data.TrialNumberInBlock == max(trials), :);
                    if trial_data_last.isSwapped(1) == 0 && trial_data_last.isRepetition(1) == 0
                        for touch = 1:height(trial_data_last) - 1
                            if trial_data_last.CurrentState(touch) == 1 && ...
                                    (strcmpi(trial_data_last.TouchCategory{touch}, 'correct') || ...
                                    strcmpi(trial_data_last.TouchCategory{touch}, 'correctSelection'))
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data_last)
                                    nextErrorObjectPosition = trial_data_last.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 2;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 4;
                                    end
                                    transitionProbs_Swap0_LastTrial(nextErrorObjectPosition) = ...
                                        transitionProbs_Swap0_LastTrial(nextErrorObjectPosition) + 1;
                                    totalTrials_Swap0_LastTrial = totalTrials_Swap0_LastTrial + 1;
                                    break;
                                end
                            end
                        end
                    end
                end

                % Collect data for Swap == 1 and first trial
                trial_data_first = block_data(block_data.TrialNumberInBlock == 7, :);
                if trial_data_first.isSwapped(1) == 1 && trial_data_first.isRepetition(1) == 0
                    for touch = 1:height(trial_data_first) - 1
                        if trial_data_first.CurrentState(touch) == 1 && ...
                                (strcmpi(trial_data_first.TouchCategory{touch}, 'correct') || ...
                                strcmpi(trial_data_first.TouchCategory{touch}, 'correctSelection'))
                            nextTouchIndex = touch + 1;
                            if nextTouchIndex <= height(trial_data_first)
                                nextErrorObjectPosition = trial_data_first.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                if nextErrorObjectPosition == 0
                                    nextErrorObjectPosition = 6;
                                elseif nextErrorObjectPosition == 2
                                    nextErrorObjectPosition = 4;
                                elseif nextErrorObjectPosition == 4
                                    nextErrorObjectPosition = 2;
                                end
                                transitionProbs_Swap1_FirstTrial(nextErrorObjectPosition) = ...
                                    transitionProbs_Swap1_FirstTrial(nextErrorObjectPosition) + 1;
                                totalTrials_Swap1_FirstTrial = totalTrials_Swap1_FirstTrial + 1;
                                break;
                            end
                        end
                    end
                end

            end

            % Normalize probabilities and store session-wise data
            if totalTrials_Swap0_LastTrial > 0
                sessionProbs_Swap0_LastTrial(end+1, :) = transitionProbs_Swap0_LastTrial / totalTrials_Swap0_LastTrial; 
            end

            if totalTrials_Swap1_FirstTrial > 0
                sessionProbs_Swap1_FirstTrial(end+1, :) = transitionProbs_Swap1_FirstTrial / totalTrials_Swap1_FirstTrial; 
            end

        end

        % Store all sessions' data across subjects for later averaging and plotting
        allSessionProbs_Swap0_LastTrial(end+1:end+size(sessionProbs_Swap0_LastTrial,1), :) = sessionProbs_Swap0_LastTrial; 
        allSessionProbs_Swap1_FirstTrial(end+1:end+size(sessionProbs_Swap1_FirstTrial,1), :) = sessionProbs_Swap1_FirstTrial; 
        
        % Store individual subject data for plotting later
        subjectData_Swap0_LastTrial{s} = sessionProbs_Swap0_LastTrial;
        subjectData_Swap1_FirstTrial{s} = sessionProbs_Swap1_FirstTrial;
		subjectData_Difference{s} = sessionProbs_Swap1_FirstTrial - sessionProbs_Swap0_LastTrial;
    end

    % Calculate averages and SEs for both conditions and their difference (Figure 3)
    meanSwap0_LastTrial = mean(allSessionProbs_Swap0_LastTrial, 1);
    seSwap0_LastTrial = std(allSessionProbs_Swap0_LastTrial, [], 1) / sqrt(size(allSessionProbs_Swap0_LastTrial, 1));

    meanSwap1_FirstTrial = mean(allSessionProbs_Swap1_FirstTrial, 1);
    seSwap1_FirstTrial = std(allSessionProbs_Swap1_FirstTrial, [], 1) / sqrt(size(allSessionProbs_Swap1_FirstTrial, 1));

    meanDifference = mean(allSessionProbs_Swap1_FirstTrial - allSessionProbs_Swap0_LastTrial, 1);
    seDifference = std(allSessionProbs_Swap1_FirstTrial - allSessionProbs_Swap0_LastTrial, [], 1) / ...
                   sqrt(size(allSessionProbs_Swap0_LastTrial, 1));

	disp(meanSwap0_LastTrial);
	disp(seSwap0_LastTrial);

	disp(meanSwap1_FirstTrial);
	disp(seSwap1_FirstTrial);

	disp(meanDifference);
	disp(seDifference);


    % Plot Figure (Swap == 0 and Last Trial)
    createFigure(meanSwap0_LastTrial, seSwap0_LastTrial, canvas_size, colors, markers, ...
                 transitionStates, 'Swap=0 & Last Trial', subjectData_Swap0_LastTrial,subjects);
    % Plot Figure (Swap == 1 and First Trial)
    createFigure(meanSwap1_FirstTrial, seSwap1_FirstTrial, canvas_size, colors, markers,...
                 transitionStates, 'Swap=1 & First Trial', subjectData_Swap1_FirstTrial,subjects);
    % Plot Figure (Difference)
    createFigure(meanDifference, seDifference, canvas_size, colors, markers, ...
                 transitionStates, 'Difference (Swap=1 First - Swap=0 Last)', subjectData_Difference,subjects);
end

function createFigure(meanData, seData, canvas_size, colors, markers, xtickLabels, titleText, individualData, subjects)
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot individual subject data with error bars and LineWidth=1
    if ~isempty(individualData)
        for iSubject = 1:length(individualData)
            subjectData = individualData{iSubject};
            if ~isempty(subjectData)
                meanSubject = mean(subjectData, 1);
                seSubject = std(subjectData, [], 1) / sqrt(size(subjectData, 1));
                % errorbar(1:length(meanSubject), meanSubject, seSubject, 'Color', colors{mod(iSubject-1, length(colors))+1}, 'Marker',markers{mod(iSubject-1, length(markers))+1}, ...
                %          'LineWidth', 1, 'LineStyle', '-');
				 plot(1:length(meanSubject), meanSubject, 'Color', colors{mod(iSubject-1, length(colors))+1}, 'Marker',markers{mod(iSubject-1, length(markers))+1}, ...
                         'LineWidth', 1, 'LineStyle', '-');
            end
        end
    end

    % Plot overall average data with thicker lines
    errorbar(1:length(meanData), meanData, seData, '-k', 'LineWidth', 2);
	disp(meanData);
	disp(seData);

    % Perform Welch's t-tests
    allDataB = [];
    allDataC = [];
    allDataD = [];
    
    for iSubject = 1:length(individualData)
        subjectData = individualData{iSubject};
        if ~isempty(subjectData)
            allDataB = [allDataB; subjectData(:,2)]; % B is index 2
            allDataC = [allDataC; subjectData(:,3)]; % C is index 3
            allDataD = [allDataD; subjectData(:,4)]; % D is index 4
        end
    end
    
    % Perform statistical tests
    [~, pBC] = performWelchTTest(allDataB, allDataC);
    [~, pBD] = performWelchTTest(allDataB, allDataD);
	[~, pCD] = performWelchTTest(allDataC, allDataD);
    
    % Add significance bars
    maxY = max(meanData + seData) + 0.1;
	minY = min(meanData + seData);
    addSignificanceBars(gca, 2, 3, maxY + 0.2, pBC);      % B vs C
    addSignificanceBars(gca, 2, 4, maxY + 0.1, pBD); % B vs D
	addSignificanceBars(gca, 3, 4, maxY + 0.2, pCD); % C vs D
	
	fprintf('B-C p-value: %.10f\n', pBC);
    fprintf('B-D p-value: %.10f\n', pBD);
	fprintf('C-D p-value: %.10f\n', pCD);
    % Adjust y-axis limit to show significance bars
	if minY < 0
		ylim([-0.5, 0.65]);
		yline(0,'--k');
	else
		% ylim([0, maxY + 0.2]);
		ylim([0, 0.9]);
	end

    ylabel('Choice Prob.', 'FontSize', 12);
    xlabel('Object', 'FontSize', 12);
    title(titleText, 'FontSize', 14);

    xticks(1:length(xtickLabels));
    xticklabels(xtickLabels);
    
    xlim([0.5 length(meanData)+0.5]);
    
    set(gca,'FontSize',12,'TickDir','out');
    
    box off;

	 % Add legend for subjects
    legendHandles = gobjects(1, length(subjects));
    for s = 1:length(subjects)
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', colors{s}, ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
	set(gcf,'renderer','Painters')
    hold off;

end


% Add these helper functions at the end of the file
function [h, p] = performWelchTTest(data1, data2)
    [h, p] = ttest2(data1, data2, 'Vartype', 'unequal');
end

function addSignificanceBars(ax, x1, x2, y, p)
    hold(ax, 'on');
    % Draw the bar
    plot(ax, [x1 x2], [y y], '-k', 'LineWidth', 1);
    
    % Add significance stars
    if p < 0.001
        sigText = '***';
    elseif p < 0.01
        sigText = '**';
    elseif p < 0.05
        sigText = '*';
    else
        sigText = 'ns';
    end
    
    text(ax, (x1+x2)/2, y+0.03, sigText, 'HorizontalAlignment', 'center', 'FontSize', 15);
end




function plotASTransitionProbPerformanceSeperated(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Define markers and colors for subjects
    markers = {'o', 's', 'd', '^', 'v'};
    colors = [39, 93, 44;
			58, 112, 175;
			107, 37, 110;
			199, 54, 55;
			233, 166, 64;
			247, 211, 76]/255;
	
    % Initialize array to store session-wise transition probabilities
allSessionProbs_Swapped = cell(numSubjects, 1);
allSessionProbs_Swapped_Good = cell(numSubjects, 1);
allSessionProbs_Swapped_Bad = cell(numSubjects, 1);
accuracy_threshold = 0.8; % Set your desired accuracy threshold

% Calculate transition probabilities for each subject and session
for s = 1:numSubjects
    subject = subjects{s};
    subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
    sessions = unique(subject_data.SessionNumber);
    subjectProbs = [];
	subjectProbsGood = [];
	subjectProbsBad = [];

    for session = sessions'
        session_data = subject_data(subject_data.SessionNumber == session, :);
        blocks = unique(session_data.BlockNumber);
        session_completed_trials = [];
        
        % Calculate learning speed for the session
        for blk = 1:numel(blocks)
            block_data = session_data(session_data.BlockNumber == blocks(blk), :);
            % Ensure the block is not swapped and not a repetition
            if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                unique_trials = unique(block_data.TrialNumberInBlock);
                numTrials = numel(unique_trials);
                completed_trials = zeros(1, numTrials);
                for trial = 1:numTrials
                    trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                    completed_trials(trial) = any(trial_data.isRewarded == 1);
                end
                % Pad completed trials to 15 if necessary
                padded_completed_trials = nan(1, 15);
                padded_completed_trials(1:length(completed_trials)) = completed_trials;
                session_completed_trials = [session_completed_trials; padded_completed_trials];
            end
		end

            transitionProbs = zeros(1, numStates);
			transitionProbsGood = zeros(1, numStates);
			transitionProbsBad = zeros(1, numStates);
            totalTrials = 0;
			totalTrialsGood = 0;
			totalTrialsBad = 0;
            
            for block = blocks'

                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);

				 if mod(block, 2) == 1  % Odd Block (Non-Swap)
    				unique_trials = unique(block_data.TrialNumberInBlock);
    				numTrials = numel(unique_trials);
    				rewardedTrials = 0;
    				
    				for trial = 1:numTrials
        				trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
        				maxState = max(trial_data.CurrentState);  % Max state in trial
        				if maxState == 5
            				rewardedTrials = rewardedTrials + 1;
        				end
    				end
    				
    				accuracy = rewardedTrials / numTrials;
    				if accuracy >= 0.8
        				blockClassification = 1;
    				else
        				blockClassification = 0;
					end
				 end

                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);
					
                    % if trial_data.isSwapped(1) == 1 && trial <= learning_speed
					if trial_data.isSwapped(1) == 1 && trial == 2
                        for touch = 1:height(trial_data) - 1
                            % A is correctly chosen
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
									if blockClassification
										transitionProbsGood(nextErrorObjectPosition) = transitionProbsGood(nextErrorObjectPosition) + 1;
										totalTrialsGood = totalTrialsGood + 1;
									else
										transitionProbsBad(nextErrorObjectPosition) = transitionProbsBad(nextErrorObjectPosition) + 1;
										totalTrialsBad = totalTrialsBad + 1;
									end
                                end
                                % only record the first correct
                                break;
                            end
                        end
                    end
                end
            end
            
            % Store transition probabilities for this session
            if totalTrials > 0
                subjectProbs = [subjectProbs; transitionProbs / totalTrials];
			end
			if totalTrialsGood > 0
                subjectProbsGood = [subjectProbsGood; transitionProbsGood / totalTrialsGood];
			end
			if totalTrials > 0
                subjectProbsBad = [subjectProbsBad; transitionProbsBad / totalTrialsBad];
			end

       
    end
    allSessionProbs_Swapped{s} = subjectProbs;
	allSessionProbs_Swapped_Good{s} = subjectProbsGood;
	allSessionProbs_Swapped_Bad{s} = subjectProbsBad;
end

    % Calculate overall mean and SEM
    allProbsCombined = vertcat(allSessionProbs_Swapped{:});
    meanProbs = mean(allProbsCombined, 1);
    semProbs = std(allProbsCombined, 0, 1) / sqrt(size(allProbsCombined, 1));
	

	allProbsCombinedGood = vertcat(allSessionProbs_Swapped_Good{:});
    meanProbsGood = mean(allProbsCombinedGood, 1);
    semProbsGood = std(allProbsCombinedGood, 0, 1) / sqrt(size(allProbsCombinedGood, 1));


	allProbsCombinedBad = vertcat(allSessionProbs_Swapped_Bad{:});
    meanProbsBad = nanmean(allProbsCombinedBad, 1);
    semProbsBad = nanstd(allProbsCombinedBad, 0, 1) / sqrt(size(~isnan(allProbsCombinedBad), 1));

    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot swarm chart for each transition state
    for i = 1:numStates
        for s = 1:numSubjects
            subjectData = allSessionProbs_Swapped{s}(:, i);
            sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
                'Marker', markers{s}, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            sm.XJitterWidth = 0.5;
        end
    end

    % Plot error bars
    errorbar(1:numStates, meanProbs, semProbs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
	
	% Plot error bars
    % errorbar(1:numStates, meanProbsGood, semProbsGood, 'r', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
	% 
	% % Plot error bars
    % errorbar(1:numStates, meanProbsBad, semProbsBad, 'b', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    % Perform Welch's t-test between B and D 
    [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 4), 'Vartype', 'unequal');

    % Add significance star if p < 0.05 and print p-value
    if p_value < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
		if p_value < 0.001
        text(3, maxY*1.35, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value < 0.01
			text(3, maxY*1.35, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value < 0.05
			text(3, maxY*1.35, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
	else
		maxY = max(allProbsCombined,[], 'all');
        plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
		text(3, maxY*1.4, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	end
    fprintf('P-value for Ordinal Position 2 vs 4 comparison: %f\n', p_value);

	    % Perform Welch's t-test between B and C 
    [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 3), 'Vartype', 'unequal');

    % Add significance star if p < 0.05 and print p-value
    if p_value < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
		if p_value < 0.001
        text(2.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value < 0.01
			text(2.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value < 0.05
			text(2.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
    else
		maxY = max(allProbsCombined,[], 'all');
        plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
		text(2.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	end
    fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);

	    % Perform Welch's t-test between C and D
    [~, p_value] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 4), 'Vartype', 'unequal');

    % Add significance star if p < 0.05 and print p-value
    if p_value < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
		if p_value < 0.001
        text(3.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value < 0.01
			text(3.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value < 0.05
			text(3.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
    else
		maxY = max(allProbsCombined,[], 'all');
        plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
		text(3.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	end
    fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);


    % Customize plot
    ylabel('Transition Probability', 'FontSize', 12);
    xlabel('Ordinal Position', 'FontSize', 12);
    title('Transition Probability After Correctly Choose A: Swapped, Before LP', 'FontSize', 14);
    xticks(1:numStates);
    xticklabels(transitionStates);
    xlim([0.5, 6.5]);
	ylim([-0.1, 1.2]);
	yline(0,'k-', HandleVisibility="off");
	yticks(0:0.2:1);
	set(gcf,'renderer','Painters')

    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
	set(gcf,'renderer','Painters')
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

function plotLastFourTrialsAccuracy(accuracy_matrices, touch_table)
    subjects = unique(touch_table.Subject);
    colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A'); [107, 37, 110]/255};
    markers = {'o', 's', 'd', '^'};
    chance_levels1 = [0.28571,0.33333, 0.4, 0.5, 0.66667];
    chance_levels2 = zeros(1,7) +0.28571;
    % First Figure - Original Plot
    figure('Position', [100 100 300 250]);
    
    % Get number of subjects
    num_subjects = size(accuracy_matrices, 3);
    all_avg_accuracy = zeros(num_subjects, 5);
    
    % For each subject
    for subj = 1:num_subjects
        % Get last 5 trials (trials 12-15)
        last_four_trials = accuracy_matrices(:, 11:15, subj);
        
        % Average across these trials
        avg_accuracy = mean(last_four_trials, 2);
        all_avg_accuracy(subj, :) = avg_accuracy';
        
        % Plot line for this subject with transparency
        plot(1:5, avg_accuracy, '-o', 'Color', [colors{subj}, 0.8], ...
             'LineWidth', 1, 'Marker', markers{subj},'MarkerSize', 6,'MarkerFaceColor', [colors{subj}],...
             'DisplayName', ['Subject ' subjects{subj}(1)]);
        hold on;
    end
    
    % Plot chance level
	plot(1:5, chance_levels1, '--','Color','[0.1,0.1,0.1]', 'DisplayName', 'Chance Performance w/o Choosing Previous Objects');
    plot(0:6, chance_levels2, 'k--', 'DisplayName', 'Chance Performance');
    
    % Calculate mean and standard error across subjects
    mean_accuracy = mean(all_avg_accuracy, 1);
    std_error = std(all_avg_accuracy, [], 1) / sqrt(num_subjects); % Standard error
    
    % Plot average of all subjects with error bars
    errorbar(1:5, mean_accuracy, std_error, '-o', 'Color', colors{5}, ...
             'LineWidth', 2, 'MarkerSize', 6, ...
             'MarkerFaceColor', colors{5}, ...
             'DisplayName', 'Average');
	disp(mean_accuracy);
	disp(std_error);

    plot(1:5, mean_accuracy, '-', 'Color', colors{5}, ...
         'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', 'Average');
    
    % Customize plot
    xlabel('Ordinal Position');
    ylabel('Accuracy');
    title('Average Accuracy for Trials (12-15)');
    xlim([0.5 5.5]);
    ylim([0 1]);
    xticks(1:5);
    legend_obj = legend('Location', 'eastoutside');
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    
    % Adjust figure sizes
    canvas_size = [300, 250];
    for fig_num = 1:1
        figure(fig_num)
        if fig_num == 1 && exist('legend_obj', 'var')
            legend_pos = get(legend_obj, 'Position');
            legend_width = legend_pos(3);
            new_fig_width = canvas_size(1) + legend_width * canvas_size(1);
            set(gcf, 'Position', [350, 250, new_fig_width, canvas_size(2)]);
		end
        set(gcf, 'renderer', 'Painters');
    end
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




function plotInitialPerformaceVsOrdinalObjectIdentityTwoConditions(block_table)
    % Sort the block table
    block_table = sortrows(block_table, {'Subject', 'SessionNumber', 'BlockNumber'});

    % Initialize arrays to store data
    initial_completion_rates = [];
    error_inferences = [];
    subjects = {};

    % Define colors and markers for each subject
    colors = {hex2rgb('#DB432C'), hex2rgb('#6C96CC'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

    % Process the data
    for i = 2:height(block_table)
        if block_table.isSwapped(i) == 1 && block_table.BlockNumber(i) > 1
            prev_block = i - 1;
            if strcmp(block_table.Subject{i}, block_table.Subject{prev_block}) && ...
                block_table.SessionNumber(i) == block_table.SessionNumber(prev_block) && ...
                block_table.isSwapped(prev_block) == 0

                initial_rate = block_table.BlockLevelAccuracy(prev_block);
                error_inference = block_table.SwappedErrorCorrection(i);

                % Check for valid data
                if isfinite(initial_rate) && isfinite(error_inference)
                    initial_completion_rates = [initial_completion_rates; initial_rate];
                    error_inferences = [error_inferences; error_inference];
                    subjects{end+1} = block_table.Subject{i};
                end
            end
        end
    end

    % Separate data into two conditions
    low_completion = initial_completion_rates < 0.7;
    high_completion = initial_completion_rates >= 0.7;

    % Calculate statistics for each condition
    low_mean = mean(error_inferences(low_completion));
    high_mean = mean(error_inferences(high_completion));
    low_sem = std(error_inferences(low_completion)) / sqrt(sum(low_completion));
    high_sem = std(error_inferences(high_completion)) / sqrt(sum(high_completion));
    low_ci = low_mean + [-1.96 * low_sem, 1.96 * low_sem];
    high_ci = high_mean + [-1.96 * high_sem, 1.96 * high_sem];

    % Print statistics
    fprintf('Low completion rate (< 80%%):\n');
    fprintf('Mean: %.4f, 95%% CI: [%.4f, %.4f]\n', low_mean, low_ci(1), low_ci(2));
    fprintf('High completion rate (>= 80%%):\n');
    fprintf('Mean: %.4f, 95%% CI: [%.4f, %.4f]\n', high_mean, high_ci(1), high_ci(2));

    % Plotting
    canvas_size = [100, 100, 300, 250];
    fig = figure('Position', canvas_size);

    hold on;

    % Create swarm chart
    	
    % Add subject-specific markers
    unique_subjects = unique(subjects);
    for i = 1:length(unique_subjects)
        subject_mask = cellfun(@(x) strcmp(x, unique_subjects{i}), subjects)';
        s1 = swarmchart(ones(sum(subject_mask&low_completion), 1), error_inferences(subject_mask&low_completion), 30, 'filled', 'MarkerFaceColor', colors{1}, 'MarkerFaceAlpha', 0.5, 'Marker', markers{i});
    	s2 = swarmchart(2 * ones(sum(subject_mask&high_completion), 1), error_inferences(subject_mask&high_completion), 30, 'filled', 'MarkerFaceColor', colors{2}, 'MarkerFaceAlpha', 0.5, 'Marker', markers{i});
		s1.XJitterWidth = 0.6;
		s2.XJitterWidth = 0.6;
    end

    % Add error bars
    errorbar([1, 2], [low_mean, high_mean], ...
        [low_ci(2)-low_mean, high_ci(2)-high_mean], ...
        [low_mean-low_ci(1), high_mean-high_ci(1)], ...
        'k', 'LineWidth', 1, 'CapSize', 10);

    % Add error bars
    errorbar([1, 2], [low_mean, high_mean], ...
        [low_ci(2)-low_mean, high_ci(2)-high_mean], ...
        [low_mean-low_ci(1), high_mean-high_ci(1)], ...
        'k', 'LineWidth', 1, 'CapSize', 10);

    % Customize the plot
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Low', 'High'}, 'TickDir', 'out', 'FontSize', 12);
    xlabel('Initial Completion Rate', 'FontSize', 12);
    ylabel('Ordinal Object identity inference', 'FontSize', 12);
%     title('Anticipated Swapping vs Initial Completion Rate', 'FontSize', 14);
	set(gcf,'renderer','Painters');

	% Set x and y limits
    xlim([0.5 2.5]);
    y_max = 1.25 * max(error_inferences);
    ylim([min(error_inferences) y_max]);

    % Perform Welch's t-test
    [~, p] = ttest2(error_inferences(low_completion), error_inferences(high_completion), 'Vartype', 'unequal');

    % Determine significance level and star representation
    
    if p < 0.01
        star = '***';
    elseif p < 0.05
        star = '**';
    elseif p < 0.1
        star = '*';
    else
        star = 'ns';
    end
	% Add significance bar and star
    y_bar = y_max * 0.9;
    plot([1, 2], [y_bar, y_bar], 'k-', 'LineWidth', 1);
    text(1.5, y_bar * 1.025, star, 'HorizontalAlignment', 'center', 'FontSize', 12);
	fprintf('p-value: %.10f\n',p)


   % Add legend for subjects
	numSubjects = 4;
	subjects = unique(block_table.Subject);
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    new_fig_size = [canvas_size(1), canvas_size(2), ...
                    canvas_size(3) + legend_pos(3), canvas_size(4)];
    set(fig, 'Position', new_fig_size);
end
function plotInitialPerformaceVsAnticipatedSwappingTwoConditions(block_table)
    % Sort the block table
    block_table = sortrows(block_table, {'Subject', 'SessionNumber', 'BlockNumber'});

    % Initialize arrays to store data
    initial_completion_rates = [];
    error_inferences = [];
    subjects = {};

    % Define colors and markers for each subject
    colors = {hex2rgb('#DB432C'), hex2rgb('#6C96CC'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

    % Process the data
    for i = 2:height(block_table)
        if block_table.isSwapped(i) == 1 && block_table.BlockNumber(i) > 1
            prev_block = i - 1;
            if strcmp(block_table.Subject{i}, block_table.Subject{prev_block}) && ...
                block_table.SessionNumber(i) == block_table.SessionNumber(prev_block) && ...
                block_table.isSwapped(prev_block) == 0

                initial_rate = block_table.BlockLevelAccuracy(prev_block);
                error_inference = block_table.AnticipatedSwapping(i);

                % Check for valid data
                if isfinite(initial_rate) && isfinite(error_inference)
                    initial_completion_rates = [initial_completion_rates; initial_rate];
                    error_inferences = [error_inferences; error_inference];
                    subjects{end+1} = block_table.Subject{i};
                end
            end
        end
    end

    % Separate data into two conditions
    low_completion = initial_completion_rates < 0.7;
    high_completion = initial_completion_rates >= 0.7;

    % Calculate statistics for each condition
    low_mean = mean(error_inferences(low_completion));
    high_mean = mean(error_inferences(high_completion));
    low_sem = std(error_inferences(low_completion)) / sqrt(sum(low_completion));
    high_sem = std(error_inferences(high_completion)) / sqrt(sum(high_completion));
    low_ci = low_mean + [-1.96 * low_sem, 1.96 * low_sem];
    high_ci = high_mean + [-1.96 * high_sem, 1.96 * high_sem];

    % Print statistics
    fprintf('Low completion rate (< 80%%):\n');
    fprintf('Mean: %.4f, 95%% CI: [%.4f, %.4f], CI: %.4f\n', low_mean, low_ci(1), low_ci(2), 1.96 * low_sem);
    fprintf('High completion rate (>= 80%%):\n');
    fprintf('Mean: %.4f, 95%% CI: [%.4f, %.4f], CI: %.4f\n', high_mean, high_ci(1), high_ci(2), 1.96 * high_sem);

    % Plotting
    canvas_size = [100, 100, 300, 250];
    fig = figure('Position', canvas_size);

    hold on;

    % Create swarm chart
    	
    % Add subject-specific markers
    unique_subjects = unique(subjects);
    for i = 1:length(unique_subjects)
        subject_mask = cellfun(@(x) strcmp(x, unique_subjects{i}), subjects)';
        s1 = swarmchart(ones(sum(subject_mask&low_completion), 1), error_inferences(subject_mask&low_completion), 50, 'filled', 'MarkerFaceColor', colors{1}, 'MarkerFaceAlpha', 0.5, 'Marker', markers{i});
    	s2 = swarmchart(2 * ones(sum(subject_mask&high_completion), 1), error_inferences(subject_mask&high_completion), 50, 'filled', 'MarkerFaceColor', colors{2}, 'MarkerFaceAlpha', 0.5, 'Marker', markers{i});
		s1.XJitterWidth = 0.6;
		s2.XJitterWidth = 0.6;
    end

    % Add error bars
    errorbar([1, 2], [low_mean, high_mean], ...
        [low_ci(2)-low_mean, high_ci(2)-high_mean], ...
        [low_mean-low_ci(1), high_mean-high_ci(1)], ...
        'k', 'LineWidth', 1, 'CapSize', 10);

    % Add error bars
    errorbar([1, 2], [low_mean, high_mean], ...
        [low_ci(2)-low_mean, high_ci(2)-high_mean], ...
        [low_mean-low_ci(1), high_mean-high_ci(1)], ...
        'k', 'LineWidth', 1, 'CapSize', 10);

    % Customize the plot
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Low', 'High'}, 'TickDir', 'out', 'FontSize', 12);
    xlabel('Initial Completion Rate', 'FontSize', 12);
    ylabel('Anticipated Swapping', 'FontSize', 12);
    title('Anticipated Swapping vs Initial Completion Rate', 'FontSize', 14);
	set(gcf,'renderer','Painters');

	% Set x and y limits
    xlim([0.5 2.5]);
    y_max = 1.25 * max(error_inferences);
    ylim([min(error_inferences) y_max]);

    % Perform Welch's t-test
    [~, p] = ttest2(error_inferences(low_completion), error_inferences(high_completion), 'Vartype', 'unequal');

    % Determine significance level and star representation
    
    if p < 0.01
        star = '***';
    elseif p < 0.05
        star = '**';
    elseif p < 0.1
        star = '*';
    else
        star = 'ns';
    end
	% Add significance bar and star
    y_bar = y_max * 0.9;
    plot([1, 2], [y_bar, y_bar], 'k-', 'LineWidth', 1);
    text(1.5, y_bar * 1.025, star, 'HorizontalAlignment', 'center', 'FontSize', 12);
	fprintf('p-value: %.10f\n',p)

   % Add legend for subjects
	numSubjects = 4;
	subjects = unique(block_table.Subject);
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    new_fig_size = [canvas_size(1), canvas_size(2), ...
                    canvas_size(3) + legend_pos(3), canvas_size(4)];
    set(fig, 'Position', new_fig_size);
end

function [updated_block_table] = calculateAnticipatedSwapping(touch_table, block_table)
    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Initialize cell array to store AnticipatedSwapping values
    anticipatedSwapping = cell(0, 4);

    % Calculate transition probabilities for each subject and session
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            blocks = unique(session_data.BlockNumber);
            
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                transitionProbs = zeros(1, numStates);
                totalTrials = 0;
                trials = unique(block_data.TrialNumberInBlock);
                
                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);
                    if trial_data.isSwapped(1) == 1
                        for touch = 1:height(trial_data) - 1
                            % A is correctly chosen
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
                                end
                                % only record the first correct
                                break;
                            end
                        end
                    end
                end
                
                % Calculate AnticipatedSwapping value
                if totalTrials > 0
                    normalizedProbs = transitionProbs / totalTrials;
                    anticipatedSwappingValue = normalizedProbs(4) - normalizedProbs(2);
                else
                    anticipatedSwappingValue = NaN;
                end
                
                % Store AnticipatedSwapping value
                anticipatedSwapping(end+1, :) = {subject, session, block, anticipatedSwappingValue};
            end
        end
    end

    % Convert anticipatedSwapping to a table
    anticipatedSwappingTable = cell2table(anticipatedSwapping, 'VariableNames', {'Subject', 'SessionNumber', 'BlockNumber', 'AnticipatedSwapping'});

    % Ensure that the data types match between anticipatedSwappingTable and block_table
    if iscell(block_table.Subject)
        anticipatedSwappingTable.Subject = anticipatedSwappingTable.Subject;
    elseif iscategorical(block_table.Subject)
        anticipatedSwappingTable.Subject = categorical(anticipatedSwappingTable.Subject);
    else
        anticipatedSwappingTable.Subject = string(anticipatedSwappingTable.Subject);
        block_table.Subject = string(block_table.Subject);
    end
    
    anticipatedSwappingTable.SessionNumber = double(anticipatedSwappingTable.SessionNumber);
    anticipatedSwappingTable.BlockNumber = double(anticipatedSwappingTable.BlockNumber);

    % Add AnticipatedSwapping to block_table
    updated_block_table = outerjoin(block_table, anticipatedSwappingTable, 'Keys', {'Subject', 'SessionNumber', 'BlockNumber'}, 'MergeKeys', true);

    % Handle any potential missing values
    updated_block_table.AnticipatedSwapping = fillmissing(updated_block_table.AnticipatedSwapping, 'constant', NaN);
end


function [updated_block_table] = calculateOrdinalObjectIdentityInference(touch_table, block_table)
    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Initialize cell array to store AnticipatedSwapping values
    anticipatedSwapping = cell(0, 4);

    % Calculate transition probabilities for each subject and session
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            blocks = unique(session_data.BlockNumber);
            
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                transitionProbs = zeros(1, numStates);
                totalTrials = 0;
                trials = unique(block_data.TrialNumberInBlock);
                
                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);
                        if trial_data.isSwapped(1) == 1
                        for touch = 1:height(trial_data) - 1
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection')) ...
                                    && trial_data.CurrentState(touch+1) == 0 && trial_data.TouchObjectCorrectOrdinalPosition(touch+1) == 4

                                nextTouchIndex = touch + 2;
                                if nextTouchIndex + 1 <= height(trial_data) && strcmpi(trial_data.TouchCategory{nextTouchIndex}, 'RetouchCorrect')
                                    nextTouchIndex = nextTouchIndex + 1;
                                end
                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
                                end
                                % only record the first correct
                                break;
                            end
                        end
                    end
                end
                
                % Calculate AnticipatedSwapping value
                if totalTrials > 0
                    normalizedProbs = transitionProbs / totalTrials;
                    anticipatedSwappingValue = normalizedProbs(4) - normalizedProbs(3);
                else
                    anticipatedSwappingValue = NaN;
                end
                
                % Store AnticipatedSwapping value
                anticipatedSwapping(end+1, :) = {subject, session, block, anticipatedSwappingValue};
            end
        end
    end
 % Convert anticipatedSwapping to a table
    anticipatedSwappingTable = cell2table(anticipatedSwapping, 'VariableNames', {'Subject', 'SessionNumber', 'BlockNumber', 'ObjectIdentityInference'});

    % Ensure that the data types match between anticipatedSwappingTable and block_table
    if iscell(block_table.Subject)
        anticipatedSwappingTable.Subject = anticipatedSwappingTable.Subject;
    elseif iscategorical(block_table.Subject)
        anticipatedSwappingTable.Subject = categorical(anticipatedSwappingTable.Subject);
    else
        anticipatedSwappingTable.Subject = string(anticipatedSwappingTable.Subject);
        block_table.Subject = string(block_table.Subject);
    end
    
    anticipatedSwappingTable.SessionNumber = double(anticipatedSwappingTable.SessionNumber);
    anticipatedSwappingTable.BlockNumber = double(anticipatedSwappingTable.BlockNumber);

    % Add ObjectIdentityInference to block_table
    updated_block_table = outerjoin(block_table, anticipatedSwappingTable, 'Keys', {'Subject', 'SessionNumber', 'BlockNumber'}, 'MergeKeys', true);

    % Check if ObjectIdentityInference column exists, if not, add it
    if ~ismember('ObjectIdentityInference', updated_block_table.Properties.VariableNames)
        updated_block_table.ObjectIdentityInference = NaN(height(updated_block_table), 1);
    end

    % Handle any potential missing values
    updated_block_table.ObjectIdentityInference = fillmissing(updated_block_table.ObjectIdentityInference, 'constant', NaN);
end

function plotSessionWiseSecondToThirdTransitionProbability(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Define markers and colors for subjects
    markers = {'o', 's', 'd', '^', 'v'};
    colors = [39, 93, 44;
			58, 112, 175;
			107, 37, 110;
			199, 54, 55;
			233, 166, 64;
			247, 211, 76]/255;
	accuracy_threshold = 0.7;
    % Initialize array to store session-wise transition probabilities
    allSessionProbs = cell(numSubjects, 1);

	% Calculate transition probabilities for each subject and session
for s = 1:numSubjects
    subject = subjects{s};
    subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
    sessions = unique(subject_data.SessionNumber);
    subjectProbs = [];

    for session = sessions'
        session_data = subject_data(subject_data.SessionNumber == session, :);
        blocks = unique(session_data.BlockNumber);
        session_completed_trials = [];
        
        % Calculate learning speed for the session
        for blk = 1:numel(blocks)
            block_data = session_data(session_data.BlockNumber == blocks(blk), :);
            % Ensure the block is not swapped and not a repetition
            if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                unique_trials = unique(block_data.TrialNumberInBlock);
                numTrials = numel(unique_trials);
                completed_trials = zeros(1, numTrials);
                for trial = 1:numTrials
                    trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                    completed_trials(trial) = any(trial_data.isRewarded == 1);
                end
                % Pad completed trials to 15 if necessary
                padded_completed_trials = nan(1, 15);
                padded_completed_trials(1:length(completed_trials)) = completed_trials;
                session_completed_trials = [session_completed_trials; padded_completed_trials];
            end
        end

        % Calculate learning speed for the session
        learning_speed = NaN;
        if ~isempty(session_completed_trials)
            mean_completion_rate = nanmean(session_completed_trials, 1);
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
        end

        % Calculate choice probabilities only if learning_speed is NaN (i.e., learning point not reached)
        if isnan(learning_speed)
			continue;
			learning_speed = 15;
		end
            transitionProbs = zeros(1, numStates);
            totalTrials = 0;

            blocks = unique(session_data.BlockNumber);
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);

                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);

                    if trial_data.isSwapped(1) == 1 && trial <= learning_speed
                        for touch = 1:height(trial_data) - 1
                            if trial_data.CurrentState(touch) == 2 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))

                                nextTouchIndex = touch + 1;

                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
								end
								break;
                            end
                        end
                    end
                end
            end

            % Store transition probabilities for this session
            if totalTrials > 0
                subjectProbs = [subjectProbs; transitionProbs / totalTrials];
            end
        end
        allSessionProbs{s} = subjectProbs;
    end

    % Calculate overall mean and SEM
    allProbsCombined = vertcat(allSessionProbs{:});
    meanProbs = mean(allProbsCombined, 1);
    semProbs = std(allProbsCombined, 0, 1) / sqrt(size(allProbsCombined, 1));
	disp(meanProbs);
	disp(semProbs);
	

    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot swarm chart for each transition state
    for i = 1:numStates
        for s = 1:numSubjects
            subjectData = allSessionProbs{s}(:, i);
            sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
                'Marker', markers{s}, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            sm.XJitterWidth = 0.5;
        end
    end

    % Plot error bars
    errorbar(1:numStates, meanProbs, semProbs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    fprintf('Transition Probabilities (mean ± 95%% CI):\n');
    for i = 1:numStates
        fprintf('%s: %.3f ± %.3f\n', transitionStates{i}, meanProbs(i), semProbs(i));
    end

    % Perform Welch's t-test between C and E (now 3 and 2)
    [~, p_value] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 5), 'Vartype', 'unequal');

    % Add significance star if p < 0.05 and print p-value
    if p_value < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([3, 5], [maxY*1.1, maxY*1.1], 'k-');
		if p_value < 0.001
        text(4, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value < 0.01
			text(4, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value < 0.05
			text(4, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
    end
    fprintf('P-value for Ordinal Position 3 vs 4 comparison: %.10f\n', p_value);

    % Customize plot
    ylabel('Transition Probability', 'FontSize', 12);
    xlabel('Ordinal Position', 'FontSize', 12);
    title('Transition Probability After Correctly Choosing D', 'FontSize', 14);
    xticks(1:numStates);
    xticklabels(transitionStates);
    xlim([0.5, 6.5]);
	ylim([-0.1, 1.2]);
	yline(0,'k-', HandleVisibility="off");
	yticks(0:0.2:1);
	set(gcf,'renderer','Painters')

    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end



function plotSessionWiseAnticipatedSwappingTransitionProbability(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Define markers and colors for subjects
    markers = {'o', 's', 'd', '^', 'v'};
    colors = [39, 93, 44;
			58, 112, 175;
			107, 37, 110;
			199, 54, 55;
			233, 166, 64;
			247, 211, 76]/255;
	
    % Initialize array to store session-wise transition probabilities
allSessionProbs_Swapped = cell(numSubjects, 1);
accuracy_threshold = 0.8; % Set your desired accuracy threshold

% Calculate transition probabilities for each subject and session
for s = 1:numSubjects
    subject = subjects{s};
    subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
    sessions = unique(subject_data.SessionNumber);
    subjectProbs = [];

    for session = sessions'
        session_data = subject_data(subject_data.SessionNumber == session, :);
        blocks = unique(session_data.BlockNumber);
        session_completed_trials = [];
        
        % Calculate learning speed for the session
        for blk = 1:numel(blocks)
            block_data = session_data(session_data.BlockNumber == blocks(blk), :);
            % Ensure the block is not swapped and not a repetition
            if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                unique_trials = unique(block_data.TrialNumberInBlock);
                numTrials = numel(unique_trials);
                completed_trials = zeros(1, numTrials);
                for trial = 1:numTrials
                    trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                    completed_trials(trial) = any(trial_data.isRewarded == 1);
                end
                % Pad completed trials to 15 if necessary
                padded_completed_trials = nan(1, 15);
                padded_completed_trials(1:length(completed_trials)) = completed_trials;
                session_completed_trials = [session_completed_trials; padded_completed_trials];
            end
        end

        % Calculate learning speed for the session
        learning_speed = NaN;
        if ~isempty(session_completed_trials)
            mean_completion_rate = nanmean(session_completed_trials, 1);
            for i = 1:length(mean_completion_rate)
				if i < length(mean_completion_rate) - 3
                	if mean(mean_completion_rate(i:i+3)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                    	learning_speed = i;
                    	break;
					end
				else
					if mean_completion_rate(i) >= accuracy_threshold
                    	learning_speed = i;
                    	break;
					end
				end
            end
        end

        % Calculate choice probabilities only if learning_speed is NaN (i.e., learning point not reached)
        if isnan(learning_speed)
			continue;
			learning_speed = 15;
		end
            transitionProbs = zeros(1, numStates);
            totalTrials = 0;
            
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);
                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);

                    if trial_data.isSwapped(1) == 1 && trial <= learning_speed
                        for touch = 1:height(trial_data) - 1
                            % A is correctly chosen
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
                                end
                                % only record the first correct
                                break;
                            end
                        end
                    end
                end
            end
            
            % Store transition probabilities for this session
            if totalTrials > 0
                subjectProbs = [subjectProbs; transitionProbs / totalTrials];
            end
       
    end
    allSessionProbs_Swapped{s} = subjectProbs;
end

    % Calculate overall mean and SEM
    allProbsCombined = vertcat(allSessionProbs_Swapped{:});
    meanProbs = mean(allProbsCombined, 1);
    semProbs = std(allProbsCombined, 0, 1) / sqrt(size(allProbsCombined, 1));
	disp('Figure 1');
	disp(meanProbs);
	disp(semProbs);
    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot swarm chart for each transition state
    for i = 1:numStates
        for s = 1:numSubjects
            subjectData = allSessionProbs_Swapped{s}(:, i);
            sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
                'Marker', markers{s}, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            sm.XJitterWidth = 0.5;
        end
    end

    % Plot error bars
    errorbar(1:numStates, meanProbs, semProbs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    % % Perform Welch's t-test between B and D 
    % [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 4), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
	% 	if p_value < 0.001
    %     text(3, maxY*1.35, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(3, maxY*1.35, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(3, maxY*1.35, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
	% else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
	% 	text(3, maxY*1.4, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 4 comparison: %f\n', p_value);
	% 
	%     % Perform Welch's t-test between B and C 
    % [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 3), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	if p_value < 0.001
    %     text(2.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(2.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(2.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
    % else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	text(2.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);
	% 
	%     % Perform Welch's t-test between C and D
    % [~, p_value] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 4), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	if p_value < 0.001
    %     text(3.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(3.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(3.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
    % else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	text(3.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);
	
		% Define which positions to compare (B=2, C=3, D=4)
    positions = [2, 3, 4];
    comparisons = nchoosek(positions, 2); % [2,3], [2,4], [3,4]
    nComp = size(comparisons,1);

    % Collect p-values for all comparisons
    pvals = zeros(nComp,1);
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        [~, pvals(i)] = ttest2(allProbsCombined(:, pos1), allProbsCombined(:, pos2), 'Vartype', 'unequal');
	end


    % --- Bonferroni correction ---
    p_bonf = min(pvals * nComp, 1); % [3]
    % --- Benjamini-Hochberg FDR correction ---
    [~, sortIdx] = sort(pvals);
    p_fdr = zeros(size(pvals));
    for k = 1:nComp
        idx = sortIdx(k);
        p_fdr(idx) = min(pvals(idx) * nComp / k, 1);
    end

    % Choose which correction to use for annotation
    % For Bonferroni: use p_bonf
    % For FDR: use p_fdr
    p_to_use = p_bonf; % or p_fdr

    % Plotting significance bars and stars
    maxY = max(allProbsCombined, [], 'all');
    y_offset = maxY * 0.2;
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        y_pos = maxY + i * y_offset;
        plot([pos1, pos2], [y_pos, y_pos], 'k-', 'LineWidth', 1);
        % Annotate with stars based on corrected p-value
        star = getStars(p_to_use(i));
        text(mean([pos1, pos2]), y_pos + y_offset*0.3, star, 'HorizontalAlignment', 'center', 'FontSize', 12);
        fprintf('Corrected P-value for Ordinal Position %d vs %d: %f\n', pos1, pos2, p_to_use(i));
	end


    % Customize plot
    ylabel('Transition Probability', 'FontSize', 12);
    xlabel('Ordinal Position', 'FontSize', 12);
    title('Transition Probability After Correctly Choose A: Swapped, Before LP', 'FontSize', 14);
    xticks(1:numStates);
    xticklabels(transitionStates);
    xlim([0.5, 6.5]);
	ylim([-0.1, 1.2]);
	yline(0,'k-', HandleVisibility="off");
	yticks(0:0.2:1);
	set(gcf,'renderer','Painters')

    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
	set(gcf,'renderer','Painters')
    hold off;



	    % Initialize array to store session-wise transition probabilities
    allSessionProbs_Initial = cell(numSubjects, 1);

% Calculate transition probabilities for each subject and session
for s = 1:numSubjects
    subject = subjects{s};
    subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
    sessions = unique(subject_data.SessionNumber);
    subjectProbs = [];

    for session = sessions'
        session_data = subject_data(subject_data.SessionNumber == session, :);
        blocks = unique(session_data.BlockNumber);
        session_completed_trials = [];
        
        % Calculate learning speed for the session
        for blk = 1:numel(blocks)
            block_data = session_data(session_data.BlockNumber == blocks(blk), :);
            % Ensure the block is not swapped and not a repetition
            if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                unique_trials = unique(block_data.TrialNumberInBlock);
                numTrials = numel(unique_trials);
                completed_trials = zeros(1, numTrials);
                for trial = 1:numTrials
                    trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                    completed_trials(trial) = any(trial_data.isRewarded == 1);
                end
                % Pad completed trials to 15 if necessary
                padded_completed_trials = nan(1, 15);
                padded_completed_trials(1:length(completed_trials)) = completed_trials;
                session_completed_trials = [session_completed_trials; padded_completed_trials];
            end
        end

        % Calculate learning speed for the session
        learning_speed = NaN;
        if ~isempty(session_completed_trials)
            mean_completion_rate = nanmean(session_completed_trials, 1);
            for i = 1:length(mean_completion_rate)
				if i < length(mean_completion_rate) - 3
                	if mean(mean_completion_rate(i:i+3)) >= accuracy_threshold && mean_completion_rate(i) >= accuracy_threshold
                    	learning_speed = i;
                    	break;
					end
				else
					if mean_completion_rate(i) >= accuracy_threshold
                    	learning_speed = i;
                    	break;
					end
				end
            end
        end

        % Calculate choice probabilities only if learning_speed is NaN (i.e., learning point not reached)
        if isnan(learning_speed)
			continue;
			learning_speed = 15;
		end
            transitionProbs = zeros(1, numStates);
            totalTrials = 0;
            
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);
                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);

                    if trial_data.isSwapped(1) == 0 && trial <= learning_speed
                        for touch = 1:height(trial_data) - 1
                            % A is correctly chosen
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection'))
                                nextTouchIndex = touch + 1;
                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 2;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 4;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
                                end
                                % only record the first correct
                                break;
                            end
                        end
                    end
                end
            end
            
            % Store transition probabilities for this session
            if totalTrials > 0
                subjectProbs = [subjectProbs; transitionProbs / totalTrials];
            end
       
    end
    allSessionProbs_Initial{s} = subjectProbs;
end
    % Calculate overall mean and SEM
    allProbsCombined = vertcat(allSessionProbs_Initial{:});
    meanProbs = mean(allProbsCombined, 1);
    semProbs = std(allProbsCombined, 0, 1) / sqrt(size(allProbsCombined, 1));
	disp('Figure 2');
	disp(meanProbs);
	disp(semProbs);
    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot swarm chart for each transition state
    for i = 1:numStates
        for s = 1:numSubjects
            subjectData = allSessionProbs_Initial{s}(:, i);
            sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
                'Marker', markers{s}, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            sm.XJitterWidth = 0.5;
        end
    end

    % Plot error bars
    errorbar(1:numStates, meanProbs, semProbs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
	% 
    % % Perform Welch's t-test between B and D 
    % [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 4), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
	% 	if p_value < 0.001
    %     text(3, maxY*1.35, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(3, maxY*1.35, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(3, maxY*1.35, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
	% else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([2, 4], [maxY*1.3, maxY*1.3], 'k-');
	% 	text(3, maxY*1.4, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 4 comparison: %f\n', p_value);
	% 
	%     % Perform Welch's t-test between B and C 
    % [~, p_value] = ttest2(allProbsCombined(:, 2), allProbsCombined(:, 3), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	if p_value < 0.001
    %     text(2.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(2.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(2.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
    % else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([2.1, 2.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	text(2.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);
	% 
	%     % Perform Welch's t-test between C and D
    % [~, p_value] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 4), 'Vartype', 'unequal');
	% 
    % % Add significance star if p < 0.05 and print p-value
    % if p_value < 0.05
    %     maxY = max(allProbsCombined,[], 'all');
    %     plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	if p_value < 0.001
    %     text(3.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif  p_value < 0.01
	% 		text(3.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	elseif p_value < 0.05
	% 		text(3.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
	% 	end
    % else
	% 	maxY = max(allProbsCombined,[], 'all');
    %     plot([3.1, 3.9], [maxY*1.1, maxY*1.1], 'k-');
	% 	text(3.5, maxY*1.2, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
	% end
    % fprintf('P-value for Ordinal Position 2 vs 3 comparison: %f\n', p_value);

	% Define which positions to compare (B=2, C=3, D=4)
    positions = [2, 3, 4];
    comparisons = nchoosek(positions, 2); % [2,3], [2,4], [3,4]
    nComp = size(comparisons,1);

    % Collect p-values for all comparisons
    pvals = zeros(nComp,1);
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        [~, pvals(i)] = ttest2(allProbsCombined(:, pos1), allProbsCombined(:, pos2), 'Vartype', 'unequal');
	end


    % --- Bonferroni correction ---
    p_bonf = min(pvals * nComp, 1); % [3]
    % --- Benjamini-Hochberg FDR correction ---
    [~, sortIdx] = sort(pvals);
    p_fdr = zeros(size(pvals));
    for k = 1:nComp
        idx = sortIdx(k);
        p_fdr(idx) = min(pvals(idx) * nComp / k, 1);
    end

    % Choose which correction to use for annotation
    % For Bonferroni: use p_bonf
    % For FDR: use p_fdr
    p_to_use = p_bonf; % or p_fdr

    % Plotting significance bars and stars
    maxY = max(allProbsCombined, [], 'all');
    y_offset = maxY * 0.2;
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        y_pos = maxY + i * y_offset;
        plot([pos1, pos2], [y_pos, y_pos], 'k-', 'LineWidth', 1);
        % Annotate with stars based on corrected p-value
        star = getStars(p_to_use(i));
        text(mean([pos1, pos2]), y_pos + y_offset*0.3, star, 'HorizontalAlignment', 'center', 'FontSize', 12);
        fprintf('Corrected P-value for Ordinal Position %d vs %d: %f\n', pos1, pos2, p_to_use(i));
	end







    % Customize plot
    ylabel('Transition Probability', 'FontSize', 12);
    xlabel('Ordinal Position', 'FontSize', 12);
    title('Transition Probability After Correctly Choose A: Initial, Before LP', 'FontSize', 14);
    xticks(1:numStates);
    xticklabels(transitionStates);
    xlim([0.5, 6.5]);
	ylim([-0.1, 1.2]);
	yline(0,'k-', HandleVisibility="off");
	yticks(0:0.2:1);
	set(gcf,'renderer','Painters')
	
    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;


	% Calculate the difference between Swapped and Initial probabilities
numSubjects = length(allSessionProbs_Swapped);
diffSessionProbs = cell(numSubjects, 1);

for s = 1:numSubjects
    swapped_probs = allSessionProbs_Swapped{s};
    initial_probs = allSessionProbs_Initial{s};
    
    % Ensure the matrices have the same size
    min_rows = min(size(swapped_probs, 1), size(initial_probs, 1));
    
    diff_probs = swapped_probs(1:min_rows, :) - initial_probs(1:min_rows, :);
    diffSessionProbs{s} = diff_probs;
end

% Plotting
canvas_size = [100, 100, 300, 250];
fig = figure('Position', canvas_size);

% Combine all differences into a single matrix
all_diffs = vertcat(diffSessionProbs{:});

% Calculate mean and standard error
mean_diff = mean(all_diffs, 1);
sem_diff = std(all_diffs, 0, 1) / sqrt(size(all_diffs, 1));
disp('Figure 3')
disp(mean_diff);
disp(sem_diff);
% Add error bars
hold on;

% Customize the plot
xlabel('Object', 'FontSize', 12);
ylabel('Difference in Choice Probability', 'FontSize', 12);
title('Difference in Choice Probability (Swapped - Initial)', 'FontSize', 14);
xticks(1:numStates);
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D', 'E','Distractor'}, 'FontSize', 12);
set(gca, 'TickDir', 'out');
xlim([0.5, 6.5]);
set(gcf,'renderer','Painters')
box off;


% Plot swarm chart for each transition state
for i = 1:numStates
    for s = 1:numSubjects
        subjectData = diffSessionProbs{s}(:, i);
        sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
            'Marker', markers{s},'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
        sm.XJitterWidth = 0.5;
    end
end


% Add legend for subjects
legendHandles = gobjects(1, numSubjects);
for s = 1:numSubjects
    legendHandles(s) = plot(NaN, NaN, markers{s}, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
end
legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
    'Location', 'eastoutside', 'FontSize', 12);
legend_obj.EdgeColor = 'none';

e = errorbar(1:length(mean_diff), mean_diff, sem_diff, 'k', 'LineStyle', 'none','LineWidth', 1.5, 'CapSize', 10);
% 
% % Perform pairwise Welch's t-tests
% positions = [2, 3, 4]; % B, C, D positions
% comparisons = nchoosek(positions, 2);
% max_y = max(all_diffs, [], 'all');
% y_offset = max_y * 0.2;
% 
% for i = 1:size(comparisons, 1)
%     pos1 = comparisons(i, 1);
%     pos2 = comparisons(i, 2);
% 
%     [~, p_value] = ttest2(all_diffs(:, pos1), all_diffs(:, pos2), 'Vartype', 'unequal');
% 
%     % Add significance bar and star
%     y_pos = max_y + i * y_offset;
%     plot([pos1, pos2], [y_pos, y_pos], 'k-', 'LineWidth', 1);
% 
%     % Add significance star and p-value
%     if p_value < 0.05
%         if p_value < 0.001
%             star = '***';
%         elseif p_value < 0.01
%             star = '**';
%         else
%             star = '*';
%         end
%         text(mean([pos1, pos2]), y_pos + y_offset*0.3, star, 'HorizontalAlignment', 'center', 'FontSize', 20);
%     else
%         text(mean([pos1, pos2]), y_pos + y_offset*0.3, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
%     end
% 
%     fprintf('P-value for Ordinal Position %d vs %d comparison: %f\n', pos1, pos2, p_value);
% end


	% Define which positions to compare (B=2, C=3, D=4)
    positions = [2, 3, 4];
    comparisons = nchoosek(positions, 2); % [2,3], [2,4], [3,4]
    nComp = size(comparisons,1);

    % Collect p-values for all comparisons
    pvals = zeros(nComp,1);
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        [~, pvals(i)] = ttest2(all_diffs(:, pos1), all_diffs(:, pos2), 'Vartype', 'unequal');
	end


    % --- Bonferroni correction ---
    p_bonf = min(pvals * nComp, 1); % [3]
    % --- Benjamini-Hochberg FDR correction ---
    [~, sortIdx] = sort(pvals);
    p_fdr = zeros(size(pvals));
    for k = 1:nComp
        idx = sortIdx(k);
        p_fdr(idx) = min(pvals(idx) * nComp / k, 1);
    end

    % Choose which correction to use for annotation
    % For Bonferroni: use p_bonf
    % For FDR: use p_fdr
    p_to_use = p_bonf; % or p_fdr

    % Plotting significance bars and stars
    maxY = max(allProbsCombined, [], 'all');
    y_offset = maxY * 0.2;
    for i = 1:nComp
        pos1 = comparisons(i,1);
        pos2 = comparisons(i,2);
        y_pos = maxY + i * y_offset;
        plot([pos1, pos2], [y_pos, y_pos], 'k-', 'LineWidth', 1, HandleVisibility='off');
        % Annotate with stars based on corrected p-value
        star = getStars(p_to_use(i));
        text(mean([pos1, pos2]), y_pos + y_offset*0.3, star, 'HorizontalAlignment', 'center', 'FontSize', 12);
        fprintf('Corrected P-value for Ordinal Position %d vs %d: %f\n', pos1, pos2, p_to_use(i));
	end


% Adjust y-axis limit to accommodate significance bars
% ylim([min(all_diffs, [], 'all') - y_offset * 1.1, max_y*1.1 + (size(comparisons, 1) + 1) * y_offset]);



% Adjust figure size to accommodate legend
legend_pos = get(legend_obj, 'Position');
legend_width = legend_pos(3);
new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
set(fig, 'Position', new_fig_size);

% Add horizontal line at y=0
yline(0, '--k');

hold off;




	function s = getStars(p)
    	if p < 0.001
        	s = '***';
    	elseif p < 0.01
        	s = '**';
    	elseif p < 0.05
        	s = '*';
    	else
        	s = 'ns';
    	end
	end

end


function plotAnticipatedSwappingOverSession(blockTable)
    % Define canvas size
    canvas_size = [100, 100, 400, 250];

    % Extract unique subjects and sessions
    subjects = unique(blockTable.Subject);
    numSubjects = numel(subjects);
    uniqueSessions = unique(blockTable.SessionNumber);

    % Initialize array for average swappedErrorCorrection
    averageSwappedErrorCorrection = nan(length(uniqueSessions), numSubjects);

    % Calculate average swappedErrorCorrection for each session and subject
    for s = 1:numSubjects
        subject = subjects{s};
        for i = 1:length(uniqueSessions)
            sessionNumber = uniqueSessions(i);
            sessionData = blockTable.AnticipatedSwapping(blockTable.SessionNumber == sessionNumber & ...
                blockTable.isSwapped == 1 & strcmp(blockTable.Subject, subject));
            averageSwappedErrorCorrection(i, s) = nanmean(sessionData);
        end
    end

    % Apply new averaging rule
	
	smoothedData = smoothDataOverDays(averageSwappedErrorCorrection);
%     smoothedData = smoothDataBySubjectCount(averageSwappedErrorCorrection);

    % Create figure
    fig = figure('Position', canvas_size);
    hold on;

    % Define colors and markers for subjects
    colors = {hex2rgb('#4A7298'), hex2rgb('#F3C846'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

    % Plot individual subject data
    for s = 1:numSubjects
        plot(uniqueSessions(1:end-5), averageSwappedErrorCorrection(1:end-5, s), ...
            ['-'], 'Color', [colors{s}, 0.5],'LineWidth', 1, 'MarkerSize', 6, ...
            'HandleVisibility', 'off');
		scatter(uniqueSessions(1:end-5), averageSwappedErrorCorrection(1:end-5, s), ...
            markers{s},'MarkerFaceColor', [colors{s}],'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none', ...
            'DisplayName', ['Subject ', subjects{s}(1)]);
    end

    % Plot the smoothed data
    plot(uniqueSessions(1:end-5), smoothedData(1:end-5), '-', 'LineWidth', 2, 'Color', 'k', ...
        'DisplayName', 'Overall Smoothed');

	% Prepare data for linear regression
    sessionNumbers = [];
    allData = [];
    for i = 1:length(uniqueSessions)
        for s = 1:numSubjects
            if ~isnan(averageSwappedErrorCorrection(i, s))
                sessionNumbers = [sessionNumbers; uniqueSessions(i)];
                allData = [allData; averageSwappedErrorCorrection(i, s)];
            end
        end
    end

    % Perform linear regression
    X = [ones(length(sessionNumbers), 1), sessionNumbers];
    [b, ~, ~, ~, stats] = regress(allData, X);
    
    slope = b(2);
    pValue = stats(3);
    rSquared = stats(1);

	% Print regression results
    fprintf('Linear regression results:\n');
    fprintf('Slope = %.4f\n', slope);
    fprintf('p-value = %.4f\n', pValue);
    fprintf('R-squared = %.4f\n', rSquared);
	R2 = rSquared;
    f2 = R2 / (1 - R2);
	fprintf('Cohen''s f2: %.4f\n', f2);

	xLine = [min(uniqueSessions), max(uniqueSessions)-4];
    yLine = b(1) + b(2) * xLine;
    plot(xLine, yLine, '--r', 'LineWidth', 2, 'DisplayName', 'Regression Line');


	 % Calculate overall mean and 95% CI
    allData = averageSwappedErrorCorrection(:);
    overallMean = mean(allData, 'omitnan');
    sem = std(allData, 'omitnan') / sqrt(sum(~isnan(allData)));
    ci95 = 1.96 * sem;
	% Perform one-sample t-test
	[~, p_value, ~, stats] = ttest(allData, 0, 'Alpha', 0.05);
	
	% Print t-test results
	fprintf('One-sample t-test results:\n');
	fprintf('t(%d) = %.4f, p = %.10f\n', stats.df, stats.tstat, p_value);
    % Add triangle with error bar
    maxSession = max(uniqueSessions);
    triangleX = maxSession -2;
    errorbar(triangleX, overallMean, ci95, 'k', 'LineWidth', 1.5, 'CapSize', 10, HandleVisibility='off');
    scatter(triangleX, overallMean, 50, 'k', 'filled', '<', 'DisplayName', 'Overall Mean');
	
	fprintf('%.2f ± %.2f\n', overallMean, ci95);

    % Plot the 0 line on y axis
    yline(0, '--k', 'LineWidth', 1, 'HandleVisibility','off');
	xlim([min(uniqueSessions), maxSession + 2]);

    % Customize plot appearance
    xlabel('Session Number', 'FontSize', 12);
    ylabel('Average Anticipated Swapping', 'FontSize', 12);
    title('Anticipated Swapping Over Sessions', 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Create legend and move it outside
    legend_obj = legend('show', 'Location', 'eastoutside');
    set(legend_obj, 'Box', 'off');
	set(gcf,'renderer','Painters')
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end



function plotOrdinalPositionInferenceOverSession(blockTable)
    % Define canvas size
    canvas_size = [100, 100, 400, 250];

    % Extract unique subjects and sessions
    subjects = unique(blockTable.Subject);
    numSubjects = numel(subjects);
    uniqueSessions = unique(blockTable.SessionNumber);

    % Initialize array for average swappedErrorCorrection
    averageSwappedErrorCorrection = nan(length(uniqueSessions), numSubjects);

    % Calculate average swappedErrorCorrection for each session and subject
    for s = 1:numSubjects
        subject = subjects{s};
        for i = 1:length(uniqueSessions)
            sessionNumber = uniqueSessions(i);
            sessionData = blockTable.SwappedErrorCorrection(blockTable.SessionNumber == sessionNumber & ...
                blockTable.isSwapped == 1 & strcmp(blockTable.Subject, subject));
            averageSwappedErrorCorrection(i, s) = nanmean(sessionData);
        end
    end

    % Apply new averaging rule
	
	smoothedData = smoothDataOverDays(averageSwappedErrorCorrection);
%     smoothedData = smoothDataBySubjectCount(averageSwappedErrorCorrection);

    % Create figure
    fig = figure('Position', canvas_size);
    hold on;

    % Define colors and markers for subjects
    colors = {hex2rgb('#4A7298'), hex2rgb('#F3C846'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

    % Plot individual subject data
    for s = 1:numSubjects
        plot(uniqueSessions(1:end-5), averageSwappedErrorCorrection(1:end-5, s), ...
            ['-'], 'Color', [colors{s}, 0.5],'LineWidth', 1, 'MarkerSize', 6, ...
            'HandleVisibility', 'off');
		scatter(uniqueSessions(1:end-5), averageSwappedErrorCorrection(1:end-5, s), ...
            markers{s},'MarkerFaceColor', [colors{s}],'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none', ...
            'DisplayName', ['Subject ', subjects{s}(1)]);
    end

    % Plot the smoothed data
    plot(uniqueSessions(1:end-5), smoothedData(1:end-5), '-', 'LineWidth', 2, 'Color', 'k', ...
        'DisplayName', 'Overall Smoothed');

	% Prepare data for linear regression
    sessionNumbers = [];
    allData = [];
    for i = 1:length(uniqueSessions)
        for s = 1:numSubjects
            if ~isnan(averageSwappedErrorCorrection(i, s))
                sessionNumbers = [sessionNumbers; uniqueSessions(i)];
                allData = [allData; averageSwappedErrorCorrection(i, s)];
            end
        end
    end

    % Perform linear regression
    X = [ones(length(sessionNumbers), 1), sessionNumbers];
    [b, ~, ~, ~, stats] = regress(allData, X);
    
    slope = b(2);
    pValue = stats(3);
    rSquared = stats(1);

	% Print regression results
    fprintf('Linear regression results:\n');
    fprintf('Slope = %.4f\n', slope);
    fprintf('p-value = %.4f\n', pValue);
    fprintf('R-squared = %.4f\n', rSquared);
	 % Calculate Cohen's f2
    R2 = rSquared;
    f2 = R2 / (1 - R2);
	fprintf('Cohen''s f2: %.4f\n', f2);

	xLine = [min(uniqueSessions), max(uniqueSessions)-4];
    yLine = b(1) + b(2) * xLine;
    plot(xLine, yLine, '--r', 'LineWidth', 2, 'DisplayName', 'Regression Line');


	 % Calculate overall mean and 95% CI
    allData = averageSwappedErrorCorrection(:);
    overallMean = mean(allData, 'omitnan');
    sem = std(allData, 'omitnan') / sqrt(sum(~isnan(allData)));
    ci95 = 1.96 * sem;
	% Perform one-sample t-test
	[~, p_value, ~, stats] = ttest(allData, 0, 'Alpha', 0.05);
	
	% Print t-test results
	fprintf('One-sample t-test results:\n');
	fprintf('t(%d) = %.4f, p = %.10f\n', stats.df, stats.tstat, p_value);
    % Add triangle with error bar
    maxSession = max(uniqueSessions);
    triangleX = maxSession -2;
    errorbar(triangleX, overallMean, ci95, 'k', 'LineWidth', 1.5, 'CapSize', 10, HandleVisibility='off');
    scatter(triangleX, overallMean, 50, 'k', 'filled', '<', 'DisplayName', 'Overall Mean');
	
	fprintf('%.2f ± %.2f\n', overallMean, ci95);

    % Plot the 0 line on y axis
    yline(0, '--k', 'LineWidth', 1, 'HandleVisibility','off');
	xlim([min(uniqueSessions), maxSession + 2]);

    % Customize plot appearance
    xlabel('Session Number', 'FontSize', 12);
    ylabel('Average Ordinal Object-identity Inference ', 'FontSize', 12);
    title('Average Ordinal Object-identity Inference  Over Sessions', 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Create legend and move it outside
    legend_obj = legend('show', 'Location', 'eastoutside');
    set(legend_obj, 'Box', 'off');
	set(gcf,'renderer','Painters')
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function smoothedData = smoothDataOverDays(data)
    [numSessions, ~] = size(data);
    smoothedData = nan(numSessions, 1);
    
    for i = 1:numSessions
        % Define the 5-day window
        startIdx = max(1, i - 2);
        endIdx = min(numSessions, i + 2);
        
        % Extract data for the current window
        windowData = data(startIdx:endIdx, :);
        
        % Calculate the mean across all subjects and days in the window
        smoothedData(i) = mean(windowData(:), 'omitnan');
    end
end

function smoothedData = smoothDataBySubjectCount(data)
    [numSessions, numSubjects] = size(data);
    smoothedData = nan(numSessions, 1);
    
    for i = 1:numSessions
        validData = data(i, ~isnan(data(i, :)));
        numValidSubjects = length(validData);
        
        if numValidSubjects >= 3
            smoothedData(i) = mean(validData);
        elseif numValidSubjects == 2
            if i < numSessions
                nextDayData = data(i+1, ~isnan(data(i+1, :)));
                smoothedData(i) = mean([validData, nextDayData]);
            else
                smoothedData(i) = mean(validData);
            end
        elseif numValidSubjects == 1
            endIdx = min(i+4, numSessions);
            nextFiveDays = data(i:endIdx, :);
            smoothedData(i) = mean(nextFiveDays(:), 'omitnan');
        end
    end
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

function plotStandardizedOrdinalInference(touch_table)
    % Extract unique subjects
    subjects = unique(touch_table.Subject);
    numSubjects = numel(subjects);
    
    % Set up the figure
    canvas_size = [100, 100, 400, 300];
    fig = figure('Position', canvas_size);
    hold on;
    
    % Define colors for subjects
    colors = lines(numSubjects);
    markers = {'o', 's', 'd', '^'};
    
    % Initialize arrays for all normalized data
    allNormalizedSessions = [];
    allSwappedErrorCorrection = [];
    
    % Initialize cell array to store individual fit results
    fitResults = cell(numSubjects, 1);
    
    % Plot individual subject data and fit linear trends
    for s = 1:numSubjects
        subject = subjects{s};
        subjectData = touch_table(strcmp(touch_table.Subject, subject), :);
        
        % Get unique sessions
        sessions = unique(subjectData.SessionNumber);
        
        % Initialize arrays for this subject
        sessionSwappedErrorCorrection = nan(length(sessions), 1);
        
        % Calculate SwappedErrorCorrection for each session
        for i = 1:length(sessions)
            sessionData = subjectData(subjectData.SessionNumber == sessions(i), :);
            
            transitionProbs = zeros(1, 6);  % A, B, C, D, E, Distractor
            totalTrials = 0;
            
            trials = unique(sessionData.TrialNumberInBlock);
            for trial = trials'
                trialData = sessionData(sessionData.TrialNumberInBlock == trial, :);
                
                if trialData.isSwapped(1) == 1
                    for touch = 1:height(trialData) - 1
                        if trialData.CurrentState(touch) == 1 && ...
                           (strcmpi(trialData.TouchCategory{touch}, 'correct') || strcmpi(trialData.TouchCategory{touch}, 'correctSelection')) && ...
                           trialData.CurrentState(touch+1) == 0 && ...
                           trialData.TouchObjectCorrectOrdinalPosition(touch+1) == 4  % Erroneously choosing B
                            
                            nextTouchIndex = touch + 2;
                            if nextTouchIndex <= height(trialData)
                                nextErrorObjectPosition = trialData.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                if nextErrorObjectPosition == 0
                                    nextErrorObjectPosition = 6;  % Distractor
                                end
                                transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                totalTrials = totalTrials + 1;
                            end
                        end
                    end
                end
            end
            
            if totalTrials > 0
                transitionProbs = transitionProbs / totalTrials;
                sessionSwappedErrorCorrection(i) = transitionProbs(4) - transitionProbs(3);  % C - D
            end
        end
        
        % Remove NaN values
        validIdx = ~isnan(sessionSwappedErrorCorrection);
        sessions = sessions(validIdx);
        sessionSwappedErrorCorrection = sessionSwappedErrorCorrection(validIdx);
        
        if isempty(sessionSwappedErrorCorrection)
            fprintf('No valid data for Monkey %s\n', subject(1));
            continue;
        end
        
        % Average over three-session windows
        windowSize = 3;
        avgSessions = [];
        avgSwappedErrorCorrection = [];
        for i = 1:windowSize:length(sessions)
            windowEnd = min(i+windowSize-1, length(sessions));
            avgSessions(end+1) = mean(sessions(i:windowEnd));
            avgSwappedErrorCorrection(end+1) = mean(sessionSwappedErrorCorrection(i:windowEnd));
        end
        
        % Normalize session numbers to 0-100% scale
        normalizedSessions = (avgSessions - min(avgSessions)) / (max(avgSessions) - min(avgSessions)) * 100;
        
        % Plot individual subject data
        plot(normalizedSessions, avgSwappedErrorCorrection, ...
            [markers{mod(s-1, numel(markers))+1}], 'Color', colors(s,:), ...
            'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', ['Monkey ' subject(1)]);
        
        % Fit linear trend for this subject
        [p, S] = polyfit(normalizedSessions, avgSwappedErrorCorrection, 1);
        x_fit = [0, 100];
        y_fit = polyval(p, x_fit);
        
        % Plot linear fit
        plot(x_fit, y_fit, '-', 'Color', colors(s,:), 'LineWidth', 1.5);
        
        % Calculate statistics for the fit
        [~, delta] = polyval(p, normalizedSessions, S);
        df = length(normalizedSessions) - 2;  % degrees of freedom
        t_stat = p(1) / (S.normr/sqrt(S.df) * sqrt(1 / sum((normalizedSessions - mean(normalizedSessions)).^2)));
        p_value = 2 * (1 - tcdf(abs(t_stat), df));
        
        % Store fit results
        fitResults{s} = struct('subject', subject, 'slope', p(1), 'intercept', p(2), 'p_value', p_value);
        
        % Store normalized data for overall trend
        allNormalizedSessions = [allNormalizedSessions; normalizedSessions'];
        allSwappedErrorCorrection = [allSwappedErrorCorrection; avgSwappedErrorCorrection'];
    end
    
    % Calculate and plot overall trend
    [sortedNormalizedSessions, sortIdx] = sort(allNormalizedSessions);
    sortedSwappedErrorCorrection = allSwappedErrorCorrection(sortIdx);
    
    % Use a moving average to smooth the overall trend
    windowSize = max(2, ceil(numel(sortedNormalizedSessions) * 0.1)); % 10% of data points, minimum 2
    smoothedSwappedErrorCorrection = movmean(sortedSwappedErrorCorrection, windowSize, 'omitnan');
    
    plot(sortedNormalizedSessions, smoothedSwappedErrorCorrection, '-', 'LineWidth', 3, ...
        'Color', [0, 0, 0], 'DisplayName', 'Overall Trend');
    
    % Perform polynomial fit for overall trend
    p_overall = polyfit(sortedNormalizedSessions, sortedSwappedErrorCorrection, 2); % 2nd degree polynomial
    x_fit_overall = linspace(0, 100, 100);
    y_fit_overall = polyval(p_overall, x_fit_overall);
    
    plot(x_fit_overall, y_fit_overall, '--', 'LineWidth', 1.5, ...
        'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Overall Polynomial Fit');
    
    % Customize plot appearance
    xlabel('Normalized Session Progress (%)', 'FontSize', 12);
    ylabel('Ordinal Object Inference', 'FontSize', 12);
    title('Standardized Ordinal Object Inference Learning', 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');
    
    % Create legend and move it outside
    legend_obj = legend('show', 'Location', 'eastoutside');
    set(legend_obj, 'Box', 'off');
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    new_fig_size = [canvas_size(1), canvas_size(2), canvas_size(3) + legend_pos(3), canvas_size(4)];
    set(fig, 'Position', new_fig_size);
    
    hold off;
    
    % Display individual fit results
    fprintf('Individual Subject Linear Fit Results (3-session averaged):\n');
    for s = 1:numSubjects
        result = fitResults{s};
        if ~isempty(result)
            fprintf('Monkey %s: Slope = %.4f, Intercept = %.4f, p-value = %.4f\n', ...
                result.subject(1), result.slope, result.intercept, result.p_value);
        end
    end
    
    % Calculate and display overall statistics
    [r, p] = corr(sortedNormalizedSessions, sortedSwappedErrorCorrection, 'Type', 'Spearman', 'Rows', 'complete');
    fprintf('\nOverall Statistics:\n');
    fprintf('Spearman correlation coefficient: %.4f\n', r);
    fprintf('p-value: %.4f\n', p);
    
    % Calculate R-squared for overall polynomial fit
    yresid = sortedSwappedErrorCorrection - polyval(p_overall, sortedNormalizedSessions);
    SSresid = sum(yresid.^2, 'omitnan');
    SStotal = (sum(~isnan(sortedSwappedErrorCorrection))-1) * var(sortedSwappedErrorCorrection, 'omitnan');
    rsq = 1 - SSresid/SStotal;
    fprintf('R-squared for overall polynomial fit: %.4f\n', rsq);
end

function glmm_analysis(block_table)
    % Get unique subjects
    subjects = unique(block_table.Subject);
    numSubjects = length(subjects);
    
    % Initialize table for modeling
    modelData = table();
    
    for subj = 1:numSubjects
        % Filter data for current subject
        subject_data = block_table(strcmp(block_table.Subject, subjects{subj}), :);
        
        % Get unique sessions for this subject
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            
            % Find repeated blocks
            repeated_blocks = session_data(session_data.isRepetition == 1 & session_data.isNewRepetition == 0, :);
            
            for i = 1:height(repeated_blocks)
                repeat_block = repeated_blocks(i, :);
                initial_block = session_data(session_data.BlockNumber == repeat_block.InitialBlockNumber, :);
                
                % Calculate accuracy difference
                acc_diff = repeat_block.BlockLevelAccuracy - initial_block.BlockLevelAccuracy;
                
                % Determine if background and color are the same
                is_same = initial_block.isColorSame == 1 && initial_block.isBackgroundSame == 1;
                
                % Add row to modelData
                newRow = table(subjects(subj), session, repeat_block.BlockNumber, ...
                               repeat_block.InitialBlockNumber, repeat_block.BlockLevelAccuracy, ...
                               initial_block.BlockLevelAccuracy, acc_diff, is_same, ...
                               'VariableNames', {'Subject', 'Session', 'RepeatBlock', ...
                               'InitialBlock', 'AccuracyRepeat', 'AccuracyInitial', ...
                               'AccuracyDifference', 'IsSame'});
                modelData = [modelData; newRow];
            end
        end
    end
    
    % Convert Subject to categorical
    modelData.Subject = categorical(modelData.Subject);
    modelData.IsSame = modelData.IsSame * 2 - 1;
    % Fit the lme
    lme = fitlme(modelData, 'AccuracyDifference ~ 1 + IsSame + (1|Subject)');
                   
    % Display the results
    disp(lme);
    
    % Perform likelihood ratio test for significance of fixed effect
    reducedModel = fitlme(modelData, 'AccuracyDifference ~ 1 + (1|Subject)');
    comparison = compare(reducedModel, glme, 'CheckNesting', true);
    
    fprintf('\nLikelihood Ratio Test:\n');
    disp(comparison);
    
    % Calculate and display effect size (Cohen's d)
    same_data = modelData.AccuracyDifference(modelData.IsSame);
    different_data = modelData.AccuracyDifference(~modelData.IsSame);
    pooled_std = sqrt((std(same_data)^2 + std(different_data)^2) / 2);
    cohens_d = (mean(same_data) - mean(different_data)) / pooled_std;
    
    fprintf('\nEffect Size:\n');
    fprintf('Cohen''s d: %.4f\n', cohens_d);
end

function plot_anticipated_swapping_difference(block_table)
    % Get unique subjects
    subjects = unique(block_table.Subject);
    numSubjects = length(subjects);
    
    % Initialize arrays to store differences
    diff_same = [];
    diff_different = [];
    subject_indices_same = [];
     subject_indices_diff = [];
    % Define shapes for subjects
    shapes = {'o', 's', 'd', '^'};
    
    % Define colors
    color_different = [20, 52, 97] / 255;
    color_same = [213, 197, 132] / 255;
    
    for subj = 1:numSubjects
        % Filter data for current subject
        subject_data = block_table(strcmp(block_table.Subject, subjects{subj}), :);
        
        % Get unique sessions for this subject
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            
            % Find repeated blocks
            repeated_blocks = session_data(session_data.isRepetition == 1 & session_data.isNewRepetition == 0 & session_data.isSwapped == 1, :);
            
            for i = 1:height(repeated_blocks)
                repeat_block = repeated_blocks(i, :);
                initial_block = session_data(session_data.BlockNumber == repeat_block.InitialBlockNumber & session_data.isRepetition == 0, :);
                
                if ~isempty(initial_block)
					if ~isnan(initial_block.isColorSame) && ~isnan(initial_block.isBackgroundSame)
                    	% Calculate SwapErrorCorrection difference
                    	sec_diff = repeat_block.AnticipatedSwapping - initial_block.AnticipatedSwapping;
                    	
                    	% Categorize based on color and background similarity
                		if initial_block.isColorSame == 1 && initial_block.isBackgroundSame == 1
                    		diff_same = [diff_same; sec_diff];
                    		subject_indices_same = [subject_indices_same; subj];
%                 		elseif initial_block.isColorSame == 0 || initial_block.isBackgroundSame == 0
						else
                    		diff_different = [diff_different; sec_diff];
                    		subject_indices_diff = [subject_indices_diff; subj];
						end
					end
                end
            end
        end
    end
    
    % Prepare data for plotting
    all_data = [diff_different; diff_same];
    group = [ones(size(diff_different)); 2*ones(size(diff_same))];
    subject_indices = [subject_indices_diff; subject_indices_same];
    % Set up the figure
    canvas_size = [100, 100, 280, 220];
    fig1 = figure('Position', canvas_size);
    set(gcf,'renderer','Painters')
    
    hold on;
    
    % Add individual data points
    for i = 1:2
        for s = 1:numSubjects
            x = ones(height(all_data(group==i & subject_indices==s)), 1) * i;
            y = all_data(group==i & subject_indices==s);
            if i == 1
            	color = color_different;
        	else
            	color = color_same;
        	end
        	sm = swarmchart(x, y, 20, shapes{s}, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5, 'Jitter', 'off');
        	sm.XJitterWidth = 0.5;
        end
    end
    
	% Create the box plot
    boxplot(all_data, group, 'Labels', {'Different', 'Same'}, 'Colors', [color_different; color_same], 'Symbol', '');
    set(findobj(gca,'type','line'),'linew',1.5);

    % Add error bars
    errorbar([1 2], [nanmean(diff_different) nanmean(diff_same)], ...
             [1.96* nanstd(diff_different)/sqrt(length(~isnan(diff_different))) 1.96* nanstd(diff_same)/sqrt(length(~isnan(diff_same)))], ...
             'k', 'LineStyle', 'none', 'LineWidth', 1);
    plot([1 2], [nanmean(diff_different) nanmean(diff_same)], 'k');
	fprintf('Diff:%.4f ± %.4f ; Same: %.4f±  %.4f \n', nanmean(diff_different), nanmean(diff_same), 1.96* nanstd(diff_different)/sqrt(length(~isnan(diff_different))), 1.96* nanstd(diff_same)/sqrt(length(~isnan(diff_same))));
    % Perform Welch's t-test
    [~, p, ~, stats] = ttest2(diff_different(~isnan(diff_different)), diff_same(~isnan(diff_same)), 'Vartype', 'unequal');
    
    % Add significance indicator
    y_max = max(ylim);
    plot([1, 2], [y_max, y_max], 'k-');
    if p < 0.001
        text(1.5, y_max*1.1, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.01
        text(1.5, y_max*1.1, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.05
        text(1.5, y_max*1.1, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.1, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    % Customize the plot
    ylabel('Anticipated Swapping Difference', 'FontSize', 12);
    title('Anticipated Swapping Difference: Different vs Same Context', 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    % Create dummy plots for legend
    subjectHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        subjectHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end
    
    % Create legend for subjects only
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    
    % Set y-axis limits (adjust as needed)
    ylim([-1, 1]);
	yticks(-1:0.5:1);
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);
    
    hold off;
    
    % Display t-test results
    fprintf('Welch''s t-test results:\n');
    fprintf('t(%0.2f) = %0.2f, p = %0.4f\n', stats.df, stats.tstat, p);
end


function plot_ordinal_identity_inference_difference(block_table)
    % Get unique subjects
    subjects = unique(block_table.Subject);
    numSubjects = length(subjects);
    
    % Initialize arrays to store differences
    diff_same = [];
    diff_different = [];
    subject_indices_same = [];
     subject_indices_diff = [];
    % Define shapes for subjects
    shapes = {'o', 's', 'd', '^'};
    
    % Define colors
    color_different = [20, 52, 97] / 255;
    color_same = [213, 197, 132] / 255;
    
    for subj = 1:numSubjects
        % Filter data for current subject
        subject_data = block_table(strcmp(block_table.Subject, subjects{subj}), :);
        
        % Get unique sessions for this subject
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            
            % Find repeated blocks
            repeated_blocks = session_data(session_data.isRepetition == 1 & session_data.isNewRepetition == 0 & session_data.isSwapped == 1, :);
            
            for i = 1:height(repeated_blocks)
                repeat_block = repeated_blocks(i, :);
                initial_block = session_data(session_data.BlockNumber == repeat_block.InitialBlockNumber & session_data.isRepetition == 0, :);
                
                if ~isempty(initial_block)
					if ~isnan(initial_block.isColorSame) && ~isnan(initial_block.isBackgroundSame)
                    	% Calculate SwapErrorCorrection difference
                    	sec_diff = repeat_block.SwappedErrorCorrection - initial_block.SwappedErrorCorrection;
                    	
                    	% Categorize based on color and background similarity
                		if initial_block.isColorSame == 1 && initial_block.isBackgroundSame == 1
                    		diff_same = [diff_same; sec_diff];
                    		subject_indices_same = [subject_indices_same; subj];
%                 		elseif initial_block.isColorSame == 0 || initial_block.isBackgroundSame == 0
						else
                    		diff_different = [diff_different; sec_diff];
                    		subject_indices_diff = [subject_indices_diff; subj];
						end
					end
                end
            end
        end
    end
    
    % Prepare data for plotting
    all_data = [diff_different; diff_same];
    group = [ones(size(diff_different)); 2*ones(size(diff_same))];
    subject_indices = [subject_indices_diff; subject_indices_same];
    % Set up the figure
    canvas_size = [100, 100, 280, 220];
    fig1 = figure('Position', canvas_size);
    
    
    hold on;
    
    % Add individual data points
    for i = 1:2
        for s = 1:numSubjects
            x = ones(height(all_data(group==i & subject_indices==s)), 1) * i;
            y = all_data(group==i & subject_indices==s);
            if i == 1
            	color = color_different;
        	else
            	color = color_same;
        	end
        	sm = swarmchart(x, y, 20, shapes{s}, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5, 'Jitter', 'off');
        	sm.XJitterWidth = 0.5;
        end
    end
    
    % Add error bars
    errorbar([1 2], [nanmean(diff_different) nanmean(diff_same)], ...
             [1.96* nanstd(diff_different)/sqrt(length(~isnan(diff_different))) 1.96* nanstd(diff_same)/sqrt(length(~isnan(diff_same)))], ...
             'k', 'LineStyle', 'none', 'LineWidth', 1);
	plot([1 2], [nanmean(diff_different) nanmean(diff_same)], 'k', LineWidth=0.5);
    fprintf('Diff:%.4f ± %.4f ; Same: %.4f±  %.4f \n', nanmean(diff_different), nanmean(diff_same), 1.96* nanstd(diff_different)/sqrt(length(~isnan(diff_different))), 1.96* nanstd(diff_same)/sqrt(length(~isnan(diff_same))));
    % Perform Welch's t-test
    [~, p, ~, stats] = ttest2(~isnan(diff_different), ~isnan(diff_same), 'Vartype', 'unequal');
    % Create the box plot
    boxplot(all_data, group, 'Labels', {'Different', 'Same'}, 'Colors', [color_different; color_same], 'Symbol', '');
    set(findobj(gca,'type','line'),'linew',1.5);
	set(gcf,'renderer','Painters')
    % Add significance indicator
    y_max = max(ylim);
    plot([1, 2], [y_max, y_max], 'k-');
    if p < 0.001
        text(1.5, y_max*1.1, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.01
        text(1.5, y_max*1.1, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.05
        text(1.5, y_max*1.1, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.1, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    % Customize the plot
    ylabel('Ordinal Identity Inference Difference', 'FontSize', 12);
    title('Ordinal Identity Inference Difference: Different vs Same Context', 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    % Create dummy plots for legend
    subjectHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        subjectHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end
    
    % Create legend for subjects only
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    
    % Set y-axis limits (adjust as needed)
    ylim([-2.2, 2.6]);
	yticks(-2:1:2);
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);
    
    hold off;
    
    % Display t-test results
    fprintf('Welch''s t-test results:\n');
    fprintf('t(%0.2f) = %0.2f, p = %0.4f\n', stats.df, stats.tstat, p);
end

function plot_block_accuracy_difference(block_table)
    % Get unique subjects
    subjects = unique(block_table.Subject);
    numSubjects = length(subjects);
    
    % Initialize arrays to store differences
    diff_same = [];
    diff_different = [];
    subject_indices = [];
	subject_indices_same = [];
	subject_indices_diff = [];

    diff_same_session = [];
    diff_different_session = [];
    subject_indices_session = [];
    
    % Define shapes for subjects
    shapes = {'o', 's', 'd', '^'};
    
    % Define colors
    color_different = [20, 52, 97] / 255;
    color_same = [213, 197, 132] / 255;
    combined_color = [color_different; color_same];
    
    for subj = 1:numSubjects
        % Filter data for current subject
        subject_data = block_table(strcmp(block_table.Subject, subjects{subj}), :);
        
        % Get unique sessions for this subject
        sessions = unique(subject_data.SessionNumber);
        
        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            
            % Find repeated blocks
            repeated_blocks = session_data(session_data.isRepetition == 1 & session_data.isNewRepetition == 0, :);
            
            session_diff_same = [];
            session_diff_different = [];
            
            for i = 1:height(repeated_blocks)
                repeat_block = repeated_blocks(i, :);
                initial_block = session_data(session_data.BlockNumber == repeat_block.InitialBlockNumber, :);
                
                % Calculate accuracy difference
                acc_diff = repeat_block.BlockLevelAccuracy - initial_block.BlockLevelAccuracy;
                
                % Categorize based on color and background similarity
                if initial_block.isColorSame == 1 && initial_block.isBackgroundSame == 1
                    diff_same = [diff_same; acc_diff];
                    subject_indices_same = [subject_indices_same; subj];
                    session_diff_same = [session_diff_same; acc_diff];
				else
%                 elseif initial_block.isColorSame == 0 || initial_block.isBackgroundSame == 0
                    diff_different = [diff_different; acc_diff];
                    subject_indices_diff = [subject_indices_diff; subj];
                    session_diff_different = [session_diff_different; acc_diff];
                end
            end
            
            diff_same_session = [diff_same_session; nanmean(session_diff_same)];
            diff_different_session = [diff_different_session; nanmean(session_diff_different)];
            subject_indices_session = [subject_indices_session; subj];
        end
	end

	   % Create Figure 1 (Trial-wise)
    create_plot(diff_different, diff_same, subject_indices_same,subject_indices_diff, 'Accuracy Difference: Repeated vs Initial Blocks');
	fprintf('Number of Blocks: Diff: %d, Same: %d, Total: %d\n', length(diff_different), length(diff_same), length(diff_different)+ length(diff_same));
%     % Create Figure 2 (Session-wise)
%     create_plot(diff_different_session, diff_same_session, subject_indices_session,subject_indices_session, 'Accuracy Difference: Repeated vs Initial Block ()');


    
	function create_plot(data_different, data_same, subj_indices_same, subj_indices_diff, fig_title)
    % Remove NaN values

	   [data_different, subj_indices_diff] = remove_outliers(data_different, subj_indices_diff);
       [data_same, subj_indices_same] = remove_outliers(data_same, subj_indices_same);


    % Prepare data for plotting
    all_data = [data_different; data_same];
    group = [ones(size(data_different)); 2*ones(size(data_same))];
    
    % Set up the figure
    canvas_size = [100, 100, 280, 220];
    fig = figure('Position', canvas_size);
    
    
    hold on;
	shapes = {'o', 's', 'd', '^'};
    
    % Add individual data points
    for i = 1:2
        if i == 1
            data = data_different;
            subj_indices_plot = subj_indices_diff;
        else
            data = data_same;
            subj_indices_plot = subj_indices_same;
        end
        for s = 1:numSubjects
            x = ones(height(data(subj_indices_plot == s)), 1) * i;
            y = data(subj_indices_plot == s);
            sm = swarmchart(x, y, 20, shapes{s},  'MarkerEdgeColor','none','MarkerFaceColor', combined_color(i,:), 'MarkerFaceAlpha', 0.5);
            sm.XJitterWidth = 0.5;
        end
    end
    
    % Add error bars
    errorbar([1 2], [mean(data_different) mean(data_same)], ...
             [1.96*std(data_different)/sqrt(length(data_different)) 1.96*std(data_same)/sqrt(length(data_same))], ...
             'k', 'LineStyle', 'none', 'LineWidth', 1, CapSize=10);
	plot([1,2], [mean(data_different) mean(data_same)], 'k', 'LineWidth', 1);
    fprintf('Diff %.4f ± %.4f  Same %.4f ± %.4f\n', mean(data_different),1.96*std(data_different)/sqrt(length(data_different)), mean(data_same), 1.96*std(data_same)/sqrt(length(data_same)))
    % Perform Welch's t-test
    [~, p, ~, stats] = ttest2(data_different, data_same, 'Vartype', 'unequal');
    

	% Create the box plot
    boxplot(all_data, group, 'Labels', {'Different', 'Same'}, 'Colors', [color_different; color_same], 'Symbol','');
    set(findobj(gca,'type','line'),'linew',1.5);
	set(gcf,'renderer','Painters')
    % Add significance indicator
    y_max = max(ylim);
    plot([1, 2], [y_max, y_max], 'k-');
    if p < 0.001
        text(1.5, y_max*1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.01
        text(1.5, y_max*1.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.05
        text(1.5, y_max*1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.05, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    % Customize the plot
    ylabel('Accuracy Difference', 'FontSize', 12);
    title(fig_title, 'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Create dummy plots for legend
    subjectHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        subjectHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end
    
    % Create legend for subjects only
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    ylim([min(all_data) * 1.2, max(all_data) * 1.2]);
    
    hold off;
    
    % Display t-test results
    fprintf('Welch''s t-test results for %s:\n', fig_title);
    fprintf('t(%0.2f) = %0.2f, p = %0.4f\n', stats.df, stats.tstat, p);

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
end
    % Function to remove outliers
    function [cleaned_data, cleaned_indices] = remove_outliers(data, indices)
        Q1 = prctile(data, 25);
        Q3 = prctile(data, 75);
        IQR = Q3 - Q1;
        lower_bound = Q1 - 1.5 * IQR;
        upper_bound = Q3 + 1.5 * IQR;
        
        valid_data = data >= lower_bound & data <= upper_bound;
        cleaned_data = data(valid_data);
        cleaned_indices = indices(valid_data);
    end


 
	
end

function plotSessionWiseTransitionProbability(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Define transition states
    transitionStates = {'A', 'B', 'C', 'D', 'E', 'Distractor'};
    numStates = length(transitionStates);

    % Define markers and colors for subjects
    markers = {'o', 's', 'd', '^', 'v'};
    colors = [111, 34, 45;
              188, 89, 111;
              133, 71, 135;
              73, 101, 61;
              64, 128, 184;
              34, 74, 137] / 255; 
	
    % Initialize array to store session-wise transition probabilities
    allSessionProbs = cell(numSubjects, 1);

    % Calculate transition probabilities for each subject and session
    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subject_data.SessionNumber);
        subjectProbs = [];

        for session = sessions'
            session_data = subject_data(subject_data.SessionNumber == session, :);
            transitionProbs = zeros(1, numStates);
            totalTrials = 0;

            blocks = unique(session_data.BlockNumber);
            for block = blocks'
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);

                for trial = trials'
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);

                    if trial_data.isSwapped(1) == 1
                        for touch = 1:height(trial_data) - 1
                            if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection')) ...
                                    && trial_data.CurrentState(touch+1) == 0 && trial_data.TouchObjectCorrectOrdinalPosition(touch+1) == 4

                                nextTouchIndex = touch + 2;
                                if nextTouchIndex + 1 <= height(trial_data) && strcmpi(trial_data.TouchCategory{nextTouchIndex}, 'RetouchCorrect')
                                    nextTouchIndex = nextTouchIndex + 1;
                                end

                                if nextTouchIndex <= height(trial_data)
                                    nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
                                    if nextErrorObjectPosition == 0
                                        nextErrorObjectPosition = 6;
                                    elseif nextErrorObjectPosition == 2
                                        nextErrorObjectPosition = 4;
                                    elseif nextErrorObjectPosition == 4
                                        nextErrorObjectPosition = 2;
                                    end
                                    transitionProbs(nextErrorObjectPosition) = transitionProbs(nextErrorObjectPosition) + 1;
                                    totalTrials = totalTrials + 1;
                                end
                            end
                        end
                    end
                end
            end

            % Store transition probabilities for this session
            if totalTrials > 0
                subjectProbs = [subjectProbs; transitionProbs / totalTrials];
            end
        end
        allSessionProbs{s} = subjectProbs;
    end

    % Calculate overall mean and SEM
    allProbsCombined = vertcat(allSessionProbs{:});
    meanProbs = mean(allProbsCombined, 1);
    semProbs = std(allProbsCombined, 0, 1) / sqrt(size(allProbsCombined, 1));


    % Create figure
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Plot swarm chart for each transition state
    for i = 1:numStates
        for s = 1:numSubjects
            subjectData = allSessionProbs{s}(:, i);
            sm = swarmchart(repmat(i, size(subjectData, 1), 1), subjectData, 30, ...
                'Marker', markers{s}, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            sm.XJitterWidth = 0.75;
        end
    end

    % Plot error bars
    errorbar(1:numStates, meanProbs, 1.96*semProbs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    fprintf('Transition Probabilities (mean ± 95%% CI):\n');
    for i = 1:numStates
        fprintf('%s: %.3f ± %.3f\n', transitionStates{i}, meanProbs(i), 1.96*semProbs(i));
    end

    % Perform Welch's t-test between C and D (now 3 and 2)
    [~, p_value] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 4), 'Vartype', 'unequal');

    % Add significance star if p < 0.05 and print p-value
    if p_value < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([3, 4], [maxY*1.1, maxY*1.1], 'k-');
		if p_value < 0.001
        text(3.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value < 0.01
			text(3.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value < 0.05
			text(3.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
    end
    fprintf('P-value for Ordinal Position 3 vs 4 comparison: %.10f\n', p_value);

    % Customize plot
    ylabel('Transition Probability', 'FontSize', 12);
    xlabel('Ordinal Position', 'FontSize', 12);
    title('Transition Probability After Erroneously Choosing B', 'FontSize', 14);
    xticks(1:numStates);
    xticklabels(transitionStates);
    xlim([0.5, 6.5]);
	ylim([-0.1, 1.2]);
	yline(0,'k-', HandleVisibility="off");
	yticks(0:0.2:1);

    % Add legend for subjects
    legendHandles = gobjects(1, numSubjects);
    for s = 1:numSubjects
        legendHandles(s) = plot(NaN, NaN, markers{s},'MarkerFaceColor','k', ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineWidth', 1);
    end
    legend_obj = legend(legendHandles, cellfun(@(x) ['Subject ', x(1)], subjects, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';

    % Adjust axis properties
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(ax, 'box', 'off');
	
	set(gcf,'renderer','Painters')
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






function plotBackgroundErrorInference(block_table)
    % Extract data for each condition
    sameBackgroundData = block_table.SwappedErrorCorrection(block_table.isBackgroundSame == 1 & block_table.isColorSame == 1 & block_table.isSwapped == 1);
    differentBackgroundData = block_table.SwappedErrorCorrection((block_table.isBackgroundSame == 0 | block_table.isColorSame == 0) & block_table.isSwapped == 1);
    % Combine data for boxplot
    combinedData = [differentBackgroundData; sameBackgroundData];
    groups = [ones(size(differentBackgroundData)); 2*ones(size(sameBackgroundData))];
    subjects = unique(block_table.Subject);
    numSubjects = length(subjects);
    % Perform Welch's t-test
    [h, p, ci, stats] = ttest2(differentBackgroundData, sameBackgroundData, 'Vartype', 'unequal');
    
    % Create figure
    figure('Position', [100, 100, 600, 400]);
    hold on;
    
    % Create box plot
   boxplot(combinedData, groups, 'Labels', {'Different', 'Same'}, ...
        'Colors', [0.3 0.3 0.3], 'Width', 0.6);
	markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
     colors = {hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A')};
    
    % Add individual data points
    gray_shades = linspace(0.3, 0.7, numSubjects);
	legendEntries = cell(1, length(subjects));
    for i = 1:numSubjects
        subject = subjects{i};
        subjectDataDiff = block_table.SwappedErrorCorrection(strcmp(block_table.Subject, subject) & ...
            block_table.isBackgroundSame == 1 & block_table.isColorSame == 1 & block_table.isSwapped == 1);
        subjectDataSame = block_table.SwappedErrorCorrection(strcmp(block_table.Subject, subject) & ...
            (block_table.isBackgroundSame == 0 | block_table.isColorSame == 0) & block_table.isSwapped == 1);
        legendEntries{i} = ['Subject ', subjects{i}(1)];

        % Convert data to numeric and remove any NaN values
    	subjectDataDiff = subjectDataDiff(~isnan(subjectDataDiff));
    	subjectDataSame = subjectDataSame(~isnan(subjectDataSame));
    	
    	% Plot both conditions with the same handle
    	if ~isempty(subjectDataDiff) || ~isempty(subjectDataSame)
        	s = scatter([ones(size(subjectDataDiff)); 2*ones(size(subjectDataSame))] + (rand(size([subjectDataDiff; subjectDataSame]))-0.5)*0.1, ...
            	[subjectDataDiff; subjectDataSame], 50, ...
            	colors{i},markers{i}, 'LineWidth', 1,'MarkerEdgeColor', colors{i}, 'DisplayName', legendEntries{i});
        	s.MarkerEdgeAlpha = 1;
    	end
	end
    
    % Add mean and SEM
    errorbar([1, 2], [mean(differentBackgroundData, 'omitnan'), mean(sameBackgroundData, 'omitnan')], ...
        [std(differentBackgroundData, 'omitnan')/sqrt(sum(~isnan(differentBackgroundData))), ...
         std(sameBackgroundData, 'omitnan')/sqrt(sum(~isnan(sameBackgroundData)))], ...
        'o', 'MarkerSize', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
        'Color', 'k', 'LineWidth', 2, 'CapSize', 15, HandleVisibility='off');
    
    % Add significance indicator
    maxY = max([differentBackgroundData; sameBackgroundData]);
    line([1, 2], [maxY + 0.1, maxY + 0.1], 'Color', 'k');
    if p < 0.001
        text(1.5, maxY + 0.07, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.01
        text(1.5, maxY + 0.07, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p < 0.05
        text(1.5, maxY + 0.07, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, maxY + 0.07, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
    
    % Formatting the plot
    title('Error Inference: Different vs Same Background', 'FontSize', 16);
    xlabel('Background Condition', 'FontSize', 14);
    ylabel('Error Inference', 'FontSize', 14);
    ylim([min([differentBackgroundData; sameBackgroundData]) - 0.1, maxY + 0.15]);
    legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 12);
    ax = gca;  % Get current axis handle
    ax.FontSize = 14;  % Increase font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    hold off;
    
    % Display statistical results
    fprintf('Welch''s t-test results:\n');
    fprintf('Different vs Same Background: t(%.2f) = %.4f, p = %.4f\n', stats.df, stats.tstat, p);
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

function plotBlockCompletionVsAnticipatedSwapping(block_table)
    block_table = sortrows(block_table, {'Subject', 'SessionNumber', 'BlockNumber'});
    
    % Initialize arrays to store data
    initial_completion_rates = [];
    error_inferences = [];
    subjects = {};
    
    % Define colors and markers for each subject
    colors = {hex2rgb('#4A7298'), hex2rgb('#F3C846'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

	
    
    % Process the data
    for i = 2:height(block_table)
        if block_table.isSwapped(i) == 1 && block_table.BlockNumber(i) > 1
            prev_block = i - 1;
            if strcmp(block_table.Subject(i), block_table.Subject(prev_block)) && ...
               block_table.SessionNumber(i) == block_table.SessionNumber(prev_block) && ...
               block_table.isSwapped(prev_block) == 0
                initial_rate = block_table.BlockLevelAccuracy(prev_block);
                error_inference = block_table.AnticipatedSwapping(i);
                
                % Check for valid data
                if isfinite(initial_rate) && isfinite(error_inference)
                    initial_completion_rates = [initial_completion_rates; initial_rate];
                    error_inferences = [error_inferences; error_inference];
                    subjects{end+1} = block_table.Subject{i};
                end
            end
        end
    end
    
    % Check if we have enough data
    if isempty(initial_completion_rates) || isempty(error_inferences)
        error('No valid data pairs found after filtering.');
    end
    
    % Define canvas size
    canvas_size = [100, 100, 270, 280];
    
    % Create the scatter plot
    fig1 = figure('Position', canvas_size);
    hold on;
    
    unique_subjects = unique(subjects);
    subject_slopes = zeros(length(unique_subjects), 1);
    subject_slope_se = zeros(length(unique_subjects), 1);
	subject_p= zeros(length(unique_subjects), 1);
    for i = 1:length(unique_subjects)
        subject_indices = strcmp(subjects, unique_subjects{i});
        scatter(initial_completion_rates(subject_indices), error_inferences(subject_indices), ...
                30, colors{i}, markers{i},'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none', 'MarkerEdgeAlpha', 0.05, ...
                'DisplayName', ['Subject ' unique_subjects{i}(1)], 'jitter', 'on','jitterAmount',0.01);

		% Add regression line for individual subject
		X = [ones(sum(subject_indices), 1), initial_completion_rates(subject_indices)];
		y = error_inferences(subject_indices);
		b = X \ y;  % Linear regression
		[b, bint,~,~,stats] = regress(y, X);
		subject_slopes(i) = b(2);
        subject_slope_se(i) = (bint(2,2) - bint(2,1)) / 2; % Approximate SE from 95% CI
		subject_p(i) = stats(3);
		x_range = 0:0.01:1;
		y_pred = b(1) + b(2) * x_range;
		plot(x_range, y_pred, 'Color', colors{i}, 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    % Add diagonal line
%     plot([0, 1], [-1, 1], 'k--', 'HandleVisibility', 'off');
    
    % Set axis limits and labels
    xlim([-0.1, 1.1]);
    ylim([-1.1, 2]);
    yticks([-1, -0.5, 0, 0.5, 1]);
    xlabel('Initial Block Completion Rate', 'FontSize', 12);
    ylabel('Anticipated Swapping', 'FontSize', 12);
%     title('Initial Block Completion Rate vs Anticipated Swapping', 'FontSize', 12);
    
    % Perform linear regression
    X = [ones(length(initial_completion_rates), 1), initial_completion_rates];
    y = error_inferences;
    [b, bint, ~, ~, stats] = regress(y, X);
    
    % Plot regression line
    x_range = 0:0.01:1;
    y_pred = b(1) + b(2) * x_range;
    plot(x_range, y_pred, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    % Calculate correlation coefficient and p-value
    [r, p] = corr(initial_completion_rates, error_inferences);
    
    % Display statistics
    text(0.05, 0.95, sprintf('R^2 = %.4f', stats(1)), 'Units', 'normalized', 'FontSize', 12);
    text(0.05, 0.87, sprintf('p-value = %.9f', stats(3)), 'Units', 'normalized', 'FontSize', 12);
%     text(0.05, 0.80, sprintf('Correlation: r = %.4f', r), 'Units', 'normalized', 'FontSize', 12);
    
    % Apply plot rules
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    % Add legend
	subjects = unique(block_table.Subject);
	subjectHandles = gobjects(1, length(subjects));
    for i = 1:length(subjects)
        subjectHandles(i) = plot(NaN, NaN, markers{mod(i-1, length(markers))+1}, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{i}(1)]);
	end

    legend_obj = legend(subjectHandles,'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    set(gcf,'renderer','Painters')
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    
    hold off;
    
    % Print summary statistics
    fprintf('\nLinear Regression Results:\n');
    fprintf('Intercept: %.4f\n', b(1));
    fprintf('Slope: %.4f\n', b(2));
    fprintf('R-squared: %.4f\n', stats(1));
    fprintf('F-statistic vs. constant model: %.4f\n', stats(2));
    fprintf('p-value: %.4f\n', stats(3));
    fprintf('\nCorrelation Analysis:\n');
    fprintf('Correlation coefficient (r): %.4f\n', r);
    fprintf('p-value: %.10f\n', p);
	   % Calculate Cohen's f2
    R2 = stats(1);
    f2 = R2 / (1 - R2);
	fprintf('Cohen''s f2: %.4f\n', f2);
    
    
    % Interpret the results
    fprintf('\nInterpretation:\n');
    if isnan(r) || isnan(p)
        fprintf('Unable to calculate correlation. Please check the debug information above for possible issues with the data.\n');
    elseif p < 0.05
        fprintf('There is a statistically significant ');
        if r > 0
            fprintf('positive');
        else
            fprintf('negative');
        end
        fprintf(' linear correlation between initial block completion rate and current error inference.\n');
    else
        fprintf('There is no statistically significant linear correlation between initial block completion rate and current error inference.\n');
    end
    
    % Additional statistics
    fprintf('\nAdditional Statistics:\n');
    fprintf('Mean Initial Completion Rate: %.4f\n', mean(initial_completion_rates));
    fprintf('Median Initial Completion Rate: %.4f\n', median(initial_completion_rates));
    fprintf('Standard Deviation of Initial Completion Rate: %.4f\n', std(initial_completion_rates));
    fprintf('Mean Error Inference: %.4f\n', mean(error_inferences));
    fprintf('Median Error Inference: %.4f\n', median(error_inferences));
    fprintf('Standard Deviation of Error Inference: %.4f\n', std(error_inferences));

	 % Create new figure for slope estimates
	 canvas_size = [100, 100, 300, 250];
    fig2 = figure('Position', canvas_size);
    hold on;

    for i = 1:length(unique_subjects)
        errorbar(i, subject_slopes(i), subject_slope_se(i), 'Color', colors{i}, 'LineStyle', 'none', 'LineWidth', 1.5, HandleVisibility='off');
        scatter(i, subject_slopes(i), 100, colors{i}, markers{i}, 'filled', ...
                'DisplayName', ['Subject ' unique_subjects{i}(1)]);
		disp(subject_p(i));
    end

    % Plot cross-subjects fit
    errorbar(length(unique_subjects) + 1, b(2), (bint(2) - bint(1)) / (2), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, HandleVisibility='off');
    scatter(length(unique_subjects) + 1, b(2), 100, 'k', 'o', 'DisplayName', 'Cross-subjects fit');

    % Customize the plot
    xlabel('Subject', 'FontSize', 12);
    ylabel('Estimated Slope', 'FontSize', 12);
    title('Estimated Slopes by Subject', 'FontSize', 14);
    xticks(1:(length(unique_subjects) + 1));
    xticklabels([cellfun(@(x) x(1), unique_subjects, 'UniformOutput', false), {'All'}]);
    
    % Add horizontal line at y = 0
    yline(0, 'k--', 'HandleVisibility', 'off');
	xlim([0.2, 5.8]);
    % Apply plot rules
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Add legend
    legend_obj = legend('Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig2, 'Position', new_fig_size);


    hold off;


end


function plotBlockCompletionVsErrorCorrection(block_table)
    block_table = sortrows(block_table, {'Subject', 'SessionNumber', 'BlockNumber'});
    
    % Initialize arrays to store data
    initial_completion_rates = [];
    error_inferences = [];
    subjects = {};
    
    % Define colors and markers for each subject
    colors = {hex2rgb('#4A7298'), hex2rgb('#F3C846'), hex2rgb('#C83E4D'), hex2rgb('#4E937A')};
    markers = {'o', 's', 'd', '^'};

	
    
    % Process the data
    for i = 2:height(block_table)
        if block_table.isSwapped(i) == 1 && block_table.BlockNumber(i) > 1
            prev_block = i - 1;
            if strcmp(block_table.Subject(i), block_table.Subject(prev_block)) && ...
               block_table.SessionNumber(i) == block_table.SessionNumber(prev_block) && ...
               block_table.isSwapped(prev_block) == 0
                initial_rate = block_table.BlockLevelAccuracy(prev_block);
                error_inference = block_table.SwappedErrorCorrection(i);
                
                % Check for valid data
                if isfinite(initial_rate) && isfinite(error_inference)
                    initial_completion_rates = [initial_completion_rates; initial_rate];
                    error_inferences = [error_inferences; error_inference];
                    subjects{end+1} = block_table.Subject{i};
                end
            end
        end
    end
    
    % Check if we have enough data
    if isempty(initial_completion_rates) || isempty(error_inferences)
        error('No valid data pairs found after filtering.');
    end
    
    % Define canvas size
    canvas_size = [100, 100, 270, 280];
    
    % Create the scatter plot
    fig1 = figure('Position', canvas_size);
    hold on;
    
    unique_subjects = unique(subjects);
    subject_slopes = zeros(length(unique_subjects), 1);
    subject_slope_se = zeros(length(unique_subjects), 1);
	subject_p= zeros(length(unique_subjects), 1);
    for i = 1:length(unique_subjects)
        subject_indices = strcmp(subjects, unique_subjects{i});
        scatter(initial_completion_rates(subject_indices), error_inferences(subject_indices), ...
                30, colors{i}, markers{i},'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none', 'MarkerEdgeAlpha', 0.1, ...
                'DisplayName', ['Subject ' unique_subjects{i}(1)], 'jitter', 'on','jitterAmount',0.01);

		% Add regression line for individual subject
		X = [ones(sum(subject_indices), 1), initial_completion_rates(subject_indices)];
		y = error_inferences(subject_indices);
		b = X \ y;  % Linear regression
		[b, bint,~,~,stats] = regress(y, X);
		subject_slopes(i) = b(2);
        subject_slope_se(i) = (bint(2,2) - bint(2,1)) / 2; % Approximate SE from 95% CI
		subject_p(i) = stats(3);
		x_range = 0:0.01:1;
		y_pred = b(1) + b(2) * x_range;
		plot(x_range, y_pred, 'Color', colors{i}, 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    % Add diagonal line
%     plot([0, 1], [-1, 1], 'k--', 'HandleVisibility', 'off');
    
    % Set axis limits and labels
    xlim([-0.1, 1.1]);
    ylim([-1.1, 2]);
    yticks([-1, -0.5, 0, 0.5, 1]);
	set(gcf,'renderer','Painters');
    xlabel('Initial Block Completion Rate', 'FontSize', 12);
    ylabel('Ordinal Obejct-identity Inference', 'FontSize', 12);
    title('Initial Block Completion Rate vs Ordinal Obejct-identity Inference', 'FontSize', 12);
    
    % Perform linear regression
    X = [ones(length(initial_completion_rates), 1), initial_completion_rates];
    y = error_inferences;
    [b, bint, ~, ~, stats] = regress(y, X);
    
    % Plot regression line
    x_range = 0:0.01:1;
    y_pred = b(1) + b(2) * x_range;
    plot(x_range, y_pred, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    % Calculate correlation coefficient and p-value
    [r, p] = corr(initial_completion_rates, error_inferences);
    
    % Display statistics
    text(0.05, 0.95, sprintf('R^2 = %.4f', stats(1)), 'Units', 'normalized', 'FontSize', 12);
    text(0.05, 0.87, sprintf('p-value = %.4f', stats(3)), 'Units', 'normalized', 'FontSize', 12);
    text(0.05, 0.80, sprintf('Correlation: r = %.4f', r), 'Units', 'normalized', 'FontSize', 12);
    
    % Apply plot rules
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    % Add legend
	subjects = unique(block_table.Subject);
	subjectHandles = gobjects(1, length(subjects));
    for i = 1:length(subjects)
        subjectHandles(i) = plot(NaN, NaN, markers{mod(i-1, length(markers))+1}, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{i}(1)]);
	end

    legend_obj = legend(subjectHandles,'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    
    hold off;
    
    % Print summary statistics
    fprintf('\nLinear Regression Results:\n');
    fprintf('Intercept: %.4f\n', b(1));
    fprintf('Slope: %.4f\n', b(2));
    fprintf('R-squared: %.4f\n', stats(1));
    fprintf('F-statistic vs. constant model: %.4f\n', stats(2));
    fprintf('p-value: %.4f\n', stats(3));
    fprintf('\nCorrelation Analysis:\n');
    fprintf('Correlation coefficient (r): %.4f\n', r);
    fprintf('p-value: %.4f\n', p);
	    % Calculate Cohen's f2
    R2 = stats(1);
    f2 = R2 / (1 - R2);
	fprintf('Cohen''s f2: %.4f\n', f2);
    
    % Interpret the results
    fprintf('\nInterpretation:\n');
    if isnan(r) || isnan(p)
        fprintf('Unable to calculate correlation. Please check the debug information above for possible issues with the data.\n');
    elseif p < 0.05
        fprintf('There is a statistically significant ');
        if r > 0
            fprintf('positive');
        else
            fprintf('negative');
        end
        fprintf(' linear correlation between initial block completion rate and current error inference.\n');
    else
        fprintf('There is no statistically significant linear correlation between initial block completion rate and current error inference.\n');
    end
    
    % Additional statistics
    fprintf('\nAdditional Statistics:\n');
    fprintf('Mean Initial Completion Rate: %.4f\n', mean(initial_completion_rates));
    fprintf('Median Initial Completion Rate: %.4f\n', median(initial_completion_rates));
    fprintf('Standard Deviation of Initial Completion Rate: %.4f\n', std(initial_completion_rates));
    fprintf('Mean Error Inference: %.4f\n', mean(error_inferences));
    fprintf('Median Error Inference: %.4f\n', median(error_inferences));
    fprintf('Standard Deviation of Error Inference: %.4f\n', std(error_inferences));

	 % Create new figure for slope estimates
	 canvas_size = [100, 100, 300, 250];
    fig2 = figure('Position', canvas_size);
    hold on;

    for i = 1:length(unique_subjects)
        errorbar(i, subject_slopes(i), subject_slope_se(i), 'Color', colors{i}, 'LineStyle', 'none', 'LineWidth', 1.5, HandleVisibility='off');
        scatter(i, subject_slopes(i), 100, colors{i}, markers{i}, 'filled', ...
                'DisplayName', ['Subject ' unique_subjects{i}(1)]);
		disp(subject_p(i));
    end

    % Plot cross-subjects fit
    errorbar(length(unique_subjects) + 1, b(2), (bint(2) - bint(1)) / (2), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, HandleVisibility='off');
    scatter(length(unique_subjects) + 1, b(2), 100, 'k', 'o', 'DisplayName', 'Cross-subjects fit');

    % Customize the plot
    xlabel('Subject', 'FontSize', 12);
    ylabel('Estimated Slope', 'FontSize', 12);
    title('Estimated Slopes by Subject', 'FontSize', 14);
    xticks(1:(length(unique_subjects) + 1));
    xticklabels([cellfun(@(x) x(1), unique_subjects, 'UniformOutput', false), {'All'}]);
    
    % Add horizontal line at y = 0
    yline(0, 'k--', 'HandleVisibility', 'off');
	xlim([0.2, 5.8]);
    % Apply plot rules
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Add legend
    legend_obj = legend('Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig2, 'Position', new_fig_size);


    hold off;


end




function analyzeSwappedErrorCorrection(block_table)

    % Calculate previous block accuracy
    prev_block_accuracy = nan(height(block_table), 1);
    for i = 2:height(block_table)
        if strcmp(block_table.Subject{i}, block_table.Subject{i-1}) && ...
           block_table.SessionNumber(i) == block_table.SessionNumber(i-1) && ...
           block_table.BlockNumber(i) == block_table.BlockNumber(i-1) + 1
            prev_block_accuracy(i) = block_table.BlockLevelAccuracy(i-1);
        end
	end

    % Filter for swapped blocks only
    swapped_data = block_table(block_table.isSwapped == 1, :);
    
    % Prepare the data for the model
    X = table();
    X.SwappedErrorCorrection = swapped_data.SwappedErrorCorrection;
    X.isRepeat = swapped_data.isRepetition == 1 & swapped_data.isNewRepetition == 0;
    X.isNewInRepeat = swapped_data.isRepetition == 1 & swapped_data.isNewRepetition == 1;
    X.Subject = categorical(swapped_data.Subject);
    X.BlockNumber = swapped_data.BlockNumber;
	X.BlockNumberSquared = swapped_data.BlockNumber.^2;
    X.PrevBlockAccuracy = prev_block_accuracy(block_table.isSwapped == 1);
    
    % Remove rows with NaN values
    X = rmmissing(X);
    
    % Fit the linear mixed-effects model: need to code dependent
	% messurement, make two plot with connected lines for each pair for NEW
	% vs REPEAT, 
    lme = fitlme(X, 'SwappedErrorCorrection ~ isRepeat + isNewInRepeat + PrevBlockAccuracy + (1|Subject)');
    
    % Display the results
    disp('Linear Mixed-Effects Model Results:');
    disp(lme);
    
    % Display the fixed effects estimates
    disp('Fixed Effects Estimates:');
    fixedEffects = lme.fixedEffects;
    disp(fixedEffects);
    
     % Display p-values for fixed effects
    disp('Fixed Effects Statistics:');
    fixedEffectTable = anova(lme);
    disp(fixedEffectTable);
    
    
%     % Plot residuals
%     figure;
%     plotResiduals(lme);
%     title('Residual Plots');
    
%     % QQ plot of residuals
%     figure;
%     qqplot(residuals(lme));
%     title('Q-Q Plot of Residuals');
%     
%     % Plot of fitted values vs residuals
%     figure;
%     scatter(fitted(lme), residuals(lme));
%     xlabel('Fitted Values');
%     ylabel('Residuals');
%     title('Fitted Values vs Residuals');
%     
%     % Histogram of residuals
%     figure;
%     histogram(residuals(lme));
%     xlabel('Residuals');
%     ylabel('Frequency');
%     title('Histogram of Residuals');


end


function plotAnticipatedSwappingThreeConditions(block_table)
    % Define colors for conditions
%     colors = {[80/255, 126/255, 186/255], [209/255, 151/255, 35/255], [179/255, 41/255, 39/255]}; % Blue, Yellow, Red
    colors = {[20, 52, 94] ./ 255,[213, 197, 131] ./ 255,  [213,174, 200]./255};
    % Define shapes for subjects
    shapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
    
    % Get unique subjects
    subjects = unique(block_table.Subject);
    
    % Extract data for each condition
    newEarlyData = block_table(block_table.isRepetition == 0 & block_table.isSwapped == 1, :);
    repeatData = block_table(block_table.isRepetition == 1 & block_table.isNewRepetition == 0 & block_table.isSwapped == 1, :);
    newLaterData = block_table(block_table.isRepetition == 1 & block_table.isNewRepetition == 1 & block_table.isSwapped == 1, :);
    
    % Combine data for boxplot
    allData = [newEarlyData.AnticipatedSwapping; repeatData.AnticipatedSwapping; newLaterData.AnticipatedSwapping];
    groups = [ones(size(newEarlyData, 1), 1); 2*ones(size(repeatData, 1), 1); 3*ones(size(newLaterData, 1), 1)];
    
    % Define canvas size
    canvas_size = [100, 100, 300, 300];
    
    % Create figure
    fig = figure('Position', canvas_size);
    hold on;
    

    
    % Add individual data points for each subject
    for i = 1:length(subjects)
        subject = subjects{i};
        shape = shapes{mod(i-1, length(shapes))+1};
        
        % New Early condition
        subjectNewEarlyData = newEarlyData(strcmp(newEarlyData.Subject, subject), :);
        x1 = ones(size(subjectNewEarlyData, 1), 1);
        s1 = swarmchart(x1, subjectNewEarlyData.AnticipatedSwapping, 30, shape,  'MarkerEdgeColor', 'none','MarkerFaceColor', colors{1}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3);
%         s1.XJitter = 'rand';
        s1.XJitterWidth = 0.75;

        % Repeat condition
        subjectRepeatData = repeatData(strcmp(repeatData.Subject, subject), :);
        x2 = 2*ones(size(subjectRepeatData, 1), 1);
        s2 = swarmchart(x2, subjectRepeatData.AnticipatedSwapping, 30, shape, 'MarkerEdgeColor', 'none','MarkerFaceColor', colors{2}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3);
        s2.XJitter = 'rand';
        s2.XJitterWidth = 0.75;

        % New Later condition
        subjectNewLaterData = newLaterData(strcmp(newLaterData.Subject, subject), :);
        x3 = 3*ones(size(subjectNewLaterData, 1), 1);
        s3 = swarmchart(x3, subjectNewLaterData.AnticipatedSwapping, 30, shape, 'MarkerEdgeColor', 'none','MarkerFaceColor', colors{3}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3);
%         s3.XJitter = 'rand';
        s3.XJitterWidth = 0.75;
	end

	    % Create boxplot
    boxplot(allData, groups, 'Colors', cell2mat(colors(:)), 'Width', 0.5, 'Symbol', '');
    set(findobj(gca,'type','line'),'linew',2);
    
    % Formatting the plot
    title('Anticipated Swapping Across Conditions', 'FontSize', 12);
    xlabel('Condition', 'FontSize', 12);
    ylabel('Anticipated Swapping', 'FontSize', 12);
    xticklabels({'New Early', 'Repeat', 'New Later'});
    xlim([0.5, 3.5]);
    yline(0, 'k--', 'HandleVisibility', 'off');
    
    ax = gca;  % Get current axis handle
    ax.FontSize = 12;  % Set font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    
    % Perform Welch's t-tests
    [~, p12, ~, stats12] = ttest2(newEarlyData.AnticipatedSwapping, repeatData.AnticipatedSwapping, 'Vartype', 'unequal');
    [~, p23, ~, stats23] = ttest2(repeatData.AnticipatedSwapping, newLaterData.AnticipatedSwapping, 'Vartype', 'unequal');
    [~, p13, ~, stats13] = ttest2(newEarlyData.AnticipatedSwapping, newLaterData.AnticipatedSwapping, 'Vartype', 'unequal');
    
    % Add significance indicators
    yMax = max(allData);
    yRange = range(allData);
    
    % Function to add significance bars
    function addSignificanceBar(x1, x2, y, p)
        plot([x1, x2], [y, y], 'k-', 'LineWidth', 1.5);
        if p < 0.001
            text(mean([x1, x2]), y, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        elseif p < 0.01
            text(mean([x1, x2]), y, '**', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        elseif p < 0.05
            text(mean([x1, x2]), y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        else
            text(mean([x1, x2]), y, 'ns', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        end
    end
    
    % Add significance bars
    addSignificanceBar(1, 2, yMax + 0.07*yRange, p12);
    addSignificanceBar(1, 3, yMax + 0.16*yRange, p13);
    addSignificanceBar(2, 3, yMax + 0.25*yRange, p23);
    
    ylim([min(allData) - 0.05*yRange, yMax + 0.35*yRange]);
    yticks([-1, -0.5, 0, 0.5, 1]);
	set(gcf,'renderer','Painters');
    % Add error bars for each condition
    conditions = {newEarlyData, repeatData, newLaterData};
    for i = 1:length(conditions)
        data = conditions{i}.AnticipatedSwapping;
        meanVal = nanmean(data);
        sem = nanstd(data) / sqrt(length(data));
        errorbar(i, meanVal, sem, 'k', 'LineWidth', 2, 'HandleVisibility', 'off', CapSize=10);
    end
    
    % Create legend for subjects
    subjectHandles = gobjects(1, length(subjects));
    for i = 1:length(subjects)
        subjectHandles(i) = plot(NaN, NaN, shapes{mod(i-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{i}(1)]);
    end
    
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
    
    hold off
    
    % Display statistical results
    fprintf('Welch''s t-test results:\n');
    fprintf('New Early vs Repeat: t(%.2f) = %.4f, p = %.4f\n', stats12.df, stats12.tstat, p12);
    fprintf('Repeat vs New Later: t(%.2f) = %.4f, p = %.4f\n', stats23.df, stats23.tstat, p23);
    fprintf('New Early vs New Later: t(%.2f) = %.4f, p = %.4f\n', stats13.df, stats13.tstat, p13);
end

function plotErrorInferenceInRepeat(block_table)
    % Define colors for conditions
%     colors = {[80/255, 126/255, 186/255],[209/255, 151/255, 35/255],  [179/255, 41/255, 39/255]}; % Yellow, Blue, Red
	colors = {[20, 52, 94] ./ 255,[213, 197, 131] ./ 255,  [213,174, 200]./255}; % 

    % Define shapes for subjects
    shapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
    
    % Get unique subjects
    subjects = unique(block_table.Subject);
    
    % Extract data for each condition
    repeatData = block_table(block_table.isRepetition == 1 & block_table.isNewRepetition == 0 & block_table.isSwapped == 1, :);
    newEarlyData = block_table(block_table.isRepetition == 0 & block_table.isSwapped == 1, :);
    newLaterData = block_table(block_table.isRepetition == 1 & block_table.isNewRepetition == 1 & block_table.isSwapped == 1, :);
    
    % Combine data for boxplot
    allData = [ newEarlyData.SwappedErrorCorrection;repeatData.SwappedErrorCorrection; newLaterData.SwappedErrorCorrection];
    groups = [ 1*ones(size(newEarlyData, 1), 1);2*ones(size(repeatData, 1), 1); 3*ones(size(newLaterData, 1), 1)];
    
    % Define canvas size
    canvas_size = [100, 100, 300, 300];
    
    % Create figure
    fig = figure('Position', canvas_size);
    hold on;
    
   
    
    % Add individual data points for each subject
    for i = 1:length(subjects)
        subject = subjects{i};
        shape = shapes{mod(i-1, length(shapes))+1};

		        % Repeat condition
        subjectRepeatData = repeatData(strcmp(repeatData.Subject, subject), :);
        x1 = 2 * ones(size(subjectRepeatData, 1), 1);
        s1 = swarmchart(x1, subjectRepeatData.SwappedErrorCorrection, 30, shape, 'MarkerFaceColor', colors{2}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3,  'MarkerEdgeColor', 'none');
%         s1.XJitter = 'rand';
        s1.XJitterWidth = 0.75;

		        % New Early condition
        subjectNewEarlyData = newEarlyData(strcmp(newEarlyData.Subject, subject), :);
        x2 = ones(size(subjectNewEarlyData, 1), 1);
        s2 = swarmchart(x2, subjectNewEarlyData.SwappedErrorCorrection, 30, shape,  'MarkerFaceColor', colors{1}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3,  'MarkerEdgeColor', 'none');
        s2.XJitter = 'rand';
        s2.XJitterWidth = 0.75;

        % New Later condition
        subjectNewLaterData = newLaterData(strcmp(newLaterData.Subject, subject), :);
        x3 = 3*ones(size(subjectNewLaterData, 1), 1);
        s3 = swarmchart(x3, subjectNewLaterData.SwappedErrorCorrection, 30,  shape, 'MarkerFaceColor', colors{3}, 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
% 		s3.XJitter = 'rand';
        s3.XJitterWidth = 0.75;
	end
    
	 % Create boxplot
    boxplot(allData, groups, 'Colors', cell2mat(colors(:)), 'Width', 0.5, 'Symbol', '');
    set(findobj(gca,'type','line'),'linew',2);

    % Formatting the plot
    title('Ordinal Object-identity Inference Across Conditions', 'FontSize', 12);
    xlabel('Condition', 'FontSize', 12);
    ylabel('Ordinal Object-identity Inference', 'FontSize', 12);
    xticklabels({ 'New Early', 'Repeat','New Later'});
    xlim([0.5, 3.5]);
    yline(0, 'k--', 'HandleVisibility', 'off');
    
    ax = gca;  % Get current axis handle
    ax.FontSize = 12;  % Set font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

	set(gcf,'renderer','Painters')

    
    % Perform Welch's t-tests
    [~, p12, ~, stats12] = ttest2(repeatData.SwappedErrorCorrection, newEarlyData.SwappedErrorCorrection, 'Vartype', 'unequal');
    [~, p23, ~, stats23] = ttest2(repeatData.SwappedErrorCorrection, newLaterData.SwappedErrorCorrection, 'Vartype', 'unequal');
    [~, p13, ~, stats13] = ttest2(newEarlyData.SwappedErrorCorrection, newLaterData.SwappedErrorCorrection, 'Vartype', 'unequal');
    
    % Add significance indicators
    yMax = max(allData);
    yRange = range(allData);
    
    % Function to add significance bars
    function addSignificanceBar(x1, x2, y, p)
        plot([x1, x2], [y, y], 'k-', 'LineWidth', 1.5);
        if p < 0.001
            text(mean([x1, x2]), y, '***', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        elseif p < 0.01
            text(mean([x1, x2]), y, '**', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        elseif p < 0.05
            text(mean([x1, x2]), y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        else
            text(mean([x1, x2]), y, 'ns', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        end
    end
    
    % Add significance bars
    addSignificanceBar(1, 2, yMax + 0.07*yRange, p12);
    addSignificanceBar(1, 3, yMax + 0.16*yRange, p13);
    addSignificanceBar(2, 3, yMax + 0.25*yRange, p23);
    
    ylim([min(allData) - 0.05*yRange, yMax + 0.35*yRange]);
    yticks(-1:0.5:1);

	    % Add error bars for each condition
    conditions = {newEarlyData,repeatData, newLaterData};
    for i = 1:length(conditions)
        data = conditions{i}.SwappedErrorCorrection;
        meanVal = nanmean(data);
        sem = nanstd(data) / sqrt(length(data));
        errorbar(i, meanVal, sem, 'k', 'LineWidth', 2, HandleVisibility='off', CapSize=10);
    end
    

    
    % Create legend for subjects
    subjectHandles = gobjects(1, length(subjects));
    for i = 1:length(subjects)
        subjectHandles(i) = plot(NaN, NaN, shapes{mod(i-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{i}(1)]);
    end
    
    legend_obj = legend(subjectHandles, 'Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
	
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
    
	hold off
    % Display statistical results
    fprintf('Welch''s t-test results:\n');
    fprintf('Repeat vs New Early: t(%.2f) = %.4f, p = %.4f\n', stats12.df, stats12.tstat, p12);
    fprintf('Repeat vs New Later: t(%.2f) = %.4f, p = %.4f\n', stats13.df, stats13.tstat, p13);
    fprintf('New Early vs New Later: t(%.2f) = %.4f, p = %.4f\n', stats23.df, stats23.tstat, p23);
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






function plotCompletionRateVsBlocks(block_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Round completion rate to 2 decimal places
    block_table.CompletionRate = round(block_table.BlockLevelAccuracy, 2);

    % Get unique subjects
    subjects = unique(block_table.Subject);

    % Create the scatter plot
    fig = figure;
    fig.Position = canvas_size;
    hold on;

    % Scatter plot for each subject
    markers = {'o', 's', 'd', '^'}; % Different markers for each subject
    colors = [hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A')];
    all_completion_rates = [];
    all_normalized_counts = [];
    subject_medians = zeros(length(subjects), 1);

    for i = 1:length(subjects)
        subject_data = block_table(strcmp(block_table.Subject, subjects{i}), :);

        % Create bins of size 0.05
        bins = 0:0.05:1;
        [count, edges] = histcounts(subject_data.CompletionRate, bins);

        % Calculate bin centers
        bin_centers = edges(2:end);

        % Normalize the count
        normalized_count = count / sum(count);

        % Scatter plot
        scatter(bin_centers, normalized_count, 50, colors(i,:), markers{i}, 'filled');

        % Fit Gaussian curve
        [fitresult, ~] = fitHalfNormal(bin_centers, normalized_count);

        % Plot fitted Gaussian curve
        x_fit = linspace(0, 1, 100);
        y_fit = feval(fitresult, x_fit);
        plot(x_fit, y_fit, 'Color', colors(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');

        all_completion_rates = [all_completion_rates; subject_data.CompletionRate];
        all_normalized_counts = [all_normalized_counts; normalized_count'];
        subject_medians(i) = median(subject_data.CompletionRate);
    end

    % Add vertical line for median of all subjects
    median_completion_rate = median(all_completion_rates);
    xline(median_completion_rate, 'k-', 'LineWidth', 2);
    y_max = max(all_normalized_counts) + 0.05;

    for i = 1:length(subjects)
        plot(subject_medians(i), y_max, 'v', 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
		fprintf('%s, %.2f\n', subjects{i}, subject_medians(i));
    end

    % Customize the plot
    subject_labels = cellfun(@(x) ['Subject ' x(1)], subjects, 'UniformOutput', false);
    xlabel('P(Completion)', 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 12);
    title('Completion Rate vs. Proportion of Blocks', 'FontSize', 14);

    % Create legend and move it outside
    legend_obj = legend([subject_labels; 'Median Overall'], 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';

    % Adjust axis limits
    xlim([0 1]);
    ylim([0 max(all_normalized_counts)*1.3]);

    % Improve plot appearance
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
	set(gcf,'renderer','Painters')
    % Adjust figure size to accommodate legend
	legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function [fitresult, gof] = fitHalfNormal(x, y)
    ft = fittype('scale * (sqrt(2/pi)/sigma) * exp(-(x-1)^2 / (2*sigma^2))', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.StartPoint = [0.1, max(y)]; % [sigma, scale]
    [fitresult, gof] = fit(x', y', ft, opts);
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



function plotErrorDecreasingRate(touch_table)
    % Set canvas size
    canvas_size = [100, 100, 300, 250];

    % Merge all subjects for analysis
    subjects = unique(touch_table.Subject);
    num_subjects = length(subjects);

    % Filter for non-repetitive trials
    touch_table = touch_table(touch_table.isRepetition == 0, :);

    % Define the error type to analyze (Distractor Errors)
    error_type = 'DistractorError';
    
    % Initialize storage for error rates
    all_error_rates = [];

    % Storage for individual subject error rates
    subject_error_rates = cell(num_subjects, 15);

    % Process data for each session
    sessions = unique(touch_table.SessionNumber);
    
    for i = 1:length(sessions)
        session = sessions(i);
        session_rows = touch_table(touch_table.SessionNumber == session, :);
        blocks = unique(session_rows.BlockNumber);
        
        for block = blocks'
            block_rows = session_rows(session_rows.BlockNumber == block, :);

            for trial = 1:15  % Assuming 15 trials per block
                % Process each subject's touches for the current trial
                for s = 1:num_subjects
                    subject_rows = block_rows(strcmp(block_rows.Subject, subjects{s}) & block_rows.TrialNumberInBlock == trial, :);
                    total_touches = height(subject_rows);

                    if total_touches == 0
                        continue  % Skip if there were no touches for this trial
                    end

                    % Calculate Distractor Errors (TouchObjectCorrectOrdinalPosition == 0)
                    distractor_touches = sum(subject_rows.TouchObjectCorrectOrdinalPosition == 0);
                    error_rate = distractor_touches / total_touches;  % Proportion of distractor errors

                    % Save error rate as {subject, trial}
                    subject_error_rates{s, trial} = [subject_error_rates{s, trial}; error_rate];
                end
            end
        end
    end

    
    % Plotting
    fig = figure;
    set(gcf,'renderer','Painters');  % Set renderer for high-quality vector graphics

    hold on;
    trials = 1:15;
    colors = [39, 93, 44;
			58, 112, 175;
			233, 166, 64;
			199, 54, 55;
			107, 37, 110;
			247, 211, 76]/255; % Color matching for subjects
    markers = {'o', 's', 'd', 'x'};  % Add more markers if needed
    
    % Plot individual subject lines
    % Plot individual subject lines
    for s = 1:num_subjects
        subject_mean = zeros(1, 15);  % Initialize mean error rates for 15 trials

        for trial = 1:15
            % Calculate the mean error rate for this subject in this trial
            if ~isempty(subject_error_rates{s, trial})
                subject_mean(trial) = mean(subject_error_rates{s, trial}, 'omitnan');
            else
                subject_mean(trial) = NaN;  % Handle missing data
            end
        end

        % Plot the subject's line
        plot(1:15, subject_mean, '-', 'Color', [colors(s,:), 0.7], 'LineWidth', 1, ...
            'DisplayName', ['Subject ', subjects{s}]);
    end

	% Calculate mean error rates across all subjects
    mean_rate = zeros(1, 15);  % Initialize mean error rate for 15 trials
	sem = zeros(1, 15);
	ci = zeros(1, 15);

    for trial = 1:15
        all_subject_errors = [];
        for s = 1:num_subjects
            if ~isempty(subject_error_rates{s, trial})
                all_subject_errors = [all_subject_errors; subject_error_rates{s, trial}];  % Collect errors from all subjects
            end
        end
        mean_rate(trial) = mean(all_subject_errors, 'omitnan');  % Calculate the mean error rate for this trial
		sem(trial) = std(all_subject_errors, 0, 1, 'omitnan') / sqrt(size(all_subject_errors, 1));
		ci(trial) = 1.96 * sem(trial);
	end

	 % Plot averaged error rate (thicker line)
	trials = 1:15;
    plot(trials, mean_rate, '-', 'Color', colors(5,:), 'LineWidth', 2, 'DisplayName', 'Average');
    fill([trials, fliplr(trials)], [mean_rate - ci, fliplr(mean_rate + ci)], colors(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Customize plot
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Proportion of Distractor Errors', 'FontSize', 12);
    title('Proportion of Distractor Errors Across Trials', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');  % Apply plotting rules
	xticks([1, 5, 10, 15]);
	xlim([0.5, 15.5]);
    % Set canvas size
    fig.Position = canvas_size;

    % Create custom legend with subject lines
    legendHandles = gobjects(1, num_subjects + 1);  % Including the average line
    for s = 1:num_subjects
        % Use NaN to avoid adding a new plot, just for the legend
        legendHandles(s) = plot(NaN, NaN, 'Color', colors(s,:), ...
                                'LineWidth', 1);
    end
    
    % Add a handle for the averaged line
    legendHandles(end) = plot(NaN, NaN, '-', 'Color', colors(5,:), 'LineWidth', 1);
    
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


function plotDistractorErrorProportionsByState(touch_table)
    % Constants
    maxState = 5;

    % Initialize storage for error counts and total observations
    distractorErrorsSwapped = zeros(maxState, 1);
    distractorErrorsNotSwapped = zeros(maxState, 1);
    trialsReachedStateSwapped = zeros(maxState, 1);
    trialsReachedStateNotSwapped = zeros(maxState, 1);

    % Unique identifiers for subjects, sessions, blocks
    subjects = unique(touch_table.Subject);
    for s = 1:length(subjects)
        subjectData = touch_table(strcmp(touch_table.Subject, subjects(s)), :);
        sessions = unique(subjectData.SessionNumber);
        for sess = 1:length(sessions)
            sessionData = subjectData(subjectData.SessionNumber == sessions(sess), :);
            blocks = unique(sessionData.BlockNumber);
            for b = 1:length(blocks)
                blockData = sessionData(sessionData.BlockNumber == blocks(b), :);
                trials = unique(blockData.TrialNumberInBlock);
                for t = 1:length(trials)
                    trialData = blockData(blockData.TrialNumberInBlock == trials(t), :);
                    if trialData.isRepetition(1) == 1
                        continue; % Skip repetition trials
                    end

                    isSwapped = trialData.isSwapped(1) == 1;

                    % Initialize array to track which positions have had errors in this trial
                    errorPositions = false(1, maxState);
                    maxReachedState = 0;

                    % Process each touch
                    for idx = 1:height(trialData)
                        currentState = trialData.CurrentState(idx);
                        errorState = 1;
                        if idx > 1
                            errorState = trialData.CurrentState(idx-1) + 1;
                            if currentState == 0 && trialData.CurrentState(idx-1) == 0
                                errorState = 1;
                            end
                        end
                        errorState = min(errorState, maxState);
                        
                        maxReachedState = max(maxReachedState, errorState);
      
                        % Check for distractor error
                        if trialData.TouchObjectCorrectOrdinalPosition(idx) == 0 && ~errorPositions(errorState)
                            if isSwapped
                                distractorErrorsSwapped(errorState) = distractorErrorsSwapped(errorState) + 1;
                            else
                                distractorErrorsNotSwapped(errorState) = distractorErrorsNotSwapped(errorState) + 1;
                            end
                            errorPositions(errorState) = true; % Mark this position as having an error
                        end
					end

					if maxReachedState < 5
						maxReachedState = maxReachedState + 1;
					end

                    % Update counts of trials that reached each state
                    if isSwapped
                        trialsReachedStateSwapped(1:maxReachedState) = trialsReachedStateSwapped(1:maxReachedState) + 1;
                    else
                        trialsReachedStateNotSwapped(1:maxReachedState) = trialsReachedStateNotSwapped(1:maxReachedState) + 1;
                    end
                end
            end
        end
    end

    % Calculate proportions and standard errors
    proportionDistrErrorsSwapped = distractorErrorsSwapped ./ trialsReachedStateSwapped;
    proportionDistrErrorsNotSwapped = distractorErrorsNotSwapped ./ trialsReachedStateNotSwapped;
    
    seSwapped = sqrt(proportionDistrErrorsSwapped .* (1 - proportionDistrErrorsSwapped) ./ trialsReachedStateSwapped);
    seNotSwapped = sqrt(proportionDistrErrorsNotSwapped .* (1 - proportionDistrErrorsNotSwapped) ./ trialsReachedStateNotSwapped);

    % Plotting
    figure; hold on;
    states = 1:maxState;

    % Plot for Swapped
    errorbar(states, proportionDistrErrorsSwapped, seSwapped, 'o-', 'Color', '#F2BB6B', 'LineWidth', 2, 'MarkerFaceColor', '#F2BB6B', 'DisplayName', 'Swapped');

    % Plot for Not Swapped
    errorbar(states, proportionDistrErrorsNotSwapped, seNotSwapped, 'o-', 'Color', '#276C9E', 'LineWidth', 2, 'MarkerFaceColor', '#276C9E', 'DisplayName', 'Not Swapped');

    % Formatting the plot
    title('Proportion of Trials Containing Distractor Error at Ordinal Position X', 'FontSize', 16);
    xlabel('Position in Sequence', 'FontSize', 14);
    ylabel('Proportion of Trials', 'FontSize', 14);
    xlim([0.5, maxState + 0.5]);
    ylim([0, 1]); % Set y-axis limit from 0 to 1 for proportion
    xticks(1:maxState);
    legend('show', 'Location', 'best', 'FontSize', 12);
    ax = gca; % Get current axis handle
    ax.FontSize = 12; % Increase font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    hold off;
end



function plotTrialCompletionRatesRepeat(touch_table)
    % Define trial range and accuracy threshold
    maxTrialNumber = 15;
    accuracy_threshold = 0.7;

    % Get unique subjects, sessions, and blocks
    subjects = unique(touch_table.Subject);
    numSubjects = length(subjects);

    % Initialize arrays to store the data
    completedTrialsRepeat = zeros(maxTrialNumber, 1);
    totalTrialsRepeat = zeros(maxTrialNumber, 1);
    completedTrialsNewEarly = zeros(maxTrialNumber, 1);
    totalTrialsNewEarly = zeros(maxTrialNumber, 1);
    completedTrialsNewLater = zeros(maxTrialNumber, 1);
    totalTrialsNewLater = zeros(maxTrialNumber, 1);

    % Initialize arrays for individual subject data
    subjectDataRepeat = zeros(numSubjects, maxTrialNumber);
    subjectDataNewEarly = zeros(numSubjects, maxTrialNumber);
    subjectDataNewLater = zeros(numSubjects, maxTrialNumber);
    subjectTrialsRepeat = zeros(numSubjects, maxTrialNumber);
    subjectTrialsNewEarly = zeros(numSubjects, maxTrialNumber);
    subjectTrialsNewLater = zeros(numSubjects, maxTrialNumber);

    % Initialize arrays for learning speeds
    learningSpeedsRepeat = [];
    learningSpeedsNewEarly = [];
    learningSpeedsNewLater = [];

    % Loop through subjects, sessions, blocks, and trials
    for s = 1:numSubjects
        subject = subjects{s};
        subjectData = touch_table(strcmp(touch_table.Subject, subject), :);
        sessions = unique(subjectData.SessionNumber);
        
        for sess = sessions'
            sessionData = subjectData(subjectData.SessionNumber == sess, :);
            blocks = unique(sessionData.BlockNumber);
            
            sessionCompletedTrialsRepeat = [];
            sessionCompletedTrialsNewEarly = [];
            sessionCompletedTrialsNewLater = [];
            
            for blk = blocks'
                blockData = sessionData(sessionData.BlockNumber == blk, :);
                
                % Check block type
                isRepetition = blockData.isRepetition(1);
                isNewRepetition = blockData.isNewRepetition(1);
                
                completedTrials = zeros(1, maxTrialNumber);
                for trial = 1:min(height(blockData), maxTrialNumber)
                    trialData = blockData(blockData.TrialNumberInBlock == trial, :);
                    
                    if isempty(trialData)
                        continue;
                    end
                    
                    isCompleted = any(trialData.isRewarded == 1);
                    completedTrials(trial) = isCompleted;
                    
                    if isRepetition == 1 && isNewRepetition == 0 % Repeat trials
                        totalTrialsRepeat(trial) = totalTrialsRepeat(trial) + 1;
                        subjectTrialsRepeat(s, trial) = subjectTrialsRepeat(s, trial) + 1;
                        if isCompleted
                            completedTrialsRepeat(trial) = completedTrialsRepeat(trial) + 1;
                            subjectDataRepeat(s, trial) = subjectDataRepeat(s, trial) + 1;
                        end
                    elseif isRepetition == 0 % New Early trials
                        totalTrialsNewEarly(trial) = totalTrialsNewEarly(trial) + 1;
                        subjectTrialsNewEarly(s, trial) = subjectTrialsNewEarly(s, trial) + 1;
                        if isCompleted
                            completedTrialsNewEarly(trial) = completedTrialsNewEarly(trial) + 1;
                            subjectDataNewEarly(s, trial) = subjectDataNewEarly(s, trial) + 1;
                        end
                    elseif isRepetition == 1 && isNewRepetition == 1 % New Later trials
                        totalTrialsNewLater(trial) = totalTrialsNewLater(trial) + 1;
                        subjectTrialsNewLater(s, trial) = subjectTrialsNewLater(s, trial) + 1;
                        if isCompleted
                            completedTrialsNewLater(trial) = completedTrialsNewLater(trial) + 1;
                            subjectDataNewLater(s, trial) = subjectDataNewLater(s, trial) + 1;
                        end
                    end
                end
                
                % Add completed trials to the appropriate session array
                if isRepetition == 1 && isNewRepetition == 0
                    sessionCompletedTrialsRepeat = [sessionCompletedTrialsRepeat; completedTrials];
                elseif isRepetition == 0
                    sessionCompletedTrialsNewEarly = [sessionCompletedTrialsNewEarly; completedTrials];
                elseif isRepetition == 1 && isNewRepetition == 1
                    sessionCompletedTrialsNewLater = [sessionCompletedTrialsNewLater; completedTrials];
                end
            end
            
            % Calculate learning speed for each condition in the session
            learningSpeedsRepeat = [learningSpeedsRepeat, calculateLearningSpeed(sessionCompletedTrialsRepeat, accuracy_threshold)];
            learningSpeedsNewEarly = [learningSpeedsNewEarly, calculateLearningSpeed(sessionCompletedTrialsNewEarly, accuracy_threshold)];
            learningSpeedsNewLater = [learningSpeedsNewLater, calculateLearningSpeed(sessionCompletedTrialsNewLater, accuracy_threshold)];
        end
    end

    % Calculate proportions and standard errors
    proportionRepeat = completedTrialsRepeat ./ totalTrialsRepeat;
    proportionNewEarly = completedTrialsNewEarly ./ totalTrialsNewEarly;
    proportionNewLater = completedTrialsNewLater ./ totalTrialsNewLater;
    seRepeat = sqrt(proportionRepeat .* (1 - proportionRepeat) ./ totalTrialsRepeat);
    seNewEarly = sqrt(proportionNewEarly .* (1 - proportionNewEarly) ./ totalTrialsNewEarly);
    seNewLater = sqrt(proportionNewLater .* (1 - proportionNewLater) ./ totalTrialsNewLater);

    % Sigmoid fit function
    sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', 'independent', 'x', 'coefficients', {'a', 'b', 'c'});
    options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, 0.1, 7.5]);

    % Fit models
    fitRepeat = fit((1:maxTrialNumber)', proportionRepeat, sigmoidModel, options);
    fitNewEarly = fit((1:maxTrialNumber)', proportionNewEarly, sigmoidModel, options);
    fitNewLater = fit((1:maxTrialNumber)', proportionNewLater, sigmoidModel, options);

    % Colors
	colors = {[20, 52, 94] ./ 255,[213, 197, 131] ./ 255,  [213,174, 200]./255};
%     colors = {'#D19723'; '#507EBA'; '#B32927'; '#728A50'}; % Updated colors
    
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Plotting Figure
    fig = figure('Position', canvas_size);
    hold on;
    trials = 1:maxTrialNumber;

    % Plot with shaded area for New Early trials
    fill([trials, fliplr(trials)], [proportionNewEarly' + 1.96 * seNewEarly', fliplr(proportionNewEarly' - 1.96 * seNewEarly')], colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(trials, feval(fitNewEarly, trials), 'Color', colors{1}, 'LineWidth', 2, 'DisplayName', 'New Early');
    plot((1:maxTrialNumber), proportionNewEarly, 'o', 'Color', colors{1}, 'MarkerFaceColor', colors{1}, 'LineStyle', 'none', 'HandleVisibility', 'off');

    % Plot with shaded area for Repeat trials
    fill([trials, fliplr(trials)], [proportionRepeat' + 1.96 * seRepeat', fliplr(proportionRepeat' - 1.96 * seRepeat')], colors{2}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(trials, feval(fitRepeat, trials), 'Color', colors{2}, 'LineWidth', 2, 'DisplayName', 'Repeat');
    plot((1:maxTrialNumber), proportionRepeat, 'o', 'Color', colors{2}, 'MarkerFaceColor', colors{2}, 'LineStyle', 'none', 'HandleVisibility', 'off');

    % Plot with shaded area for New Later trials
    fill([trials, fliplr(trials)], [proportionNewLater' + 1.96 * seNewLater', fliplr(proportionNewLater' - 1.96 * seNewLater')], colors{3}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(trials, feval(fitNewLater, trials), 'Color', colors{3}, 'LineWidth', 2, 'DisplayName', 'New Later');
    plot((1:maxTrialNumber), proportionNewLater, 'o', 'Color', colors{3}, 'MarkerFaceColor', colors{3}, 'LineStyle', 'none', 'HandleVisibility', 'off');

    % Add triangles for mean learning speeds with error bars
    plotLearningSpeed(learningSpeedsNewEarly, 0.94, colors{1},'New Early');
    plotLearningSpeed(learningSpeedsRepeat, 0.95, colors{2}, 'Repeat');
    plotLearningSpeed(learningSpeedsNewLater, 0.96, colors{3}, 'New Late');

    % Formatting the plot
    title('Trial Completion Rates: Repeat vs. New Early vs. New Later', 'FontSize', 12);
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Completion Rate', 'FontSize', 12);
    xticks(1:2:15);
    ylim([0.5, 1]);
    xlim([0.5, 15.5]);
    yline(0.7, '--', 'HandleVisibility', 'off');

    % Apply plot rules
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
	set(gcf,'renderer','Painters')

    % Add legend
    legend_obj = legend('Location', 'eastoutside', 'FontSize', 12);
    set(legend_obj, 'Box', 'off');
    
    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end

function learning_speed = calculateLearningSpeed(session_completed_trials, accuracy_threshold)
    if isempty(session_completed_trials)
        learning_speed = NaN;
        return;
    end
    
    mean_completion_rate = nanmean(session_completed_trials, 1);
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
end

function learning_speed = calculateLearningSpeed2(session_completed_trials, accuracy_threshold)
if isempty(session_completed_trials)
	learning_speed = NaN;
	return;
end
mean_completion_rate = nanmean(session_completed_trials, 1);
mean_completion_rate(end) = mean(mean_completion_rate(end-3:end));
mean_completion_rate(end-1) = mean(mean_completion_rate(end-3:end-1));
half_sigmoid = @(b,x) b(1)./(1 + exp(-b(2)*(x)));
x_data = 1:length(mean_completion_rate);
initial_guess = [1, 0.1];
[fitted_params, ~] = nlinfit(x_data, mean_completion_rate, half_sigmoid, initial_guess);
a = fitted_params(1);
b = fitted_params(2);
learning_speed = NaN;
if a >= accuracy_threshold
	% Find where sigmoid = 0.8 (accuracy_threshold)
	crossing_point = -log((a/accuracy_threshold) - 1)/b;
	if ~isnan(crossing_point) && crossing_point >= 0  && crossing_point < 16
		learning_speed = crossing_point;
	end
end

end

function plotLearningSpeed(learning_speeds, x_offset, color, name)
    valid_speeds = learning_speeds(~isnan(learning_speeds));
    mean_speed = nanmean(valid_speeds);
    ci = 1.96 * nanstd(valid_speeds) / sqrt(length(valid_speeds));
    fprintf('%s: %.2f ± %.2f\n',name, mean_speed, ci);
	fprintf('valid sessions: %d, all sessions %d \n', length(~isnan(learning_speeds)), length(learning_speeds));
    plot(mean_speed, x_offset, 'v', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 8, 'HandleVisibility', 'off');
    errorbar( mean_speed, x_offset, ci, 'horizontal', 'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
end



function plotDistractorErrorRates(touch_table)
    % Extract unique subjects
    subjects = unique(touch_table.Subject);
    maxTrialNumber = 15; % Define the maximum trial number for plotting

    % Initialize arrays to store counts of trials with and without distractor errors
    trialsWithErrorSwapped = zeros(maxTrialNumber, 1);
    trialsWithErrorNotSwapped = zeros(maxTrialNumber, 1);
    totalTrialsPerNumberSwapped = zeros(maxTrialNumber, 1);
    totalTrialsPerNumberNotSwapped = zeros(maxTrialNumber, 1);

    % Iterate through each subject
    for s = 1:numel(subjects)
        subjectData = touch_table(strcmp(touch_table.Subject, subjects{s}), :);
        % Extract sessions for current subject
        sessions = unique(subjectData.SessionNumber);

        % Iterate through each session
        for sess = 1:numel(sessions)
            sessionData = subjectData(subjectData.SessionNumber == sessions(sess), :);
            % Extract blocks for current session
            blocks = unique(sessionData.BlockNumber);

            % Iterate through each block
            for b = 1:numel(blocks)
                blockData = sessionData(sessionData.BlockNumber == blocks(b), :);
                % Extract trials for current block
                trials = unique(blockData.TrialNumberInBlock);

                % Iterate through each trial
                for t = 1:numel(trials)
                    trialData = blockData(blockData.TrialNumberInBlock == trials(t), :);

                    if any(trialData.isRepetition == 1) % Skip if trial is a repetition
                        continue;
                    end

                    trialNumber = trials(t);
                    if trialData.isSwapped(1) % Check if swapped
                        totalTrialsPerNumberSwapped(trialNumber) = totalTrialsPerNumberSwapped(trialNumber) + 1;
                        if any(trialData.TouchObjectCorrectOrdinalPosition == 0)
                            trialsWithErrorSwapped(trialNumber) = trialsWithErrorSwapped(trialNumber) + 1;
                        end
                    else
                        totalTrialsPerNumberNotSwapped(trialNumber) = totalTrialsPerNumberNotSwapped(trialNumber) + 1;
                        if any(trialData.TouchObjectCorrectOrdinalPosition == 0)
                            trialsWithErrorNotSwapped(trialNumber) = trialsWithErrorNotSwapped(trialNumber) + 1;
                        end
                    end
                end
            end
        end
    end

    % Calculate the proportion of trials with at least one distractor error
    distractorErrorProportionsSwapped = trialsWithErrorSwapped ./ totalTrialsPerNumberSwapped;
    distractorErrorProportionsNotSwapped = trialsWithErrorNotSwapped ./ totalTrialsPerNumberNotSwapped;
    distractorErrorSEsSwapped = sqrt(distractorErrorProportionsSwapped .* (1 - distractorErrorProportionsSwapped) ./ totalTrialsPerNumberSwapped);
    distractorErrorSEsNotSwapped = sqrt(distractorErrorProportionsNotSwapped .* (1 - distractorErrorProportionsNotSwapped) ./ totalTrialsPerNumberNotSwapped);

    % Plotting
    figure; hold on;
    trials = 1:maxTrialNumber;
    colors = {'#DE7833', '#329845'}; % Colors for swapped and not swapped

    % Sigmoid model fitting
    sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', 'independent', 'x', 'coefficients', {'a', 'b', 'c'});
    initialGuess = [1, 0.1, maxTrialNumber / 2];

    % Plot and fit for Swapped
    [fitResultSwapped, ~] = fit(trials', distractorErrorProportionsSwapped, sigmoidModel, 'StartPoint', initialGuess);
    plot(trials, feval(fitResultSwapped, trials), '-', 'Color', colors{1}, 'LineWidth', 2, 'DisplayName', 'Swapped');
    fill([trials, fliplr(trials)], [distractorErrorProportionsSwapped' + 1.96 * distractorErrorSEsSwapped', fliplr(distractorErrorProportionsSwapped' - 1.96 * distractorErrorSEsSwapped')],hex2rgb(colors{1}), 'FaceAlpha', 0.1, 'EdgeColor', 'none', HandleVisibility='off');

    % Plot and fit for Not Swapped
	
    [fitResultNotSwapped, ~] = fit(trials', distractorErrorProportionsNotSwapped, sigmoidModel, 'StartPoint', initialGuess);
    plot(trials, feval(fitResultNotSwapped, trials), '-', 'Color', colors{2}, 'LineWidth', 2, 'DisplayName', 'Not Swapped');
    fill([trials, fliplr(trials)], [distractorErrorProportionsNotSwapped' + 1.96 * distractorErrorSEsNotSwapped', fliplr(distractorErrorProportionsNotSwapped' - 1.96 * distractorErrorSEsNotSwapped')], hex2rgb(colors{2}), 'FaceAlpha', 0.1, 'EdgeColor', 'none', HandleVisibility='off');
	
	plot(trials, distractorErrorProportionsSwapped, 'o', 'Color', colors{1}, 'MarkerFaceColor', colors{1}, HandleVisibility='off');
    plot(trials, distractorErrorProportionsNotSwapped, 'o', 'Color', colors{2}, 'MarkerFaceColor', colors{2}, HandleVisibility='off');

	% Formatting the plot
    title('Proportion of Trials with Distractor Errors by Trial Number', 'FontSize', 12);
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Proportion of Errors', 'FontSize', 12);
    xlim([0.5, maxTrialNumber + 0.5]);
    xticks(1:maxTrialNumber);
    legend('Location', 'northeast', 'FontSize', 12);    ax = gca;  % Get current axis handle
    ax.FontSize = 14;  % Increase font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
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

function calculateSessionLearningSpeed(touch_table, block_table, accuracy_threshold)
    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    
    % Initialize containers for storing data
    subject_session_learning_speeds = containers.Map;
    all_completed_trials = [];
    all_subject_ids = {};
    subject_session_counts = zeros(1, numel(subjects));
    subject_block_counts = zeros(1, numel(subjects));
    all_learning_speeds = [];

    % Loop through each subject
    for subj = 1:numel(subjects)
        subject_data = touch_table(strcmp(touch_table.Subject, subjects{subj}), :);
        sessions = unique(subject_data.SessionNumber);
        subject_session_counts(subj) = numel(sessions);
        
        session_learning_speeds = [];
        block_count = 0;

        % Loop through each session
        for sess = 1:numel(sessions)
            session_data = subject_data(subject_data.SessionNumber == sessions(sess), :);
            blocks = unique(session_data.BlockNumber);
            
            session_completed_trials = [];

            % Loop through each block
            for blk = 1:numel(blocks)
                block_data = session_data(session_data.BlockNumber == blocks(blk), :);
                % Ensure the block is not swapped and not a repetition
                if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                    block_count = block_count + 1;
                    unique_trials = unique(block_data.TrialNumberInBlock);
                    numTrials = numel(unique_trials);
                    completed_trials = zeros(1, numTrials);
                    for trial = 1:numTrials
                        trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                        completed_trials(trial) = any(trial_data.isRewarded == 1);
                    end
                    
                    % Pad completed trials to 15 if necessary
                    padded_completed_trials = nan(1, 15);
                    padded_completed_trials(1:length(completed_trials)) = completed_trials;
                    session_completed_trials = [session_completed_trials; padded_completed_trials];
                end
            end
            
            % Calculate learning speed for the session
            if ~isempty(session_completed_trials)
                mean_completion_rate = nanmean(session_completed_trials, 1);
                learning_speed = NaN;

% 				half_sigmoid = @(b,x) b(1)./(1 + exp(-b(2)*(x) ));
% 				x_data = 1:length(mean_completion_rate);
% 				initial_guess = [1, 0.1, 7.5];
% 				[fitted_params, ~] = nlinfit(x_data, mean_completion_rate, half_sigmoid, initial_guess);
%             	a = fitted_params(1);
%             	b = fitted_params(2);
% 				
% 
% 				if a > accuracy_threshold
%                 % Find where sigmoid = 0.8 (accuracy_threshold)
%                 crossing_point = -log((a/accuracy_threshold) - 1)/b;
%                 	if ~isnan(crossing_point) && crossing_point >= 0  && crossing_point < 16
%                     	session_learning_speeds = [session_learning_speeds, crossing_point];
%                     	all_learning_speeds = [all_learning_speeds, crossing_point];
%                 	end
% 				end

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
                    session_learning_speeds = [session_learning_speeds, learning_speed];
                    all_learning_speeds = [all_learning_speeds, learning_speed];
				end
            end
            
            % Accumulate all completed trials for overall plot
% 			if ~isnan(crossing_point) && crossing_point >= 0  && crossing_point < 16
			if ~isnan(learning_speed)
    		all_completed_trials = [all_completed_trials; session_completed_trials];
    		all_subject_ids = [all_subject_ids; repmat(subjects(subj), size(session_completed_trials, 1), 1)];
		
			end
		end
        
        % Store session learning speeds for the subject
        subject_session_learning_speeds(subjects{subj}) = session_learning_speeds;
        subject_block_counts(subj) = block_count;
    end

    % Calculate average learning speed and SEM for each subject
    subject_avg_speeds = zeros(1, numel(subjects));
    subject_sem_speeds = zeros(1, numel(subjects));
    for subj = 1:numel(subjects)
        speeds = subject_session_learning_speeds(subjects{subj});
        subject_avg_speeds(subj) = mean(speeds);
        subject_sem_speeds(subj) = std(speeds) / sqrt(length(speeds));
        
        % Print out subject information
        fprintf('Subject %s:\n', subjects{subj});
        fprintf('  Average Learning Speed = %.2f ± %.2f\n', subject_avg_speeds(subj), subject_sem_speeds(subj));
        fprintf('  Number of Sessions: %d\n', subject_session_counts(subj));
        fprintf('  Number of Blocks: %d\n', subject_block_counts(subj));
        fprintf('  Number of Sessions with Valid Learning Speed: %d\n\n', length(speeds));
    end

    % Calculate and print overall average learning speed and SEM
    overall_avg_speed = mean(all_learning_speeds);
    overall_sem_speed = std(all_learning_speeds) / sqrt(length(all_learning_speeds));
    fprintf('Overall Average Learning Speed (across all subjects and sessions):\n');
    fprintf('  %.2f ± %.2f\n', overall_avg_speed, overall_sem_speed);
    fprintf('Total Number of Sessions with Valid Learning Speed: %d\n\n', length(all_learning_speeds));

    % Plotting
    canvas_size = [100, 100, 300, 250];
    fig = figure('Position', canvas_size);
    hold on;
	
    % Define colors for each subject
    colors = [hex2rgb('#4A7298'); hex2rgb('#F3C846'); hex2rgb('#C83E4D'); hex2rgb('#4E937A')];

    trials = 1:15;

    % Plot average completion rate
    meanCompletionRateAll = nanmean(all_completed_trials, 1);
    plot(trials, meanCompletionRateAll, '.', 'MarkerSize', 12, 'Color', '#EDAE92', 'LineWidth', 3, 'DisplayName', 'All Subjects', 'HandleVisibility', 'off');

    % Add shaded error bars for the average data
    sem_all = std(all_completed_trials, [], 1, 'omitnan') / sqrt(size(all_completed_trials, 1));
    ci_all = 1.96 * sem_all; % 95% CI
    lowerBoundAll = meanCompletionRateAll - ci_all;
    upperBoundAll = meanCompletionRateAll + ci_all;
    fill([trials fliplr(trials)], [lowerBoundAll fliplr(upperBoundAll)], hex2rgb('#EDAE92'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Plot each subject's completion rates and average learning speed marker
    for subj = 1:numel(subjects)
        subj_trials = all_completed_trials(strcmp(all_subject_ids, subjects{subj}), :);
        meanCompletionRateSubj = nanmean(subj_trials, 1);
        
        % Plot subject's completion rate with alpha 0.5
        plot(trials, meanCompletionRateSubj, 'Color', [colors(subj,:), 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        
        % Plot average learning speed with error bar
        errorbar(subject_avg_speeds(subj), 0.95 + 0.015 * subj, subject_sem_speeds(subj), 'horizontal', ...
                 'Color', colors(subj,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(subject_avg_speeds(subj), 0.95+ 0.015 * subj, 'v', 'Color', colors(subj,:), ...
             'MarkerFaceColor', colors(subj,:), ...
             'MarkerSize', 8, 'DisplayName', ['Subject ' subjects{subj}(1)]);
    end

    % Plot overall average learning speed
    errorbar(overall_avg_speed, 1.1, overall_sem_speed, 'horizontal', ...
             'Color', hex2rgb('#EDAE92'), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(overall_avg_speed, 1.1, 'v', 'Color', hex2rgb('#EDAE92'), ...
         'MarkerFaceColor', hex2rgb('#EDAE92'), ...
         'MarkerSize', 8, 'DisplayName', 'Overall Average');

    % Fit sigmoid model to all data
    sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', ...
                           'independent', 'x', ...
                           'coefficients', {'a', 'b', 'c'});
    initialGuess = [1, 0.1, 7.5];
    [fitResultAll, ~] = fit(trials', meanCompletionRateAll', sigmoidModel, 'StartPoint', initialGuess);
    plot(trials, feval(fitResultAll, trials), '-', 'Color', hex2rgb('#EDAE92'), 'LineWidth', 2, 'DisplayName', 'Overall Fit');

    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Completion Rate', 'FontSize', 12);
    title('Individual Subject Completion Rates', 'FontSize', 12);
    legend_obj = legend('Location', 'eastoutside', 'FontSize', 12, 'Box', 'off');
    xticks(1:2:15);
    yticks([0.2, 0.5, 0.8, 1]);
    yline(0.8, '--', 'HandleVisibility', 'off');
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
    ylim([0.2, 1.2]);
	set(gcf,'renderer','Painters')
    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);
end




function calculateSessionBasedLearningSpeed(touch_table, block_table, accuracy_threshold)
    % Extract unique subjects and sessions
    subjects = unique(touch_table.Subject);
    
    % Initialize containers for storing data
    all_completed_trials = [];
    all_block_performances = [];
	all_block_error_rate = [];
    session_learning_speeds_good = [];
    session_learning_speeds_poor = [];

    % Loop through each subject
    for subj = 1:numel(subjects)
        subject_data = touch_table(strcmp(touch_table.Subject, subjects{subj}), :);
        sessions = unique(subject_data.SessionNumber);
        
        % Loop through each session
        for sess = 1:numel(sessions)
            session_data = subject_data(subject_data.SessionNumber == sessions(sess), :);
            blocks = unique(session_data.BlockNumber);
            
            session_completed_trials = [];
            session_block_performances = [];
			session_block_error_rate = [];
            
            % Loop through each block
            for blk = 1:numel(blocks)
                block_data = session_data(session_data.BlockNumber == blocks(blk), :);
                % Ensure the block is not swapped and not a repetition
                if all(block_data.isSwapped == 0) && all(block_data.isRepetition == 0)
                    unique_trials = unique(block_data.TrialNumberInBlock);
                    numTrials = numel(unique_trials);
                    completed_trials = zeros(1, numTrials);
					total_number_errors = 0;
                    for trial = 1:numTrials
                        trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
                        completed_trials(trial) = any(trial_data.isRewarded == 1);
						total_number_errors = total_number_errors + sum(contains(trial_data.TouchCategory, 'Error'));
                    end
                    
                    % Pad completed trials to 15 if necessary
                    padded_completed_trials = nan(1, 15);
                    padded_completed_trials(1:length(completed_trials)) = completed_trials;
                    session_completed_trials = [session_completed_trials; padded_completed_trials];
                    session_block_error_rate = [session_block_error_rate; total_number_errors / 15];
                    % Calculate and store block performance
                    block_performance = sum(completed_trials) / numTrials;
                    session_block_performances = [session_block_performances; block_performance];
                end
            end
            
            % Separate good and poor blocks within the session
            % good_block_indices = session_block_performances >= accuracy_threshold;
            % poor_block_indices = ~good_block_indices;
            good_block_indices = session_block_error_rate < 5.33;
            poor_block_indices = session_block_error_rate > 6.6;

            % Calculate mean completion rates for good and poor blocks
            meanCompletionRateGood = nanmean(session_completed_trials(good_block_indices,:), 1);
            meanCompletionRatePoor = nanmean(session_completed_trials(poor_block_indices,:), 1);
            
            % Calculate learning speed for good blocks
            good_learning_speed = NaN;
            for i = 1:length(meanCompletionRateGood)
                if i < length(meanCompletionRateGood) - 3
                	if mean(meanCompletionRateGood(i:i+3)) >= accuracy_threshold && meanCompletionRateGood(i) >= accuracy_threshold
                    	good_learning_speed = i;
                    	break;
					end
				else
					if mean(meanCompletionRateGood(i:end)) >= accuracy_threshold && meanCompletionRateGood(i) >= accuracy_threshold
                    	good_learning_speed = i;
                    	break;
					end
				end
			end

            session_learning_speeds_good = [session_learning_speeds_good; good_learning_speed];
            
            % Calculate learning speed for poor blocks
            poor_learning_speed = NaN;
            for i = 1:length(meanCompletionRatePoor)
                % if all(meanCompletionRatePoor(i:end) >= accuracy_threshold)
                %     poor_learning_speed = i;
                %     break;
                % end
				if i < length(meanCompletionRatePoor) - 3
                	if mean(meanCompletionRatePoor(i:i+3)) >= accuracy_threshold && meanCompletionRatePoor(i) >= accuracy_threshold
                    	poor_learning_speed = i;
                    	break;
					end
				else
					if mean(meanCompletionRatePoor(i:end)) >= accuracy_threshold && meanCompletionRatePoor(i) >= accuracy_threshold
                    	poor_learning_speed = i;
                    	break;
					end
				end
			end

            session_learning_speeds_poor = [session_learning_speeds_poor; poor_learning_speed];
            
            % Accumulate all completed trials and block performances
            all_completed_trials = [all_completed_trials; session_completed_trials];
            all_block_performances = [all_block_performances; session_block_performances];
			all_block_error_rate = [all_block_error_rate; session_block_error_rate];
        end
    end
	
	[~, sortedIndices] = sort(all_block_error_rate, 'descend');
	oneThird = floor(length(all_block_error_rate) / 4);
	% Indices for the top 1/3
	topIndices = sortedIndices(1:oneThird);
	% Indices for the bottom 1/3
	bottomIndices = sortedIndices(end-oneThird+1:end);
	bad_threshold = min(all_block_error_rate(topIndices));
	good_threshold = max(all_block_error_rate(bottomIndices));


    % Calculate overall mean completion rates
    % good_block_indices = all_block_performances >= accuracy_threshold;
    % poor_block_indices = ~good_block_indices;
	good_block_indices = all_block_error_rate < 5.33;
    poor_block_indices = all_block_error_rate > 6.6;

    meanCompletionRateGood = nanmean(all_completed_trials(good_block_indices,:), 1);
    meanCompletionRatePoor = nanmean(all_completed_trials(poor_block_indices,:), 1);

    % Plotting
    canvas_size = [100, 100, 300, 250];
    fig1 = figure('Position', canvas_size);
    hold on;

    trials = 1:size(all_completed_trials, 2);

    % Sigmoid model
    sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', ...
                           'independent', 'x', ...
                           'coefficients', {'a', 'b', 'c'});

    initialGuess = [1, 0.1, trials(round(end / 2))];

    % Fit sigmoid models
    [fitResultGood, ~] = fit(trials', meanCompletionRateGood', sigmoidModel, 'StartPoint', initialGuess);
    [fitResultPoor, ~] = fit(trials', meanCompletionRatePoor', sigmoidModel, 'StartPoint', initialGuess);

    % Calculate confidence intervals
    sem_good = std(all_completed_trials(good_block_indices,:), [], 1, 'omitnan') / sqrt(sum(good_block_indices));
    ci_good = 1.96 * sem_good;
    sem_poor = std(all_completed_trials(poor_block_indices,:), [], 1, 'omitnan') / sqrt(sum(poor_block_indices));
    ci_poor = 1.96 * sem_poor;

    % Plot shaded areas for confidence intervals
    fill([trials fliplr(trials)], [meanCompletionRateGood-ci_good fliplr(meanCompletionRateGood+ci_good)], ...
         hex2rgb('#6C96CC'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([trials fliplr(trials)], [meanCompletionRatePoor-ci_poor fliplr(meanCompletionRatePoor+ci_poor)], ...
         hex2rgb('#DB432C'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Plot data points and fitted curves
    plot(trials, meanCompletionRateGood, '.', 'MarkerSize', 12, 'Color', '#6C96CC', 'HandleVisibility', 'off');
    plot(trials, feval(fitResultGood, trials), 'Color', '#6C96CC', 'LineWidth', 2, 'DisplayName', 'Well Performed Blocks Fit');

    plot(trials, meanCompletionRatePoor, '.', 'MarkerSize', 12, 'Color', '#DB432C', 'HandleVisibility', 'off');
    plot(trials, feval(fitResultPoor, trials), 'Color', '#DB432C', 'LineWidth', 2, 'DisplayName', 'Poorly Performed Blocks Fit');
	
    % Plot average learning speeds
    avg_good_speed = nanmean(session_learning_speeds_good);
	disp('good');
	disp(length(session_learning_speeds_good));
	disp(sum(~isnan(session_learning_speeds_good)));

    sem_good_speed = nanstd(session_learning_speeds_good) / sqrt(sum(~isnan(session_learning_speeds_good)));
    avg_poor_speed = nanmean(session_learning_speeds_poor);
	disp(length(session_learning_speeds_poor));
	disp(sum(~isnan(session_learning_speeds_poor)));
    sem_poor_speed = nanstd(session_learning_speeds_poor) / sqrt(sum(~isnan(session_learning_speeds_poor)));

    errorbar(avg_good_speed, 1.05, sem_good_speed, 'horizontal', 'Color', '#6C96CC', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(avg_good_speed, 1.05, 'v', 'Color', '#6C96CC', 'MarkerFaceColor', '#6C96CC', 'MarkerSize', 8, 'DisplayName','Well Performed Block Learning Speed');

    errorbar(avg_poor_speed, 1.05, sem_poor_speed, 'horizontal', 'Color', '#DB432C', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(avg_poor_speed, 1.05, 'v', 'Color', '#DB432C', 'MarkerFaceColor', '#DB432C', 'MarkerSize', 8, 'DisplayName','Poorly Performed Block Learning Speed');
	
	set(gcf,'renderer','Painters');
    xlabel('Trial Number', 'FontSize', 12);
    ylabel('Completion Rate', 'FontSize', 12);
    title('Trial Completion Rates by Performance', 'FontSize', 12);
    xticks(1:2:15);
    yticks([0.2, 0.5, 0.8, 1]);
    yline(0.8, '--', 'HandleVisibility', 'off');
    legend_obj = legend('Location', 'eastoutside', 'FontSize', 12, 'Box', 'off');
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
	set(gcf,'renderer','Painters')
    ylim([0.2, 1.1]);
    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    % Print out the average learning speeds
    fprintf('Average learning speed for well-performed blocks: %.2f ± %.2f\n', avg_good_speed, sem_good_speed);
    fprintf('Average learning speed for poorly-performed blocks: %.2f ± %.2f\n', avg_poor_speed, sem_poor_speed);
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

function plotAvgSearchTimeByOrdinalPosition(touch_table)
    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Constants for the plot
    ordinalPositions = 1:5; % Assuming 5 as the maximum ordinal position
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;
    correctTouches = touch_table((strcmp(touch_table.TouchCategory, 'Correct') | strcmp(touch_table.TouchCategory, 'correctSelection')) & touch_table.isRepetition == 0 & touch_table.isSwapped == 0, :);
    subjects = unique(touch_table.Subject);
    sessions = unique(touch_table.SessionNumber);

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
    allData = cell(1, length(ordinalPositions));
    maxSearchTime = 0;  % Variable to store the maximum search time
    for position = ordinalPositions
        positionData = [];
        for s = 1:length(subjects)
            for sess = sessions'
                sessionData = correctTouches(correctTouches.TouchObjectCorrectOrdinalPosition == position & ...
                                             strcmp(correctTouches.Subject, subjects{s}) & ...
                                             correctTouches.SessionNumber == sess, :);
                if ~isempty(sessionData)
                    meanSearchTime = mean(sessionData.SearchTime, 'omitnan');
                    positionData = [positionData; position, meanSearchTime, s];
                    maxSearchTime = max(maxSearchTime, meanSearchTime);  % Update max search time
                end
            end
        end
        allData{position} = positionData;
        
        % Plot swarm chart
%         swarmchart(positionData(:,1), positionData(:,2), 36, colors(position,:), 'filled', ...
%                    'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', colors(position,:), 'MarkerEdgeAlpha', 0.2);
        
        % Add shapes for each subject
        for s = 1:length(subjects)
            subjectData = positionData(positionData(:,3) == s, :);
            if ~isempty(subjectData)
                swarmchart(subjectData(:,1), subjectData(:,2), 36, colors(position,:), shapes{mod(s-1, length(shapes))+1}, ...
                        'MarkerEdgeColor', colors(position,:), 'MarkerEdgeAlpha', 0.3, 'LineWidth', 1.5);
            end
        end
    end

    % Add error bars
    for position = ordinalPositions
        positionData = cell2mat(allData(position));
        meanValue = mean(positionData(:,2));
        stdError = std(positionData(:,2)) / sqrt(size(positionData, 1));
		ciError = 1.96*stdError;
		fprintf('%d, %.3f ± %.3f\n', position, meanValue, ciError);
        errorbar(position, meanValue, ciError, 'k', 'LineWidth', 2, 'CapSize', 10);
    end

    % Customize the plot
    title('Search Time Distribution by Ordinal Position', 'FontSize', 14);
    xlabel('Ordinal Position', 'FontSize', 12);
    ylabel('Average Search Time (s)', 'FontSize', 12);
    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
    xlim([0.25, 5.75]);
    ylim([0, maxSearchTime * 1.15]);  % Set y-axis limit to 1.15 times the maximum value
	set(gcf,'renderer','Painters')
    % Add legend with "Subject [initial]" format
    legendLabels = cellfun(@(x) ['Subject ' x(1)], subjects, 'UniformOutput', false);
    legendHandles = gobjects(1, length(subjects));
    for s = 1:length(subjects)
        legendHandles(s) = plot(NaN, NaN, shapes{mod(s-1, length(shapes))+1}, 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 8);
    end
    legend_obj = legend(legendHandles, legendLabels, 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Location = 'eastoutside';

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig, 'Position', new_fig_size);

    hold off;
end


function plotAvgSearchTimeByTrialAndPerformance(touch_table, block_table, accuracy_threshold)
    % Get unique trial numbers
    uniqueTrials = unique(touch_table.TrialNumberInBlock);

    % Initialize arrays for storing mean search times and standard errors
    meanSearchTimes = nan(length(uniqueTrials), 2); % Column 1: Poor, Column 2: Good
    seSearchTimes = nan(length(uniqueTrials), 2); % Standard Error

    % Extract unique subjects
    subjects = unique(touch_table.Subject);
    numSubjects = numel(subjects);
    gray_shades = linspace(0.3, 0.7, numSubjects);

    % Loop through each trial number to calculate mean and SE for good and poor performance
    for i = 1:length(uniqueTrials)
        for perf = [0, 1] % 0 for poor, 1 for good
            if perf == 1
                perfCondition = block_table.BlockLevelAccuracy >= accuracy_threshold;
            else
                perfCondition = block_table.BlockLevelAccuracy < accuracy_threshold;
            end

            % Find the corresponding blocks with the performance condition
            validBlocks = block_table(perfCondition, {'Subject', 'SessionNumber', 'BlockNumber'});

            % Find the indices in touch_table that match the valid blocks and trial number
            trialIndices = ismember(touch_table(:, {'Subject', 'SessionNumber', 'BlockNumber'}), validBlocks) & ...
                touch_table.TrialNumberInBlock == uniqueTrials(i);

            trialData = touch_table(trialIndices, :);

            % Calculate mean search time for the trial
            if ~isempty(trialData)
                meanSearchTimes(i, perf+1) = nanmean(trialData.SearchTime);
                % Calculate standard error
                N = sum(~isnan(trialData.SearchTime)); % Number of non-NaN entries
                seSearchTimes(i, perf+1) = nanstd(trialData.SearchTime) / sqrt(N);
            end
        end
    end

    % Plotting the first figure
    figure; hold on;
    colors = {hex2rgb('#DB432C'), hex2rgb('#6C96CC')}; % Red for poor, Blue for good
    labels = {'Poorly Performed Block', 'Well Performed Block'};

    for perf = [2, 1] % Good performance first for legend ordering
        x = 1:length(uniqueTrials);
        y = meanSearchTimes(:, perf);
        err = seSearchTimes(:, perf) * 1.96; % For 95% CI

        % Main line
        plot(x, y, 'Color',colors{perf}, 'DisplayName', labels{perf}, 'Marker', 'o', 'LineWidth', 2);

        % Shaded area for CI
        fill([x fliplr(x)], [y' + err' fliplr(y' - err')], colors{perf}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    title('Average Search Time by Trial Number with 95% CI', 'FontSize', 16);
    xlabel('Trial Number', 'FontSize', 14);
    ylabel('Average Search Time (s)', 'FontSize', 14);
    legend('Location', 'southeast', 'FontSize', 14);
    grid off;
    hold off;

    % Plotting the second figure
    figure; hold on;
    x = 1:length(uniqueTrials);
    y = mean(meanSearchTimes, 2, 'omitnan'); % Overall mean
    err = mean(seSearchTimes, 2, 'omitnan') * 1.96; % Overall SE

    plot(x, y, 'DisplayName', 'Overall Mean', 'Color', [0.8500, 0.3250, 0.0980], 'Marker', 'o', 'LineWidth', 2);
    fill([x fliplr(x)], [y' + err' fliplr(y' - err')], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    for s = 1:numSubjects
        subject = subjects{s};
        subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
        subject_mean = arrayfun(@(t) nanmean(subject_data.SearchTime(subject_data.TrialNumberInBlock == t)), uniqueTrials);

        plot(x, subject_mean, '-o', 'Color', [gray_shades(s), gray_shades(s), gray_shades(s)], 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', sprintf('%s.', subject(1)));
    end

    title('Overall Average Search Time by Trial Number', 'FontSize', 16);
    xlabel('Trial Number', 'FontSize', 14);
    ylabel('Average Search Time (s)', 'FontSize', 14);
    legend('FontSize', 14, 'Location', 'best');
    grid off;
    hold off;

    figure; hold on;
    maxTrialNumber = 15; % Maximum trial number for plotting

    % Calculate mean and standard error for correct and error touches
    correctIndices = contains(touch_table.TouchCategory, 'Correct', 'IgnoreCase', true);
    errorIndices = contains(touch_table.TouchCategory, 'Error', 'IgnoreCase', true);

    % Arrays to hold means and standard errors
    correctSearchTimes = zeros(1, maxTrialNumber);
    correctSearchSEs = zeros(1, maxTrialNumber);
    errorSearchTimes = zeros(1, maxTrialNumber);
    errorSearchSEs = zeros(1, maxTrialNumber);

    % Calculate means and SEs
    for t = 1:maxTrialNumber
        correctTrialData = touch_table.SearchTime(correctIndices & touch_table.TrialNumberInBlock == t);
        errorTrialData = touch_table.SearchTime(errorIndices & touch_table.TrialNumberInBlock == t);

        correctSearchTimes(t) = nanmean(correctTrialData);
        correctSearchSEs(t) = nanstd(correctTrialData) / sqrt(sum(~isnan(correctTrialData)));

        errorSearchTimes(t) = nanmean(errorTrialData);
        errorSearchSEs(t) = nanstd(errorTrialData) / sqrt(sum(~isnan(errorTrialData)));
    end

    % Plot with shaded error bars for correct touches
    x = 1:maxTrialNumber;
    fill([x, fliplr(x)], [correctSearchTimes + correctSearchSEs * 1.96, fliplr(correctSearchTimes - correctSearchSEs * 1.96)], ...
         'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', HandleVisibility='off');
    plot(x, correctSearchTimes, '-bo', 'LineWidth', 2, 'DisplayName', 'Correct Touches');

    % Plot with shaded error bars for error touches
    fill([x, fliplr(x)], [errorSearchTimes + errorSearchSEs * 1.96, fliplr(errorSearchTimes - errorSearchSEs * 1.96)], ...
         'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', HandleVisibility='off');
    plot(x, errorSearchTimes, '-ro', 'LineWidth', 2, 'DisplayName', 'Error Touches');

    title('Average Search Time by Touch Correctness', 'FontSize', 16);
    xlabel('Trial Number', 'FontSize', 14);
    ylabel('Average Search Time (s)', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 14);
    grid off;
    hold off;

end




function [dataWithPerf, avgSearchTimesByBlock] = classifyPerformance(touch_table, performance_threshold)
% Initialize a table to store block performance
blockPerformance = []; % Stores the performance of each block

uniqueSubject = unique(touch_table.Subject);
for m = 1:length(uniqueSubject)
	subjectdata = touch_table(strcmp(touch_table.Subject, uniqueSubject(m)), :);

	% Extract unique session numbers
	uniqueSessions = unique(subjectdata.SessionNumber);

	for i = 1:length(uniqueSessions)
		sessionData = subjectdata(subjectdata.SessionNumber == uniqueSessions(i), :);
		uniqueBlocks = unique(sessionData.BlockNumber);

		for j = 1:length(uniqueBlocks)
			blockData = sessionData(sessionData.BlockNumber == uniqueBlocks(j), :);
			uniqueTrials = unique(blockData.TrialNumberInBlock);

			rewardedTrialsCount = 0;
			for k = 1:length(uniqueTrials)
				trialData = blockData(blockData.TrialNumberInBlock == uniqueTrials(k), :);
				% Check if the trial is rewarded (considering the trial rewarded if any touch in the trial was rewarded)
				if any(trialData.isRewarded == 1)
					rewardedTrialsCount = rewardedTrialsCount + 1;
				end
			end

			% Calculate the performance metric for the block
			totalTrials = length(uniqueTrials);
			perfMetric = rewardedTrialsCount / totalTrials;

			isGoodPerf = perfMetric >= performance_threshold;

			% Append performance data
			blockPerformance = [blockPerformance; table(uniqueSubject(m), uniqueSessions(i), uniqueBlocks(j), perfMetric, isGoodPerf, 'VariableNames', {'Subject','SessionNumber', 'BlockNumber', 'PerformanceMetric', 'IsGoodPerformance'})];
		end
	end
end

% Join the performance data back to the original data
dataWithPerf = join(touch_table, blockPerformance, 'Keys', {'Subject','SessionNumber', 'BlockNumber'});

% Calculate average search times by block, considering all touches
avgSearchTimesByBlock = varfun(@mean, dataWithPerf, 'InputVariables', 'SearchTime', ...
	'GroupingVariables', {'Subject','SessionNumber', 'BlockNumber', 'IsGoodPerformance'});
end


function plotNextErrorFrequency(touch_table, accuracy_threshold)
% Define the maximum object number
maxObjectNum = 6;
% Extract unique subjects
subjects = unique(touch_table.Subject);
numSubjects = numel(subjects);

markers = {'o', 's', 'd', '^'};
% Initialize containers to count the occurrences of next error objects for good and bad performances
nextErrorFrequencyGood = zeros(maxObjectNum, numSubjects);
nextErrorFrequencyBad = zeros(maxObjectNum, numSubjects);


% Initialize arrays to store the number of trials for good and bad performances
totalTrialsGood = zeros(numSubjects, 1);
totalTrialsBad = zeros(numSubjects, 1);

% Define gray shades for subjects
gray_shades = linspace(0.3, 0.7, numSubjects);

% Loop through each subject
for s = 1:numSubjects
	subject = subjects{s};
	subject_data = touch_table(strcmp(touch_table.Subject, subject), :);
	sessions = unique(subject_data.SessionNumber);

	% Loop through each session for the current subject
	for session = 1:numel(sessions)
		session_data = subject_data(subject_data.SessionNumber == sessions(session), :);
		blocks = unique(session_data.BlockNumber);

		% Loop through each block in the current session
		for block = 1:numel(blocks)
			block_data = session_data(session_data.BlockNumber == blocks(block), :);
			trials = unique(block_data.TrialNumberInBlock);
			% Calculate BlockLevelAccuracy
			unique_trials = unique(block_data.TrialNumberInBlock);
			numRewardedTrials = sum(arrayfun(@(trial) any(block_data.isRewarded(block_data.TrialNumberInBlock == trial)), unique_trials));
			numTrials = numel(unique_trials);
			blockLevelAccuracy = numRewardedTrials / numTrials;
			
			if block_data.isSwapped(1) == 0
				total_number_errors = 0;
            	for trial = 1:numTrials
                	trial_data = block_data(block_data.TrialNumberInBlock == unique_trials(trial), :);
					total_number_errors = total_number_errors + sum(contains(trial_data.TouchCategory, 'Error'));
				end
				block_level_error_rate = total_number_errors/15;
			end

			% Loop through each trial in the current block
			for trial = 1:numel(trials)
				trial_data = block_data(block_data.TrialNumberInBlock == trials(trial), :);

				% Check if the trial is swapped
				if trial_data.isSwapped(1) == 1
					% Loop through each touch in the current trial
					for touch = 1:height(trial_data) - 1

						if trial_data.CurrentState(touch) == 1 && (strcmpi(trial_data.TouchCategory{touch}, 'correct') || strcmpi(trial_data.TouchCategory{touch}, 'correctSelection')) ...
								&& trial_data.CurrentState(touch+1) == 0 && trial_data.TouchObjectCorrectOrdinalPosition(touch+1) == 4
							% Find the touch after making an error (touch = correct at 1, touch + 1 = error at slot 2, touch + 2 = error correction)
							nextTouchIndex = touch + 2;

							if  nextTouchIndex + 1 <= height(trial_data) && strcmpi(trial_data.TouchCategory{nextTouchIndex}, 'RetouchCorrect')
								nextTouchIndex = nextTouchIndex + 1;
							end


							if nextTouchIndex <= height(trial_data)
								nextErrorObjectPosition = trial_data.TouchObjectCorrectOrdinalPosition(nextTouchIndex);
								if nextErrorObjectPosition == 0
									nextErrorObjectPosition = 6;
								end
								% if blockLevelAccuracy >= accuracy_threshold
								% 	nextErrorFrequencyGood(nextErrorObjectPosition, s) = nextErrorFrequencyGood(nextErrorObjectPosition,s) + 1;
								% 	totalTrialsGood(s,1) = totalTrialsGood(s,1) + 1;
								% else
								% 	nextErrorFrequencyBad(nextErrorObjectPosition, s) = nextErrorFrequencyBad(nextErrorObjectPosition,s) + 1;
								% 	totalTrialsBad(s,1) = totalTrialsBad(s,1) + 1;
								% end

								if block_level_error_rate <= 5.33
									nextErrorFrequencyGood(nextErrorObjectPosition, s) = nextErrorFrequencyGood(nextErrorObjectPosition,s) + 1;
									totalTrialsGood(s,1) = totalTrialsGood(s,1) + 1;
								elseif block_level_error_rate >= 6.6
									nextErrorFrequencyBad(nextErrorObjectPosition, s) = nextErrorFrequencyBad(nextErrorObjectPosition,s) + 1;
									totalTrialsBad(s,1) = totalTrialsBad(s,1) + 1;
								end
							end
						end
					end
				end
			end
		end
	end
end

temp = nextErrorFrequencyGood(4,:);
nextErrorFrequencyGood(4,:) = nextErrorFrequencyGood(2,:);
nextErrorFrequencyGood(2,:) = temp;

temp = nextErrorFrequencyBad(4,:);
nextErrorFrequencyBad(4,:) = nextErrorFrequencyBad(2,:);
nextErrorFrequencyBad(2,:) = temp;

% Calculate the proportions and standard errors for good and bad performances
proportionGood = sum(nextErrorFrequencyGood, 2) ./ sum(totalTrialsGood);
proportionBad = sum(nextErrorFrequencyBad, 2) ./ sum(totalTrialsBad);
SE_good = sqrt(proportionGood .* (1 - proportionGood) ./ sum(totalTrialsGood));
SE_bad = sqrt(proportionBad .* (1 - proportionBad) ./ sum(totalTrialsBad));

    canvas_size = [100, 100, 300, 250];

	gray_color = [0.7,0.7,0.7];
    % Plot 1: Proportions for good and bad performance
    fig1 = figure;
    fig1.Position = canvas_size;
    hold on;

    for s = 1:numSubjects
        plot(1:maxObjectNum, nextErrorFrequencyGood(:, s) ./ totalTrialsGood(s), ['-', markers{s}], 'Color', gray_color, 'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
        plot(1:maxObjectNum, nextErrorFrequencyBad(:, s) ./ totalTrialsBad(s), ['--', markers{s}], 'Color', gray_color, 'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end
    % Plot overall results
    errorbar(1:maxObjectNum, proportionGood, SE_good, '-', 'DisplayName', 'Well Performed Block', 'LineWidth', 2, 'Color', hex2rgb('#6C96CC'));
    errorbar(1:maxObjectNum, proportionBad, SE_bad, '-', 'DisplayName', 'Poorly Performed Block', 'LineWidth', 2, 'Color',  hex2rgb('#DB432C'));
	disp('Good');
	disp(proportionGood');
	disp(SE_good');
	disp('Poor');
	disp(proportionBad');
	disp(SE_bad');
    legend_obj1 = legend('Location', 'eastoutside', 'FontSize', 12);
    legend_obj1.EdgeColor = 'none';
    title('Proportion of Errors by Object', 'FontSize', 14);
    xlabel('Object ID', 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 12);
    xticks(1:maxObjectNum);
    ylim_values = [min(min(proportionGood - SE_good), min(proportionBad - SE_bad)) - 0.05, ...
                   max(max(proportionGood + SE_good), max(proportionBad + SE_bad)) * 1.5];
    ylim(ylim_values);
	xlim([0.5,6.5]);
    xticklabels({'A', 'B', 'C', 'D', 'E', 'Distractor'});

    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos1 = get(legend_obj1, 'Position');
    legend_width1 = legend_pos1(3);
    new_fig_width1 = canvas_size(3) + legend_width1 * canvas_size(3);
    new_fig_size1 = [canvas_size(1), canvas_size(2), new_fig_width1, canvas_size(4)];
    set(fig1, 'Position', new_fig_size1);

    hold off;

    % Plot 2: Differences in proportions
    fig2 = figure;
    fig2.Position = canvas_size;
    hold on;
	
	differenceProportions = proportionGood - proportionBad;
	SE_diff = sqrt(SE_good.^2 + SE_bad.^2);

    for s = 1:numSubjects
        diffGood = nextErrorFrequencyGood(:, s) / totalTrialsGood(s);
        diffBad = nextErrorFrequencyBad(:, s) / totalTrialsBad(s);
        diffProportion = diffGood - diffBad;
        plot(1:maxObjectNum, diffProportion, ['-', markers{s}], 'Color', gray_color, 'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
    end

    errorbar(1:maxObjectNum, differenceProportions, SE_diff, '-s', 'DisplayName', 'Difference in Proportions', 'LineWidth', 2, 'Color', 'k');
	yline(0, '--', 'HandleVisibility', 'off');

    ylim_diff = [min(min(differenceProportions - 2*SE_diff), -0.2), ...
                 max(max(differenceProportions + 3*SE_diff), 0.3)];
    ylim(ylim_diff);
    legend_obj2 = legend('Location', 'eastoutside', 'FontSize', 12);
    legend_obj2.EdgeColor = 'none';

    title('Difference in Proportions: Well Performed - Poorly Performed', 'FontSize', 14);
    xlabel('Object ID', 'FontSize', 12);
    ylabel('Difference in Proportion', 'FontSize', 12);
    
	xticks(1:maxObjectNum);
    xticklabels({'A', 'B', 'C', 'D', 'E', 'Distractor'});

    ax = gca;
    ax.FontSize = 12;
    ax.TickDir = 'out';
    set(gca, 'box', 'off');

    % Adjust figure size to accommodate legend
    legend_pos2 = get(legend_obj2, 'Position');
    legend_width2 = legend_pos2(3);
    new_fig_width2 = canvas_size(3) + legend_width2 * canvas_size(3);
    new_fig_size2 = [canvas_size(1), canvas_size(2), new_fig_width2, canvas_size(4)];
    set(fig2, 'Position', new_fig_size2);

    hold off;
end

function plotOrdinalPositionInference(blockTable)
% Extract unique subjects
subjects = unique(blockTable.Subject);
numSubjects = numel(subjects);

% Extract unique session numbers
uniqueSessions = unique(blockTable.SessionNumber);

% Initialize an array to store the average swappedErrorCorrection for each session for each subject
averageSwappedErrorCorrection = nan(length(uniqueSessions), numSubjects);

% Loop through each subject
for s = 1:numSubjects
	subject = subjects{s};

	% Loop through each session
	for i = 1:length(uniqueSessions)
		sessionNumber = uniqueSessions(i);

		% Get the swappedErrorCorrection values for the current session and subject
		sessionData = blockTable.SwappedErrorCorrection(blockTable.SessionNumber == sessionNumber & ...
			blockTable.isSwapped == 1 & strcmp(blockTable.Subject, subject));

		% Calculate the nanmean of swappedErrorCorrection values for the current session and subject
		averageSwappedErrorCorrection(i, s) = nanmean(sessionData);
	end
end

% Calculate the overall average and smoothed version using a moving average of 3 sessions
overallAverageSwappedErrorCorrection = nanmean(averageSwappedErrorCorrection, 2);
smoothedSwappedErrorCorrection = movmean(overallAverageSwappedErrorCorrection, [0 2], 'omitnan');

% Perform linear regression for the overall data
[p, S] = polyfit(uniqueSessions, overallAverageSwappedErrorCorrection', 1);

% Calculate the regression line
yfit = polyval(p, uniqueSessions);

% Calculate R-squared value
yresid = overallAverageSwappedErrorCorrection' - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(overallAverageSwappedErrorCorrection) - 1) * var(overallAverageSwappedErrorCorrection, 'omitnan');
rsq = 1 - SSresid / SStotal;

% Calculate regression statistics using regstats
stats = regstats(overallAverageSwappedErrorCorrection', uniqueSessions, 'linear');

% Get p-values for regression coefficients
p_value = stats.tstat.pval;

% Create the figure
figure;
hold on;

% Define gray shades for subjects
gray_shades = linspace(0, 0.9, numSubjects);
dot_shape = {'o','x','*','^'};

% Plot individual subject data
for s = 1:numSubjects
	% 		if ~strcmp(subjects{s},'Jotun') && ~strcmp(subjects{s},'Bard')
	% 			continue
	% 		end
	plot(uniqueSessions, averageSwappedErrorCorrection(:, s), ['-', dot_shape{s}], 'Color', [gray_shades(s), gray_shades(s), gray_shades(s)], ...
		'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', ['Subject ', subjects{s}(1)]);
	hold on
end

% Plot the overall smoothed data
plot(uniqueSessions, smoothedSwappedErrorCorrection, '-x', 'LineWidth', 3, 'Color', 'b', 'MarkerSize', 8, 'DisplayName', 'Overall Smoothed');

% Plot the regression line
plot(uniqueSessions, yfit, '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Linear Fit');

% Display regression statistics
fprintf('Overall Regression:\n');
fprintf('Coefficients: %.4f (slope), %.4f (intercept)\n', p);
fprintf('R-squared: %.4f\n', rsq);
fprintf('p-value: %.4f\n', p_value(2));

xlabel('Session Number');
ylabel('Average Ordinal Object Inference');
title('Average Ordinal Object Inference Over Sessions');
legend('show');
grid on;

% Customize plot appearance
ax = gca;  % Get current axis handle
ax.FontSize = 14;  % Increase font size of axis labels
ax.TickDir = 'out';
set(gca, 'box', 'off');

hold off;
end


function plotContextConditionAccuracy(blockTable)
% Helper function to calculate mean and SEM
	function [avgAccRepeat, avgAccNonRepeat, semRepeat, semNonRepeat] = calcStats(data)
		avgAccRepeat = mean(data.BlockLevelAccuracy(data.isRepetition == true));
		avgAccNonRepeat = mean(data.BlockLevelAccuracy(data.isRepetition == false));
		semRepeat = std(data.BlockLevelAccuracy(data.isRepetition == true)) / sqrt(sum(data.isRepetition == true));
		semNonRepeat = std(data.BlockLevelAccuracy(data.isRepetition == false)) / sqrt(sum(data.isRepetition == false));
	end

% Helper function to calculate differences and SEM
	function [diff, semDiff] = calcDiffAndSEM(avgRepeat, avgNonRepeat, semRepeat, semNonRepeat)
		diff = avgRepeat - avgNonRepeat;
		semDiff = sqrt(semRepeat^2 + semNonRepeat^2);
	end

% Calculate statistics for isColorSame conditions
[avgAccColorSame0Repeat, avgAccColorSame0NonRepeat, semColorSame0Repeat, semColorSame0NonRepeat] = calcStats(blockTable(blockTable.isColorSame == 0, :));
[avgAccColorSame1Repeat, avgAccColorSame1NonRepeat, semColorSame1Repeat, semColorSame1NonRepeat] = calcStats(blockTable(blockTable.isColorSame == 1, :));

% Calculate statistics for isBackgroundSame conditions
[avgAccBackgroundSame0Repeat, avgAccBackgroundSame0NonRepeat, semBackgroundSame0Repeat, semBackgroundSame0NonRepeat] = calcStats(blockTable(blockTable.isBackgroundSame == 0, :));
[avgAccBackgroundSame1Repeat, avgAccBackgroundSame1NonRepeat, semBackgroundSame1Repeat, semBackgroundSame1NonRepeat] = calcStats(blockTable(blockTable.isBackgroundSame == 1, :));

% Calculate differences and standard errors
[diffColorSame0, semDiffColorSame0] = calcDiffAndSEM(avgAccColorSame0Repeat, avgAccColorSame0NonRepeat, semColorSame0Repeat, semColorSame0NonRepeat);
[diffColorSame1, semDiffColorSame1] = calcDiffAndSEM(avgAccColorSame1Repeat, avgAccColorSame1NonRepeat, semColorSame1Repeat, semColorSame1NonRepeat);
[diffBackgroundSame0, semDiffBackgroundSame0] = calcDiffAndSEM(avgAccBackgroundSame0Repeat, avgAccBackgroundSame0NonRepeat, semBackgroundSame0Repeat, semBackgroundSame0NonRepeat);
[diffBackgroundSame1, semDiffBackgroundSame1] = calcDiffAndSEM(avgAccBackgroundSame1Repeat, avgAccBackgroundSame1NonRepeat, semBackgroundSame1Repeat, semBackgroundSame1NonRepeat);

% Calculate statistics for both conditions same and any condition different
bothSame = blockTable(blockTable.isColorSame == 1 & blockTable.isBackgroundSame == 1, :);
[avgAccBothSameRepeat, avgAccBothSameNonRepeat, semBothSameRepeat, semBothSameNonRepeat] = calcStats(bothSame);
anyDifferent = blockTable(blockTable.isColorSame ~= 1 | blockTable.isBackgroundSame ~= 1, :);
[avgAccAnyDifferentRepeat, avgAccAnyDifferentNonRepeat, semAnyDifferentRepeat, semAnyDifferentNonRepeat] = calcStats(anyDifferent);

[diffBothSame, semDiffBothSame] = calcDiffAndSEM(avgAccBothSameRepeat, avgAccBothSameNonRepeat, semBothSameRepeat, semBothSameNonRepeat);
[diffAnyDifferent, semDiffAnyDifferent] = calcDiffAndSEM(avgAccAnyDifferentRepeat, avgAccAnyDifferentNonRepeat, semAnyDifferentRepeat, semAnyDifferentNonRepeat);

% Prepare data for box plot
subjects = unique(blockTable.Subject);
gray_shades = linspace(0.3, 0.7, numel(subjects));
subjectDiffBothSame = NaN(numel(subjects), 1);
subjectDiffAnyDifferent = NaN(numel(subjects), 1);

for i = 1:numel(subjects)
	subject = subjects{i};
	subjDataBothSame = bothSame(strcmp(bothSame.Subject, subject), :);
	subjDataAnyDifferent = anyDifferent(strcmp(anyDifferent.Subject, subject), :);
	[subjAvgAccBothSameRepeat, subjAvgAccBothSameNonRepeat] = deal(mean(subjDataBothSame.BlockLevelAccuracy(subjDataBothSame.isRepetition == true)), ...
		mean(subjDataBothSame.BlockLevelAccuracy(subjDataBothSame.isRepetition == false)));
	[subjAvgAccAnyDifferentRepeat, subjAvgAccAnyDifferentNonRepeat] = deal(mean(subjDataAnyDifferent.BlockLevelAccuracy(subjDataAnyDifferent.isRepetition == true)), ...
		mean(subjDataAnyDifferent.BlockLevelAccuracy(subjDataAnyDifferent.isRepetition == false)));
	subjectDiffBothSame(i) = subjAvgAccBothSameRepeat - subjAvgAccBothSameNonRepeat;
	subjectDiffAnyDifferent(i) = subjAvgAccAnyDifferentRepeat - subjAvgAccAnyDifferentNonRepeat;
end

% Create the first figure with two horizontal subplots
fig1 = figure;
set(fig1, 'Position', [100, 100, 1200, 400]);

% Plotting function
	function createSubplot(position, diffs, sems, titleText, xlabelText, xticksLabels)
		subplot(1, 2, position);
		hold on;
		bar(diffs, 'FaceColor', 'flat');
		errorbar(diffs, sems, 'k', 'LineStyle', 'none');
		title(titleText);
		xlabel(xlabelText);
		ylabel('Difference in Avg Accuracy');
		xticks(1:length(xticksLabels));
		xticklabels(xticksLabels);
		customizePlot();
		hold off;
	end

% First subplot: Average Accuracy for Repeat - Non-Repeat for isColorSame
createSubplot(1, [diffColorSame0, diffColorSame1], [semDiffColorSame0, semDiffColorSame1], 'Avg Accuracy for Repeat - Non-Repeat (isColorSame)', 'isColorSame', {'0', '1'});

% Second subplot: Average Accuracy for Repeat - Non-Repeat for isBackgroundSame
createSubplot(2, [diffBackgroundSame0, diffBackgroundSame1], [semDiffBackgroundSame0, semDiffBackgroundSame1], 'Avg Accuracy for Repeat - Non-Repeat (isBackgroundSame)', 'isBackgroundSame', {'0', '1'});

sgtitle('Average Accuracy Comparison for isColorSame and isBackgroundSame');

% Create a new figure for both conditions same vs any condition different (box plot)
fig2 = figure;
set(fig2, 'Position', [200, 200, 600, 400]);
hold on;

% Create the box plots
boxplot([subjectDiffAnyDifferent, subjectDiffBothSame], 'Labels', {'Different', 'Same'});

% Plot subject dots and error bars in gray scale
for i = 1:numel(subjects)
    plot(1, subjectDiffAnyDifferent(i), '*', 'Color', [gray_shades(i), gray_shades(i), gray_shades(i)], 'MarkerSize', 8);
    plot(2, subjectDiffBothSame(i), 'o', 'Color', [gray_shades(i), gray_shades(i), gray_shades(i)], 'MarkerSize', 8);
end

% Plot overall results
errorbar(1, diffAnyDifferent, semDiffAnyDifferent, 'o', 'MarkerSize', 10, 'Color', 'b', 'LineWidth', 2, 'CapSize', 15);
errorbar(2, diffBothSame, semDiffBothSame, 'o', 'MarkerSize', 10, 'Color', 'r', 'LineWidth', 2, 'CapSize', 15);

% Perform statistical tests
[h_any_different, p_any_different, ~, stats_any_different] = ttest(subjectDiffAnyDifferent);
[h_both_same, p_both_same, ~, stats_both_same] = ttest(subjectDiffBothSame);
[h_comparison, p_comparison, ~, stats_comparison] = ttest(subjectDiffAnyDifferent, subjectDiffBothSame);

% Function to get significance stars
function stars = get_stars(p)
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end

% Add significance stars and lines
y_max = max([subjectDiffAnyDifferent; subjectDiffBothSame]);
y_text = y_max + 0.03;

% Any Different vs 0
text(1, y_text, get_stars(p_any_different), 'HorizontalAlignment', 'center', 'FontSize', 14);
line([1 1], [y_max+0.01 y_text-0.005], 'Color', 'k');

% Both Same vs 0
text(2, y_text, get_stars(p_both_same), 'HorizontalAlignment', 'center', 'FontSize', 14);
line([2 2], [y_max+0.01 y_text-0.005], 'Color', 'k');

% Any Different vs Both Same
y_comparison = y_text + 0.03;
line([1 2], [y_comparison y_comparison], 'Color', 'k');
text(1.5, y_comparison + 0.01, get_stars(p_comparison), 'HorizontalAlignment', 'center', 'FontSize', 14);

% % Add statistical results to the plot
% text(0.7, 0.18, sprintf('Any Different vs 0:\nt(%.0f) = %.2f, p = %.3f', stats_any_different.df, stats_any_different.tstat, p_any_different), 'FontSize', 10);
% text(0.7, 0.15, sprintf('Both Same vs 0:\nt(%.0f) = %.2f, p = %.3f', stats_both_same.df, stats_both_same.tstat, p_both_same), 'FontSize', 10);
% text(0.7, 0.12, sprintf('Any Different vs Both Same:\nt(%.0f) = %.2f, p = %.3f', stats_comparison.df, stats_comparison.tstat, p_comparison), 'FontSize', 10);

title('Avg Accuracy for Repeat - Non-Repeat');
xlabel('Condition');
ylabel('Difference in Avg Accuracy');
ylim([-0.05, 0.25]);
xlim([0.5, 2.5]);
customizePlot();
hold off;


	function customizePlot()
		ax = gca;  % Get current axis handle
		ax.FontSize = 14;  % Increase font size of axis labels
		ax.TickDir = "out";
		set(gca, 'box', 'off');
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


function [avg_speed_learning, se_speed_learning, avg_good_blocks, se_good_blocks, avg_poor_blocks, se_poor_blocks] = calculateLearningPerformance(touch_table)
% Keep unique rows based on specific columns
unique_rows = unique(touch_table(:, {'Subject', 'SessionNumber', 'BlockNumber', 'TrialNumberInBlock', 'isRepetition', 'isSwapped', 'isRewarded'}));

all_completed_trials = [];
all_block_performances = [];
subjects = unique(unique_rows.Subject);

% Initialize a structure to store subject-wise data
subject_data = struct();

for i = 1:numel(subjects)
	subject = subjects{i};
	subject_rows = unique_rows(strcmp(unique_rows.Subject, subject), :);
	sessions = unique(subject_rows.SessionNumber);

	for j = 1:numel(sessions)
		session = sessions(j);
		session_rows = subject_rows(subject_rows.SessionNumber == session, :);
		blocks = unique(session_rows.BlockNumber);

		for block = blocks'
			block_indices = session_rows.BlockNumber == block;
			block_rows = session_rows(block_indices, :);

			if any(block_rows.isSwapped) || any(block_rows.isRepetition)
				continue
			end

			completed_trials = block_rows.isRewarded;
			padded_completed_trials = nan(1, 15);
			padded_completed_trials(1:length(completed_trials)) = completed_trials';
			all_completed_trials = [all_completed_trials; padded_completed_trials];

			block_performance = mean(completed_trials);
			all_block_performances = [all_block_performances; block_performance];

			% Store subject-wise data
			if ~isfield(subject_data, subject)
				subject_data.(subject) = [];
			end
			subject_data.(subject) = [subject_data.(subject); padded_completed_trials];
		end
	end
end

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

good_block_indices = all_block_performances >= 0.7;
poor_block_indices = ~good_block_indices;

avg_good_blocks = mean(learning_trial(good_block_indices), 'omitnan');
se_good_blocks = std(learning_trial(good_block_indices), 'omitnan') / sqrt(sum(good_block_indices));

avg_poor_blocks = mean(learning_trial(poor_block_indices), 'omitnan');
se_poor_blocks = std(learning_trial(poor_block_indices), 'omitnan') / sqrt(sum(poor_block_indices));

meanCompletionRateGood = nanmean(all_completed_trials(good_block_indices,:), 1);
meanCompletionRateBad = nanmean(all_completed_trials(poor_block_indices,:), 1);
meanCompletionRateAll = nanmean(all_completed_trials, 1);

trials = 1:size(all_completed_trials, 2);

% Sigmoid model
sigmoidModel = fittype('a / (1 + exp(-b * (x - c)))', ...
	'independent', 'x', ...
	'coefficients', {'a', 'b', 'c'});

% Fit the model to data for good, poor, and all blocks
initialGuess = [1, 0.1, trials(round(end / 2))];

% Good blocks
[fitResultGood, ~] = fit(trials', meanCompletionRateGood', sigmoidModel, 'StartPoint', initialGuess);

% Poor blocks
[fitResultBad, ~] = fit(trials', meanCompletionRateBad', sigmoidModel, 'StartPoint', initialGuess);

% All blocks
[fitResultAll, ~] = fit(trials', meanCompletionRateAll', sigmoidModel, 'StartPoint', initialGuess);

% Plotting
figure;
hold on;

% For good blocks
sem_good = std(all_completed_trials(good_block_indices,:), [], 1, 'omitnan') / sqrt(sum(good_block_indices));
ci_good = 1.96 * sem_good; % 95% CI

% For poor blocks
sem_poor = std(all_completed_trials(poor_block_indices,:), [], 1, 'omitnan') / sqrt(sum(poor_block_indices));
ci_poor = 1.96 * sem_poor; % 95% CI

% For all blocks
sem_all = std(all_completed_trials, [], 1, 'omitnan') / sqrt(size(all_completed_trials, 1));
ci_all = 1.96 * sem_all; % 95% CI

% For good blocks
lowerBoundGood = meanCompletionRateGood - ci_good;
upperBoundGood = meanCompletionRateGood + ci_good;
fill([trials fliplr(trials)], [lowerBoundGood fliplr(upperBoundGood)], hex2rgb('#6C96CC'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% For poor blocks
lowerBoundPoor = meanCompletionRateBad - ci_poor;
upperBoundPoor = meanCompletionRateBad + ci_poor;
fill([trials fliplr(trials)], [lowerBoundPoor fliplr(upperBoundPoor)], hex2rgb('#DB432C'), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot individual subject lines
subject_names = fieldnames(subject_data);
num_subjects = numel(subject_names);
gray_shades = linspace(0, 0.8, num_subjects);
subject_labels = cell(num_subjects, 1);

for i = 1:num_subjects
	subj_trials = subject_data.(subject_names{i});
	subj_mean_trials = nanmean(subj_trials, 1);
	plot(trials, subj_mean_trials, '-', 'Color', [gray_shades(i), gray_shades(i), gray_shades(i)], 'LineWidth', 1, 'DisplayName',subject_names{i});
	subject_labels{i} = ['Monkey ', subject_names{i}(1)];
end



% Plot the overall data fits
plot(trials, meanCompletionRateGood, '.', 'MarkerSize', 12, 'Color', '#6C96CC', 'HandleVisibility', 'off');
plot(trials, feval(fitResultGood, trials), 'Color', '#6C96CC', 'LineWidth', 2);

plot(trials, meanCompletionRateBad, '.', 'MarkerSize', 12, 'Color', '#DB432C', 'HandleVisibility', 'off');
plot(trials, feval(fitResultBad, trials), 'Color', '#DB432C', 'LineWidth', 2);

xlabel('Trial Number', 'FontSize', 20);
ylabel('Completion Rate', 'FontSize', 14);
title('Trial Completion Rates', 'FontSize', 24);
xticks(trials);
yticks([0, 0.5, 0.8, 1]);
yline(0.8, '--', 'HandleVisibility', 'off');
legend([subject_labels; 'Well Performed Blocks'; 'Poorly Performed Blocks'], 'Location', 'southeast', 'FontSize', 14);
set(gca, 'FontSize', 20);
ylim([0, 1]);
    ax = gca; % Get current axis handle
    ax.FontSize = 12; % Increase font size of axis labels
    ax.TickDir = 'out';
    set(gca, 'box', 'off');
hold off;
end

function rgb = hex2rgb(hexStr)
if hexStr(1) ~= '#'
	error('Hex string must start with "#"');
end

hexStr(1) = [];

if numel(hexStr) == 3
	hexStr = repmat(hexStr, 2, 1);
end

r = hex2dec(hexStr(1:2));
g = hex2dec(hexStr(3:4));
b = hex2dec(hexStr(5:6));

rgb = [r g b] / 255;
end

function fixed_touch_table = fixRuleBreakingErrors(touch_table)
    % Sort the table by Subject, SessionNumber, BlockNumber, TrialNumberInBlock, and TouchNumber
    sorted_table = sortrows(touch_table, {'Subject', 'SessionNumber', 'BlockNumber', 'TrialNumberInBlock', 'TouchNumber'});
    
    % Initialize variables to track the state
    current_trial = [0, 0, 0, 0]; % [Subject, SessionNumber, BlockNumber, TrialNumberInBlock]
    correct_selection_made = false;
    
    % Iterate through all rows
    for i = 1:height(sorted_table)
        row = sorted_table(i, :);
        
        % Check if this is a new trial
        if ~isequal(current_trial, [row.Subject, row.SessionNumber, row.BlockNumber, row.TrialNumberInBlock])
            current_trial = [row.Subject, row.SessionNumber, row.BlockNumber, row.TrialNumberInBlock];
            correct_selection_made = false;
        end
        
        % If correct selection hasn't been made yet, fix the error type if needed
        if ~correct_selection_made
            if strcmp(row.TouchCategory, 'ruleBreakingError')
                if row.TouchObjectCorrectOrdinalPosition == 0
                    sorted_table.TouchCategory{i} = 'distractorRuleAbidingError';
                else
                    sorted_table.TouchCategory{i} = 'ruleAbidingError';
                end
            elseif strcmp(row.TouchCategory, 'correctSelection')
                correct_selection_made = true;
            end
        end
    end
    
    % Return the fixed table
    fixed_touch_table = sorted_table;
    
    % Report the number of fixes made
    num_fixes = sum(~strcmp(touch_table.TouchCategory, fixed_touch_table.TouchCategory));
    fprintf('Number of error types fixed: %d\n', num_fixes);
end
