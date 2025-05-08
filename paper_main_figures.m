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


%% Plotting

%% Figure 1

% Figure 1C
% plotLastFourTrialsAccuracyDots(touch_table);

% Figure 1D
% calculateSessionLearningSpeed(touch_table, block_table, 0.8);

% Figure 1E
% calculateSessionBasedLearningSpeed(touch_table, block_table, 0.8);


%% Figure 2

% Figure 2B
% plotSessionWiseTransitionProbability(touch_table);

% Figure 2C
% plotNextErrorFrequency(touch_table);

% Figure 2F
% plotSessionWiseAnticipatedSwappingTransitionProbability(touch_table);

% Figure 2G
% plotSessionWiseSecondToThirdTransitionProbability(touch_table);



%% Figure 4

% Figure 4ABC
% plotTransitionProbabilitiesWithAccuracyFilter(touch_table, 0)

% Figure 4DEF
% plotTransitionProbabilitiesWithAccuracyFilter(touch_table, 1)

% Figure 4JHI
% plotTransitionProbabilitiesWithAccuracyFilter(touch_table, 2)

% Figure 4J
% plotTransitionProbabilitiesFirstTrial(touch_table)

% Figure 4I
% plotTransitionDifferencesByTrial(touch_table)



%% Figure 5
% Figure 5A
% plotTrialCompletionRatesRepeat(touch_table);
% Figure 5B
% plotAnticipatedSwappingThreeConditions(block_table)
% Figure 5C
% plotErrorInferenceInRepeat(block_table);

% Figure 5D
% plotInitialPerformaceVsAnticipatedSwappingTwoConditions(block_table)
% Figure 5E
% plotBlockCompletionVsAnticipatedSwapping(block_table);
% Figure 5H
% plotAnticipatedSwappingOverSession(block_table);

% Figure 5F
% plotInitialPerformaceVsOrdinalObjectIdentityTwoConditions(block_table);
% Figure 5J
% plotBlockCompletionVsErrorCorrection(block_table);
% Figure 5I
% plotOrdinalPositionInferenceOverSession(block_table);


%% Figure 6

% Figure 6B
% WM_Acc_WWW_ACC_MixedModel(WM_WWW_Interaction)
% Figure 6C
% WM_Acc_AnticipatedSwapping_MixedModel(WM_WWW_Interaction)

% Figure 6D
% plot_block_accuracy_difference(block_table);
% % Figure 6E
% plot_ordinal_identity_inference_difference(block_table);
% % Figure 6F
% plot_anticipated_swapping_difference(block_table);














function WM_Acc_AnticipatedSwapping_MixedModel(WM_WWW_Interaction)
    % Ensure the required toolbox is available
    if ~license('test', 'Statistics_Toolbox')
        error('This function requires the Statistics and Machine Learning Toolbox.');
    end

    % Prepare data
    subjects = WM_WWW_Interaction.Subject;
    wm_acc = WM_WWW_Interaction.OverallAccuracy;
    anticipated_swapping = WM_WWW_Interaction.AvgAnticipatedSwapping;
    session_number = WM_WWW_Interaction.SessionNumber;

    % Remove any rows with NaN values
    valid_idx = ~isnan(wm_acc) & ~isnan(anticipated_swapping) & ~isnan(session_number);
    subjects = subjects(valid_idx);
    wm_acc = wm_acc(valid_idx);
    anticipated_swapping = anticipated_swapping(valid_idx);
    session_number = session_number(valid_idx);

    % Create table for mixed-effects model
    tbl = table(subjects, wm_acc, anticipated_swapping, session_number, 'VariableNames', {'Subject', 'WM_Accuracy', 'Anticipated_Swapping', 'SessionNumber'});

    % Fit linear mixed-effects model as specified
    lme = fitlme(tbl, 'Anticipated_Swapping ~ WM_Accuracy + SessionNumber + WM_Accuracy*SessionNumber + (1 + SessionNumber|Subject)');

    % Display model results
    disp(lme);

    % Extract fixed effects and their p-values
    [beta, names, stats] = fixedEffects(lme);
    overall_intercept = beta(1);
    overall_wm_slope = beta(2);
    overall_session_slope = beta(3);
    overall_interaction = beta(4);
    overall_wm_p = stats.pValue(2);
    overall_session_p = stats.pValue(3);
    overall_interaction_p = stats.pValue(4);

    % Calculate individual correlations and slopes
    unique_subjects = unique(subjects);
    individual_correlations = zeros(length(unique_subjects), 1);
    individual_slopes = zeros(length(unique_subjects), 1);
    individual_p_values = zeros(length(unique_subjects), 1);

    fprintf('\nIndividual Subject Correlations and Slopes:\n');
    fprintf('--------------------------------------------\n');
    fprintf('Subject\tCorrelation\tSlope\t\tp-value\n');

    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_ant_swap = anticipated_swapping(subj_idx);

        % Calculate correlation
        [r, p] = corr(subj_wm_acc, subj_ant_swap);
        individual_correlations(i) = r;
        individual_p_values(i) = p;

        % Calculate slope
        coeffs = polyfit(subj_wm_acc, subj_ant_swap, 1);
        individual_slopes(i) = coeffs(1);

        fprintf('%s\t%.4f\t\t%.4f\t\t%.4f\n', subj, r, coeffs(1), p);
    end

    % Define canvas size
    canvas_size = [100, 100, 300, 250];

    % Create scatter plot
    fig1 = figure('Position', canvas_size);
    set(gcf, 'renderer', 'Painters');
    hold on;

    % Define colors for each subject
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;

    % Define markers for each subject
    markers = {'o', 's', 'd', '^', 'v'};

    % Plot data points and subject-specific lines
    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_ant_swap = anticipated_swapping(subj_idx);

        % Plot scatter for this subject
        scatter(subj_wm_acc, subj_ant_swap, 50, 'filled', 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'none', 'Marker', markers{i}, 'LineWidth', 1);

        % Plot subject-specific line (using individual slope)
        x_range = linspace(min(subj_wm_acc), max(subj_wm_acc), 100);
        y_fit = individual_slopes(i) * x_range + mean(subj_ant_swap) - individual_slopes(i) * mean(subj_wm_acc);
        plot(x_range, y_fit, 'Color', colors(i,:), 'LineWidth', 2);
    end

    % Plot overall fixed effect line
    x_range_all = linspace(min(wm_acc), max(wm_acc), 100);
    y_fit_all = overall_intercept + overall_wm_slope * x_range_all + overall_session_slope * mean(session_number) + overall_interaction * x_range_all * mean(session_number);
    plot(x_range_all, y_fit_all, 'k--', 'LineWidth', 2);

    % Finalize plot
    xlabel('WM Accuracy', 'FontSize', 12);
    ylabel('Anticipated Swapping', 'FontSize', 12);
    title('WM Accuracy vs Anticipated Swapping', 'FontSize', 12);

    % Add legend for subjects
    num_subjects = length(unique_subjects);
    legendHandles = gobjects(1, num_subjects + 1);
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, [ markers{s} '-'], 'MarkerFaceColor', colors(s,:), 'MarkerEdgeColor', 'none','Color', colors(s,:), 'MarkerSize', 8, 'LineWidth', 1);
    end
    legendHandles(end) = plot(NaN, NaN, 'k--', 'LineWidth', 2);
    legend_obj = legend([legendHandles], [cellfun(@(x) ['Subject ', x(1)], unique_subjects, 'UniformOutput', false); {'Overall Trend'}], ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Box = 'off';

    % Adjust axis limits
    xlim([min(wm_acc) - 0.05, max(wm_acc) + 0.05]);
    ylim([min(anticipated_swapping) - 0.05, max(anticipated_swapping) + 0.05]);

    % Adjust font size and tick direction
    set(gca, 'FontSize', 12, 'TickDir', 'out');

    % Remove top and right spines
    box off;
    ax = gca;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.LineWidth = 1;

    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    % Display additional information
    fprintf('\nFixed Effects:\n');
    fprintf('Intercept: %.4f\n', overall_intercept);
    fprintf('WM Accuracy: %.4f (p=%.4f)\n', overall_wm_slope, overall_wm_p);
    fprintf('Session Number: %.4f (p=%.4f)\n', overall_session_slope, overall_session_p);
    fprintf('Interaction (WM Accuracy * Session Number): %.4f (p=%.4f)\n', overall_interaction, overall_interaction_p);

end

function WM_Acc_WWW_ER_MixedModel(WM_WWW_Interaction)

    % Prepare data
    subjects = WM_WWW_Interaction.Subject;
    wm_acc = WM_WWW_Interaction.OverallAccuracy;
    WWW_Accuracy = WM_WWW_Interaction.BlockLevelErrorRate;
    session_number = WM_WWW_Interaction.SessionNumber;

    % Remove any rows with NaN values
    valid_idx = ~isnan(wm_acc) & ~isnan(WWW_Accuracy) & ~isnan(session_number);
    subjects = subjects(valid_idx);
    wm_acc = wm_acc(valid_idx);
    WWW_Accuracy = WWW_Accuracy(valid_idx);
    session_number = session_number(valid_idx);

    % Create table for mixed-effects model
    tbl = table(subjects, wm_acc, WWW_Accuracy, session_number, 'VariableNames', {'Subject', 'WM_Accuracy', 'WWW_Accuracy', 'SessionNumber'});

    % Fit linear mixed-effects model with session number as a fixed effect and random effect
    lme = fitlme(tbl, 'WWW_Accuracy ~ WM_Accuracy + SessionNumber + WM_Accuracy * SessionNumber + (SessionNumber|Subject)');

    % Display model results
    disp(lme);

    % Extract fixed effects and their p-values
    [beta, names, stats] = fixedEffects(lme);
    overall_intercept = beta(1);
    overall_wm_slope = beta(2);  % The slope for WM_Accuracy
    overall_session_slope = beta(3);  % The slope for SessionNumber
    overall_wm_p = stats.pValue(2);
    overall_session_p = stats.pValue(3);

    % Calculate individual correlations and slopes
    unique_subjects = unique(subjects);
    individual_correlations = zeros(length(unique_subjects), 1);
    individual_wm_slopes = zeros(length(unique_subjects), 1);
    individual_session_slopes = zeros(length(unique_subjects), 1);
    individual_p_values = zeros(length(unique_subjects), 1);
    individual_intercepts = zeros(length(unique_subjects), 1);

    fprintf('\nIndividual Subject Correlations and Slopes:\n');
    fprintf('--------------------------------------------\n');
    fprintf('Subject\tCorrelation\tWM Slope\tSession Slope\tp-value\n');

    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_www_acc = WWW_Accuracy(subj_idx);
        subj_session = session_number(subj_idx);

        % Calculate partial correlation (WM_Accuracy and WWW_Accuracy, controlling for SessionNumber)
        [r, p] = partialcorr(subj_wm_acc, subj_www_acc, subj_session);
        individual_correlations(i) = r;
        individual_p_values(i) = p;

        % Calculate slopes and intercept using multiple regression
        X = [ones(size(subj_wm_acc)), subj_wm_acc, subj_session];
        coeffs = X \ subj_www_acc;
        individual_intercepts(i) = coeffs(1);
        individual_wm_slopes(i) = coeffs(2);
        individual_session_slopes(i) = coeffs(3);

        fprintf('%s\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', subj, r, coeffs(2), coeffs(3), p);
    end

    % Create scatter plot
    canvas_size = [100, 100, 300, 250];
    fig1 = figure('Position', canvas_size);
    set(gcf, 'renderer', 'Painters');
    hold on;

    % Define colors for each subject
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;

    % Define markers for each subject
    markers = {'o', 's', 'd', '^', 'v'};

    % Plot data points and individual correlation lines
    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_www_acc = WWW_Accuracy(subj_idx);

        % Plot scatter for this subject
        scatter(subj_wm_acc, subj_www_acc, 50, 'filled', 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'none', 'Marker', markers{i}, 'LineWidth', 1.5);
        
        % Plot individual correlation line
        x_range = linspace(min(subj_wm_acc), max(subj_wm_acc), 100);
        y_fit = individual_intercepts(i) + individual_wm_slopes(i) * x_range + ...
                individual_session_slopes(i) * mean(session_number(subj_idx));
        plot(x_range, y_fit, '-', 'Color', colors(i,:), 'LineWidth', 1);
    end

    % Plot overall fixed effect line for WM_Accuracy (assuming mean session number)
    x_range_all = linspace(min(wm_acc), max(wm_acc), 100);
    y_fit_all = overall_intercept + overall_wm_slope * x_range_all + overall_session_slope * mean(session_number);
    plot(x_range_all, y_fit_all, 'k--', 'LineWidth', 2);

    % Finalize plot
    xlabel('WM Accuracy', 'FontSize', 12);
    ylabel('WWW Error Rate', 'FontSize', 12);
    title('WM Accuracy vs WWW Error Rate (Controlling for Session)', 'FontSize', 12);
    
    % Add legend for subjects
    num_subjects = length(unique_subjects);
    legendHandles = gobjects(1, num_subjects + 1);
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, [markers{s}, '-'], 'MarkerFaceColor', colors(s,:), 'MarkerEdgeColor', 'none', 'Color', colors(s,:), 'MarkerSize', 8, 'LineWidth', 1);
    end
    legendHandles(end) = plot(NaN, NaN, 'k--', 'LineWidth', 2);
    legend_obj = legend(legendHandles, [cellfun(@(x) ['Subject ', x(1)], unique_subjects, 'UniformOutput', false); {'Overall Trend'}], ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Box = 'off';

    % Adjust axis limits
    xlim([min(wm_acc) - 0.05, max(wm_acc) + 0.05]);
    ylim([min(WWW_Accuracy) - 0.05, max(WWW_Accuracy) + 0.05]);

    % Adjust font size and tick direction
    set(gca, 'FontSize', 12, 'TickDir', 'out');

    % Remove top and right spines
    box off;
    ax = gca;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.LineWidth = 1;

    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    % Add overall effect information to the plot
%     annotation('textbox', [0.15, 0.95, 0.7, 0.05], ...
%                'String', sprintf('WM effect: %.4f (p=%.4f), Session effect: %.4f (p=%.4f)', ...
%                overall_wm_slope, overall_wm_p, overall_session_slope, overall_session_p), ...
%                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

function WM_Acc_WWW_ACC_MixedModel(WM_WWW_Interaction)

    % Prepare data
    subjects = WM_WWW_Interaction.Subject;
    wm_acc = WM_WWW_Interaction.OverallAccuracy;
    WWW_Accuracy = WM_WWW_Interaction.AvgBlockLevelAccuracy;
    session_number = WM_WWW_Interaction.SessionNumber;

    % Remove any rows with NaN values
    valid_idx = ~isnan(wm_acc) & ~isnan(WWW_Accuracy) & ~isnan(session_number);
    subjects = subjects(valid_idx);
    wm_acc = wm_acc(valid_idx);
    WWW_Accuracy = WWW_Accuracy(valid_idx);
    session_number = session_number(valid_idx);

    % Create table for mixed-effects model
    tbl = table(subjects, wm_acc, WWW_Accuracy, session_number, 'VariableNames', {'Subject', 'WM_Accuracy', 'WWW_Accuracy', 'SessionNumber'});

    % Fit linear mixed-effects model with session number as a fixed effect and random effect
    lme = fitlme(tbl, 'WWW_Accuracy ~ WM_Accuracy + SessionNumber + WM_Accuracy * SessionNumber + (SessionNumber|Subject)');

    % Display model results
    disp(lme);

    % Extract fixed effects and their p-values
    [beta, names, stats] = fixedEffects(lme);
    overall_intercept = beta(1);
    overall_wm_slope = beta(2);  % The slope for WM_Accuracy
    overall_session_slope = beta(3);  % The slope for SessionNumber
    overall_wm_p = stats.pValue(2);
    overall_session_p = stats.pValue(3);

    % Calculate individual correlations and slopes
    unique_subjects = unique(subjects);
    individual_correlations = zeros(length(unique_subjects), 1);
    individual_wm_slopes = zeros(length(unique_subjects), 1);
    individual_session_slopes = zeros(length(unique_subjects), 1);
    individual_p_values = zeros(length(unique_subjects), 1);
    individual_intercepts = zeros(length(unique_subjects), 1);

    fprintf('\nIndividual Subject Correlations and Slopes:\n');
    fprintf('--------------------------------------------\n');
    fprintf('Subject\tCorrelation\tWM Slope\tSession Slope\tp-value\n');

    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_www_acc = WWW_Accuracy(subj_idx);
        subj_session = session_number(subj_idx);

        % Calculate partial correlation (WM_Accuracy and WWW_Accuracy, controlling for SessionNumber)
        [r, p] = partialcorr(subj_wm_acc, subj_www_acc, subj_session);
        individual_correlations(i) = r;
        individual_p_values(i) = p;

        % Calculate slopes and intercept using multiple regression
        X = [ones(size(subj_wm_acc)), subj_wm_acc, subj_session];
        coeffs = X \ subj_www_acc;
        individual_intercepts(i) = coeffs(1);
        individual_wm_slopes(i) = coeffs(2);
        individual_session_slopes(i) = coeffs(3);

        fprintf('%s\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', subj, r, coeffs(2), coeffs(3), p);
    end

    % Create scatter plot
    canvas_size = [100, 100, 300, 250];
    fig1 = figure('Position', canvas_size);
    set(gcf, 'renderer', 'Painters');
    hold on;

    % Define colors for each subject
    colors = [162, 34, 91;
              218, 114, 46;
              152, 173, 54;
              79, 174, 226;
              37, 90, 164] / 255;

    % Define markers for each subject
    markers = {'o', 's', 'd', '^', 'v'};

    % Plot data points and individual correlation lines
    for i = 1:length(unique_subjects)
        subj = unique_subjects{i};
        subj_idx = strcmp(subjects, subj);
        subj_wm_acc = wm_acc(subj_idx);
        subj_www_acc = WWW_Accuracy(subj_idx);

        % Plot scatter for this subject
        scatter(subj_wm_acc, subj_www_acc, 50, 'filled', 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'none', 'Marker', markers{i}, 'LineWidth', 1.5);
        
        % Plot individual correlation line
        x_range = linspace(min(subj_wm_acc), max(subj_wm_acc), 100);
        y_fit = individual_intercepts(i) + individual_wm_slopes(i) * x_range + ...
                individual_session_slopes(i) * mean(session_number(subj_idx));
        plot(x_range, y_fit, '-', 'Color', colors(i,:), 'LineWidth', 1);
    end

    % Plot overall fixed effect line for WM_Accuracy (assuming mean session number)
    x_range_all = linspace(min(wm_acc), max(wm_acc), 100);
    y_fit_all = overall_intercept + overall_wm_slope * x_range_all + overall_session_slope * mean(session_number);
    plot(x_range_all, y_fit_all, 'k--', 'LineWidth', 2);

    % Finalize plot
    xlabel('WM Accuracy', 'FontSize', 12);
    ylabel('WWW Accuracy', 'FontSize', 12);
    title('WM Accuracy vs WWW Accuracy (Controlling for Session)', 'FontSize', 12);
    
    % Add legend for subjects
    num_subjects = length(unique_subjects);
    legendHandles = gobjects(1, num_subjects + 1);
    for s = 1:num_subjects
        legendHandles(s) = plot(NaN, NaN, [markers{s}, '-'], 'MarkerFaceColor', colors(s,:), 'MarkerEdgeColor', 'none', 'Color', colors(s,:), 'MarkerSize', 8, 'LineWidth', 1);
    end
    legendHandles(end) = plot(NaN, NaN, 'k--', 'LineWidth', 2);
    legend_obj = legend(legendHandles, [cellfun(@(x) ['Subject ', x(1)], unique_subjects, 'UniformOutput', false); {'Overall Trend'}], ...
        'Location', 'eastoutside', 'FontSize', 12);
    legend_obj.EdgeColor = 'none';
    legend_obj.Box = 'off';

    % Adjust axis limits
    xlim([min(wm_acc) - 0.05, max(wm_acc) + 0.05]);
    ylim([min(WWW_Accuracy) - 0.05, max(WWW_Accuracy) + 0.05]);

    % Adjust font size and tick direction
    set(gca, 'FontSize', 12, 'TickDir', 'out');

    % Remove top and right spines
    box off;
    ax = gca;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.LineWidth = 1;

    hold off;

    % Adjust figure size to accommodate legend
    legend_pos = get(legend_obj, 'Position');
    legend_width = legend_pos(3);
    new_fig_width = canvas_size(3) + legend_width * canvas_size(3);
    new_fig_size = [canvas_size(1), canvas_size(2), new_fig_width, canvas_size(4)];
    set(fig1, 'Position', new_fig_size);

    % Add overall effect information to the plot
%     annotation('textbox', [0.15, 0.95, 0.7, 0.05], ...
%                'String', sprintf('WM effect: %.4f (p=%.4f), Session effect: %.4f (p=%.4f)', ...
%                overall_wm_slope, overall_wm_p, overall_session_slope, overall_session_p), ...
%                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
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

function plotTransitionProbabilitiesWithAccuracyFilter(touch_table, figure_require)
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
            		 
					if figure_require == 1
            			if block_level_error_rate > 6.6% or > 6.6
							includeBlock(i) = 1;
                			includeBlock(i + 1) = 1;  % Exclude next block if it exists
						end
					elseif figure_require == 2
						if block_level_error_rate < 5.3% or > 6.6
							includeBlock(i) = 1;
                			includeBlock(i + 1) = 1;  % Exclude next block if it exists
						end
					else
						
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


function plotNextErrorFrequency(touch_table)
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
    [~, p_value34] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 4), 'Vartype', 'unequal');
	[~, p_value35] = ttest2(allProbsCombined(:, 3), allProbsCombined(:, 5), 'Vartype', 'unequal');
	[~, p_value45] = ttest2(allProbsCombined(:, 4), allProbsCombined(:, 5), 'Vartype', 'unequal');
	
	n_compare = 3;
	p_value34 = p_value34 * n_compare;
	p_value35 = p_value35 * n_compare;
	p_value45 = p_value45 * n_compare;

	fprintf('P-value for Ordinal Position 3 vs 4 comparison: %.10f\n', p_value34);
	fprintf('P-value for Ordinal Position 3 vs 5 comparison: %.10f\n', p_value35);
	fprintf('P-value for Ordinal Position 4 vs 5 comparison: %.10f\n', p_value45);

    % Add significance star if p < 0.05 and print p-value
    if p_value34 < 0.05
        maxY = max(allProbsCombined,[], 'all');
        plot([3, 4], [maxY*1.1, maxY*1.1], 'k-');
		if p_value34 < 0.001
        text(3.5, maxY*1.15, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif  p_value34 < 0.01
			text(3.5, maxY*1.15, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
		elseif p_value34 < 0.05
			text(3.5, maxY*1.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
		end
    end
    
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
