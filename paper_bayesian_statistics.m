analyze_touch_table_HB_hierarchical(touch_table)


function analyze_touch_table_HB_hierarchical(touch_table)
% Main entry point: builds the data structure, specifies a Beta prior on θ,
% then calls the hierarchical sampler which uses MH for θ.
   %% Define outcome set and constants
    % We work only with outcomes after the correct A touch.
    % Outcome categories: B, C, D, E, Distractor.
    outcomes = {'B', 'C', 'D', 'E', 'Distractor'};  % 5 outcomes
    % In the raw data the TouchObjectCorrectOrdinalPosition is coded as:
    %   1: A, 2: B, 3: C, 4: D, 5: E, 0: Distractor -> recoded as 6.
    % For swapped trials, we exchange positions 2 and 4.
    
    %% Build session-wise data structure for the hierarchical model
    subjects = unique(touch_table.Subject);
    nSubjects = numel(subjects);
    
    % We will store data per subject, per session in the structure "data"
    data.subject = [];
    % Also, for later regression we store observation-level info:
    obs_Subject = [];  % subject id (categorical)
    obs_SessionOrder = [];  % session number (e.g., as provided by SessionNumber)
    obs_delta = []; % will hold delta = theta_D/theta_B (posterior estimate)
    accuracy_threshold = 0.8;

	
    for si = 1:nSubjects
        subj = subjects{si};
        subj_data = touch_table(strcmp(touch_table.Subject, subj), :);
        sessions = unique(subj_data.SessionNumber);
        % Sort sessions if necessary
        sessions = sort(sessions);
        % Initialize subject-level session structure array
        subjSessions = [];
        
        for s_idx = 1:numel(sessions)
            sessionNum = sessions(s_idx);
            session_data = subj_data(subj_data.SessionNumber == sessionNum, :);
            blocks = unique(session_data.BlockNumber);

			 % Calculate learning speed for the session
			 session_completed_trials = [];
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


            
            % Initialize count vector (length 6: positions 1 to 6)
            % We will later drop index 1 because outcome A is not included.
            session_counts_full = zeros(1,6);
            total_valid_trials = 0;
            
            % Loop over blocks and then trials
            for b = 1:numel(blocks)
                block = blocks(b);
                block_data = session_data(session_data.BlockNumber == block, :);
                trials = unique(block_data.TrialNumberInBlock);
                for t = 1:min(numel(trials), learning_speed)
                    trial = trials(t);
                    trial_data = block_data(block_data.TrialNumberInBlock == trial, :);
                    % Process only swapped trials (experimental phase where sequence changed)
                    if trial_data.isSwapped(1) == 1
                        % Look for the first correct A touch within the trial
                        for touchIdx = 1:(height(trial_data)-1)
                            if (trial_data.CurrentState(touchIdx) == 1) && ...
                               (strcmpi(trial_data.TouchCategory{touchIdx}, 'correct') || ...
                                strcmpi(trial_data.TouchCategory{touchIdx}, 'correctSelection'))
                            
                                nextIdx = touchIdx + 1;
                                if nextIdx <= height(trial_data)
                                    pos = trial_data.TouchObjectCorrectOrdinalPosition(nextIdx);
                                    % Recode: if pos==0 then set to 6 (Distractor)
                                    if pos == 0
                                        pos = 6;
                                    end
                                    % For swapped trials: exchange positions 2 and 4
                                    if pos == 2
                                        pos = 4;
                                    elseif pos == 4
                                        pos = 2;
                                    end
                                    % Increment count at position "pos"
                                    % (the count vector indices run 1:6)
                                    session_counts_full(pos) = session_counts_full(pos) + 1;
                                    total_valid_trials = total_valid_trials + 1;
                                end
                                % Only record the first correct A transition in the trial
                                break;
                            end
                        end
                    end % if isSwapped
                end % trials loop
            end % blocks loop
            
            % Only consider session if at least one valid trial was observed.
            if total_valid_trials > 0
                % Remove the outcome "A" (position 1) since we need outcomes for B, C, D, E, Distractor.
                session_counts = session_counts_full(2:6);  % indices correspond to [B, C, D, E, Distractor]
                % Save session-level info in the data structure.
                sess_struct.counts = session_counts;
                sess_struct.total = total_valid_trials;
                sess_struct.session_order = sessionNum;  % or use s_idx for sequential order if preferred
                
                subjSessions = [subjSessions, sess_struct];
            end
        end % sessions loop
        
        % Save subject-level data only if there is at least one session.
        if ~isempty(subjSessions)
            data.subject(si).id = subj;
            % Each subject holds a structure array "sessions"
            data.subject(si).sessions = subjSessions;
        end
    end % subjects loop
    

    %% 2) Prior & sampling settings
    a_beta     = 2;    % Beta prior α
    b_beta     = 5;    % Beta prior β
    numSamples = 200000;
    burnIn     = 100000;
    
    model = struct(...
      'likelihood',   @multinomial_logLikelihood, ...
      'prior_params', struct('a',a_beta,'b',b_beta), ...
      'hierarchy',    1  ...  % no hyper‐alpha any more
    );
	% 
	% % Modified model specification
	% alpha = [2 2 2 2 2]; % Dirichlet prior parameters (symmetric)
	% model = struct(...
    % 	'likelihood', @multinomial_logLikelihood, ...
    % 	'prior_params', struct('alpha', alpha), ...
    % 	'hierarchy', 1 ...
	% );
	

    %% 3) Run MCMC
    [posterior, diagnostics] = hbm_sampler_aggregated(data, numSamples, burnIn, model);
    
	% At the end of analyze_touch_table_HB_hierarchical function:
	visualize_mcmc_diagnostics(diagnostics, posterior);

    % %% 4) Visualize
    visualize_aggregated_posterior(posterior, diagnostics);

	
end



function [posterior, diagnostics] = hbm_sampler_aggregated(data, numSamples, burnIn, model)
% sample a Multinomial–Beta posterior.
rng(5)
nSubj = numel(data.subject);
K     = numel(data.subject(1).sessions(1).counts);
posterior.burnIn = burnIn;

for i = 1:nSubj
    % aggregate counts
    totalCounts = zeros(1,K);
    for s = data.subject(i).sessions
        totalCounts = totalCounts + s.counts;
    end
    
    % % define log‐posterior for θ
    logpost_theta = @(theta) ...
      multinomial_logLikelihood(theta, totalCounts) + ...
      sum((model.prior_params.a-1)*log(theta) + (model.prior_params.b-1)*log(1-theta));
	
    
    % run MH chain for θ
    thetaChain = zeros(numSamples, K);
    % initialize at MLE
    thetaChain(1,:) = totalCounts / sum(totalCounts);
    
    for t = 2:numSamples
        % Gaussian RW proposal in each component
        prop = thetaChain(t-1,:) + 0.05*randn(1,K);
        % project back to (ε,1–ε)
        prop = max(min(prop,1-1e-6),1e-6);
        % renormalize to sum=1 (simplex constraint)
        prop = prop / sum(prop);
        
        logA = logpost_theta(prop) - logpost_theta(thetaChain(t-1,:));
        if log(rand) < logA
            thetaChain(t,:) = prop;
        else
            thetaChain(t,:) = thetaChain(t-1,:);
        end
    end
    
    % store results
    samples = thetaChain(burnIn+1:end,:);
    posterior.subject(i).id           = data.subject(i).id;
    posterior.subject(i).thetaSamples = samples;
    posterior.subject(i).theta        = mean(samples,1);
    diagnostics.theta_chain{i}        = thetaChain;
end

end

function ll = multinomial_logLikelihood(theta, counts)
% MULTINOMIAL_LOGLIKELIHOOD  log P(counts | theta) under Multinomial
%   theta  – 1×K probability vector (sum(theta)=1, theta>0)
%   counts – 1×K nonnegative integer counts
%
% Returns -Inf if any theta<=0 to enforce support.

    if any(theta <= 0)
        ll = -Inf;
        return;
    end

    % Pure multinomial log-likelihood (ignoring the constant term)
    ll = sum(counts .* log(theta));
end



function visualize_mcmc_diagnostics(diagnostics, posterior)
    % VISUALIZE_MCMC_DIAGNOSTICS Visualize MCMC convergence diagnostics
    %
    % Inputs:
    %   diagnostics - struct containing theta_chain for each subject
    %   posterior   - struct from hbm_sampler_aggregated containing burnIn
    %   subjLabels  - cell array of subject labels
    
    subjLabels  = {'Subject B','Subject J','Subject K','Subject S'};
    
    % Define outcome labels
    outcomes = {'B', 'C', 'D', 'E', 'Distractor'};
    nOutcomes = length(outcomes);
    
    % Define colors for different parameters
   	colors = [
    	182, 118, 108;  % B
    	231, 188, 198;  % C
    	138, 140, 191;  % D
    	78,  101, 155;  % E
    	95,   95,  95   % Distractor
	] / 255;
    
    burnIn = posterior.burnIn;
    nSubjects = numel(diagnostics.theta_chain);
    
        
    % 2. Calculate Gelman-Rubin diagnostic (R-hat)
    fprintf('\nGelman-Rubin Diagnostic (R-hat):\n');
    fprintf('-------------------------------\n');
    fprintf('Values close to 1 indicate convergence. R-hat > 1.1 suggests lack of convergence.\n\n');
    % Initialize storage for R-hat values across all subjects
all_rhat = zeros(nSubjects, nOutcomes);

% For each subject
for i = 1:nSubjects
    fprintf('%s:\n', subjLabels{i});
    
    % Get post burn-in samples
    chain = diagnostics.theta_chain{i}(burnIn+1:end, :);
    
    % Split the chain into multiple segments for Gelman-Rubin calculation
    chainLength = size(chain, 1);
    numSegments = min(4, max(2, floor(chainLength / 1000))); 
    segmentLength = floor(chainLength / numSegments);
    
    if segmentLength > 50  % Only compute if we have reasonable segment lengths
        
        % Initialize storage for R-hat values
        rhat = zeros(1, nOutcomes);
        
        for j = 1:nOutcomes
            % Create matrix of chains by splitting the single chain
            chains = zeros(segmentLength, numSegments);
            
            for seg = 1:numSegments
                startIdx = (seg - 1) * segmentLength + 1;
                endIdx = startIdx + segmentLength - 1;
                chains(:, seg) = chain(startIdx:endIdx, j);
            end
            
            % Compute Gelman-Rubin statistic
            rhat(j) = compute_gelman_rubin(chains);
            all_rhat(i, j) = rhat(j);
            
            % Display with interpretations
            if rhat(j) > 1.1
                fprintf('  θ_%s: R-hat = %.4f (Potential convergence issue)\n', outcomes{j}, rhat(j));
            elseif rhat(j) > 1.05
                fprintf('  θ_%s: R-hat = %.4f (Marginal convergence)\n', outcomes{j}, rhat(j));
            else
                fprintf('  θ_%s: R-hat = %.4f (Good convergence)\n', outcomes{j}, rhat(j));
            end
        end
    else
        fprintf('  Insufficient samples to compute reliable Gelman-Rubin diagnostic\n');
        all_rhat(i, :) = NaN;
    end
    fprintf('\n');
end


% Create figure with specified size
figure('Position', [100, 100, 1000, 150]);

% Create 1x4 subplots
for i = 1:nSubjects
    subplot(1, 4, i);
    
    % Create bar graph with the specified colors
    b = bar(1:nOutcomes, all_rhat(i,:));
    
    % Assign colors to each bar
    b.FaceColor = 'flat';
    for j = 1:nOutcomes
        b.CData(j,:) = colors(j,:);
    end
    
    % Add reference lines
    yline(1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'Perfect convergence');
    yline(1.1, 'r--', 'LineWidth', 1, 'DisplayName', 'Convergence threshold');
    
    % Labels and formatting
    xlabel('Outcome', 'FontSize', 12);
    if i == 1
        ylabel('R-hat', 'FontSize', 12);
    end
    title(['Subject'], 'FontSize', 14);
    set(gca, 'XTick', 1:nOutcomes, 'XTickLabel', outcomes, 'FontSize', 10);
    ylim([0.9, max(1.15, max(all_rhat(:)))]);
    grid on;
end

% Add main title
sgtitle('Gelman-Rubin Diagnostic (R-hat)', 'FontSize', 14);

% Create a legend for the entire figure
% legendLabels = [outcomes, {'Perfect convergence', 'Convergence threshold'}];
% legend('show', 'Location', 'eastoutside');

    
    % 3. Running mean plots
    figure('Name', 'MCMC Running Means', 'Units', 'pixels', 'Position', [300, 300, 1000, 800]);
    
    for i = 1:min(nSubjects, 4)
        for j = 1:nOutcomes
            subplot(nOutcomes, min(nSubjects, 4), (j-1)*min(nSubjects, 4) + i);
            
            % Get post burn-in samples
            chain = diagnostics.theta_chain{i}(burnIn+1:end, j);
            
            % Calculate running means
            running_means = cumsum(chain) ./ (1:length(chain))';
            
            % Plot running means
            plot(running_means, 'Color', colors(j,:), 'LineWidth', 1.5);
            
            if j == 1
                title(sprintf('%s', subjLabels{i}), 'FontSize', 12);
            end
            
            if i == 1
                ylabel(sprintf('Running Mean \\theta_{%s}', outcomes{j}), 'FontSize', 12);
            end
            
            if j == nOutcomes
                xlabel('Iteration (post burn-in)', 'FontSize', 12);
            end
            
            grid on;
        end
    end
    
    % 4. Compute and display diagnostic statistics
    disp('MCMC Diagnostics Summary:');
    disp('-------------------------');
    for i = 1:nSubjects
        fprintf('\n%s:\n', subjLabels{i});
        
        post_samples = diagnostics.theta_chain{i}(burnIn+1:end, :);
        
        for j = 1:nOutcomes
            % Calculate effective sample size (ESS) using autocorrelation
            chain = post_samples(:, j);
            [acf, ~] = xcorr(chain, min(50, length(chain)-1));
            
            % Compute ESS using initial positive sequence estimator
            pos_idx = find(acf(2:end) < 0, 1);
            if isempty(pos_idx)
                pos_idx = length(acf);
            end
            
            rho_sum = 2 * sum(acf(2:pos_idx));
            ess = length(chain) / (1 + rho_sum);
            
            fprintf('  θ_%s: Mean = %.4f, SD = %.4f', ...
                outcomes{j}, mean(chain), std(chain));
        end
    end
    
    % Save figures
    saveas(figure(2), 'mcmc_running_means.svg');
	saveas(figure(1), 'gelman_rubin_diagnostic.svg');
end

function rhat = compute_gelman_rubin(chains)
    % COMPUTE_GELMAN_RUBIN Calculate Gelman-Rubin statistic
    % chains: matrix where each column is a chain segment
    
    [n, m] = size(chains);  % n = chain length, m = number of chain segments
    
    % Calculate chain means
    chain_means = mean(chains, 1);
    
    % Calculate grand mean
    grand_mean = mean(chain_means);
    
    % Calculate between-chain variance B
    B = (n/(m-1)) * sum((chain_means - grand_mean).^2);
    
    % Calculate within-chain variance W
    % Use var() without the 'omitnan' flag since we shouldn't have NaNs
    W = mean(var(chains, 0, 1));
    
    % Handle numerical issues
    if W < eps
        if B < eps
            rhat = 1.0;  % Both variances are effectively zero
        else
            rhat = 10.0;  % Large between-chain variance, near-zero within-chain variance
        end
    else
        % Calculate variance estimate
        var_est = ((n-1)/n) * W + (1/n) * B;
        
        % Calculate potential scale reduction factor
        rhat = sqrt(var_est / W);
    end
    
    % Ensure rhat is not NaN or Inf
    if isnan(rhat) || isinf(rhat)
        rhat = 10.0;  % Set to a high value to indicate a problem
    end
end


function visualize_aggregated_posterior(posterior, diagnostics)
    % VISUALIZE_AGGREGATED_POSTERIOR  Wrapper to plot all summaries
    % Inputs:
    %   posterior   – struct from hbm_sampler_aggregated
    %   diagnostics – struct (unused here)

    % Define outcome labels
    outcomes    = {'B','C','D','E','Distractor'};
    % Define subject labels for indices 1–4
    subjLabels  = {'Subject B','Subject J','Subject K','Subject S'};

    % Collect theta samples per subject
    nSubjects = numel(posterior.subject);
    theta_samples = cell(nSubjects,1);
    for i = 1:nSubjects
        theta_samples{i} = posterior.subject(i).thetaSamples;
    end

    % Call each plotting routine
    plot_theta_distributions(theta_samples, outcomes, subjLabels);
    plot_hypothesis_testing(theta_samples, subjLabels);
    plot_model_comparison(theta_samples);
end
%% ------------------------------------------------------------------------
function plot_theta_distributions(theta_samples, outcomes, subjLabels)
    % Compute mean±std for B,C,D for each subject
    n = numel(theta_samples);
    meanVals = zeros(n,3);
    stdVals  = zeros(n,3);
    for i = 1:n
        % only columns 1–3 correspond to B,C,D
        meanVals(i,:) = mean(theta_samples{i}(:,1:3));
        stdVals(i,:)  = std (theta_samples{i}(:,1:3));
    end

    % Print summary lines
    fprintf('\nPosterior θ Summary:\n');
    for i = 1:n
        fprintf('%s: B = %.2f ± %.2f,  C = %.2f ± %.2f,  D = %.2f ± %.2f\n', ...
            subjLabels{i}, ...
            meanVals(i,1), stdVals(i,1), ...
            meanVals(i,2), stdVals(i,2), ...
            meanVals(i,3), stdVals(i,3));
    end
    fprintf('\n');

    % Now proceed with plotting
    colB = [190,  20,  32]/255;  % #be1420
    colC = [ 95,  95,  95]/255;  % #5f5f5f
    colD = [102, 154, 186]/255;  % #669aba
    cols = [colB; colC; colD];

    figure('Name','Posterior Distributions of θ', ...
           'Units','pixels','Position',[100,100,600,400]);
    for i = 1:min(n,4)
        ax = subplot(2,2,i); hold(ax,'on');
        for k = 1:3
            % Example using a fixed number of bins (e.g. 20 bins):
			histogram(theta_samples{i}(:,k), ...
    			'Normalization','probability', ...
    			'NumBins', 15, ...
    			'FaceColor', cols(k,:), ...
    			'EdgeColor', 'none', ...
    			'FaceAlpha', 0.5);
        	xline(meanVals(i,k), ...
              	'Color',cols(k,:), ...
              	'LineWidth',3, ...
              	'HandleVisibility','off');
        end
        title(subjLabels{i}, 'FontSize',12);
        xlabel('Posterior Distribution','FontSize',12);
        ylabel('Density','FontSize',12);
        % legend(outcomes{1:3}, 'FontSize',10,'Location','best');
        set(ax,'FontSize',12);
        hold(ax,'off');
    end

    saveas(gcf,'posterior_distributions.svg');
end


%% ------------------------------------------------------------------------
function plot_hypothesis_testing(theta_samples, subjLabels)
    n   = numel(theta_samples);
    pDB = zeros(n,1);
    pDC = zeros(n,1);
    sdDB = zeros(n,1);
    sdDC = zeros(n,1);

    % Define new colors
    colDB = [255, 127,   0]/255;   % #ff7f00
    colDC = [152,  78, 163]/255;   % #984ea3

    % Compute P(D>B), P(D>C) and their std-dev (across the Bernoulli indicator)
    for i = 1:n
        diffDB = theta_samples{i}(:,3) - theta_samples{i}(:,1);
        diffDC = theta_samples{i}(:,3) - theta_samples{i}(:,2);
        % probability estimates
        pDB(i) = mean(diffDB>0);
        pDC(i) = mean(diffDC>0);
        % std of the 0/1 indicator
        sdDB(i) = std(double(diffDB>0));
        sdDC(i) = std(double(diffDC>0));
    end

    % Print summary lines
    fprintf('\nHypothesis Summary:\n');
    for i = 1:n
        fprintf('  %s: P(D > B) = %.2f ± %.2f,  P(D > C) = %.2f ± %.2f\n', ...
                subjLabels{i}, pDB(i), sdDB(i), pDC(i), sdDC(i));
    end
    fprintf('\n');

    % --------------------------------------------------------------------
    % Now plot the distributions of differences
    figure('Name','Hypothesis Testing', ...
           'Units','pixels','Position',[200,200,600,400]);
    for i = 1:min(n,4)
        diffDB = theta_samples{i}(:,3) - theta_samples{i}(:,1);
        diffDC = theta_samples{i}(:,3) - theta_samples{i}(:,2);

        ax = subplot(2,2,i); hold(ax,'on');
        histogram(diffDB, ...
                  'Normalization','probability', ...
				  'NumBins', 15, ...
                  'FaceColor',colDB, ...
				  'EdgeColor','none',...
                  'FaceAlpha',0.5);
        histogram(diffDC, ...
                  'Normalization','probability', ...
				  'NumBins', 15, ...
                  'FaceColor',colDC, ...
				  'EdgeColor','none',...
                  'FaceAlpha',0.5);

        xline(0,      'k--','LineWidth',1);
        xline(mean(diffDB), 'Color',colDB,'LineWidth',2,'HandleVisibility','off');
        xline(mean(diffDC), 'Color',colDC,'LineWidth',2,'HandleVisibility','off');

        title(sprintf('%s', subjLabels{i}), 'FontSize',12);
        xlabel('\theta Difference','FontSize',12);
        ylabel('Density','FontSize',12);
        % legend({'diffDB','diffDC'},'FontSize',10,'Location','best');
        set(ax,'FontSize',12);
        hold(ax,'off');
    end
    saveas(gcf,'hypothesis_testing.svg');
end

%% ------------------------------------------------------------------------
function plot_model_comparison(theta_samples)
    n = numel(theta_samples);
    models = {'Positional','Interference','Chaining','Mixed'};
    probs = zeros(n,4);

    for i = 1:n
        samp = theta_samples{i};
        assign = arrayfun(@classify_model, ...
                          mat2cell(samp,ones(size(samp,1),1),5));
        for m = 1:4
            probs(i,m) = mean(assign==m);
        end
    end
	% Add data reporting for each subject
    fprintf('\n--- Model Classification Report ---\n\n');
    for i = 1:n
        fprintf('Subject %d: position coding: %.1f%%, memory interference: %.1f%%, serial inference: %.1f%%, and mixed strategy: %.1f%%\n', ...
                i, probs(i,1)*100, probs(i,2)*100, probs(i,3)*100, probs(i,4)*100);
    end
    fprintf('\n');

    figure('Name','Model Comparison','Units','pixels', ...
           'Position',[300,300,600,400]);
    % Top: per‐subject stacked
    ax1 = subplot(2,1,1);
    bar(probs,'stacked');
    xlabel('Subject','FontSize',12);
    ylabel('P','FontSize',12);
    legend(models,'FontSize',10,'Location','eastoutside');
    title('Subject‐level Model Posterior','FontSize',12);
    set(ax1,'FontSize',12);

    % Bottom: overall
    ax2 = subplot(2,1,2);
    overall = mean(probs,1);
    bar(overall);
    set(ax2,'XTickLabel',models,'FontSize',12);
    ylabel('P','FontSize',12);
    title('Overall Model Posterior','FontSize',12);
    ylim([0 1]);
    % no grid on either
    saveas(gcf,'model_comparison.svg');
end

%% ------------------------------------------------------------------------
function model = classify_model(theta)
    % CLASSIFY_MODEL  Returns 1–4 based on θ_B,θ_C,θ_D
    theta = theta{1};  % from mat2cell
    tB = theta(1);
    tC = theta(2);
    tD = theta(3);
    if tB <= tD && tC <= tD
        model = 1;
    elseif tB >= tD && tB >= tC && tC >= tD
        model = 2;
    elseif tC >= tB && tC >= tD
        model = 3;
    else
        model = 4;
    end
end
