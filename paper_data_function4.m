
analyze_touch_table_HB_hierarchical(touch_table)

function analyze_touch_table_HB_hierarchical(touch_table)
% analyze_touch_table_HB_hierarchical performs a hierarchical Bayesian analysis
% for sequence learning data stored in touch_table. For each subject and session:
%
% -- The observed outcome (the touch immediately following a correct A)
%    is assumed to be one of {B, C, D, E, Distractor}. The outcome is obtained
%    from the TouchObjectCorrectOrdinalPosition field (with recoding: 0->Distractor,
%    swap positions 2 and 4 for swapped trials).
%
% -- Likelihood: y_ij ~ Multinomial(n_ij, theta_ij) where theta_ij =
%    [theta_ij^B, theta_ij^C, theta_ij^D, theta_ij^E, theta_ij^Distractor].
%
% -- Hierarchical priors: theta_ij ~ Dirichlet(alpha_i + 1/2) with
%    alpha_i ~ Gamma(mu_alpha, sigma_alpha^2).
%
% -- Contrast: delta_ij = theta_ij^D / theta_ij^B.
%    A population regression model is assumed:
%         log(delta_ij) = beta0 + beta1 * X_ij + u_i,
%    where X_ij is the session order (or trial number within session)
%    and u_i ~ N(0, sigma_u^2) is a subject-specific random intercept.
%
% Input:
%   touch_table - a table with the columns:
%       'Subject', 'SessionNumber', 'BlockNumber', 'TrialNumberInBlock',
%       'isSwapped', 'isRepetition', 'isRewarded', 'CurrentState',
%       'TouchCategory', 'TouchObjectCorrectOrdinalPosition'
%
% Note: This function assumes that the hierarchical Bayesian sampler
% (e.g., hbm_sampler from the hierarchical-matlab toolbox) is available
% on the MATLAB path. If not, the observed proportions will be used.
%

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
    
    %% Hierarchical Bayesian Inference
	
	mu_alpha = 2;      % (example value; adjust as needed)
    sigma_alpha = 2;   % (example value; adjust as needed)
    
    model = struct(...
        'likelihood', @dirichlet_multinomial_logLikelihood, ...
        'prior_params', struct('mu', mu_alpha, 'sigma', sigma_alpha), ...
        'hierarchy', 2);
    
    numSamples = 10000;
    burnIn = 5000;
    
    % Call the custom hierarchical Bayesian sampler.
    [posterior, diagnostics] = hbm_sampler_aggregated(data, numSamples, burnIn, model);
    

	visualize_aggregated_posterior(posterior, diagnostics);

end








function [posterior, diagnostics] = hbm_sampler_aggregated(data, numSamples, burnIn, model)
% HBM_SAMPLER_AGGREGATED  Hierarchical sampler using only subject‐total counts
%
% Inputs:
%   data           – struct with fields .subject(i).id and .subject(i).sessions(j).counts
%   numSamples     – total MCMC draws per chain
%   burnIn         – number of initial draws to discard
%   model.prior_params.mu, .sigma – Gamma(shape,scale) hyperparameters
%
% Outputs:
%   posterior      – struct with fields
%       .burnIn
%       .subject(i).id
%       .subject(i).alpha     = posterior mean of α_i
%       .subject(i).theta     = posterior mean of θ_i
%       .subject(i).thetaSamples (numSamples–burnIn × K)
%   diagnostics   – struct with .alpha_samples{i}

nSubjects = numel(data.subject);
K = numel(data.subject(1).sessions(1).counts);

posterior.burnIn = burnIn;
diagnostics.alpha_samples = cell(nSubjects,1);

for i = 1:nSubjects
    % 1) aggregate counts across sessions
    sessions = data.subject(i).sessions;
    totalCounts = zeros(1,K);
    for j = 1:numel(sessions)
        totalCounts = totalCounts + sessions(j).counts;
    end

    % 2) define log-posterior for alpha_i using only totalCounts
    logpost_agg = @(alpha) logpost_alpha_agg(alpha, totalCounts, ...
                                             model.prior_params.mu, ...
                                             model.prior_params.sigma);

    % 3) run Metropolis–Hastings on log(alpha)
    logalpha0 = log(1);
    propStd   = 0.3;
    logChain  = mhsample(logalpha0, numSamples, ...
                         'logpdf', @(la) logpost_agg(exp(la)) + la, ...
                         'proprnd', @(la) normrnd(la, propStd), ...
                         'symmetric', true);

    alpha_samples = exp(logChain);
    diagnostics.alpha_samples{i} = alpha_samples;

    % posterior mean alpha
    alpha_hat = mean(alpha_samples(burnIn+1:end));

    % 4) sample θ_i from Dirichlet(totalCounts + (α+0.5))
    a_vec = totalCounts + (alpha_samples(burnIn+1:end));
    nPost = numSamples - burnIn;
    thetaSamples = zeros(nPost, K);
    for t = 1:nPost
        thetaSamples(t,:) = dirichlet_sample(a_vec(t,:));
    end

    % store results
    posterior.subject(i).id            = data.subject(i).id;
    posterior.subject(i).alpha         = alpha_hat;
    posterior.subject(i).thetaSamples  = thetaSamples;
    posterior.subject(i).theta         = mean(thetaSamples,1);
end

end

%% Log‐posterior for aggregated counts
function lp = logpost_alpha_agg(alpha, counts, mu, sigma)
    % Gamma prior
    if alpha <= 0
        lp = -Inf; return
    end
    lp = (mu-1)*log(alpha) - alpha/sigma - mu*log(sigma) - gammaln(mu);
    % one Dirichlet–multinomial term
    a = alpha + 0.5;
    K = numel(counts);
    N = sum(counts);
    lp = lp + ...
         gammaln(K*a) - gammaln(N + K*a) + ...
         sum(gammaln(counts + a) - gammaln(a) - gammaln(counts+1));
end

%% Dirichlet sampler (as before)
function x = dirichlet_sample(alphaVec)
    y = gamrnd(alphaVec,1);
    x = y ./ sum(y);
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
    % plot_theta_distributions(theta_samples, outcomes, subjLabels);
    plot_hypothesis_testing(theta_samples, subjLabels);
    % plot_model_comparison(theta_samples);
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
            histogram(theta_samples{i}(:,k), ...
                      'Normalization','probability', ...
                      'FaceColor',cols(k,:), ...
					  'EdgeColor','none', ...
                      'FaceAlpha',0.5);
            xline(meanVals(i,k), ...
                  'Color',cols(k,:), ...
                  'LineWidth',3, ...
                  'HandleVisibility','off');
        end
        title(subjLabels{i}, 'FontSize',12);
        xlabel('Posterior Distribution','FontSize',12);
        ylabel('Density','FontSize',12);
        legend(outcomes{1:3}, 'FontSize',10,'Location','best');
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
                  'FaceColor',colDB, ...
				  'EdgeColor','none',...
                  'FaceAlpha',0.5);
        histogram(diffDC, ...
                  'Normalization','probability', ...
                  'FaceColor',colDC, ...
				  'EdgeColor','none',...
                  'FaceAlpha',0.5);

        xline(0,      'k--','LineWidth',1);
        xline(mean(diffDB), 'Color',colDB,'LineWidth',2,'HandleVisibility','off');
        xline(mean(diffDC), 'Color',colDC,'LineWidth',2,'HandleVisibility','off');

        title(sprintf('%s: P(D>B)=%.2f, P(D>C)=%.2f', ...
                      subjLabels{i}, pDB(i), pDC(i)), ...
              'FontSize',12);
        xlabel('\theta Difference','FontSize',12);
        ylabel('Density','FontSize',12);
        legend({'diffDB','diffDC'},'FontSize',10,'Location','best');
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
    elseif tC >= tB && tC >= tD && tB >= tD
        model = 3;
    else
        model = 4;
    end
end
