% Calculates basic measures of accuracy, RT and learning, LR - alpha, beta
% and WSLSbeta

% Alicia Rybicki - July 2022

clearvars; close all;clc


%% specify directories

origdir  = '/Users/rybickia-admin/Documents/Projects_analysis/SLT/Social_Learning_Course/SLS_analysis/';
datadir = ([origdir 'data/gorilla']);
%datadir = ([origdir 'data/yourdata']);
addpath(genpath(origdir));
cd(datadir)

% % Rename files 
% for f = 1:length(files)
%     oldName = files(f).name;
% 
%     [~, nameNoExt, ext] = fileparts(oldName);
% 
%     newSuffix = nameNoExt(end-1:end);
% 
%     newName = [newSuffix '_2' ext];
% 
%     movefile(fullfile(datadir, oldName), fullfile(datadir, newName));
% 
%     fprintf('Renamed: %s -> %s\n', oldName, newName);
% end


%% data - Get subject list from data folder
filePattern = fullfile(datadir, '*.csv');
files = dir(filePattern);

for f = 1:length(files)

    filename = [files(f).name];
    dash = strfind(filename, '_');
    loc = dash(1);
    names{f} = filename(1:loc-1);
    subjects = unique(names);
end

%% ======================== ANALYSIS LOOP ===========================

out = table();

% Initialize cell arrays -  parameter vectors for each subject
all_p_prc_p = {};  % Cell array to store est.p_prc.p for each subject
all_p_obs_p = {};  % Cell array to store est.p_obs.p for each subject
all_p_prc_p_sls = {};  % Cell array to store est_sls.p_prc.p for each subject  
all_p_obs_p_sls = {};  % Cell array to store est_sls.p_obs.p for each subject


for subj = 1:length(subjects)

    day = 2;
    cd(datadir)

    Filename = [char(subjects(:,subj)) '_' num2str(day) '.csv'];

    if exist(Filename,"file")
        d = readtable(Filename);
    else
    end

    fprintf(1, 'Now analysing %s\n', Filename);

    %% clean up input data file

    out.sID(subj,1) = d.ParticipantPrivateID(1);
    d.Properties.VariableNames{27} = 'randomiser';
    rand = d.randomiser(1);

    group = d.Spreadsheet(1);
    A = {'A'};
    B = {'B'};

    if contains(group,A)
        schedule = 1; % A
    else
        schedule = 2; % B
    end


    data = d(find(contains(d.ScreenName,'offer')),1:end); %#ok<FNDSB>

    data.block = double(strcmp(data.condition,'expert')); % 1 = expert, 0 = inexpert

    switch schedule
        case 1

            data.choices = double(strcmp(data.Response,'redBox.png')); % comment for pilot
            data.actual_corr_choice = double(strcmp(data.ANSWER,'redBox.png')); % comment for pilot
        case 2
            data.choices = double(strcmp(data.Response,'yellowBox.png')); % comment for pilot
            data.actual_corr_choice = double(strcmp(data.ANSWER,'yellowBox.png')); % comment for pilot
    end

    data.accuracy = data.Correct;
    data.correct = data.Correct;
    data.groupAccuracy = data.group_corr;
    data.conformity = 1-(abs(data.group_choice-data.choices));

    side_left = 0;
    for l = 1:length (data.choices)
        if (data.blue_pos(l) == 1 && data.choices(l) == 1) || (data.blue_pos(l) == 0 && data.choices(l) == 0)
            side_left = side_left + 1;
        else
        end
    end

    %% Volatility

    IndivProb  = data.prob_blue_corr;
    SocialProb  = data.prob_group_corr;

    volatile{1} = [ones(1,30) zeros(1,(60-30)) ones(1,(90-60)) zeros(1,(120-90))];  %Individual % 1 if volatile, 0 if stable
    volatile{2} = [zeros(1,30) ones(1,(60-30)) zeros(1,(90-60)) ones(1,(120-90))]; % Social

    stableIndiv = [31:60, 91:120];
    volIndiv = [1:30, 61:90];
    stableSocial = [1:30, 61:90];
    volSocial = [31:60, 91:120];

    switch schedule
        case 1
            expert = 1:60;
            inexpert = 61:120;

        case 2

            expert = 61:120;
            inexpert = 1:60;
    end


    %%
    cd ([origdir '/tapas'])  % volatile social

    inputs_reward = data.actual_corr_choice;
    response = data.choices;
    response_group2 = data.conformity;
    inputs_advice = data.group_choice;
    inputs_groupcorrectness = data.groupAccuracy;
    condition = data.block;

    est = tapas_fitModel_vol_social(response, [inputs_reward inputs_groupcorrectness ], volatile, inputs_advice, condition); %

    all_p_prc_p{subj} = est.p_prc.p;
    all_p_obs_p{subj} = est.p_obs.p;

    %% To simulate the data - uncomment the next BLOCK :

    sim = tapas_simModel(est.u, 'tapas_rw_social_reward_vol', [est.p_prc.p], volatile, inputs_advice, 'rw_softmax_constant_weight_social_reward', [est.p_obs.p]);

    % est.y  = responses;
    % est.u  = inputs (blue/group);
    % est.p_prc.p - estimates of perceptual paramters
    % est.p_obs.p - estimates of observation paramters
    % outputs: sim.y  = simulated data

   
    %%
    % separate for expert and inexpert, two LR

    volatile{1} = data.block'; % 1 if expert, 0 if inexpert
    volatile{2} = data.block';

    est_sls = tapas_fitModel_vol_social(response, [inputs_reward inputs_groupcorrectness ], volatile, inputs_advice, condition); %
    
    all_p_prc_p_sls{subj} = est_sls.p_prc.p;
    all_p_obs_p_sls{subj} = est_sls.p_obs.p;
    %% save output
    cd(datadir)
    out.subj(subj,1) = subj;


    out.schedule(subj,1) = schedule;
    out.acc(subj,1) = mean(data.accuracy);
    out.conf(subj,1) = mean(data.conformity);
    out.alien_acc(subj,1) = mean(data.groupAccuracy);
    out.side_check(subj,1) = sum(side_left)/120;
    out.col_check(subj,1) = sum(data.choices)/120;
    out.acc_expert(subj,1)  = mean(data.correct(data.block == 1));
    out.acc_inexpert(subj,1)  = mean(data.correct(data.block == 0));
    out.rt(subj,1)  = mean(data.ReactionTime);
    out.rt_expert(subj,1)  = mean(data.ReactionTime(data.block == 1));
    out.rt_inexpert(subj,1)  = mean(data.ReactionTime(data.block == 0));
    out.acc_vol(subj,1)  = mean(data.correct(volIndiv));
    out.acc_stable(subj,1) = mean(data.correct(stableIndiv));
    out.conf_vol(subj,1)  = mean(data.conformity(volSocial));
    out.conf_stable(subj,1)  = mean(data.conformity(stableSocial));
    out.conf_expert(subj,1)  = mean(data.conformity(data.block == 1));
    out.conf_inexpert(subj,1)  = mean(data.conformity(data.block == 0));

    % RW model - vol social

    out.indiv_LR_vol(subj,1) = est.p_prc.al_v_r;
    out.indiv_LR_stable(subj,1) = est.p_prc.al_s_r;
    out.soc_LR_vol(subj,1) = est.p_prc.al_v_a;
    out.soc_LR_stable(subj,1) = est.p_prc.al_s_a;

    out.v_0_r(subj,1) = est.p_prc.vr_0;
    out.v_0_a(subj,1) = est.p_prc.va_0;
    out.beta(subj,1) = est.p_obs.ze2;
    out.zeta(subj,1) = est.p_obs.ze1;

    BIC(1,subj,1) = est.optim.BIC;
    AIC(1,subj,1) = est.optim.AIC;
    Fvalues(1,subj,1) = est.optim.LME;

    % RW model - expert/inexpert social

    out.indiv_LR_expert(subj,1) = est_sls.p_prc.al_v_r;
    out.indiv_LR_inexpert(subj,1) = est_sls.p_prc.al_s_r;
    out.soc_LR_expert(subj,1) = est_sls.p_prc.al_v_a;
    out.soc_LR_inexpert(subj,1) = est_sls.p_prc.al_s_a;

    out.v_0_r_sls(subj,1) = est_sls.p_prc.vr_0;
    out.v_0_a_sls(subj,1) = est_sls.p_prc.va_0;
    out.beta_sls(subj,1) = est_sls.p_obs.ze2;
    out.zeta_sls(subj,1) = est_sls.p_obs.ze1;
    BIC(2,subj,1) = est_sls.optim.BIC;
    AIC(2,subj,1) = est_sls.optim.AIC;
    Fvalues(2,subj,1) = est_sls.optim.LME;

    all_choices(subj, :) = data.choices'; 

end

% Calculate average parameters across all subjects
avg_p_prc_p = mean(cell2mat(all_p_prc_p'), 1);
avg_p_obs_p = mean(cell2mat(all_p_obs_p'), 1);
avg_p_prc_p_sls = mean(cell2mat(all_p_prc_p_sls'), 1);
avg_p_obs_p_sls = mean(cell2mat(all_p_obs_p_sls'), 1);
real_choices = mean(all_choices, 1);
remove = [];
for w = 1:height(out)

    if out.sID(w) == 0
        remove = [remove, w];
    else
    end
end
out(remove,:) = [];

%sanity_figs
cd(origdir)
save('data_csc','out')
writetable(out, 'data_csc.txt')
save('model_comparisons','BIC', 'AIC','Fvalues')

est.p_prc.p = avg_p_prc_p;  % Average perceptual parameters
est.p_obs.p = avg_p_obs_p;  % Average observation parameters
save('averaged_parameters', 'est')
save('inputs_advice', 'inputs_advice');
save('real_choices', 'real_choices'); 
%% plots

%clearvars
close all

origdir  = '/Users/rybickia-admin/Documents/Projects_analysis/SLT/Social_Learning_Course/SLS_analysis/';
% workdir = ([origdir '/all_data_gorilla/data_2']);% Prolific data + RPS from V33 onwards
datadir = ([origdir 'gorilla']);
addpath(genpath(origdir));
cd(origdir)
load data_csc.mat

% accuracy
figdata{1} = [out.acc_expert;out.acc_inexpert]; % PLA
figcond{1} = [ones(1,height(out))';2*ones(1,height(out))']; %
%rt
figdata{2} = [out.conf_expert;out.conf_inexpert]; % PLA
figcond{2} = [ones(1,height(out))';2*ones(1,height(out))']; %

niceManyGroupsPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {'Accuracy', 'Conformity'}, 'condlab', {'expert', 'inexpert'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
%set(gca, 'YLim', [0 0.3]);
title('SLS')
ylabel({'Average'})



%% Learning rate alpha by volatility
figure
clear figdata figcond

figdata{1} = [out.indiv_LR_vol;out.soc_LR_vol]; % action
figcond{1} = [ones(1,height(out))';2*ones(1,height(out))']; %

figdata{2} = [out.indiv_LR_stable;out.soc_LR_stable]; % PLA]; % colour
figcond{2} = [ones(1,height(out))';2*ones(1,height(out))']; %

niceManyGroupsPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {'volatile', 'stable'}, 'condlab', {'Individual', 'Social'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
set(gca, 'YLim', [0 1]);
title('SLS')
ylabel({'Learning rate \alpha)'})


%% Learning rate alpha by expert status

figure
clear figdata figcond

figdata{1} = [out.indiv_LR_expert;out.soc_LR_expert]; % action
figcond{1} = [ones(1,height(out))';2*ones(1,height(out))']; %

figdata{2} = [out.indiv_LR_inexpert;out.soc_LR_inexpert]; % PLA]; % colour
figcond{2} = [ones(1,height(out))';2*ones(1,height(out))']; %

niceManyGroupsPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {'Expert', 'Inexpert'}, 'condlab', {'Individual', 'Social'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
set(gca, 'YLim', [0 1]);
title('SLS')
ylabel({'Learning rate \alpha)'})

%% 
load('model_comparisons.mat')

Fvalues = [Fvalues(1:2,:,:)];
[ep,out] = VBA_groupBMC(Fvalues);


AICgroup = -AIC(1:2,:,:);
[h, p] = VBA_groupBMC(AIC);
% 
BICgroup = -BIC(1:2,:,:);
[h, p] = VBA_groupBMC(BIC);


%%
% Figures  

% model comparison

% Posterior prob
figure
subplot(1,2,1)    
figs = out.Ef;

set(gca, 'YLim', [0 2]);

set(gca, 'Xdir', 'reverse');
b = bar([1 2 ],[figs(1); figs(2)]);

view([90 90]);
ylabel('p(y|m)'); 
xticklabels({'model 1','model 2'});

b.FaceColor = [.9 .6 .5];
set(gca,'FontSize',15)
yticks([0 0.25 0.5 1 1.25])

% Exceedance probability
subplot(1,2,2)
figs = out.ep';
set(gca, 'YLim', [0 1]);

set(gca, 'Xdir', 'reverse');
b = bar([1 2 ],[figs(1,1); figs(2,1)]);

view([90 90]);
ylabel('ϕ'); 
xticklabels({'model 1','model 2'});

b.FaceColor = [.8 .5 .4];
yticks([0 1])
set(gca,'FontSize',15)


%% Model Validation
cd ([origdir '/tapas'])  % volatile social
% parameter recovery
% correlation between actual and simulated response data


% Simulate response data using estimated model parameter values (tapas_simModel.m). 

% Number of simulations
n_sims = 500;

% Pre-allocate structure to store simulation results
sim_results = struct();


volatile{1} = [ones(1,30) zeros(1,(60-30)) ones(1,(90-60)) zeros(1,(120-90))];
volatile{2} = [zeros(1,30) ones(1,(60-30)) zeros(1,(90-60)) ones(1,(120-90))];

load('averaged_parameters.mat')
load('inputs_advice.mat')
% Initialize arrays to store simulated data
sim_results.y = cell(n_sims, 1);  % Simulated responses for each run
sim_results.parameters.p_prc = est.p_prc.p;  % original perceptual parameters
sim_results.parameters.p_obs = est.p_obs.p;  % original observation parameters
sim_results.inputs.u = est.u;  % Store original inputs
sim_results.inputs.volatile = volatile;  % Store volatility inputs
sim_results.inputs.advice = inputs_advice;  % Store advice inputs

% Run simulation loop
fprintf('Running %d simulations...\n', n_sims);
for i = 1:n_sims
    if mod(i, 20) == 0  % Progress indicator every 20 simulations
        fprintf('Completed %d/%d simulations\n', i, n_sims);
    end
    
    % Run single simulation
    sim = tapas_simModel(est.u, 'tapas_rw_social_reward_vol', ...
                        [est.p_prc.p], volatile, inputs_advice, ...
                        'rw_softmax_constant_weight_social_reward', ...
                        [est.p_obs.p]);

    % Store simulated responses
    sim_results.y{i} = sim.y;
    
    % store other simulation outputs
    % sim_results.other_field{i} = sim.other_field;
end

% Save results to .mat file
cd (origdir )  % volatile social

save('tapas_simulation_results.mat', 'sim_results');

fprintf('Simulation completed! Results saved to tapas_simulation_results.mat\n');

% summary statistics
fprintf('\nSummary:\n');
fprintf('Number of simulations: %d\n', n_sims);
fprintf('Number of trials per simulation: %d\n', length(sim_results.y{1}));

% Calculate mean response across all simulations
all_responses = cell2mat(sim_results.y);  % Convert cell array to matrix
mean_response_per_trial = mean(all_responses, 2);
fprintf('Mean response across all simulations: %.3f\n', mean(mean_response_per_trial));


%% 
% STATS ON SIMUALTED DATA - Learning rates 

%% Plot real vs simulated choice data over time
cd(origdir)

% Real data choices (averaged across subjects)
load('real_choices.mat'); real_choices = real_choices;

% Simulated data (averaged across n_sims)
% Convert simulations into matrix: [n_sims x n_trials]
all_responses = cellfun(@(x) x(:)', sim_results.y, 'UniformOutput', false);
all_responses = vertcat(all_responses{:});   % [n_sims x 120]

% average across simulations-  mean choice probability per trial
mean_sim_choices = mean(all_responses, 1);   

figure; hold on
plot(real_choices, 'b--', 'LineWidth', 2); % Real choices
plot(mean_sim_choices, 'k--', 'LineWidth', 2); % Simulated mean choices
xlabel('Trial'); ylabel('Choice (0/1)');
legend({'Real data','Simulated data'}, 'Location','best');
title('Real vs Simulated Choices over Time');


% fit model and (re)estimate parameters

%% Re-fit model to simulated data to estimate learning rates
n_sims = length(sim_results.y);
recovered_params = struct();
recovered_params.indiv_LR_vol = nan(n_sims,1);
recovered_params.indiv_LR_stable = nan(n_sims,1);
recovered_params.soc_LR_vol = nan(n_sims,1);
recovered_params.soc_LR_stable = nan(n_sims,1);
recovered_params.beta = nan(n_sims,1);
recovered_params.zeta = nan(n_sims,1);

for i = 1:n_sims
    fprintf('Fitting model to simulation %d/%d...\n', i, n_sims);

    % Simulated responses
    sim_response = sim_results.y{i};

    % Fit learning model as above.. 
    est_sim = tapas_fitModel_vol_social(sim_response, ...
        [inputs_reward inputs_groupcorrectness], ...
        sim_results.inputs.volatile, ...
        sim_results.inputs.advice, ...
        condition);

    % recovered learning rates
    recovered_params.indiv_LR_vol(i)    = est_sim.p_prc.al_v_r;
    recovered_params.indiv_LR_stable(i) = est_sim.p_prc.al_s_r;
    recovered_params.soc_LR_vol(i)      = est_sim.p_prc.al_v_a;
    recovered_params.soc_LR_stable(i)   = est_sim.p_prc.al_s_a;
    recovered_params.beta(i) = est_sim.p_obs.ze2;
    recovered_params.zeta(i) = est_sim.p_obs.ze1;
end

save('recovered_params.mat', 'recovered_params');

%% Plot recovered learning rates

figure
clear figdata figcond

figdata{1} = [recovered_params.indiv_LR_vol;recovered_params.soc_LR_vol]; % action
figcond{1} = [ones(1,500)';2*ones(1,500)']; %

figdata{2} = [recovered_params.indiv_LR_stable;recovered_params.soc_LR_stable]; % PLA]; % colour
figcond{2} = [ones(1,500)';2*ones(1,500)']; %

niceManyGroupsPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {'volatile', 'stable'}, 'condlab', {'Individual', 'Social'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
set(gca, 'YLim', [0 1]);
title('SLS')
ylabel({'Recovered learning rate \alpha)'})

% Plot recovered beta and zeta 
load('data_csc.mat')
figure
clear figdata figcond

figdata{1} = [recovered_params.beta;out.beta]; % action
n = height(out)
figcond{1} = [ones(1,500)';2*ones(1,n)']; %

niceSoloPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {''},'condlab', {'Recovered', 'Actual'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
%set(gca, 'YLim', [0 1]);
title('Recovered vs actual \beta')
ylabel({'\beta'})
legend off

figure
clear figdata figcond

figdata{1} = [recovered_params.zeta;out.zeta]; % action
n = height(out)
figcond{1} = [ones(1,500)';2*ones(1,n)']; %

niceSoloPlot(figdata,figcond, 'dotsize',5,'col1', 9, 'col2', 7, ...
    'gplab', {''},'condlab', {'Recovered', 'Actual'}, ...
    'transp', 0.3,'whatplot', 3);
set(gca, 'FontSize', 26);
set(gca,'FontName', 'Arial')
%set(gca, 'YLim', [0 1]);
title('Recovered vs actual \zeta')
ylabel({'\zeta'})
legend off

