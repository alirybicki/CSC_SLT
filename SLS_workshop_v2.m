% Rename gorilla files

% Alicia Rybicki - August 2025

clearvars; close all;clc


origdir  = '/Users/rybickia-admin/Documents/Projects_analysis/SLT/Social_Learning_Course/SLS_analysis/';
datadir = ([origdir 'data/yourdata']);

cd(datadir)

%%
filePattern = fullfile(datadir, '*.csv');
files = dir(filePattern);

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

% Calculates basic measures of accuracy, RT and learning, LR - alpha, beta
% and WSLSbeta


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

for i = 1:length(subjects)
    subj = subjects{i};  % subj is now a string like '68'
    day = 2;
    cd(datadir)

    Filename = [subj '_' num2str(day) '.csv'];  % <- use subj directly

    if exist(Filename, "file")
        d = readtable(Filename);
       
    else
        fprintf(1, 'File %s not found\n', Filename);
    end

    fprintf(1, 'Now analysing %s\n', Filename);

    %% clean up input data file

    out.sID(subj,1) = d.ParticipantPrivateID(1);
    d.Properties.VariableNames{27} = 'randomiser';
    rand = d.randomiser(1);

    group = d.SpreadsheetName{1};
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
    data.Correct = str2double(string(data.Correct)); % converts 'NaN' -> NaN
    data.Correct(isnan(data.Correct)) = 0;


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


    %% To simulate the data - uncomment the next line :

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





    %% save output
    cd(datadir)
    
    out.subj{i,1} = subj;  


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






end


%%

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
save('gorilla_model_comparisons','BIC', 'AIC','Fvalues')

%% plots

clearvars
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
load('gorilla_model_comparisons.mat')

Fvalues = [Fvalues(1:2,:,:)];
[ep,out] = VBA_groupBMC(Fvalues);


AICgroup = -AIC(1:2,:,:);
[h, p] = VBA_groupBMC(AIC);
% 
BICgroup = -BIC(1:2,:,:);
[h, p] = VBA_groupBMC(BIC);


%%
% Suppl. Figures  

% model comparison

% Posterior prob
figure
subplot(1,2,1)    
figs = out.Ef;

set(gca, 'YLim', [0 2]);

% flip plot so 'x' axis is vertical and 'y' is horizontal

set(gca, 'Xdir', 'reverse');
b = bar([1 2 ],[figs(1); figs(2)]);

view([90 90]);
ylabel('p(y|m)'); 
xticklabels({'model 1','model 2'});

b.FaceColor = [.9 .6 .5];
set(gca,'FontSize',15)
yticks([0 0.25 0.5 1 1.25])

%%Exceedance probability
subplot(1,2,2)
figs = out.ep';
set(gca, 'YLim', [0 1]);

% flip plot so 'x' axis is vertical and 'y' is horizontal

set(gca, 'Xdir', 'reverse');
b = bar([1 2 ],[figs(1,1); figs(2,1)]);

view([90 90]);
ylabel('EP'); 
xticklabels({'model 1','model 2'});

b.FaceColor = [.8 .5 .4];
yticks([0 1])
set(gca,'FontSize',15)
