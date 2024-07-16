% Calculates basic measures of accuracy, RT and learning, LR - alpha, beta
% and WSLSbeta

% Alicia Rybicki - July 2022

clearvars; close all;clc


%% specify directories

origdir  = '/Users/rybickia-admin/Documents/Projects_analysis/SLT/Social_Learning_Course/SLS_analysis/';
datadir = ([origdir 'data/sls']);
yourdata = ([origdir 'data/gorilla']);
addpath(genpath(origdir));

cd(datadir)
cd(yourdata)
%% ______  Rename files if using data directly from gorilla 
% 
% filePattern = fullfile(yourdata, '*.csv');
% subjects = dir(filePattern);
% for subj = 1:length(subjects) %loop through data folder
%     Filename = [subjects(subj).name];
%     if exist(Filename,"file")
%         d = readtable(Filename);
%     else
%     end
%     if height(d) > 0
%         sID = d.ParticipantPrivateID(1);
%         d.Properties.VariableNames{27} = 'randomiser';
%         rand = d.randomiser(1);
% 
%         group = d.Spreadsheet(1);
%         A = {'A'};
%         B = {'B'};
%         if contains(group,A)
%             schedule = 1; % A
%         else
%             schedule = 2; % B
%         end
%         if contains(rand,A)
%             if schedule == 1
%                 day = 1
%             else
%                 day = 2
%             end
%         else
%             if schedule == 2
%                 day = 1
%             else
%                 day = 2
% 
%             end
%         end
%         movefile(Filename,[num2str(sID) '_' num2str(day) '.csv']);
% 
% 
% 
%     end
% end
% 



%% Or just use sample data - Get subject list from data folder
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

sess = 1;

out = table();

for subj = 1:length(subjects)

    for j = 1
        day = sess(j);
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
                data.choices = double(strcmp(data.Response,'blueBox.png')); % comment for pilot
                data.actual_corr_choice = double(strcmp(data.ANSWER,'blueBox.png')); % comment for pilot
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
        % Initilaise variables for simulations

        % Initialise some variables
        model_acc = [];
        est_acc_vol =[];
        est_acc_stable =[];

        responses_sched1 = [];
        ModPrediction_sched1 = [];

        model_acc_vol = [];
        model_acc_stable = [];

        vDirect = [];
        vIndirect = [];
        peDirect = [];
        peIndirect = [];
        count = 1;



        %% RESCORLA-WAGNER ANALYSIS
        %
        cd ([origdir '/tapas'])  %Alicia

        inputs_reward = data.actual_corr_choice;
        response = data.choices;
        response_group2 = data.conformity;
        inputs_advice = data.group_choice;
        inputs_groupcorrectness = data.groupAccuracy;
        condition = data.block;

        est = tapas_fitModel_vol_social(response, [inputs_reward inputs_groupcorrectness ], volatile, inputs_advice, condition); %
   
        %% To simulate the data - uncomment the next line :

         sim = tapas_simModel(est.u, 'tapas_rw_social_reward_vol', [est.p_prc.p], volatile, inputs_advice, 'rw_softmax_constant_weight_social_reward', [est.p_obs.p]);
         
         % % est.y  = responses;
         % % est.u  = inputs (blue/group);
         % % est.p_prc.p - estimates of perceptual paramters
         % % est.p_obs.p - estimates of observation paramters

         % % outputs: sim.y  = simulated data 


         %% Parameter recovery  - uncomment the next line :
         
         % % estimate paramters from simulated data (sim.y) 
         est_recovered = tapas_fitModel_vol_social([sim.y], [inputs_groupcorrectness inputs_reward ], volatile, inputs_advice, group); %


        %%
        % separate for expert and inexpert
        volatile{1} = data.block'; % 1 if expert, 0 if inexpert
        volatile{2} = data.block';

        inputs_reward = data.actual_corr_choice;
        response = data.choices;
        response_group2 = data.conformity;
        inputs_advice = data.group_choice;
        inputs_groupcorrectness = data.groupAccuracy;
        condition = data.block;

        est_sls = tapas_fitModel_vol_social(response, [inputs_reward inputs_groupcorrectness ], volatile, inputs_advice, condition); %


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

        % RW Learning rate overall
        out.indiv_LR_vol(subj,1) = est.p_prc.al_v_r;
        out.indiv_LR_stable(subj,1) = est.p_prc.al_s_r;
        out.soc_LR_vol(subj,1) = est.p_prc.al_v_a;
        out.soc_LR_stable(subj,1) = est.p_prc.al_s_a;

        out.v_0_r(subj,1) = est.p_prc.vr_0;
        out.v_0_a(subj,1) = est.p_prc.va_0;
        out.beta(subj,1) = est.p_obs.ze2;
        out.zeta(subj,1) = est.p_obs.ze1;
        BIC(subj,1) = est.optim.BIC;
        AIC(subj,1) = est.optim.AIC;
        LME(subj,1) = est.optim.LME;

        % expert/inexpert instrad of vol/stable

        out.indiv_LR_expert(subj,1) = est_sls.p_prc.al_v_r;
        out.indiv_LR_inexpert(subj,1) = est_sls.p_prc.al_s_r;
        out.soc_LR_expert(subj,1) = est_sls.p_prc.al_v_a;
        out.soc_LR_inexpert(subj,1) = est_sls.p_prc.al_s_a;

        out.v_0_r_sls(subj,1) = est_sls.p_prc.vr_0;
        out.v_0_a_sls(subj,1) = est_sls.p_prc.va_0;
        out.beta_sls(subj,1) = est_sls.p_obs.ze2;
        out.zeta_sls(subj,1) = est_sls.p_obs.ze1;
        %AIC_sls(subj,1) = est_sls.optim.AIC;
        %LME_sls(subj,1) = est_sls.optim.LME;


        out.sim_LR_vol_indiv(subj,1) = est_recovered.p_prc.al_v_r;
        out.sim_LR_stable_indiv(subj,1) = est_recovered.p_prc.al_s_r;
        out.sim_LR_vol_soc(subj,1) = est_recovered.p_prc.al_v_a;
        out.sim_LR_stable_soc(subj,1) = est_recovered.p_prc.al_s_a;
        out.sim_LR_vr_0(subj,1) = est_recovered.p_prc.vr_0;
        out.sim_LR_va_0(subj,1) = est_recovered.p_prc.va_0;
        out.sim_ze1(subj,1) = est_recovered.p_obs.ze1;
        out.sim_beta(subj,1) = est_recovered.p_obs.ze2;
                


    end
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

%% plots

clearvars
close all

origdir  = '/Users/rybickia-admin/Documents/Projects_analysis/SLT/Social_Learning_Course/SLS_analysis/';
% workdir = ([origdir '/all_data_gorilla/data_2']);% Prolific data + RPS from V33 onwards
datadir = ([origdir 'sls']);
addpath(genpath(origdir));

load data_sls_adults_wide.mat

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
