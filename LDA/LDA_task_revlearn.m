%% LDA Script
% Converts filtered EEG data into LDA format (X) and runs LDA
clear all
close all

%% Set-Up

% set filepath from main project folder (e.g. 'DATA/EEG/')
filepath = '/DATA/derivatives/';
analyses_path='/ANALYSIS/';

subFind = dir(fullfile(filepath,'sub-*')); %find subjects to loop
sub = {subFind.name};

% get_X
offset = [-100:10:850]; % sliding window of analysis in samples (1000Hz)
tbase = 100; % number of samples subtracted in baseline correction during epoching

% single_trial_analysis 
dur = 60; % size of window in samples (centred around offset)
skipLOO = 0; % whether to skip the leave-one-out validation (1 = skip)
perm = [0 500 0.01]; % permutation test for Az sig: [flag(1|0), nperms, sig_thresh]
eigratio = 0; % input for logistpca(), minimum ratio of eigval of each component to 
              % largest eigval. Components with ratio < eigratio are excluded.
gamma = 0; % gamma = 1: discard covariance values, adding some degree of regularisation: 
           % gamma=0.01/0.05 (near 0)
method = 0; % 0: logistic regression, 1: Regularised Fisher's 

%% Run Single Trial Analysis

%define contrasts
cope={'rewPosVSrewNeg' 'punNegVSpunPos'};

% Loop through subjects
for sidx = 1:length(sub)
    
    % set current subject
    sj = sub{sidx};
    
    fprintf('Processing subject: %s\n',sj);
    
    %load events data
    filename = [sj, '/eeg/', sj, '_task-revlearn'];
    load([filepath, filename, '_events']);
    
    % load epoched EEG data
    load([filepath, filename, '_epoched']);

    % set parameters for get_X
    % cond_id sets CONDITION A = 1 and CONDITION B = 2, where A and B are
    % being discriminated (e.g. negative vs positive feedback, as below)
    for copeID=1:length(cope)
    
        fprintf('Running contrast: %s\n',cope{copeID}); 

         if copeID==1
             cond_id =((2*events.fdb_win.*events.valence)+(events.fdb_neutral.*events.valence)); 
         elseif copeID==2               
             cond_id = ((events.fdb_loss.*~events.valence)+(2*events.fdb_neutral.*~events.valence));     
         end

           cond_discr = [1,2]; 
    
    % Get X and run single_trial_analysis on cond_discr
    fprintf('Offset: ')
        for i = 1:length(offset);

            fprintf('%d, ', offset(i))

            % get X 
            [X,truth] = get_X(allData, cond_id, cond_discr, offset(i)+tbase, dur);

            % Single Trial Analysis (Logistic)
            [Azloo(i),~,Y(i),a(i,:),v(i,:),D(i)] = single_trial_analysis(...
                X, truth, dur, skipLOO, perm, eigratio, gamma, method);   
        end

    file_name=sprintf('%s_LDA.mat',cope{copeID});
       
    save([analyses_path, sj, '/eeg/', file_name],...
             'Azloo', 'Y', 'a', 'v', 'D', 'tbase', 'offset', 'method');

    clear X truth cond_id cond_discr Azloo Y a v D 
    
    end
end

fprintf('Done!')

%% to plot - optional
hold on
plot(offset,Azloo,'k');
xlabel('Time (ms)')
ylabel ('Az')
title('Stim1 vs. Stim2', 'FontSize',12);
legend boxoff
