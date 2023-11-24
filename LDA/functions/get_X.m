function [X, truth] = get_X(eegdata, cond_id, cond_discr, offset, dur)

% eegdata:      epoched data in 3D format - CHANNELS x TIME x TRIALS
% cond_id:      vector of condition ID's for each trial [1 X TRIALS]
% cond_discr:   two element vector with condition IDs for discrimination - e.g. [1 4]
% offset:       onset time of analysis window in samples (relative to trial onset)
% duration:     duration of analysis window in samples

% get trial indices for condition A
idx_A = find(cond_id==cond_discr(1));
% get trial indices for condition B
idx_B = find(cond_id==cond_discr(2));

% build a logical (binary) truth labels vector (0: condition A, 1: condition B)
truth = [zeros(1,length(idx_A)) ones(1,length(idx_B))];

% extract all desired data and build an analysis matrix with dimensions
% Channels x [Analysis Duration x Trials]
X = [reshape(eegdata(:,offset+1:offset+dur,idx_A),size(eegdata,1),dur*length(idx_A)) ...
     reshape(eegdata(:,offset+1:offset+dur,idx_B),size(eegdata,1),dur*length(idx_B))];

