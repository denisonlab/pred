function [s] = pred_trials(p, s, devNum, data)
% Makes trial structure for hot adapt 

%% Generate trials
% Construct trials matrix
s.trials_headers = {'sOrientation','sPhase','sContrast','tOrientation','tPhase, tContrast'}; 

% Make sure column indices match trials headers
s.standOrientation = strcmp(s.trials_headers,'sOrientation'); % 1 (-45) ot 2 (45)
s.standPhase = strcmp(s.trials_headers,'sPhase'); % standard phase
s.standContrast = strcmp(s.trials_headers,'sContrast'); % standard contrast
s.testOrientation = strcmp(s.trials_headers,'tOrientation'); % 0 (low) or 1 (high)
s.testPhase = strcmp(s.trials_headers,'testPhase'); % test phase value
s.testContrast = strcmp(s.trials_headers,'tContrast'); % test contrast

%target orientations same in HOT as in prediction: -45, 45,
% shu wade3 contrast threshold multiplier

s.trials = fullfact([numel(p.gratingOrientations),... 
    numel(p.targetState),...
    numel(p.gratingContrasts) ...
    numel(p.testPhase)]); 
s.trials = repmat(s.trials, p.nRepsPerBlock, 1); 

