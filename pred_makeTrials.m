function [s] = pred_makeTrials(p, s, devNum, data)


s.trials_headers = {'precueValidity','tOrientation','tPhase', 'tContrast','tStatus','precue','responseKey','response','accuracy','rt'};
% make sure column indices match trials headers
s.precueValidityIdx = strcmp(trials_headers,'precueValidity');
s.testOrientationIdx = strcmp(trials_headers,'tOrientation'); % 0 (low) or 1 (high)
s.testPhaseIdx = strcmp(trials_headers,'tPhase'); % test phase value
s.testContrastIdx = strcmp(trials_headers,'tContrast'); % test contrast
s.testStatusIdx = strcmp(trials_headers,'tStatus'); % test contrast
s.precueIdx=strcmp(trials_headers,'precue'); % precue
s.responseKeyIdx = strcmp(trials_headers,'responseKey');
s.responseIdx = strcmp(trials_headers,'response');
s.accuracyIdx = strcmp(trials_headers,'accuracy');
s.rtIdx = strcmp(trials_headers,'rt');


s.trials = fullfact([numel(p.precueValidities),...
    numel(p.gratingOrientations),...
    numel(p.testPhases),...
    numel(p.gratingContrasts)...
    numel(p.testStatus)]);

s.trials = repmat(s.trials, p.reps, 1); 


