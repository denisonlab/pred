function predPlaid
% ___________________________________________________________________
% pred.m runs Waffle experiment
% Melissa Allouche
% January 2024
%
% p = parameters
% s = stimuli, rects, trials
% window = window
% data = data
%
% Code adapted from Angus Chapman triad_attn, Karen Tian hot-adapt,
% Juneau Wang ta-auditory, and Rachel Denison ta-demo, ta-template
% __________________________________________________________________
clear; close all;
PsychPortAudio('Close');
rng("shuffle");
% if ~exist(pahandle,'var')
%     PsychPortAudio('Stop', pahandle);
%     PsychPortAudio('Close', pahandle);
%     Screen('CloseAll');
% end
%% Input
% Subject and session info
p.subjectID = input('Enter subject ID:  ','s');
p.sessionNum = input('Enter session number (1,2...) ');
p.counter=input('Choose 1 or 2 for participant \n'); %counterbalance 
p.eyeTracking=input('Eyetracking (0/1)? ');
p.debug=input('Debug (0/1)? ');

%% Setup
% Add paths
directory = pwd; % get project directory path, set to prediction folder parent level
addpath(genpath(directory))
% Get parameters
p = predParamsPlaid(p);

if p.counter==2
    p.toneFreqs=flip(p.toneFreqs); 
end

data.dataDir = sprintf('%s/data',pwd);
if ~exist(data.dataDir, 'dir')
    mkdir(data.dataDir)
end

data.dataDir_subjectID = sprintf('%s/data/%s',pwd,p.subjectID);
if ~exist(data.dataDir_subjectID , 'dir')
    mkdir(data.dataDir_subjectID )
end


data.dataDir_sub = sprintf('%s/data/%s/session_%d',pwd,p.subjectID,p.sessionNum);
if ~exist(data.dataDir_sub, 'dir')
   mkdir(data.dataDir_sub)
end

% data.subDir = sprintf('%s/%s',data.dataDir,p.subjectID);
% if ~exist(data.subDir, 'dir')
%     mkdir(data.subDir)
% end
% are we overwriting data
date = datestr(now, 'yymmdd_HHMM'); 
if exist(sprintf('%s/%s_s%d_predplaid_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'file')
    error('This subject/session already has saved data');
end

% Running on PTB-3? Abort otherwise
AssertOpenGL;

%% Setup Eyetracking
if p.eyeTracking==1
    eyeDataDir = 'eyedata';
    eyeFile = sprintf('%s%s', p.subjectID(1:2), datestr(now, 'mmdd'));
    data.eyeDataDir = sprintf('%s/eyedata',pwd);
    if ~exist(data.eyeDataDir, 'dir')
        mkdir(data.eyeDataDir)
    end
    % Check to see if this eye file already exists
    if exist(sprintf('%s/%s_s%d_eye.edf',data.eyeDataDir,p.subjectID,p.sessionNum),'file')
        error('This subject/session already has saved eye data');
    end
end
%% Display key settings to the experimenter
fprintf('\nExperiment settings:\n')
fprintf('subject = %s\n', p.subjectID)
fprintf('aperture = %s\n', p.aperture)

%% Screen

% Skip screen tests - ONLY for demo, not for real
Screen('Preference', 'SkipSyncTests', 1);
screenNumber = max(Screen('Screens')); %screen to display on

% Get window size
%[window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);

[window, rect] = Screen('OpenWindow', screenNumber);

[screenWidthPx, screenHeight] = Screen('WindowSize', window);
flipInterval = Screen('GetFlipInterval', window); % % Get refresh rate, frame duration (s)

% We will request the screen to flip half a refresh (the "slack" time) before we
% actually want the screen to change. This helps to avoid late screen
% flips. So let's define this "slack" variable for convenience.
slack = flipInterval/2;
p.slack=slack;

% Get x and y coordinates for the center of the window
[cx, cy] = RectCenter(rect);
white = WhiteIndex(window);
black = BlackIndex(window);

%% Calculate stimulus dimensions (px) and position
pixelsPerDegree = ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);
fixSize = p.fixSize*pixelsPerDegree;
gaborSD = round(p.gaborSD * p.ppd); 

%% Making sounds ...
% 10^0.5 for every 10dB
%% Make a pure tone for each cue frequency
cueTones = [];
for iF = 1:numel(p.toneFreqs)
    tone0 = MakeBeep(p.toneFreqs(iF), p.toneDur, p.Fs);
    % Apply an envelope so the sound doesn't click at the beginning and end
    tone = applyEnvelope(tone0, p.Fs);
    cueTones(iF,:) = tone;
end
%cueTones(iF+1,:) = mean(cueTones,1); % neutral precue, both tones together, DO WE NEED THIS?
p.cueTones=cueTones;

%% Keyboard
% Check all "devices" (keyboards, mice) for response input
devNum = -1;

KbName('UnifyKeyNames');
validKeys = KbName({'9(','0)'});



%% MAKE PLAID REFERENCE STIMULI FOR VERSION 2 HERE
for iC = 1:numel(p.plaidContrasts1) 
    contrast1 = p.plaidContrasts1(iC);
    for iP=1:numel(p.testPhases)
        phase = p.testPhases(iP);
        for iS=1:numel(p.plaidSPF)
            SPF=p.plaidSPF(iS);
            grating1 = rd_grating(pixelsPerDegree, p.imSize(1), ...
                SPF, p.plaidAxes(1), phase, contrast1); % 0 to 1
            grating2 = rd_grating(pixelsPerDegree, p.imSize(1), ...
                 p.plaidSPFconstant, p.plaidAxes(2), phase, p.plaidContrasts2); % 0 to 1
            plaid = (grating1 + grating2) - 0.5; 
            
            [plaid, aps] = rd_aperture(plaid, 'gaussian', gaborSD);
            tex_plaid{iC,iP} = Screen('MakeTexture', window, plaid*white);
        end
    end
end

imRectL = CenterRectOnPoint([0 0 imSize], cx+imPos(1)-2*p.plaidEcc*p.ppd, cy+imPos(2));
imRectR = CenterRectOnPoint([0 0 imSize], cx+imPos(1)+2*p.plaidEcc*p.ppd, cy+imPos(2));

%% Define fixation box for eyetracker
fixBoxWidth=5*pixelsPerDegree; %width of fixation box in deg
fixRect = [cx-.5*fixBoxWidth, cy-.5*fixBoxWidth, cx+.5*fixBoxWidth, cy+.5*fixBoxWidth];
%rad = 70; 
rad=round(ang2pix(p.eyerad,p.screenWidthCm, screenWidthPx, p.viewDistCm,'central')); % radius of allowable eye movement in pixels
%rad=50;
%% RUN expt
%HideCursor(window,-1);

%% Plaid (waffles) trials
trials_headers2 = {'precueValidity','plaidOrientation','plaidPhase' 'plaidContrast'};
% make sure column indices match trials headers
precueValidityIdx2 = strcmp(trials_headers2,'precueValidity'); %precue validity for plaid
plaidOrientationIdx = strcmp(trials_headers2,'plaidOrientation'); % plaid orientation
plaidPhaseIdx = strcmp(trials_headers2,'plaidPhase'); % plaid phase value
plaidContrastIdx = strcmp(trials_headers2,'plaidContrast'); % plaid contrast

trials2 = fullfact([ numel(p.precueValiditiesWaffle)... % 1 precue validity
     numel(p.plaidOrientations),... % 2 plaid ori
     numel(p.testPhases),... % 3 plaid phase 
     numel(p.plaidContrasts1)]); % 4 plaid contrast
%% Merge trial count
nTrials = size(trials2,1);  % total trials = number of grating trials + number of waffle trials
nBlocks=nTrials/p.BlockTrials; 

trialOrder = randperm(nTrials); %randomize trial order

%% Eyetracker
if p.eyeTracking==1
    % Initialize eye tracker
    [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
    if exitFlag
        return
    end

    % Write subject ID into the edf file
    Eyelink('message', 'BEGIN DESCRIPTIONS');
    Eyelink('message', 'Subject code: %s', p.subjectID);
    Eyelink('message', 'END DESCRIPTIONS');

    el.drift_correction_target_beep = [0 0 0];
    el.drift_correction_failed_beep = [0 0 0];
    el.drift_correction_success_beep = [0 0 0];

    % Accept input from all keyboards
    el.devicenumber = -1; %see KbCheck for details of this value

    % Update with custom settings
    EyelinkUpdateDefaults(el);

    % Calibrate eye tracker
    [cal, exitFlag] = rd_eyeLink('calibrate', window, el);
    if exitFlag
        return
    end
    fixations=[];
end

%% Sound
% Initialize the sound driver
InitializePsychSound(1); % 1 for precise timing
PsychPortAudio('Close');

% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel

if p.eyeTracking
    rd_eyeLink('startrecording',window,{el, fixRect});
end

%% Show instruction screen and wait for a button press
instructions = 'This is the plaid version of the experiment\n\n';
Screen('FillRect', window, white*p.backgroundColor);
instructions1 = sprintf('%s\n\nThere will be a tone that signals the start of the trial followed by a plaid. \n\n Press to continue!', instructions);
DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

KbWait(devNum);

instructions2 ='\n\n You must report which of the two orientations in the plaid is stronger in contrast. \n\n For -45 degrees (press 9).  \n\n For +45 degrees (press 0)! \n\n  Press to continue!';
DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white); DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack);
WaitSecs(1);
KbWait(devNum);

timeStart = GetSecs;
correct = [];
testStronger=[];
block=1;
eyeSkip = zeros(size(trials1,1)+size(trials2,1),1); % trials skipped due to an eye movement, same size as trials matrix
skippedTrials = [];
iTrial=1;
completedTrials=0;

disp('reached')
stairIdxExp=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
lastFewAccExp=[]; % no previous grating trials have occured so this is empty
stairIdxUn=length(p.stairs);
lastFewAccUnexp=[];  % no previous grating trials have occured so this is empty

while iTrial<=nTrials
    trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
    stopThisTrial=0;
    %% Initialize for eye tracking trial breaks
     if iTrial>1
         eyeSkip(iTrial-1) = stopThisTrial; % this is for the previous trial
     end
     stopThisTrial = 0;

    %% Get condition information for this trial

    plaidStatus=2; % if the trial id is a value > the number of grating trials, this trial will be a waffle trial
    precueValidity = p.precueValiditiesWaffle(trials2(trialIdx-size(trials1,1), precueValidityIdx2));
    plaidOrientation = trials2(trialIdx-size(trials1,1), plaidOrientationIdx);
    plaidPhase = trials2(trialIdx-size(trials1,1), plaidPhaseIdx);
    plaidContrast = trials2(trialIdx-size(trials1,1), plaidContrastIdx);
    pOrientation=p.plaidOrientations(plaidOrientation); %get the orientation value (+/-45) using the plaid orientation index

    %% tone related to overall orientation
    toneName = p.precueNames{precueValidity};
    switch toneName
        case 'valid'
            tone=cueTones(plaidOrientation,:);
            toneVersion=plaidOrientation;
        case 'invalid'
            if  plaidOrientation==1
                tone = cueTones(2,:);
                toneVersion=2;
            elseif  plaidOrientation==2
                tone = cueTones(1,:);
                toneVersion=1;
            end
        otherwise
            error('precueName not recognized')
    end
    %% Store stimulus information in trials matrix
    %% %%%% Play the trial %%%%
    
    d.plaidOrientation(iTrial) = plaidOrientation; %store orientation of plaid stimuli
    d.plaidPhase(iTrial) = plaidPhase; %store phase of plaid stimuli
    d.plaidContrast(iTrial) = plaidContrast; %store contrast of plaid stimuli
    d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
    d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
    d.plaidContrastValue(iTrial)=p.plaidContrasts1(plaidContrast);

    d.block(iTrial)=block;
    d.timeEnd(iTrial) = NaN;
    d.timeElapsed(iTrial) = NaN;
    d.timeTargetResponse(iTrial) = NaN;
    d.targetRT(iTrial) = NaN;
    d.targetResponseKey(iTrial) = NaN;
    d.timeSpentSaving(iTrial) = NaN;

   
    
    % Present fixation rest (grey)
    drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
    timePreStart = Screen('Flip', window);
    % Present fixation active (white)
    drawFixation(window, cx, cy, fixSize, p.fixColor*white);
    timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
    if p.eyeTracking
        Eyelink('Message', 'FixOn')
    end
    %% Check fixation hold
    if p.eyeTracking
        driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad,fixRect});
        if driftCorrected
            % restart trial
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
        end
    end
    %% Present predictive tone
    PsychPortAudio('FillBuffer', pahandle, tone);
    timeTone = PsychPortAudio('Start', pahandle, [], timeFix + p.signalStart, 1); % waitForStart = 1 in order to return a timestamp of playback
    
    %% EYE TRACKING

    if p.eyeTracking
        Eyelink('Message', 'EVENT_CUE');
    end

    if p.eyeTracking
        while GetSecs < timeTone + p.toneSOA - p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            fixations = [fixations fixation];

            if fixation==0
                stopThisTrial = 1;
                WaitSecs(1);
            
                if iTrial==1
                    firstNonWaffle=0;
                end
                % redo this trial at the end of the experiment
                % this can be easily done by appending the trial number to the end of
                % trialOrder
                trialOrder(end+1) = trialOrder(iTrial);
                d.stopThisTrial(iTrial) = stopThisTrial;
                skippedTrials(end+1) = trialOrder(iTrial);
                d.timeTargetResponse(iTrial) = NaN;
                d.correct(iTrial) = NaN;
                d.testStronger(iTrial)=NaN;

                d.block(iTrial)=block;
                d.timeStart(iTrial)=timeStart;
                d.timeTone(iTrial)=timeTone;
                d.timePreStart(iTrial)=timePreStart;
                d.timeFix(iTrial)=timeFix;
                d.timeS(iTrial) = NaN;
                d.timeBlank1(iTrial) = NaN;

                nTrials = nTrials + 1;
                iTrial = iTrial + 1;

                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize,p.fixColor*p.dimFactor*white);
                instructionsSkip = 'Fixation lost!';
                DrawFormattedText(window, instructionsSkip, 'center', cy-imRect(1)/2, [1 1 1]*white);
                Screen('Flip', window,p.demoInstructDur-slack);
                WaitSecs(1);
            else 
                stopThisTrial = 0;
            end
        end

        if stopThisTrial
            continue
        end
    end
    %% Present stimuli
    Screen('DrawTexture', window, tex_plaid{plaidContrast,plaidPhase}, [], imRect, pOrientation);
    drawFixation(window, cx, cy, fixSize, p.fixColor*white)
    
    timeS = Screen('Flip', window, timeTone + p.toneSOA - slack); % timeFix+ how much i want to wait from white (active) to standard
        
    if p.eyeTracking
        Eyelink('Message', 'EVENT_CUE');
    end

    if p.eyeTracking
        while GetSecs < timeS + p.imDur- p.eyeSlack && ~stopThisTrial
            WaitSecs(.01);
            fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
            fixations = [fixations fixation];

            if fixation==0
                stopThisTrial = 1;
                WaitSecs(1);

                if iTrial==1
                    firstNonWaffle=0;
                end
            
                % redo this trial at the end of the experiment
                % this can be easily done by appending the trial number to the end of
                % trialOrder
                trialOrder(end+1) = trialOrder(iTrial);
                d.stopThisTrial(iTrial) = stopThisTrial;
                skippedTrials(end+1) = trialOrder(iTrial);
                d.timeTargetResponse(iTrial) = NaN;
                d.correct(iTrial) = NaN;
                d.testStronger(iTrial)=NaN;

                d.timeStart(iTrial)=timeStart;
                d.timeTone(iTrial)=timeTone;
                d.timePreStart(iTrial)=timePreStart;
                d.timeFix(iTrial)=timeFix;
                d.timeS(iTrial)=timeS;
               
                d.timeT(iTrial) = NaN;
                d.timeBlank1(iTrial) = NaN;

                nTrials = nTrials + 1;
                iTrial = iTrial + 1;

                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize,p.fixColor*p.dimFactor*white);
                instructionsSkip = 'Fixation lost!';
                DrawFormattedText(window, instructionsSkip, 'center', cy-imRect(1)/2, [1 1 1]*white);
                Screen('Flip', window,p.demoInstructDur-slack);
                WaitSecs(1);
            else 
                stopThisTrial = 0;
            end
        end

        if stopThisTrial
            continue
        end
    end

    % blank
    if p.debug
        debugCue = sprintf('Precue= %d !',precueValidity);
        DrawFormattedText(window, debugCue, 'center', cy-imRect(1)/2, [1 1 1]*white);

    end
    drawFixation(window, cx, cy, fixSize, p.fixColor*white);
    timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
 
    %% Eyetracking 
    if p.eyeTracking
        Eyelink('Message', 'EVENT_RESPCUE');
    end
     
    %% Wait for response
    % check only valid response keys
    targetResponseKey = []; targetResponseKeyName = []; 
    while isempty(targetResponseKey) 
        [timeTargetResponse, keyCode] = KbWait(devNum);
        targetRT = timeTargetResponse - timeS;
        targetResponseKey = find(ismember(validKeys,find(keyCode)));
        targetResponseKeyName = KbName(validKeys(targetResponseKey));
    
        if iscell(targetResponseKeyName)
            targetResponseKey = []; 
            targetResponseKeyName = []; 
        end
    end
    if isempty(targetResponseKey)
        targetRT = NaN; % timeout
        targetResponseKey = NaN;
        targetResponseKeyName = NaN;
    end
    if p.eyeTracking
        Eyelink('Message', 'TRIAL_END');
    end

   %% Percent correct for orientation & test stronger
   if plaidStatus==2 && strcmp('0)',targetResponseKeyName) % if they respond with +45
       if plaidOrientation==2 % and orientation is +45
           testStronger=1; % they reported test as stronger
       elseif plaidOrientation==1 % and orientation is -45
           testStronger=0; % they did not report test as stronger
       end
   elseif plaidStatus==2 && strcmp('9(',targetResponseKeyName) % if they respond with -45
       if plaidOrientation==2 % and orientation is +45
           testStronger=0; % they did not report test as stronger
       elseif plaidOrientation==1 % and orientation is -45
           testStronger=1;  % they reported test as stronger
       end
   end
    %% Response collected (blue fixation)
    if targetResponseKey
        drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
        Screen('Flip',window);
    
    end
    
    timeEnd=GetSecs(); 

    d.plaidOrientation(iTrial) = plaidOrientation; %store orientation of plaid stimuli
    d.plaidPhase(iTrial) = plaidPhase; %store phase of plaid stimuli
    d.plaidContrast(iTrial) = plaidContrast; %store contrast of plaid stimuli
    d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
    d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
   
    %% Store trial info
    p.trials2=trials2;
    tic
    % save data
    d.block(iTrial) = block; %block trial is in
    d.timeStart(iTrial)=timeStart; % time at which trial starts
    d.timeEnd(iTrial)=timeEnd; % time at which trial ends
    d.timeElapsed(iTrial)=timeEnd-timeStart; % trial duration
    d.timeTargetResponse(iTrial)=timeTargetResponse; % time at which response was made
    d.targetRT(iTrial) = targetRT; %response time duration
    d.targetResponseKey(iTrial) = targetResponseKey; % response key for trial
    d.timePreStart(iTrial)=timePreStart; % time before standard stimuli presentation
    d.timeFix(iTrial)=timeFix; % time for fixation
    d.timeS(iTrial)=timeS; % time for standard
    d.timeBlank1(iTrial)=timeBlank1; % time after standard
    d.timeTone(iTrial)=timeTone; % time at tone
    d.stopThisTrial(iTrial) = stopThisTrial;
    d.testStronger(iTrial)=testStronger;
    d.trialOrder(iTrial)=trialOrder;
    save(sprintf('%s/%s_s%d_predplaid_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

    %time between trials after response
    d.timeSpentSaving(iTrial)=toc;
    WaitSecs(p.ITI-d.timeSpentSaving(iTrial));
    
    completedTrials = completedTrials + 1; 
    %% SANTIY CHECK
    fprintf('Trial %d/%d in block %d, trial %d of %d total \n', completedTrials, p.BlockTrials, block,iTrial,nTrials);
    %% Blocking
    if (mod(completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block)  
        % Calculate block accuracy
        blockStartTrial = ((iTrial-completedTrials)/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
        if blockStartTrial < 0 % we are doing less than one block
            blockStartTrial = 1;
        end
        fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
       

     
        blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
        if iTrial==nTrials
            keyMessage = 'All done! Thank you for participating!';
        else
            keyMessage = 'Press any key to go on.';
        end

        breakMessage = sprintf('%s\n%s', blockMessage, keyMessage);
        DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window);
        WaitSecs(1);
        if iTrial < nTrials
            KbWait(-1);
        end

        block = block+1; % keep track of block for block message only
        completedTrials=0;
    end
    iTrial = iTrial + 1;
end

if p.eyeTracking
    p.eyefixation = fixations;
    p.eyefixT1 = timeFix;
    d.eyeSkip=eyeSkip;
end

%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});

    % rename eye file
    eyeFileFull = sprintf(sprintf('%s/%s_s%d_plaid_eye.edf',data.eyeDataDir,p.subjectID,p.sessionNum),'file');
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
end

%% Clean Up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')
ShowCursor;
end
