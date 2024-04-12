function pred
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
p.debug = input('Debug? (Y/N) ','s');
p.task = input(['Prediction task run:\n' ...
    '1 - Task original\n'...
    '2 - Demo\n'...
    '3 - Task version 2 43008 trials\n'...
    '4 - Task version 2 800 trials\n'...
    '5 - Task version Kok + Waffles\n' ...
    '6 - Demo for Kok + Waffles\n']);
p.reps = input ('How many reps(1/2)? ');
p.fullScreen = input ('Full screen(0/1)? ');
p.eyeTracking=input('Eyetracking (0/1)? ');
p.maskType=input('Masktype (none/color)? ','s');

%% Setup
% Add paths
directory = pwd; % get project directory path, set to prediction folder parent level
addpath(genpath(directory))
% Get parameters
p = predParams(p);

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
if exist(sprintf('%s/%s_s%d_pred_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'file')
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
    if exist(sprintf('%s/%s.edf',eyeDataDir, eyeFile),'file')
        error('This subject/session already has saved eye data');
    end
end
%% Display key settings to the experimenter
fprintf('\nExperiment settings:\n')
fprintf('subject = %s\n', p.subjectID)
fprintf('aperture = %s\n', p.aperture)
fprintf('expt stage = %d\n', p.task)

%% Screen

% Skip screen tests - ONLY for demo, not for real
Screen('Preference', 'SkipSyncTests', 1);
screenNumber = max(Screen('Screens')); %screen to display on


% Get window size
%[window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);

if p.fullScreen==1
    [window, rect] = Screen('OpenWindow', screenNumber);
elseif p.fullScreen==0
    [window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);
end

% Debug settings
if p.debug=="Y" || p.debug=="y"
    nTrials=p.debugTrials;
    nBlocks=nTrials/p.BlockTrialsDebug;
    % elseif p.debug=="N" || p.debug=="n"
    %     [window, rect] = Screen('OpenWindow', screenNumber);
end

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
if p.task==1 || p.task==2
    validKeys = KbName({'1!','2@','9(','0)'});
elseif p.task==3 || p.task==4 ||p.task==5
    validKeys = KbName({'9(','0)'});
end
%% Make TEST stimuli

imPos = round(p.imPos*pixelsPerDegree); % image position (pixels)
% gratingRadius = round(p.gratingDiameter/2*pixelsPerDegree); % radius (pixels)

tex = cell(numel(p.gratingPhases), numel(p.gratingSPF));
tex_stand= cell(numel(p.testPhases), numel(p.gratingSPF));

% for iP=1:numel(p.testPhases)
%         phase = p.testPhases(iP);
%         for iS=1:numel(p.gratingSPF)
%             SPF = p.gratingSPF(iS);
%             grating = rd_grating(pixelsPerDegree, p.imSize(1), ...
%                 SPF, 0, phase, p.standardContrast2); % 0 to 1
%             [gabor, aps] = rd_aperture(grating, 'gaussian', gaborSD);
%             tex{iP,iS} = Screen('MakeTexture', window, gabor*white);
%         end
%  end

if p.task ~=5
    for iC = 1:numel(p.gratingContrasts)  %HOW TO LOOP THROUGH PHASE OR VARY ACROSS PHASES
        contrast = p.gratingContrasts(iC);
        for iP=1:numel(p.testPhases)
            phase = p.testPhases(iP);
            grating = rd_grating(pixelsPerDegree, p.imSize(1), ...
                p.gratingSF, 0, phase, contrast); % 0 to 1
            [gabor, aps] = rd_aperture(grating, 'gaussian', gaborSD);
            tex{iC,iP} = Screen('MakeTexture', window, gabor*white);
        end
    end
elseif p.task==5
 
    for iP=1:numel(p.gratingPhases) %TEST
        phase = p.testPhases(iP);
        for iS=1:numel(p.gratingSPF)
            SPF = p.gratingSPF(iS);
            grating = rd_grating(pixelsPerDegree, p.imSize(1), ...
               SPF, 0, phase, p.gratingContrast2 ); % 0 to 1
             [gabor, aps] = rd_aperture(grating, 'gaussian', gaborSD);
            tex{iP,iS} = Screen('MakeTexture', window, gabor*white);
        end
    end
    

    for iP=1:numel(p.testPhases) %STANDARD
        phase = p.testPhases(iP);
        for iS=1:numel(p.gratingSPF)
            SPF = p.gratingSPF(iS);
            grating = rd_grating(pixelsPerDegree, p.imSize(1), ...
                SPF, 0, phase, p.standardContrast2); % 0 to 1
            [gabor, aps] = rd_aperture(grating, 'gaussian', gaborSD);
            tex_stand{iP,iS} = Screen('MakeTexture', window, gabor*white);
        end
    end
end
 
%Make rects for placing image
imSize = size(grating);
imRect = CenterRectOnPoint([0 0 imSize], cx+imPos(1), cy+imPos(2));


%% Make STANDARD stimuli

gratingS = rd_grating(pixelsPerDegree, p.gratingSize, p.gratingSF, 0, 0, p.standardContrast); % First make a grating image
imPosS = round(p.imPos*pixelsPerDegree); % standard image position (pixels)
% gratingRadiusS = round(p.gratingDiameter/2*pixelsPerDegree); % standard radius (pixels)

% Place an Gaussian aperture on the image to turn it into a Gabor
gaborS = rd_aperture(gratingS, 'gaussian', gaborSD);


% Make an image "texture"
% Choose which image you want to make into a texture
imS = gaborS;

% Make the texture
imTexStandard = Screen('MakeTexture', window, imS*white); % multiply by "white" to scale from 0-255
imSizeS = size(imS); % Make the rects for placing the images in the window
imRectS = CenterRectOnPoint([0 0 imSizeS(1) imSizeS(2)], cx+imPosS(1), cy+imPosS(2));

%% MAKE PLAID REFERENCE STIMULI FOR VERSION 2 HERE
for iC = 1:numel(p.plaidContrasts1) 
    contrast1 = p.plaidContrasts1(iC);
    for iP=1:numel(p.testPhases)
        phase = p.testPhases(iP);
        grating1 = rd_grating(pixelsPerDegree, p.imSize(1), ...
            p.plaidSPF, p.plaidAxes(1), phase, contrast1); % 0 to 1
        grating2 = rd_grating(pixelsPerDegree, p.imSize(1), ...
             p.plaidSPF, p.plaidAxes(2), phase, p.plaidContrasts2); % 0 to 1
        plaid = (grating1 + grating2) - 0.5; 
        
        [plaid, aps] = rd_aperture(plaid, 'gaussian', gaborSD);
        tex_plaid{iC,iP} = Screen('MakeTexture', window, plaid*white);
        
    end
end

imRectL = CenterRectOnPoint([0 0 imSize], cx+imPos(1)-2*p.plaidEcc*p.ppd, cy+imPos(2));
imRectR = CenterRectOnPoint([0 0 imSize], cx+imPos(1)+2*p.plaidEcc*p.ppd, cy+imPos(2));

%% Define fixation box for eyetracker
fixBoxWidth=5*pixelsPerDegree; %width of fixation box in deg
fixRect = [cx-.5*fixBoxWidth, cy-.5*fixBoxWidth, cx+.5*fixBoxWidth, cy+.5*fixBoxWidth];
rad = 70; % radius of allowable eye movement in pixels

%% RUN expt
%HideCursor(window,-1);


switch p.task % task and demo
    case 1
        %% 1 TASK
        % Make trials structure
        %% %%%% Generate trials in different conditions %%%%

        trials_headers = {'precueValidity','tOrientation','tPhase', 'tContrast','tStatus','precue','responseKey','response','accuracy','rt'};
        % make sure column indices match trials headers
        precueValidityIdx = strcmp(trials_headers,'precueValidity');
        testOrientationIdx = strcmp(trials_headers,'tOrientation'); % 0 (low) or 1 (high)
        testPhaseIdx = strcmp(trials_headers,'tPhase'); % test phase value
        testContrastIdx = strcmp(trials_headers,'tContrast'); % test contrast
        testStatusIdx = strcmp(trials_headers,'tStatus'); % test status absent or present
        precueIdx=strcmp(trials_headers,'precue'); % precue
        responseKeyIdx = strcmp(trials_headers,'responseKey');
        responseIdx = strcmp(trials_headers,'response');
        accuracyIdx = strcmp(trials_headers,'accuracy');
        rtIdx = strcmp(trials_headers,'rt');


        trials = fullfact([numel(p.precueValidities),...
            numel(p.gratingOrientations),...
            numel(p.testPhases),...
            numel(p.gratingContrasts)...
            numel(p.testStatus)]);
        %if p.debug=="Y" || p.debug=="y"
        %nTrials=p.debugTrials;
        %nBlocks=nTrials/p.BlockTrialsDebug; %60 trials/2 blocks
        if (p.debug=="N" || p.debug=="n") && p.reps==1
            trials=repmat(trials,p.reps*p.repScale1,1);
            nTrials = size(trials,1); %448*2= 896, each trial twice
            nBlocks=nTrials/p.BlockTrials; %((448)*2)/64 = 14 blocks
        elseif (p.debug=="N" || p.debug=="n") && p.reps==2
            trials=repmat(trials,p.reps*p.repScale2,1); %((448)*3)=
            nTrials = size(trials,1);
            nTrialsNoSkip=nTrials;
            nBlocks=p.trialsNeeded/nTrials; %21 blocks 3 mins each
        end

        trialOrder = randperm(nTrials);

        %% Eyetracker
        if p.eyeTracking
            % Initialize eye tracker
            [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
            if exitFlag
                return
            end

            % Write subject ID into the edf file
            Eyelink('message', 'BEGIN DESCRIPTIONS');
            Eyelink('message', 'Subject code: %s', p.subjectID);
            Eyelink('message', 'END DESCRIPTIONS');

            % No sounds indicating success of calibration
            %     el.targetbeep = false;
            %     el.calibration_failed_beep = [0 0 0];
            %     el.calibration_success_beep = [0 0 0];
            el.drift_correction_target_beep = [0 0 0];
            el.drift_correction_failed_beep = [0 0 0];
            el.drift_correction_success_beep = [0 0 0];

            % Accept input from all keyboards
            el.devicenumber = -1; %see KbCheck for details of this value

            % Update with custom settings
            EyelinkUpdateDefaults(el);

            % Calibrate eye tracker
            [~, exitFlag] = rd_eyeLink('calibrate', window, el);
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

        %% generate mask for every trial and generate standard mask for every trial too, use rd_aperture
        tic
        maskTexs=[]; 
        masks = []; maskTexs = []; 
        for i= 1:nTrials
            for iJ = 1:2 % Standard then test
                [mask,~,~]=kt_makeFilteredNoise(p.gratingSize, p.maskContrast,0, 0, p.gratingSF/2, p.gratingSF*2,pixelsPerDegree, 0,'allOrientations');
                masks1=mask;
                % gratingRadiusM = round(p.gratingDiameter/2*pixelsPerDegree); 
                mask=rd_aperture(mask, 'gaussian', gaborSD);
                maskTex = Screen('MakeTexture', window, mask*white);

                masks(i,iJ,:,:) = mask; 
                maskTexs(i,iJ)=maskTex;
            end
        end
        toc

        %% Show instruction screen and wait for a button press
        instructions = 'This is the main experiment\n\n';
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a reference patch followed by a predictive tone and a second patch.\n\nA low tone predicts a CCW orientation.\n\nA high tone predicts a CW orientation. \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 ='Your goal is to determine whether the second patch is: \n\n higher contrast and counterclockwise (press 1), \n\nlower contrast and counterclockwise (press 2), \n\nlower contrast and clockwise (press 9), \n\nor higher contrast and clockwise (press 0) with reference to the first patch!\n\nPress to start!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);

        %pred_instructions(window, p, devNum);
        KbWait(devNum);
        timeStart = GetSecs;
        correct = [];
        block=1;
        %eyeSkip = zeros(size(trials,1),1); % trials skipped due to an eye movement, same size as trials matrix
        % Start recording
        %rd_eyeLink('startrecording',window,{el,fixRect});
        skippedTrials = [];
        iTrial=1;
        completedTrials=0;

        disp('reached')

        %for iTrial = 1:nTrials % the iteration in the trial loop
        while iTrial<=nTrials
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
            stopThisTrial = 0;
            % %% Initialize for eye tracking trial breaks
            % if iTrial>1
            %     eyeSkip(trialIdx) = stopThisTrial; % this is for the previous trial
            % end

            %% %%%% Present one trial %%%
            %% Get condition information for this trial
            precueValidity = p.precueValidities(trials(trialIdx, precueValidityIdx));
            testOrientation = trials(trialIdx, testOrientationIdx);
            testPhase = trials(trialIdx, testPhaseIdx);
            testContrast = trials(trialIdx, testContrastIdx);
            testStatus=p.testStatus(trials(trialIdx, testStatusIdx));

            % orientation
            orientation=p.gratingOrientations(testOrientation);
            % tone
            toneName = p.precueNames{precueValidity};
            switch toneName
                case 'valid'
                    tone=cueTones(testOrientation,:);
                    toneVersion=testOrientation;
                case 'invalid'
                    if testOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif testOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
                    end
                otherwise
                    error('precueName not recognized')
            end

            %% Set up test mask for this trial

            %% Store trial information
            d.testStatus(iTrial) = testStatus; %absent or present
            d.testOrientation(iTrial) = testOrientation; %orientation of test stimuli
            d.testPhase(iTrial) = testPhase; %phase of test stimuli
            d.testContrast(iTrial) = testContrast; %contrast of test stimuli
            d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)

            %% Store stimulus information in trials matrix

            trials(trialIdx, precueIdx) = toneVersion;
            %% %%%% Play the trial %%%%
            %% Present fixation rest (grey)
            drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
            timePreStart = Screen('Flip', window);
            %% Present fixation active (white)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
            
            if p.eyeTracking
                Eyelink('Message', 'FixOn')
            end
            %% SHOULD I SWITCH THESEEEEE SO THAT IT TURNS WHITE WHEN YOU ARE FIXATING !!!
            %% Check fixation hold
            if p.eyeTracking
                driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});

                if driftCorrected
                    % restart trial
                    drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                    timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
                end
            end

            %% Present STANDARD image
            Screen('DrawTexture', window, imTexStandard, [], imRectS, 0);
            drawFixation(window, cx, cy, fixSize, p.fixColor*white)
            timeS = Screen('Flip', window, timeFix+p.signalStart - slack); %timeFix+ how much i want to wait from white(active) to standard
            %sound(click, p.Fs)  %play click
            
            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
            %% PRESENT STANDARD MASK
            if ~strcmp(p.maskType,'none')
                Screen('DrawTexture', window, maskTexs(trialIdx,1), [], imRect);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeMask1 = Screen('Flip', window, timeS + p.maskSOA - slack);

                % blank
                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)                
                timeMaskBlank1 = Screen('Flip', window, timeMask1 + p.maskDur - slack);
            else
                timeMask1 = NaN;
                timeMaskBlank1 = NaN;
            end
            
            %% Present predictive tone
            PsychPortAudio('FillBuffer', pahandle, tone);
            timeTone = PsychPortAudio('Start', pahandle, [], timeS + p.standSOA, 1); % waitForStart = 1 in order to return a timestamp of playback
            
            timeTone = timeS + p.standSOA; 

            if p.eyeTracking
                Eyelink('Message', 'EVENT_CUE');
            end
            
            if p.eyeTracking
                while GetSecs < timeTone + p.standSOA - p.eyeSlack && ~stopThisTrial
                    WaitSecs(.01);
                    %             fixation = mod(iTrial,10); %%% for testing
                    fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                   
                    fixations = [fixations fixation];

                    if fixation==0
                        stopThisTrial = 1;
                        WaitSecs(1);
                        
                        trialOrder(end+1) = trialOrder(iTrial);
                        d.stopThisTrial(iTrial) = stopThisTrial;
                        skippedTrials(end+1) = trialOrder(iTrial);
                        d.timeTargetResponse(iTrial) = NaN;
                        d.correct(iTrial) = NaN;
                        nTrials = nTrials + 1;
                        iTrial = iTrial + 1;

                        % fixations = [fixations fixation];
                        %DrawFormattedText(window, 'Fixation lost. Please press space when ready to fixate.', 'center', 'center', [1 1 1]*white);
                        Screen('FillRect', window, white*p.backgroundColor);
                        Screen('Flip', window);
                        WaitSecs(1);
                        %Screen('Flip', window);
                        %KbWait(devNum);

                    else 
                        stopThisTrial = 0;
                    end
                end


                if stopThisTrial
                    continue
                end
            end

            %% Present TEST image
            %drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            if testStatus==1
                Screen('DrawTexture', window, tex{testContrast, testPhase}, [], imRect, orientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            elseif testStatus==0
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            end
            timeT = Screen('Flip', window, timeTone + p.toneSOA - slack); % is it p.toneSOA, DIDN'T GET LOGIC
           
            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            Screen('Flip', window, timeT + p.imDur - slack);

            %% PRESENT TEST MASK
            if ~strcmp(p.maskType,'none')
                Screen('DrawTexture', window, maskTexs(trialIdx,2), [], imRect);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeMask2 = Screen('Flip', window, timeT + p.maskSOA - slack);

                % blank
                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)                
                timeMaskBlank2 = Screen('Flip', window, timeMask2 + p.maskDur - slack);
            else
                timeMask2 = NaN;
                timeMaskBlank2 = NaN;
            end
            
            %% Eyetracking 
            if p.eyeTracking
                Eyelink('Message', 'EVENT_RESPCUE');
            end
      
            %% Wait for response
            % check only valid response keys
            targetResponseKey = [];
            while isempty(targetResponseKey)
                [timeTargetResponse, keyCode] = KbWait(devNum);
                targetRT = timeTargetResponse - timeMask2;
                targetResponseKey = find(ismember(validKeys,find(keyCode)));
                targetResponseKeyName = KbName(validKeys(targetResponseKey));
                correct = NaN;
            end
            if isempty(targetResponseKey)
                targetRT = NaN; % timeout
                targetResponseKey = NaN;
                targetResponseKeyName = NaN;
                correct = NaN;
            end
            if p.eyeTracking
                Eyelink('Message', 'TRIAL_END');
            end

            %% Response collected (blue fixation)
            if targetResponseKey
                drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
                Screen('Flip',window);
               
            end
            
            %% Percent correct for orientation

            if testStatus~=0 && (strcmp('9(',targetResponseKeyName) || strcmp('0)',targetResponseKeyName)) % response CCW
                if testOrientation==2 % stimuli CW
                    correct = 1;  % stimuli = response
                elseif testOrientation==1 % stimuli CCW
                    correct = 0; % stimuli != response
                end
            elseif testStatus~=0 && (strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName)) % response CCW
                if testOrientation==1 %stimuli CCW
                    correct = 1; % stimuli = response
                elseif testOrientation==2 %stimuli CW
                    correct = 0; % stimuli != response
                end
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            timeEnd=GetSecs();

            %% Store trial info
            trials(trialIdx, rtIdx) = targetRT;
            trials(trialIdx, responseKeyIdx) = targetResponseKey;
            p.trials=trials;

            tic
            % save data
            d.block(iTrial) = block; %block trial is in
            d.timeStart(iTrial)=timeStart; % time at which trial starts
            d.timeEnd(iTrial)=timeEnd; % time at which trial ends
            d.timeElapsed(iTrial)=timeEnd-timeStart; % trial duration
            
            d.correct(iTrial) = correct; % response correctness (0/1)
            d.timeTargetResponse(iTrial)=timeTargetResponse; % time at which response was made
            d.targetRT(iTrial) = targetRT; %response time duration
            d.targetResponseKey(iTrial) = targetResponseKey; % response key for trial
            d.timePreStart(iTrial)=timePreStart; % time before standard stimuli presentation
            d.timeFix(iTrial)=timeFix; % time for fixation
            d.timeS(iTrial)=timeS; % time for standard
            d.timeBlank1(iTrial)=timeBlank1; % time after standard
            d.timeTone(iTrial)=timeTone; % time at tone
            d.timeT(iTrial)=timeT; % time at test
            d.trialOrder=trialOrder;
            d.testMask(iTrial)=masks(trialIdx,2);
            d.standardMask(iTrial)=masks(trialIdx,1);
            %d.trials(iTrial)=trials(trialIdx, : );
            save(sprintf('%s/%s_s%d_pred_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

            %time between trials after response
            d.timeSpentSaving(iTrial)=toc;
            WaitSecs(p.ITI-d.timeSpentSaving(iTrial));
            
            completedTrials = completedTrials + 1; 
            %% SANTIY CHECK
            fprintf('Trial %d/%d in block %d, trial %d of %d total \n', completedTrials, p.BlockTrials, block,iTrial,nTrials);
            %% Blocking
            %if completedTrials==p.BlockTrials && (mod(iTrial-completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block) 
            if (mod(completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block)  
                % Calculate block accuracy
                %blockStartTrial = (iTrial/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                blockStartTrial = ((iTrial-completedTrials)/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                if blockStartTrial < 0 % we are doing less than one block
                    blockStartTrial = 1;
                end

                current=d.correct(blockStartTrial:iTrial);
                nan_correctnessIdx= isnan(current);
                %non_nan_correctness=d.correct(non_nan_correctnessIdx);
                non_nan_correctness=current(~nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);

                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
                DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
                %         DrawFormattedText(window, cueText, 'center', cy-3*textOffset, white);
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

        %% CASE 2
    case 2 % DEMO

        %% Sound
        % Initialize the sound driver
        InitializePsychSound(1); % 1 for precise timing
        % Open audio device for low-latency output
        reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
        pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel
        
        %% Generate trials in different conditions %%%%

        trials_headers = {'precueValidity','tOrientation','tPhase', 'tContrast','tStatus','precue','responseKey','response','accuracy','rt'};
        % make sure column indices match trials headers
        precueValidityIdx = strcmp(trials_headers,'precueValidity');
        testOrientationIdx = strcmp(trials_headers,'tOrientation'); % 0 (low) or 1 (high)
        testPhaseIdx = strcmp(trials_headers,'tPhase'); % test phase value
        testContrastIdx = strcmp(trials_headers,'tContrast'); % test contrast
        testStatusIdx = strcmp(trials_headers,'tStatus'); % test status absent or present
        precueIdx=strcmp(trials_headers,'precue'); % precue
        responseKeyIdx = strcmp(trials_headers,'responseKey');
        responseIdx = strcmp(trials_headers,'response');
        accuracyIdx = strcmp(trials_headers,'accuracy');
        rtIdx = strcmp(trials_headers,'rt');


        trials = fullfact([numel(p.precueValidities),...
            numel(p.gratingOrientations),...
            numel(p.testPhases),...
            numel(p.gratingContrasts)...
            numel(p.testStatus)]);

        repmat(trials,p.reps*p.repScaleDemo,1); % demo: 448 trials
        nTrials = size(trials,1); %448, each trial once
        nBlocks=nTrials/p.BlockTrials; %((448)/64 = 7 blocks


        trialOrder = randperm(nTrials);

          %% generate mask for every trial and generate standard mask for every trial 
        tic
        maskTexs=[]; 
        masks = []; maskTexs = []; 
        for i= 1:nTrials
            for iJ = 1:2 % Standard then test
                [mask,~,~]=kt_makeFilteredNoise(p.gratingSize, p.maskContrast,0, 0, p.gratingSF/2, p.gratingSF*2,pixelsPerDegree, 0,'allOrientations');
                masks1=mask;

                mask=rd_aperture(mask, 'gaussian', gaborSD);
                maskTex = Screen('MakeTexture', window, mask*white);

                masks(i,iJ,:,:) = mask; 
                maskTexs(i,iJ)=maskTex;
            end
        end
        toc
        instructions = 'This is a demo of the main experiment\n\n';
        %% Show instruction screen and wait for a button press
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a reference patch followed by a predictive tone and a second patch.\n\nA low tone predicts a CCW orientation.\n\nA high tone predicts a CW orientation. \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 ='Your goal is to determine whether the second patch is: \n\n higher contrast and counterclockwise (press 1), \n\nlower contrast and counterclockwise (press 2), \n\nlower contrast and clockwise (press 9), \n\nor higher contrast and clockwise (press 0) with reference to the first patch!\n\nPress to continue!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        timeInstruct2=Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);
        KbWait(devNum);

        instructionsToneLow ='This is what the low tone will sound like.\n\nPress to continue!';
        DrawFormattedText(window, instructionsToneLow, 'center', 'center', [1 1 1]*white);
        timeToneLowText=Screen('Flip', window,timeInstruct2-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'

        toneLow=cueTones(1,:);
        PsychPortAudio('FillBuffer', pahandle, toneLow);
        timeToneLow = PsychPortAudio('Start', pahandle, [], timeToneLowText+p.standSOA, 1); % waitForStart = 1 in order to return a timestamp of playback
        WaitSecs(1);

        KbWait(devNum);

        instructionsToneHigh ='This is what the high tone will sound like. Press to continue!';
        DrawFormattedText(window, instructionsToneHigh, 'center', 'center', [1 1 1]*white);
        timeToneHighText=Screen('Flip', window,timeToneLow+p.toneSOA-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        
        WaitSecs(1);
        toneHigh=cueTones(2,:);
        PsychPortAudio('FillBuffer', pahandle, toneHigh);
        timeToneHigh = PsychPortAudio('Start', pahandle, [], timeToneHighText+p.standSOA-slack, 1); % waitForStart = 1 in order to return a timestamp of playback
        
        toc
        WaitSecs(1);
        KbWait(devNum);

        instructionsStartDemo ='Press to start demo trials!';
        DrawFormattedText(window, instructionsStartDemo, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeToneHigh+p.toneSOA-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);
        KbWait(devNum);

        timeStart = GetSecs;
        correct = [];
        block=1;

        for iTrial = 1:nTrials % the iteration in the trial loop
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix

            %% Santiy Check
            fprintf('Trial %d/%d in block %d, trial %d of %d total \n\n', ceil(iTrial/p.BlockTrials), p.BlockTrials,iTrial,nTrials);

            %% %%%% Present one trial %%%
            %% Get condition information for this trial
            precueValidity = p.precueValidities(trials(trialIdx, precueValidityIdx));
            testOrientation = trials(trialIdx, testOrientationIdx);
            testPhase = trials(trialIdx, testPhaseIdx);
            testContrast = trials(trialIdx, testContrastIdx);
            testStatus=p.testStatus(trials(trialIdx, testStatusIdx));

            % orientation
            orientation=p.gratingOrientations(testOrientation);
            % tone
            toneName = p.precueNames{precueValidity};
            switch toneName
                case 'valid'
                    tone=cueTones(testOrientation,:);
                    toneVersion=testOrientation;
                case 'invalid'
                    if testOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif testOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
                    end
                otherwise
                    error('precueName not recognized')
            end

            %% Store stimulus information in trials matrix

            trials(trialIdx, precueIdx) = toneVersion;

            %% PLAY TRIAL
            %% Present fixation rest (grey)
            drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
            timePreStart = Screen('Flip', window);
            %% Present fixation active (white)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
            
            %% Present STANDARD image
            Screen('DrawTexture', window, imTexStandard, [], imRectS, 0);
            drawFixation(window, cx, cy, fixSize, p.fixColor*white)
            timeS = Screen('Flip', window, timeFix+p.signalStart - slack); %timeFix+ how much i want to wait from white(active) to standard
            %sound(click, p.Fs)  %play click
            
            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
            %% PRESENT STANDARD MASK
            if ~strcmp(p.maskType,'none')
                Screen('DrawTexture', window, maskTexs(trialIdx,1), [], imRect);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeMask1 = Screen('Flip', window, timeS + p.maskSOA - slack);

                % blank
                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)                
                timeMaskBlank1 = Screen('Flip', window, timeMask1 + p.maskDur - slack);
            else
                timeMask1 = NaN;
                timeMaskBlank1 = NaN;
            end
            
            %% Present predictive tone
            PsychPortAudio('FillBuffer', pahandle, tone);
            timeTone = PsychPortAudio('Start', pahandle, [], timeS + p.standSOA, 1); % waitForStart = 1 in order to return a timestamp of playback

            %% Present TEST image
            if testStatus==1
                Screen('DrawTexture', window, tex{testContrast, testPhase}, [], imRect, orientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);

            elseif testStatus==0
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            end
            timeT = Screen('Flip', window, timeTone + p.toneSOA - slack); % is it p.toneSOA, DIDN'T GET LOGIC

            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            Screen('Flip', window, timeT + p.imDur - slack);

            %% PRESENT TEST MASK
            if ~strcmp(p.maskType,'none')
                Screen('DrawTexture', window, maskTexs(trialIdx,2), [], imRect);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeMask2 = Screen('Flip', window, timeT + p.maskSOA - slack);

                % blank
                Screen('FillRect', window, white*p.backgroundColor);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)                
                timeMaskBlank2 = Screen('Flip', window, timeMask2 + p.maskDur - slack);
            else
                timeMask2 = NaN;
                timeMaskBlank2 = NaN;
            end
            
            %% Wait for response
            % check only valid response keys
            targetResponseKey = [];
            while isempty(targetResponseKey)
                [timeTargetResponse, keyCode] = KbWait(devNum);
                targetRT = timeTargetResponse - timeT;
                targetResponseKey = find(ismember(validKeys,find(keyCode)));
                targetResponseKeyName = KbName(validKeys(targetResponseKey));
                correct = NaN;
            end
            if isempty(targetResponseKey)
                targetRT = NaN; % timeout
                targetResponseKey = NaN;
                targetResponseKeyName = NaN;
                correct = NaN;
            end

            if targetResponseKey
                drawFixation(window, cx, cy, fixSize,p.dimFactor*p.fixColor*white);
                Screen('Flip',window);
            end

            %% Percent correct for orientation

            if testStatus~=0 && (strcmp('9(',targetResponseKeyName) || strcmp('0)',targetResponseKeyName)) % response CCW
                if testOrientation==2 % stimuli CW
                    correct = 1;  % stimuli = response
                elseif testOrientation==1 % stimuli CCW
                    correct = 0; % stimuli != response
                end
            elseif testStatus~=0 && (strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName)) % response CCW

                if testOrientation==1 %stimuli CCW
                    correct = 1; % stimuli = response
                elseif testOrientation==2 %stimuli CW
                    correct = 0; % stimuli != response
                end
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            %% Response Acceptance Marker based on Correctness of Response
            if targetResponseKey && correct==1
                drawFixation(window, cx, cy, fixSize,[0 1 0]*p.fixColor*white);    
                Screen('Flip',window);
     
            elseif targetResponseKey && correct==0
                drawFixation(window, cx, cy, fixSize,[1 0 0]*p.fixColor*white);
                Screen('Flip',window);
            end

            timeEnd=GetSecs();

            %% Store trial info
            trials(trialIdx, rtIdx) = targetRT;
            trials(trialIdx, responseKeyIdx) = targetResponseKey;
            p.trials=trials;

            tic
            % save data
            d.block(iTrial) = block; %block trial is in
            d.timeStart(iTrial)=timeStart; % time at which trial starts
            d.timeEnd(iTrial)=timeEnd; % time at which trial ends
            d.timeElapsed(iTrial)=timeEnd-timeStart; % trial duration
            d.testStatus(iTrial) = testStatus; %absent or present
            d.testOrientation(iTrial) = testOrientation; %orientation of test stimuli
            d.testPhase(iTrial) = testPhase; %phase of test stimuli
            d.testContrast(iTrial) = testContrast; %contrast of test stimuli
            d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
            d.correct(iTrial) = correct; % response correctness (0/1)
            d.timeTargetResponse(iTrial)=timeTargetResponse; % time at which response was made
            d.targetRT(iTrial) = targetRT; %response time duration
            d.targetResponseKey(iTrial) = targetResponseKey; % response key for trial
            d.timePreStart(iTrial)=timePreStart; % time before standard stimuli presentation
            d.timeFix(iTrial)=timeFix; % time for fixation
            d.timeS(iTrial)=timeS; % time for standard
            d.timeBlank1(iTrial)=timeBlank1; % time after standard
            d.timeTone(iTrial)=timeTone; % time at tone
            d.timeT(iTrial)=timeT; % time at test

            save(sprintf('%s/%s_demo_s%d_pred_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,datestr(now, 'mmdd')),'d','p');

            %time between trials after response
            d.timeSpentSaving(iTrial)=toc;
            WaitSecs(p.ITI-d.timeSpentSaving(iTrial));

            if mod(iTrial,p.BlockTrials)==0 || iTrial>p.BlockTrials*block
                % Calculate block accuracy
                blockStartTrial = (iTrial/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                if blockStartTrial < 0 % we are doing less than one block
                    blockStartTrial = 1;
                end

                current=d.correct(blockStartTrial:iTrial);
                nan_correctnessIdx= isnan(current);
                non_nan_correctness=current(~nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f / %.2f (block/total)', ...
                    100*blockCorrectness, 100*mean(d.correct));

                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, ceil(nTrials/p.BlockTrials));
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
                DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
                Screen('Flip', window);
                WaitSecs(1);
                if iTrial < nTrials
                    KbWait(-1);
                end

                block = block+1;
            end
        end
    case 3
        %% 3 TASK version 2
        %% %%%% Generate trials in different conditions %%%%
        % 43008 trial version 
        trials_headers = {'precueValidity','gratingOrientation','gratingPhase', 'gratingContrast','plaidStatus','plaidOrientation','plaidPhase' 'plaidContrast','precue','responseKey','response','accuracy','rt'};
        %make sure column indices match trials headers
        precueValidityIdx = strcmp(trials_headers,'precueValidity');
        gratingOrientationIdx = strcmp(trials_headers,'gratingOrientation'); % 0 (low) or 1 (high)
        gratingPhaseIdx = strcmp(trials_headers,'gratingPhase'); % test phase value
        gratingContrastIdx = strcmp(trials_headers,'gratingContrast'); % test contrast

        plaidStatusIdx = strcmp(trials_headers,'plaidStatus'); % plaid status absent (1) or present (2)
        plaidOrientationIdx = strcmp(trials_headers,'plaidOrientation'); % 0 (low) or 1 (high)
        plaidPhaseIdx = strcmp(trials_headers,'plaidPhase'); % test phase value
        plaidContrastIdx = strcmp(trials_headers,'plaidContrast'); % test contrast

        testStatusIdx = strcmp(trials_headers,'tStatus'); % test status absent or present

        precueIdx=strcmp(trials_headers,'precue'); % precue
        responseKeyIdx = strcmp(trials_headers,'responseKey');
        responseIdx = strcmp(trials_headers,'response');
        accuracyIdx = strcmp(trials_headers,'accuracy');
        rtIdx = strcmp(trials_headers,'rt');


        trials = fullfact([numel(p.precueValidities),... % 1
            numel(p.gratingOrientations),... % 2 grating ori
            numel(p.testPhases),... % 3 grating phase 
            numel(p.gratingContrasts)... % 4 grating contrast
            numel(p.plaidStatus)... % 5 plaid status
            numel(p.plaidOrientations),... % 6 plaid ori
            numel(p.testPhases),... % 7 plaid phase 
            numel(p.plaidContrasts1)]); % 8 plaid contrast
       
        if (p.debug=="N" || p.debug=="n") && p.reps==1
            %trials=repmat(trials,p.reps*p.repScale1,1);
            nTrials = size(trials,1); 
            nBlocks=nTrials/p.BlockTrials; 
          
        elseif (p.debug=="N" || p.debug=="n") && p.reps==2
            trials=repmat(trials,p.reps*p.repScale2,1); 
            nTrials = size(trials,1);
            nTrialsNoSkip=nTrials;
            nBlocks=p.trialsNeeded/nTrials; 
        end

        trialOrder = randperm(nTrials);

        %% Eyetracker
        if p.eyeTracking
            % Initialize eye tracker
            [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
            if exitFlag
                return
            end

            % Write subject ID into the edf file
            Eyelink('message', 'BEGIN DESCRIPTIONS');
            Eyelink('message', 'Subject code: %s', p.subjectID);
            Eyelink('message', 'END DESCRIPTIONS');

            % No sounds indicating success of calibration
            %     el.targetbeep = false;
            %     el.calibration_failed_beep = [0 0 0];
            %     el.calibration_success_beep = [0 0 0];
            el.drift_correction_target_beep = [0 0 0];
            el.drift_correction_failed_beep = [0 0 0];
            el.drift_correction_success_beep = [0 0 0];

            % Accept input from all keyboards
            el.devicenumber = -1; %see KbCheck for details of this value

            % Update with custom settings
            EyelinkUpdateDefaults(el);

            % Calibrate eye tracker
            [~, exitFlag] = rd_eyeLink('calibrate', window, el);
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

        %% Show instruction screen and wait for a button press
        instructions = 'This is the second version of the experiment\n\n';
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a tone followed by a patch.\n\n You will have to report the orientation of the patch.\n\n Your goal is to determine whether the patch is: \n\n counterclockwise (press 9)  \n\n or clockwise (press 0)! \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 =' \n\n At random, a plaid patch will appear instead. \n\n For those trials you must report which of the two orientations in the plaid is stronger in contrast. \n\n For counterclockwise (press 9).  \n\n For clockwise (press 0)!  Press to start!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);

        KbWait(devNum);
        timeStart = GetSecs;
        correct = [];
        block=1;

        skippedTrials = [];
        iTrial=1;
        completedTrials=0;

        disp('reached')

        while iTrial<=nTrials
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
            stopThisTrial = 0;

            %% %%%% PRESENT A TRIAL %%%
            %% Get condition information for this trial
            precueValidity = p.precueValidities(trials(trialIdx, precueValidityIdx));
            gratingOrientation = trials(trialIdx, gratingOrientationIdx);
            gratingPhase = trials(trialIdx, gratingPhaseIdx);
            gratingContrast = trials(trialIdx, gratingContrastIdx);
            
            plaidStatus=p.plaidStatus(trials(trialIdx, plaidStatusIdx));

            plaidOrientation = trials(trialIdx, plaidOrientationIdx);
            plaidPhase = trials(trialIdx, plaidPhaseIdx);
            plaidContrast = trials(trialIdx, plaidContrastIdx);
            
            
            % orientation
            gOrientation=p.gratingOrientations(gratingOrientation);
            pOrientation=p.plaidOrientations(plaidOrientation);
            % tone
            toneName = p.precueNames{precueValidity};
            switch toneName
                case 'valid'
                    if plaidStatus==1
                        tone=cueTones(gratingOrientation,:);
                        toneVersion=gratingOrientation;
                    elseif plaidStatus==2
                        tone=cueTones(plaidOrientation,:);
                        toneVersion=plaidOrientation;
                    end
                case 'invalid'
                    if plaidStatus==1 && gratingOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==1 && gratingOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
 
                    elseif plaidStatus==2 && plaidOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==2 && plaidOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
                    end
                otherwise
                    error('precueName not recognized')
            end

            %% Store trial information
            
            d.gratingOrientation(iTrial) = gratingOrientation; %orientation of test stimuli
            d.gratingPhase(iTrial) = gratingPhase; %phase of test stimuli
            d.gratingContrast(iTrial) = gratingContrast; %contrast of test stimuli
            d.plaidPhase(iTrial) = plaidPhase; %phase of test stimuli
            d.plaidContrast(iTrial) = plaidContrast; %contrast of test stimuli
            d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
            d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)

            %% Store stimulus information in trials matrix

            trials(trialIdx, precueIdx) = toneVersion;
            %% %%%% Play the trial %%%%
            %% Present fixation rest (grey)
            drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
            timePreStart = Screen('Flip', window);
            %% Present fixation active (white)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
            if p.eyeTracking
                Eyelink('Message', 'FixOn')
            end
            %% Check fixation hold
            if p.eyeTracking
                driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
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
                while GetSecs < timeTone + p.standSOA - p.eyeSlack && ~stopThisTrial
                    WaitSecs(.01);
                    fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                    fixations = [fixations fixation];

                    if fixation==0
                        stopThisTrial = 1;
                        WaitSecs(1);
                    
                        % redo this trial at the end of the experiment
                        % this can be easily done by appending the trial number to the end of
                        % trialOrder
                        trialOrder(end+1) = trialOrder(iTrial);
                        d.stopThisTrial(iTrial) = stopThisTrial;
                        skippedTrials(end+1) = trialOrder(iTrial);
                        d.timeTargetResponse(iTrial) = NaN;
                        d.correct(iTrial) = NaN;
                        nTrials = nTrials + 1;
                        iTrial = iTrial + 1;

                        Screen('FillRect', window, white*p.backgroundColor);
                        Screen('Flip', window);
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
            
            if plaidStatus==2
                Screen('DrawTexture', window, tex_plaid{plaidContrast,plaidPhase}, [], imRect, pOrientation);
            elseif plaidStatus==1
                Screen('DrawTexture', window, tex{gratingContrast,gratingPhase }, [], imRect, gOrientation);
            end
            
            drawFixation(window, cx, cy, fixSize, p.fixColor*white)
            
            timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard

            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
         
            %% Eyetracking 
            if p.eyeTracking
                Eyelink('Message', 'EVENT_RESPCUE');
            end
          
            %% Wait for response
            % check only valid response keys
            targetResponseKey = [];
            while isempty(targetResponseKey)
                [timeTargetResponse, keyCode] = KbWait(devNum);
                targetRT = timeTargetResponse - timeS;
                targetResponseKey = find(ismember(validKeys,find(keyCode)));
                targetResponseKeyName = KbName(validKeys(targetResponseKey));
                correct = NaN;
            end
            if isempty(targetResponseKey)
                targetRT = NaN; % timeout
                targetResponseKey = NaN;
                targetResponseKeyName = NaN;
                correct = NaN;
            end
            if p.eyeTracking
                Eyelink('Message', 'TRIAL_END');
            end

            %% Response collected (blue fixation)
            if targetResponseKey
                drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
                Screen('Flip',window);
            
            end

            %% Percent correct for orientation

            if plaidStatus~=2 && (strcmp('9(',targetResponseKeyName) || strcmp('0)',targetResponseKeyName)) % response CCW
                if gratingOrientation==2 % stimuli CW
                    correct = 1;  % stimuli = response
                elseif gratingOrientation==1 % stimuli CCW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus~=2 &&(strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName)) % response CCW
                if gratingOrientation==1 %stimuli CCW
                    correct = 1; % stimuli = response
                elseif gratingOrientation==2 %stimuli CW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus==2
                correct=NaN;
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            timeEnd=GetSecs();

            %% Store trial info
            trials(trialIdx, rtIdx) = targetRT;
            trials(trialIdx, responseKeyIdx) = targetResponseKey;
            p.trials=trials;

            tic
            % save data
            d.block(iTrial) = block; %block trial is in
            d.timeStart(iTrial)=timeStart; % time at which trial starts
            d.timeEnd(iTrial)=timeEnd; % time at which trial ends
            d.timeElapsed(iTrial)=timeEnd-timeStart; % trial duration
            d.correct(iTrial) = correct; % response correctness (0/1)
            d.timeTargetResponse(iTrial)=timeTargetResponse; % time at which response was made
            d.targetRT(iTrial) = targetRT; %response time duration
            d.targetResponseKey(iTrial) = targetResponseKey; % response key for trial
            d.timePreStart(iTrial)=timePreStart; % time before standard stimuli presentation
            d.timeFix(iTrial)=timeFix; % time for fixation
            d.timeS(iTrial)=timeS; % time for standard
            d.timeBlank1(iTrial)=timeBlank1; % time after standard
            d.timeTone(iTrial)=timeTone; % time at tone
            d.trialOrder=trialOrder;
            save(sprintf('%s/%s_s%d_predv2_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

            %time between trials after response
            d.timeSpentSaving(iTrial)=toc;
            WaitSecs(p.ITI-d.timeSpentSaving(iTrial));
            
            completedTrials = completedTrials + 1; 
            %% SANTIY CHECK
            fprintf('Trial %d/%d in block %d, trial %d of %d total \n', completedTrials, p.BlockTrials, block,iTrial,nTrials);
            %% Blocking
            %if completedTrials==p.BlockTrials && (mod(iTrial-completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block) 
            if (mod(completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block)  
                % Calculate block accuracy
                %blockStartTrial = (iTrial/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                blockStartTrial = ((iTrial-completedTrials)/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                if blockStartTrial < 0 % we are doing less than one block
                    blockStartTrial = 1;
                end

                current=d.correct(blockStartTrial:iTrial);
                nan_correctnessIdx= isnan(current);
                %non_nan_correctness=d.correct(non_nan_correctnessIdx);
                non_nan_correctness=current(~nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);

                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
                DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
                %         DrawFormattedText(window, cueText, 'center', cy-3*textOffset, white);
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

    case 4 
        % 600-ish trial version

        %% Grating trials
        trials_headers1 = {'precueValidity','gratingOrientation','gratingPhase', 'gratingContrast','precue','responseKey','response','accuracy','rt'};
        % make sure column indices match trials headers
        precueValidityIdx1 = strcmp(trials_headers1,'precueValidity'); %precue validity index 1 2 3 (for 1) 4 (for 2)
        gratingOrientationIdx = strcmp(trials_headers1,'gratingOrientation'); % grating orientation index 1 (-45) or 2 (+45)
        gratingPhaseIdx = strcmp(trials_headers1,'gratingPhase'); % test phase value 
        gratingContrastIdx = strcmp(trials_headers1,'gratingContrast'); % test contrast
        
        precueIdx=strcmp(trials_headers1,'precue'); % precue
        responseKeyIdx = strcmp(trials_headers1,'responseKey');
        responseIdx = strcmp(trials_headers1,'response');
        accuracyIdx = strcmp(trials_headers1,'accuracy');
        rtIdx = strcmp(trials_headers1,'rt');
        
        trials1 = fullfact([numel(p.precueValidities),... % 1 precue valdities
             numel(p.gratingOrientations),... % 2 grating ori
             numel(p.testPhases),... % 3 grating phase 
             numel(p.gratingContrasts) ]); % 4 grating contrast
       
        %% Plaid (waffle) Trials
        trials_headers2 = {'precueValidity','plaidOrientation','plaidPhase' 'plaidContrast'};
        % make sure column indices match trials headers
        precueValidityIdx2 = strcmp(trials_headers2,'precueValidity');
        plaidOrientationIdx = strcmp(trials_headers2,'plaidOrientation'); % 0 (low) or 1 (high)
        plaidPhaseIdx = strcmp(trials_headers2,'plaidPhase'); % test phase value
        plaidContrastIdx = strcmp(trials_headers2,'plaidContrast'); % test contrast
       
        trials2 = fullfact([ numel(p.precueValidities)...
             numel(p.plaidOrientations),... % 1 plaid ori
             numel(p.testPhases),... % 2 plaid phase 
             numel(p.plaidContrasts1)]); % 3 plaid contrast

        if (p.debug=="N" || p.debug=="n") % && p.reps==1
            %trials=repmat(trials,p.reps*p.repScale1,1);
            trials1=repmat(trials1,3,1);
            nTrials = size(trials1,1)+size(trials2,1); 
            nBlocks=nTrials/p.BlockTrials; 
        end

        trialOrder = randperm(nTrials);

        %% Eyetracker
        if p.eyeTracking
            % Initialize eye tracker
            [el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
            if exitFlag
                return
            end

            % Write subject ID into the edf file
            Eyelink('message', 'BEGIN DESCRIPTIONS');
            Eyelink('message', 'Subject code: %s', p.subjectID);
            Eyelink('message', 'END DESCRIPTIONS');

            % No sounds indicating success of calibration
            %     el.targetbeep = false;
            %     el.calibration_failed_beep = [0 0 0];
            %     el.calibration_success_beep = [0 0 0];
            el.drift_correction_target_beep = [0 0 0];
            el.drift_correction_failed_beep = [0 0 0];
            el.drift_correction_success_beep = [0 0 0];

            % Accept input from all keyboards
            el.devicenumber = -1; %see KbCheck for details of this value

            % Update with custom settings
            EyelinkUpdateDefaults(el);

            % Calibrate eye tracker
            [~, exitFlag] = rd_eyeLink('calibrate', window, el);
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

        %% Show instruction screen and wait for a button press
        instructions = 'This is the second version of the experiment\n\n';
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a tone followed by a patch.\n\n You will have to report the orientation of the patch.\n\n Your goal is to determine whether the patch is: \n\n counterclockwise (press 9)  \n\n or clockwise (press 0)! \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 =' \n\n At random, a plaid patch will appear instead. \n\n For those trials you must report which of the two orientations in the plaid is stronger in contrast. \n\n For counterclockwise (press 9).  \n\n For clockwise (press 0)!  Press to start!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);

        KbWait(devNum);
        timeStart = GetSecs;
        correct = [];
        block=1;

        skippedTrials = [];
        iTrial=1;
        completedTrials=0;

        disp('reached')

        %for iTrial = 1:nTrials % the iteration in the trial loop
        while iTrial<=nTrials
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
            stopThisTrial = 0;

            %% %%%% PRESENT A TRIAL %%%
            %% Get condition information for this trial

            if trialIdx<=size(trials1,1)
                plaidStatus=1;
                precueValidity = p.precueValidities(trials1(trialIdx, precueValidityIdx1));
                gratingOrientation = trials1(trialIdx, gratingOrientationIdx);
                gratingPhase = trials1(trialIdx, gratingPhaseIdx);
                gratingContrast = trials1(trialIdx, gratingContrastIdx);
                gOrientation=p.gratingOrientations(gratingOrientation);
            elseif trialIdx>size(trials1,1)
                plaidStatus=2;
                precueValidity = p.precueValidities(trials2(trialIdx-size(trials1,1), precueValidityIdx2));
                plaidOrientation = trials2(trialIdx-size(trials1,1), plaidOrientationIdx);
                plaidPhase = trials2(trialIdx-size(trials1,1), plaidPhaseIdx);
                plaidContrast = trials2(trialIdx-size(trials1,1), plaidContrastIdx);
                pOrientation=p.plaidOrientations(plaidOrientation);
            end     
            
            % tone
            toneName = p.precueNames{precueValidity};
            switch toneName
                case 'valid'
                    if plaidStatus==1
                        tone=cueTones(gratingOrientation,:);
                        toneVersion=gratingOrientation;
                    elseif plaidStatus==2
                        tone=cueTones(plaidOrientation,:);
                        toneVersion=plaidOrientation;
                    end
                case 'invalid'
                    if plaidStatus==1 && gratingOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==1 && gratingOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
 
                    elseif plaidStatus==2 && plaidOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==2 && plaidOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
                    end
                otherwise
                    error('precueName not recognized')
            end

            %% Store trial information
             if plaidStatus==1
                d.gratingOrientation(iTrial) = gratingOrientation; %orientation of test stimuli
                d.gratingPhase(iTrial) = gratingPhase; %phase of test stimuli
                d.gratingContrast(iTrial) = gratingContrast; %contrast of test stimuli
                d.plaidOrientation(iTrial) = NaN; %orientation of plaid stimuli
                d.plaidPhase(iTrial) = NaN; %phase of test stimuli
                d.plaidContrast(iTrial) = NaN; %contrast of test stimuli
                d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
                d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
            elseif plaidStatus==2
                d.plaidOrientation(iTrial) = plaidOrientation; %orientation of plaid stimuli
                d.plaidPhase(iTrial) = plaidPhase; %phase of test stimuli
                d.plaidContrast(iTrial) = plaidContrast; %contrast of test stimuli
                d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
                d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
                d.gratingOrientation(iTrial) = NaN; %orientation of test stimuli
                d.gratingPhase(iTrial) = NaN; %phase of test stimuli
                d.gratingContrast(iTrial) = NaN; %contrast of test stimuli
            end
            %% Store stimulus information in trials matrix

            %trials(trialIdx, precueIdx) = toneVersion;
            %% %%%% Play the trial %%%%
            %% Present fixation rest (grey)
            drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
            timePreStart = Screen('Flip', window);
            %% Present fixation active (white)
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
            if p.eyeTracking
                Eyelink('Message', 'FixOn')
            end
            %% Check fixation hold
            if p.eyeTracking
                driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
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
                while GetSecs < timeTone + p.standSOA - p.eyeSlack && ~stopThisTrial
                    WaitSecs(.01);
                    fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                    fixations = [fixations fixation];

                    if fixation==0
                        stopThisTrial = 1;
                        WaitSecs(1);
                    
                        % redo this trial at the end of the experiment
                        % this can be easily done by appending the trial number to the end of
                        % trialOrder
                        trialOrder(end+1) = trialOrder(iTrial);
                        d.stopThisTrial(iTrial) = stopThisTrial;
                        skippedTrials(end+1) = trialOrder(iTrial);
                        d.timeTargetResponse(iTrial) = NaN;
                        d.correct(iTrial) = NaN;
                        nTrials = nTrials + 1;
                        iTrial = iTrial + 1;

                        Screen('FillRect', window, white*p.backgroundColor);
                        Screen('Flip', window);
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
            
            if plaidStatus==2
                Screen('DrawTexture', window, tex_plaid{plaidContrast,plaidPhase}, [], imRect, pOrientation);
            elseif plaidStatus==1
                Screen('DrawTexture', window, tex{gratingContrast,gratingPhase }, [], imRect, gOrientation);
            end
            
            drawFixation(window, cx, cy, fixSize, p.fixColor*white)
            
            timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard

            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
         
            %% Eyetracking 
            if p.eyeTracking
                Eyelink('Message', 'EVENT_RESPCUE');
            end
          
            %% Wait for response
            % check only valid response keys
            targetResponseKey = [];
            while isempty(targetResponseKey)
                [timeTargetResponse, keyCode] = KbWait(devNum);
                targetRT = timeTargetResponse - timeS;
                targetResponseKey = find(ismember(validKeys,find(keyCode)));
                targetResponseKeyName = KbName(validKeys(targetResponseKey));
                correct = NaN;
            end
            if isempty(targetResponseKey)
                targetRT = NaN; % timeout
                targetResponseKey = NaN;
                targetResponseKeyName = NaN;
                correct = NaN;
            end
            if p.eyeTracking
                Eyelink('Message', 'TRIAL_END');
            end

            %% Response collected (blue fixation)
            if targetResponseKey
                drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
                Screen('Flip',window);
            
            end

            %% Percent correct for orientation

            if plaidStatus~=2 && (strcmp('9(',targetResponseKeyName) || strcmp('0)',targetResponseKeyName)) % response CCW
                if gratingOrientation==2 % stimuli CW
                    correct = 1;  % stimuli = response
                elseif gratingOrientation==1 % stimuli CCW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus~=2 &&(strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName)) % response CCW
                if gratingOrientation==1 %stimuli CCW
                    correct = 1; % stimuli = response
                elseif gratingOrientation==2 %stimuli CW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus==2
                correct=NaN;
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            timeEnd=GetSecs();

            %% Store trial info
            % trials(trialIdx, rtIdx) = targetRT;
            % trials(trialIdx, responseKeyIdx) = targetResponseKey;
            % p.trials=trials;

            tic
            % save data
            d.block(iTrial) = block; %block trial is in
            d.timeStart(iTrial)=timeStart; % time at which trial starts
            d.timeEnd(iTrial)=timeEnd; % time at which trial ends
            d.timeElapsed(iTrial)=timeEnd-timeStart; % trial duration
            d.correct(iTrial) = correct; % response correctness (0/1)
            d.timeTargetResponse(iTrial)=timeTargetResponse; % time at which response was made
            d.targetRT(iTrial) = targetRT; %response time duration
            d.targetResponseKey(iTrial) = targetResponseKey; % response key for trial
            d.timePreStart(iTrial)=timePreStart; % time before standard stimuli presentation
            d.timeFix(iTrial)=timeFix; % time for fixation
            d.timeS(iTrial)=timeS; % time for standard
            d.timeBlank1(iTrial)=timeBlank1; % time after standard
            d.timeTone(iTrial)=timeTone; % time at tone
            d.trialOrder=trialOrder;
            save(sprintf('%s/%s_s%d_predv2_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

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

                current=d.correct(blockStartTrial:iTrial);
                nan_correctnessIdx= isnan(current);
                non_nan_correctness=current(~nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);

                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
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

    case 5 
        %% Kok 2012 + Waffles interleaved

        %% Grating Trials
        trials_headers1 = {'precueValidity','staticGrating','Orientation','standardPhase','standardSPF','gratingSPF','differenceMultiplier','precue','responseKey','response','accuracy','rt'};
        % make sure column indices match trials headers
        precueValidityIdx1 = strcmp(trials_headers1,'precueValidity'); %precue validity
        staticGratingIdx=strcmp(trials_headers1,'staticGrating'); % which grating has the constant orientation of +/- 45
        OrientationIdx=strcmp(trials_headers1,'Orientation'); % orientation
        standardPhaseIdx=strcmp(trials_headers1,'standardPhase'); %stand grating phase
        standardSPFIdx=strcmp(trials_headers1,'standardSPF'); % standard grating SPF
        gratingSPFIdx=strcmp(trials_headers1,'gratingSPF'); % test grating SPF
        differenceMultiplierIdx=strcmp(trials_headers1,'differenceMultiplier'); % -1 or +1, sign of difference between 2 grating orientations

        precueIdx=strcmp(trials_headers1,'precue'); % precue
        responseKeyIdx = strcmp(trials_headers1,'responseKey');
        responseIdx = strcmp(trials_headers1,'response');
        accuracyIdx = strcmp(trials_headers1,'accuracy');
        rtIdx = strcmp(trials_headers1,'rt');
        
        trials1 = fullfact([numel(p.precueValidities),... %1 precue
             numel(p.staticGrating),... % 2 static grating
             numel(p.gratingOrientations),... % 3 standard orientation
             numel(p.testPhases),... % 4 standard phase 
             numel(p.gratingSPF)... % 5 standard SPF
             numel(p.gratingSPF)...% 6 grating SPF
             numel(p.differenceMultiplier)]); % 7 grating diff multiplier
        
        %% Plaid (waffles) trials
        trials_headers2 = {'precueValidity','plaidOrientation','plaidPhase' 'plaidContrast'};
        % make sure column indices match trials headers
        precueValidityIdx2 = strcmp(trials_headers2,'precueValidity'); %precue validity for plaid
        plaidOrientationIdx = strcmp(trials_headers2,'plaidOrientation'); % plaid orientation
        plaidPhaseIdx = strcmp(trials_headers2,'plaidPhase'); % plaid phase value
        plaidContrastIdx = strcmp(trials_headers2,'plaidContrast'); % plaid contrast
       
        trials2 = fullfact([ numel(p.precueValidities)... % 1 precue validity
             numel(p.plaidOrientations),... % 2 plaid ori
             numel(p.testPhases),... % 3 plaid phase 
             numel(p.plaidContrasts1)]); % 4 plaid contrast
        %% Merge trial count
        if (p.debug=="N" || p.debug=="n") % && p.reps==1
            nTrials = size(trials1,1)+size(trials2,1);  % total trials = number of grating trials + number of waffle trials
            %nBlocks=nTrials/p.BlockTrials; 
            nBlocks=nTrials/p.BlockTrials; 
        end

        trialOrder = randperm(nTrials); %randomize trial order
        
        %% Eyetracker
        if p.eyeTracking
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
            [~, exitFlag] = rd_eyeLink('calibrate', window, el);
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

        %% Show instruction screen and wait for a button press
        instructions = 'This is the second version of the experiment\n\n';
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a tone followed by 2 gratings.\n\n You will have to report the orientation of the second grating relative to the first.\n\n Your goal is to determine whether the grating is: \n\n counterclockwise (press 9)  \n\n or clockwise (press 0)! \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 =' \n\n At random, a plaid patch will appear instead. \n\n For those trials you must report which of the two orientations in the plaid is stronger in contrast. \n\n For counterclockwise (press 9).  \n\n For clockwise (press 0)! \n\n  Press to start!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack);
        WaitSecs(1);

        KbWait(devNum);
        timeStart = GetSecs;
        correct = [];
        block=1;

        skippedTrials = [];
        iTrial=1;
        completedTrials=0;

        disp('reached')
        firstNonWaffle=0; % tracking the first trial that is a non-waffle trial

        while iTrial<=nTrials
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
            stopThisTrial = 0;

            %% Get condition information for this trial

            if trialIdx<=size(trials1,1)

                plaidStatus=1; % if the trial id is a value <= the number of grating trials, this trial will be a grating trial
                
                differenceMultiplier=p.differenceMultiplier(trials1(trialIdx, differenceMultiplierIdx)); % get difference multiplier
                precueValidity = p.precueValidities(trials1(trialIdx, precueValidityIdx1)); % get precue validity (1: valid, 2: invalid)
                staticGrating=trials1(trialIdx, staticGratingIdx); % find out which of the gratings has the +/-45 "static" orientation
                Orientation = trials1(trialIdx, OrientationIdx); % get the orientation index
                standardPhase = trials1(trialIdx, standardPhaseIdx); % get the index for the phase of the standard grating
                gratingPhase=standardPhase; % the index of the phase of the test is the same as the index of the standard
                standardSPF= trials1(trialIdx, standardSPFIdx); % get index of standard SPF
                gratingSPF=trials1(trialIdx, gratingSPFIdx); % get index of test SPF
        
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial
                    stairIdx=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
                    lastFewAcc=[]; % no previous grating trials have occured so this is empty
                    firstNonWaffle=1; % firstNonWaffle becomes 1 because this is officially the first grating trial
                elseif firstNonWaffle==1 % if this variable is 1, this trial is not the first grating trial
                    corrects=d.correct; %get all corrects
                    skipRowsCorrect=isnan(corrects); %find the NaN values
                    corrects(skipRowsCorrect)=[];% delete the NaN values associated with the waffle trials
                    [stairIdx lastFewAcc]=updateStaircase(p.stairs, stairIdx, corrects, correct); % get new stair index value
                end

                if staticGrating==1 % standard is +/-45 and test is slightly different
                    sOrientation=p.gratingOrientations(Orientation);
                    gOrientation=sOrientation+p.stairs(stairIdx)*differenceMultiplier;
                elseif staticGrating==2 %test is +/-45 and standard is slightly different
                    gOrientation=p.gratingOrientations(Orientation);
                    sOrientation=gOrientation+p.stairs(stairIdx)*differenceMultiplier;
                end

            elseif trialIdx>size(trials1,1) 
                plaidStatus=2; % if the trial id is a value > the number of grating trials, this trial will be a waffle trial
                precueValidity = p.precueValidities(trials2(trialIdx-size(trials1,1), precueValidityIdx2));
                plaidOrientation = trials2(trialIdx-size(trials1,1), plaidOrientationIdx);
                plaidPhase = trials2(trialIdx-size(trials1,1), plaidPhaseIdx);
                plaidContrast = trials2(trialIdx-size(trials1,1), plaidContrastIdx);
                pOrientation=p.plaidOrientations(plaidOrientation); %get the orientation value (+/-45) using the plaid orientation index
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial
                    stairIdx=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
            
                end
            end
            % tone
            toneName = p.precueNames{precueValidity};
            switch toneName
                case 'valid'
                    if plaidStatus==1
                        tone=cueTones(Orientation,:);
                        toneVersion=Orientation;
                    elseif plaidStatus==2
                        tone=cueTones(plaidOrientation,:);
                        toneVersion=plaidOrientation;
                    end
                case 'invalid'
                    if plaidStatus==1 && Orientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==1 && Orientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
 
                    elseif plaidStatus==2 && plaidOrientation==1
                        tone = cueTones(2,:);
                        toneVersion=2;
                    elseif plaidStatus==2 && plaidOrientation==2
                        tone = cueTones(1,:);
                        toneVersion=1;
                    end
                otherwise
                    error('precueName not recognized')
            end

            %% Set up test mask for this trial

            %% Store trial information
             if plaidStatus==1 % if this is a grating trial
                d.stairIdx(iTrial) = stairIdx; %store stair index stimuli
                d.stair(iTrial)= p.stairs(stairIdx); %store stair value
                d.differenceMultiplier(iTrial)=differenceMultiplier; %store difference multiplier
                d.staticGrating(iTrial)=staticGrating; %store the static grating index
                d.gratingContrast(iTrial) = p.gratingContrast2; %store contrast of test stimuli
                d.gratingSPF(iTrial) = gratingSPF; %store SPF of test stimuli

                d.orientation(iTrial) = Orientation; %store orientation of trial
                d.standardPhase(iTrial) = standardPhase; %store phase of standard stimuli
                d.standardContrast(iTrial) = p.standardContrast2; % store contrast of standard stimuli
                d.standardSPF(iTrial) = standardSPF; %store SPF of standard stimuli
                d.precueValidity(iTrial) = precueValidity; % store cue validity (1:valid, 2: invalid)
                
                d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
                d.plaidOrientation(iTrial) = NaN; % keep orientation of plaid stimuli as NaN
                d.plaidPhase(iTrial) = NaN; %keep phase of plaid stimuli as NaN
                d.plaidContrast(iTrial) = NaN; %keep contrast of plaid stimuli as NaN

             elseif plaidStatus==2 %if this is a waffle trial
                d.plaidOrientation(iTrial) = plaidOrientation; %store orientation of plaid stimuli
                d.plaidPhase(iTrial) = plaidPhase; %store phase of plaid stimuli
                d.plaidContrast(iTrial) = plaidContrast; %store contrast of plaid stimuli
                d.plaidStatus(iTrial) = plaidStatus; %is it a plaid this time?
                d.precueValidity(iTrial) = precueValidity; % cue validity (valid/invalid)
                d.stairIdx(iTrial) = stairIdx;
                d.gratingOrientation(iTrial) = NaN; %orientation of test stimuli
                d.gratingPhase(iTrial) = NaN; %phase of test stimuli
                d.gratingContrast(iTrial) = NaN; %contrast of test stimuli
                d.standardOrientation(iTrial) = NaN; %orientation of test stimuli
                d.standardPhase(iTrial) = NaN; %phase of test stimuli
                d.standardContrast(iTrial) = NaN; %contrast of test stimuli
                d.standardSPF(iTrial) = NaN; %contrast of test stimuli
            end
            %% Store stimulus information in trials matrix

            %trials(trialIdx, precueIdx) = toneVersion;
            %% %%%% Play the trial %%%%

            if plaidStatus==1
                %% Present fixation rest (grey)
                drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
                timePreStart = Screen('Flip', window);
                %% Present fixation active (white)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
                if p.eyeTracking
                    Eyelink('Message', 'FixOn')
                end
                %% Check fixation hold
                if p.eyeTracking
                    driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
                    if driftCorrected
                        % restart trial
                        drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                        timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
                    end
                end

                %% Present predictive tone
                PsychPortAudio('FillBuffer', pahandle, tone);
                timeTone = PsychPortAudio('Start', pahandle, [], timeFix+p.signalStart, 1); % waitForStart = 1 in order to return a timestamp of playback
                
               
                %% Present STANDARD
                Screen('DrawTexture', window, tex_stand{standardPhase,standardSPF}, [], imRect, sOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
                
                % blank
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
                
                %% EYE TRACKING
    
                if p.eyeTracking
                    Eyelink('Message', 'EVENT_CUE');
                end
    
                if p.eyeTracking
                    while GetSecs < timeTone + p.standSOA - p.eyeSlack && ~stopThisTrial
                        WaitSecs(.01);
                        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                        fixations = [fixations fixation];
    
                        if fixation==0
                            stopThisTrial = 1;
                            WaitSecs(1);
                        
                            % redo this trial at the end of the experiment
                            % this can be easily done by appending the trial number to the end of
                            % trialOrder
                            trialOrder(end+1) = trialOrder(iTrial);
                            d.stopThisTrial(iTrial) = stopThisTrial;
                            skippedTrials(end+1) = trialOrder(iTrial);
                            d.timeTargetResponse(iTrial) = NaN;
                            d.correct(iTrial) = NaN;
                            nTrials = nTrials + 1;
                            iTrial = iTrial + 1;
    
                            Screen('FillRect', window, white*p.backgroundColor);
                            Screen('Flip', window);
                            WaitSecs(1);
                        else 
                            stopThisTrial = 0;
                        end
                    end

                    if stopThisTrial
                        continue
                    end
                end
                %% Present TEST
                Screen('DrawTexture', window, tex{gratingPhase,gratingSPF}, [], imRect, gOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                
                timeT = Screen('Flip', window, timeS + p.standSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
    
                % blank
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank2 = Screen('Flip', window, timeT + p.imDur - slack); %timeS + how long i want stimulus to be presented for
             
                %% Eyetracking 
                if p.eyeTracking
                    Eyelink('Message', 'EVENT_RESPCUE');
                end
            elseif plaidStatus==2
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
                    driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
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
                    while GetSecs < timeTone + p.standSOA - p.eyeSlack && ~stopThisTrial
                        WaitSecs(.01);
                        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
                        fixations = [fixations fixation];
    
                        if fixation==0
                            stopThisTrial = 1;
                            WaitSecs(1);
                        
                            % redo this trial at the end of the experiment
                            % this can be easily done by appending the trial number to the end of
                            % trialOrder
                            trialOrder(end+1) = trialOrder(iTrial);
                            d.stopThisTrial(iTrial) = stopThisTrial;
                            skippedTrials(end+1) = trialOrder(iTrial);
                            d.timeTargetResponse(iTrial) = NaN;
                            d.correct(iTrial) = NaN;
                            nTrials = nTrials + 1;
                            iTrial = iTrial + 1;
    
                            Screen('FillRect', window, white*p.backgroundColor);
                            Screen('Flip', window);
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
                
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
    
                % blank
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
             
                %% Eyetracking 
                if p.eyeTracking
                    Eyelink('Message', 'EVENT_RESPCUE');
                end
            end

            %% Wait for response
            % check only valid response keys
            targetResponseKey = [];
            while isempty(targetResponseKey)
                [timeTargetResponse, keyCode] = KbWait(devNum);
                if plaidStatus==2
                    targetRT = timeTargetResponse - timeS;
                elseif plaidStatus==1
                    targetRT = timeTargetResponse - timeT;
                end
                targetResponseKey = find(ismember(validKeys,find(keyCode)));
                targetResponseKeyName = KbName(validKeys(targetResponseKey));
                correct = NaN;
            end
            if isempty(targetResponseKey)
                targetRT = NaN; % timeout
                targetResponseKey = NaN;
                targetResponseKeyName = NaN;
                correct = NaN;
            end
            if p.eyeTracking
                Eyelink('Message', 'TRIAL_END');
            end

            %% Response collected (blue fixation)
            if targetResponseKey
                drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
                Screen('Flip',window);
            
            end

            %% Percent correct for orientation

            if plaidStatus==1 && strcmp('0)',targetResponseKeyName) % response CW
                if (differenceMultiplier==-1 && staticGrating==1)||(differenceMultiplier==1 && staticGrating==2) % stimuli CCW
                    correct = 0;  % stimuli != response
                elseif (differenceMultiplier==1 && staticGrating==1) ||(differenceMultiplier==-1 && staticGrating==2) % stimuli CW
                    correct = 1; % stimuli = response
                end
            elseif plaidStatus==1 && strcmp('9(',targetResponseKeyName)  % response CCW
                if (differenceMultiplier==-1 && staticGrating==1) ||(differenceMultiplier==1 && staticGrating==2)%stimuli CCW
                    correct = 1; % stimuli = response
                elseif (differenceMultiplier==1 && staticGrating==1) ||(differenceMultiplier==-1 && staticGrating==2)  9 %stimuli CW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus==2
                correct_p=NaN;
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            timeEnd=GetSecs();

            %% Keep track of last 3 trials
            if plaidStatus==1  
                lastFewAcc=[lastFewAcc correct];
            end

            %% Store trial info
            %trials(trialIdx, rtIdx) = targetRT;
            %trials(trialIdx, responseKeyIdx) = targetResponseKey;
            p.trials1=trials1;
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

            if plaidStatus==1
                d.timeT(iTrial)=timeT;
                d.timeBlank2(iTrial)=timeBlank2;
                d.correct(iTrial) = correct; % response correctness (0/1)
            elseif plaidStatus==2
                d.timeT(iTrial)=NaN;
                d.timeBlank2(iTrial)=NaN;
                d.correct(iTrial) = correct_p;
            end

            d.trialOrder=trialOrder;
            save(sprintf('%s/%s_s%d_predv2_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

            %time between trials after response
            d.timeSpentSaving(iTrial)=toc;
            WaitSecs(p.ITI-d.timeSpentSaving(iTrial));
            
            completedTrials = completedTrials + 1; 
            %% SANTIY CHECK
            fprintf('Trial %d/%d in block %d, trial %d of %d total \n', completedTrials, p.BlockTrials, block,iTrial,nTrials);
            %% Blocking
            %if completedTrials==p.BlockTrials && (mod(iTrial-completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block) 
            if (mod(completedTrials,p.BlockTrials)==0 || (iTrial-completedTrials)>p.BlockTrials*block)  
                % Calculate block accuracy
                blockStartTrial = ((iTrial-completedTrials)/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                if blockStartTrial < 0 % we are doing less than one block
                    blockStartTrial = 1;
                end

                current=d.correct(blockStartTrial:iTrial);
                nan_correctnessIdx= isnan(current);
                non_nan_correctness=current(~nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);

                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s', blockMessage, accMessage, keyMessage);
                DrawFormattedText(window, breakMessage, 'center', 'center', [1 1 1]*white);
                %         DrawFormattedText(window, cueText, 'center', cy-3*textOffset, white);
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

    otherwise
        error('Expt stage not found.')
end

if p.eyeTracking
    p.eyefixation = fixations;
    p.eyefixT1 = timeFix;
end

%% Questionnaire
% adapted from PTB3 GetEchoString
% maxNumChar==200;
% vLineSpacing=2;
%Screen('FillRect', window, white*p.backgroundColor);
%instructions = 'Did you notice anything about the experiment?';
%DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
% while true
%     if useKbCheck
%         char = GetKbChar(varargin{:});
%     else
%         char = GetChar;
%     end
%     if isempty(char)
%         string = '';
%         break;
%     end
%     switch (abs(char))
%         case {13, 3, 10}
%             % ctrl-C, enter, or return
%             break;
%         case 8
%             % backspace
%             if ~isempty(string)
%                 % Remove last character from string:
%                 string = string(1:length(string)-1);
%             end
%         otherwise
%             string = [string, char];
%     end
% 
%     output = [msg, ' ', string];
%     output=WrapString(output,maxNumChar);
%     DrawFormattedText(windowPtr,output,x,y,textColor,[],0,0,vLineSpacing);
%     Screen('Flip',windowPtr);
% end





%% Save eye data and shut down the eye tracker
if p.eyeTracking
    rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});

    % rename eye file
    eyeFileFull = sprintf(sprintf('%s/%s_s%d_eye.edf',data.eyeDataDir,p.subjectID,p.sessionNum),'file');
    copyfile(sprintf('%s/%s.edf', eyeDataDir, eyeFile), eyeFileFull)
end

%% Clean Up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')
ShowCursor;
end

