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
p.task = input(['Prediction task run:\n' ...
    '1 - Demo\n'...
    '2 - Task version Kok + Waffles\n']);
p.counter=input('Choose 1 or 2 for participant \n'); %counterbalance 
p.eyeTracking=input('Eyetracking (0/1)? ');

%% Setup
% Add paths
directory = pwd; % get project directory path, set to prediction folder parent level
addpath(genpath(directory))
% Get parameters
p = predParams(p);

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
if p.task==1 || p.task==2
    validKeys = KbName({'1!','2@','9(','0)'});
elseif p.task==3 || p.task==4 ||p.task==5
    validKeys = KbName({'9(','0)'});
end
%% Make TEST stimuli

imPos = round(p.imPos*pixelsPerDegree); % image position (pixels)

tex = cell(numel(p.gratingPhases), numel(p.gratingSPF));
tex_stand= cell(numel(p.testPhases), numel(p.gratingSPF));

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

 
%Make rects for placing image
imSize = size(grating);
imRect = CenterRectOnPoint([0 0 imSize], cx+imPos(1), cy+imPos(2));

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
%rad = 70; 
rad=round(ang2pix(p.eyerad,p.screenWidthCm, screenWidthPx, p.viewDistCm,'central')); % radius of allowable eye movement in pixels
%rad=50;
%% RUN expt
%HideCursor(window,-1);


  switch p.task % task and demo
    case 1 % DEMO

        %% Sound
        % Initialize the sound driver
        InitializePsychSound(1); % 1 for precise timing
        % Open audio device for low-latency output
        reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
        pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel
        
        %% Generate trials in different conditions %%%%
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
       
        trials2 = fullfact([ numel(p.precueValiditiesWaffle)... % 1 precue validity
             numel(p.plaidOrientations),... % 2 plaid ori
             numel(p.testPhases),... % 3 plaid phase 
             numel(p.plaidContrasts1)]); % 4 plaid contrast
        %% Merge trial count
        nTrials = size(trials1,1)+size(trials2,1);  % total trials = number of grating trials + number of waffle trials
        %nBlocks=nTrials/p.BlockTrials; 
        nBlocks=nTrials/p.BlockTrials; 
        

        trialOrder = randperm(nTrials); %randomize trial order

        instructions = 'This is a demo of the main experiment\n\n';
        %% Show instruction screen and wait for a button press
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a tone that signals the start of the trial followed by either 2 gratings or one plaid.\n\n For the gratings, you will have to report the tilt of the second grating relative to the first.\n\n Your goal is to determine whether the grating is: \n\n counterclockwise (press 9)  \n\n or clockwise (press 0)! \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 ='\n\n For the plaid trials, you must report which of the two orientations in the plaid is stronger in contrast. \n\n For -45 degrees (press 9).  \n\n For +45 degrees (press 0)! \n\n  Press to continue!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        timeInstruct2=Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);
        KbWait(devNum);
        

        instructionsGrating='This is an example of a grating.\n\n Press to continue!';
        DrawFormattedText(window, instructionsGrating, 'center', cy-imRect(1)/2, [1 1 1]*white);
        Screen('DrawTexture', window, tex{2,1}, [], imRect, 45);
        timeGrating = Screen('Flip', window, timeInstruct2 - slack);
        WaitSecs(1);
        KbWait(devNum); 

        instructionsPlaid='This is an example of a plaid patch.\n\n Press to continue!';
        DrawFormattedText(window, instructionsPlaid, 'center', cy-imRect(1)/2, [1 1 1]*white);
        Screen('DrawTexture', window, tex_plaid{3,2}, [], imRect, 45);
        timePlaid = Screen('Flip', window, timeGrating - slack);
        WaitSecs(1);
        KbWait(devNum);
       
        instructionsStartDemo ='After each block, you will be able to see your mean tilt score from 1 to 20. \n\n Your goal is to keep this score as low as possible. \n\n Press to start demo trials!';
        DrawFormattedText(window, instructionsStartDemo, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timePlaid-slack); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        WaitSecs(1);
        KbWait(devNum);

        timeStart = GetSecs;
        correct = [];
        block=1;
        firstNonWaffle=0; % tracking the first trial that is a non-waffle trial
        completedTrials=0;
        for iTrial = 1:nTrials % the iteration in the trial loop
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
        
                % separate staircasing for expected and unexpected
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial
                    if precueValidity==1 %if expected first non-waffle, stair exp is hightest
                        stairIdxExp=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
                        lastFewAccExp=[]; % no previous grating trials have occured so this is empty
                    elseif precueValidity==2 %if unexpected first non-waffle, stair unexp is highest
                        stairIdxUn=length(p.stairs);
                        lastFewAccUnexp=[];  % no previous grating trials have occured so this is empty
                    end
                    firstNonWaffle=1; % firstNonWaffle becomes 1 because this is officially the first grating trial
                
                elseif firstNonWaffle==1 % if this variable is 1, this trial is not the first grating trial
                    corrects=d.correct; %get all corrects
                    skipRowsCorrect=isnan(corrects); %find the NaN values
                    corrects(skipRowsCorrect)=[];% delete the NaN values associated with the waffle trials
                    
                    if precueValidity==1
                        expTrials=d.precueValidity==1;
                        expTrials(skipRowsCorrect)=[];
                        if (expTrials==0)
                            stairIdxExp=length(p.stairs);
                            lastFewAccExp=[];
                        else
                            stairIdxExp_all=d.stairIdxExp;
                            skipUnexp=isnan(stairIdxExp_all);
                            stairIdxExp_all(skipUnexp)=[];
                            stairIdxExp_last=stairIdxExp_all(length(stairIdxExp_all));
                            whichExp=find(expTrials);
                            corrects_exp= (corrects(whichExp));
                            [stairIdxExp lastFewAccExp]=updateStaircase(p.stairs, stairIdxExp_last, corrects_exp, correct); % get new stair index value
                        end
                    

                    elseif precueValidity==2
                        unexpTrials=d.precueValidity==2;
                        unexpTrials(skipRowsCorrect)=[];
                        if (unexpTrials==0) 
                            stairIdxUn=length(p.stairs);
                            lastFewAccUnexp=[]; 
                        else
                            stairIdxUnexp_all=d.stairIdxUn;
                            skipExp=isnan(stairIdxUnexp_all);
                            stairIdxUnexp_all(skipExp)=[];
                            stairIdxUn_last=stairIdxUnexp_all(length(stairIdxUnexp_all));
                            whichUnex=find(unexpTrials);
                            corrects_unexp= (corrects(whichUnex));  
                            [stairIdxUn lastFewAccUnexp]=updateStaircase(p.stairs, stairIdxUn_last, corrects_unexp, correct); % get new stair index value
                    
                        end
                    end
                end

                if staticGrating==1 % standard is +/-45 and test is slightly different
                    if precueValidity==1
                        stairIdx=stairIdxExp;
                    elseif precueValidity==2
                        stairIdx=stairIdxUn;
                    end
                    sOrientation=p.gratingOrientations(Orientation);
                    gOrientation=sOrientation+p.stairs(stairIdx)*differenceMultiplier;
                elseif staticGrating==2 %test is +/-45 and standard is slightly different
                    if precueValidity==1
                        stairIdx=stairIdxExp;
                    elseif precueValidity==2
                        stairIdx=stairIdxUn;
                    end
                    gOrientation=p.gratingOrientations(Orientation);
                    sOrientation=gOrientation+p.stairs(stairIdx)*differenceMultiplier;
                end

            elseif trialIdx>size(trials1,1) 
                plaidStatus=2; % if the trial id is a value > the number of grating trials, this trial will be a waffle trial
                precueValidity = p.precueValiditiesWaffle(trials2(trialIdx-size(trials1,1), precueValidityIdx2));
                plaidOrientation = trials2(trialIdx-size(trials1,1), plaidOrientationIdx);
                plaidPhase = trials2(trialIdx-size(trials1,1), plaidPhaseIdx);
                plaidContrast = trials2(trialIdx-size(trials1,1), plaidContrastIdx);
                pOrientation=p.plaidOrientations(plaidOrientation); %get the orientation value (+/-45) using the plaid orientation index
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial hasn't happened yet
                    stairIdx=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
                    stairIdxExp=length(p.stairs);
                    stairIdxUn=length(p.stairs);
                end
            end
            % tone related to overall orientation
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


            %% Store trial information
             if plaidStatus==1 % if this is a grating trial
                d.stairIdx(iTrial) = stairIdx; %store stair index stimuli
                if precueValidity==1
                    d.stairIdxExp(iTrial) = stairIdxExp; %store stair index stimuli expected
                    d.stairIdxUn(iTrial)=NaN;
                elseif precueValidity==2
                    d.stairIdxExp(iTrial) = NaN; %store stair index stimuli expected
                    d.stairIdxUn(iTrial)=stairIdxUn;
                end

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
                d.stairIdx(iTrial) = NaN;
                d.stairIdxExp(iTrial) = NaN; %store stair index stimuli expected
                d.stairIdxUn(iTrial)=NaN;
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
                timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want resT

                %% Present predictive tone
                PsychPortAudio('FillBuffer', pahandle, tone);
                timeTone = PsychPortAudio('Start', pahandle, [], timeFix+p.signalStart, 1); % waitForStart = 1 in order to return a timestamp of playback
           
                %% Present STANDARD
                Screen('DrawTexture', window, tex_stand{standardPhase,standardSPF}, [], imRect, sOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
                
                %% Fixation
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
                
                %% Present TEST
                Screen('DrawTexture', window, tex{gratingPhase,gratingSPF}, [], imRect, gOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                
                timeT = Screen('Flip', window, timeBlank1 + p.standSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
    
                %% Fixation
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank2 = Screen('Flip', window, timeT + p.imDur - slack); %timeS + how long i want stimulus to be presented for
             
            elseif plaidStatus==2
                % Present fixation rest (grey)
                drawFixation(window, cx, cy, fixSize, p.fixColor*p.dimFactor*white);
                timePreStart = Screen('Flip', window);
                % Present fixation active (white)
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
                
                %% Present predictive tone
                PsychPortAudio('FillBuffer', pahandle, tone);
                timeTone = PsychPortAudio('Start', pahandle, [], timeFix + p.signalStart, 1); % waitForStart = 1 in order to return a timestamp of playback
                
                %% Present stimuli
                Screen('DrawTexture', window, tex_plaid{plaidContrast,plaidPhase}, [], imRect, pOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
    
                % blank
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
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
                elseif (differenceMultiplier==1 && staticGrating==1) ||(differenceMultiplier==-1 && staticGrating==2)  %stimuli CW
                    correct = 0; % stimuli != response
                end
            elseif plaidStatus==2
                correct_p=NaN;
            else
                responseText = 'Please push one of the buttons: 9 (CW), 0(CW)';
            end

            %% Response Acceptance Marker based on Correctness of Response
            if targetResponseKey && plaidStatus==1 && correct==1
                drawFixation(window, cx, cy, fixSize,[0 1 0]*p.fixColor*white);    
                Screen('Flip',window);
     
            elseif targetResponseKey && plaidStatus==1 && correct==0
                drawFixation(window, cx, cy, fixSize,[1 0 0]*p.fixColor*white);
                Screen('Flip',window);
            elseif targetResponseKey && plaidStatus==2
                drawFixation(window, cx, cy, fixSize,[0 0 1]*p.fixColor*white);
                Screen('Flip',window);
            end

            timeEnd=GetSecs();

            %% Keep track of last 3 trials
            if plaidStatus==1  
                if precueValidity==1
                    lastFewAccExp=[lastFewAccExp correct];
                elseif precueValidity==2
                    lastFewAccUnexp=[lastFewAccUnexp correct];
                end
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
            save(sprintf('%s/%s_s%d_demo_predv2_s%s.mat',data.dataDir_sub,p.subjectID,p.sessionNum,date),'d','p');

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
                tiltID=d.stairIdx(blockStartTrial:iTrial);
                nan_tiltID=isnan(tiltID);
                non_nan_tiltID=tiltID(~nan_tiltID);
                blockTilt=mean(non_nan_tiltID);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);
                fprintf('Mean tilt in block: %d.\n-----\n', blockTilt);


                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                tiltMessage= sprintf('Mean Tilt in Block: %d%', round(blockTilt));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s\n\n%s', blockMessage, accMessage, tiltMessage, keyMessage);
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
        end
    case 2
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
        nTrials = size(trials1,1)+size(trials2,1);  % total trials = number of grating trials + number of waffle trials
        %nBlocks=nTrials/p.BlockTrials; 
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
        instructions = 'This is the main version of the experiment\n\n';
        Screen('FillRect', window, white*p.backgroundColor);
        instructions1 = sprintf('%s\n\nThere will be a tone that signals the start of the trial followed by either 2 gratings or one plaid.\n\n For the gratings, you will have to report the tilt of the second grating relative to the first.\n\n Your goal is to determine whether the grating is: \n\n counterclockwise (press 9)  \n\n or clockwise (press 0)! \n\n Press to continue!', instructions);
        DrawFormattedText(window, instructions1, 'center', 'center', [1 1 1]*white);
        timeInstruct1=Screen('Flip', window,p.demoInstructDur-slack);

        KbWait(devNum);

        instructions2 ='\n\n For the plaid trials, you must report which of the two orientations in the plaid is stronger in contrast. \n\n For -45 degrees (press 9).  \n\n For +45 degrees (press 0)! \n\n  Press to continue!';
        DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white); DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window,timeInstruct1+p.demoInstructDur-slack);
        WaitSecs(1);
        KbWait(devNum);
        
        timeStart = GetSecs;
        correct = [];
        block=1;
        eyeSkip = zeros(size(trials1,1)+size(trials2,1),1); % trials skipped due to an eye movement, same size as trials matrix
        skippedTrials = [];
        iTrial=1;
        completedTrials=0;

        disp('reached')
        firstNonWaffle=0; % tracking the first trial that is a non-waffle trial
        
        
        
       
        while iTrial<=nTrials
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix
            stopThisTrial=0;
            % %% Initialize for eye tracking trial breaks
            %  if iTrial>1
            %      eyeSkip(iTrial-1) = stopThisTrial; % this is for the previous trial
            %  end
            %  stopThisTrial = 0;
            % 
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
        
                % separate staircasing for expected and unexpected
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial
                    if precueValidity==1 %if expected first non-waffle, stair exp is hightest
                        stairIdxExp=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
                        lastFewAccExp=[]; % no previous grating trials have occured so this is empty
                    elseif precueValidity==2 %if unexpected first non-waffle, stair unexp is highest
                        stairIdxUn=length(p.stairs);
                        lastFewAccUnexp=[];  % no previous grating trials have occured so this is empty
                    end
                    firstNonWaffle=1; % firstNonWaffle becomes 1 because this is officially the first grating trial
                
                elseif firstNonWaffle==1 % if this variable is 1, this trial is not the first grating trial
                    corrects=d.correct; %get all corrects
                    skipRowsCorrect=isnan(corrects); %find the NaN values
                    corrects(skipRowsCorrect)=[];% delete the NaN values associated with the waffle trials
                    
                    if precueValidity==1
                        expTrials=d.precueValidity==1;
                        expTrials(skipRowsCorrect)=[];
                        if (expTrials==0)
                            stairIdxExp=length(p.stairs);
                            lastFewAccExp=[];
                        else
                            stairIdxExp_all=d.stairIdxExp;
                            skipUnexp=isnan(stairIdxExp_all);
                            stairIdxExp_all(skipUnexp)=[];
                            stairIdxExp_last=stairIdxExp_all(length(stairIdxExp_all));
                            whichExp=find(expTrials);
                            corrects_exp= (corrects(whichExp));
                            [stairIdxExp lastFewAccExp]=updateStaircase(p.stairs, stairIdxExp_last, corrects_exp, correct); % get new stair index value
                        end
                    

                    elseif precueValidity==2
                        unexpTrials=d.precueValidity==2;
                        unexpTrials(skipRowsCorrect)=[];
                        if (unexpTrials==0) 
                            stairIdxUn=length(p.stairs);
                            lastFewAccUnexp=[]; 
                        else
                            stairIdxUnexp_all=d.stairIdxUn;
                            skipExp=isnan(stairIdxUnexp_all);
                            stairIdxUnexp_all(skipExp)=[];
                            stairIdxUn_last=stairIdxUnexp_all(length(stairIdxUnexp_all));
                            whichUnex=find(unexpTrials);
                            corrects_unexp= (corrects(whichUnex));  
                            [stairIdxUn lastFewAccUnexp]=updateStaircase(p.stairs, stairIdxUn_last, corrects_unexp, correct); % get new stair index value
                    
                        end
                    end
                end

                if staticGrating==1 % standard is +/-45 and test is slightly different
                    if precueValidity==1
                        stairIdx=stairIdxExp;
                    elseif precueValidity==2
                        stairIdx=stairIdxUn;
                    end
                    sOrientation=p.gratingOrientations(Orientation);
                    gOrientation=sOrientation+p.stairs(stairIdx)*differenceMultiplier;
                elseif staticGrating==2 %test is +/-45 and standard is slightly different
                    if precueValidity==1
                        stairIdx=stairIdxExp;
                    elseif precueValidity==2
                        stairIdx=stairIdxUn;
                    end
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
                if firstNonWaffle==0 % if this variable is still 0 it means this is the first trial that is a grating trial hasn't happened yet
                    stairIdx=length(p.stairs); % the stair index should be the easiest (i.e. the last) value
                    stairIdxExp=length(p.stairs);
                    stairIdxUn=length(p.stairs);
                end
            end
            
            %% tone related to overall orientation
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

            %% Store trial information
             if plaidStatus==1 % if this is a grating trial
                d.stairIdx(iTrial) = stairIdx; %store stair index stimuli
                if precueValidity==1
                    d.stairIdxExp(iTrial) = stairIdxExp; %store stair index stimuli expected
                    d.stairIdxUn(iTrial)=NaN;
                elseif precueValidity==2
                    d.stairIdxExp(iTrial) = NaN; %store stair index stimuli expected
                    d.stairIdxUn(iTrial)=stairIdxUn;
                end

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
                d.stairIdx(iTrial) = NaN;
                d.stairIdxExp(iTrial) = NaN; %store stair index stimuli expected
                d.stairIdxUn(iTrial)=NaN;
                d.gratingOrientation(iTrial) = NaN; %orientation of test stimuli
                d.gratingPhase(iTrial) = NaN; %phase of test stimuli
                d.gratingContrast(iTrial) = NaN; %contrast of test stimuli
                d.standardOrientation(iTrial) = NaN; %orientation of test stimuli
                d.standardPhase(iTrial) = NaN; %phase of test stimuli
                d.standardContrast(iTrial) = NaN; %contrast of test stimuli
                d.standardSPF(iTrial) = NaN; %contrast of test stimuli
            end
            %% Store stimulus information in trials matrix
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
                    driftCorrected = rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad,fixRect});
                    if driftCorrected
                        % restart trial
                        drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                        timeFix = Screen('Flip', window, timePreStart+p.signalRestDur-slack); %timePreStart + how long i want rest
                    end
                end

                %% Present predictive tone
                PsychPortAudio('FillBuffer', pahandle, tone);
                timeTone = PsychPortAudio('Start', pahandle, [], timeFix+p.signalStart, 1); % waitForStart = 1 in order to return a timestamp of playback
                
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
               
                %% Present STANDARD
                Screen('DrawTexture', window, tex_stand{standardPhase,standardSPF}, [], imRect, sOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
                
                if p.eyeTracking
                    Eyelink('Message', 'EVENT_CUE');
                end
    
                if p.eyeTracking
                    while GetSecs < timeS + p.imDur - p.eyeSlack && ~stopThisTrial
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

                %% Fixation
                drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack); %timeS + how long i want stimulus to be presented for
                
                if p.eyeTracking
                    Eyelink('Message', 'EVENT_CUE');
                end
    
                if p.eyeTracking
                    while GetSecs < timeBlank1 + p.standSOA - p.eyeSlack && ~stopThisTrial
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

                
                %% Present TEST
                Screen('DrawTexture', window, tex{gratingPhase,gratingSPF}, [], imRect, gOrientation);
                drawFixation(window, cx, cy, fixSize, p.fixColor*white)
                
                timeT = Screen('Flip', window, timeBlank1 + p.standSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
    
                 if p.eyeTracking
                    Eyelink('Message', 'EVENT_CUE');
                end
    
                if p.eyeTracking
                    while GetSecs < timeT + p.imDur- p.eyeSlack && ~stopThisTrial
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

                %% Fixation
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
                
                timeS = Screen('Flip', window, timeTone+p.toneSOA - slack); %timeFix+ how much i want to wait from white(active) to standard
                
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
                responseText = 'Please push one of the buttons: 9 (CW), 0 (CW)';
            end

            timeEnd=GetSecs();

            %% Keep track of last 3 trials
            if plaidStatus==1  
                if precueValidity==1
                    lastFewAccExp=[lastFewAccExp correct];
                elseif precueValidity==2
                    lastFewAccUnexp=[lastFewAccUnexp correct];
                end
            end

            %% Store trial info
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
            d.stopThisTrial(iTrial) = stopThisTrial;

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
                tiltID=d.stairIdx(blockStartTrial:iTrial);
                nan_tiltID=isnan(tiltID);
                non_nan_tiltID=tiltID(~nan_tiltID);
                blockTilt=mean(non_nan_tiltID);

                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f (block)', ...
                    100*blockCorrectness);
                fprintf('Mean tilt in block: %d.\n-----\n', blockTilt);


                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                tiltMessage= sprintf('Mean Tilt in Block: %d%', round(blockTilt));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, nBlocks);
                if iTrial==nTrials
                    keyMessage = 'All done! Thank you for participating!';
                else
                    keyMessage = 'Press any key to go on.';
                end

                breakMessage = sprintf('%s\n%s\n\n%s\n\n%s', blockMessage, accMessage, tiltMessage, keyMessage);
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

