function pred_fake
% ___________________________________________________________________
% pred.m runs experiment

% p = parameters
% s = stimuli, rects, trials
% window = window
% data = data
% ___________________________________________________________________
clear; close all;
%% Input
% Subject and session info
p.subjectID = input('Enter subject ID:  ','s');
p.sessionNum = input('Enter session number (1)');
%p.debug = input('Debug? (Y/N)','s');
p.task = input(['Prediction task run:\n' ...
    '1 - Demo\n' ...
    '2 - Task\n' ...
    '3 - Debug']);
p.reps = input ('How many reps(1/2)?');
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

% data.subDir = sprintf('%s/%s',data.dataDir,p.subjectID);
% if ~exist(data.subDir, 'dir')
%     mkdir(data.subDir)
% end

% are we overwriting data
if exist(sprintf('%s/%s_s%d_pred.mat',data.dataDir,p.subjectID,p.sessionNum),'file')
    error('This subject/session already has saved data');
end

% Running on PTB-3? Abort otherwise
AssertOpenGL;

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
[window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);
%[window, rect] = Screen('OpenWindow', screenNumber);

[screenWidthPx, screenHeight] = Screen('WindowSize', window);
flipInterval = Screen('GetFlipInterval', window); % % Get refresh rate, frame duration (s)

% We will request the screen to flip half a refresh (the "slack" time) before we
% actually want the screen to change. This helps to avoid late screen
% flips. So let's define this "slack" variable for convenience.
slack = flipInterval/2;

% Get x and y coordinates for the center of the window
[cx, cy] = RectCenter(rect);
white = WhiteIndex(window);
black = BlackIndex(window);

%% Calculate stimulus dimensions (px) and position
pixelsPerDegree = ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);
fixSize = p.fixSize*pixelsPerDegree;

%% Sound
% Initialize the sound driver
InitializePsychSound(1); % 1 for precise timing
% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel

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

%% Keyboard
% Check all "devices" (keyboards, mice) for response input
devNum = -1;

KbName('UnifyKeyNames');
validKeys = KbName({'1!','2@','9(','0)'});

%% Make TEST stimuli

[s,tex] = pred_makeTest(window, p, s, devNum); 

%% Make STANDARD stimuli
gratingS = rd_grating(pixelsPerDegree, p.gratingSize, p.gratingSF, 0, 0, p.standardContrast); % First make a grating image
imPosS = round(p.imPos*pixelsPerDegree); % standard image position (pixels)
gratingRadiusS = round(p.gratingDiameter/2*pixelsPerDegree); % standard radius (pixels)

% Place an Gaussian aperture on the image to turn it into a Gabor
gaborS = rd_aperture(gratingS, 'gaussian', gratingRadiusS(1)/4);


% Make an image "texture"
% Choose which image you want to make into a texture
imS = gaborS;

% Make the texture
imTexStandard = Screen('MakeTexture', window, imS*white); % multiply by "white" to scale from 0-255
imSizeS = size(imS); % Make the rects for placing the images in the window
imRectS = CenterRectOnPoint([0 0 imSizeS(1) imSizeS(2)], cx+imPosS(1), cy+imPosS(2));

%% SANNITY CHECK: SHOW ALL TEXTURES ONE BY ONE
% % Loop through each texture and display it
% for iC = 1:numel(p.gratingContrasts)
%     for iO = 1:numel(p.gratingOrientations)
%         % Draw the texture on the screen
%         Screen('DrawTexture', window, tex{iC, iO}, [], [], [], [], [], [], [], kPsychDontDoRotation, [0, 0, 1, 1]);
%
%         % Flip the screen to display the texture
%         Screen('Flip', window);
%
%         % Wait for a certain duration (adjust as needed)
%         WaitSecs(1); % Adjust the duration as needed
%     end
% end
%% RUN expt
%HideCursor(window,-1);
%
%
switch p.task % one scenario for now, add demo for teaching purposes! add feedback in demo ????
    case 1
        
    case 2
        %% 1 TASK
        % Make trials structure
        %% %%%% Generate trials in different conditions %%%%

        trials_headers = {'precueValidity','tOrientation','tPhase', 'tContrast','tStatus','precue','responseKey','response','accuracy','rt'};
        % make sure column indices match trials headers
        precueValidityIdx = strcmp(trials_headers,'precueValidity');
        testOrientationIdx = strcmp(trials_headers,'tOrientation'); % 0 (low) or 1 (high)
        testPhaseIdx = strcmp(trials_headers,'tPhase'); % test phase value
        testContrastIdx = strcmp(trials_headers,'tContrast'); % test contrast
        testStatusIdx = strcmp(trials_headers,'tStatus'); % test contrast
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
        if p.debug=="Y" || p.debug=="y"
            nTrials=p.debugTrials;
        elseif p.debug=="N" || p.debug=="n" && p.reps==1
            nTrials = size(trials,1); %448
            nBlocks=p.trialsNeeded/nTrials; %1400/448 = 3.25
            %nTrials = p.nTrialsPerContrast*size(trials,1);
            %nBlocks=1400/
        elseif
        end

        trialOrder = randperm(nTrials);

        instructions = 'This is the main experiment\n\n';
        %% Show instruction screen and wait for a button press
        Screen('FillRect', window, white*p.backgroundColor);
        instructionsFull = sprintf('%s\n\nThere will be a reference patch followed by a predictive tone and a second patch. \n\nYour goal is to determine whether the second patch is: \n\n higher contrast and counterclockwise (press 1), \n\n lower contrast and counterclockwise (press 2), \n\lower contrast and clockwise (press 9), \n\nor higher contrast and clockwise (press 0) \n\nwith reference to the first patch!', instructions);
        DrawFormattedText(window, instructionsFull, 'center', 'center', [1 1 1]*white);
        Screen('Flip', window); %command to change what's on the screen. we're drawing on a behind (hidden) window, and at the moment we screenflip, we flip that window to the front position'
        %pred_instructions(window, p, devNum);
        KbWait(devNum);
        timeStart = GetSecs;
        correct = [];
        block=1;

        for iTrial = 1:nTrials % the iteration in the trial loop
            trialIdx = trialOrder(iTrial); % the trial number in the trials matrix

            %%SANTIY CHECK
            fprintf('Trial %d/%d\n', iTrial, nTrials);

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

            %% Set up stimuli for this trial
            % tone
            % if testOrientation==1
            %     tone=cueTones(1,:);
            % elseif testOrientation==2
            %     tone=cueTones(2,:);
            % end

            %tone = cueTones(toneVersion,:);


            %% Store stimulus information in trials matrix
            
            trials(trialIdx, precueIdx) = toneVersion;


            %% %%%% Play the trial %%%%
            %% Present fixationwhite
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timePreStart = Screen('Flip', window);
            drawFixation(window, cx, cy, fixSize, p.fixColor*white*p.dimFactor);
            timeFix = Screen('Flip', window, p.fixSOA+timePreStart-slack);

            %% Present STANDARD image
            drawFixation(window, cx, cy, fixSize, p.fixColor*white) % DO I EVEN NEED THIS????
            Screen('DrawTexture', window, imTexStandard, [], imRectS, 0);
            timeS = Screen('Flip', window, timeFix+p.standSOA - slack);

            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            timeBlank1 = Screen('Flip', window, timeS + p.imDur - slack);

            %% Present predictive tone
            PsychPortAudio('FillBuffer', pahandle, tone);
            timeTone = PsychPortAudio('Start', pahandle, [], timeBlank1 + p.toneSOA, 1); % waitForStart = 1 in order to return a timestamp of playback

            %% Present TEST image
            %drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            if testStatus==1
                Screen('DrawTexture', window, tex{testContrast, testPhase}, [], imRect, orientation);
                %timeT = Screen('Flip', window, timeTone + p.toneSOA - slack); % is it p.toneSOA, DIDN'T GET LOGIC

                %Screen('DrawTexture', windowPointer, texturePointer [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);
            elseif testStatus==0
                %drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                drawFixation(window, cx, cy, fixSize, white*p.backgroundColor);
            end
            timeT = Screen('Flip', window, timeTone + p.toneSOA - slack); % is it p.toneSOA, DIDN'T GET LOGIC
            % blank
            drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            Screen('Flip', window, timeT + p.imDur - slack);
            %timeBlank2 = Screen('Flip', window, timeT + p.imDur - slack);

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
                %drawFixation(window, cx, cy, fixSize,[0 1 0]*p.fixColor*white);
                %Screen('Flip', window, timeTargetResponse+p.ITI - slack);
                Screen('Flip',window);
                % blank
                %drawFixation(window, cx, cy, fixSize, p.fixColor*white);
                %Screen('Flip', window, timeEnd - slack);
            end

            

            %% MEASURE ISSUES
            %  SHOULD I DO SEPARATE MEASURES FOR CONTRAST AND
            % ORIENTATION

            %CONTRAST AND ORIENTATION: THIS WOULDN'T MAKE SENSE BECAUSE
            %WHAT IF THEY WERE THE SAME, YOU WOULD NEVER BE ABLE TO REPORT
            %YOUR PERCEPTION OF THAT AND SO YOU CAN'T MEASURE ACCURACY OF
            %AN INCOMPLETE THING
            % if testContrast>p.standardContrast && testOrientation==1 % higher 45 CW
            %     if strcmp('0)',targetResponseKeyName)% higher 45 CW
            %         correct = 1;
            %     elseif strcmp('1!',targetResponseKeyName) || strcmp('9)',targetResponseKeyName)|| strcmp('2@',targetResponseKeyName)
            %         correct = 0;
            %     end
            % elseif testContrast>p.standardContrast && testOrientation==2 %higher -45 CCW
            %     if strcmp('1!',targetResponseKeyName)% higher -45 CCW
            %         correct = 1;
            %     elseif strcmp('0)',targetResponseKeyName) || strcmp('9)',targetResponseKeyName)|| strcmp('2@',targetResponseKeyName)
            %         correct = 0;
            %     end
            % elseif testContrast<p.standardContrast && testOrientation==2 %higher  CCW
            %
            %     targetPresent && strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName) % -45, 1 or 2
            %     if testOrientation==-45
            %         correct = 1;
            %     elseif targetOri==45
            %         correct = 0;
            %     end
            % else
            %     responseText = 'Please push valid key [1,2] or [9,0]';
            % end

            %% Percent correct for orientation
     
            if testStatus~=0 && (strcmp('9(',targetResponseKeyName) || strcmp('0)',targetResponseKeyName)) % response CCW
                if testOrientation==2 % stimuli CW
                    correct = 1;  % stimuli = response
                elseif testOrientation==1 % stimuli CCW
                    correct = 0; % stimuli != response
                end
            elseif testStatus~=0 && (strcmp('1!',targetResponseKeyName) || strcmp('2@',targetResponseKeyName)) % response CCW
            %elseif testStatus~=0 && (targetResponseKey==1 || targetResponseKey==2) % response CCW
                if testOrientation==1 %stimuli CCW
                    correct = 1; % stimuli = response
                elseif testOrientation==2 %stimuli CW
                    correct = 0; % stimuli != response
                end
            else
                responseText = 'Please push one of the buttons: 1 (higher CCW), 2 (lower CCW), 9 (lower CW), 0(higher CW)';
            end

            %% Response Acceptance Marker based on Correctness of Response
            % if targetResponseKey && correct==1
            %     drawFixation(window, cx, cy, fixSize,[0 1 0]*p.fixColor*white);
            %     %drawFixation(window, cx, cy, fixSize,*p.fixColor*white);
            %     %Screen('Flip', window, timeTargetResponse+p.ITI - slack);
            %     Screen('Flip',window);
            %     % blank
            %     %drawFixation(window, cx, cy, fixSize, p.fixColor*white);
            %     %Screen('Flip', window, timeEnd - slack);
            % 
            % elseif targetResponseKey && correct==0
            %     drawFixation(window, cx, cy, fixSize,[1 0 0]*p.fixColor*white);
            %     Screen('Flip',window);
            % end
            % 
            timeEnd=GetSecs();

            %% Store trial info
            trials(trialIdx, rtIdx) = targetRT;
            trials(trialIdx, responseKeyIdx) = targetResponseKey;
            p.trials=trials;

            tic
            % save data
            d.block(iTrial) = block;
            d.timeStart(iTrial)=timeStart;
            d.timeEnd(iTrial)=timeEnd;
            d.timeElapsed(iTrial)=timeEnd-timeStart;
            d.testStatus(iTrial) = testStatus;
            d.testOrientation(iTrial) = testOrientation;
            d.testPhase(iTrial) = testPhase;
            d.testContrast(iTrial) = testContrast;
            d.precueValidity(iTrial) = precueValidity;
            d.correct(iTrial) = correct;
            d.timeTargetResponse(iTrial)=timeTargetResponse;
            d.targetRT(iTrial) = targetRT;
            d.targetResponseKey(iTrial) = targetResponseKey;
            %d.trials(iTrial)=trials(trialIdx, : );
            save(sprintf('%s/%s_s%d_pred.mat',data.dataDir,p.subjectID,p.sessionNum),'d','p');

            %time between trials after response
            timeSpent=toc;
            WaitSecs(p.ITI-timeSpent);

            if mod(iTrial,p.BlockTrials)==0 || iTrial>p.BlockTrials*block
                % Calculate block accuracy
                blockStartTrial = (iTrial/p.BlockTrials)*p.BlockTrials - p.BlockTrials + 1;
                if blockStartTrial < 0 % we are doing less than one block
                    blockStartTrial = 1;
                end
               
                non_nan_correctnessIdx= ~isnan(d.correct(blockStartTrial:iTrial));
                non_nan_correctness=d.correct(non_nan_correctnessIdx);
                blockCorrectness = mean(non_nan_correctness);
                
                fprintf('Block %d of %d complete.\n-----\n', block, nBlocks);
                fprintf('Percent Correct: %.2f / %.2f (block/total)', ...
                    100*blockCorrectness, 100*mean(d.correct));
           
                accMessage = sprintf('Percent Correct: %d%', round(blockCorrectness*100));
                blockMessage = sprintf('%d of %d blocks completed. Great Job! Keep Going!', block, ceil(nTrials/p.BlockTrials));
                if iTrial==nTrials
                    keyMessage = '';
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
            end
        end

    otherwise
        error('Expt stage not found.')
end

PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')
ShowCursor;
end

