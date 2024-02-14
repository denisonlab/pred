%% PTB setup
% Check we are running PTB-3
AssertOpenGL;

% Skip screen tests - ONLY for demo, not for real experiments
Screen('Preference', 'SkipSyncTests', 1);

%Priority(MaxPriority(w)); %max priority for ptb
%% Basic info
% Name the subject
subjectID = 'test';

% Load parameters for this experiment
p = predParams;
%% Screens
screenNumber = max(Screen('Screens')); %screen to display on
% Open a PTB window
% [window, rect] = Screen('OpenWindow', screenNumber); % defaults to full screen
[window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);

%% Getting many useful properties of the window ...
% Get x and y coordinates for the center of the window
[cx, cy] = RectCenter(rect);
% Get the color code for white
white = WhiteIndex(window);  

% Get window size
[p.screenWidthPx, p.screenHeightPx] = Screen('WindowSize', window);

% Get refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)

% We will request the screen to flip half a refresh (the "slack" time) before we 
% actually want the screen to change. This helps to avoid late screen
% flips. So let's define this "slack" variable for convenience.
slack = flipInterval/2;

% 
% %% Keyboard
% % Check all "devices" (keyboards, mice) for response input
% devNum = -1;
% Make vertical gratings at different contrasts
grating = rd_grating(p.ppd, p.gratingSize, ...
    p.gratingSF, p.gratingOrientations, p.phase, 1); % 0 to 1 contrastThreshold = 1
[gabor, ~] = rd_aperture(grating, 'gaussian', s.gratingRadius(1)/4);
% Make texture
targetTex = Screen('MakeTexture', window, gabor*s.white);

% 
% %% Exercises with the square image
% % 1. Turn the horizontal stripe gray
% % 2. Make a black square image with a white vertical stripe
% % 3. Make a black square image witha smaller white square inside
% 
% pixelsPerDegree = ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);
% fixSize = p.fixSize*pixelsPerDegree;
% %% Make a gabor image
% % First make a grating image
% grating = rd_grating(pixelsPerDegree, p.gratingSize, p.gratingSF, 0, 0, 1);
% 
% % Place an Gaussian aperture on the image to turn it into a Gabor
% gabor = rd_aperture(grating, 'gaussian', p.gaborSD*pixelsPerDegree);
% im = gabor;
% 
% % Make the texture
% imTex = Screen('MakeTexture', window, im*white); % multiply by "white" to scale from 0-255
% 
% %% Make the rects for placing the images in the window
% imSize = size(im);
% imRect = CenterRectOnPoint([0 0 imSize(1) imSize(2)], cx+p.imPos(1), cy+p.imPos(2));

%% Making sounds ...
% 10^0.5 for every 10dB
%% Sound
% Initialize the sound driver
InitializePsychSound(1); % 1 for precise timing
% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel
% %% Make a pure tone
tone0 = MakeBeep(p.toneFreq, p.toneDur, p.Fs);




%% Apply an envelope so the sound doesn't click at the beginning and end
tone = applyEnvelope(tone0, p.Fs);

PsychPortAudio('FillBuffer', pahandle, tone);
timeTone = PsychPortAudio('Start', pahandle, [], [], 1); % waitForStart = 1 in order to return a timestamp of playback

PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll') % CLOSE SCREEN