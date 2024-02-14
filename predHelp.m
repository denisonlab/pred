% rd_PTBDemo.m

% This script shows how to do some common things we do in the lab with 
% Psychtoolbox.
%
% Rachel Denison
% October 2021

%% Handy commands
% "clear all"   Clear all variables from the workspace
% "close all"   Close all figures
% "sca"         (or "ScreenCloseAll") Close the PTB window
% Ctrl+c        Kill any processes that are running (if the program is hanging)
% Cmd+0         Put the cursor in the command window
% Tab           Press tab when typing to show options for completing a
%               function name or variable (aka "tab complete")
% Cmd+/         Comment lines of code
% Cmd+t         Uncomment lines of code
% Cmd+i         Magic indent
% Shift+Fn+F7   Execute lines of code
% "help [function name]"  Show the documentation for that function

%% PTB setup
% Check we are running PTB-3
AssertOpenGL;

% Skip screen tests - ONLY for demo, not for real experiments
Screen('Preference', 'SkipSyncTests', 1);

%% Basic info
% Name the subject
subjectID = 'test';

% Load parameters for this experiment
p = predParams;

%% %%%% Setup: window, sound, and keyboard %%%%
%% Window
% Here we open a PTB window and see how to get several properties of the
% window.

%% Setting up the window ...
% Pick the screen on which to display the window
screenNumber = max(Screen('Screens'));

% Open a PTB window
% [window, rect] = Screen('OpenWindow', screenNumber); % defaults to full screen
[window, rect] = Screen('OpenWindow', screenNumber, [], [0 0 600 400]);

%% Getting many useful properties of the window ...
% Get x and y coordinates for the center of the window
[cx, cy] = RectCenter(rect);

% Get the color code for white
white = WhiteIndex(window);  

% Get window size
[screenWidthPx, screenHeightPx] = Screen('WindowSize', window);

% Get refresh rate
flipInterval = Screen('GetFlipInterval', window); % frame duration (s)

% We will request the screen to flip half a refresh (the "slack" time) before we 
% actually want the screen to change. This helps to avoid late screen
% flips. So let's define this "slack" variable for convenience.
slack = flipInterval/2;

%% Sound
% Initialize the sound driver
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, p.Fs, 1); % 1 = single-channel

%% Keyboard
% Check all "devices" (keyboards, mice) for response input
devNum = -1;

%% %%%% Make stimuli %%%%
%% Making images ...
%% Make a black square with a white stripe
% We make images scaled 0-1. 0 = black, 1 = white.
sz = 100;
square = zeros(sz,sz); 
square(30:50,:) = 1; % set a horizontal stripe to 1

% Two ways to view the image
figure
imagesc(square)
colorbar

figure
imshow(square)

%% Exercises with the square image
% 1. Turn the horizontal stripe gray
% 2. Make a black square image with a white vertical stripe
% 3. Make a black square image witha smaller white square inside

%% Calculate stimulus dimensions (px) and position
pixelsPerDegree = ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);

fixSize = p.fixSize*pixelsPerDegree;

%% Make a gabor image
% First make a grating image
grating = rd_grating(pixelsPerDegree, p.gratingSize, p.gratingSF, 0, 0, 1);

% View the grating
figure
imshow(grating)

% Place an Gaussian aperture on the image to turn it into a Gabor
gabor = rd_aperture(grating, 'gaussian', p.gaborSD*pixelsPerDegree);

% View the gabor
figure
imshow(gabor)

%% Exercises with the gabor
% 1. Change the spatial frequency of the gabor
% 2. Change the orientation of the gabor
% 3. Change the contrast of the gabor

%% Make an image "texture"
% Choose which image you want to make into a texture
im = gabor;

% Make the texture
imTex = Screen('MakeTexture', window, im*white); % multiply by "white" to scale from 0-255

%% Make the rects for placing the images in the window
imSize = size(im);
imRect = CenterRectOnPoint([0 0 imSize(1) imSize(2)], cx+p.imPos(1), cy+p.imPos(2));

%% Making sounds ...
% 10^0.5 for every 10dB
%% Make a pure tone
tone0 = MakeBeep(p.toneFreq, p.toneDur, p.Fs);

% View tone
figure
t = 0:1/p.Fs:p.toneDur;
plot(t, tone0)
xlabel('Time (s)')
ylabel('Amplitude')

% Listen to tone
sound(tone0, p.Fs)

%% Apply an envelope so the sound doesn't click at the beginning and end
tone = applyEnvelope(tone0, p.Fs);

% View tone
figure
plot(t, tone)
xlabel('Time (s)')
ylabel('Amplitude')

% Listen to tone
sound(tone, p.Fs)

%% Exercises with tones
% 1. Make the tone 500 ms long
% 2. Change the tone frequency

%% %%%% Present one trial %%%%
%% Show instruction screen and wait for a button press
Screen('FillRect', window, white*p.backgroundColor);
DrawFormattedText(window, 'Press a key as soon as you see the image\n\nPress any key to begin', 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
KbWait(devNum);
timeStart = GetSecs;

%% Present fixation
drawFixation(window, cx, cy, fixSize, p.fixColor*white);
timeFix = Screen('Flip', window);

%% Present tone
PsychPortAudio('FillBuffer', pahandle, tone);
timeTone = PsychPortAudio('Start', pahandle, [], [], 1); % waitForStart = 1 in order to return a timestamp of playback

%% Present image
drawFixation(window, cx, cy, fixSize, p.fixColor*white);
Screen('DrawTexture', window, imTex, [], imRect);
timeIm = Screen('Flip', window);

% blank
drawFixation(window, cx, cy, fixSize, p.fixColor*white);
timeBlank = Screen('Flip', window, timeIm + p.imDur - slack);

%% Collect response 
[secs, keyCode] = KbWait(devNum);
rt = secs - timeIm;
responseKey = find(keyCode);

%% Exercises with stimulus presentation and trial sequence
% 1. Change the location where the image is presented to somewhere in the
% upper right quadrant of the screen
% 2. Set the image to be presented exactly 1 s after the tone
% 3. Switch the order of the tone and the image. Now the image should be 1
% s after the tone.
% 4. Update the RT so that it tells us how fast the observer responded to
% the tone. (You can change the task instructions as well.)

%% Store trial info
trialsHeaders = {'RT','ResponseKey'};
trialIdx = 1;
trials(trialIdx, 1) = rt;
trials(trialIdx, 2) = responseKey;

DrawFormattedText(window, sprintf('Your reaction time was %.2f s!', rt), 'center', 'center', [1 1 1]*white);
Screen('Flip', window);
WaitSecs(2);

%% Store expt info
expt.subjectID = subjectID;
expt.p = p;
expt.trialsHeaders = trialsHeaders;
expt.trials = trials;

%% Clean up
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
Screen('CloseAll')