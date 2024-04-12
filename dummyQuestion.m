%dummy
function dummyQuestion

clear; close all;
PsychPortAudio('Close');

% Running on PTB-3? Abort otherwise
AssertOpenGL;


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

% Get x and y coordinates for the center of the window
[cx, cy] = RectCenter(rect);
white = WhiteIndex(window);
black = BlackIndex(window);
bgColor=0.5;
%% Calculate stimulus dimensions (px) and position
pixelsPerDegree = ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);
fixSize = p.fixSize*pixelsPerDegree;
gaborSD = round(p.gaborSD * p.ppd); 

%% Keyboard
% Check all "devices" (keyboards, mice) for response input
devNum = -1;

KbName('UnifyKeyNames');

maxNumChar==200;
vLineSpacing=2;
Screen('FillRect', window, white*bgColor);
instructions = 'Did you notice anything about the experiment?';
DrawFormattedText(window, instructions, 'center', 'center', [1 1 1]*white);

GetEchoStringPreserve(windowPtr, msg, x, y, textColor, bgColor, useKbCheck, varargin,maxNumChar,vLineSpacing);

Screen('FillRect', window, white*bgColor);
instructions2 = 'Did you notice any relationship between tone and stimuli?';
DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);

GetEchoStringPreserve(windowPtr, msg, x, y, textColor, bgColor, useKbCheck, varargin,maxNumChar,vLineSpacing);

Screen('FillRect', window, white*bgColor);
instructions2 = 'There was a relationship between tone and stimuli orientation!\n\n If you said yes, you are correct. Tell us what you think it is. \n\n If you said no, that is okay. Try to guess what the relationship is. ';
DrawFormattedText(window, instructions2, 'center', 'center', [1 1 1]*white);

GetEchoStringPreserve(windowPtr, msg, x, y, textColor, bgColor, useKbCheck, varargin,maxNumChar,vLineSpacing);

end