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

%% Keyboard
% Check all "devices" (keyboards, mice) for response input
devNum = -1;

KbName('UnifyKeyNames');

maxNumChar=80;
vLineSpacing=1;

x=cx;
y=cy;

textColor=white;
useKbCheck=1;
bgcolor=black;

Screen('FillRect', window, white*bgColor);
instructions1='Did you notice anything about the experiment? Please type your answer';
[response1,temp1]=GetEchoStringPreserve(window,instructions1, x, y, white, white*bgColor,length(instructions1),vLineSpacing,1);

Screen('FillRect', window, white*bgColor);
instructions2 = 'Did you notice any relationship between the tone and stimuli?';

[response2,temp2]=GetEchoStringPreserve(window, instructions2, x, y, textColor, bgColor,length(instructions2),vLineSpacing,1);

Screen('FillRect', window, white*bgColor);
instructions3 = 'There was a relationship between tone and stimuli orientation! If you said yes, you are correct. Tell us what you think it is or make a guess. ';

[response3,temp3]=GetEchoStringPreserve(window, instructions3, x, y, textColor, bgColor,length(instructions3),vLineSpacing,1);

end