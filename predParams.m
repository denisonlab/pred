function p = predParams(p)

%% Screen
p.screenWidthCm = 20;
p.viewDistCm = 57;
%% Fixation
p.fixSize = 0.15; % degrees visual angle
p.fixColor = 1; % white
p.dimFactor=0.5;
%% Debug
p.debugTrials=30*2;
p.BlockTrialsDebug=2;

%% Demo
p.demoInstructDur=3;

%% Images
p.backgroundColor = 0.5; % gray
p.imPos = [0 0];
p.gratingSize = 2; % degrees visual angle
p.gratingSF = 3; % cycles per degree90020012200110001190292
p.gratingOrientations = [-45 45]; % vertical 
p.phase = 0;
p.gaborSD = 0.5; % about 4 SDs will be visible at full contrast
p.aperture = 'gaussian';
p.gratingDiameter = [2 0];
% 0.125 0.25 0.5 1
%p.gratingContrasts = [0 0.14678 0.215443 0.316228 0.464159 0.681292 1];
% 2.5 to 50 percent, 6-79% 
%p.gratingContrasts = [0 logspace(log10(0.025),log10(0.5),6)]; % 7 contrast levels and 0 (absent stimuli) 
p.gratingContrasts = [0 logspace(log10(0.2),log10(0.8),6)];
p.maskContrast=0.8;
%p.gratingContrasts=[1];
p.testStatus=[1]; %0 (absent stimuli)
p.standardContrast = p.gratingContrasts(4); %WHAT SHOULD CONTRAST BE???
p.standardOrientation=0;
%p.imSize=[2 2]; %p.imSize=2;
p.imSize=2;
p.testPhases = linspace(0,2*pi,4); %%PHASE WHAT AM I GONNA DO WITH IT

%% Calculate pixels per degree
%pix= ang2pix(1, p.screenWidthCm, screenWidthPx, p.viewDistCm);
%[p.screenWidthPx, p.screenHeightPx] = Screen('WindowSize', window);
p.ppd = ang2pix(1,  p.screenWidthCm, 1920, p.viewDistCm); %

%% Sounds
p.Fs = 44100; % samples per second
p.toneFreqs = [440 1760]; % Hz % try 990 Hz instead ?
p.toneDur = 0.2; % s TONE DURATION 200 ms
p.deviceName = 'Scarlett';
%% Clicks
% p.clickDur=0.0005;
% p.clickRampDur=p.clickDur/2;
% p.clickFreq=8000;
%% Timing
p.fixSOA=0.5; %precue fixation time 
p.standSOA = 1; % standard to tone
p.toneSOA = 0.5; % tone to test
p.imDur = 0.03; % image duration 0.03
p.ITI=0.2;
p.signalStart=1;
p.signalRestDur=1;
p.maskDur=0.02;
p.maskSOA = 0; % should be 0
%% Validity
p.precueNames = {'valid','invalid'}; %p.precueNames = {'valid','neutral','invalid'};
p.precueValidities = [1 1 1 2]; %tone validity 75% 1 in 4 is invalid

%% Trials & Blocks
p.trialsNeeded=7*50*2*2; %7 contrast levels x 50 trials per contrast level x 2 absent/present x 2 validity
p.BlockTrials=64; %64 trials 3 mins, 128 trials approx 6 mins
%p.BlockTrials=5;%eyetracker testing
p.repScale1=2; % 42 mins for 896 trials
p.repScale2=3; % 63 mins for 1344 trials
p.repScaleDemo=1;
%% Eyetracking
p.eyeSlack=0.12;
