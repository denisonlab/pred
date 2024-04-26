function p = predParams(p)

%% Screen
p.screenWidthCm = 20;
p.viewDistCm = 57;
p.eyerad=2.5;
%% Fixation
p.fixSize = 0.1; % degrees visual angle  0.15 used before making smaller
p.fixColor = 1; % white
p.dimFactor=0.5;

%% Debug
p.debugTrials=30*2;
p.BlockTrialsDebug=2;

%% Demo
p.demoInstructDur=3;

%% Gratings
p.backgroundColor = 0.5; % gray
p.imPos = [0 0];
p.gratingSize = 2; % degrees visual angle

p.gratingSPF=[1.0 1.5];
p.gratingOrientations = [-45 45]; % vertical 
p.gaborSD = 0.5; % about 4 SDs will be visible at full contrast
p.aperture = 'gaussian';
p.gratingDiameter = [2 0];
p.testStatus = [1]; %0 (absent stimuli)

p.imSize = 4;
p.testPhases = linspace(0,2*pi,4);
p.gratingPhases=p.testPhases+pi;

%% Plaid parameters
% p.plaidContrasts1 = [logspace(log10(0.21),log10(0.6),7)];
% p.plaidContrasts2 = 0.2;
p.plaidContrastsRange=[logspace(log10(.1333),log10(.3),5)];
p.plaidContrasts1=[p.plaidContrastsRange(1) p.plaidContrastsRange(2) p.plaidContrastsRange(3) p.plaidContrastsRange(3) p.plaidContrastsRange(3) p.plaidContrastsRange(3) p.plaidContrastsRange(4) p.plaidContrastsRange(5)];
p.plaidContrasts2 = 0.2;
p.plaidAxes = [0 90]; % base orientations 
p.plaidEcc = 1; % dva from center for L and R stimultaneous stimuli 

p.plaidOrientations= [-45 45];
p.plaidStatus=[1 1 1 2]; 
p.plaidSPF=1.5;
p.precueValiditiesWaffle=[1 2];

%% Waffles+Kok
p.standardContrast2 = 0.8; 
p.gratingContrast2 = 0.8;
p.gratingOrientationDiff=[-5 5];
p.staticGrating=[1 2];
%% Calculate pixels per degree
p.ppd = ang2pix(1,  p.screenWidthCm, 1920, p.viewDistCm); %

%% Sounds
p.Fs = 44100; % samples per second
p.toneFreqs = [450 1000]; % Hz
p.toneDur = 0.2; % s tone duration 200 ms
p.deviceName = 'Scarlett';

%% Timing
p.fixSOA=0.5; %precue fixation time 
p.standSOA = 0.1; % standard to tone
p.toneSOA = 0.55; % tone to standard 550
%p.toneSOA = 2; % tone to standard 550
p.imDur = 0.5; % image duration 0.03
%p.imDur = 2; % image duration 0.03
p.signalStart=1;
p.signalRestDur=1;
p.maskDur=0.02;
p.maskSOA = 0; % should be 0
p.ITI=0.2;
%% Stairs
p.stairs=[logspace(log10(0.1),log10(10),20)];
p.differenceMultiplier=[-1 1];
%% Validity
p.precueNames = {'valid','invalid'}; %p.precueNames = {'valid','neutral','invalid'};
p.precueValidities = [1 1 1 2]; %tone validity 75% 1 in 4 is invalid

%% Trials & Blocks
p.BlockTrials=48;
p.repScale1=2; % 42 mins for 896 trials
p.repScale2=3; % 63 mins for 1344 trials
p.repScaleDemo=1;

%% Eyetracking
p.eyeSlack=0.12;
