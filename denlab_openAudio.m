function pahandle = denlab_openAudio(deviceName, reqlatencyclass, Fs)
% pahandle = denlab_openAudio(deviceName)
% Opens PsychPortAudio and Snd to the same pahandle for desired sound device. 
% Eyetracker uses Snd; Psychtoolbox uses PsychPortAudio. 
% Inputs: 
%   deviceName: (string) desired device name for presenting audio stimuli 
%       Options: 'Scarlett' 'sysdefault' 'default' 
%       If no deviceName specified, opens 'Scarlett' if detected, 'default' otherwise 
%   reqlatencyclass: (double) elects how aggressive PsychPortAudio should be about
%       minimizing sound latency and getting good deterministic timing 
%       1 (default) = Try to get the lowest latency that is possible under the constraint of reliable playback
%   Fs: (double) Requested playback/capture rate in samples per second (Hz)
%       44100 (default) 
% Outputs: 
%   pahandle: (double) device handle for the device

%% Close audio, then initialize 
PsychPortAudio('Close');

% Initialize the sound driver
InitializePsychSound(1);

%% Find all sound devices 
PsychPortAudioDevices = PsychPortAudio('GetDevices');
PsychPortAudioDeviceNames = {PsychPortAudioDevices.DeviceName};
PsychPortAudioDeviceIndex = [PsychPortAudioDevices.DeviceIndex];

%% Check inputs 
if nargin < 3
    Fs = 44100; 
end
if nargin < 2
    reqlatencyclass = 1; 
end
if nargin < 1 
    if any(strcmp(PsychPortAudioDeviceNames,'Scarlett'))
        deviceName = 'Scarlett'; 
        disp('Scarlett sound device found...')
    elseif any(strcmp(PsychPortAudioDeviceNames,'default'))
        deviceName = 'default'; 
        disp('Default sound device found...')
    else 
        error('Scarlett and default sound devices not found. Please check sound setup.')
    end
end 

%% Find desired sound device
deviceNameIdxs = find( contains(PsychPortAudioDeviceNames, deviceName) );
deviceIdx =  PsychPortAudioDeviceIndex( deviceNameIdxs(:,1) );
if isempty(deviceNameIdxs)
    error('Desired sound device %s not found. Try again.', deviceName);
else
    fprintf('Desired sound device %s found. Initializing...', deviceName);
end

%% Open PsychPortAudio and Snd to same pahandle on desired sound device 
pahandle = PsychPortAudio('Open', deviceIdx(1), 1, reqlatencyclass, Fs, 2); %scarlett, mode, latency, Fs, stereo
Snd('Open', pahandle, 1); % links eyetracker sound output to psychportaudio
fprintf('%s initialized to pahandle %d ...', deviceName, pahandle)