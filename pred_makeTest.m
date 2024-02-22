function [s,tex] = pred_makeTest(window, p, s, devNum, data)

s.imPos = round(p.imPos*pixelsPerDegree); % image position (pixels)
s.gratingRadius = round(p.gratingDiameter/2*pixelsPerDegree); % radius (pixels)

tex = cell(numel(p.gratingContrasts), numel(p.gratingOrientations));

for iC = 1:numel(p.gratingContrasts)  %HOW TO LOOP THROUGH PHASE OR VARY ACROSS PHASES
    contrast = p.gratingContrasts(iC);
    for iP=1:numel(p.testPhases)
        phase = p.testPhases(iP);
        grating = rd_grating(pixelsPerDegree, p.imSize, ...
            p.gratingSF, 0, phase, contrast); % 0 to 1
        [gabor, aps] = rd_aperture(grating, 'gaussian', gratingRadius(1)/4);
        tex{iC,iP} = Screen('MakeTexture', window, gabor*white);
    end
end

%Make rects for placing image
s.imSize = size(grating);
s.imRect = CenterRectOnPoint([0 0 imSize], cx+imPos(1), cy+imPos(2));