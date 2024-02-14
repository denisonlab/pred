function [s,tex] = pred_stim(window, p, s, devNum, data)

grating = rd_grating(p.ppd, p.gratingSize, p.gratingSF, 0, 0, 1); %makes grating

gabor = rd_aperture(grating, 'gaussian', p.gaborSD*p.ppd); %makes aperature
im = gabor;
imSize = size(im); %size of gabor

s.imPos = round(p.imPos*p.ppd); % image position (pixels) 
s.gratingRadius = round(p.gratingDiameter/2*p.ppd); % radius (pixels) 

tex = cell(numel(p.gratingContrasts), numel(p.gratingOrientations));

for iC = 1:numel(p.gratingContrasts)
    contrast = p.gratingContrasts(iC);
    for iO=1:numel(p.gratingOrientations)
        orientation = p.gratingOrientations(iO);
        grating = rd_grating(p.ppd, p.imSize, ...
            p.gratingSF, orientation, p.phase, contrast); % 0 to 1 
        [gabor, aps] = rd_aperture(grating, 'gaussian', s.gratingRadius(1)/4);
        tex{iC,iO} = Screen('MakeTexture', window, gabor*s.white);
    end
end

s.imSize = size(grating);
s.imRect = CenterRectOnPoint([0 0 s.imSize], s.cx+s.imPos(1), s.cy+s.imPos(2));

