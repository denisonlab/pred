function tintedimage = imtint(im, rgbtints)

% function tintedimage = imtint(im, rgbtints)
%
% <im> is the input image (grayscale or RGB), <rgbtints> is a [R G B] 1x3 
% array of scalar tint values, which will be multiplied by image values in 
% each color channel.

if size(im,3) == 1 % if 2d grayscale image
    
    tintedimage = cat(3, im*rgbtints(1), im*rgbtints(2), im*rgbtints(3));
    
elseif size(im,3) == 3 % if 3d RGB image

    tintedimage(:,:,1) = im(:,:,1) * rgbtints(1);
    tintedimage(:,:,2) = im(:,:,2) * rgbtints(2);
    tintedimage(:,:,3) = im(:,:,3) * rgbtints(3);

else
    
    error('The input image has unexpected dimensions.')
    
end

% color image may need to be uint8 to be displayed properly
tintedimage = uint8(tintedimage);