function [im] = readimg(imgfile)

im = im2double(flipud(fitsread(imgfile)));
im(im<0) = 0;

end

