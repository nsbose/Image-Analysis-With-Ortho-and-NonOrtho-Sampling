function [ OrthPSNR NonOrthPSNR OrthMSE NonOrthMSE ] = grayImageCompare( img_input, str_title )
%GRAYIMAGECOMPARE 
%   Checks the shape of the image spectrum
%   Analyzes gray scale images in both orthogonal and non-orthogonal manner
%   Outputs their respective MSE and PSNR values

img_gray = rgb2gray(im2double(img_input));
isoCheck(img_gray,str_title);
[NonOrthPSNR NonOrthMSE NonOrthImage] = NonOrthoCheck(img_gray,str_title);
[OrthPSNR OrthMSE OrthImage] = OrthoCheck(img_gray,str_title);

str_Orth = '_Orth.tif';
str_NonOrth = '_NonOrth.tif';
str_Orth_C = [str_title str_Orth];
str_NonOrth_C = [str_title str_NonOrth];
str_gray = '_gray.tif';
str_gray_C = [str_title str_gray];

imwrite(OrthImage, str_Orth_C); %Saves the orthogonal image
imwrite(NonOrthImage, str_NonOrth_C); %Saves the non-orthogonal image
imwrite(img_gray,str_gray_C); %Saves the grayscale image

end

