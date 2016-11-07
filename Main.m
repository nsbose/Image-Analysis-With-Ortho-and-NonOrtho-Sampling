%Smart Camera: Scene Adaptive Image Capture
%EC 520. Spring 2015
%Neladri Bose and Shivani Sheopory 
%Professor: Janusz Konrad

%Read images from folder. From the Test Images set, add all the images you
%want to test to the folder containing these programs. The imagefiles
%variable will cread all these images and the program will run 5 figures
%for each image being tested.

imagefiles = dir('*.jpg'); 	 
nfiles = length(imagefiles);	% Number of files found

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   img_load = imread(currentfilename);

%grayImageCompare.m will output the PSNR and MSE value for orthogonal and
%non-orthogonal cases. To use the function, enter the name of the image
%being tested without the file type
[OrthPSNR NonOrthPSNR OrthMSE NonOrthMSE] = grayImageCompare(img_load, currentfilename); 


end
