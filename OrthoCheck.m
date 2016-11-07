function [ OrthPSNR OrthMSE OrthImage ] = OrthoCheck( img_input,img_title )
%% ORTHOCHECK 
%   Subsamples quantized image on orthogonal square storage lattice
%   and returns restored image and MSE and PSNR values

[N M] = size(img_input); %Determines the matrix size

%because 1/sqrt(2) = irrational number, we are using (2/3) (computationally
%faster)
L = round(N*(2/3));
a = round((N/2)-(L/2)+1); b = a+L-1;
%Generates the square filter
s_filter = zeros(N);
for i = a:b
	for j = a:b
    	s_filter(i,j) = 1;
	end
end

%Gaussian to that filter
fN = 15; sigma = 20;
[X Y] = meshgrid(round(-fN/2):round(fN/2), round(-fN/2):round(fN/2));
gfil = exp(-X.^2/(2*sigma^2)-Y.^2/(2*sigma^2));
gfil = gfil./sum(gfil(:));

s_filter = imfilter(s_filter,gfil,'same'); %Adds the Gaussian to the filter
%Filter applied to the image
ff_img = fft2(img_input); sf_img = fftshift(ff_img).*s_filter;
s_img = ifft2(ifftshift(sf_img));

figure;
str_OrigTitle= sprintf('Original Image: ');
subplot(2,2,1); imshow(img_input); title({str_OrigTitle,img_title});
subplot(2,2,2); imshow(fftshift(log(abs(ff_img)+1)),[]); title('Magnitude of DFT of the original image');
subplot(2,2,3); imshow(s_img); title('Image with Square Filter');
subplot(2,2,4); imshow(fftshift(log(abs(fft2(s_img))+1)),[]); title('Magnitude of DFT of the modified Image');

%% Downsampling by 1.5 to display the image
%Upsampling the image orthogonally by 2
uf = [0 ((1/2)-0.1) ((1/2)+0.1) 1]; ua = [2.0 2.0 0 0];
ub = firpm(19,uf,ua); mat_ub = ub.'*ub;
u_img = upsample(s_img,2); u_img = upsample(u_img',2); u_img = u_img';
uf_img = imfilter(u_img,mat_ub);
%Downsampling the image orthogonally by 3
df = [0 ((1/3)-0.01) ((1/3)+0.01) 1]; da = [1.0 1.0 0.0 0.0];
db = firpm(31,df,da); mat_db = db.'*db; d_img = imfilter(uf_img,mat_db);
d_img2 = downsample(uf_img,3); d_img2 = downsample(d_img2',3);
d_img2 = d_img2';
t_hold = d_img2;
%% Upsampling by 1.5 to display the image
% Up sampling by a factor 3
uf = [0 ((1/3)-0.1) ((1/3)+0.1) 1]; ua = [3.0 3.0 0 0]; 
ub = firpm(20,uf,ua); mat_ub = ub.'*ub;
u_img = upsample(d_img2,3); u_img = upsample(u_img',3); u_img = u_img';
uf_img = imfilter(u_img,mat_ub);
% Up sampling by a factor 2
df = [0 ((1/2)-0.05) ((1/2)+0.05) 1]; da = [1.0 1.0 0.0 0.0];
db = firpm(30,df,da); mat_db = db.'*db; d_img = imfilter(uf_img,mat_db);
d_img2 = downsample(uf_img,2); d_img2 = downsample(d_img2',2);
d_img2 = d_img2';

%Adjusts the rows and cols due to rounding
[Rows Cols] = size(img_input);
d_img2 = d_img2([1:Rows],[1:Cols]);
%Applies the square filter
ff_dmg = fftshift(fft2(d_img2)); ff_new = ff_dmg.*s_filter;
d_img2 = ifft2(ifftshift(ff_new));
%% 
% PSNR Calculations
figure;
str_OrigTitle= sprintf('Original Image: ');
subplot(2,2,1); imshow(img_input); title({str_OrigTitle,img_title});
subplot(2,2,2); imshow(t_hold); title('Subsampled Image');
OrthPSNR = psnr(d_img2,img_input);%PSNR Calculation
% MSE Calculation
E = (img_input - d_img2);
OrthMSE = abs((sum(sum(E.^2))) / (M * N));

str_Title = sprintf('Restored Image');
str_psnr = sprintf('PSNR value: %0.2f', OrthPSNR);
str_mse = sprintf('MSE value: %0.2e', OrthMSE);
subplot(2,2,3); imshow(d_img2); title({str_Title,str_psnr,str_mse});
subplot(2,2,4); imshow(fftshift(log(abs(fft2(d_img2))+1)),[]); title('DFT Magnitude of Restored Image'); 

OrthImage = d_img2;
end

