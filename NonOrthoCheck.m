function [ NonOrthPSNR NonOrthMSE NonOrthImage ] = NonOrthoCheck( img_input,img_title )
%% NONORTHOCHECK 
%   Sub-samples quantized image on non-orthogonal diamond storage lattice
%   and returns restored image and MSE and PSNR values

[N M] = size(img_input); %Determines the matrix size
 
%Creation of the diamond filter
%Builds top half of the diamond
n = N/2;
d_filter = zeros(N);
for i = 1:n
	for j = ((n+1)-i):round(n+i)
    	d_filter(i,j) = 1.0;
	end
end
%Builds bottom half of the diamond
c = 0;
for i = n+1:N
	for j = (i-n):N-c
    	d_filter(i,j) = 1.0;
	end
	c = c+1;
end

%Gaussian 

fN = 15; sigma = 20;
[X Y] = meshgrid(round(-fN/2):round(fN/2), round(-fN/2):round(fN/2));
gfil = exp(-X.^2/(2*sigma^2)-Y.^2/(2*sigma^2));
gfil = gfil./sum(gfil(:));

d_filter = imfilter(d_filter,gfil,'same');
ff_img_input = fft2(img_input); df_img = (fftshift(ff_img_input)).*d_filter;
t_img = ifft2(ifftshift(df_img)); %diamond filter applied coeff

figure;
str_OrigTitle= sprintf('Original Image: ');
subplot(2,2,1); imshow(img_input); title({str_OrigTitle,img_title});
subplot(2,2,2); imshow(fftshift(log(abs(ff_img_input)+1)),[]); title('Magnitude of DFT of the original image');
subplot(2,2,3); imshow(t_img); title('Image with Diamond Filter');
subplot(2,2,4); imshow(fftshift(log(abs(fft2(t_img))+1)),[]); title('Magnitude of DFT of the Image with Diamond Filter');

%Manual subsampling
for i = 1:2:N
    for j = 2:2:N
        t_img(i,j) = 0;
    end
end
for h = 2:2:N
    for k = 1:2:N
        t_img(h,k) = 0;
    end
end

figure; 
subplot(2,2,2); imshow(t_img); title('Subsampled Image [2X 0], [X X]');
ff_t_img = fft2(t_img); 
d_filter = d_filter.*2; %Upsample by a factor of 2 because we are down sampling similar
ff_output = fftshift(ff_t_img).*d_filter;
d_img = ifft2(ifftshift(ff_output));

%Display
str_Title= sprintf('Original Image: ');
subplot(2,2,1); imshow(img_input); title({str_Title,img_title});
NonOrthPSNR = psnr(d_img,img_input); %PSNR Calculations

%MSE Calculations
E = (img_input - d_img);
NonOrthMSE = abs((sum(sum(E.^2))) / (M * N));

str_Title = sprintf('Restored Image');
str_psnr = sprintf('PSNR value: %0.2f', NonOrthPSNR);
str_mse = sprintf('MSE value: %0.2e', NonOrthMSE);
subplot(2,2,3); imshow(d_img); title({str_Title,str_psnr,str_mse});

subplot(2,2,4); imshow(fftshift(log(abs(fft2(d_img))+1)),[]); title('DFT Magnitude of Restored Image'); 


NonOrthImage = d_img;
end

