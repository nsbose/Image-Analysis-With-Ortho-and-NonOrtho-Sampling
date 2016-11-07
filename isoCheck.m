function isoCheck( img_input,img_title )
%ISOCHECK 
%   Converts grayscale image into frequency domain and retains 1/4
%   highest value coefficients

t_img_gray = img_input;
figure; 
str_OrigTitle= sprintf('Original Image: ');
subplot(2,2,1); imshow(img_input); title({str_OrigTitle,img_title});

ff_img_gray = fft2(t_img_gray);
subplot(2,2,3); imshow(fftshift(log(abs(ff_img_gray)+1)),[]);
title('DFT Magnitude of Original Image');

ff_uimg = ff_img_gray;
[uRows uCols] = size(ff_uimg);
vect_ff_uimg = ff_uimg(:); vec_ff_uE = vect_ff_uimg.^2; %finds the energy
[Bu,Iu] = sort(vec_ff_uE,'descend'); t_thres = round((1/4)*length(Bu));
%retains 1/4 amount of coeffs
for i = 1:length(Bu)
    if (i > t_thres)
        vect_ff_uimg(Iu(i)) = 0;
    end
end

newMatu = reshape(vect_ff_uimg,[uRows,uCols]);
subplot(2,2,2); imshow(ifft2(newMatu)); title('Image with coeff removed');
subplot(2,2,4); imshow(fftshift(log(abs(newMatu)+1)),[]);
title('DFT Magnitude of Quantized Image');

end

