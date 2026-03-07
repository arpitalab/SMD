function [filtered, stdWavelet1] = wavelet_filter(image)

%Applies a wavelet filter to the input image as described in:
%
%Izeddin et al. "Wavelet analysis for single molecule localization
%microscopy", Opt. Express 20, 2081-2095 (2012) 
%
%Input:
%   image       -   Single input image given by 2d array of pixel values
%
%Output:
%   filtered    -   Wavelet filtered image
%
%
%   stdWavelet1 -   Standard deviation of the first wavelet map.
%                   This can interpreted as a measure of the background
%                   noise and is later used to calculate the detection
%                   threshold.



bordObj = 4; % the amount of pixels on the border which need to be set to 0.

% initialize wavelet coefficients (undecimated wavelet transform using cubic
% B-spline function)
H0 = 3/8;
H1 = 1/4;
H2 = 1/16;

% the filter matrices are obtained from wavelet coefficients
filter_1 = [H2,H1,H0, H1, H2]' * [H2,H1,H0, H1, H2];
filter_2 = [H2,0,H1, 0, H0,0,H1, 0, H2]' * [H2,0,H1,0,H0,0 H1,0, H2];

if ~isa(image,'double') && ~isa(image,'single')
    image = single(image);
end

% Filter using the first matrix
Coef1 = conv2(image,filter_1,'same');   

% Filter using the second filter matrix
Coef2 = conv2(Coef1,filter_2,'same');   % the second coefficient map

%Create background noise image by subtracting the first coefficient map
%from the original image
Wavelet1 = image - Coef1;

%Get standard deviaion of the first wavelet filtered plan 
%(estimation of the background noise)
stdWavelet1 = std(Wavelet1(:));

% the second wavelet map is the filtered image and is obtained as the
% difference of coefficient maps.
Wavelet2 = Coef1 - Coef2;


% the borders need to be set to zeros to exclude the resulting weird
% effects near image borders.
Wavelet2(1:bordObj,:) = 0;
Wavelet2((end - bordObj + 1):end,:) = 0;
Wavelet2(:,1:bordObj) = 0;
Wavelet2(:,(end - bordObj + 1):end) = 0;
Wavelet2(Wavelet2 < 0) = 0;
filtered = Wavelet2;

end
