function varargout = llr(I, k, w, t, return_filt)
%LLR  Log-likelihood ratio test for Gaussian spots

    if nargin < 2, k = 1.0; end
    if nargin < 3, w = 9; end
    if nargin < 4, t = 20.0; end
    if nargin < 5, return_filt = false; end

    [H, W] = size(I);
    I = double(I);  % Ensure double precision for accuracy

    [G_rft, Sgc2] = mle_amp_setup(H, W, k, w);

    n_pixels = w^2;

    % Local mean and mean of squares using uniform filter
    uniform_kernel = ones(w, w) / n_pixels;
    A = imfilter(I, uniform_kernel, 'conv', 'same');
    B = imfilter(I.^2, uniform_kernel, 'conv', 'same');

    % FFT convolution with centered Gaussian
    I_rft = fft2(I);
    full_conv = ifft2(I_rft .* G_rft);  % Sizes now match perfectly
    C = real(fftshift(full_conv));
    C(C < 0) = 0;  % Only convex spots

    % Log-likelihood ratio
    variance = B - A.^2;
    denom = n_pixels * Sgc2 * variance;

    ratio = zeros(size(C));
    valid = denom > eps;  % Avoid division by zero or tiny values
    ratio(valid) = C(valid).^2 ./ denom(valid);
    ratio(~valid) = 0.001;

    L = 1 - ratio;

    % Zero out near edges
    hw = floor(w/2);
    L(1:hw, :) = 1;
    L(:, 1:hw) = 1;
    L(end-hw+1:end, :) = 1;
    L(:, end-hw+1:end) = 1;

    score_img = -(n_pixels / 2) * log(max(L, eps));  % Avoid log(0)

    % Detect spots
    binary = score_img > t;
    [y, x] = find(binary);
    coords = [y, x];  % [y x] = row/col

    if return_filt
        varargout{1} = score_img;
        varargout{2} = binary;
        varargout{3} = coords;
    else
        varargout{1} = coords;
    end

end


