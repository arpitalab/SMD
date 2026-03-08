function varargout = llr(I, k, w, t, return_filt)
%LLR  Log-likelihood ratio test for Gaussian spots
%
%   coords = llr(I, k, w, t)
%   [score_img, binary, coords] = llr(I, k, w, t, true)
%
%   Tests the likelihood of a Gaussian spot centered in each w x w
%   sub-window versus flat background with Gaussian noise.
%
%   Reference implementation: alecheckert/quot (findSpots.py)
%
%   Args:
%       I           : 2D image (will be cast to double)
%       k           : Gaussian kernel sigma (default 1.0)
%       w           : window size, odd integer (default 9)
%       t           : score threshold (default 20.0)
%       return_filt : if true, return [score_img, binary, coords]
%
%   Returns:
%       coords      : Nx2 [row, col] coordinates of detected spots
%                     (one per connected component, at max-score position)

    if nargin < 2, k = 1.0; end
    if nargin < 3, w = 9; end
    if nargin < 4, t = 20.0; end
    if nargin < 5, return_filt = false; end

    [H, W] = size(I);
    I = double(I);

    % Kernel setup (centered Gaussian in frequency domain)
    [G_rft, Sgc2] = mle_amp_setup(H, W, k, w);
    n_pixels = w^2;

    % --- Local statistics via uniform filtering ---------------------------
    % Use 'symmetric' boundary to match scipy.ndimage.uniform_filter('reflect')
    uniform_kernel = ones(w, w) / n_pixels;
    A = imfilter(I, uniform_kernel, 'conv', 'same', 'symmetric');
    B = imfilter(I.^2, uniform_kernel, 'conv', 'same', 'symmetric');

    % --- Centered-Gaussian convolution via FFT ----------------------------
    C = real(fftshift(ifft2(fft2(I) .* G_rft)));
    C(C < 0) = 0;  % only convex (bright) spots

    % --- Log-likelihood ratio ---------------------------------------------
    % L = 1 - C^2 / (n * Sgc2 * local_variance)
    % Match quot's stable_divide_array: clip denominator to >= 0.001
    local_var = B - A.^2;
    denom = n_pixels * Sgc2 * local_var;
    denom = max(denom, 0.001);     % stable division (quot: zero=0.001)

    L = 1.0 - (C.^2) ./ denom;

    % Zero out edges (same as quot)
    hw = floor(w / 2);
    L(1:hw, :)       = 1.0;
    L(:, 1:hw)       = 1.0;
    L(end-hw+1:end, :) = 1.0;
    L(:, end-hw+1:end) = 1.0;

    % Score image: -(n/2)*log(L)
    % Clamp L to small positive value to avoid log(<=0) and preserve
    % dynamic range among strong spots (eps is too small — gives a flat
    % plateau that breaks peak finding)
    score_img = -(n_pixels / 2) * log(max(L, 1e-6));

    % --- Spot extraction via connected-component labeling -----------------
    % Match quot's threshold_image + label_spots(mode='max'):
    % one detection per connected component, at the max-score position
    binary = score_img > t;
    [labels, n_labels] = bwlabel(binary, 8);
    coords = zeros(n_labels, 2);
    for ii = 1:n_labels
        mask = (labels == ii);
        masked_score = score_img;
        masked_score(~mask) = -Inf;
        [~, lin_idx] = max(masked_score(:));
        [r, c] = ind2sub([H, W], lin_idx);
        coords(ii, :) = [r, c];
    end

    % --- Output -----------------------------------------------------------
    if return_filt
        varargout{1} = score_img;
        varargout{2} = binary;
        varargout{3} = coords;
    else
        varargout{1} = coords;
    end
end
