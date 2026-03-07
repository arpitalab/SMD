function [G_rft, Sgc2] = mle_amp_setup(H, W, k, w)
%MLE_AMP_SETUP  Prepare centered Gaussian kernel in frequency domain (full  H x W shape)

    S = 2 * (k^2);

    % Symmetric coordinate grid centered at zero
    half = (w-1)/2;
    [X, Y] = meshgrid(-half:half, -half:half);
    dist_sq = X.^2 + Y.^2;

    g = exp(-dist_sq / S);
    g = g / sum(g(:));

    gc = g - mean(g(:));  % mean-subtracted (centered)

    Sgc2 = sum(gc(:).^2);

    % Create exactly H x W zero-padded array
    gc_padded = zeros(H, W);

    % Center position (1-based indexing)
    cy = floor(H/2) + 1;
    cx = floor(W/2) + 1;
    hw = floor(w/2);

    % Place the w x w kernel exactly at the center
    y_idx = cy - hw : cy - hw + w - 1;
    x_idx = cx - hw : cx - hw + w - 1;
    gc_padded(y_idx, x_idx) = gc;

    % Full FFT2 — this is H x W complex array
    G_rft = fft2(gc_padded);

end