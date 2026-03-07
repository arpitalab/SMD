function results = ls_int_gaussian(I, sigma, ridge, max_iter, damp, convergence, divergence)
%LS_INT_GAUSSIAN  Maximum-likelihood fit of symmetric 2D integrated Gaussian PSF
%
%   Uses radial symmetry for initial guess, then Levenberg-Marquardt refinement.
%
%   Args:
%       I           : 2D array, the observed spot image
%       sigma       : float, fixed Gaussian width (default 1.0)
%       ridge       : float, initial Hessian regularization (default 1e-4)
%       max_iter    : int, max iterations (default 20)
%       damp        : float, LM damping factor (default 0.3)
%       convergence : float, convergence threshold for y/x updates (default 1e-4)
%       divergence  : float, divergence threshold (default 1.0, currently unused)
%
%   Returns:
%       results     : struct with fields y, x, I0, bg, errors, H_det, rmse, etc.

    if nargin < 2, sigma       = 1.0;      end
    if nargin < 3, ridge       = 0.0001;   end
    if nargin < 4, max_iter    = 20;       end
    if nargin < 5, damp        = 0.3;      end
    if nargin < 6, convergence = 1e-4;     end
    if nargin < 7, divergence  = 1.0;      end

    % ---- Initial guess -------------------------------------------------
    [y0,x0,sig]=radialcenter(I);                     % you must provide rs.m or implement it
    bg_guess  = ring_mean(I);

    I0_guess = estimate_I0(I, y0, x0, bg_guess, sigma);

    guess = [y0, x0, I0_guess, bg_guess];   % [y, x, I0, bg]

    % ---- Core fitting --------------------------------------------------
    [pars, err, H_det, rmse, n_iter] = fit_ls_int_gaussian( ...
        I, guess, sigma, ridge, max_iter, damp, convergence, divergence);

    % ---- Quality control and SNR ---------------------------------------
    % Replace check_2d_gauss_fit with your own sanity check if desired
    error_flag = 0;  % placeholder

    snr = estimate_snr(I, pars(3), sigma);   % now passes sigma for amp_from_I

    % ---- Output struct -------------------------------------------------
    results = struct( ...
        'y',        pars(1), ...
        'x',        pars(2), ...
        'I0',       pars(3), ...
        'bg',       pars(4), ...
        'y_err',    err(1), ...
        'x_err',    err(2), ...
        'I0_err',   err(3), ...
        'bg_err',   err(4), ...
        'H_det',    H_det, ...
        'error_flag', error_flag, ...
        'snr',      snr, ...
        'rmse',     rmse, ...
        'n_iter',   n_iter);

% =========================================================================
% =========================== LOCAL FUNCTIONS =============================
% =========================================================================

    function m = ring_mean(I)
        if numel(I) < 4
            m = mean(I(:));
            return;
        end
        top    = mean(I(1, 1:end-1));
        bottom = mean(I(end, 2:end));
        left   = mean(I(2:end-1, 1));
        right  = mean(I(1:end-1, end));
        m = (top + bottom + left + right) / 4;
    end

    function v = ring_var(I, ddof)
        if nargin < 2, ddof = 0; end
        
        [Ny, Nx] = size(I);
        
        % Handle very small images
        if Ny <= 2 || Nx <= 2
            v = var(I(:), ddof);
            return;
        end
        
        % Extract all unique border pixels and concatenate as column vectors
        top    = I(1,   :)';          % Nx × 1  (full top row including corners)
        right  = I(1:end-1, end);     % (Ny-1) × 1  (right column excluding top corner)
        bottom = I(end, end-1:-1:1)'; % Nx × 1  (bottom row reversed, excluding right corner)
        left   = I(end-1:-1:2, 1);    % (Ny-1) × 1  (left column reversed, excluding bottom/top)
        
        % Vertically concatenate all into one long column vector
        border = [top; right; bottom; left];
        
        v = var(border, ddof);  % ddof = 0 (population), 1 (sample)
    end
    function I0_est = estimate_I0(I, y, x, bg, sigma)
        if nargin < 5, sigma = 1.0; end
        [~, lin_idx] = max(I(:));
        [ym, xm] = ind2sub(size(I), lin_idx);
        psf_y = int_gauss_psf_1d(ym-1, y-1, sigma);  % convert to 0-based temporarily if needed
        psf_x = int_gauss_psf_1d(xm-1, x-1, sigma);
        denom = psf_y * psf_x;
        if abs(denom) < eps
            I0_est = NaN;
        else
            I0_est = (I(ym,xm) - bg) / denom;
        end
    end

    function snr = estimate_snr(I, I0, sigma)
        if nargin < 3, sigma = 1.0; end
        amplitude_sq = amp_from_I(I0, sigma)^2;
        bg_var = ring_var(I, 1);  % sample variance
        if bg_var <= 0
            snr = Inf;
        else
            snr = amplitude_sq / bg_var;
        end
    end

    function amp = amp_from_I(I0, sigma)
        if nargin < 2, sigma = 1.0; end
        amp = I0 / (2 * pi * sigma^2);
    end

    function proj = int_gauss_psf_1d(Y, yc, sigma)
        if nargin < 3, sigma = 1.0; end
        S = sigma * sqrt(2);
        arg1 = (Y + 0.5 - yc) / S;
        arg2 = (Y - 0.5 - yc) / S;
        proj = 0.5 * (erf(arg1) - erf(arg2));
    end

    function dproj = int_gauss_psf_deriv_1d(Y, yc, sigma)
        if nargin < 3, sigma = 1.0; end
        A = 2 * sigma^2;
        B = sigma * sqrt(2 * pi);
        term1 = exp( -(Y - 0.5 - yc).^2 / A );
        term2 = exp( -(Y + 0.5 - yc).^2 / A );
        dproj = (term1 - term2) / B;
    end

    function [pars, err, H_det, rmse, n_iter] = fit_ls_int_gaussian(img, guess, sigma, ridge, max_iter, damp, convergence, ~)
        if nargin < 3, sigma = 1.0; end
        if nargin < 4, ridge = 0.0001; end
        if nargin < 5, max_iter = 20; end
        if nargin < 6, damp = 0.3; end
        if nargin < 7, convergence = 1e-4; end

        [Ny, Nx] = size(img);
        n_pixels = Ny * Nx;

        index_y = (0:Ny-1)';
        index_x = (0:Nx-1);

        pars   = guess(:)';
        update = ones(1,4);

        iter_idx = 0;
        while iter_idx < max_iter && any(abs(update(1:2)) > convergence)
            proj_y = int_gauss_psf_1d(index_y, pars(1), sigma);
            proj_x = int_gauss_psf_1d(index_x, pars(2), sigma);

            dproj_y_dy = int_gauss_psf_deriv_1d(index_y, pars(1), sigma);
            dproj_x_dx = int_gauss_psf_deriv_1d(index_x, pars(2), sigma);

            J_col1 = pars(3) * (dproj_y_dy * proj_x);
            J_col2 = pars(3) * (proj_y * dproj_x_dx);
            J_col3 = (proj_y * proj_x);
            J_col4 = ones(Ny, Nx);

            J = [J_col1(:), J_col2(:), J_col3(:), J_col4(:)];

            model = pars(3) * J_col3(:) + pars(4);

            residuals = img(:) - model;
            sig2 = mean(residuals.^2);

            grad = (J' * residuals) / sig2;

            H = -(J' * J) / sig2;
            H_inv = invert_hessian_local(H, ridge);

            update = - (H_inv * grad);   % grad is 4x1, H_inv is 4x4
            pars = pars + damp * update';

            iter_idx = iter_idx + 1;
        end

        err = sqrt(diag(-H_inv))';
        H_det = det(H_inv);
        rmse = sqrt(sig2);
        n_iter = iter_idx;
    end

    function H_inv = invert_hessian_local(H, ridge)
        H_reg = H;
        H_reg(1:5:end) = H_reg(1:5:end) + ridge;
        H_inv = inv(H_reg);
    end

end