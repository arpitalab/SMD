function localize(self)
%LOCALIZE  Detect and fit single-molecule spots in an image stack
%
%   self.localize()
%
%   Dispatches to one of two pipelines based on self.detection_method:
%     'wavelet' - a-trous wavelet + Crocker-Grier + Poisson MLE Gaussian
%     'llr'     - log-likelihood ratio + radial symmetry / integrated Gaussian
%
%   Results are stored in self.spots as a cell array (one entry per frame).
%   Each entry is an Nx9 matrix:
%     [x, y, A, sigma, frame, SNR, iframe, negL, variance]
%
%   See also: SMD, track

    % --- Load image stack -------------------------------------------------
    if ~endsWith(self.dirname, '/') && ~endsWith(self.dirname, '\')
        self.dirname = [self.dirname '/'];
    end
    fullpath = [self.dirname, self.fname];
    [~,~,ext] = fileparts(which(fullpath));
    if isempty(ext)
        [~,~,ext] = fileparts(fullpath);
    end

    if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
        stack = tiffread2(fullpath);
    elseif strcmpi(ext, '.nd2')
        stack = mynd2read(fullpath);
    else
        error('SMD:localize', 'Unsupported file format: %s', ext);
    end
    stack = data2photons(stack);

    num_frames = min(1200, length(stack));
    self.num_frames = num_frames;
    fprintf('Localizing %d frames (%s pipeline)\n', num_frames, self.detection_method);

    % --- Dispatch ---------------------------------------------------------
    switch lower(self.detection_method)
        case 'wavelet'
            [spots, stop_frame, timestamps, maxint_img] = ...
                detect_fit_wavelet(self, stack, num_frames);
        case 'llr'
            [spots, stop_frame, timestamps, maxint_img] = ...
                detect_fit_llr(self, stack, num_frames);
        otherwise
            error('SMD:localize', 'Unknown detection_method: %s', self.detection_method);
    end

    % --- Store results ----------------------------------------------------
    self.spots      = spots;
    self.stop_frame = stop_frame;
    self.maxint     = maxint_img;
    if ~isempty(timestamps)
        self.timestamps = timestamps;
    end
end

% =========================================================================
%  WAVELET PIPELINE (existing method, refactored)
% =========================================================================
function [spots, stop_frame, timestamps, maxint_img] = detect_fit_wavelet(self, stack, num_frames)

    thresh = self.localization_threshold;
    dx     = self.localization_box;  % half-width of fitting box

    spots       = cell(1, num_frames);
    timestamps  = [];
    stop_frame  = num_frames;
    total_spots = 0;
    discarded   = 0;

    [~,~,ext] = fileparts(self.fname);

    [H, W] = size(stack(1).data);
    pixel_values = uint16(zeros(H, W, num_frames));

    se = strel('disk', 8);

    h = waitbar(0, 'Wavelet detection...');
    for iframe = 1:num_frames
        waitbar(iframe / num_frames, h);

        if strcmpi(ext, '.nd2')
            timestamps(iframe) = stack(iframe).ts;
        end
        pixel_values(:,:,iframe) = uint16(stack(iframe).data);

        [bg_mean, ~] = calcBG(double(stack(iframe).data), H, W, 21);
        mean_sub = double(stack(iframe).data) - bg_mean;
        mean_sub(mean_sub < 0) = 0;

        % A-trous wavelet filter
        [filtered, stdWavelet1] = wavelet_filter(mean_sub);

        % Denoised image for SNR
        im1 = imtophat(mean_sub, se);
        im2 = wiener2(im1, [3 3]);
        im3 = imgaussfilt(im2, 2);

        % Crocker-Grier peak finding
        out = pkfnd(filtered, thresh * stdWavelet1, 7);

        if isempty(out) || size(out,1) == 0
            spots{iframe} = [1 1 -1 1 1 1 iframe 0 0];
            stop_frame = min(stop_frame, iframe);
            continue;
        end

        loc = cntrd(filtered, out, 7, 0);
        total_spots = total_spots + size(loc, 1);
        loc2 = zeros(size(loc,1), 9);

        for ii = 1:size(loc, 1)
            px = round(loc(ii,2));
            py = round(loc(ii,1));

            % Bounds check
            if px-dx < 1 || px+dx > H || py-dx < 1 || py+dx > W
                loc2(ii,:) = [loc(ii,1) loc(ii,2) 0 0 iframe 0 iframe 0 0];
                continue;
            end

            [A, x0, y0, sigma, ~, negL, variance] = ...
                gaussfit2DMLE(double(filtered(px-dx:px+dx, py-dx:py+dx)));
            loc2(ii,1) = x0 + py - (dx+1);
            loc2(ii,2) = y0 + px - (dx+1);
            loc2(ii,3) = A;
            loc2(ii,4) = sigma;
            loc2(ii,5) = iframe;
            loc2(ii,8) = negL;
            loc2(ii,9) = variance;

            try
                [A1, ~] = calc_snr(loc2(1:ii,:), im3);
            catch
                A1 = 3;
            end
        end
        loc2(:,6) = A1(:);
        loc2(:,7) = iframe;

        % Keep spots with SNR > 2
        idx = loc2(:,6) > 2;
        discarded = discarded + (size(out,1) - sum(idx));
        if any(idx)
            spots{iframe} = loc2(idx, :);
        else
            spots{iframe} = [1 1 -1 1 1 1 iframe 0 0];
            stop_frame = min(stop_frame, iframe);
        end
    end
    close(h);
    fprintf('Wavelet: discarded %d / %d spots\n', discarded, total_spots);
    maxint_img = max(pixel_values, [], 3);
end

% =========================================================================
%  LLR PIPELINE (new)
% =========================================================================
function [spots, stop_frame, timestamps, maxint_img] = detect_fit_llr(self, stack, num_frames)

    k = self.llr_params.k;
    w = self.llr_params.w;
    t = self.llr_params.t;
    dx = self.localization_box;
    psf_sig = self.psf_sigma;
    min_sep = 2 * dx + 1;  % minimum separation between spots (pixels)

    spots       = cell(1, num_frames);
    timestamps  = [];
    stop_frame  = num_frames;
    total_spots = 0;
    discarded   = 0;

    [~,~,ext] = fileparts(self.fname);

    [H, W] = size(stack(1).data);
    pixel_values = uint16(zeros(H, W, num_frames));

    se = strel('disk', 8);

    h = waitbar(0, 'LLR detection...');
    for iframe = 1:num_frames
        waitbar(iframe / num_frames, h);

        if strcmpi(ext, '.nd2')
            timestamps(iframe) = stack(iframe).ts;
        end
        pixel_values(:,:,iframe) = uint16(stack(iframe).data);

        % Background subtraction
        [bg_mean, ~] = calcBG(double(stack(iframe).data), H, W, 21);
        mean_sub = double(stack(iframe).data) - bg_mean;
        mean_sub(mean_sub < 0) = 0;

        % LLR detection: score image + binary + spot coords
        % llr() now returns one detection per connected component
        % at the max-score position (matching quot's label_spots)
        [score_img, ~, det_coords] = llr(mean_sub, k, w, t, true);

        if isempty(det_coords)
            spots{iframe} = [1 1 -1 1 1 1 iframe 0 0];
            stop_frame = min(stop_frame, iframe);
            continue;
        end

        py_all = det_coords(:, 1);  % row
        px_all = det_coords(:, 2);  % col

        npeaks = numel(py_all);
        total_spots = total_spots + npeaks;
        loc2 = zeros(npeaks, 9);

        % Denoised image for SNR (same as wavelet path)
        im1 = imtophat(mean_sub, se);
        im2 = wiener2(im1, [3 3]);
        im3 = imgaussfilt(im2, 2);

        for ii = 1:npeaks
            py = py_all(ii);
            px = px_all(ii);

            % Bounds check
            if py-dx < 1 || py+dx > H || px-dx < 1 || px+dx > W
                loc2(ii,:) = [px py 0 0 iframe 0 iframe 0 0];
                continue;
            end

            subimg = mean_sub(py-dx:py+dx, px-dx:px+dx);

            switch lower(self.noise_model)
                case 'gaussian'
                    % Integrated Gaussian fit (least-squares, fast)
                    res = ls_int_gaussian(subimg, psf_sig);
                    x0 = res.x + px - (dx+1);
                    y0 = res.y + py - (dx+1);
                    loc2(ii,1) = x0;
                    loc2(ii,2) = y0;
                    loc2(ii,3) = res.I0;
                    loc2(ii,4) = psf_sig;
                    loc2(ii,5) = iframe;
                    loc2(ii,8) = res.rmse;
                    loc2(ii,9) = 0;

                case 'poisson'
                    % Radial symmetry initial guess -> Poisson MLE refinement
                    [xc, yc, ~] = radialcenter(subimg);
                    params0 = [dx+1, xc, yc, psf_sig, max(subimg(:))-min(subimg(:))];
                    [A, x0, y0, sigma, ~, negL, variance] = ...
                        gaussfit2DMLE(subimg, [], params0);
                    loc2(ii,1) = x0 + px - (dx+1);
                    loc2(ii,2) = y0 + py - (dx+1);
                    loc2(ii,3) = A;
                    loc2(ii,4) = sigma;
                    loc2(ii,5) = iframe;
                    loc2(ii,8) = negL;
                    loc2(ii,9) = variance;

                otherwise
                    error('SMD:localize', 'Unknown noise_model: %s', self.noise_model);
            end
        end

        % SNR on denoised image (consistent with wavelet pipeline)
        try
            [A1, ~] = calc_snr(loc2, im3);
        catch
            A1 = repmat(3, 1, npeaks);
        end
        loc2(:,6) = A1(:);
        loc2(:,7) = iframe;

        % Keep spots with SNR > 2
        idx = loc2(:,6) > 2;
        discarded = discarded + (npeaks - sum(idx));
        if any(idx)
            spots{iframe} = loc2(idx, :);
        else
            spots{iframe} = [1 1 -1 1 1 1 iframe 0 0];
            stop_frame = min(stop_frame, iframe);
        end
    end
    close(h);
    fprintf('LLR: discarded %d / %d spots\n', discarded, total_spots);
    maxint_img = max(pixel_values, [], 3);
end

% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function [SNRspots, intPeak] = calc_snr(spots, originalIm)
%CALC_SNR  Signal-to-noise ratio of detected spots

    spotRadius = 5;
    halfWindowSize = 8;
    originalIm = double(originalIm);

    nSpots = size(spots, 1);
    imageSize = size(originalIm);
    [columnsInImage, rowsInImage] = meshgrid(1:imageSize(2), 1:imageSize(1));

    spotMask = false(imageSize);
    for spotIdx = 1:nSpots
        curSpotY = round(spots(spotIdx, 2));
        curSpotX = round(spots(spotIdx, 1));
        spotMask = spotMask | ...
            (rowsInImage - curSpotY).^2 + (columnsInImage - curSpotX).^2 <= spotRadius^2;
    end

    bgMask = ~spotMask;
    bgIm   = originalIm .* bgMask;

    intPeak  = zeros(1, nSpots);
    SNRspots = zeros(1, nSpots);

    for spotIdx = 1:nSpots
        curSpotY = round(spots(spotIdx, 2));
        curSpotX = round(spots(spotIdx, 1));

        xMin = max(curSpotX - halfWindowSize, 1);
        xMax = min(curSpotX + halfWindowSize, imageSize(2));
        yMin = max(curSpotY - halfWindowSize, 1);
        yMax = min(curSpotY + halfWindowSize, imageSize(1));

        curSpotBgIm = bgIm(yMin:yMax, xMin:xMax);

        sr = 1;
        [columnsInSpotImage, rowsInSpotImage] = ...
            meshgrid(1:xMax-xMin+1, 1:yMax-yMin+1);
        smask = (rowsInSpotImage - (yMax-yMin)/2 - 1).^2 + ...
                (columnsInSpotImage - (xMax-xMin)/2 - 1).^2 <= sr^2;

        curSpotIm = smask .* originalIm(yMin:yMax, xMin:xMax);

        bgPixelValues  = curSpotBgIm(curSpotBgIm ~= 0);
        meanBgI        = mean(bgPixelValues);
        stdBg          = std(bgPixelValues);

        spotPixelValues = curSpotIm(curSpotIm ~= 0);
        meanSpotI       = mean(spotPixelValues);
        intPeak(spotIdx)  = max(spotPixelValues);
        SNRspots(spotIdx) = round((meanSpotI - meanBgI) / stdBg, 2);
    end
end

function tmp = mynd2read(fname)
    bfr = BioformatsImage(fname);
    tmp = struct([]);
    for ii = 1:bfr.sizeT
        [I, ts, ~] = getPlane(bfr, 1, 1, ii);
        tmp(ii).data = I;
        tmp(ii).ts   = ts(ii);
    end
end

function stack = data2photons(stack)
%DATA2PHOTONS  Convert EMCCD ADU to photon counts
    gain_factor = 1.8;
    offset = 500;
    num_frames = min(1200, length(stack));
    for iframe = 1:num_frames
        I = stack(iframe).data;
        I = (I - offset) ./ gain_factor;
        stack(iframe).data = uint16(I);
    end
end
