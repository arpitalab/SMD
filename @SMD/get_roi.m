function get_roi(self)
%GET_ROI  Segment nuclei from fluorescence image using gradient flow
%
%   self.get_roi()
%
%   Detects nuclear ROIs using a Cellpose-inspired gradient flow algorithm:
%     1. Otsu thresholding for robust foreground detection
%     2. Distance transform gradient field
%     3. Vectorized gradient flow simulation (pixels flow toward centers)
%     4. Convergence-point clustering for instance segmentation
%
%   Falls back to manual polygon drawing if the user rejects the result.
%   ROIs are stored as cell array of Nx2 [x, y] polygons in self.ROIs,
%   compatible with inpolygon() in cull_tracks.
%
%   See also: SMD, cull_tracks

    % --- Load image if needed ---------------------------------------------
    if isempty(self.ROI_image)
        I = load_roi_image(self);
        self.ROI_image = I;
    else
        I = self.ROI_image;
    end

    % --- Segment nuclei via gradient flow ---------------------------------
    labels = gradient_flow_segmentation(double(I));

    % --- Extract boundaries per nucleus -----------------------------------
    bb = {};
    region_ids = setdiff(unique(labels(:)), 0);
    for k = 1:numel(region_ids)
        mask_k = labels == region_ids(k);
        bdry = bwboundaries(mask_k, 8, 'noholes');
        if ~isempty(bdry)
            % bwboundaries returns [row, col] = [y, x]; convert to [x, y]
            bb{end+1} = [bdry{1}(:,2), bdry{1}(:,1)]; %#ok<AGROW>
        end
    end

    % --- Display and confirm ----------------------------------------------
    figure;
    imagesc(I); axis equal tight; colormap gray; hold on;
    for ii = 1:numel(bb)
        plot(bb{ii}(:,1), bb{ii}(:,2), 'r', 'LineWidth', 2);
    end

    query = false;
    while ~query
        response = input('Is the contour satisfactory (y=1 or n=0)? ');
        if response == 1
            query = true;
            self.ROIs = bb;
            close(gcf);
        else
            % Manual fallback
            bb = {};
            done = false;
            k = 1;
            disp('Draw contour by hand');
            while ~done
                h = drawpolygon;
                bb{k} = h.Position;
                response2 = input('More polygons to draw (y=1 or n=0)? ');
                if response2 == 1
                    k = k + 1;
                else
                    done = true;
                end
            end
            query = true;
            self.ROIs = bb;
            close(gcf);
        end
    end
end

% =========================================================================
%  GRADIENT FLOW INSTANCE SEGMENTATION
% =========================================================================
function labels = gradient_flow_segmentation(I)
%GRADIENT_FLOW_SEGMENTATION  Cellpose-inspired nucleus segmentation
%
%   labels = gradient_flow_segmentation(I)
%
%   Returns a label matrix where each nucleus has a unique integer ID.
%   Background pixels are 0.

    [H, W] = size(I);

    % --- 1. Robust foreground detection -----------------------------------
    % Normalize to [0, 1] for Otsu
    I_norm = I - min(I(:));
    mx = max(I_norm(:));
    if mx > 0
        I_norm = I_norm / mx;
    end

    % Gaussian blur to suppress noise
    I_blur = imgaussfilt(I_norm, 2);

    % Otsu threshold
    level = graythresh(I_blur);
    foreground = imbinarize(I_blur, level);

    % Morphological cleanup
    foreground = imopen(foreground, strel('disk', 3));       % remove speckle
    foreground = imclose(foreground, strel('disk', 5));       % close gaps
    foreground = imfill(foreground, 'holes');                 % fill holes
    foreground = bwareaopen(foreground, 500);                 % remove small debris

    % --- 2. Distance transform + normalized gradient field ----------------
    D = bwdist(~foreground);
    D = imgaussfilt(D, 2);  % smooth to avoid noisy gradients at boundary

    [Gy, Gx] = gradient(D);

    % Normalize to unit vectors (only inside foreground)
    mag = sqrt(Gx.^2 + Gy.^2) + eps;
    Gx = Gx ./ mag;
    Gy = Gy ./ mag;

    % Zero out gradient outside foreground
    Gx(~foreground) = 0;
    Gy(~foreground) = 0;

    % --- 3. Gradient flow simulation (vectorized) -------------------------
    [row_all, col_all] = find(foreground);
    n_fg = numel(row_all);

    if n_fg == 0
        labels = zeros(H, W);
        return;
    end

    % Current positions (floating point)
    py = double(row_all);
    px = double(col_all);

    step_size = 0.5;
    max_steps = 200;
    converge_tol = 0.1;

    for step = 1:max_steps
        % Interpolate gradient at current positions
        dx = interp2(Gx, px, py, 'linear', 0);
        dy = interp2(Gy, px, py, 'linear', 0);

        % Update positions
        px = px + step_size * dx;
        py = py + step_size * dy;

        % Clamp to image bounds
        px = max(1, min(W, px));
        py = max(1, min(H, py));

        % Check convergence
        if max(abs(dx * step_size)) < converge_tol && ...
           max(abs(dy * step_size)) < converge_tol
            break;
        end
    end

    % --- 4. Cluster converged positions → instance labels -----------------
    % Round to nearest pixel
    py_r = round(py);
    px_r = round(px);
    py_r = max(1, min(H, py_r));
    px_r = max(1, min(W, px_r));

    % Build a label image from converged positions
    converge_img = zeros(H, W);
    idx_lin = sub2ind([H, W], py_r, px_r);
    converge_img(idx_lin) = 1;  % mark converged-to pixels

    % Label connected clusters of convergence points
    CC_conv = bwconncomp(converge_img, 8);

    % Map each convergence cluster back to the original foreground pixels
    % Build lookup: converged linear index → cluster ID
    cluster_map = zeros(H, W);
    for k = 1:CC_conv.NumObjects
        cluster_map(CC_conv.PixelIdxList{k}) = k;
    end

    % Assign each foreground pixel to the cluster its flow converged to
    labels = zeros(H, W);
    for ii = 1:n_fg
        cid = cluster_map(py_r(ii), px_r(ii));
        if cid > 0
            labels(row_all(ii), col_all(ii)) = cid;
        end
    end

    % --- 5. Post-processing -----------------------------------------------
    % Re-label connected components (a cluster might split spatially)
    labels_clean = zeros(H, W);
    new_id = 0;
    unique_ids = setdiff(unique(labels(:)), 0);
    for k = 1:numel(unique_ids)
        mask_k = labels == unique_ids(k);
        CC_k = bwconncomp(mask_k, 8);
        for j = 1:CC_k.NumObjects
            if numel(CC_k.PixelIdxList{j}) >= 500  % min nucleus area
                new_id = new_id + 1;
                labels_clean(CC_k.PixelIdxList{j}) = new_id;
            end
        end
    end

    labels = labels_clean;
end

% =========================================================================
%  IMAGE LOADING (preserved from original)
% =========================================================================
function I = load_roi_image(self)
%LOAD_ROI_IMAGE  Load the fluorescence image for ROI detection

    if ~endsWith(self.dirname, '/') && ~endsWith(self.dirname, '\')
        dn = [self.dirname '/'];
    else
        dn = self.dirname;
    end

    tmp = split(self.fname, '.');
    fnamesplit = split(tmp{1}, '_');
    fileident = str2double(fnamesplit{end});
    ROIfile = sprintf('%03d', fileident + 1);
    fnamesplit{end} = ROIfile;
    xtmp = join(fnamesplit, '_');
    file = fullfile(dn, [xtmp{1}, '.', tmp{2}]);
    file_exists = isfile(file);

    if file_exists
        d1 = dir(file);
    end

    if file_exists && d1.bytes < 1e7
        if strcmpi(tmp{2}, 'nd2')
            bfr = BioformatsImage(file);
            [I, ~, ~] = getPlane(bfr, 1, 1, 1);
        else
            I = imread(file);
        end
    else
        parts = split(self.fname, '.');
        f = msgbox({self.fname; 'Choose the appropriate nucleus image'});
        [file, ~] = uigetfile(['*.', parts{2}]);
        close(f);
        if strcmpi(parts{2}, 'nd2')
            bfr = BioformatsImage(file);
            [I, ~, ~] = getPlane(bfr, 1, 1, 1);
        else
            I = imread(file);
        end
    end
end
