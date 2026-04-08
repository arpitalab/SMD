function [tracks, adjacency_tracks] = quot_tracker(points, varargin)
%QUOT_TRACKER  Hungarian tracker with augmented cost matrix (quot-style)
%
%   [tracks, adjacency_tracks] = QUOT_TRACKER(points)
%   [tracks, adjacency_tracks] = QUOT_TRACKER(points, 'Key', Value, ...)
%
%   Drop-in replacement for simpletracker with three improvements over the
%   simpletracker default:
%
%     1. Augmented cost matrix: the assignment problem is formulated as a
%        square (n_traj+n_locs) x (n_traj+n_locs) matrix.  Dropping a
%        trajectory or starting a new one costs InitCost, so the solver
%        explicitly weighs linking against fragmentation.
%
%     2. Inline gap handling: unmatched trajectories stay active for up to
%        MaxGapClosing additional frames using their last known position,
%        rather than a separate nearest-neighbor gap-closing pass.
%
%     3. No over-linking: trajectories beyond MaxLinkingDistance are
%        forbidden (cost = Inf), and the dummy-column / dummy-row
%        mechanism prevents forced assignments.
%
%   INPUTS
%     points  - cell array (1 x n_frames); each cell is [n x 2] (or [n x d])
%               containing [x y ...] localisations for that frame.
%               Empty cells / zero-row cells are fine.
%
%   KEY/VALUE PARAMETERS
%     'MaxLinkingDistance'  max jump in pixels (default Inf). Callers
%                           should convert from physical units before
%                           passing (e.g. radius_um / pixelsize).
%     'MaxGapClosing'       max consecutive missed frames before a
%                           trajectory is terminated (default 3)
%     'InitCost'            penalty (in squared pixels) for dropping an
%                           existing trajectory or starting a new one.
%                           Defaults to MaxLinkingDistance^2 (or 1e8 if
%                           MaxLinkingDistance is Inf).  A link is preferred
%                           over a drop/start whenever its squared distance
%                           is less than InitCost.
%     'Debug'               print per-frame progress (default false)
%
%   OUTPUTS
%     tracks          - cell (n_tracks x 1); tracks{i} is [n_frames x 1]
%                       with the 1-based in-frame particle index at each
%                       frame, or NaN where no detection was linked.
%                       Format is identical to simpletracker output.
%     adjacency_tracks - cell (n_tracks x 1); adjacency_tracks{i} is a
%                       column vector of 1-based global indices into
%                       vertcat(points{:}), with NaN frames omitted.
%                       Format is identical to simpletracker output.
%
%   ALGORITHM
%     Frame-to-frame reconnection follows quot/track.py (Heckert lab,
%     github.com/alecheckert/quot).  At each frame transition the cost
%     matrix is:
%
%       W = [ D_sq          | init_cost * 1 ]   <- real trajs
%           [ init_cost * 1 |       0       ]   <- dummy rows (new trajs)
%
%     where D_sq(i,j) is the squared distance from trajectory i's last
%     known position to localisation j (Inf if beyond MaxLinkingDistance).
%     The lower-right block is 0 so that dummy-dummy assignments are free.
%     Munkres minimises the total cost over the square matrix.
%
%     Trajectory i is linked to loc j iff assignment(i) <= n_locs.
%     Otherwise the trajectory's blink counter is incremented; it is
%     terminated when the counter exceeds MaxGapClosing.
%     A dummy row assigned to column j (j <= n_locs) starts a new track.
%
%   DEPENDENCIES
%     munkres.m  (Yi Cao, available on MATLAB File Exchange)
%
%   EXAMPLE
%     coords = { rand(5,2); rand(4,2); rand(6,2) };
%     [T, A] = quot_tracker(coords, 'MaxLinkingDistance', 3, ...
%                                   'MaxGapClosing', 2, 'InitCost', 9);
%
%   See also: simpletracker, hungarianlinker, nearestneighborlinker

    %% -----------------------------------------------------------------------
    %  Parse arguments
    %% -----------------------------------------------------------------------
    ip = inputParser;
    ip.addParamValue('MaxLinkingDistance', Inf,  @isnumeric);   %#ok<NVREPL>
    ip.addParamValue('MaxGapClosing',      3,    @isnumeric);   %#ok<NVREPL>
    ip.addParamValue('InitCost',           [],   @(x) isempty(x) || isnumeric(x)); %#ok<NVREPL>
    ip.addParamValue('Debug',              false, @islogical);  %#ok<NVREPL>
    ip.parse(varargin{:});

    max_dist   = ip.Results.MaxLinkingDistance;
    max_blinks = ip.Results.MaxGapClosing;
    init_cost  = ip.Results.InitCost;
    debug      = ip.Results.Debug;

    % Default InitCost: one squared-pixel more expensive than the worst
    % allowed link, so every in-radius link is preferred over a drop/start.
    if isempty(init_cost)
        if isinf(max_dist)
            init_cost = 1e8;
        else
            init_cost = max_dist ^ 2;
        end
    end

    %% -----------------------------------------------------------------------
    %  Pre-compute global index offsets
    %% -----------------------------------------------------------------------
    n_slices  = numel(points);
    n_cells   = cellfun(@(x) size(x, 1), points);
    cum_cells = [0, cumsum(n_cells(:)')];   % cum_cells(f+1) = last global idx in frame f
    n_total   = cum_cells(end);

    if n_total == 0
        tracks           = {};
        adjacency_tracks = {};
        return;
    end

    % Determine coordinate dimensionality from first non-empty frame
    first_ne  = find(n_cells > 0, 1);
    n_dim_pos = size(points{first_ne}, 2);

    % Stack all positions for O(1) lookup by global index
    all_pos = zeros(n_total, n_dim_pos);
    for fi = 1:n_slices
        if n_cells(fi) > 0
            all_pos(cum_cells(fi)+1 : cum_cells(fi+1), :) = points{fi};
        end
    end

    max_dist_sq = max_dist ^ 2;   % precompute once

    %% -----------------------------------------------------------------------
    %  Initialise trajectories from frame 1
    %  Each trajectory: struct with fields
    %    .indices   - [1 x L] row vector of 1-based global indices
    %    .n_blinks  - consecutive missed frames since last successful link
    %% -----------------------------------------------------------------------
    active    = cell(n_cells(1), 1);
    for ii = 1:n_cells(1)
        active{ii} = struct('indices', cum_cells(1)+ii, 'n_blinks', 0);
    end
    completed = {};

    %% -----------------------------------------------------------------------
    %  Main loop: frame 2 .. n_slices
    %% -----------------------------------------------------------------------
    for fi = 2:n_slices

        n_locs = n_cells(fi);
        n_traj = numel(active);

        if debug
            fprintf('Frame %3d/%d : %3d active trajs  %3d locs\n', ...
                fi, n_slices, n_traj, n_locs);
        end

        % 1-based global indices of localisations in this frame
        loc_gidx = cum_cells(fi) + (1:n_locs);

        % ------------------------------------------------------------------
        % Case A: no localisations – blink all active trajectories
        % ------------------------------------------------------------------
        if n_locs == 0
            surv = {};
            for ti = 1:n_traj
                active{ti}.n_blinks = active{ti}.n_blinks + 1;
                if active{ti}.n_blinks <= max_blinks
                    surv{end+1} = active{ti}; %#ok<AGROW>
                else
                    completed{end+1} = active{ti}; %#ok<AGROW>
                end
            end
            active = surv;
            continue;
        end

        % ------------------------------------------------------------------
        % Case B: no active trajectories – seed new ones from all locs
        % ------------------------------------------------------------------
        if n_traj == 0
            for li = 1:n_locs
                active{end+1} = struct('indices', loc_gidx(li), 'n_blinks', 0); %#ok<AGROW>
            end
            continue;
        end

        % ------------------------------------------------------------------
        % General case: build and solve the augmented assignment problem
        % ------------------------------------------------------------------

        % Last known positions of active trajectories (n_traj x 2)
        last_pos = zeros(n_traj, 2);
        for ti = 1:n_traj
            last_pos(ti, :) = all_pos(active{ti}.indices(end), 1:2);
        end

        % Positions of current-frame localisations (n_locs x 2)
        curr_pos = all_pos(loc_gidx, 1:2);

        % Squared-Euclidean distance matrix D  (n_traj x n_locs)
        D = zeros(n_traj, n_locs);
        for ti = 1:n_traj
            diffs      = curr_pos - repmat(last_pos(ti, :), n_locs, 1);
            D(ti, :)   = sum(diffs .^ 2, 2)';
        end

        % Forbid links beyond the search radius
        D(D > max_dist_sq) = Inf;

        % Build augmented cost matrix W  size (n_traj+n_locs) x (n_traj+n_locs)
        %
        %   W = [ D            | init_cost ]   rows = real trajs
        %       [ init_cost    |     0     ]   rows = dummy (potential new trajs)
        %
        % Lower-right block is 0: dummy-dummy assignments are free and
        % never force spurious links.
        n_aug = n_traj + n_locs;
        W     = zeros(n_aug, n_aug);
        W(1:n_traj,       1:n_locs)       = D;
        W(1:n_traj,       n_locs+1:n_aug) = init_cost;
        W(n_traj+1:n_aug, 1:n_locs)       = init_cost;
        % W(n_traj+1:n_aug, n_locs+1:n_aug) stays 0

        % Solve with Munkres (minimise total cost)
        assignment = munkres(W);   % length n_aug; assignment(i) = col, 0 = unassigned

        % ------------------------------------------------------------------
        % Parse assignments for real trajectory rows (1 .. n_traj)
        % ------------------------------------------------------------------
        surv = {};
        for ti = 1:n_traj
            j = assignment(ti);
            if j >= 1 && j <= n_locs
                % Trajectory i linked to localisation j
                active{ti}.indices(end+1) = loc_gidx(j);
                active{ti}.n_blinks       = 0;
                surv{end+1}               = active{ti}; %#ok<AGROW>
            else
                % Trajectory i dropped (assigned to dummy column)
                active{ti}.n_blinks = active{ti}.n_blinks + 1;
                if active{ti}.n_blinks <= max_blinks
                    surv{end+1}       = active{ti}; %#ok<AGROW>
                else
                    completed{end+1}  = active{ti}; %#ok<AGROW>
                end
            end
        end

        % ------------------------------------------------------------------
        % Parse assignments for dummy rows (n_traj+1 .. n_aug)
        % A dummy row assigned to a real loc column starts a new trajectory.
        % ------------------------------------------------------------------
        for di = n_traj+1:n_aug
            j = assignment(di);
            if j >= 1 && j <= n_locs
                surv{end+1} = struct('indices', loc_gidx(j), 'n_blinks', 0); %#ok<AGROW>
            end
        end

        active = surv;
    end

    % Collect remaining active trajectories
    completed = [completed, active(:)'];

    %% -----------------------------------------------------------------------
    %  Build output in simpletracker format
    %% -----------------------------------------------------------------------
    n_tracks         = numel(completed);
    adjacency_tracks = cell(n_tracks, 1);
    tracks           = cell(n_tracks, 1);

    for ti = 1:n_tracks
        gidx = completed{ti}.indices(:);    % column vector, no NaN

        % adjacency_tracks: global indices, NaN-free (same as simpletracker)
        adjacency_tracks{ti} = gidx;

        % tracks: per-frame in-frame index (1-based), NaN for missed frames
        t_vec = NaN(n_slices, 1);
        for k = 1:numel(gidx)
            g  = gidx(k);
            fi = find(cum_cells(2:end) >= g, 1, 'first');
            t_vec(fi) = g - cum_cells(fi);
        end
        tracks{ti} = t_vec;
    end

end
