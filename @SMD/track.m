function track(self)
%TRACK  Link detected spots into trajectories
%
%   self.track()
%
%   Connects detections stored in self.spots into trajectories using the
%   linker selected by self.tracking_params.tracker:
%
%     'simpletracker' (default) - frame-to-frame Hungarian assignment
%                                  followed by a nearest-neighbour
%                                  gap-closing pass (Tinevez 2011).
%
%     'quot'                    - augmented-cost-matrix Hungarian tracker
%                                  inspired by quot/track.py (Heckert lab).
%                                  Drops/starts cost InitCost = R^2 px^2;
%                                  gaps are handled inline (no second pass).
%
%   Tracking parameters (self.tracking_params fields):
%     trackingRadius      - max frame-to-frame jump in pixels (default 2)
%     gapFrames           - max gap length in frames         (default 3)
%     minLengthBeforeGap  - min track length to keep         (default 4)
%     tracker             - linker to use (see above)        (default 'simpletracker')
%
%   Results are stored in self.tracks.
%
%   See also: SMD, localize, cull_tracks, quot_tracker, simpletracker

    spots  = self.spots;
    coords = cell(1, self.stop_frame - 1);

    for iframe = 1:self.stop_frame - 1
        if ~isempty(spots{iframe})
            coords{iframe} = spots{iframe}(:, 1:2);
        else
            coords{iframe} = zeros(0, 2);
        end
    end

    % Resolve tracker choice (default to simpletracker if field absent)
    if isfield(self.tracking_params, 'tracker')
        tracker_name = self.tracking_params.tracker;
    else
        tracker_name = 'simpletracker';
    end

    % Convert trackingRadius from µm to pixels for both trackers
    radius_px = self.tracking_params.trackingRadius / self.pixelsize;

    switch lower(tracker_name)

        case 'simpletracker'
            [~, adjacency_tracks] = simpletracker(coords, ...
                'MaxLinkingDistance', radius_px, ...
                'MaxGapClosing',      self.tracking_params.gapFrames, ...
                'Debug',              false);

        case 'quot'
            [~, adjacency_tracks] = quot_tracker(coords, ...
                'MaxLinkingDistance', radius_px, ...
                'MaxGapClosing',      self.tracking_params.gapFrames, ...
                'Debug',              false);

        otherwise
            error('SMD:track', ...
                'Unknown tracker ''%s''. Use ''simpletracker'' or ''quot''.', ...
                tracker_name);
    end

    all_points = vertcat(spots{1:self.stop_frame - 1});
    n_tracks   = numel(adjacency_tracks);
    result     = {};
    k = 1;

    for ii = 1:n_tracks
        track_points = all_points(adjacency_tracks{ii}, :);
        if size(track_points, 1) >= self.tracking_params.minLengthBeforeGap
            result{k} = track_points;
            result{k}(:, 1:2) = result{k}(:, 1:2) * self.pixelsize;
            k = k + 1;
        end
    end

    self.tracks = result;
end
