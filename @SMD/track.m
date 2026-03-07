function track(self)
%TRACK  Link detected spots into trajectories using simpletracker
%
%   self.track()
%
%   Uses the spots stored in self.spots and the parameters in
%   self.tracking_params to connect detections across frames into
%   trajectories. Results are stored in self.tracks.
%
%   See also: SMD, localize, cull_tracks

    spots  = self.spots;
    coords = cell(1, self.stop_frame - 1);

    for iframe = 1:self.stop_frame - 1
        if ~isempty(spots{iframe})
            coords{iframe} = spots{iframe}(:, 1:2);
        else
            coords{iframe} = [];
        end
    end

    [tracks_raw, adjacency_tracks] = simpletracker(coords, ...
        'MaxLinkingDistance', self.tracking_params.trackingRadius, ...
        'MaxGapClosing',     self.tracking_params.gapFrames, ...
        'Debug',             false);

    all_points = vertcat(spots{1:self.stop_frame - 1});
    n_tracks   = numel(tracks_raw);
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
