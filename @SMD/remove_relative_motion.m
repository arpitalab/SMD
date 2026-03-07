function remove_relative_motion(self)
    % make sure that column 10 exists. If not then rerun to identify the
    % specific ROI that the track belongs to and populate it
    data={};
    tracks = {};
    for itrack=1:numel(self.tracks)
        data{itrack}=self.tracks{itrack}(:,[1:2,7,10]);
    end
    for iroi=1:length(self.ROIs)
        k=0;
        ix=[];
        for itrack=1:numel(data)
            if(data{itrack}(1,end)==iroi)
                k=k+1;
                ix(k)=itrack;
            end
        end
        if(length(ix)>10)
            try
                [forest, ~, forestInfo] = maxSpanningForest(data(ix), 10, 2);
                [relTracks, ~] = relativeTracksFromForest(data(ix), forest);
                if forestInfo.numComponents > 1
                    fprintf('  ROI %d: %d components\n', iroi, forestInfo.numComponents);
                end
                if(numel(relTracks)>5)
                    tracks=[tracks,relTracks];
                end
            catch
                disp('cannot do correction');
                continue;
            end
        end
    end
    self.relative_tracks = tracks;
end