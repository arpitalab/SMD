function cull_tracks(self)
% culls tracks by discarding any that are not within any ROIs 
% ROIs are in self.ROIs (from the flood fill image)
% identity of ROI (cell) 
% added 11/03/2025
% Also checks if each trajectory has no gaps in the first M frames
% if there is a gap, all frames upto the first gap are discarded and
% recursively checked to only keep fragments that satisfy above criterion
% added 10/31/2025

    if(isempty(self.ROIs))
        disp('calculating ROIs')
        self.get_roi;
    end
    ntracks = length(self.tracks);
    culled_tracks = {};
    minLengthBeforeGap = self.tracking_params.minLengthBeforeGap;
    k=1;
    for ii=1:ntracks
        tmp=self.tracks{ii}(:,1:2)/self.pixelsize;
        in=zeros(size(tmp,1),length(self.ROIs));
        % check to see if points in the track are in any of the ROIs
        % (nuclei)
        for jj=1:length(self.ROIs)
            in(:,jj)=inpolygon(tmp(:,1),tmp(:,2),self.ROIs{jj}(:,1),self.ROIs{jj}(:,2));
        end
        % if most of the track is within the ROI, keep it
        maximum = 0;
        ROI_belong = 0;
        for jj=1:length(self.ROIs)
            if(nnz(in(:,jj))>0.5*size(tmp,1))
                maximum = max(nnz(in(:,jj)),maximum);
                ROI_belong = jj;
            end
        end
        if(ROI_belong > 0)
          culled_tracks{k}=self.tracks{ii};
          culled_tracks{k}(:,10) = ROI_belong;
          k=k+1; 
        end
    end
    % Use shared utility (frame numbers are in column 5 for SMD tracks)
    validTracks = TrackUtils.filterFragmentsWithInitialGaps(culled_tracks, 5, minLengthBeforeGap);
    for ii=1:length(validTracks)
        validTracks{ii}(:,5) = validTracks{ii}(:,5)-validTracks{ii}(1,5)+1;
    end
    self.tracks = validTracks;
end