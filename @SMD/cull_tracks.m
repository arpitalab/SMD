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
    validTracks = filterFragmentsWithInitialGaps(culled_tracks, 3);
    for ii=1:length(validTracks)
        validTracks{ii}(:,5) = validTracks{ii}(:,5)-validTracks{ii}(1,5)+1;
    end
    self.tracks = validTracks;
end

function validFragments = filterFragmentsWithInitialGaps(fragments, M)
    % fragments: cell array of Nx3 matrices [x, y, frame]
    % M: number of initial frames that must be consecutive
    % Returns: cell array of valid fragments (no gaps in first M frames)
    
    validFragments = {};   
    for i = 1:length(fragments)
        frag = fragments{i};
        if size(frag, 1) < M
            continue; % skip invalid
        end
        
        % Sort by frame number just in case
        frag = sortrows(frag, 5);
        frames = frag(:, 5);
        
        % Find the longest valid suffix that has no gap in its first M frames
        startIdx = findValidStart(frames, M);
        if ~isempty(startIdx)
            validFrag = frag(startIdx:end, :);
            validFragments{end+1} = validFrag;
        end
    end
end

function startIdx = findValidStart(frames, M)
    % Recursively find the earliest index where the next M frames are consecutive
    n = length(frames);
    if n < M
        startIdx = [];
        return;
    end
    
    % Check if first M frames are consecutive
    if all(diff(frames(1:M)) == 1) && frames(1) == frames(1) % already sorted
        startIdx = 1;
        return;
    end
    
    % Find first gap in the sequence
    diffs = diff(frames);
    gapIdx = find(diffs > 1, 1, 'first'); % first position where gap occurs
    
    if isempty(gapIdx)
        % No gap at all — check if first M are consecutive from start
        if all(diff(frames(1:M)) == 1)
            startIdx = 1;
        else
            startIdx = [];
        end
    else
        % Gap found between gapIdx and gapIdx+1
        % Recursively check the fragment starting at gapIdx+1
        startIdx = findValidStart(frames(gapIdx+1:end), M);
        if ~isempty(startIdx)
            startIdx = gapIdx + startIdx; % adjust index
        end
    end
end