function obj = track(obj)

%This function is used to create a utrack-compatible spots structure which
%is then passed on to u-track to connect detected spots in a single-molecule 
%movie into tracks.
%
%
%tracks = uTrack_wrapper(spots, trackingRadius, gapFrames,minLengthBeforeGap, motionType)
%
%
%
%Input:
%   spots               -   Cell array that has as many entries as there are
%                           frames in the image stack. Each cell contains an 
%                           array with at least 2 columns with the x-coordinates 
%                           in the 1st column and y-coordinates in the 2nd column.
%                           Eg. spots{10} = [5.1,10.2; 20.5,30.1] indicates that
%                           there are two spots, one at (x=5.1, y=10.2) and one
%                           at(x=20.5, y=30.1) in the 10th frame of the stack.
%   trackingRadius      -   Maximum allowed distance to connect spots into tracks
%   gapFrames           -   Amount of frames a molecule is allowed to dissapear to  
%                           still be connected into a track.
%   minLengthBeforeGap  -   Minimum amount of frames a track has to exist 
%                           before a connection over a gap frame is allowed.
%   motionType          -   Set to 1 for using linear motion + brownian
%                           motion model or to 0 for using only a brownian
%                           motion model
%                           
%Output:
%   tracks              -   Cell array that has as many entries as there
%                           are tracks. Each cell contains an array with at
%                           least 3 columns [frame,xpos,ypos]. If spot
%                           positions were refined with the fit_spots
%                           function each cell contains an array with 9 columns:
%                           [frame,xpos,ypos,A,BG,sigma_x,sigma_y,angle,exitflag]
%                               1: frame in which spot appears
%                               2: x-position of the fitted spots
%                               3: y-position of the fitted spots
%                               4: A - maximum of the fitted gaussian
%
%
%Example:
%   tracks{1} = [1, 5.1, 10.2, 1044.3;
%                2, 7.2,  9.2, 1030.1;
%                4, 5.5, 30.1, 1050.9]
%   This would be the first track that has been found. The spots connected 
%   to the track appear in the 1st, 2nd and 4th frame. The 3rd frame is a
%   "gap frame".




nFramesAnalyzed = length(spots);
spotsUTrack = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),nFramesAnalyzed,1);

for m = 1:nFramesAnalyzed
    if~isempty(spots{m})
        nCand = size(spots{m},1);
        spotsUTrack(m).xCoord = [spots{m}(:,1) zeros(nCand,1)];     % x-coordinate without using standard deviation of fit
        %         spotsUTrack(m).xCoord = [spots{m}(:,1) spots{m}(:,5)];        % x-coordinate plus standard deviation of fit
        spotsUTrack(m).yCoord = [spots{m}(:,2) zeros(nCand,1)];     % y-coordinate without using standard deviation of fit
        %         spotsUTrack(m).yCoord = [spots{m}(:,2) spots{m}(:,6)];      % y-coordinate plus standard deviation of fit
        spotsUTrack(m).amp    = [spots{m}(:,3) zeros(nCand,1)];     % Fitted spot intensity
    end
end

%% ----------Adopted from u-tracks scriptTrackGeneral.m function-----------

% Jaqaman, K., Loerke, D., Mettlen, M., Kuwata, H., Grinstein, S., Schmid, S. L., & Danuser, G. (2008).
% Robust single-particle tracking in live-cell time-lapse sequences. Nature Methods, 5(8), 695–702. https://doi.org/10.1038/nmeth.1237

gapCloseParam.timeWindow = gapFrames+1; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = minLengthBeforeGap; %minimum length of track segments from linking to be used in gap closing.
gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%parameters

parameters.linearMotion = linearMotion; %use linear motion Kalman filter.
parameters.minSearchRadius = 0; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = trackingRadius; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = linearMotion; %use linear motion Kalman filter.

parameters.minSearchRadius = 0; %minimum allowed search radius.
parameters.maxSearchRadius = trackingRadius; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0.25 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

parameters.ampRatioLimit = []; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^(n-1)).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

%NEW PARAMETER
parameters.gapExcludeMS = 1; %flag to allow gaps to exclude merges and splits

%NEW PARAMETER
parameters.strategyBD = -1; %strategy to calculate birth and death cost

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%verbose state
verbose = 1;

%problem dimension
probDim = 2;

[tracksFromUTrack,~,~] = trackCloseGapsKalmanSparse(spotsUTrack,costMatrices,gapCloseParam,kalmanFunctions,probDim,0,verbose);
%% ------------------------------------------------------------------------

%% Restructure tracking data to obtain compatible results with nearest neighbor algorithm and software output


nTracks = numel(tracksFromUTrack);
runningTrackId = 1;
tracks = cell(nTracks,1);

for trackIdx = 1:nTracks
    
    firstFrame = tracksFromUTrack(trackIdx).seqOfEvents(1,1); %Get frame where track appears first
    
    %Get current track from u-track output structure
    curTrack = [tracksFromUTrack(trackIdx).tracksCoordAmpCG(1:8:end)',tracksFromUTrack(trackIdx).tracksCoordAmpCG(2:8:end)',tracksFromUTrack(trackIdx).tracksCoordAmpCG(4:8:end)'];
    
    trackLength = size(curTrack,1);
    frameIndices = firstFrame:firstFrame+trackLength-1;
    curTrackWithIdx = [frameIndices', curTrack];
    
    %Delete gap frames (rows with NaN)
    curTrackWithIdx = curTrackWithIdx(~isnan(curTrackWithIdx(:,end)),:);
    
    %Sometimes utrack returns tracks with zeros or nans where x or y
    %coordinates should be. Discard these tracks and save all other tracks
    %in tracks structure
    spotsUTrack = curTrackWithIdx(:,2:3);
    if ~any(isnan(spotsUTrack(:))) && ~any(spotsUTrack(:)==0)
        tracks(runningTrackId) = {curTrackWithIdx}; %[frame,x,y,z,amp]
        runningTrackId = runningTrackId + 1;
    end
end

%Delete empty cells
tracks = tracks(~cellfun('isempty', tracks));

%% Get non-linked spots


allSpots = vertcat(spots{:});
frameNumsOfSpots = ones(size(allSpots,1),1);
runningIdx = 1;

for frameIdx = 1:length(spots)
    if isempty(spots{frameIdx})
        continue
    end
    
    nSpotsInCurFrame = size(spots{frameIdx},1);
    frameNumsOfSpots(runningIdx:runningIdx+nSpotsInCurFrame-1) = frameIdx;
    runningIdx = runningIdx+nSpotsInCurFrame;
end

allSpots = [frameNumsOfSpots, allSpots];
allTracks = vertcat(tracks{:});

inSpotsAndTracksIdx = ismember(allSpots(:,2:3),allTracks(:,2:3));
inSpotsAndTracksIdx = inSpotsAndTracksIdx(:,1)&inSpotsAndTracksIdx(:,2);

nonLinkedSpots = num2cell(allSpots(~inSpotsAndTracksIdx,1:4),2);
tracks = [tracks', nonLinkedSpots'];

end