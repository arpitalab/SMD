classdef SMD < handle
%SMD  Single-Molecule Data class for localization and tracking
%
%   obj = SMD(dirname, fname, exposure_time, frame_rate)
%
%   The SMD class provides a complete pipeline for single-molecule
%   microscopy analysis: spot detection, sub-pixel localization, tracking,
%   and visualization.
%
%   Two detection/fitting pipelines are available, selected via the
%   detection_method property:
%
%     'llr'      - (default) Log-likelihood ratio test + radial symmetry or
%                  integrated Gaussian fitting
%     'wavelet'  - a-trous wavelet filter + Crocker-Grier peak finding +
%                  Poisson MLE Gaussian fitting
%
%   Example:
%       smd = SMD('/data/experiment/', 'movie.tif', 0.02, 50);
%       smd.localize();
%       smd.track();
%       smd.all_plot();
%
%   See also: localize, track, all_plot

    % ----- User-configurable settings (public) ----------------------------
    properties
        dirname         char                        % parent directory for raw data
        fname           char                        % TIFF / ND2 filename
        camera          char   = 'EMCCD'            % camera type ('EMCCD' or 'CMOS')
        pixelsize       double = 0.16               % pixel size in microns
        exposure_time   double                      % exposure time in seconds
        frame_rate      double                      % acquisition frame rate (Hz)
        localization_box       double = 3           % half-width of fitting box (pixels)
        localization_threshold double = 1.75        % wavelet detection threshold (1-2)
        tracking_params struct = struct( ...        % tracking parameters
            'trackingRadius',    0.32, ...          %   max jump distance in µm
            'gapFrames',         3, ...
            'minLengthBeforeGap',4, ...
            'linearMotion',      0, ...
            'tracker',           'quot')            % 'simpletracker' or 'quot'
        detection_method char   = 'llr'             % 'wavelet' or 'llr'
        noise_model      char   = 'gaussian'        % 'gaussian' (ls_int_gaussian) or 'poisson' (radialcenter + gaussfit2DMLE)
        psf_sigma        double = 1.3               % PSF std dev in pixels (used by LLR and integrated Gaussian)
        llr_params       struct = struct( ...       % LLR detection parameters
            'k', 1.0, ...                           %   Gaussian kernel width
            'w', 9,   ...                           %   window size (odd integer)
            't', 20.0)                              %   score threshold
    end

    % ----- Result properties (private set, public read via getters) -------
    properties (SetAccess = private)
        spots           cell   = {}                 % cell array of localized spots per frame
        tracks          cell   = {}                 % cell array of linked trajectories
        relative_tracks cell   = {}                 % drift-corrected trajectories
        stop_frame      double = Inf                % last analysable frame
        timestamps      double = []                 % per-frame timestamps (ND2 only)
        maxint                                      % maximum intensity projection
        num_frames      double = 0                  % number of analyzed frames
        ROIs            cell   = {}                 % cell array of ROI polygons
        ROI_image                                   % fluorescence image for ROI detection
        BF_image                                    % bright-field image
        BF_zoom                                     % zoomed bright-field image
    end

    methods
        % ----- Constructor ------------------------------------------------
        function obj = SMD(dirname, fname, exposure_time, frame_rate)
        %SMD  Construct an SMD object
        %   obj = SMD(dirname, fname, exposure_time, frame_rate)
            arguments
                dirname       char
                fname         char
                exposure_time double
                frame_rate    double
            end
            obj.dirname        = dirname;
            obj.fname          = fname;
            obj.exposure_time  = exposure_time;
            obj.frame_rate     = frame_rate;
        end

        % ----- Getter methods for private properties ----------------------
        function v = get_spots(self),           v = self.spots;           end
        function v = get_tracks(self),          v = self.tracks;          end
        function v = get_relative_tracks(self), v = self.relative_tracks; end
        function v = get_stop_frame(self),      v = self.stop_frame;      end
        function v = get_timestamps(self),      v = self.timestamps;      end
        function v = get_maxint(self),          v = self.maxint;          end
        function v = get_num_frames(self),      v = self.num_frames;      end
        function v = get_ROIs(self),            v = self.ROIs;            end
        function v = get_ROI_image(self),       v = self.ROI_image;       end
        function v = get_BF_image(self),        v = self.BF_image;        end
        function v = get_BF_zoom(self),         v = self.BF_zoom;         end

        % ----- External method declarations -------------------------------
        localize(self);
        get_roi(self);
        track(self);
        cull_tracks(self);
        image_regions(self);
        remove_relative_motion(self);
        all_plot(self);
    end
end
