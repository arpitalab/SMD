function all_plot(self)
%ALL_PLOT  Quick visualization of tracks and associated images
%
%   self.all_plot()
%
%   Creates a 2x2 figure with:
%     (1) Tracks overlaid on ROI image
%     (2) Bright-field image
%     (3) Zoomed bright-field image
%     (4) Maximum intensity projection (inverted)
%
%   See also: SMD, localize, track

    figure;

    % --- Panel 1: tracks on ROI image ------------------------------------
    subplot(2,2,1);
    if ~isempty(self.tracks) && ~isempty(self.ROI_image)
        imagesc(self.ROI_image); axis equal tight; colormap gray; hold on;
        for ii = 1:numel(self.tracks)
            plot(self.tracks{ii}(:,1) / self.pixelsize, ...
                 self.tracks{ii}(:,2) / self.pixelsize);
        end
    else
        text(0.5, 0.5, 'No ROI image or tracks', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end

    % --- Panel 2: bright-field image -------------------------------------
    subplot(2,2,2);
    if ~isempty(self.BF_image)
        imagesc(self.BF_image); axis equal tight; colormap gray;
    else
        text(0.5, 0.5, 'No bright-field image', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end

    % --- Panel 3: zoomed bright-field ------------------------------------
    subplot(2,2,3);
    if ~isempty(self.BF_zoom)
        imagesc(self.BF_zoom); axis equal tight; colormap gray;
    else
        text(0.5, 0.5, 'No zoomed bright-field image', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end

    % --- Panel 4: max intensity projection --------------------------------
    subplot(2,2,4);
    if ~isempty(self.maxint)
        imagesc(imcomplement(self.maxint)); axis equal tight; colormap gray;
    else
        text(0.5, 0.5, 'No max intensity projection', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
end
