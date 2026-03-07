function get_roi(self)
    if(isempty(self.ROI_image))
        cd(self.dirname);
        tmp = split(self.fname,'.');
        fnamesplit = split(tmp{1},'_');
        fileident = str2num(fnamesplit{end});
        ROIfile = sprintf('%03d',fileident+1);%num2str(fileident+1);
        fnamesplit{end} = ROIfile;
        xtmp=join(fnamesplit,'_');
        file = [xtmp{1},'.',tmp{2}];
        file_exists = isfile(file);
        if(file_exists) 
            d1 = dir(file);
        end
        if(file_exists && d1.bytes<1e7)
            if(tmp{2} == 'nd2')
                bfr=BioformatsImage(file);
                [I, ~, ~] = getPlane(bfr, 1,1,1);
            else
               I=imread(file);
            end
        else
            parts = split(self.fname,'.');
            query = false;
            f = msgbox({self.fname;'Choose the appropriate JF549 nucleus image'});
            [file,~] = uigetfile(['*.',parts{2}]);
            close(f);
            figure;
            if(parts{2} == 'nd2')
                bfr=BioformatsImage(file);
                [I, ~, ~] = getPlane(bfr, 1,1,1);
            else
                I=imread(file);
            end
        end
        self.ROI_image = I;
    else
        I = self.ROI_image;
        query = false;
    end
    imagesc(I);axis equal;colormap gray;
    hold on;
    M=[1 -2 1; -2 4 -2; 1 -2 1]; 
    [H,W]=size(I);
    Sigma=sum(sum(abs(conv2(I, M))));
    Sigma=Sigma*sqrt(0.5*pi)./(6*(W-2)*(H-2));
    thresh = mean(I(:))+1.2*Sigma;
    I(I<thresh)=0;
    BW1=bwmorph(I,'majority',2);
    %BW1=imdilate(BW1,se1);
    BW1=imfill(BW1,'holes');
    D = -bwdist(~BW1);
    D = imgaussfilt(D, 1); % Smooth the distance transform slightly

    % Modify the distance transform to avoid over-segmentation
    D = imhmin(D, 1); % Suppress minima shallower than 1 to prevent over-segmentation

    % Apply watershed transform
    L = watershed(D);
    BW1(L == 0) = 0; 
    outI = bwareaopen(BW1,2000); 
    bb=bwboundaries(outI);
    for ii=1:length(bb)
        plot(bb{ii}(:,2),bb{ii}(:,1),'r','linewidth',2);
    end
    query = false;
    while ~query
        response = input('is the contour satisfactory (y=1 or n=0)? ');
        if(response == 1)
             query = true;
             for ii=1:length(bb)
                 bb_1{ii}(:,1) = bb{ii}(:,2);
                 bb_1{ii}(:,2) = bb{ii}(:,1);
             end
             self.ROIs = bb_1;
             % need to flip 1 and 2
             close(gcf);
        else
            bb = {};
            done = false;
            k=1;
            disp('draw contour by hand')
            while ~done
                h = drawpolygon;
                bb{k} = h.Position;
                response2 = input('more polygons to draw (y=1 or n=0)?');
                if(response2 == 1)
                    done = false;
                    k=k+1;
                else
                    done = true;
                end
            end
            query = true;
            close(gcf);
            % hand drawn ROIs are correct
            self.ROIs = bb;
        end
    end
end