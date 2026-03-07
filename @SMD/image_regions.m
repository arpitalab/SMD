function image_regions(self)
    cd(self.dirname);

    tmp = split(self.fname,'.');
    fnamesplit = split(tmp{1},'_');
    fileident = str2num(fnamesplit{end});
    BFfile = sprintf('%03d',fileident+2);%num2str(fileident+2);
    BFzoomfile = sprintf('%03d',fileident+3);%num2str(fileident+3);
    fnamesplit{end} = BFfile;
    xtmp=join(fnamesplit,'_');
    bffilename = [xtmp{1},'.',tmp{2}];
    fnamesplit{end} = BFzoomfile;
    xtmp=join(fnamesplit,'_');
    bfzoomfilename = [xtmp{1},'.',tmp{2}];
    bffile_exists = isfile(bffilename);
    bfzoom_exists = isfile(bfzoomfilename);
    if(bffile_exists) 
        d1 = dir(bffilename);
    end
    if(bfzoom_exists) 
        d2 = dir(bfzoomfilename);
    end
    if(bffile_exists && bfzoom_exists && d1.bytes<1e7 && d2.bytes<1e7)
        if(tmp{2} == 'nd2')
            bfr=BioformatsImage(bffilename);
            [I, ~, ~] = getPlane(bfr, 1,1,1);
            self.BF_image = I;
            bfr=BioformatsImage(bfzoomfilename);
            [I, ~, ~] = getPlane(bfr, 1,1,1);
            self.BF_zoom = I;
        else
           I=imread(bffilename);
           self.BF_image = I;
           I=imread(bfzoomfilename);
           self.BF_zoom = I;
        end       
    else
        parts = split(self.fname,'.');
        query = false;
        f = msgbox({self.fname;'Choose the appropriate brightfield image'});
        [file,~] = uigetfile(['*.',parts{2}]);
        close(f);
        %figure;
        if(parts{2} == 'nd2')
            bfr=BioformatsImage(file);
            [I, ~, ~] = getPlane(bfr, 1,1,1);
        else
            I=imread(file);
        end
        self.BF_image = I;
        %delete(f);
        f = msgbox({self.fname;'Choose the appropriate brightfield zoom image'});
        [file,~] = uigetfile(['*.',parts{2}]);
        close(f);
        %figure;
        if(parts{2} == 'nd2')
            bfr=BioformatsImage(file);
            [I, ~, ~] = getPlane(bfr, 1,1,1);
        else
            I=imread(file);
        end
        self.BF_zoom = I;
    end

end