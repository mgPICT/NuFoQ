function [ signal ] = segment3Dfcm( img3D, dims, spotsizedef, resXY, resZ, Sig_qual)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        output = double(img3D);
        fprintf('3D segmentation in progress\n');
        
        spotsize = round(sqrt(2)*spotsizedef);
        
        %% Removal of spots
%         disk = false(spotsize*2+1, spotsize*2+1, spotsize*2+1);
%         [xx, yy, zz] = meshgrid(1:spotsize*2+1, 1:spotsize*2+1, 1:spotsize*2+1);
%         disk = (xx - spotsize-1).^2 + (yy - spotsize-1).^2 + (zz - spotsize-1).^2 <= spotsize.^2;
%         output = imopen(output, disk);
        
        %% 2D localisation of nuclei
        sumproj = sum(output, 3);
%         if Sig_qual
%             sumproj = (sumproj - min(sumproj(:)))/(max(sumproj(:)) - min(sumproj(:)));
%             sumproj = imadjust(sumproj, [], [], .5);
%         end
%         bckgproj = imopen(sumproj, strel('disk', spotsize, 0));
%         bckgproj = imfilter(sumproj, fspecial('gaussian', spotsize*5, spotsize));
        %bckgproj = imtophat(sumproj, strel('disk', spotsize*2, 0));
        if Sig_qual
            bckgproj = sumproj;
        else
            bckgproj = imopen(sumproj, strel('disk', spotsize, 0));
        end
        
        imfilt = zeros(dims(1:2));

        sigmax = 10:5:20;
        sigmay = 10:5:20;
        theta = 0:pi/4:pi/4;
        
        reps=min(10, length(sigmax));
        fprintf(['\t[',repmat(' ',1,reps),']']);
        did = 0;

        for i_sigx = 1:length(sigmax)
            for i_sigy = 1:length(sigmay)
               for i_theta = 1:length(theta)
                % Generate 3D gaussian kernel
                    cellsize = max(max(sigmax), max(sigmay));

                    wxyz = ceil(cellsize*7);
                    cxyz = wxyz/2+0.5;
                    [xx, yy] = meshgrid(1:wxyz);

                    a = cos(theta(i_theta))^2/(2*sigmax(i_sigx)^2) + sin(theta(i_theta))^2/(2*sigmay(i_sigy)^2);
                    b = -sin(2*theta(i_theta))/(4*sigmax(i_sigx)^2) + sin(2*theta(i_theta))^2/(4*sigmay(i_sigy)^2);
                    c = sin(theta(i_theta))^2/(2*sigmax(i_sigx)^2) + cos(theta(i_theta))^2/(2*sigmay(i_sigy)^2);

                    filter = exp(-(a.*(xx - cxyz).^2 + 2*b.*(xx - cxyz).*(yy - cxyz) + c.*(yy - cxyz).^2));
        %             filter = exp(-( ((xx-cxyz).^2)/(2*sigmax(i_sigx)^2) + ((yy-cxyz).^2)/(2*sigmay(i_sigy)^2)));
                    filter = filter./sum(filter(:));
                    filter = del2(filter);


                    imfilt = imfilt + imfilter(double(bckgproj), filter, 'symmetric').*(1+log(sigmax(i_sigx))).*(1+log(sigmay(i_sigy)));
                end
            end

            if sum(i_sigx == round(1:(length(sigmax)/reps):length(sigmax))) == 1
                fprintf(repmat('\b',1,reps+1-did)) %send the cursor back to the start
                did = did +1;
                fprintf('-');
                fprintf(repmat(' ',1,reps-did));
                fprintf(']');
            end
        end

        seg = 1-mat2gray(imfilt);
        
        locs = imregionalmax(seg); %medfilt2(seg, [cellsize cellsize]));
        if Sig_qual
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 2) = 0;
            segseg = segseg > 1;
        else
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 3) = 0;
            segseg = segseg > 1;
        end
        
        % First separation step
        ws = watershed(-seg);
        segseg(ws == 0) = 0;
        
        % Check if noisy segmentation (big connected objects)
        CC = bwconncomp(segseg);
        LL = regionprops(CC, 'Area', 'Solidity');
        if corr(cat(1,LL.Area), cat(1,LL.Solidity)) < -0.5
            % Possibly bad segmentation - try to separate
            fm = FastFCMeans(im2uint8(seg), 3, 1.1);
            ws = watershed(3-fm);
            segseg(ws == 0) = false;
            
            CC = bwconncomp(segseg);
            LL = regionprops(CC, 'Area', 'Solidity');
            if corr(cat(1,LL.Area), cat(1,LL.Solidity)) < -0.5
                % Gnaan
                locs = imregionalmax(seg);
                segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
                locs(segseg < 3) = 0;
                segseg = segseg == 3;
            end
        end
        clear CC LL
        
        % Check if noisy segmentation (background kept)
        if sum(segseg(:))/double(dims(1))/double(dims(2)) > .8 % 80% of image (in 2D) segmented
            locs = imregionalmax(seg);
            segseg = FastFCMeans(im2uint8(seg), 3, 1.1);
            locs(segseg < 3) = 0;
            segseg = segseg == 3;
        end
        
        
%         gmmobj = gmdistribution.fit(seg(:), 2);
%         [~, id] = max(gmmobj.mu);
%         seg2 = reshape(gmmobj.cluster(seg(:)) == id, size(seg));
%         tset = imreconstruct(locs, seg2, 4);

        [px, py] = find(locs);
        candidates = [px, py];
        T = clusterdata(candidates, 'criterion', 'distance', 'cutoff', spotsize);
        nbLocs = max(T);
        newlocs = zeros(nbLocs, 2);
        for i = 1:nbLocs
            newlocs(i,:) = round(mean(candidates(T == i, :), 1));
        end

        locs = zeros(dims(1), dims(2));
        locs(sub2ind([dims(1), dims(2)], newlocs(:,1), newlocs(:,2))) = 1;
        
        ridges = watershed(~segseg);
        map = double(ridges>0);
        map(locs(:)>0) = 2;
        
        labels = watershed(-map);

%         dmap = bwdist(locs)
        
    %     dmap(segseg == 0) = Inf;
%         labels = watershed(dmap);%locs==0);

        fprintf('3D segmentation of nuclei\n');
        clear imfilt locs seg

        wz = ceil(spotsize*resZ/resXY*2);

        wxyz = ceil(spotsize*2);
        cxyz = wxyz/2+0.5;
        cz = wz/2+0.5;
        [xx, yy, zz] = meshgrid(1:wxyz, 1:wxyz, 1:wz);
        filter = exp(-( ((xx-cxyz).^2)/(spotsize^2) + ...
                        ((yy-cxyz).^2)/(spotsize^2) + ...
                        ((zz-cz).^2)/((spotsize*resZ/resXY)^2) ));
        filter = filter./sum(filter(:));

    %     alfred = img3D;
        output = imfilter(output, filter, 'symmetric');
        output(~repmat(imdilate(segseg, strel('disk', 10, 0)), 1, 1, dims(3))) = nan;

        labs = double(repmat(labels, 1, 1, dims(3)));
        
        outputadj = im2uint8((output - nanmin(output(:)))./(nanmax(output(:)) - nanmin(output(:))));
        pixlist = outputadj(~isnan(output));
        if Sig_qual
            classified = FastFCMeans(pixlist, 3) == 3;
            bobi = false(dims);
            bobi(~isnan(output)) = classified;
        else
            classified = FastFCMeans(pixlist, 3) > 1;
            bobi = false(dims);
            bobi(~isnan(output)) = classified;
        end
        CC = bwconncomp(bobi);
        if CC.NumObjects < max(labels(:))*0.05
            bobi = imdilate(bobi, strel('disk', 5, 0));
            output(bobi) = nan;
        end

        bobi = nan(dims);
        for ilabs = 1:max(labels(:))
            subpixs = find(labs == ilabs);
            vals = output(subpixs);
            subpixs(isnan(vals)) = [];
            vals(isnan(vals)) = [];
 
            vals = ((vals-min(vals(:)))./(max(vals(:))-min(vals(:))));
            try
                valseg = otsu(vals, 3);
%                 valseg = imquantize(vals, multithresh(vals, 3));
            catch err
                if strcmp(err.message, 'Edge vector must be monotonically non-decreasing.')
                    continue
                end
            end
            
            if Sig_qual % [mod 06/11/19]
                if sum(valseg == 3) < 4/3*pi*spotsize^3*2
                    bobi(subpixs) = valseg > 1;
                else
                    bobi(subpixs) = valseg > 2;
                end
            else % [/mod 06/11/19]
                bobi(subpixs) = valseg;
            end
        end
        
        if Sig_qual
            bobi = bwareaopen(bobi >0, 25); %bobi = bwareaopen(bobi >2, 25); [mod 06/11/19] % >0 since bwareaopen seems to consider NaNs as true...
        else
            bobi = bwareaopen(bobi >2, 25);
        end

    %     nucsbw = double(bobi).*labs;

        conts = zeros(dims);
        for i_z = 1:dims(3)
            conts(:,:,i_z) = edge(bobi(:,:,i_z));
        end
        conts = conts.*labs;

        nbNucs = max(labs(:));
        fprintf('\t Nuclei located and labeled: %d found !\n', nbNucs);

        if nbNucs > 255
            signal = zeros(dims, 'uint16');
        else
            signal = zeros(dims, 'uint8');
        end
        for i_nuc = 1:nbNucs
            pointMatrix = find(conts == i_nuc);
            if length(pointMatrix) < 4
                continue
            end
            [px, py, pz] = ind2sub(dims, pointMatrix);
            dt = DelaunayTri([px, py, pz]);  %# Create a Delaunay triangulation
            if isempty(dt.Triangulation)
                continue
            end

            minx = min(px); maxx = max(px);
            miny = min(py); maxy = max(py);
            minz = min(pz); maxz = max(pz);
            [YY,XX,ZZ] = meshgrid(miny:maxy, minx:maxx, minz:maxz);   %# Create a mesh of coordinates for your volume
            simplexIndex = pointLocation(dt,XX(:),YY(:),ZZ(:));  %# Find index of simplex that
                                                              %#   each point is inside
            maskloc = ~isnan(simplexIndex);    %# Points outside the convex hull have a
                                            %#   simplex index of NaN
            maskloc = double(reshape(maskloc,size(XX)));

            locs = signal(minx:maxx, miny:maxy,minz:maxz) == 0;
            maskloc(locs) = maskloc(locs).*i_nuc;
            
            locs_sig = signal(minx:maxx, miny:maxy,minz:maxz);
            maskloc(~locs) = locs_sig(~locs);

            signal(minx:maxx, miny:maxy,minz:maxz) = maskloc;
        end
        
        signal = imclearborder(signal, 8);
        
        fprintf('\t Nuclei fully segmented\n');
        
        % Check for empty image
        valid = zeros(nbNucs, 1);
        bg = median(img3D(signal == 0));
        for i = 1:nbNucs
            valid(i) = median(img3D(signal == i));
        end
        
        if sum(valid(~isnan(valid)) > bg*1.05) == 0
            fprintf('\tThe image does not contain any nuclei\n');
            signal = zeros(dims);
        end
        
        %  Relabel if needed [mod 12/11/19]
        gnaan = unique(signal(signal>0));
        if length(gnaan) < nbNucs
            nbNucs = length(gnaan);
            for i = 1:nbNucs
                signal(signal == gnaan(i)) = i;
            end
            fprintf('Nuclei re-labeled, correct number of nuclei segmented: %d\n', nbNucs);
        end % [/mod 12/11/19]

end
