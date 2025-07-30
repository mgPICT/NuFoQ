function [spots, signal, spheresout] = spotslocator_HoG_LoG(imgDec, GFP, DEC, signal, dims, spotsize, resXY, resZ, sigQual, minSigTh)
%rrr
% Spots
% 1: X
% 2: Y
% 3: ZminSigTh
% 4: Volume
% 5: Intensity
% 6: Nucleus ID
% 7: Position (distance to membrane)

    fprintf('3D spots localisation\n');
    imDec = double(imgDec);
    GFP = double(GFP);
    DEC = double(DEC);
    
    clear imgDec
    
    sizeCrop = 3*max(spotsize)+1;
    
    sigLarge = imdilate(signal, strel('disk', sizeCrop));
    
    % Big Muff
    transformed = Anscombe_forward(imDec);
    transformed = medfilt3(transformed, [3, 3, 3]);
    if ~sigQual
        for iz = 1:dims(3)
            transformed(:,:,iz) = wiener2(transformed(:,:,iz), [sizeCrop sizeCrop]);
        end
    end
    transformed = Anscombe_inverse_exact_unbiased(transformed);

    if ~sigQual
        transformed = imtophat(transformed, offsetstrel('ball', ceil(max(spotsize)*sqrt(3)), ceil(max(spotsize)*sqrt(3)*resZ/resXY)));
    end

	transformed = (transformed - min(transformed(:)))./(max(transformed(:))-min(transformed(:)));

    % Frangi filter - Blob detection
	options3d = struct('HoGScaleRange', [min(spotsize)/sqrt(3), max(spotsize)],...
                   'HoGScaleRatio', sqrt(3)/3, ...
                   'HoGAlpha', 0.5, ...
                   'HoGBeta', 0.5, ...
                   'HoGC', 500, ...
                   'BlackWhite', false, ...
                   'verbose', true, ...
                   'Output', 'max');
               
    [J] = HoG_blobs3D(transformed, options3d);
    J = (J-min(J(:)))./(max(J(:))-min(J(:)));
    J(~sigLarge) = min(J(:));
    
    % MSER CV toolbox
    r2dareas = [floor(spotsize(1)), ceil(pi*3*spotsize(end))];
    regions = cell([dims(3), 1]);
    cc = cell([dims(3), 1]);
    for ii = 1:dims(3)
        [regions{ii},cc{ii}] = detectMSERFeatures(J(:,:,ii), 'ThresholdDelta', 0.8, 'RegionAreaRange', r2dareas, 'MaxAreaVariation', 1); %.5
    end
    M = false(dims);
    for iz = 1:dims(3)
        Mz = false(dims(1:2));
        pl = cat(2, cc{iz}.PixelIdxList);
        Mz(cell2mat(pl')) = true;
        M(:,:,iz) = Mz;
    end
    
    % Check if local maxima
    maxloc = imregionalmax(imhmax(J, (1+minSigTh).*std(J(:))), 26);
    M = imreconstruct(maxloc, M);
    
    clear maxloc
    
    % Add connected spots to nuclei segmentation (for good measure)
    switch class(signal)
        case 'uint8'
            M = im2uint8(M);
        case 'uint16'
            M = im2uint16(M);
        case 'double'
            M = max(signal(:)).*im2double(M);
            
    end
    
    signal = imreconstruct(signal, M+signal);
    M(signal == 0) = 0;
    
    % Watershed the spots since matlab mser tend to return non convex
    % regions [mod 08/11/19]
    ws = watershed(-J);
    M(ws == 0) = 0;
    
    clear dmap ws;
    
    % Get non bg ones (areas and HoG response) - mod [190924] th per nuc
    signalbg = signal;
    signalbg(M>0) = 0;
    CN0 = bwconncomp(signalbg);
    CN = CN0;
    CN.NumObjects = max(signalbg(:));
    CN.PixelIdxList = cell([1 CN.NumObjects]);
    for i = 1:CN.NumObjects
        CN.PixelIdxList{i} = find(signalbg(:) == i); % [mod 22/11/19]
    end % [/mod 12/11/19]
    ths = cellfun(@(x) mean(DEC(x)) + minSigTh*std(DEC(x)), CN.PixelIdxList);

    CC = bwconncomp(M, 6);
    
    % Clear ponctual spot detections
    vol = cat(1, cellfun(@length, CC.PixelIdxList));
    tokill = vol<9;
    for icc = find(tokill)
        M(CC.PixelIdxList{icc}) = 0;
    end
    CC.PixelIdxList(tokill) = [];
    CC.NumObjects = CC.NumObjects - sum(tokill);
    
    % Return to spots intensity filtering on dec (to avoid out of focus light)
    idnucs = cellfun(@(x) median(signal(x)), CC.PixelIdxList);
    respI = cellfun(@(x) max(DEC(x)), CC.PixelIdxList);
    
    tokill = respI<ths(idnucs);
    
    for icc = find(tokill)
        M(CC.PixelIdxList{icc}) = 0;
    end
    
    idnucs(tokill) = [];
    CC.PixelIdxList(tokill) = [];
    CC.NumObjects = CC.NumObjects - sum(tokill);
    
    % Filter again on ori image
    signalbg = signal;
    signalbg(M>0) = 0;
    
    CN0 = bwconncomp(signalbg);
    CN = CN0;
    CN.NumObjects = max(signalbg(:));
    CN.PixelIdxList = cell([1 CN.NumObjects]);
    for i = 1:CN.NumObjects
        CN.PixelIdxList{i} = find(signalbg(:) == i); % [mod 22/11/19]
    end % [/mod 12/11/19]
    
    ths = cellfun(@(x) mean(GFP(x)) + minSigTh*std(GFP(x)), CN.PixelIdxList);
    respI = cellfun(@(x) mean(GFP(x)), CC.PixelIdxList);
    
    tokill = respI<ths(idnucs);
    
    for icc = find(tokill)
        M(CC.PixelIdxList{icc}) = 0;
    end
    CC.PixelIdxList(tokill) = [];
    CC.NumObjects = CC.NumObjects - sum(tokill);
    
    % Extract spots information
    spheres = bwlabeln(M, 6);
    spheresout = zeros(size(spheres));
    dmap = bwdist(~signal);
    spotsCC = table2struct(regionprops3(spheres, GFP, 'Centroid', 'Volume', 'VoxelValues', 'Solidity'));%, 'PrincipalAxisLength'));
    
    ske = bwskel(spheres>0);
    ske_lab = ske.*spheres;
    skeCC = table2struct(regionprops3(ske_lab, 'Volume'));
    tau = ceil(4/3*pi*(max(spotsize)/2/sqrt(2))^3);
    elong = cat(1,skeCC.Volume);
    elong = 1-exp(-(elong-1)/tau);
    
    nbPosi = size(spotsCC, 1);
    spots = zeros(nbPosi, 8);
    
    for ispot = 1:nbPosi
        posi = round(spotsCC(ispot).Centroid);
        spots(ispot,1) = posi(2);
        spots(ispot,2) = posi(1);
        spots(ispot,3) = posi(3);
        spots(ispot,4) = spotsCC(ispot).Volume;
        spots(ispot,5) = sum(spotsCC(ispot).VoxelValues);
        spots(ispot,6) = signal(posi(2), posi(1), posi(3));
        spots(ispot,7) = dmap(posi(2), posi(1), posi(3));
        spots(ispot,8) = elong(ispot);
    end
    
    % Remove outsiders
    outsiders = spots(:, 6) == 0;
    spots(outsiders,:) = [];
    
    % Rebuilt spheres map
    nbPosi = size(spots, 1);
    for ispot = 1:nbPosi
        sphereslab = spheres(spots(ispot,1), spots(ispot,2), spots(ispot, 3));
        spheresout(spheres == sphereslab) = ispot;
    end
    fprintf('%d spots detected\n', nbPosi);
end
