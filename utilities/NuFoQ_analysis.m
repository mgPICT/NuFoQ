function [error] = NuFoQ_analysis(analysisFolder, resultsFolder, saveImg, params, onylOne, filename)
% Main NuFoQ function
% Manages all the steps of the script, from denoising to quantification, by
% calling the relevant functions and count the number of files returning
% errors during processing.

%% Initialisation
% Get parameters and prepare folders
if iscell(analysisFolder)
    Loc_data = analysisFolder{1};
    Loc_data2 = analysisFolder{2};
    if onylOne
        filename = filename{1};
    end
    two_sets = true;
else
    Loc_data = analysisFolder;
    two_sets = false;
end
Loc_resu = resultsFolder;

resXY = str2double(params{1});
resZ =  str2double(params{2});
spotsize =  eval(params{3});
toSkip = str2double(params{4});
Sig_info = params{5};
Sig_qual = eval(params{6});
minSigTh = eval(params{7});

if two_sets
    fprintf('Data locations: %s\t%s\nResults location: %s\n\n', Loc_data, Loc_data2, Loc_resu);
else
    fprintf('Data location: %s\nResults location: %s\n\n', Loc_data, Loc_resu);
end
fprintf('\t---\n\tParameters:\n');
fprintf('\tPixel size: %.3fx%.3fx%.3f\n', resXY, resXY, resZ);
fprintf('\tSpot sizes:');
for ii = 1:length(spotsize)
    fprintf(' %.2f', spotsize(ii));
end
fprintf('\n\tMinimum variations (in stds): %.2f\n', minSigTh);
fprintf('\tSignal tag: %s\n', Sig_info);
fprintf('\tLow signal?: %d\n\t---\n', Sig_qual);


if(~isfolder(Loc_resu))
	mkdir(Loc_resu);
end

% Set inter parameters and get files
error = 0;

cd(Loc_data);

files = dir;
files = files(3:end);

nbFiles = length(files);
tokill = zeros(nbFiles, 1);
for i_f = nbFiles:-1:1
    if files(i_f).isdir
        tokill(i_f) = 1;
    else
        try
            imfinfo(files(i_f).name);
        catch err
            fprintf('File %s: skipped !\n%s\n', files(i_f).name, err.message);
            tokill(i_f) = 1;
        end
    end
end
files(tokill == 1) = [];
nbFiles = length(files);

if two_sets
    cd(Loc_data2);
    files2 = dir;
    files2 = files2(3:end);
    
    nbFiles2 = length(files2);
    tokill = zeros(nbFiles2, 1);
    for i_f = nbFiles2:-1:1
        if files2(i_f).isdir
            tokill(i_f) = 1;
        else
            try
                imfinfo(files2(i_f).name);
            catch err
                fprintf('File %s: skipped !\n%s\n', files2(i_f).name, err.message);
                tokill(i_f) = 1;
            end
        end
    end
    files2(tokill == 1) = [];
    nbFiles2 = length(files2);
    
    if nbFiles ~= nbFiles2
        warning('Warnings:nbFilesnbFiles2', 'Numbers of files in both directory are not the same! Images with missing correspondance will be skipped!');
    end
end

%% Associate Segmentation help images to quantification ones
gfpimgs = zeros(nbFiles, 1);
for i_f = 1:nbFiles
    % Check for gfp
    if strfind(files(i_f).name, Sig_info)
        gfpimgs(i_f) = 1;
    elseif strfind(files(i_f).name, upper(Sig_info))
        gfpimgs(i_f) = 1;
    end
end
nbGFP = sum(gfpimgs);
gfpimgs = find(gfpimgs);
if two_sets
    decimgs = zeros(nbFiles2, 1);
    for i_f = 1:nbFiles2
        % Check for dec
        if strfind(files2(i_f).name, Sig_info)
            decimgs(i_f) = 1;
        elseif strfind(files2(i_f).name, upper(Sig_info))
            decimgs(i_f) = 1;
        end
    end
    nbDec = sum(decimgs);
    if nbDec ~= nbGFP
        warning('Warnings:DecGfpNb', 'Numbers of quantification and segmentation images are not the same! Images with missing correspondance will be skipped!');
    end
    decimgs = find(decimgs);
end

%% Analysis loop on every files
for i_file = 1:nbGFP
    tic;
    cd(Loc_data);
    i_gfp = gfpimgs(i_file);
    if strcmp(files(i_gfp).name(end-3:end), '.stk')
        continue
    end
    
    if onylOne
        if ~strcmp(files(i_gfp).name, filename)
            continue
        end
    end
    
    if two_sets
        dotpos = strfind(files(i_gfp).name, '.');
        basename = files(i_gfp).name(1:dotpos(end));
        % Get dec image
        for i_file2 = 1:nbDec
            if strfind(files2(i_file2).name, basename)
                i_dec = decimgs(i_file2);
                break;
            end
        end
    end

    if ~two_sets
        fprintf('\nAnalysing files:\t%s\t:\t%s\n', Sig_info, files(i_gfp).name);
    else
        fprintf('\nAnalysing files:\t%s\t:\t%s\n\t\t\tSpots\t:\t%s\n', Sig_info, files(i_gfp).name, files2(i_dec).name);
    end
    
    %% File reading
    fprintf('\t> Reading data information and preparing output files. \n');

    info = imfinfo(files(i_gfp).name);
    dims = uint16([info(1).Height info(1).Width length(info)]);
    
    switch info(1).BitDepth
        case 16
            type = 'uint16';
        case 8
            type = 'uint8';
        otherwise
            type = 'double';
    end
    
    dims(3) = dims(3) - toSkip;
    GFP = zeros(dims(1), dims(2), dims(3), type);
    if two_sets
        DEC = zeros(dims(1), dims(2), dims(3), type);
    end
    for iz = (1+toSkip):dims(3)
        GFP(:,:,iz) = imread(files(i_gfp).name, iz);
        if two_sets
            DEC(:,:,iz) = imread([Loc_data2, filesep, files2(i_dec).name], iz);
        end
    end
    
    %% STEP 0: Correct for background variations and denoise
    transformed = double(GFP);
    bckg_med = zeros(dims(3), 1);
    for i_z = 1:dims(3)
        slice = transformed(:,:,i_z);
        bckg_med(i_z) = median(slice(:));
    end

    for i_z = 1:dims(3)
        transformed(:,:,i_z) = transformed(:,:,i_z)./bckg_med(i_z);
    end

    transformed = Anscombe_forward(transformed);
    wienersize = max(min(spotsize), 3)*2;
    for i_z = 1:dims(3)
        transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
    end
    GFPmod = Anscombe_inverse_exact_unbiased(transformed);
    
    if two_sets
        transformed = double(DEC);
        bckg_med = zeros(dims(3), 1);
        for i_z = 1:dims(3)
            slice = transformed(:,:,i_z);
            bckg_med(i_z) = median(slice(:));
        end

        for i_z = 1:dims(3)
            transformed(:,:,i_z) = transformed(:,:,i_z)./bckg_med(i_z);
        end

        transformed = Anscombe_forward(transformed);
        wienersize = max(min(spotsize), 3)*2;
        for i_z = 1:dims(3)
            transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
        end
        DECmod = Anscombe_inverse_exact_unbiased(transformed);
    end
    
    %% STEP 1: NUCLEI SEGMENTATION (in 3D)
    try
        [seg_gfp] = segment3Dfcm(GFPmod, dims, max(spotsize)+1, resXY, resZ, Sig_qual);
    catch err
        fprintf(err.message);
        error = error+1;
        continue
    end

    if max(seg_gfp(:)) == 0
        fprintf('No Nuclei have been found\n');
        outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
        fid = fopen(outputfile, 'w+');
        fprintf(fid, 'No Nuclei were found');
        fclose(fid);
        continue
    end
    
    %% STEP 2: SPOTS SEGMENTATION (in 3D)
    try
        if two_sets
            [spots, seg_gfp, seg_spots] = spotslocator_HoG_LoG(DECmod, GFP, DEC, seg_gfp, dims, spotsize, resXY, resZ, Sig_qual, minSigTh);
        else
            [spots, seg_gfp, seg_spots] = spotslocator_HoG_LoG(GFPmod, GFP, GFP, seg_gfp, dims, spotsize, resXY, resZ, Sig_qual, minSigTh);
        end

    catch err
        fprintf(err.message);
        error = error+1;
        continue
    end
    if size(spots,1) == 0
        fprintf('No Foci have been found\n');
        outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
        fid = fopen(outputfile, 'w+');
        fprintf(fid, 'No Foci have been found');
        fclose(fid);
        toc;
        continue
    end
    

    nbSpots = size(spots,1);
    if saveImg
        seg_name = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '_segSpots.tif');
        imwrite(uint16(seg_spots(:,:,1)), seg_name, 'writemode', 'overwrite');
        for i = 2:dims(3)
            imwrite(uint16(seg_spots(:,:,i)), seg_name, 'writemode', 'append');
        end
        seg_name = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '_segNuc.tif');
        imwrite(uint16(seg_gfp(:,:,1)), seg_name, 'writemode', 'overwrite');
        for i = 2:dims(3)
            imwrite(uint16(seg_gfp(:,:,i)), seg_name, 'writemode', 'append');
        end
    end

    %% STEP 3: NUCLEI STATISTCS (in 3D)
    try
        [nuclei] = measureInfosNuclei(GFP, GFPmod, seg_gfp, seg_spots, spots, saveImg, strcat(Loc_resu, '/',files(i_file).name(1:end-4), '_segSpotsSpheres.tif'));
    catch err
        fprintf(err.message);
        error = error+1;
        continue
    end

    %% WRITE OUTPUT  NUCLEI % --------------------------------------------------------- %
	% Nucleus: 1> Intensity					Output: 1> Nucleus Volume                   %
	% Nucleus: 2> Intensity background				2> Nucleus intensity                %
	% Nucleus: 3> Number of foci					3> Nucleus intensity background     %
	% Nucleus: 4> Volume							4> Nucleus number of spots      	%
	%                       						5> ToKeep                       	%
	% --------------------------------------------------------------------------------- %
	
    outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
    fid = fopen(outputfile, 'w+');
    fprintf(fid, 'Nucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\n');

    bckg_int = median(double(GFP(:)));
	nbNucs = size(nuclei,1);
    
    tokeepnuc = nuclei(:,4) > 999; % Hard-coded to avoir smaller chunk of signal (dust and whatnot)
    
    for i = 1:nbNucs
		%              I Ibg   V  nbSpots  keep
        fprintf(fid, '%d\t%d\t%.0f\t%.0f\t%.0f', ...
            nuclei(i, 1)-bckg_int*nuclei(i, 4), nuclei(i, 2)-bckg_int*nuclei(i, 4), nuclei(i, 4), nuclei(i, 3), tokeepnuc(i));
        if i < nbNucs
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    %% WRITE OUTPUT  SPOTS % --------------------------------------------------------------------------------------------------------------------------------- %
	% Spots: 1> X						Nucleus: 1> Intensity					Output: 1> X spot							8> Relative distance to border		%
	%		 2> Y						Nucleus: 2> Intensity background				2> Y spot							9> Nucleus ID						%
	%		 3> Z						Nucleus: 3> Number of foci						3> Z spot							A> Nucleus intensity				%
	%		 4> Volume					Nucleus: 4> Volume								4> Volume spot						B> Nucleus intensity background		%
	%		 5> Intensity																5> Intensity spot					C> Nucleus Volume					%
	%		 6> Nucleus ID																6> Relative spot intensity			D> Nucleus number of spots			%
	%		 7> distance to border		8> Sphericity (if low)							7> Distance to nucleus border		E> ToKeep		F> Sphericity		%
	% --------------------------------------------------------------------------------------------------------------------------------------------------------- %
	
    outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.spots');
    fid = fopen(outputfile, 'w+');
    fprintf(fid, 'X\tY\tZ\tVolume\tIntensity\tRelative intensity\tDistance to border\tRelative distance to border\tNucleus ID\tNucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\n');

    bckg_int = median(double(GFP(:)));
	nuclei_radii = (3/(4*pi)*nuclei(:,4)).^(1/3);
    
    for i = 1:nbSpots
        relInt = (spots(i,5) -bckg_int*spots(i,4))/(nuclei(spots(i,6), 1) - nuclei(spots(i,6), 4) * bckg_int);
		%              X     Y     Z     V      I     Ir    D     Dr   nuc  I Ibg   V  nbSpots  keep   sphericity
        fprintf(fid, '%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.4f\t%.0f\t%.2f\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f\t%d', ...
            spots(i,1), spots(i,2), spots(i,3), spots(i,4), spots(i,5)-spots(i,4)*bckg_int, relInt, spots(i,7), spots(i,7)/nuclei_radii(spots(i,6)),...
            spots(i,6), nuclei(spots(i,6), 1)-bckg_int*nuclei(spots(i,6), 4), nuclei(spots(i,6), 2)-bckg_int*nuclei(spots(i,6), 4), nuclei(spots(i,6), 4), nuclei(spots(i,6), 3), tokeepnuc(spots(i,6)), spots(i,8));
        if i < nbSpots
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    toc;
end
fprintf('\nqFOCI analysis finished !\n-----------------------------\n');
