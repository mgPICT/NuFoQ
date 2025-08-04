function supervised_filtering_results(analysisFolder, resultsFolder,save_imgs, ~, onlyOne, filename, loc)
%% supervised_filtering_results
%  Read NoFoQ .nucs files along with images for a supervised filtering of
%  the results
% 
%% Output
%  Figures designed with Antoine Hocher on R displaying:
%      > the histogram of the number of foci per cell
%      > the pdfs of the percentage of intensity contained in the foci
%      > the pdfs of the total intensity contained in the foci
%      > the pdfs of the foci position in the cell (distance to border)
%  for each cell line (WT and mutants)
%
%  Bonus figures:
%      > Boxplots per yeastLines of intensities (absolutes and relatives)
%      > Pdfs of relatives and absolutes intensities of spots
%
% Micka�l Garnier - UMR3664 - PICT@PASTEUR - 2018 - 2025

% Initiliaze colors
Ckeep = [.25 .5 .25];
Ckill = [.5 .25 .25];
Ccurr = [1 1 0];
CcurrSpots = [1 0 1];

% Get parameters and prepare folders
if iscell(analysisFolder)
    Loc_data = analysisFolder{1}; % Display original image
else
    Loc_data = analysisFolder;
end
Loc_resu = resultsFolder;

if(~isfolder(Loc_resu))
	error('The selected result folder does not exist, you need to run the analysis first !');
end

fprintf('Initializing interface data\n');

% Manually select the files to analyse
if onlyOne
    baseName = filename(1:strfind(filename, '.'));
    nucsFiles = [baseName, 'spots'];
    nucsPath = Loc_resu;
else
    [nucsFiles, nucsPath] = uigetfile([Loc_resu, filesep, '*.spots'], 'Select the .spots files to check', 'MultiSelect', 'off');
end
nbNucsFiles = size(nucsFiles, 2);
if ~iscell(nucsFiles)
    nucsFiles = {nucsFiles};
    nbNucsFiles = 1;
end

% Retrieve corresponding images
imgList = dir(Loc_data);
nbImgs = size(imgList, 1);
imgFiles = cell(size(nucsFiles));
transFiles = cell(size(nucsFiles));
for inucs = 1:nbNucsFiles
    baseName = nucsFiles{inucs}(1:strfind(nucsFiles{inucs}, '.spots'));
    
    for iim = 1:nbImgs
        if strfind(imgList(iim).name, baseName) == 1
            % Check if trans or signal
            if strfind(imgList(iim).name, 'Trans')
                transFiles{inucs} = imgList(iim).name;
            elseif strfind(imgList(iim).name, 'TRANS')
                 transFiles{inucs} = imgList(iim).name;   
            else
                imgFiles{inucs} = imgList(iim).name;
            end
        end
    end
end

%% Read Data
for inucs = 1:nbNucsFiles
    fprintf('Supervising: %s\n', imgFiles{inucs});
    
    % Read image
    fprintf('Reading images\n');
    imginfo = imfinfo([Loc_data, filesep, imgFiles{inucs}]);
    dims = [imginfo(1).Height imginfo(1).Width length(imginfo)];
    imageRatio = dims(1)/dims(2);
    
    img = zeros(dims);
    for iz = 1:dims(3)
        img(:,:,iz) = imread([Loc_data, filesep, imgFiles{inucs}], iz);
    end
    
    % Look for trans image (for multiple display)
    if ~isempty(transFiles{inucs})
        trans = imread([Loc_data, filesep, transFiles{inucs}]);
        figure, imshow(trans,[]);
    end
    
    % Pre-compute projections
    fprintf('Pre-computing projections\n');
    maxProj = max(img, [], 3);
    sumProj = sum(img, 3);
    stdProj = std(img, 1, 3);
    meanProj = mean(img, 3);
    medProj = median(img, 3);
    minProj = min(img, [], 3);

    minInt = min(img(:));
    maxInt = max(img(:));

    % Read analysis outpus
    fprintf('Loading analysis data\n');
    [X,Y,Z,Volume,~,~,~,~,NucleusID,~,~,NucleusVolume,~,Tokeep] = importSpots([nucsPath, filesep, nucsFiles{inucs}]);
    nucsfile = replace(nucsFiles{inucs}, '.spots', '.nucs');
    [~,~,NucleusVolume_nucs,~,ToKeepNucs] = importNucs([nucsPath, filesep, nucsfile]);
    
    % Check files integrity !
    for ii = 1:size(NucleusID, 1)
        if NucleusVolume(ii) == NucleusVolume_nucs(NucleusID(ii))
            continue
        else
            error("Missmatch between nucleus intensities in *.spots and *.nucs files !");
        end
    end
    
    % Get nucleus lists
    nbNucs = length(ToKeepNucs);
    nucList = (1:nbNucs)';
    
%     [nucList2, ~, ispots2nucs] = unique(NucleusID);
    spotsPerNucs = histc(NucleusID, 1:max(nucList));
    if sum(ToKeepNucs==0) > 0
        listNucsKept = mat2cell(num2str(nucList), ones(length(nucList), 1))';
        
        ids = find(ToKeepNucs == 0);
        listNucsKill = listNucsKept(ids);
        listNucsKept(ids) = [];
    else
        listNucsKept = mat2cell(num2str(nucList), ones(length(nucList), 1))';
        listNucsKill = {};
    end
    
    listNucsKill = cellfun(@strtrim, listNucsKill, 'UniformOutput', false);
    listNucsKept = cellfun(@strtrim, listNucsKept, 'UniformOutput', false);
    
    % Get spot lists
    nbSpots = size(X,1);
    if sum(Tokeep) > 0
        listSpots = find(Tokeep);
        listSpotsKept = num2cell(listSpots, length(listSpots))';
        listSpotsKept = cellfun(@num2str, listSpotsKept, 'UniformOutput', false);
    else 
        listSpotsKept = {};
    end
    if sum(~Tokeep) > 0
        listSpots = find(~Tokeep);
        listSpotsKill = num2cell(listSpots, length(listSpots))';
        listSpotsKill = cellfun(@num2str, listSpotsKill, 'UniformOutput', false);
     else 
        listSpotsKill = {};
    end
    
    
    % Prepare nuclei map
    fprintf('Pre-rendering contour maps\n');
    C_nucs = cell([1 max(nucList)]);
    if (save_imgs)
        spheres_nucs = zeros(dims);
        for iz = 1:dims(3)
            spheres_nucs(:,:,iz) = imread([nucsPath, filesep, nucsFiles{inucs}(1:end-6), '_segNuc.tif'], iz);
        end
    else
        disp('No images saved, generating spheres for nuclei and spots based on positions and volume');
        spheres_nucs = zeros(dims);
        for ispot = 1:nbSpots
            nucid = NucleusID(ispot);
            radii = ceil((NucleusVolume(ispot)*3/4/pi)^(1/3));
           % Create 3D neighborhood
            voi_R = round(X(ispot)) + (-ceil(radii) : ceil(radii));
            voi_C = round(Y(ispot)) + (-ceil(radii) : ceil(radii));
            voi_Z = round(Z(ispot)) + (-ceil(radii) : ceil(radii));

            voi_R = voi_R(voi_R<=dims(1));
            voi_R = voi_R(voi_R>0);
            voi_C = voi_C(voi_C<=dims(2));
            voi_C = voi_C(voi_C>0);
            voi_Z = voi_Z(voi_Z<=dims(3));
            voi_Z = voi_Z(voi_Z>0);


            % Create local sphere
            [YY, XX, ZZ] = meshgrid(voi_C, voi_R, voi_Z);

            sphere = (XX - X(ispot)).^2 + ...
                     (YY - Y(ispot)).^2 + ...
                     (ZZ - Z(ispot)).^2 <= radii^2;

            spheres_nucs(voi_R, voi_C, voi_Z) = max(cat(4,spheres_nucs(voi_R, voi_C, voi_Z), double(sphere).*nucid), [], 4);
        end
    end
        
    % Contour map (for app) and contour image (for volumeViewer)
    tmp = max(spheres_nucs, [], 3);

    for ival = nucList'
        if sum(tmp(:) == ival) == 0  % [mod 12/11/19]
            continue
        end                          % [/mod 12/11/19]
        structtmp = contourdata(contourc(double(tmp == ival), 1));
        C_nucs{ival}{1} = cat(1,structtmp.xdata);
        C_nucs{ival}{2} = cat(1,structtmp.ydata);
    end

%     nucs_edges = zeros(dims);
%     for iz = 1:dims(3)
%         nucs_edges(:,:,iz) = edge(spheres_nucs(:,:,iz), 'canny').*spheres_nucs(:,:,iz);
%     end
    
    
    % Prepare spots map
    C_spots = cell([1 max(nucList)]);
    if (save_imgs)
        spheres_spots = zeros(dims);
        for iz = 1:dims(3)
            spheres_spots(:,:,iz) = imread([nucsPath, filesep, nucsFiles{inucs}(1:end-6), '_segSpots.tif'], iz);
        end
    else
        spheres_spots = zeros(dims);
        for ispot = 1:nbSpots
            radii = ceil((Volume(ispot)*3/4/pi)^(1/3));
           % Create 3D neighborhood
            voi_R = round(X(ispot)) + (-ceil(radii) : ceil(radii));
            voi_C = round(Y(ispot)) + (-ceil(radii) : ceil(radii));
            voi_Z = round(Z(ispot)) + (-ceil(radii) : ceil(radii));

            voi_R = voi_R(voi_R<=dims(1));
            voi_R = voi_R(voi_R>0);
            voi_C = voi_C(voi_C<=dims(2));
            voi_C = voi_C(voi_C>0);
            voi_Z = voi_Z(voi_Z<=dims(3));
            voi_Z = voi_Z(voi_Z>0);


            % Create local sphere
            [YY, XX, ZZ] = meshgrid(voi_C, voi_R, voi_Z);

            sphere = (XX - X(ispot)).^2 + ...
                     (YY - Y(ispot)).^2 + ...
                     (ZZ - Z(ispot)).^2 <= radii^2;

            spheres_spots(voi_R, voi_C, voi_Z) = max(cat(4,spheres_spots(voi_R, voi_C, voi_Z), double(sphere).*ispot), [], 4);
        end
    end
    
    % Contour map (for app) and contour image (for volumeViewer)
    tmp = max(spheres_spots, [], 3);

    spots2fill = ones(1, length(spotsPerNucs));
    for inuc = nucList'
        C_spots{inuc} = cell([spotsPerNucs(inuc) 3]);
    end
    for ispot = 1:nbSpots
        inuc = NucleusID(ispot);
        
        map = contourc(double(tmp == ispot), 1);
        if isempty(map)
            map = contourc(double(max(spheres_spots == ispot, [], 3)), 1);
        end
        structtmp = contourdata(map);
        C_spots{inuc}{spots2fill(inuc), 1} = structtmp.xdata;
        C_spots{inuc}{spots2fill(inuc), 2} = structtmp.ydata;
        C_spots{inuc}{spots2fill(inuc), 3} = ispot;
        spots2fill(inuc) = spots2fill(inuc)+1;
    end

%     spots_edges = zeros(dims);
%     for iz = 1:dims(3)
%         spots_edges(:,:,iz) = edge(spheres_spots(:,:,iz), 'canny').*spheres_spots(:,:,iz);
%     end
    
    % Create volumeViewer image and launch viewer
%     volume = img./max(img(:));
%     volume(nucs_edges>0) = 1;
%     volume(spots_edges>0) = 1;
%     volumeViewer(volume);

    %% Load interface
    fprintf('Loading interface\n');
    % Basic panel
    screenSize = get(0,'ScreenSize');
    figureHeight = screenSize(4)*.95;
%     figureLength = min(screenSize(3:4))*.99;
    figureWidth =  screenSize(3)*.95; %max(screenSize(3), dims(2)+600)*.95;
    figureSize = [1, 1, figureWidth, figureHeight];%(figureLength/imageRatio)+600, figureLength];
    filterFig = uifigure('Visible','on','Position', figureSize,...
    'Color',[0.85 0.85 0.87], 'Resize','on',...
        'Name', 'qFOCI 2.0 - Institut Curie - Supervised filtering of analysis outputs',...  % Title figure
        'NumberTitle', 'off',... % Do not show figure number
        'Tag', [nucsPath, filesep, nucsFiles{inucs}], ...
        'MenuBar', 'none');
    movegui(filterFig, 'center');

    % Image display
    axesPos = zeros(1,4)+25;
    %axesPos(2) = 200;
    axesPos(3) = figureSize(3)-675;%(figureSize(4)-75)/imageRatio; %(figureSize(4)-550)*imageRatio;
    axesPos(4) = figureSize(4)-75;%figureSize(4)-550;

    image_panel=uiaxes(filterFig,...
        'Visible', 'on', ...
        'Position', axesPos,...
        'XLim', [1 dims(2)], 'YLim', [1 dims(1)], ...
        'XTickLabelMode', 'manual', 'YTickLabelMode', 'manual', ...
        'XTickMode', 'manual', 'YTickMode', 'manual', ...
        'DataAspectRatioMode', 'manual', ...
        'PlotBoxAspectRatioMode', 'manual', ...
        'BackgroundColor',[.15 .15 .15]);
    
    H_nucs = cell(size(C_nucs));
    H_spots = cell(size(C_nucs));
    
    H_map = imagesc(image_panel, maxProj);
    H_map.CDataMapping = 'scaled';
    colormap(image_panel, gray);
    
    % Add overlayed spots and nucleus and add to plot handles !
    % Kept
    for ic = 1:length(listNucsKept)
        inuc = str2double(listNucsKept{ic});
        hold(image_panel, 'on');
        H_nucs{inuc} = plot(image_panel, C_nucs{inuc}{1}, C_nucs{inuc}{2}, 'Color', Ckeep, 'LineWidth', 2);
        H_spots{inuc} = cell([2 spotsPerNucs(inuc)]);
        for ics = 1:spotsPerNucs(inuc)
            hold(image_panel, 'on');
            H_spots{inuc}{1, ics} = plot(image_panel, C_spots{inuc}{ics, 1}, C_spots{inuc}{ics, 2}, 'Color', min([1, 1, 1], Ckeep.*2), 'LineWidth', 1);
            H_spots{inuc}{2, ics} = C_spots{inuc}{ics, 3};
        end
    end
    % Killed
    for ic = 1:length(listNucsKill)
        inuc = str2double(listNucsKill{ic});
        hold(image_panel, 'on');
        H_nucs{inuc} = plot(image_panel, C_nucs{inuc}{1}, C_nucs{inuc}{2}, 'Color', Ckill, 'LineWidth', 2);
        H_spots{inuc} = cell([2 spotsPerNucs(inuc)]);
        for ics = 1:spotsPerNucs(inuc)
            hold(image_panel, 'on');
            H_spots{inuc}{1, ics} = plot(image_panel, C_spots{inuc}{ics, 1}, C_spots{inuc}{ics, 2}, 'Color', min([1, 1, 1], Ckill.*2), 'LineWidth', 1);
            H_spots{inuc}{2, ics} = C_spots{inuc}{ics, 3};
        end
    end
    
    % Correct killed spots color
    if ~isempty(listSpotsKill)
        nucswithspots = find(cellfun(@length, H_spots));
        spotspernucs = cellfun(@(x) size(x,2), H_spots(nucswithspots));
        nucswithspotspernucs = repelem(nucswithspots, spotspernucs);
        spotsorder = cat(2, H_spots{:});
        [~, spotsorder] = sort(cat(2,spotsorder{2,:}));
        spots2nucs = nucswithspotspernucs(spotsorder);
        for ic = 1:length(listSpotsKill)
            ispot = str2double(listSpotsKill(ic));
            inuc = spots2nucs(ispot);
            ics = sum(spots2nucs(1:ispot)==inuc);
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
        end
    end
    

    % Buttons - display
    display_panel = uipanel(filterFig, 'Title','Display settings',...
        'Position', [figureSize(3)-625, figureSize(4)-500, 600, 450], ...%[25, 25, (figureSize(4)-550)*imageRatio, 100],...
        'FontWeight', 'bold', ...
        'BackgroundColor',[0.85 0.85 0.87],...
        'ForegroundColor',[0 0 0]);

    uilabel(display_panel, ...
                     'Position', [400, 50, 150, 25], ... % 325 50 225 25
                     'HorizontalAlignment', 'center', ...
                     'Text', 'Colormap');
    uidropdown(display_panel, ...
                        'Items', {'gray', 'red', 'blue', 'green', 'parula', 'jet', 'hsv', 'hot', 'cool', 'sping', 'summer', 'autumn', 'winter', 'bone', 'copper', 'pink'}, ...
                        'Value', 'gray', ...
                        'Position', [400, 25, 150, 25] ,...
                        'ValueChangedFcn', @(cmap_choice,event) updateCmap(cmap_choice,image_panel));
                    
    uilabel(display_panel, ...
                     'Position', [210, 50, 150, 25], ...
                     'HorizontalAlignment', 'center', ...
                     'Text', 'Overlay');
	uiswitch(display_panel, 'slider', ...
                     'Value', 'On', ...
                     'Position', [260, 25, 100, 25], ...
                     'ValueChangedFcn', @(overlay_switch, event) updateOverlay(overlay_switch, H_spots, H_nucs));

    uilabel(display_panel, ...
                     'Position', [25, 50, round(((figureSize(4)-550)*imageRatio)/4), 25], ...
                     'HorizontalAlignment', 'center', ...
                     'Text', 'Contrast');

    low_slider = uislider(display_panel);
    high_slider = uislider(display_panel);
        low_slider.Limits = [0, maxInt*2];
        low_slider.Value = minInt;
        low_slider.MajorTicks = []; low_slider.MinorTicks = [];
        low_slider.Position = [60, 40, min(round(((figureSize(4)-550)*imageRatio)/4), 150), 3]; %[25, 40, round(((figureSize(4)-550)*imageRatio)/4), 3];
        low_slider.ValueChangedFcn = @(low_slider,event) updateInt(low_slider,high_slider,image_panel);
        high_slider.Limits = [0, maxInt*2];
        high_slider.Value = maxInt;
        high_slider.MajorTicks = []; high_slider.MinorTicks = []; 
        high_slider.Position = [60, 20, min(round(((figureSize(4)-550)*imageRatio)/4), 150), 3]; %[25, 20, round(((figureSize(4)-550)*imageRatio)/4), 3];
        high_slider.ValueChangedFcn = @(high_slider,event) updateInt(low_slider,high_slider,image_panel);
    uilabel(display_panel, ...
                     'Position', [25, 40, 30, 15], ...
                     'HorizontalAlignment', 'center', ...
                     'Text', 'Min');
    uilabel(display_panel, ...
                     'Position', [25, 20, 30, 15], ...
                     'HorizontalAlignment', 'center', ...
                     'Text', 'Max');
        
    uiknob(display_panel, 'discrete', ...
                   'Items', {'Std deviation', 'Maximum', 'Minimum', 'Sum', 'Mean', 'Median'}, ...
                   'Value', 'Maximum', ...
                   'Position', [225, 200, 200, 150],... 
                   'ValueChangedFcn',@(proj_knob,event) updateProj(proj_knob, image_panel, H_map, stdProj, maxProj, sumProj, meanProj, medProj, minProj, low_slider, high_slider));

    % Buttons - data
    dps = figureSize(4)-550;
    data_panel = uipanel(filterFig,'Title','Data',...
        'Position', [figureSize(3)-625, 25, 600, dps], ...%[(figureSize(4)-550)*imageRatio+50, 25, figureLength+400-(figureSize(4)-550)*imageRatio-75, figureSize(4)-425],...
        'FontWeight', 'bold', ...
        'BackgroundColor',[0.85 0.85 0.87],...
        'ForegroundColor',[0 0 0]);

    uilabel(data_panel, ...
                     'Position', [25, 5*dps/6+50, 250, 50], ...
                     'HorizontalAlignment', 'left', ...
                     'Text', 'Nuclei kept', ...
                     'FontWeight', 'bold');
                     
    uilabel(data_panel, ...
                     'Position', [325, 5*dps/6+50, 250, 50], ...
                     'HorizontalAlignment', 'left', ...
                     'Text', 'Nuclei killed', ...
                     'FontWeight', 'bold');
                 
	nucs_box_kept = uilistbox(data_panel);
	nucs_box_killed = uilistbox(data_panel);
    spots_box_kept = uilistbox(data_panel);
    spots_box_killed = uilistbox(data_panel);
        nucs_box_kept.Value = {};
        nucs_box_kept.Position = [25, dps/2+50, 250, dps/3];
        nucs_box_kept.Multiselect = 'on';
        nucs_box_kept.Tag = 'kept';
        nucs_box_kept.Items = listNucsKept;
        nucs_box_kept.FontColor = Ckeep.*2;
        nucs_box_kept.BackgroundColor = [.15 .15 .15];
        nucs_box_killed.Value = {};
        nucs_box_killed.Position = [300, dps/2+50, 250, dps/3];
        nucs_box_killed.Multiselect = 'on';
        nucs_box_killed.Tag = 'killed';
        nucs_box_killed.Items = listNucsKill;
        nucs_box_killed.FontColor = Ckill*2;
        nucs_box_killed.BackgroundColor = [.15 .15 .15];        
        spots_box_kept.Value = {};
        spots_box_kept.Position = [25, 125, 250, dps/3];
        spots_box_kept.Multiselect = 'on';
        spots_box_kept.Enable = 'off';
        spots_box_kept.Tag = 'kept';
        spots_box_kept.Items = listSpotsKept;
        spots_box_kept.FontColor = Ckeep.*2;
        spots_box_kept.BackgroundColor = [.15 .15 .15];
        spots_box_killed.Value = {};
        spots_box_killed.Position = [300, 125, 250, dps/3];
        spots_box_killed.Multiselect = 'on';
        spots_box_killed.Enable = 'off';
        spots_box_killed.Tag = 'killed';
        spots_box_killed.Items = listSpotsKill;
        spots_box_killed.FontColor = Ckill*2;
        spots_box_killed.BackgroundColor = [.15 .15 .15];
        nucs_box_kept.ValueChangedFcn = @(nucs_box_kept,event) updateOtherList(nucs_box_killed, nucs_box_kept, image_panel, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr);
        nucs_box_killed.ValueChangedFcn = @(nucs_box_killed,event) updateOtherList(nucs_box_kept, nucs_box_killed, image_panel, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr);
        spots_box_kept.ValueChangedFcn = @(spots_box_kept,event) updateOtherListSpots(spots_box_killed, spots_box_kept, image_panel, H_nucs, H_spots, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots);
        spots_box_killed.ValueChangedFcn = @(spots_box_killed,event) updateOtherListSpots(spots_box_kept, spots_box_killed, image_panel, H_nucs, H_spots, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots);
        
	uilabel(data_panel, ...
                     'Position', [25, dps/3+125, 250, 50], ...
                     'HorizontalAlignment', 'left', ...
                     'Text', 'Spots kept', ...
                     'FontWeight', 'bold');
                     
    uilabel(data_panel, ...
                     'Position', [325, dps/3+125, 250, 50], ...
                     'HorizontalAlignment', 'left', ...
                     'Text', 'Spots killed', ...
                     'FontWeight', 'bold');

    nucsOrSpots_button = uibutton(data_panel, 'state', ...
                          'Text', 'Nuclei', ...
                          'FontWeight', 'bold', ...
                          'Position', [25, 25, 75, 75], ...
                          'ValueChangedFcn', @(nucsOrSpots_button, event) nucsOrSpots(nucsOrSpots_button, image_panel, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, H_nucs, H_spots, Ckeep, Ckill));
        
    uibutton(data_panel, 'push', ...
                          'Text', 'Keep & Kill', ...
                          'FontWeight', 'bold', ...
                          'Position', [125, 25, 75, 75], ...
                          'ButtonPushedFcn', @(killswitch, event) updateNucsOrSpots(nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, killswitch, image_panel, nucsOrSpots_button, H_nucs, H_spots));
                      
	uibutton(data_panel, 'push', ...
                         'Position', [475, 25, 75, 75] , ...
                         'Text', 'Next / End', 'FontWeight', 'Bold', ...
                         'ButtonPushedFcn', @(ok_button, event) fcnOK(spots_box_killed, spots_box_kept, nucs_box_killed, nucs_box_kept, filterFig, inucs, nbNucsFiles));
                     
	[imageLogo, ~, alpha] = imread([loc, filesep, 'institut_curie.png']);
    logdims = size(imageLogo);
    logratio = logdims(1)/logdims(2);
    logo_axes = uiaxes(data_panel, ...
                       'Visible', 'on', ...
                       'Position', [275, 15, 100/logratio, 100],...
                       'XTickLabelMode', 'manual', 'YTickLabelMode', 'manual', ...
                       'XTickMode', 'manual', 'YTickMode', 'manual', ...
                       'Color', [0.85 0.85 0.87], ...
                       'GridLineStyle', 'none', ...
                       'XColor', [0.85 0.85 0.87], 'YColor', [0.85 0.85 0.87], ...
                       'BackgroundColor',[0.85 0.85 0.87]);
	
    f=image(imageLogo,'parent',logo_axes);
    axis off;
    set(f, 'AlphaData', alpha);
    
    % Click invisible overlay !
    imagesc(image_panel, maxProj, 'AlphaData', 0, 'ButtonDownFcn', @(h_map, event) click(h_map, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots, nucsOrSpots_button));
    
    % Zoom panel (if release >= 2019)
    verstr = version;
    vers = str2double(verstr(strfind(verstr, '(R')+2:strfind(verstr,')')-2));
    if vers > 2018
        image_panel.Toolbar.Visible = 'on';
    end
    
end
end
%% %%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%% %%
%% %%%%%% %%% APP FUNCTIONS %%% %%%%%% %%
%% %%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%% %%

%% Display functions
% function updateZ(z_slider,image_panel,img)
%     imagesc(image_panel, img(:,:,round(z_slider.Value)));
%     drawnow;
% end

function updateCmap(cmap_choice,image_panel)
    switch cmap_choice.Value
        case 'red'
            cmap = [0:1/255:1; zeros(2, 256)]';
        case 'green'
            cmap = [zeros(1, 256); 0:1/255:1; zeros(1, 256)]';
        case 'blue'
            cmap = [zeros(2, 256); 0:1/255:1]';
        otherwise
            cmap = cmap_choice.Value;
    end
    colormap(image_panel, cmap);
    drawnow;
end

function updateOverlay(overlay_switch, H_spots, H_nucs)
    switch overlay_switch.Value
        case 'Off'
            new = 'off';
        case 'On'
            new = 'on';
    end
    
    nbNucs = size(H_nucs,2);
    for inuc = 1:nbNucs
        H_nucs{inuc}.Visible = new;
        
        nbSpots = size(H_spots{inuc}, 2);
        if nbSpots > 0
            for isp = 1:nbSpots
                H_spots{inuc}{1, isp}.Visible = new;
            end
        end
    end

end

function updateInt(low_slider,high_slider,image_panel)
    if low_slider.Value > high_slider.Value
        high_range = diff(high_slider.Limits);
        high_slider.Value = low_slider.Value+high_range/100;
    end
    if high_slider.Value < low_slider.Value
        low_range = diff(low_slider.Limits);
        low_slider.Value = high_slider.Value-low_range/100;
    end
    image_panel.CLim = [low_slider.Value, high_slider.Value];
    drawnow;
end

function updateProj(proj_knob, image_panel, H_map, stdProj, maxProj, sumProj, meanProj, medProj, minProj, low_slider, high_slider)
    
    switch proj_knob.Value
        case 'Std deviation'
            H_map.CData = stdProj;
        case 'Maximum'
            H_map.CData = maxProj;
        case 'Sum'
            H_map.CData = sumProj;
        case 'Minimum'
            H_map.CData = minProj;
        case 'Mean'
            H_map.CData = meanProj;
        case 'Median'
            H_map.CData = medProj;
    end
    
    mincur = min(H_map.CData(:));
    maxcur = max(H_map.CData(:));
    low_slider.Limits = [mincur maxcur];
    high_slider.Limits = [mincur maxcur];
    low_slider.Value = mincur;
    high_slider.Value = maxcur;
    image_panel.CLim = [mincur, maxcur];
    H_map.CDataMapping = 'scaled';
    drawnow;
end

%% Data functions
function click(h_map, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots, nucsOrSpots_button)
    % Check if Nuclei or Spots
    selected = nucsOrSpots_button.Value;
    if selected
        H = cat(2, H_spots{:});
        [~, order] = sort(cat(2, H{2,:}));
        Hkill = cellfun(@isempty, H(2,:));
        H = H(1,~Hkill);
        H = H(order);
        enlarge = 2;
    else
        H = H_nucs;
        enlarge = 1;
    end
    polys = cellfun(@(cc) [cc.XData; cc.YData], H, 'UniformOutput', false);

    % Extract position inside image panel
    pos_abs = get(0, 'PointerLocation');
    pos_app = get(gcbf, 'Position');
    
    % Check for dual screen misconfiguration
    MonitorPositions = get(0,'MonitorPosition');
    if size(MonitorPositions,1) == 2
        if MonitorPositions(2,1) < 0
            pos_app(1) = pos_app(1) + MonitorPositions(2,1);
            MonitorPositions(1,1) = MonitorPositions(2,3)+1;
            MonitorPositions(2,1) = 1;
        end
    end
    
    pos_axes = h_map.Parent.Position;%pos_app + [25, 0, -625, -75];%h_map.Parent.Position; only if didn't move
    pos_axes(1:2) = pos_axes(1:2) + pos_app(1:2);
    
    % Correct for aspect ratio difference
    dims = size(h_map.CData);
    
    ratio_app = pos_axes(3)/pos_axes(4);%h_map.Parent.Position(3)/h_map.Parent.Position(4);
    ratio_img = dims(2)/dims(1);

    if ratio_img/ratio_app > 1
        pos_im = pos_axes;
        pos_im(4) = pos_im(4)*ratio_app/ratio_img;
        pos_im(2) = pos_im(2) + (pos_axes(4) - pos_im(4))/2;
        zoom = pos_im(3)/dims(2);
    else
        pos_im = pos_axes;
        pos_im(3) = pos_im(3)*ratio_img/ratio_app;
        pos_im(1) = pos_im(1) + (pos_axes(3) - pos_im(3))/2;
        zoom = pos_im(4)/dims(1);
    end
    
    pos_inside = pos_abs - pos_im(1:2);
    pos_inside(2) = pos_im(4) - pos_inside(2);
    
    % Extract position as image pixel coordinates
    posi = pos_inside./zoom;
    
    % Get associated object
    empty = cellfun(@isempty, H);
    polys = polys(~empty);
    polys = cellfun(@(cc) dilatepoly(cc, enlarge), polys, 'UniformOutput', false);
    isinside = cellfun(@(cc) inpolygon(posi(1), posi(2), cc(1,:), cc(2,:)), polys);
    
    [val, sel_temp_id] = max(isinside);
    if ~val
        if selected
            disp('Click inside a spot to select it');
        else
            disp('Click inside a nucleus to select it');
        end
        return;
    end
    
    sel_id = find(cumsum(~empty) == sel_temp_id, 1, 'first');
    
    % Update lists
    if selected
        ListKept = cellfun(@str2double, spots_box_kept.Items);

        isKept = ismember(sel_id, ListKept);
        if isKept
            spots_box_kept.Value = num2str(sel_id);
            updateOtherListSpots(spots_box_killed, spots_box_kept, h_map.Parent, H_nucs, H_spots, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots);
        else
            spots_box_killed.Value = num2str(sel_id);
            updateOtherListSpots(spots_box_kept, spots_box_killed, h_map.Parent, H_nucs, H_spots, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots);
        end
        
    else
        ListKept = cellfun(@str2double, nucs_box_kept.Items);
    
        isKept = ismember(sel_id, ListKept);
        if isKept
            nucs_box_kept.Value = num2str(sel_id);
            updateOtherList(nucs_box_killed, nucs_box_kept, h_map.Parent, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr);
        else
            nucs_box_killed.Value = num2str(sel_id);
            updateOtherList(nucs_box_kept, nucs_box_killed, h_map.Parent, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr);
        end
    end
    
end

function outs = dilatepoly(lines, scale)
    if isstruct(lines)
        outs = lines;
        outs.XData = scale*lines.XData+(1-scale)*mean(lines.XData(1:end-1));
        outs.YData = scale*lines.YData+(1-scale)*mean(lines.YData(1:end-1));
    else
        outs(1,:) = scale*lines(1,:)+(1-scale)*mean(lines(1,1:end-1));
        outs(2,:) = scale*lines(2,:)+(1-scale)*mean(lines(2,1:end-1));
    end
end

function nucsOrSpots(nucsOrSpots_button, ~, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, H_nucs, H_spots, Ckeep, Ckill)
    % Change selection type
    val = nucsOrSpots_button.Value;
    
    if val
        nucsOrSpots_button.Text = 'Spots';
        
        spots_box_kept.Enable = 'on';
        spots_box_killed.Enable = 'on';
        nucs_box_kept.Enable = 'off';
        nucs_box_killed.Enable = 'off';
    else
        nucsOrSpots_button.Text = 'Nuclei';
        
        spots_box_kept.Enable = 'off';
        spots_box_killed.Enable = 'off';
        nucs_box_kept.Enable = 'on';
        nucs_box_killed.Enable = 'on';
    end
    
    % Deselect everything
    nucs_box_kept.Value = {};
    nucs_box_killed.Value = {};
    spots_box_kept.Value = {};
    spots_box_killed.Value = {};
    
    % De-color everything
    ListKept = cellfun(@str2double, nucs_box_kept.Items);
    ListKilled = cellfun(@str2double, nucs_box_killed.Items);
    for ic = 1:length(ListKept)
        inuc = ListKept(ic);
        H_nucs{inuc}.Color = Ckeep;
        for ics = 1:size(H_spots{inuc}, 2)
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckeep.*2);
        end
    end
    for ic = 1:length(ListKilled)
        inuc = ListKilled(ic);
        H_nucs{inuc}.Color = Ckill;
        for ics = 1:size(H_spots{inuc}, 2)
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
        end
    end
    ListKilled = cellfun(@str2double, spots_box_killed.Items);
    nucswithspots = find(cellfun(@length, H_spots));
    spotspernucs = cellfun(@(x) size(x,2), H_spots(nucswithspots));
    nucswithspotspernucs = repelem(nucswithspots, spotspernucs);
    spotsorder = cat(2, H_spots{:});
    [~, spotsorder] = sort(cat(2,spotsorder{2,:}));
    spots2nucs = nucswithspotspernucs(spotsorder);
    for ic = 1:length(ListKilled)
        ispot = ListKilled(ic);
        inuc = spots2nucs(ispot);
        ics = sum(spots2nucs(1:ispot)==inuc);
        H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
    end
    
end

function updateOtherListSpots(list2mod, listmod, image_panel, H_nucs, H_spots, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr, CcurrSpots)
    list2mod.Value = {};
    
    % Re-Color previous guy
    ListKept = cellfun(@str2double, spots_box_kept.Items);
    ListKilled = cellfun(@str2double, spots_box_killed.Items);
    
    nucswithspots = find(cellfun(@length, H_spots));
    spotspernucs = cellfun(@(x) size(x,2), H_spots(nucswithspots));
    nucswithspotspernucs = repelem(nucswithspots, spotspernucs);
    spotsorder = cat(2, H_spots{:});
    [~, spotsorder] = sort(cat(2,spotsorder{2,:}));
    spots2nucs = nucswithspotspernucs(spotsorder);
    for ic = 1:length(ListKept)
        ispot = ListKept(ic);
        inuc = spots2nucs(ispot);
        ics = sum(spots2nucs(1:ispot)==inuc);
        H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckeep.*2);
    end
    for ic = 1:length(ListKilled)
        ispot = ListKilled(ic);
        inuc = spots2nucs(ispot);
        ics = sum(spots2nucs(1:ispot)==inuc);
        H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
    end
    
    % Color new guys
    selectedVals = listmod.Value;
    hold(image_panel, 'on');
    for ic = 1:length(selectedVals)
        ispot = str2double(selectedVals{ic});
        inuc = spots2nucs(ispot);
        ics = sum(spots2nucs(1:ispot)==inuc);
        H_spots{inuc}{1, ics}.Color = min([1, 1, 1], CcurrSpots.*2);  
    end
    drawnow;
    
end

function updateOtherList(list2mod, listmod, image_panel, H_nucs, H_spots, nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, Ckeep, Ckill, Ccurr)
    list2mod.Value = {};
    
    % Re-Color previous guy
    ListKept = cellfun(@str2double, nucs_box_kept.Items);
    ListKilled = cellfun(@str2double, nucs_box_killed.Items);
    for ic = 1:length(ListKept)
        inuc = ListKept(ic);
        H_nucs{inuc}.Color = Ckeep;
        for ics = 1:size(H_spots{inuc}, 2)
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckeep.*2);
        end
    end
    for ic = 1:length(ListKilled)
        inuc = ListKilled(ic);
        H_nucs{inuc}.Color = Ckill;
        for ics = 1:size(H_spots{inuc}, 2)
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
        end
    end
    
    % Color new guys
    selectedVals = listmod.Value;
    hold(image_panel, 'on');
    for ic = 1:length(selectedVals)
        inuc = str2double(selectedVals{ic});
        H_nucs{inuc}.Color = Ccurr;
        for ics = 1:size(H_spots{inuc}, 2)
            H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ccurr.*2);
        end
    end
    
    % Re-color killed spots according to spot lists
    ListKilled = cellfun(@str2double, spots_box_killed.Items);
    nucswithspots = find(cellfun(@length, H_spots));
    spotspernucs = cellfun(@(x) size(x,2), H_spots(nucswithspots));
    nucswithspotspernucs = repelem(nucswithspots, spotspernucs);
    spotsorder = cat(2, H_spots{:});
    [~, spotsorder] = sort(cat(2,spotsorder{2,:}));
    spots2nucs = nucswithspotspernucs(spotsorder);
    for ic = 1:length(ListKilled)
        ispot = ListKilled(ic);
        inuc = spots2nucs(ispot);
        ics = sum(spots2nucs(1:ispot)==inuc);
        H_spots{inuc}{1, ics}.Color = min([1, 1, 1], Ckill.*2);
    end
    
    drawnow;
    
end

function updateNucsOrSpots(nucs_box_kept, nucs_box_killed, spots_box_kept, spots_box_killed, ~, ~, nucsOrSpots_button, H_nucs, H_spots)
    
    val = nucsOrSpots_button.Value;
    if val
        ListKept = cellfun(@str2double, spots_box_kept.Items);
        ListKilled = cellfun(@str2double, spots_box_killed.Items);
        
        if isempty(spots_box_killed.Value) % Kill items
            listmod = spots_box_kept.Value;
            v2mod = cellfun(@str2double, listmod);
            lia = ismember(ListKept, v2mod);
            ListKept(lia) = []; % Remove from kept list
            ListKilled = sort(cat(2, ListKilled, v2mod)); % Add to killed list

            % Update list boxes
            spots_box_kept.Items = mat2cell(num2str(ListKept'), ones(length(ListKept), 1))';
            spots_box_kept.Items = cellfun(@strtrim, spots_box_kept.Items, 'UniformOutput', false);
            spots_box_killed.Items = mat2cell(num2str(ListKilled'), ones(length(ListKilled), 1))';
            spots_box_killed.Items = cellfun(@strtrim, spots_box_killed.Items, 'UniformOutput', false);
            spots_box_killed.Value = cellfun(@strtrim, listmod, 'UniformOutput', false);
            spots_box_kept.Value = {};
        else % Save items                                                         
            listmod = spots_box_killed.Value;
            v2mod = cellfun(@str2double, listmod);
            lia = ismember(ListKilled, v2mod);
            ListKilled(lia) = []; % Remove from killed list
            ListKept = sort(cat(2, ListKept, v2mod)); % Add to kept list

            % Update list boxes
            spots_box_killed.Items = mat2cell(num2str(ListKilled'), ones(length(ListKilled), 1))';
            spots_box_killed.Items = cellfun(@strtrim, spots_box_killed.Items, 'UniformOutput', false);
            spots_box_kept.Items = mat2cell(num2str(ListKept'), ones(length(ListKept), 1))';
            spots_box_kept.Items = cellfun(@strtrim, spots_box_kept.Items, 'UniformOutput', false);
            spots_box_kept.Value = cellfun(@strtrim, listmod, 'UniformOutput', false);
            spots_box_killed.Value = {};
        end
        
    else
        ListKept = cellfun(@str2double, nucs_box_kept.Items);
        ListKilled = cellfun(@str2double, nucs_box_killed.Items);
        
        if isempty(nucs_box_killed.Value) % Kill items
            listmod = nucs_box_kept.Value;
            v2mod = cellfun(@str2double, listmod);
            lia = ismember(ListKept, v2mod);
            ListKept(lia) = []; % Remove from kept list
            ListKilled = sort(cat(2, ListKilled, v2mod)); % Add to killed list

            % Update list boxes
            nucs_box_kept.Items = mat2cell(num2str(ListKept'), ones(length(ListKept), 1))';
            nucs_box_kept.Items = cellfun(@strtrim, nucs_box_kept.Items, 'UniformOutput', false);
            nucs_box_killed.Items = mat2cell(num2str(ListKilled'), ones(length(ListKilled), 1))';
            nucs_box_killed.Items = cellfun(@strtrim, nucs_box_killed.Items, 'UniformOutput', false);
            nucs_box_killed.Value = cellfun(@strtrim, listmod, 'UniformOutput', false);
            nucs_box_kept.Value = {};

            % Kill the corresponding spots
            listKeptSpots = cellfun(@str2double, spots_box_kept.Items);
            listKilledSpots = cellfun(@str2double, spots_box_killed.Items);
            nucswithspots = find(cellfun(@length, H_spots));
            spotspernucs = cellfun(@(x) size(x,2), H_spots(nucswithspots));
            nucswithspotspernucs = repelem(nucswithspots, spotspernucs);
            spotsorder = cat(2, H_spots{:});
            [~, spotsorder] = sort(cat(2,spotsorder{2,:}));
            spots2nucs = nucswithspotspernucs(spotsorder);
            
            ics = find(spots2nucs == v2mod);
            if ~isempty(ics)
                lia = ismember(listKeptSpots, ics);
                listKeptSpots(lia) = [];
                listKilledSpots = unique(cat(2, listKilledSpots,  ics));

                spots_box_kept.Items = mat2cell(num2str(listKeptSpots'), ones(length(listKeptSpots), 1))';
                spots_box_kept.Items = cellfun(@strtrim, spots_box_kept.Items, 'UniformOutput', false);
                spots_box_killed.Items = mat2cell(num2str(listKilledSpots'), ones(length(listKilledSpots), 1))';
                spots_box_killed.Items = cellfun(@strtrim, spots_box_killed.Items, 'UniformOutput', false);
                spots_box_killed.Value = {};
                spots_box_kept.Value = {};
            end
            
        else % Save items                                                         
            listmod = nucs_box_killed.Value;
            v2mod = cellfun(@str2double, listmod);
            lia = ismember(ListKilled, v2mod);
            ListKilled(lia) = []; % Remove from killed list
            ListKept = sort(cat(2, ListKept, v2mod)); % Add to kept list

            % Update list boxes
            nucs_box_killed.Items = mat2cell(num2str(ListKilled'), ones(length(ListKilled), 1))';
            nucs_box_killed.Items = cellfun(@strtrim, nucs_box_killed.Items, 'UniformOutput', false);
            nucs_box_kept.Items = mat2cell(num2str(ListKept'), ones(length(ListKept), 1))';
            nucs_box_kept.Items = cellfun(@strtrim, nucs_box_kept.Items, 'UniformOutput', false);
            nucs_box_kept.Value = cellfun(@strtrim, listmod, 'UniformOutput', false);
            nucs_box_killed.Value = {};
        end
    end
    
    drawnow;
end

function fcnOK(spots_box_killed, spots_box_kept, nucs_box_killed, nucs_box_kept, filterFig, ~, ~)
    % Get nuclear lists
    ListKept = cellfun(@str2double, spots_box_kept.Items);
    ListKilled = cellfun(@str2double, spots_box_killed.Items);
    
    % Check for conflicts
    if ~isempty(intersect(ListKept, ListKilled))
        warning('Conflicts between the list of killed and kept spots ! The kept list will be privileged');
    end
    
    % Get file name and open data
    fname = filterFig.Tag;
    [X,Y,Z,Volume,Intensity,Relativeintensity,Distancetoborder,Relativedistancetoborder,NucleusID,Nucleusintensity,Nucleusbackground,NucleusVolume,Nucleusspotnumber,~] = importSpots(fname);
    nbSpots = size(X, 1);
    
    % Find to spots to keep/kill
%     toKill = double(ismember(NucleusID, ListKilled));
    toKeep = double(ismember(1:nbSpots, ListKept)); 
    
    % Export modified *.spots file
    fid = fopen(fname, 'w');
    fprintf(fid, 'X\tY\tZ\tVolume\tIntensity\tRelative intensity\tDistance to border\tRelative distance to border\tNucleus ID\tNucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\n');
    for ii = 1:nbSpots
        fprintf(fid, '%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.4f\t%.0f\t%.2f\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f', ...
            X(ii),Y(ii),Z(ii),Volume(ii),Intensity(ii),Relativeintensity(ii),Distancetoborder(ii),Relativedistancetoborder(ii),NucleusID(ii),Nucleusintensity(ii),Nucleusbackground(ii),NucleusVolume(ii),Nucleusspotnumber(ii),toKeep(ii));
        if ii < nbSpots
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    % Export modified *.nucs file
    fname = replace(fname, '.spots', '.nucs');
    
    ListKept = cellfun(@str2double, nucs_box_kept.Items);
    ListKilled = cellfun(@str2double, nucs_box_killed.Items);
    
    % Check for conflicts
    if ~isempty(intersect(ListKept, ListKilled))
        warning('Conflicts between the list of killed and kept nuclei ! The kept list will be privileged');
    end
    
    [Nucleusintensity_nucs,Nucleusbackground_nucs,NucleusVolume_nucs,Nucleusspotnumber_nucs,~] = importNucs(fname);
    nbNucs = size(Nucleusintensity_nucs, 1);
    
    toKeep = double(ismember(1:nbNucs, ListKept)); 
    
    fid = fopen(fname, 'w');
    fprintf(fid, 'Nucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\n');
    for ii = 1:nbNucs
        fprintf(fid, '%d\t%d\t%.0f\t%.0f\t%.0f', ...
            Nucleusintensity_nucs(ii),Nucleusbackground_nucs(ii),NucleusVolume_nucs(ii),Nucleusspotnumber_nucs(ii),toKeep(ii));
        if ii < nbNucs
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    close(gcbf);
end
