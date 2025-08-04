function NuFoQ()

clear variables

filepath = fileparts(mfilename('fullpath'));
cd(filepath);

addpath(filepath);
addpath(genpath([filepath, filesep,'utilities']));

%% ================================================= %%
disp('============= NuFoQ / Institut Curie - UMR 3664 / Mickaël GARNIER =============');

%% ================================================= %%
% User interface default colors
colorFgd=[0 0 0];
colorBgd=[0.85 0.85 0.87]; %% "Natural" color background
colorBut=[0.88 0.87 0.85];
% colorBgd=[1 1 1]; %% White color background
%text_size = 0.65;
%edit_text_size = 0.25;

% Create and hide the GUI as it is being constructed.
frontpanel = figure('Visible','on','Position',[00,000,500,400],...
    'Color',colorBgd,'Resize','on',...
    'Name', 'NuFoQ - Institut Curie',...  % Title figure
    'NumberTitle', 'off',... % Do not show figure number
    'MenuBar', 'none');
movegui(frontpanel, 'center');

% Prepare parameter variables

START_PATH=filepath;
ANALYSE_FOLDER=0;
ANALYZE_PATHNAME='';
ANALYZE_FILENAME='';
RESULTS_FILENAME='';
SAVE_IMGS=false;
TWO_SETS=false;
PARAMS={'0.065';'0.2';'2:3'; '0'; 'gfp'; 'false'; '1'};

setappdata(gcf,'START_PATH',START_PATH);
setappdata(gcf,'ANALYSE_FOLDER',ANALYSE_FOLDER);
setappdata(gcf,'ANALYZE_PATHNAME',ANALYZE_PATHNAME);
setappdata(gcf,'ANALYZE_FILENAME',ANALYZE_FILENAME);
setappdata(gcf,'RESULTS_FILENAME',RESULTS_FILENAME);
setappdata(gcf,'SAVE_IMGS', SAVE_IMGS);
setappdata(gcf,'TWO_SETS',TWO_SETS);
setappdata(gcf,'PARAMS',PARAMS);

%% ================================================= %%
% The command line below defined the spatial structure of the user
% interface. Their related functions are defined at the end of this file.

%% ================================================= %%
%
%                       Panel qFOCI
%
%% ================================================= %%

%% ================================================= %%

%% ================================================= %%
LoadData_panel = uipanel('Parent',frontpanel,'Title','Define path to folder/file to analyze',...
    'Position',[.05 .05 .5 .85],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);
uicontrol('Parent',frontpanel,'Style','text',...
    'String','NuFoQ','Visible','on',...
    'FontSize',20', ...
    'Units','normalized','Position',[0.1,0.9,0.90,0.1],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',[1 0.5216 0]);

%% ================================================= %%
LoadData_File_or_Folder_buttonGroup = uibuttongroup('Parent',LoadData_panel,...
    'Visible','on',...
    'Units','normalized','Position',[0.05,0.45,0.9,0.5],...
    'SelectionChangeFcn',@LoadData_File_or_Folder_buttonGroup_Callback, ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);

uicontrol('Parent',LoadData_File_or_Folder_buttonGroup,'Style','text',...
    'String','Parameters','Visible','on',...
    'Units','normalized','Position',[0.05,0.70,0.90,0.25],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);

% Create two radio buttons in the button group.
uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','Radio','String',' Single file analysis',...
    'Units','normalized','Position',[0.05,0.6,0.90,0.25],...
    'HandleVisibility','on', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);
uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','Radio','String',' Folder analysis',...
    'Units','normalized','Position',[0.05,0.35,0.90,0.25],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);
uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','CheckBox','String',' Use deconvolution images',...
    'Units','normalized','Position',[0.05,0.2,0.90,0.15],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Callback',@checkbox_useTwoSets_callback);
uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','CheckBox','String',' Save segmented images',...
    'Units','normalized','Position',[0.05,0.03,0.90,0.15],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Callback',@checkbox_saveimgs_callback);


%% ================================================= %%
uicontrol('Parent',LoadData_panel,'Style','pushbutton',...
    'String','Choose data location',...
    'Units','normalized','Position',[0.05,0.25,0.9,0.15],...
    'TooltipString', 'Folder containing the data or file to analysis',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Folder_pushbutton_Callback});

%% ================================================= %%
uicontrol('Parent',LoadData_panel,'Style','pushbutton',...
    'String','Choose results location',...
    'Units','normalized','Position',[0.05,0.05,0.9,0.15],...
    'TooltipString', 'Folder in which the analysis results will be written',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Result_pushbutton_Callback});

%% ================================================= %%
% To easily find the location of the icon image without considering O/S and
% display it in the GUI
pathImgUI=strcat([filepath,filesep,'institut_curie.png']);
[imageLogo, ~, alpha] = imread(pathImgUI);
Logo_panel=uipanel('Parent',frontpanel,...
    'Position',[0.60,0.63,0.35,0.25],...
    'BackgroundColor',[1 1 1],...
    'ForegroundColor',colorFgd);
ax = axes('parent',Logo_panel,'units','normalized','position',[0 0 1 1]);
f=image(imageLogo,'parent',ax);
axis off;
set(f, 'AlphaData', alpha);

%% ================================================= %%
uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','Set parameters',...
    'Units','normalized','Position',[0.60,0.51,0.35,0.1],...
    'TooltipString', 'Open a window to enter the analysis parameters',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Params_pushbutton_Callback});

%% ================================================= %%
uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','ANALYSIS', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.35,0.350,0.1215],...
    'TooltipString', 'Run the analysis on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@RunAnalysis_pushbutton_Callback});

uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','CHECK', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.225,0.350,0.1215],...
    'TooltipString', 'Supervised result filtering on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@FilterResults_pushbutton_Callback});

uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','GRAPHS', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.095,0.350,0.1215],...
    'TooltipString', 'Generate figures from the analysis conducted on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@DrawPlots_pushbutton_Callback});
% End of the qFOCIpanel
%% ===========  End of the qFOCIpanel   ============ %%
%%================================================== %%
%%================================================== %%
%%================================================== %%

end%function

%% function LoadData_File_or_Folder_buttonGroup_Callback(source,eventdata)
% Set the king of analysis: only on a single file (ANALYSE_FOLDER=0) or
% on a whole folder (ANALYSE_FOLDER=1)
function LoadData_File_or_Folder_buttonGroup_Callback(source,~)

str=get(get(source,'SelectedObject'),'String');

ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
if strcmp(str,' Single file analysis')
    ANALYSE_FOLDER=0;
    disp(' Single file analysis')
end

if strcmp(str,' Folder analysis')
    ANALYSE_FOLDER=1;
    disp(' Folder analysis')
end

setappdata(gcbf,'ANALYSE_FOLDER',ANALYSE_FOLDER);

end

%% function checkbox_useTwoSets_callback(source,eventdata)
% Set wether the segmentation is done using two sets of images or not
function checkbox_useTwoSets_callback(source,~)
    val = get(source,'Value');
    
    TWO_SETS = val;
    if val
        disp('Spots and nuclei segmentation will be done on different sets of images');
    else
        disp('All the segmentations and quantifications will be done on a single set of images');
    end
    
    setappdata(gcbf, 'TWO_SETS', TWO_SETS);
end


%% function checkbox_saveimgs_callback(source,eventdata)
% Set the saving of images: checked (true) unchecked (false)
function checkbox_saveimgs_callback(source,~)

val=get(source,'Value');

SAVE_IMGS = val;
if val
    disp(' Segmented images will be written in the result location.')
else
    disp(' Segmented images will not be written.')
end

setappdata(gcbf,'SAVE_IMGS',SAVE_IMGS);

end

%% function Folder_pushbutton_Callback(hObject, eventdata, handles)
% Open a dialog box to specify where the files to analyzed are located. 
% Depending the kind of analysis (single file or all file in a folder), the
% program ask to specify the image or the folder to be analyzed.
function Folder_pushbutton_Callback(~, ~, ~)

ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
START_PATH=getappdata(gcbf,'START_PATH');
ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
TWO_SETS = getappdata(gcbf, 'TWO_SETS');
ext = {'*.tif','TIF image (.tif)';'*.TIF','TIF image (.TIF)'};

if TWO_SETS
    if (ANALYSE_FOLDER)
        [ANALYZE_PATHNAME] = uigetdir(START_PATH,'Select the folder containing the quantification set of images.');
        [ANALYZE_PATHNAME2] = uigetdir(START_PATH,'Select the folder containing the spots segmentation set of images.');
    else
        [ANALYZE_FILENAME,ANALYZE_PATHNAME] = uigetfile(ext,'Select the quantification image file to analyse.');
        [ANALYZE_FILENAME2,ANALYZE_PATHNAME2] = uigetfile(ext,'Select the spots segmentation image file to analyse.');
        ANALYZE_FILENAME = {ANALYZE_FILENAME, ANALYZE_FILENAME2};
    end
    ANALYZE_PATHNAME = {ANALYZE_PATHNAME, ANALYZE_PATHNAME2};
else
    if (ANALYSE_FOLDER)
        [ANALYZE_PATHNAME] = uigetdir(START_PATH,'Select the folder containing the images.');
    else
        [ANALYZE_FILENAME,ANALYZE_PATHNAME] = uigetfile(ext,'Select the image file to analyse.');
    end
end

setappdata(gcbf,'ANALYZE_PATHNAME',ANALYZE_PATHNAME);
setappdata(gcbf,'ANALYZE_FILENAME',ANALYZE_FILENAME);
end

function Result_pushbutton_Callback(~, ~, ~)

START_PATH=getappdata(gcbf,'START_PATH');

[ANALYZE_PATHNAME] = uigetdir(START_PATH,'Select (or create) the output folder.');

cd(char(ANALYZE_PATHNAME));
setappdata(gcbf,'RESULTS_PATHNAME',ANALYZE_PATHNAME);
end

%% function Params_pushbutton_Callback(hObject, eventdata, handles)
function Params_pushbutton_Callback(~, ~, ~)
PARAMS=getappdata(gcbf,'PARAMS');

[params]=gui_params(PARAMS);

setappdata(gcbf,'PARAMS',params);
end


%% function RunAnalysis_pushbutton_Callback(hObject, eventdata, handles)
function RunAnalysis_pushbutton_Callback(~, ~, ~)

ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
SAVE_IMGS = getappdata(gcbf,'SAVE_IMGS');
PARAMS=getappdata(gcbf,'PARAMS');

error = NuFoQ_analysis(ANALYZE_PATHNAME, RESULTS_PATHNAME, SAVE_IMGS, PARAMS, ~ANALYSE_FOLDER, ANALYZE_FILENAME);
if error > 0
    fprintf('%d images not processed', error);
end
end

function FilterResults_pushbutton_Callback(~, ~, ~)

START_PATH = getappdata(gcbf,'START_PATH');
ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
SAVE_IMGS = getappdata(gcbf,'SAVE_IMGS');
PARAMS=getappdata(gcbf,'PARAMS');

supervised_filtering_results(ANALYZE_PATHNAME, RESULTS_PATHNAME, SAVE_IMGS, PARAMS, ~ANALYSE_FOLDER, ANALYZE_FILENAME, START_PATH);
end

function DrawPlots_pushbutton_Callback(~, ~, ~)

START_PATH = getappdata(gcbf,'START_PATH');
ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
PARAMS=getappdata(gcbf,'PARAMS');

generate_figures(ANALYZE_PATHNAME, RESULTS_PATHNAME, PARAMS, ~ANALYSE_FOLDER, ANALYZE_FILENAME, START_PATH);
end
