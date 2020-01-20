%% @authors: Sara Garbarino & Mara Scussolini
% email": garabrino@dima.unige.it; scussolini@dima.unige.it

function  main(action,varargin)
clc
warning('off')

if nargin<1,
    
    global analysis_gui

    if ispc==1, analysis_gui.slash='\'; else analysis_gui.slash='/'; end
    
    analysis_gui.start_path = pwd;
    analysis_gui.path = [];
    
    pathfile=[analysis_gui.start_path analysis_gui.slash 'CONFIG' analysis_gui.slash 'path.txt'];
    if exist(pathfile,'file')
        fo=fopen(pathfile,'r');
        gui_pyr.path=fgetl(fo);
        fclose(fo);
    else
        SetPath_Callback;
    end
    
    if exist(pathfile,'file')
        comp_g;
    end
    
else
    feval(action,varargin{:});
end

end

function analysis_gui = comp_g()

global analysis_gui

analysis_gui.OUTPUTfolder = [];
analysis_gui.DATAfolder = [];

analysis_gui.fig = figure('Visible','off','Position',[350,350,450,350],'MenuBar','none'); %,'Toolbar','figure');

analysis_gui.setplot = uicontrol('Style','pushbutton','String','Model',...
    'Position',[175,285,90,30], ...
    'HorizontalAlignment','Right',...
    'Callback',{@setPLOT_Callback});

analysis_gui.setoutput = uicontrol('Style','pushbutton','String','set OUTPUT',...
    'Position',[110,241.5,225,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setOUTPUT_Callback});

analysis_gui.setdata = uicontrol('Style','pushbutton','String','Load Target Tissue Data',...
    'Position',[110,195,225,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setDATA_Callback});

analysis_gui.setRT = uicontrol('Style','pushbutton','String','Load Reference Tissue Data',...
    'Position',[110,150,225,25], ...
    'HorizontalAlignment','Right',...
    'Callback',{@setRT_Callback});

analysis_gui.setif = uicontrol('Style','pushbutton','String','Load asymptotic blood data',...
    'Position',[110,105,225,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setIF_Callback});

analysis_gui.setstart = uicontrol('Style','pushbutton','String','START',...
    'Position',[110,30,100,50],...
    'HorizontalAlignment','Right',...
    'Callback',{@setSTART_Callback});

analysis_gui.setexit = uicontrol('Style','pushbutton','String','EXIT',...
    'Position',[235,30,100,50],...
    'HorizontalAlignment','Right',...
    'Callback',{@setEXIT_Callback});

%-------------------------------------------------------------------------%
set([analysis_gui.fig,analysis_gui.setoutput,analysis_gui.setplot,... 
    analysis_gui.setdata,analysis_gui.setif,analysis_gui.setRT,...
    analysis_gui.setstart, analysis_gui.setexit],...
    'Units','normalized');

set(analysis_gui.fig,'Name','Compartmental Analysis - Reference Tissue Model GUI');
set(analysis_gui.fig,'NumberTitle','off');
movegui(analysis_gui.fig,'center');
%-------------------------------------------------------------------------%
set(analysis_gui.fig,'Visible','on');

end

%% Callback

%-------------------------------------------------------------------------%
function SetPath_Callback(hObject, evendata, handles)

global analysis_gui;

if isempty(analysis_gui.path),
    str='No path selected. Click Yes to procede...';
else
    str=sprintf('Current path:\n %s\n\n Click Yes to change it...',analysis_gui.path);
end

choice = questdlg(str,'path setting:', ...
    'Yes','No','Yes');
switch choice
    case 'Yes'
        path=uigetdir(analysis_gui.start_path, 'Choose path');
        
        if path == 0
        else
            analysis_gui.path=path;
            config_dir=[analysis_gui.start_path analysis_gui.slash 'CONFIG'];
            if ~exist(config_dir,'dir'), mkdir(config_dir); end
            
            fo=fopen([config_dir,analysis_gui.slash 'path.txt'],'w');
            fprintf(fo,'%s',path);
            fclose(fo);
        end
    case 'No'
end

end

%-------------------------------------------------------------------------%
function setOUTPUT_Callback(hObject, evendata, handles)

global analysis_gui

OUTPUTfolder=uigetdir(analysis_gui.start_path,'Choose OUTPUT folder');

set([analysis_gui.setoutput],'String',OUTPUTfolder,'Units','normalized');

if OUTPUTfolder == 0
else
    analysis_gui.OUTPUTfolder=OUTPUTfolder;    
    if ~exist(analysis_gui.OUTPUTfolder,'dir') 
        mkdir(analysis_gui.OUTPUTfolder);
    end
end

end

%-------------------------------------------------------------------------%
function setPLOT_Callback(hObject, evendata, handles)

global analysis_gui

pathfigure=[analysis_gui.start_path analysis_gui.slash 'model' analysis_gui.slash];

I = imread([pathfigure 'RTM_model_eposter.jpg']);
figure('numbertitle', 'off');
imshow(I);

end

%-------------------------------------------------------------------------%
function setDATA_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_TT,location_TT]=uigetfile('*.VOISTAT','Select Data');

set([analysis_gui.setdata],'String',namefile_TT,'Units','normalized');
if namefile_TT ~= 0
    if ~iscell(namefile_TT)
        if namefile_TT==0
        else
            if isfield(analysis_gui,'namefile_TT'), analysis_gui=rmfield(analysis_gui,'namefile_TT'); end
            if isfield(analysis_gui,'location_TT'), analysis_gui=rmfield(analysis_gui,'location_TT');  end
        end
    end
    
    analysis_gui.DATA = strcat(location_TT,namefile_TT);
    analysis_gui.DATA = struct('model','RTM');
    
    delimiter = '\t';
    startRow = 4;
    
    formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
    fileID = fopen(strcat(location_TT,namefile_TT),'r');
    dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
    fclose(fileID);
    
    analysis_gui.DATA.ROI.Study_Name = dataArray{:, 1};
    analysis_gui.DATA.ROI.Study_Name = analysis_gui.DATA.ROI.Study_Name(1,:);
    
    analysis_gui.DATA.ROI.Time_ROI = dataArray{:, 4};  
    analysis_gui.DATA.ROI.Averaged_ROI = dataArray{:, 7};    
    analysis_gui.DATA.ROI.Volume_ROI = dataArray{:, 11};
    analysis_gui.DATA.ROI.Volume_ROI = analysis_gui.DATA.ROI.Volume_ROI(1,:);
    
else  set([analysis_gui.setdata],'String','Choose Data','Units','normalized');
end

end

%-------------------------------------------------------------------------%
function setRT_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_RT,location_RT]=uigetfile('*.VOISTAT','Select Data');

if namefile_RT ~= 0
set([analysis_gui.setRT],'String',namefile_RT,'Units','normalized');

if ~iscell(namefile_RT)
    if namefile_RT==0
    else
        if isfield(analysis_gui,'namefile_RT'), analysis_gui=rmfield(analysis_gui,'namefile_RT'); end
        if isfield(analysis_gui,'location_RT'), analysis_gui=rmfield(analysis_gui,'location_RT');  end
    end
end

analysis_gui.RT = strcat(location_RT,namefile_RT);

delimiter = '\t';
startRow = 4;

formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
fileID = fopen(strcat(location_RT,namefile_RT),'r');
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
fclose(fileID);

analysis_gui.RT = struct('Averaged_RT',dataArray{:, 7});
analysis_gui.RT.Volume_RT = dataArray{:, 11};
analysis_gui.RT.Volume_RT = analysis_gui.RT.Volume_RT(1,:);

else 
    set([analysis_gui.setRT],'String','Choose RT Data','Units','normalized');
end

end

%-------------------------------------------------------------------------%
function setIF_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_if,location_if]=uigetfile('*.VOISTAT','Select Data');
if namefile_if ~= 0
set([analysis_gui.setif],'String',namefile_if, ...
    'Units','normalized');

if ~iscell(namefile_if)
    if namefile_if==0
    else
        if isfield(analysis_gui,'namefile_if'), analysis_gui=rmfield(analysis_gui,'namefile_if'); end
        if isfield(analysis_gui,'location_if'), analysis_gui=rmfield(analysis_gui,'location_if');  end
    end
end

analysis_gui.IF = strcat(location_if,namefile_if);

delimiter = '\t';
startRow = 4;

formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
fileID = fopen(strcat(location_if,namefile_if),'r');
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
fclose(fileID);

analysis_gui.IF = struct('Averaged_IF',dataArray{:, 7});
analysis_gui.IF.Volume_IF = dataArray{:, 11};
analysis_gui.IF.Volume_IF = analysis_gui.IF.Volume_IF(1,:);

else  set([analysis_gui.setif],'String','Choose IF','Units','normalized');
end

end

%-------------------------------------------------------------------------%
function setSTART_Callback(hObject, evendata, handles)

global analysis_gui

if ~exist(analysis_gui.OUTPUTfolder,'dir')
    h = msgbox('Warning: You have to select OUTPUT folder');
end

if ~isfield(analysis_gui,'DATA')
    h = msgbox('Warning: You have to choose TT data');
end

if  ~isfield(analysis_gui,'RT')
    hh = msgbox('Warning: You have to choose RT data');
end

if  ~isfield(analysis_gui,'IF')
    hh = msgbox('Warning: You have to choose IF');
end

if exist(analysis_gui.OUTPUTfolder,'dir') && isfield(analysis_gui,'DATA') && isfield(analysis_gui,'RT') && isfield(analysis_gui,'IF')
    RTM_gui;
end

end

%-------------------------------------------------------------------------%
function setEXIT_Callback(hObject, evendata, handles)
clc
clear all
close all
% exit
end