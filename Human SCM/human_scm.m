function varargout = human_scm(varargin)
% HUMAN_SCM MATLAB code for human_scm.fig
%      HUMAN_SCM, by itself, creates a new HUMAN_SCM or raises the existing
%      singleton*.
%
%      H = HUMAN_SCM returns the handle to a new HUMAN_SCM or the handle to
%      the existing singleton*.
%
%      HUMAN_SCM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HUMAN_SCM.M with the given input arguments.
%
%      HUMAN_SCM('Property','Value',...) creates a new HUMAN_SCM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before human_scm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to human_scm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help human_scm

% Last Modified by GUIDE v2.5 03-Mar-2015 16:02:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @human_scm_OpeningFcn, ...
    'gui_OutputFcn',  @human_scm_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before human_scm is made visible.
function human_scm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to human_scm (see VARARGIN)

% Choose default command line output for human_scm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes human_scm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = human_scm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Load Image Sets
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir('R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\',...
    'Select Image folder');

%if file selection is cancelled, pathname should be zero and nothing should happen
if folder_name == 0
    return
end

UserData.folder_name = folder_name; %save to handle so can grab from other objects
last_image = 0; %no images have been loaded yet

%---get the name of images in the foldler---%
image_list = ls([folder_name '\','*.bmp']);
order = NaN(2,36);
for il = 1:size(image_list,1);
    num  = sscanf(image_list(il,:), '%d');
    if num < 10
        str  = sprintf(image_list(il,2:end), '%s');
        if strcmpi(str(1),'.')
            row = 1;%novel image
        else
            row = 2;%manipulated or repeat image
        end
    else
        str  = sprintf(image_list(il,3:end), '%s');
        if strcmpi(str(1),'.')
            row = 1;%novel image
        else
            row = 2;%manipulated or repeat image
        end
    end
    order(row,num) = il;
end

%if organizing back-to-back
image_names = [];
for img = 1:36;
    %strtrim removes extra spaces if they exist in the string names
    image_names = [image_names {strtrim(image_list(order(1,img),:))} ...
        {strtrim(image_list(order(2,img),:))}];
end
UserData.response = NaN(1,72); %where we're going to store the data, row1 type, row 2 guess type
UserData.response(1:2:72) = 1;%1 for novel presentation

UserData.image_names = image_names;%save to can acess across objects
UserData.last_image = last_image; %no images have been loaded yet

set(hObject,'UserData',UserData);%store structure array to pushbutton object

display_next_image(folder_name,image_names,last_image,handles);

%% Button for Stating it was a repeat presentation
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
UserData  = get(pub1_h,'UserData');
image_names = UserData.image_names;
last_image =UserData.last_image;
folder_name = UserData.folder_name;

UserData.response(last_image) = 2;%type 2 for repeat
set(pub1_h,'UserData',UserData)

if last_image == 72
    %save the data
    save_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\B2B_responses\';
    pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
    UserData  = get(pub1_h,'UserData');
    Responses = UserData.response;
    image_names = UserData.image_names;
    folder_name = UserData.folder_name;
    slashes = strfind(folder_name,'\');
    Set = folder_name(slashes(end)+1:end); %get the set name
    Setnum  = num2str(Set(4:end));
    str_h = findobj('Tag','edit1');
    initial = get(str_h,'String');
    save([save_dir 'B2B_Set' num2str(Setnum) '_' initial],'Responses','image_names','Setnum')
        close all
else
    display_next_image(folder_name,image_names,last_image,handles)
end

%% Button for Stating it was a Replaced presentation
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
UserData  = get(pub1_h,'UserData');
image_names = UserData.image_names;
last_image =UserData.last_image;
folder_name = UserData.folder_name;

UserData.response(last_image) = 3;%type 2 for replaced
set(pub1_h,'UserData',UserData)

if last_image == 72
    %save the data
    save_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\B2B_responses\';
    pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
    UserData  = get(pub1_h,'UserData');
    Responses = UserData.response;
    image_names = UserData.image_names;
    folder_name = UserData.folder_name;
    slashes = strfind(folder_name,'\');
    Set = folder_name(slashes(end)+1:end); %get the set name
    Setnum  = num2str(Set(4:end));
    str_h = findobj('Tag','edit1');
    initial = get(str_h,'String');
    save([save_dir 'B2B_Set' num2str(Setnum) '_' initial],'Responses','image_names','Setnum')
        close all
else
    display_next_image(folder_name,image_names,last_image,handles)
end

%% Button for Stating it was a Moved presentation
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
UserData  = get(pub1_h,'UserData');
image_names = UserData.image_names;
last_image =UserData.last_image;
folder_name = UserData.folder_name;

UserData.response(last_image) = 4;%type 2 for moved
set(pub1_h,'UserData',UserData)

if last_image == 72
    %save the data
    save_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\B2B_responses\';
    pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
    UserData  = get(pub1_h,'UserData');
    Responses = UserData.response;
    image_names = UserData.image_names;
    folder_name = UserData.folder_name;
    slashes = strfind(folder_name,'\');
    Set = folder_name(slashes(end)+1:end); %get the set name
    Setnum  = num2str(Set(4:end));
    str_h = findobj('Tag','edit1');
    initial = get(str_h,'String');
    save([save_dir 'B2B_Set' num2str(Setnum) '_' initial],'Responses','image_names','Setnum')
        close all
else
    display_next_image(folder_name,image_names,last_image,handles)
end

%% Button for stating if not sure
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
UserData  = get(pub1_h,'UserData');
image_names = UserData.image_names;
last_image =UserData.last_image;
folder_name = UserData.folder_name;

UserData.response(last_image) = 5;%type 5 for unsure
set(pub1_h,'UserData',UserData)

if last_image == 72
    %save the data
    save_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\B2B_responses\';
    pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
    UserData  = get(pub1_h,'UserData');
    Responses = UserData.response;
    image_names = UserData.image_names;
    folder_name = UserData.folder_name;
    slashes = strfind(folder_name,'\');
    Set = folder_name(slashes(end)+1:end); %get the set name
    Setnum  = num2str(Set(4:end));
    str_h = findobj('Tag','edit1');
    initial = get(str_h,'String');
    save([save_dir 'B2B_Set' num2str(Setnum) '_' initial],'Responses','image_names','Setnum')
        close all
else
    display_next_image(folder_name,image_names,last_image,handles)
end

function display_next_image(folder_name,image_names,last_image,handles)
pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
UserData  = get(pub1_h,'UserData');
last_image = UserData.last_image;
last_image = last_image + 1;


UserData.last_image = last_image;
pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
set(pub1_h,'UserData',UserData);%store structure array to pushbutton object

if  ~isempty(strfind(image_names{last_image}(1:3),'p'))
    type = 2; %familiar image
elseif ~isempty(strfind(image_names{last_image}(1:3),'r'))
    type = 3; %replaced
elseif   ~isempty(strfind(image_names{last_image}(1:3),'m'))
    type = 4; %moved image
else
    type = 1;%novel image
end

%---load both images so computer does lag at a bad time
img1 = imread('Crosshair.jpg');
img2 = imread([folder_name '\' image_names{last_image}]);
img3 = imread('AskResponse.jpg');

%---display images
axes(handles.axes1);
cla
handles.h1=image(img1);
axis off
pause(2) %inter-trial-interval in seconds

axes(handles.axes1);
cla
handles.h1=image(img2);
axis off
pause(5)  %image presentation duration in seconds

if type ~= 1
    axes(handles.axes1);
    cla
    handles.h1=image(img3);
    axis off
    return
else
    if last_image == 72
        %save the data
        save_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\B2B_responses\';
        pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton3 that has ROI data
        UserData  = get(pub1_h,'UserData');
        Responses = UserData.response;
        image_names = UserData.image_names;
        folder_name = UserData.folder_name;
        slashes = strfind(folder_name,'\');
        Set = folder_name(slashes(end)+1:end); %get the set name
        Setnum  = num2str(Set(4:end));
        str_h = findobj('Tag','edit1');
        initial = get(str_h,'String');
        save([save_dir 'B2B_Set' num2str(Setnum) '_' initial],'Responses','image_names','Setnum')
            close all
    else
        display_next_image(folder_name,image_names,last_image,handles)
    end
end