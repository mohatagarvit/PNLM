function varargout = actual(varargin)
% =========================================================================
%% Code for GUI of Pruned Non-Local Means
%
%  x{1}                : input grayscale image
%  x{2}                : noisy grayscale image
%  x_hat               : output denoised grayscale image
%  diff_x_hat_with_y   : sum of the partial derivative terms of x_hat with 
%                        respect to y(i) for every i
%  PSNR_noisy          : Peak Signal to Noise Ratio of noisy image
%  PSNR_pnlm           : Peak Signal to Noise Ratio PNLM-processed output 
%                        denoised grayscale image
%
%%
%  Authors     : Garvit Mohata, Sanjay Ghosh, and Kunal Narayan Chaudhury.
%
%  Reference   : S. Ghosh, A. K. Mandal, and K. N. Chaudhury, "Pruned Non-Local Means", 
%               IET Image Processing, vol. 11, no. 5, pp. 317-323, April 2017.
%
% =========================================================================
%%

% ACTUAL MATLAB code for actual.fig
%      ACTUAL, by itself, creates a new ACTUAL or raises the existing
%      singleton*.
%
%      H = ACTUAL returns the handle to a new ACTUAL or the handle to
%      the existing singleton*.
%
%      ACTUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACTUAL.M with the given input arguments.
%
%      ACTUAL('Property','Value',...) creates a new ACTUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before actual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to actual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help actual

% Last Modified by GUIDE v2.5 26-Jun-2017 15:24:13

% Begin initialization code - DO NOT EDIT
clc;
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @actual_OpeningFcn, ...
                   'gui_OutputFcn',  @actual_OutputFcn, ...
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


% --- Executes just before actual is made visible.
function actual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to actual (see VARARGIN)

% Choose default command line output for actual
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes actual wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = actual_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Retrieving value of S from user
function S_Callback(hObject, eventdata, handles)
% hObject    handle to S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of S as text
%        str2double(get(hObject,'String')) returns contents of S as a double

S = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function S_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Retrieving value of K from user
function K_Callback(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K as text
%        str2double(get(hObject,'String')) returns contents of K as a double

K = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Retrieving value of sigma from user
function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double

sigma = str2double(get(hObject,'String'));

        
% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Retrieving value of h from user
function h_Callback(hObject, eventdata, handles)
% hObject    handle to h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h as text
%        str2double(get(hObject,'String')) returns contents of h as a double

h = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Taking the imagepath entered by the user
function imagepath_Callback(hObject, eventdata, handles)
% hObject    handle to imagepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imagepath as text
%        str2double(get(hObject,'String')) returns contents of imagepath as a double


% --- Executes during object creation, after setting all properties.
function imagepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Browsing and Loading the Image
% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, user_canceled] = imgetfile;
if (~user_canceled)
    imagepath_Handle = findobj(gcbf,'Tag','imagepath');
    set(imagepath_Handle,'String',filename);
    mread = imread(get(imagepath_Handle,'String'));
    if ndims(mread)>2
        error('error : The input given is not a grayscale image. The desired input is a grayscale image.');
    end
    axes(findobj(gcbf,'Tag','clean_image'));
    cla; hold on;
    imshow(mread); axis('image', 'off');
    x{1} = double(mread);
    set(gcbf,'UserData',x);
end


%% Adding the noise and displaying the noisy image
% --- Executes on button press in add_noise.
function add_noise_Callback(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        
x                =  get(gcbf,'UserData');
editsigma_Handle =  findobj(gcbf,'Tag','sigma');
sigma            =  eval(get(editsigma_Handle,'String'));
        
x{2}             =  x{1}+sigma.*randn(size(x{1}));
[m,n]            =  size(x{1});
peak             =  255;

axes(findobj(gcbf,'Tag','noisy_image'));
cla; hold on;
imshow(uint8(x{2})); axis('image', 'off');

PSNR_noisy       =  10 * log10(m * n * peak^2 / sum(sum((x{2} - x{1}).^2)) );
set(handles.psnr_noisy,'String',PSNR_noisy);
set(gcbf,'UserData',x);


%% Processing for pnlm denoised image & displaying it & its method noise
% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x                =  get(gcbf,'UserData');

editS_Handle     =  findobj(gcbf,'Tag','S');
S                =  str2double(get(editS_Handle,'String'));
editK_Handle     =  findobj(gcbf,'Tag','K');
K                =  str2double(get(editK_Handle,'String'));
editsigma_Handle =  findobj(gcbf,'Tag','sigma');
sigma            =  eval(get(editsigma_Handle,'String'));
edith_Handle     =  findobj(gcbf,'Tag','h');
h                =  str2double(get(edith_Handle,'String'));

h                =  h*sigma;
alpha            =  100;
[m,n]            =  size(x{1});

y_padded         =  padarray(x{2},[S+K S+K],'symmetric');
peak             =  255;

[x_hat,diff_x_hat_with_y]  =  GUI_mex(y_padded,sigma,S,K,alpha,h);
PSNR_pnlm                  =  10 * log10(m * n * peak^2 / sum(sum((x_hat - x{1}).^2)) );
set(handles.psnr_pnlm,'String',PSNR_pnlm);
%SURE       =  sum(sum((x_hat - y).^2))/(m*n) - sigma^2 + 2*sigma^2*diff_x_hat_with_y/(m*n);
%MSE_pnlm   =  sum(sum((x_hat-x).^2))/(m*n);
%SSIM_pnlm  =  ssim(x,x_hat);  

axes(findobj(gcbf,'Tag','pnlm'));
cla; hold on;
imshow(uint8(x_hat)); axis('image', 'off');
axes(findobj(gcbf,'Tag','method_noise'));
cla; hold on;
imshow(uint8(x{2}-x_hat)); axis('image', 'off');

%% Saving Images for future use
% --- Executes on selection change in saveimage.
function saveimage_Callback(hObject, eventdata, handles)
% hObject    handle to saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns saveimage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from saveimage
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val};
case 'Clean Image' % User selects Clean Image to save.
   imsave(findobj(gcbf,'Tag','clean_image'));
case 'Noisy Image' 
   imsave(findobj(gcbf,'Tag','noisy_image'));
case 'PNLM' 
   imsave(findobj(gcbf,'Tag','pnlm'));
case 'Method Noise' 
   imsave(findobj(gcbf,'Tag','method_noise'));
case 'All'
   imsave(findobj(gcbf,'Tag','clean_image'));
   imsave(findobj(gcbf,'Tag','noisy_image'));
   imsave(findobj(gcbf,'Tag','pnlm'));
   imsave(findobj(gcbf,'Tag','method_noise'));
end

% --- Executes during object creation, after setting all properties.
function saveimage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Closing the GUI window
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
