function [cam,prv] = initCamera
% INITCAMERA initialize a camera object and opens a preview.
%
%   cam = INITCAMERA returns the camera object. This is required for using
%   getsnapshot.m:
%       -> im = getsnapshot(cam);
%
%   [cam,prv] = INITCAMERA returns both the camera object and the preview
%   handle. Note that "prv" can be used to get images faster than 
%   getsnapshot.m using:
%       -> im = get(prv,'CData');
%
%   NOTE: This requires an installed version of the "USB Webcams" package
%       >> supportPackageInstaller
%           -> select "Install from Internet"
%           -> select “OS Generic Video Interface”
%           -> login to mathworks using email and password
%           -> Install
%
%   M. Kutzer, 02Mar2016, USNA

% Updates
%   18Jan2017 - Updated documentation
%   08Jan2019 - Updated documentation

%% Declare persistent variable to declare new camera names
% TODO - remove persistent and replace with device selection
persistent callNum

if isempty(callNum)
    callNum = 1;
end

%% Check for installed adapters
goodAdaptor = false;
info = imaqhwinfo;
for i = 1:numel(info.InstalledAdaptors)
    switch lower(info.InstalledAdaptors{i})
        case 'winvideo'
            goodAdaptor = true;
            break
    end
end

if ~goodAdaptor
    error('initCam:BadAdaptor',...
        ['The "winvideo" adaptor is not detected.\n',...
        ' -> Run "supportPackageInstaller"\n',...
        ' -> Select and install "OS Generic Video Interface".\n'])
end

%% Check for cameras
devices = imaqhwinfo('winvideo');
if isempty(devices.DeviceIDs)
    error('No connected camera found');
end

n = numel(devices.DeviceIDs);
if n > 1
    for i = 1:n
        camList{i} = devices.DeviceInfo(i).DeviceName;
    end
    [camIdx,OK] = listdlg('PromptString','Select camera:',...
                      'SelectionMode','single',...
                      'ListString',camList);
    if ~OK
        error('No camera selected.');
    end
else
    camIdx = 1;
end

%% Check for existing devices
% out1 = imaqfind('Type', 'videoinput')
vids = imaqfind;
m = size(vids,2);
if m > 0
    for i = 1:m
        switch lower(vids(i).Type)
            case 'videoinput'
                if vids(i).DeviceID == camIdx
                    cam = vids(i);
                    if nargout > 1
                        prv = preview(cam);
                    end
                    return
                end
        end
    end
end

%% Check for available formats
formatIDX = 1; % Use default format if good one is unavailable
m = numel(devices.DeviceInfo(camIdx).SupportedFormats);
for i = 1:m
    switch devices.DeviceInfo(camIdx).SupportedFormats{i}
        case 'YUY2_640x480'
            formatIDX = i;
            break
    end
end
         
%% Create video input object
cam = videoinput('winvideo',camIdx,...
    devices.DeviceInfo(camIdx).SupportedFormats{formatIDX});

%% Setup camera parameters
set(cam,'ReturnedColorSpace','rgb');
set(cam,'Name',sprintf('camera%d',callNum));
callNum = callNum + 1;

%% Update camera properties
src_obj = getselectedsource(cam); 
try
    set(src_obj, 'ExposureMode', 'manual');
    set(src_obj, 'Exposure', -4);
catch
    warning('Unable to set "ExposureMode" to manual.');
end

%% Start camera and create preview
triggerconfig(cam,'manual'); 
start(cam);

if nargout > 1
    prv = preview(cam);
end