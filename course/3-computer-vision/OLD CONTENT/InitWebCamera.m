function usbCam = InitWebCamera()
% Initializes iCube Webcam and 
% Returns an image aquisition object webCam
% Can now get frames with im = getsnapshot(usbCam);
% NOTES:
% If you have multiple cameras installed you need to change the "1" on line 11
% If you want to do this for a different model of camera use imaqtool to
% see the options for the format in line 11 (ex. YUY2_640X480) 
% Esposito 3/31/2012
disp('Initializing Camera...');
imaqreset;
usbCam = videoinput('winvideo',1, 'YUY2_640X480'); %YOU MAY NEED TO CHANGE 
% THE 1 to a 2 IF YOUR LAPTOP HAS A BUILT IN CAMERA


set(usbCam, 'SelectedSourceName', 'input1');  
set(usbCam, 'ReturnedColorSpace', 'rgb');

%% These lines provide an example of how to change a camera property like 
% auto white balance for example.  Use imaqtool to figure out a list of
% proerties and their acceptable values
src_obj = getselectedsource(usbCam) 
set(src_obj, 'HueMode', 'manual')
get(src_obj)


%% this is important to get a fast frame rate
triggerconfig(usbCam, 'manual'); 
start(usbCam);
preview(usbCam);
disp('Camera initialization done.  Preview can be closed if desired.')