%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: G:\My Drive\Teaching\EW Courses, Spring\EW309\EW309 Spring 2020\EW309 AY20\10. Turret System Modeling\Turrett Data\EncoderForMike_PWM_Sine.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 13-Mar-2020 12:31:00

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "Sheet1";
%opts.DataRange = "A1:D1182";
opts.DataRange = "A1:D501";

% Specify column names and types
opts.VariableNames = ["PWM", "Time", "Counts", "Radians"];
opts.SelectedVariableNames = ["PWM", "Time", "Counts", "Radians"];%["VarName1", "VarName2", "VarName3", "VarName4"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Import the data
%data = readtable("EncoderForMike_PWM_45.xlsx", opts, "UseExcel", false);
data = readtable("EncoderForMike_PWM_Sinev2.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

%% Plot results
figure; plot(data.Time,data.PWM);
xlabel('Time (s)');
ylabel('PWM');

figure; plot(data.Time,data.Counts);
xlabel('Time (s)');
ylabel('Encoder Counts');

figure; plot(data.Time,data.Radians);
xlabel('Time (s)');
ylabel('Radians');