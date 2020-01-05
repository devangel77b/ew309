%% Circular Error Probable Computation
% ES309
% MIDN 2/C 
% MIDN 2/C
% Date:
%  This file takes data at ONE range and computes the bias, CEP, and 95%
%  confidence intervals
%

%% Prep MATLAB workspace
clc;
clear   all;
close   all;
warning off;
%
%% Load data from excel file 

[filename,pathname] = uigetfile('*.xlsx','Select appropriate data file');
%parse the numerical, text, and raw data 
[num,txt,raw] = xlsread(filename);

%you can also load data from a MATLAB file, for example
%load('rg_13ft.mat');
%%
%Data - recall this is for a single distance but easily extractable for
%multiple distances

x = num(:,1);
y = num(:,2);
range = txt(1,1); %string that lets you know what distance this data is taken
n = length(x); %number of samples

%% Basic Sample Statistics

%Sample mean
xbar = 
ybar = 

%Sample Standard Deviation
Sx = 
Sy = 

% 95% confidence on mean
Sxbar = 
Sybar = 

% Bias corrected data
bias_corr = 
xB = 
yB = 
SxB = 
SyB = 


% bias and precision standard deviations
sigma_P = 
sigma_B = 
sigma = 

% CEP (bias corrected)
cep = 

%95% confiednce on CEP
Sx_low   = 
Sx_high  = 
Sy_low   = 
Sy_high  = 

errlow  = 
errhigh = 

CEP_lower  = 
CEP_upper  = 

%% Plotting

f1 = figure(1);
p1 = plot(x,y,'o','markerfacecolor','g');
grid on; 
hold on;

% Label axes
xlabel('(m)');
ylabel('(m)');
%title_txt = ['CEP = ',num2str(cep,'%.3f'),'(m)','range =',range];
title_txt = ['CEP = ',num2str(cep,'%.3f'),'(m)'];
h = title(title_txt);

CEP_cir = viscircles([0 0],cep,'EdgeColor','b');

% Draw 95% confidence boundaries
CEP_upp = viscircles([0 0],cep+CEP_upper,'EdgeColor','g');

% A radius less than zero does not make sense.
if(cep-CEP_lower>0)
  CEP_low = viscircles([0 0],cep-CEP_lower,'EdgeColor','g');
end

% Make axes limits same
curr_axes = gca;
x_lim = get(curr_axes,'XLim');
y_lim = get(curr_axes,'YLim');
min_pt = min(x_lim(1),y_lim(1));
axis([min_pt -min_pt min_pt -min_pt]);
axis square;

% Legend
legend('POI');