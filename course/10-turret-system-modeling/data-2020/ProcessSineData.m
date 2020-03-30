%% Import data from spreadsheet
clear all
clc

fname = 'EncoderForMike_PWM_Sinev2.xlsx';

dat = readtable(fname);

time = dat{:,2};
dc = dat{:,1};
enc = dat{:,3};
pos = dat{:,4};

figure(1); clf
% subplot(2,1,1)
plot(time,pos)
axis([0 15 -inf inf])
xlabel('Time (s)')
ylabel('Orientation (rad)')
title('Experimental Data')
% subplot(2,1,2)
% plot(time,dc)
% axis([0 15 -inf inf])

figure(3); clf
plot(enc,pos)

%%
fname = 'EncoderForMike_PWM_45.xlsx';

dat = readtable(fname);

time = dat{:,2};
dc = 0.45*ones(size(time));
dc(time<=1) = 0;
enc = dat{:,3};
pos = dat{:,4};

figure(1); clf
subplot(2,1,1)
plot(time,pos)
subplot(2,1,2)
plot(time,dc)


