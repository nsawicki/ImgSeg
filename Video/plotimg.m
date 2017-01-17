%%  plotimg.m
%   Plot and save desirable figures used in the report
%
%   Daniel Vial, Yang Xiao

clear; close all;

load('man_fullpipe_gamma=0.1_xThresh=0.0002_mult=500.mat')

%% TODO: Plot the images you want to put into the report
imagesc(f(:,:,:,1));

%%  Display figures in true size without extra white areas
truesize;
set(gca,'position',[0,0,1,1],'units','normalized');
set(gca,'XTick',''); set(gca,'YTick','');
saveas(gcf,'man/xTract-white-unseg.png');
