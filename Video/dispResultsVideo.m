%%  dispResultsVideo.m
%   Display the segmented video frames and performance measures
%   
%   Anish Lahiri, Yang Xiao

clear; close all;

%%  Init
load('../man_fullpipe_gamma=0.1_xThresh=0.0002_mult=500.mat');

extractTimes1 = [naiveTimes(1); extractTimes];
combTimes = extractTimes1 + smoothTimes;
totTimes_naive = cumsum(naiveTimes(1:end));
totTimes_ext = cumsum(combTimes);

totTimes_naive = round(totTimes_naive);
totTimes_ext = round(totTimes_ext);

%%  Display video frames
fig =1;
for idx = 1:6
    figure(fig)
    imagesc(uNaive(:,:,:,idx));
    axis off;
    title(['Segmentation-Naive. Finished @ ',num2str(totTimes_naive(idx)),'s'],'FontSize',16);
    fig = fig+1;
    
    figure(fig)
    imagesc(uSmooth(:,:,:,idx));
    axis off;
    title(['Segmentation-Accel. Finished @ ',num2str(totTimes_ext(idx)),'s'],'FontSize',16);
    fig = fig+1;
    
end

%%  Display runnning times
fprintf('Naive Times  |  Accelerated Times')
[totTimes_naive totTimes_ext]

fprintf('Naive ObjFunction | Accelerated ObjFunction')
[ objNaive objSmooth ]