%%  segmentExtracted.m
%   To test the performance of difference extraction function.
%   
%   Daniel Vial
clear; clc;

%%  Init
numFrames = 6;
extractThres = 0.005;
tau = 2; admmThres = 10;
gamma = 0.05; mu = 0.01*gamma;

tmp = im2double(imread('man/1.tiff'));
[m,n,s] = size(tmp);
f = zeros(m,n,s,numFrames);
% Load frames
f(:,:,:,1) = tmp(:,:,1:s);
for i=2:10
    tmp = im2double(imread(['man/' num2str(i) '.tiff']));
    f(:,:,:,i) = tmp(:,:,1:s);
end

%%  Extract difference rectangles
[fExtract,boundary] = extractSimilarVideo(f,extractThres);

uExtract = cell(size(fExtract));
extractTimes = zeros(size(fExtract));
for i=1:length(uExtract)
    disp(['running image ' num2str(i)]);
    tic;
    uExtract{i} = segmentImage(fExtract{i},gamma,tau,mu,admmThres);
    extractTimes(i) = toc;
    save('tmp.mat','fExtract','boundary','uExtract','extractTimes');
end

%%  Draw boundaries of these difference rectangles
for i = 1 : numFrames
    figure(i);
    
    subplot(2,1,1); hold on;
    imagesc(f(:,:,:,i));
    line([boundary(i,3),boundary(i,3)],[boundary(i,1),boundary(i,2)],...
        'color','r','linewidth',3);
    line([boundary(i,4),boundary(i,4)],[boundary(i,1),boundary(i,2)],...
        'color','r','linewidth',3);
    line([boundary(i,3),boundary(i,4)],[boundary(i,1),boundary(i,1)],...
        'color','r','linewidth',3);
    line([boundary(i,3),boundary(i,4)],[boundary(i,2),boundary(i,2)],...
        'color','r','linewidth',3);
    
    subplot(2,1,2); hold on;
    imagesc(f(:,:,:,i+1));
    line([boundary(i,3),boundary(i,3)],[boundary(i,1),boundary(i,2)],...
        'color','r','linewidth',3);
    line([boundary(i,4),boundary(i,4)],[boundary(i,1),boundary(i,2)],...
        'color','r','linewidth',3);
    line([boundary(i,3),boundary(i,4)],[boundary(i,1),boundary(i,1)],...
        'color','r','linewidth',3);
    line([boundary(i,3),boundary(i,4)],[boundary(i,2),boundary(i,2)],...
        'color','r','linewidth',3);
end