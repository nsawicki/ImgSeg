%%  gammaSelect.m
%
%   To generate and compare segmentation results with different gammas.
%
%   Daniel Vial, Yang Xiao, Anish Lahiri

clear; clc;

%%  Init
numFrames = 6;
gamma = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5];
tau = 2; 
admmThres = 10;

tmp = im2double(imread('man/1.tiff'));
[m,n,s] = size(tmp);
f = zeros(m,n,s,numFrames);
f(:,:,:,1) = tmp(:,:,1:s);

%%  Read images
for i = 2 : numFrames
    tmp = im2double(imread(['man/' num2str(i) '.tiff']));
    f(:,:,:,i) = tmp(:,:,1:s);
end

%%  Run segmentation
naiveTimes = zeros(size(f,4),length(gamma));
uNaive = zeros([size(f),length(gamma)]);
for i=1:length(gamma)
    for j=1:size(f,4)
        disp(['running gamma = ' num2str(gamma(i)) ', image = ' num2str(j)]);
        tic;
        uNaive(:,:,:,j,i) = segmentImage(f(:,:,:,j),gamma(i),tau,0.01*gamma(i),admmThres);
        naiveTimes(j,i) = toc;    
    end
end

%% Save results
save('tmpGammaSel.mat','uNaive','naiveTimes');




