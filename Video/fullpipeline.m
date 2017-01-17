%%  fullpipeline.m
%
%   Main program for video segmentation. Output segmented video frames and
%   performance mesures.
%
%   Daniel Vial

clear; clc;

%%  Init
% parameters
gamma = 0.1; mu = 0.01*gamma; mult = 500; tau = 2; 
admmThres = 10; 
extractThres = 0.0002;
numFrames = 6;

% objective functions
objNaive = zeros(1,numFrames);
objSmooth = zeros(1,numFrames);

% load images
tmp = im2double(imread('man/1.tiff'));
[m,n,s] = size(tmp);
f = zeros(m,n,s,numFrames);
f(:,:,:,1) = tmp(:,:,1:s);
for i = 2 : numFrames
    tmp = im2double(imread(['man/' num2str(i) '.tiff']));
    f(:,:,:,i) = tmp(:,:,1:s);
end

%%  Naive method
naiveTimes = zeros(size(f,4),1);
uNaive = zeros(size(f));
vNaive = zeros(size(f));
wNaive = zeros(size(f));
zNaive = zeros(size(f));

% Run segmentation for every frame successively
for i=1:size(f,4)
    disp(['naive method ' num2str(i)]);
    tic;
    [uNaive(:,:,:,i),vNaive(:,:,:,i),wNaive(:,:,:,i),zNaive(:,:,:,i)] ...
        = segmentImage(f(:,:,:,i),gamma,tau,mu,admmThres);
    naiveTimes(i) = toc;
    % obj
    objNaive(i) = objfun(uNaive(:,:,:,i),vNaive(:,:,:,i),...
        wNaive(:,:,:,i),zNaive(:,:,:,i),f(:,:,:,i),gamma);
end


%%  Accelerated method
% Extract rectangles
[fExtract,boundary] = extractSimilarVideo(f,extractThres);

% segment extracted
uExtract = cell(size(fExtract));
extractTimes = zeros(size(fExtract));
for i=1:length(uExtract)
    disp(['extracted ' num2str(i)]);
    tic;
    uExtract{i} = segmentImage(fExtract{i},gamma,tau,mu,admmThres);
    extractTimes(i) = toc;
end

% combine and smooth images
uSmooth = zeros(size(uNaive));
vSmooth = zeros(size(vNaive));
wSmooth = zeros(size(wNaive));
zSmooth = zeros(size(zNaive));
uCombined = zeros(size(uNaive));

smoothTimes = zeros(size(f,4),1);
uSmooth(:,:,:,1) = uNaive(:,:,:,1);
vSmooth(:,:,:,1) = vNaive(:,:,:,1);
wSmooth(:,:,:,1) = wNaive(:,:,:,1);
zSmooth(:,:,:,1) = zNaive(:,:,:,1);

objSmooth(1) = objNaive(1);

% Run segmentation for every frame successively
for i = 2 : size(uSmooth,4)
    disp(['smoothing ' num2str(i)]);
    uCombined(:,:,:,i) = uSmooth(:,:,:,i-1);
    uCombined(boundary(i-1,1):boundary(i-1,2),boundary(i-1,3):boundary(i-1,4),:,i) = uExtract{i-1};
    tic;
    [uSmooth(:,:,:,i),vSmooth(:,:,:,i),wSmooth(:,:,:,i),zSmooth(:,:,:,i)] ...
        = segmentImage(uCombined(:,:,:,i),gamma,tau,mu*mult,admmThres);
    smoothTimes(i) = toc;
    % compute objective function
    objSmooth(i) = objfun(uSmooth(:,:,:,i),vSmooth(:,:,:,i),...
        wSmooth(:,:,:,i),zSmooth(:,:,:,i),f(:,:,:,i),gamma);
end

%% Save results
save(['man_fullpipe_gamma=',num2str(gamma),'_xThresh=',num2str(extractThres),'_mult=',num2str(mult),'.mat'],...
    'uNaive','naiveTimes','fExtract','boundary',...
    'uExtract','extractTimes','uSmooth','uCombined',...
    'smoothTimes','objNaive','objSmooth');