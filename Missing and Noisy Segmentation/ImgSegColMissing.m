%%  Image Partitioning using ADMM
%
% Yang Xiao <yangshaw@umich.edu>
% 25 NOV 2016

%%  Init
clear all; close all;
tic;
%load Image
imageRGB = im2double(imread('pebbles4.jpg'));
f = imageRGB;
[m,n,s] = size(f);

%Set Parameters tau, gamma
tau = 2;
gamma = .5;
wh = ones(m,1);
wv = ones(1,n);

%Reshape the image for processing
ftmp = reshape(f,m*n,3);

%Set how much of the image is kept/Masked
total_samps = m*n;
fraction_keep = .7;
numsamps_kept = floor(total_samps.*fraction_keep);
samps_keep = randperm(m*n,numsamps_kept);

%Set weights to zero for pixels that are missing or corrupted
f_out = zeros(m*n,3);
weights = zeros(m*n,1);
f_out(samps_keep,:) = ftmp(samps_keep,:);
weights(samps_keep) = 1;

%Reshape back to the original dimension of f
f_masked = reshape(f_out,[size(f)]);
weights = reshape(weights,m,n);

figure()
imagesc(f_masked);
figure()
imagesc(weights)



mu = 0.01*gamma;



%%  ADMM
[u,objectiveFuncion] = ADMM8colormissing(f_masked,gamma,tau,mu,10^-8,weights,f);


figure()
imshow(f_masked);
title('Original')
figure()
imshow(u);
title('Partitioned')
toc;