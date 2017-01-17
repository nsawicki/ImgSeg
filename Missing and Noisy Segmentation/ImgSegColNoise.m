%%  Image Partitioning using ADMM
%
% Yang Xiao <yangshaw@umich.edu>
% Adapted by Nathan Sawicki <nsawicki@umich.edu>
% 12 Dec 2016

%%  Init
clear all; close all;
tic;

%load image
imageRGB = im2double(imread('pebbles4.jpg'));
f = imageRGB;


[m,n,s] = size(f);

%Set Parameters
tau = 2;
gamma = .5;
wh = ones(m,1); % weights are 1 when data is not missing
wv = ones(1,n);

%Apply Gaussian Noise to f
f_noise = imnoise(f,'gaussian',0,.1)


figure()
imagesc(f_noise);

%Set parameter Mu
mu = 0.01*gamma;



%%  ADMM
[u] = ADMM8color(f_noise,gamma,tau,mu,10^-8,wh,wv);

figure(1)
subplot(1,2,1)
imshow(f_noise);
title('Original')
subplot(1,2,2)
imshow(u);
title('Partitioned')
toc;