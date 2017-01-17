%%  dispGammaSelimgs.m
%   Display segmentation results for different gammas
%
%   Anish Larihi

load('tmpGammaSel.mat');

gamma = [0.0100    0.0200    0.0500    0.1000    0.2000    0.3000    0.4000    0.5000];
for idx =1 : length(gamma)
    for sbidx = 1:6
        imagesc(uNaive(:,:,:,sbidx,idx));
        title(['gamma = ',num2str(gamma(idx))]);
        pause
        disp 'hit key'
    end
end