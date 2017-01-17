function [fExtract,boundary] = extractSimilarVideo(f,thres)
%   Output the extacted rectangles and corresponding boundaries for all
%   frames
%   
%   Daniel Vial

    fExtract = cell(size(f,4)-1,1);
    boundary = zeros(size(f,4)-1,4);
    for i=1:length(fExtract)
        [r1,r2,c1,c2,fExtract{i}] = xtract_similar(f(:,:,:,i),f(:,:,:,i+1),thres);
        boundary(i,:) = [r1,r2,c1,c2];        
    end

end


function [ r1, r2, c1, c2, ximg ] = xtract_similar( img1, img2, thres )
%   Compare two images and output a rectangle containing their difference 
%   and the horizontal/vertical boundary values. 
%
%   Yang Xiao, Anish Lahiri

    [m,n,~] = size(img1);

    diffRGB = abs(img1 - img2).^3;
    diff = sum(diffRGB,3)*100;
    
    % Horizontal distribution
    hDist = sum(diff,1);
    hDist = hDist/sum(hDist);
    % Vertical distribution
    vDist = sum(diff,2);
    vDist = vDist/sum(vDist);
    
    % The significance of the difference is determined by the threshold
    hDist(hDist < thres) = 0;
    vDist(vDist < thres) = 0;
    hDist(hDist>0)=1;
    vDist(vDist>0)=1;

    % Output boundaries
    c1 = find(cumsum(hDist)==1);
    c2 = n-find(cumsum(fliplr(hDist))==1);
    r1 = find(cumsum(vDist)==1);
    r2 = m-find(cumsum(flipud(vDist))==1);
    r1=r1(1);
    r2=r2(1);
    c1=c1(1);
    c2=c2(1);

    ximg = img2(r1:r2,c1:c2,:);
end

