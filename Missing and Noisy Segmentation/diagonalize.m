%%  Split the Diagonals of an Image into a Cells
%
% Yang Xiao <yangshaw@umich.edu>
% 25 NOV 2016

function [ fd ] = diagonalize( f, direction )
%
%   Generate m+n-1 cells, each each corresponds to a 1-2 diagonal of f

[m,n,s] = size(f);
fd = cell(m+n-1,1);

%Updward Diagonal
if strcmp(direction,'12')
    for j = 1 : n
        for i = 1 : m
            fd{i+j-1} = [fd{i+j-1},f(i,j,:)];
        end
    end
%Downward Diagonal
elseif strcmp(direction,'21')
    f = flipud(f);
    for j = 1 : n
        for i = 1 : m
            fd{i+j-1} = [fd{i+j-1},f(i,j,:)];
        end
    end
else
    error('Direction must be 12 or 21');
end

end

