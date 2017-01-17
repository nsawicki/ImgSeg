%%  Concatenate the Diagonals of an Image
%
% Yang Xiao <yangshaw@umich.edu>
% 25 NOV 2016

function [ f ] = rectangularize( fd , m , n , s, direction )
%
%   Convert cells of diagonal elements into m by n matrix

f = zeros(m,n,s);

%Updward Diagonal
if strcmp(direction,'12')
    for j = n : -1 : 1
        for i = m : -1 : 1
            f(i,j,:) = fd{i+j-1}(1,end,:);
            fd{i+j-1}(:,end,:) = [];
        end
    end
    
%Downward Diagonal
elseif strcmp(direction,'21')
    for j = n : -1 : 1
        for i = m : -1 : 1
            f(i,j,:) = fd{i+j-1}(1,end,:);
            fd{i+j-1}(:,end,:) = [];
        end
    end
    f = flipud(f);
else
    error('Direction must be 12 or 21');
end

end

