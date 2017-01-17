%%  Compute the Potts Objective Function
%
% Nathan Sawicki <nsawicki@umich.edu>
% 12 Dec 2016

%% 
function [cost] = objfun(u,v,w,z,f,gamma)

omg_c = sqrt(2) - 1;
omg_d = 1 - sqrt(2)/2;


[m,n,s] = size(f);

%count the number of row-wise color jumps
u_count = 0;
for i = 1:n
    pix_prev = u(1,i,:);
   for j = 2:m
       if(norm(squeeze(pix_prev(1,1,1:s) - u(j,i,1:s)),'fro') > 10^-8)
           u_count = u_count + 1;
       end
       pix_prev = u(j,i,:);
   end
end

%count the number of column-wise color jumps
v_count = 0;
for i = 1:m
    pix_prev = v(i,1,:);
   for j = 2:n
       if(norm(squeeze(pix_prev - v(i,j,:)),'fro') > 10^-8)
           v_count = v_count + 1;
       end
       pix_prev = v(i,j,:);
   end
end

wd = diagonalize(w,'12');

%count the number of upper-diagonal-wise color jumps
w_count = 0;
for i = 1:length(wd)
    w_curr = wd{i};
    pix_prev = w_curr(1,1,1:s);
    for j = 2:size(w_curr,2)
        if(norm(squeeze(pix_prev - w_curr(1,j,:)),'fro') > 10^-8)
           w_count = w_count + 1;
        end
        pix_prev = w_curr(1,j,1:s);
    end
end

%count the number of downward-diagonal-wise color jumps
z_count = 0;
zd = diagonalize(z,'21');
for i = 1:length(zd)
    z_curr = zd{i};
    pix_prev = z_curr(1,1,1:s);
    for j = 2:size(z_curr,2)
        if(norm(squeeze(pix_prev - z_curr(1,j,:)),'fro') > 10^-8)
           z_count = z_count + 1;
        end
        pix_prev = z_curr(1,j,:);
    end
end
       

%compute the total Potts cost
cost = gamma.*omg_c.*(u_count + v_count) + gamma.*omg_d.*(w_count + z_count)...
        + .25.* (  norm(u(:) - f(:),'fro').^2 + norm(v(:) - f(:),'fro').^2 + norm(w(:) - f(:),'fro').^2 + norm(z(:) - f(:),'fro').^2)