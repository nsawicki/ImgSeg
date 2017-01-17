function [ obj ] = objfun( u,v,w,z,f, gamma )
%	Compute the objective function of ADMM8 Potts problem
%   
%   Anish Lahiri, Nathan Sawicki, Yang Xiao

[m,n,~] = size(f);

wc = sqrt(2) - 1;
wd = 1 - sqrt(2)/2;

sumu=0;
sumv=0;
sumw=0;
sumz=0;

for i = 1 : m
    for j = 1 : n
        if i <= m-1 && norm(squeeze(u(i,j,:))-squeeze(u(i+1,j,:))) == 0
            sumu = sumu + 1;
        end
        if j <= n-1 && norm(squeeze(v(i,j,:))-squeeze(v(i,j+1,:))) == 0
            sumv = sumv + 1;
        end
        if i > 1 && j <= n-1 && norm(squeeze(w(i,j,:))-squeeze(w(i-1,j+1,:))) == 0
            sumw = sumw + 1;
        end
        if i <= m-1 && j <= n-1 && norm(squeeze(z(i,j,:))-squeeze(z(i+1,j+1,:))) == 0
            sumz = sumz + 1;
        end
    end
end

% output
obj =  gamma*wc*(sumu+sumv) + gamma*wd*(sumw+sumz) ...
    + 0.25*(norm(u(:)-f(:),'fro')^2+norm(v(:)-f(:),'fro')^2 ...
    + norm(w(:)-f(:),'fro')^2+norm(z(:)-f(:),'fro')^2);

end
