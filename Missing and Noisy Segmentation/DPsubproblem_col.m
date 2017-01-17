%   Anish Lahiri, Daniel Vial, Nathan Sawicki, Yang Xiao
%   Missing data additions by Nathan Sawicki

function h = DPsubproblem_col( f , gamma,  w )
%
%   Accelerated Dynamic Programming for the Univariate Potts Problem

% Init
[n,s] = size(f);
h = zeros(n,s);
P = zeros(n,1);
M = zeros(n,s);
S = zeros(n,1);
W = zeros(n,1);
J = zeros(n,1);

% Find the optimal jump locations
% for r = 1
r = 1; 
M(r,:) = 0 + w(r)*f(r,:);
S(r) = 0 + w(r)*norm(f(r,:))^2;
W(r) = 0 + w(r);
P(r) = S(r) - norm(M(r,:))^2/W(r);
J(r) = 0;
%{
for l = r : 2
    d = S(r)- M(r)^2/W(r);
    if P(r) < d - gamma
        break;
    end
    p = gamma + d;
    if p <= P(r)
        P(r) = p;
        J(r) = l - 1;
    end     
end
%}
% for r = 2 : n
for r = 2 : n
    M(r,:) = M(r-1,:) + w(r)*f(r,:);
    S(r) = S(r-1) + w(r)*norm(f(r,:))^2;
    W(r) = W(r-1) + w(r);
    P(r) = S(r) - norm(M(r,:))^2/W(r);
    J(r) = 0;
    for l = r : -1 : 2
        d = S(r) - S(l-1) - norm(M(r,:)-M(l-1,:))^2/(W(r)-W(l-1));
        if P(r) < d - gamma
            break;
        end
        p = P(l-1) + gamma + d;
        if p <= P(r)
            P(r) = p;
            J(r) = l - 1;
        end
        
    end
end

% Reconstruct the minimizer h from the optimal jump locations
r = n;
l = J(r);
while r > 0

    h(l+1:r,:) = repmat((w(l+1:r)')*f(l+1:r,:)/sum(w(l+1:r)),[r-l,1]);
    r = l;
    if r == 0
        break;
%         keyboard;
    end
    l = J(r);
end

end

