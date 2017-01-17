%   Anish Lahiri, Daniel Vial, Nathan Sawicki, Yang Xiao
%   Missing data additions by Nathan Sawicki

function [u] = ADMM8color(f,gamma,tau,mu,threshold,wh,wv)


objective = 0;
delta = 1;
[m,n,s] = size(f);

%Set compass and diagonal weights
omg_c = sqrt(2) - 1;
omg_d = 1 - 1/sqrt(2);

%Init Lambda
lambda1 = zeros(m,n,s); lambda2 = zeros(m,n,s); lambda3 = zeros(m,n,s);
lambda4 = zeros(m,n,s); lambda5 = zeros(m,n,s); lambda6 = zeros(m,n,s);

%Initialize u,v,w,z = f
u = f; v = f; w = f; z = f;

wd = cell(m+n-1,1); zd = cell(m+n-1,1);

%break when the norm difference of u,v,w,z are lower than the threshold
while delta > threshold
    imshow(u)
    drawnow;
    
    % u update
    a = (f+2*mu*(v+w+z)+2*(-lambda1-lambda2-lambda3))/(1+6*mu);
    for i = 1 : m
        ai = squeeze(a(i,:,:));
        u(i,:,:) = reshape(DPsubproblem_col(ai,4*gamma*omg_c/(1+6*mu),wv'),[1,n,s]); % transpose weights?
    end
    
    
    % w update
    b = (f+2*mu*(u+v+z)+2*(lambda2+lambda4-lambda6))/(1+6*mu);
    b12 = diagonalize(b,'12');
    for k = 1 : m + n - 1
        if numel(b12{k}) == 3
            wd{k} = reshape(DPsubproblem_col(squeeze(b12{k})',4*gamma*omg_d/(1+6*mu),ones(size(squeeze(b12{k}),2),1)),[1,1,3]);
        else
            wd{k} = reshape(DPsubproblem_col(squeeze(b12{k}),...
                4*gamma*omg_d/(1+6*mu),ones(size(squeeze(b12{k}),1),1)),1,[],3);
        end
    end
    w = rectangularize(wd,m,n,s,'12');
    
    % v update
    c = (f+2*mu*(u+w+z)+2*(lambda1-lambda4-lambda5))/(1+6*mu);
    for j = 1 : n
        cj = squeeze(c(:,j,:));
        v(:,j,:) = reshape(DPsubproblem_col(cj,4*gamma*omg_c/(1+6*mu),wh),[m,1,s]); % transpose weights?
    end
    
    % z update
    d = (f+2*mu*(u+v+w)+2*(lambda3+lambda5+lambda6))/(1+6*mu);
    d21 = diagonalize(d,'21');
    for l = 1 : m + n - 1
        if numel(d21{l}) == 3
            zd{l} = reshape(DPsubproblem_col(squeeze(d21{l})',4*gamma*omg_d/(1+6*mu),ones(size(squeeze(d21{l}),2),1)),[1,1,3]);
        else
            zd{l} = reshape(DPsubproblem_col(squeeze(d21{l}),...
                4*gamma*omg_d/(1+6*mu),ones(size(squeeze(d21{l}),1),1)),1,[],3);        
        end
    end
    z = rectangularize(zd,m,n,s,'21');
    
    % update dual variables
    lambda1 = lambda1 + mu*(u-v);
    lambda2 = lambda2 + mu*(u-w);
    lambda3 = lambda3 + mu*(u-z);
    lambda4 = lambda4 + mu*(v-w);
    lambda5 = lambda5 + mu*(v-z);
    lambda6 = lambda6 + mu*(w-z);
    
    
    % Compute Objective Function
    objective = [objective objfun(u,v,w,z,f,gamma)];
    
    % Compute norm difference of consensus variables
    mu = tau*mu;
    delta = norm(u(:)-v(:),'fro')+norm(u(:)-w(:),'fro')+norm(u(:)-z(:),'fro')...
        +norm(v(:)-w(:),'fro')+norm(v(:)-z(:),'fro')+norm(w(:)-z(:),'fro')
    
end

figure()
plot(objective(2:end))