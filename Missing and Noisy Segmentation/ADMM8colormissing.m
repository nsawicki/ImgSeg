%   Anish Lahiri, Daniel Vial, Nathan Sawicki, Yang Xiao
%   Missing data additions by Nathan Sawicki

function [u,obj] = ADMM8colormissing(f,gamma,tau,mu,threshold,weights,f_unmasked)

obj = [];
delta = 1;
delta_prev = 2;
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
    
    if(abs(delta_prev - delta) <= .0001)
        weights = ones(size(weights));
    end
    delta_prev = delta;
    
    subplot(1,4,1)
    imshow(u)
    subplot(1,4,2)
    imshow(v)
    subplot(1,4,3)
    imshow(w)
    subplot(1,4,4)
    imshow(z)
    drawnow;
    
    % u update
    a = (f+2*mu*(v+w+z)+2*(-lambda1-lambda2-lambda3))/(1+6*mu);
    for i = 1 : m
        ai = squeeze(a(i,:,:));
        wh = reshape(weights(i,:),n,1);
        u(i,:,:) = reshape(DPsubproblem_colMiss(ai,4*gamma*omg_c/(1+6*mu),wh),[1,n,s]); % transpose weights?
    end
    
    % w update
    b = (f+2*mu*(u+v+z)+2*(lambda2+lambda4-lambda6))/(1+6*mu);
    b12 = diagonalize(b,'12');
    w12 = diagonalize(weights,'12');
    for k = 1 : m + n - 1
        if numel(b12{k}) == 3
            weight_diag = reshape(squeeze(w12{k}),size(squeeze(b12{k}),2),1);
            wd{k} = reshape(DPsubproblem_colMiss(squeeze(b12{k})',4*gamma*omg_d/(1+6*mu),weight_diag),[1,1,3]);
        else
            weight_diag = reshape(squeeze(w12{k}),size(squeeze(b12{k}),1),1);
            wd{k} = reshape(DPsubproblem_colMiss(squeeze(b12{k}),...
                4*gamma*omg_d/(1+6*mu),weight_diag),1,[],3);
        end
    end
    w = rectangularize(wd,m,n,s,'12');
    
    % v update
    c = (f+2*mu*(u+w+z)+2*(lambda1-lambda4-lambda5))/(1+6*mu);
    for j = 1 : n
        cj = squeeze(c(:,j,:));
        wv = reshape(weights(:,j),m,1);
        v(:,j,:) = reshape(DPsubproblem_colMiss(cj,4*gamma*omg_c/(1+6*mu),wv),[m,1,s]); % transpose weights?
    end
    
    % z update
    d = (f+2*mu*(u+v+w)+2*(lambda3+lambda5+lambda6))/(1+6*mu);
  
   
    d21 = diagonalize(d,'21');
    w21 = diagonalize(weights,'21');
    for l = 1 : m + n - 1
        if numel(d21{l}) == 3
            weight_diag = reshape(squeeze(w21{l}),size(squeeze(d21{l}),2),1);
            zd{l} = reshape(DPsubproblem_colMiss(squeeze(d21{l})',4*gamma*omg_d/(1+6*mu),weight_diag),[1,1,3]);
        else
            weight_diag = reshape(squeeze(w21{l}),size(squeeze(d21{l}),1),1);
            zd{l} = reshape(DPsubproblem_colMiss(squeeze(d21{l}),...
                4*gamma*omg_d/(1+6*mu),weight_diag),1,[],3);        
        end
    end
    z = rectangularize(zd,m,n,s,'21');
%      figure(9)
%      subplot(1,2,1)
%     imshow(z)
%     subplot(1,2,2)
%     imshow(d);
%     drawnow
%       pause(5)
    
    % update dual variables
    lambda1 = lambda1 + mu*(u-v);
    lambda2 = lambda2 + mu*(u-w);
    lambda3 = lambda3 + mu*(u-z);
    lambda4 = lambda4 + mu*(v-w);
    lambda5 = lambda5 + mu*(v-z);
    lambda6 = lambda6 + mu*(w-z);
    
    mu = tau*mu;
    
    %Compute Objective function for this iteration
    obj = [obj objfun(u,v,w,z,f,gamma)];
    
    % Compute norm difference between consensus variables
    delta = norm(u(:)-v(:),'fro')+norm(u(:)-w(:),'fro')+norm(u(:)-z(:),'fro')...
        +norm(v(:)-w(:),'fro')+norm(v(:)-z(:),'fro')+norm(w(:)-z(:),'fro')
    
%     if(max(isnan(delta)))
%         figure()
%         subplot(1,4,1)
%         imshow(u)
%         subplot(1,4,2)
%         imshow(v)
%         subplot(1,4,3)
%         imshow(w)
%         subplot(1,4,4)
%         imshow(z)
%         drawnow;
%         delta = delta_prev - .01
%     end
    
end

figure()
plot(obj)
