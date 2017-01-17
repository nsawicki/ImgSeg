function [u,v,w,z,obj] = segmentImage(f,gamma,tau,mu,thres)
%   The essence of this project: Using ADMM to split segmentation task into
%   subproblems. The subproblems are solved with DPsubproblem.
%
%   Anish Lahiri, Daniel Vial, Nathan Sawicki, Yang Xiao

    % Init
    obj = [];

    delta = Inf;
    rotate = 0;
    [m,n,s] = size(f);
    if m == n
        disp('error - no square images');
    elseif m < n
        f = rot90(f);
        [m,n,s] = size(f);
        rotate = 1; % ensures segmented image rotated back at end
    else
        % everything works as is
    end
    
    omg_c = sqrt(2)-1;
    omg_d = 1-1/sqrt(2);

    lambda = zeros(m,n,s,6);

    % Initialize segmentation for four directions
    u = zeros(size(f)); v = f; w = f; z = f;

    % Run iterations until convergence
    while delta > thres
        % row updates
        a = (f+2*mu*(v+w+z)+2*(-lambda(:,:,:,1)-...
            lambda(:,:,:,2)-lambda(:,:,:,3)))/(1+6*mu);
        for i = 1 : m
            ai = squeeze(a(i,:,:));
            u(i,:,:) = reshape(DPsubproblem(ai,...
                4*gamma*omg_c/(1+6*mu)),[1,n,s]);
        end

        % 12 diagonal updates
        b = (f+2*mu*(u+v+z)+2*(lambda(:,:,:,2)+lambda(:,:,:,4)-...
            lambda(:,:,:,6)))/(1+6*mu);
        b12 = diagonalize(b,'12');
        wd = zeros(min(m,n),m+n-1,s);
        for k = 1 : m + n - 1
            idx1 = max(1,k+1-max(m,n)); idx2 = min(k,min(m,n));
            bk = squeeze(b12(idx1:idx2,k,:));
            if numel(bk) == s
                bk = bk';
            end
            wd(idx1:idx2,k,:) = reshape(DPsubproblem(bk,...
                4*gamma*omg_d/(1+6*mu)),[1,size(bk,1),s]);
        end
        w = rectangularize(wd,m,n,s,'12');

        % column updates
        c = (f+2*mu*(u+w+z)+2*(lambda(:,:,:,1)-lambda(:,:,:,4)-...
            lambda(:,:,:,5)))/(1+6*mu);
        for j = 1 : n
            cj = squeeze(c(:,j,:));
            v(:,j,:) = reshape(DPsubproblem(cj,...
                4*gamma*omg_c/(1+6*mu)),[m,1,s]);
        end

        % 21 diagonal updates
        d = (f+2*mu*(u+v+w)+2*(lambda(:,:,:,3)+lambda(:,:,:,5)+...
            lambda(:,:,:,6)))/(1+6*mu);
        d21 = diagonalize(d,'21');
        zd = zeros(min(m,n),m+n-1,s);
        for l = 1 : m + n - 1
            idx1 = max(1,l+1-max(m,n));
            idx2 = min(l,min(m,n));
            dl = squeeze(d21(idx1:idx2,l,:));
            if numel(dl) == s
                dl = dl';
            end
            zd(idx1:idx2,l,:) = reshape(DPsubproblem(dl,...
                4*gamma*omg_d/(1+6*mu)),[1,size(dl,1),s]); 
        end
        z = rectangularize(zd,m,n,s,'21');

        % update dual variables
        lambda(:,:,:,1) = lambda(:,:,:,1) + mu*(u-v);
        lambda(:,:,:,2) = lambda(:,:,:,2) + mu*(u-w);
        lambda(:,:,:,3) = lambda(:,:,:,3) + mu*(u-z);
        lambda(:,:,:,4) = lambda(:,:,:,4) + mu*(v-w);
        lambda(:,:,:,5) = lambda(:,:,:,5) + mu*(v-z);
        lambda(:,:,:,6) = lambda(:,:,:,6) + mu*(w-z);

        % update mu and delta
        mu = tau*mu;
        delta = norm(u(:)-v(:),'fro') + norm(u(:)-w(:),'fro') + ...
            norm(u(:)-z(:),'fro') + norm(v(:)-w(:),'fro') + ...
            norm(v(:)-z(:),'fro') + norm(w(:)-z(:),'fro');
        
        % output obj
        obj = [obj,objfun(u,v,w,z,f,gamma)];
        
    end
    
    % rotate image back if originally rotated
    if rotate
        u = rot90(u,-1);
        v = rot90(v,-1);
        w = rot90(w,-1);
        z = rot90(z,-1);
    end
    
end

function h = DPsubproblem( f , gamma )
%   Use dynamic programming to solve Potts subproblems. f can be a row,
%   column, 12 diagonal, or 21 diagonal of the input image. The detail
%   procedure is presented as Algorithm 1 in our report.
%
%   Anish Lahiri, Daniel Vial, Nathan Sawicki, Yang Xiao

    % Init
    [n,s] = size(f);
    M = cumsum(f);
    S = sum(cumsum(f.*f),2);
    h = zeros(n,s);
    P = zeros(n,1);
    J = zeros(n,1);

    % Find optimal jump locations
    for r = 1 : n
        P(r) = S(r) - norm(M(r,:))^2/r;
        J(r) = 0;
        for l = r : -1 : 2
            d = S(r) - S(l-1) - norm(M(r,:)-M(l-1,:))^2/(r-l+1);
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
        h(l+1:r,:) = repmat(mean(f(l+1:r,:),1),[r-l,1]);
        r = l;
        if r == 0
            break;
        end
        l = J(r);
    end
end

function fd = diagonalize(f,direction )
%   Convert a square image matrix into an array of its diagonals, either 12
%   direction or 21 direction
%
%   Daniel Vial, Yang Xiao

    [m,n,s] = size(f);
    fd = zeros(min(m,n),m+n-1,s);

    if strcmp(direction,'12')
        frot = rot90(f);
        for i=1:s
            fd(:,:,i) = flipud(spdiags(frot(:,:,i)));
        end
    elseif strcmp(direction,'21')
        for i=1:s
            fd(:,:,i) = spdiags(f(:,:,i));
        end
    else
        error('Direction must be 12 or 21');
    end

end

function f = rectangularize(fd,m,n,s,direction)
%   Convert an array of diagonal items back into a square image matrix.
%
%   Daniel Vial, Yang Xiao

    f = zeros(m,n,s);
    if strcmp(direction,'12')
        for i=1:s
            f(:,:,i) = spdiags(fd(:,:,i))';
        end
    elseif strcmp(direction,'21')
        for i=1:s
            f(:,:,i) = rot90(spdiags(fd(:,:,i)));
        end        
    else
        error('Direction must be 12 or 21');
    end
end


