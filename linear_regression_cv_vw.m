%% Training voxel-wise encoding models % modified line 7 'floor!'
function [W, Rmat, Lambda] = linear_regression_cv_vw(Y, X, lambda, nfold)

    % validation: nfold cross validation    
    Nt = size(Y,1); % number of total time points
    Nc = size(Y,2); % number of components
    dT = floor(Nt/nfold); % number of time points in a fold
    T = Nt - dT; % number of time points for training
    Nv = size(X,2); % number of voxels
    
    Rmat = zeros(Nv,length(lambda), nfold,'single');
    
    disp('Validating regularization parameters ..');
    for nf = 1 : nfold
        idx_t = 1:Nt; % 't' means training
        idx_v = (nf-1)*dT+1:nf*dT; % 'v' means validation
        idx_t(idx_v) = [];
        
        Y_t = Y(idx_t,:);
        Y_v = Y(idx_v,:);
        YTY = Y_t'*Y_t;
        X_v = X(idx_v,:);
        X_t = X(idx_t,:);
        
        for k = 1 : length(lambda)
            lmb = lambda(k);
            disp(['fold: ', num2str(nf),' lambda: ', num2str(lmb)]);
            M = (YTY + lmb*T*eye(Nc))\Y_t'; 

            R = zeros(Nv,1);
            v1 = 1;
            while (v1 <= Nv)
                v2 = min(v1+4999,Nv);    
                W = M*X_t(:,v1:v2);
                X_vp = Y_v * W; % predicted
                X_vg = X_v(:,v1:v2); % ground truth
                R(v1:v2) = amri_sig_corr(X_vg, X_vp, 'mode', 'auto');
                v1 = v2+1;
            end
            Rmat(:,k,nf) = R;
        end
    end
    
    % choose optimal regularization parameters
    [~, Lambda] = max(mean(Rmat,3),[],2);
    Lambda = lambda(Lambda);
    
    % training with optimal regulariztion paramters
    disp('Training encoidng models with optimal regularization parameters ..');
    YTY = Y'*Y;
    W = zeros(Nc, Nv, 'single');
    
    for k = 1 : length(lambda)
        disp(['Progress: ',num2str(k/length(lambda)*100),'%']);
        lmb = lambda(k);
        M = (YTY + lmb*Nt*eye(Nc))\Y'; 
        M = (YTY + lmb*Nt*eye(Nc))\Y'; 
        voxels = (Lambda==lmb);
        W(:,voxels) = M*X(:,voxels);
     end
    
end
