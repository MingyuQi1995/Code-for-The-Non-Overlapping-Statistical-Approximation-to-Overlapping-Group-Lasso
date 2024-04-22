%Funtion of cross-validation for overlapping group lasso
%Input: X,y,lambda_seq (regularization path), opts, beta_real (the true coefficient vector)
%Output: lambda_LL (the lambda with the samllest loss),  lambda_auc (the lambda with the best loss). Other three are unused.





function [lambda_LL, lambda_ber, lambda_auc, cv_max,cv_min,auc_max] = cvreal(X, y, lambda_seq, opts, K,foldinf)
%cross validation

%if isempty(nlambda)
%    nlambda = 100;
%end

if isempty(K)
    K = 5;
end

[n,p] = size(X);
%if n<=p
%    lambda_min_ratio = 0.01;
%else
%    lambda_min_ratio = 0.0001;
%end

%lambda_min = lambda_max * lambda_min_ratio;

%C = (log(lambda_max) - log(lambda_min))/(nlambda - 1);

%lambda_seq = lambda_min*exp(((1:nlambda) - 1)*C);

nlambda = length(lambda_seq);


cv_LL = zeros(K, nlambda);
cv_ber = zeros(K, nlambda);
cv_auc = zeros(1, nlambda);
pred_y = zeros(n,nlambda);
y_test0 = zeros(n,nlambda);

for i = 1:K
    
    
    test_idx = foldinf(i,:);
    X_train = X;
    X_train(test_idx,:) = [];
    X_test = X(test_idx,:);
    y_train = y;
    y_train(test_idx,:) = [];
    y_test = y(test_idx,:);
    
for j = 1:nlambda
        
        z = [0, lambda_seq(j)];
        [beta_est, c, ~, ~]= overlapping_LogisticR(X_train, y_train, z, opts);
        y_prob = 1./(1+exp(-X_test*beta_est - c));
        
        %deviance holdout likelihood
        y_test1 = 0.5*y_test + 0.5;
        y_test0(test_idx, j) = y_test1;
        pred_y(test_idx,j) = y_prob;
        y_LL = y_test1.*log(y_prob) + (1- y_test1).*log(1 - y_prob);
        
        cv_LL(i,j) = sum(y_LL);
        
        conf1 = round(y_prob) + y_test1;
        d = length(find(conf1 == 2));
        a = length(find(conf1 == 0));
        conf2 = round(y_prob) - y_test1;
        b = length(find(conf2 == 1));
        m = length(find(conf2 == -1));



        blanceacc_ggl = 0.5*(b/(a+b) + m/(m+d));
        cv_ber(i,j) = blanceacc_ggl;
    end
    
end
cv_m = mean(cv_LL,1);
cv_max = max(cv_m);
idx_max = find(cv_m == cv_max,1,'last');
lambda_LL = lambda_seq(idx_max);

cv_berm = mean(cv_ber,1);
cv_min = min(cv_berm);
idx_min = find(cv_berm == cv_min,1,'last');
lambda_ber = lambda_seq(idx_min);

for k = 1:nlambda
    [~,~,~,auc_tmp] = perfcurve(y_test0(:,k), pred_y(:,k), 1);
    cv_auc(k) = auc_tmp;
end

auc_max = max(cv_auc);
idx_auc = find(cv_auc == auc_max,1,'last');
lambda_auc = lambda_seq(idx_auc);

end






