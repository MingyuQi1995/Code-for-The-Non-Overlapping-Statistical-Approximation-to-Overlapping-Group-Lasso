%Function to use a line search to find the regularization path for overlapping group lasso
%Additionally, this function generates on a log scale similar to glm package in R.
%Input: X,y, nlambda (the number of elements in the path), opts (see SLEP: Sparse Learning with Efficient Projections)
%Output: A sequence of lambda 


function [lambda_seq] = lambda_gen_log_real(X,y,nlambda,opts,orderedVector)

lambda_max_init = 1e10;
lambda_min_init =  1e-10;
beta1 = 0.9;
beta2 = 1.1;
%----------------------- Set optional items -----------------------

orderedVector = orderedVector;




flag  = 1;
while(flag)
    
   z = lambda_max_init;
   [beta_est,~, ~, ~]=  glLogisticR(X, y, z, opts);
 %  beta_est_ori = accumarray(orderedVector', beta_est');
   if(sum(abs(beta_est) > 1e-50) < 1e-20)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end


flag  = 1;
while(flag)
    z =  lambda_min_init;
   [beta_est, c,~, ~]=  glLogisticR(X, y, z, opts);

   beta_est_ori = accumarray(orderedVector', beta_est');
   
   if(sum(abs(beta_est_ori) < 1e-50) < 1e-20)
       lambda_min_init = lambda_min_init*beta2;
   else
       flag = 0;
   end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);


end