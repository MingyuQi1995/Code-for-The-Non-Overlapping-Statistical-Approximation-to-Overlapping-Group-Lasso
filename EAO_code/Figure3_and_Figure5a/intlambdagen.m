%Function to use line search to find the regularization path for overlapping group lasso
%Additionally, this function generates on a log scale similar to glm package in R.
%Input: X,y, nlambda (the number of elements in the path), opts (see SLEP: Sparse Learning with Efficient Projections),
%type ("0" for overlapping group lasso,"1" for lasso, "2" for proposed method).
%Output: lambda_seq (a sequence of lambda ).


function [lambda_seq] = intlambdagen(X,y,nlambda,opts,type)

lambda_max_init =  1e8;

lambda_min_init =  1e-8;
beta1 = 0.9;
beta2 = 1.1;


    if strcmp(type, '0')
 
       flag  = 1;

       while(flag)
       z = [0, lambda_max_init];
       [beta_est, ~, ~]= overlapping_LeastR(X, y, z, opts);
   if(sum(abs(beta_est) > 1e-8) < 1e-10)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end

  flag  = 1;
   opts.init = 1;
   opts.x0 = beta_est;
   while(flag)
     z = [0, lambda_min_init];
     [beta_est, ~, ~]= overlapping_LeastR(X, y, z, opts);
     if(sum(abs(beta_est) < 1e-8) < 1e-10)
        lambda_min_init = lambda_min_init*beta2;
     else
       flag = 0;
     end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);


    elseif strcmp(type, '1')
    
       flag  = 1;
     while(flag)
    
       rho = lambda_max_init;
       [beta_est, ~, ~]= LeastR(X, y, rho, opts);
      if(sum(abs(beta_est) > 1e-8) < 1e-10)
        lambda_max_init = lambda_max_init*beta1;
       else
        flag = 0;
      end
    
     end
     
% Change:del
lambda_min_init =  1e-5;

flag  = 1;

   while(flag)
      opts.init = 1;
      opts.x0 = beta_est;
      rho = lambda_min_init*beta2;
     [beta_est, ~, ~]= LeastR(X, y, rho, opts);
     if(sum(abs(beta_est) < 1e-8) < 1e-10)
        lambda_min_init = lambda_min_init*beta2;
     else
       flag = 0;
     end
   end

 C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
 lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);
 
    elseif strcmp(type, '2')
    
      flag  = 1;
      while(flag)
    
   rho = lambda_max_init;
   [beta_est, ~, ~]= glLeastR(X, y, rho, opts);
   if(sum(abs(beta_est) > 1e-8) < 1e-9)
       lambda_max_init = lambda_max_init*beta1;
   else
       flag = 0;
   end
    
end
flag  = 1;
while(flag)
    rho =  lambda_min_init;
    opts.init = 1;
      opts.x0 = beta_est;
    [beta_est, ~, ~]= glLeastR(X, y, rho, opts);
     if(sum(abs(beta_est) < 1e-8) < 1e-10)
        lambda_min_init = lambda_min_init*beta2;
     else
       flag = 0;
     end
end

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);
    else
        % Handle the case when 'type' is not "0", "1", or "2"
        disp('Invalid type. No action performed.');
    end

end