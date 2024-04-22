%Funtion of cross-validation for overlapping group lasso
%Input: X,y,lambda_seq (regularization path), opts, beta_real (the true coefficient vector),type
%Output: Lambda_min (the lambda with the samllest MSE), lambda_min_pattern (the lambda with the samllest support discrepancy)



function [lambda_min,lambda_min_pattern] = realcv(X,y,lambda_seq, opts, beta_real, type)

beta_realpattern = abs(beta_real);

threshold = 1e-10;

beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
nlambda = length(lambda_seq);

Delta = [];
Deltapattern = [];
  
    if strcmp(type, '0')

         opts.init=2;   
     
       opts.maxIter2 = 1e+5;
       opts.flag2 = 2 ;

       z = [0, lambda_seq(nlambda)];

       [beta_est, ~, ~]= overlapping_LeastR(X, y, z, opts);
 
       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(positions);
 
    for j = 1:(nlambda-1)
        opts.init = 1;
        opts.x0 = beta_est;
        z = [0, lambda_seq(nlambda - j)];
        [beta_est, ~, ~]= overlapping_LeastR(X, y, z, opts);
   
       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(positions);

    end
    Delta_min = min(Delta);
    idx_min = find(Delta == Delta_min,1,"last");
    lambda_min= lambda_seq(idx_min);
 

    Delta_minpattern = min(Deltapattern);
    idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
    lambda_min_pattern = lambda_seq(idx_minpattern );
    

    elseif strcmp(type, '1')

       opts.init=2;   

       rho = lambda_seq(nlambda);
       [beta_est, ~, ~]= LeastR(X, y, rho, opts);
       beta_estpattern = abs(beta_est);
 
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(positions);
 

 
    for j = 1:(nlambda-1)
        opts.init = 1;
        opts.x0 = beta_est;
        rho = lambda_seq(nlambda - j);
        [beta_est, ~, ~]= LeastR(X, y, rho, opts);
        Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
        beta_estpattern = abs(beta_est);
        beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
        beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
        Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(positions);
 
    end


   Delta_min = min(Delta);
   idx_min = find(Delta == Delta_min,1,"last");
   lambda_min= lambda_seq(idx_min);
   Delta_minpattern = min(Deltapattern);
   idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
   lambda_min_pattern = lambda_seq(idx_minpattern );
   
    elseif strcmp(type, '2')

        
        z = lambda_seq(nlambda);


  [beta_est, ~, ~]= glLeastR(X, y, z, opts);
  beta_estpattern = abs(beta_est);
  beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
  beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
  Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
  Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(positions);


 
    for j = 1:(nlambda-1)
        opts.init = 1;
        opts.x0 = beta_est;
        z =  lambda_seq(nlambda - j);
      
        [beta_est, ~, ~]= glLeastR(X, y, z, opts);
  
       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       beta_estpattern = abs(beta_est);
      beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
      beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(positions);
 
    end

Delta_min = min(Delta);
idx_min = find(Delta == Delta_min,1,"last");
lambda_min= lambda_seq(idx_min);

 Delta_minpattern = min(Deltapattern);
idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
lambda_min_pattern = lambda_seq(idx_minpattern );
    else
        % Handle the case when 'type' is not "0", "1", or "2"
        disp('Invalid type. No action performed.');
    end
    
        opts.init = 2;   
end