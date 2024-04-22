%Funtion of cross-validation
%Input: X,y,W0, set1, set2, lambda_seq (regularization path), beta_real (the true coefficient vector), package, type
%When package = slep, W0 is useless, set1 is the function values, set2 is the opts, type "0" is overlapping group lasso type "1" is lasso, and type "2" is proposed group lasso.
%When package = spam, W0 is the intial value, set1 is useless, set2 is param, type "0" is overlapping group lasso, type "1" is lasso, and type "2" is proposed group lasso.


%Output: fval (function value), lambda_min (the lambda with the samllest estimation error),lambda_min_pattern (the lambda with the samllest support discrepancy)




function[fval,lambda_min,lambda_min_pattern] = treecv(data,y,W0,set1,set2,lambda_seq,beta_real, package, type)

      beta_realpattern = abs(beta_real);
      threshold = 1e-10;
      beta_realpattern(beta_realpattern < threshold) = 0;
      beta_realpattern(beta_realpattern >= threshold) = 1;
      nlambda = length(lambda_seq);
    

    if strcmp(package, 'slep')
   
       if strcmp(type, '0')

       X =data;

       fval = set1;

       opts = set2;

       Delta = [];

       Deltapattern = [];
 
       
       opts.init=2;   
    
       opts.tFlag = 2;
       
       opts.tol = fval(end);
       
        
       z = lambda_seq(nlambda);

       [beta_est, ~, ~]= tree_LeastR(X, y, z, opts);

       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);
       
       nRows = size(X,1); 
      

    for j = 1:(nlambda-1)
        
         opts.x0 = beta_est;

         opts.init=1;   
         
         opts.tFlag = 2;

         opts.tol = fval(nlambda - j);
         
    
         
        z = lambda_seq(nlambda - j);

       [beta_est, ~, ~]= tree_LeastR(X, y, z, opts);

       % z = [0,lambda_seq(nlambda - j)];
       % [beta_est, ~, ~] = overlapping_LeastR(X, y, z, opts);

  
       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end


    Delta_min = min(Delta);

    idx_min = find(Delta == Delta_min,1,"last");

    lambda_min= lambda_seq(idx_min);
 

    Delta_minpattern = min(Deltapattern);
    idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
    lambda_min_pattern = lambda_seq(idx_minpattern );
    
 
      fval = 1;

      elseif strcmp(type, '3')

       A = data;

       fval = set1;

       opts = set2;

       Delta = [];

       Deltapattern = [];
 
       opts.init=2;   
    
       opts.tFlag = 2;
       
       opts.tol = fval(end);
       
       opts.method = 1;
  
       z = [0,lambda_seq(nlambda)];
       
       [beta_est, ~, ~] = overlapping_LeastR(A, y, z, opts);

       beta_estpattern = abs(beta_est);
       
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);
       
        
      n  = size(A,1); % Number of rows in X
      

    for j = 1:(nlambda-1)

        
         opts.init=1;  

         opts.x0 = beta_est;
         
         opts.tFlag = 2;

         opts.maxIter2= 3000;
        
        
         opts.tol = 1.0001*fval(nlambda - j);
         
  
        
     

         opts.flag2 = 2;
         
         opts.tol2 = fval(nlambda - j)/10;

 
      
         z = [0,lambda_seq(nlambda - j)];
         
       
        [beta_est, ~, ~] = overlapping_LeastR(A, y, z, opts);

  
        Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
        beta_estpattern = abs(beta_est);
        beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
        beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
        Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end


    Delta_min = min(Delta);

    idx_min = find(Delta == Delta_min,1,"last");

    lambda_min= lambda_seq(idx_min);
 

    Delta_minpattern = min(Deltapattern);
    idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
    lambda_min_pattern = lambda_seq(idx_minpattern );
    
 
      fval = 1;

      elseif strcmp(type, '1')

       Xwlas = data;
       fval = set1;
       opts = set2;
       Delta = [];
       Deltapattern = [];
       
       opts.init=2;   
       opts.tFlag = 2;
       opts.tol = fval(end);
       rho = lambda_seq(nlambda);
       [beta_est, ~, ~]=  LeastR(Xwlas, y, rho, opts);
       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

          
      nRows = size(Xwlas,1); % Number of rows in X

    for j = 1:(nlambda-1)
    
        opts.tFlag = 2;
        opts.x0 = beta_est;  
        opts.init=1;   
        
      opts.tol = 1.0001*fval(nlambda - j);
       
       
               
        
   
        rho = lambda_seq(nlambda - j);
        
        

       [beta_est, ~, ~]=  LeastR(Xwlas, y, rho, opts);
       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       beta_estpattern = abs(beta_est);
       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
       Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end

    fval = 1;
    Delta_min = min(Delta);
    idx_min = find(Delta == Delta_min,1,"last");
    lambda_min= lambda_seq(idx_min);
 

    Delta_minpattern = min(Deltapattern);
    idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
    lambda_min_pattern = lambda_seq(idx_minpattern );
    
       elseif strcmp(type, '2')
           
       Xwlas = data;
       fval = set1;
       opts = set2;
       Delta = [];
       Deltapattern = [];
       
       opts.init=2;   
       opts.tFlag = 2;
       
       opts.tol = fval(end);

       rho = lambda_seq(nlambda);

       [beta_est, ~, ~]=  glLeastR(Xwlas, y, rho, opts);

       beta_estpattern = abs(beta_est);

       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;

       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;

       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));

       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);
 
    

    for j = 1:(nlambda-1)

    
         opts.x0 = beta_est;  
         opts.init=1;   
         opts.tFlag = 2;

         opts.tol = 1.0001*fval(nlambda - j);
           

        rho = lambda_seq(nlambda - j);

        [beta_est, ~, ~]=  glLeastR(Xwlas, y, rho, opts);

       Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));

       beta_estpattern = abs(beta_est);

       beta_estpattern (abs(beta_estpattern)  < threshold) = 0;

       beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;

       Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end
    Delta_min = min(Delta);
    idx_min = find(Delta == Delta_min,1,"last");
    lambda_min= lambda_seq(idx_min);
 

    Delta_minpattern = min(Deltapattern);
    idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
    lambda_min_pattern = lambda_seq(idx_minpattern );
    
    fval = 1;
       end


    elseif strcmp(package, 'spam')


        if strcmp(type, '0')

        B = data;

        tree = set1;

        param = set2;
        
        nRows = size(B,1); 
        param.tol= param.tol * nRows;
        
        fval = [];
        Delta = [];
        Deltapattern = [];

        param.lambda = lambda_seq(nlambda);

        [beta_est,fv]=mexFistaTree(y,B,W0,tree,param);

        beta_est = flip(beta_est);

        fval(nlambda) = fv(1);

        beta_estpattern = abs(beta_est);
        beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
        beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
        Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
        Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    for j = 1:(nlambda-1)

        param.lambda = lambda_seq(nlambda - j);
        
        W0 = beta_est;
        
       [beta_est,fv]=mexFistaTree(y,B,W0,tree,param);
   
       fval(nlambda - j) =  fv(1);

      beta_est = flip(beta_est);

      Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));

      beta_estpattern = abs(beta_est);

      beta_estpattern (abs(beta_estpattern)  < threshold) = 0;

      beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;

      Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end

Delta_min = min(Delta);
idx_min = find(Delta == Delta_min,1,"last");
lambda_min= lambda_seq(idx_min);

 Delta_minpattern = min(Deltapattern);
idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
lambda_min_pattern = lambda_seq(idx_minpattern );


   elseif strcmp(type, '1')

        X = data;
        param = set2;
        param.lambda2 = 0;
        param.mode = 2;
        nRows = size(X,1); 
        k =  sqrt(nRows);
        param.tol= param.tol * nRows;
        fval = [];
        Delta = [];
        Deltapattern = [];

        param.lambda = lambda_seq(nlambda);
        [beta_est, ~] = mexLasso(y,X,param);
        beta_est = full(beta_est);
        r1 = y - X*beta_est;
       fval(nlambda) = 0.5*(sum(r1.^2)) + param.lambda*sum(abs(beta_est));

       beta_estpattern = abs(beta_est);
       beta_estpattern (beta_estpattern  < threshold) = 0;
       beta_estpattern (beta_estpattern  >= threshold) = 1;
       Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
       Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);
 
       for j = 1:(nlambda-1)
       
        param.lambda = lambda_seq(nlambda - j);
       [beta_est, ~] = mexLasso(y,X,param);
        beta_est = full(beta_est);
        r1 = y - X*beta_est;

        fval(nlambda - j) = 0.5*(sum(r1.^2)) + param.lambda*sum(abs(beta_est));
       
      Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
      beta_estpattern = abs(beta_est);
      beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
      beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
      Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

       end

Delta_min = min(Delta);
idx_min = find(Delta == Delta_min,1,"last");
lambda_min= lambda_seq(idx_min);

 Delta_minpattern = min(Deltapattern);
idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
lambda_min_pattern = lambda_seq(idx_minpattern );



        elseif strcmp(type, '2')

        D = data;
      
      
        paramg = set2;
        nRows = size(D,1); 
        paramg.tol= paramg.tol*nRows;
        fval = [];

        Delta = [];

        Deltapattern = [];

        paramg.lambda = lambda_seq(nlambda);
        [beta_est, fv]=mexFistaFlat(y,D,W0,paramg);

        fval(nlambda) = fv(1);

  beta_estpattern = abs(beta_est);
  beta_estpattern (beta_estpattern  < threshold) = 0;
  beta_estpattern (beta_estpattern  >= threshold) = 1;
  Delta(nlambda) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
  Deltapattern(nlambda) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);
 
    for j = 1:(nlambda-1)

       % disp(j);
       
       W0 = beta_est;

        paramg.lambda = lambda_seq(nlambda - j);
       [beta_est, fv]=mexFistaFlat(y,D,W0,paramg);

      fval(nlambda - j) = fv(1);
       
      
    
      Delta(nlambda - j) =  sum((beta_est - beta_real).^2)./(sum(beta_real.^2));
      beta_estpattern = abs(beta_est);
      beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
      beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
      Deltapattern(nlambda-j) = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

    end

Delta_min = min(Delta);
idx_min = find(Delta == Delta_min,1,"last");
lambda_min= lambda_seq(idx_min);

 Delta_minpattern = min(Deltapattern);
idx_minpattern = find(Deltapattern == Delta_minpattern ,1,"first");
lambda_min_pattern = lambda_seq(idx_minpattern );

         
        end
    else
        error('Invalid package. Supported packages are "slep" and "spam".');
    end
end