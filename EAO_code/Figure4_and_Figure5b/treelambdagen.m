%Function to use a line search to find the regularization path.
%Input: X,y, W0,set1,set2, nlambda (the number of elements in the path),  package, type
%Output: lambda_seq (a sequence of lambda)



function[lambda_seq] = treelambdagen(data,y,W0,set1,set2,nlambda, package, type)
           
    lambda_max_init = 1e8;
    lambda_min_init =  1e-8;
    beta1 = 0.9;
    beta2 = 1.1;
    

    if ~ismember(type, {'0', '1', '2'})
        error('Invalid type. Type must be "0", "1", or "2".');
    end

    % Check which package is specified and perform actions accordingly
    if strcmp(package, 'slep')
       
       if strcmp(type, '0') 

       B = data;
          
       opts = set2;
        
       opts.method = 0;
       opts.mFlag=0;      
       opts.lFlag=0;  
        
       opts.tFlag = 3;
        
       opts.tol = 1e-4;
        
        flag  = 1;
        
           while(flag)
              z = lambda_max_init;
                [beta_est, ~, ~]= tree_LeastR(B, y, z, opts);
              if(sum(abs(beta_est) > 1e-8) < 1e-9)
                lambda_max_init = lambda_max_init*beta1;
              else
              flag = 0;
             end
           end

          opts.init=2;  

          z = lambda_min_init;
          [beta_est, ~, ~]= tree_LeastR(B, y, z, opts);
          lambda_min_init = lambda_min_init*beta2;
          flag  = 1;

           while(flag)
            opts.x0 = beta_est;
            opts.init=1;   
         
            z = lambda_min_init;
            %disp(z);
        
            [beta_est, ~, ~]= tree_LeastR(B, y, z, opts);

             if(sum(abs(beta_est) < 1e-8) < 1e-9)
                lambda_min_init = lambda_min_init*beta2;
             else
                  flag = 0;
             end
            end

         C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
         lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);

      elseif strcmp(type, '1')

       Xwlas = data;
       opts = set2;
       opts.tFlag = 3;
       opts.tol = 1e-4;
       
       flag  = 1;
        
       while(flag)
         rho = lambda_max_init;

        [beta_est, ~, ~]=  LeastR(Xwlas, y, rho, opts);
       if(sum(abs(beta_est) > 1e-8) < 1e-9)
         lambda_max_init = lambda_max_init*beta1;
       else
        flag = 0;
      end
     end


   rho = lambda_min_init;
   [beta_est, ~, ~]=  LeastR(Xwlas, y, rho, opts);
   lambda_min_init = lambda_min_init*beta2;

     flag  = 1;

   while(flag)
       opts.x0 = beta_est;
     opts.init =1;
     rho = lambda_min_init;
    % disp(rho);
     [beta_est, ~, ~]=  LeastR(Xwlas, y, rho, opts);
    if(sum(abs(beta_est) < 1e-8) < 1e-9)
               lambda_min_init = lambda_min_init*beta2;
    else
          flag = 0;
    end
   end
   
 C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
 lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);


        elseif strcmp(type, '2')

       Xwlas = data;
       opts = set2;
       opts.tFlag = 3;
       opts.tol = 1e-4;
       flag  = 1;
        
    while(flag)
        rho = lambda_max_init;
     [beta_est, ~, ~]=  glLeastR(Xwlas, y, rho, opts);
      if(sum(abs(beta_est) > 1e-8) < 1e-9)
        lambda_max_init = lambda_max_init*beta1;
       else
        flag = 0;
      end
    end

    rho = lambda_min_init;
     [beta_est, ~, ~]=  glLeastR(Xwlas, y, rho, opts);
       lambda_min_init = lambda_min_init*beta2;

     flag  = 1;
   while(flag)
       opts.x0 = beta_est;
       opts.init = 1;
     rho = lambda_min_init;
     [beta_est, ~, ~]=  glLeastR(Xwlas, y, rho, opts);
    if(sum(abs(beta_est) < 1e-8) < 1e-9)
               lambda_min_init = lambda_min_init*beta2;
       else
          flag = 0;
       end
    end
    
    
 C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
 lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);

 fval =1;

        end


    elseif strcmp(package, 'spam')
       
       if strcmp(type, '0')

            B = data;

            tree = set1;

            param = set2;

            param.tol=1e-3;

            flag  = 1;

         while(flag)

            param.lambda = lambda_max_init;

           [beta_est, ~]=mexFistaTree(y,B,W0,tree,param);


          if(sum(abs(beta_est) > 1e-8) < 1e-9)

            lambda_max_init = lambda_max_init*beta1;

          else

             flag = 0;

          end

        end
     
       param.lambda = lambda_min_init;

       [beta_est,~]=mexFistaTree(y,B,W0,tree,param);

        lambda_min_init = lambda_min_init*beta2;

      flag  = 1;


       while(flag)

            W0 = beta_est;
          param.lambda = lambda_min_init;

         [beta_est,~]=mexFistaTree(y,B,W0,tree,param);
     
       if(sum(abs(beta_est) < 1e-8) < 1e-9)

               lambda_min_init = lambda_min_init*beta2;
       else

          flag = 0;

       end

      end
   

C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);

lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);

% disp("lambda_seq generated")

 elseif strcmp(type, '1')

      Xwlas = data;

      param = set2;
   
      flag  = 1;
        
    while(flag)
        param.lambda = lambda_max_init;

       [beta_est, ~]=mexLasso(y,Xwlas,param);

        beta_est = full(beta_est);

      if(sum(abs(beta_est) > 1e-8) < 1e-9)
        lambda_max_init = lambda_max_init*beta1;
       else
        flag = 0;
      end
    end

     flag  = 1;

   while(flag)
     param.lambda = lambda_min_init;
    %disp(lambda_min_init);
       [beta_est, ~]=mexLasso(y,Xwlas,param);
        beta_est = full(beta_est);

       if (sum(abs(beta_est) < 1e-8) < 1e-9)
               lambda_min_init = lambda_min_init*beta2;
       else
          flag = 0;
       end
    end
    
 % lambda_min_init = max(lambda_min_init, 0.1);
    
  

 C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);

 lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);

 % disp("lambda_seq generated")


      elseif strcmp(type, '2')

      Xwlas = data;

      paramg = set2;

      paramg.tol=1e-4;

      flag  = 1;
        
    while(flag)

        paramg.lambda = lambda_max_init;


       [beta_est, ~]=mexFistaFlat(y,Xwlas,W0,paramg);

      if(sum(abs(beta_est) > 1e-8) < 1e-9)

        lambda_max_init = lambda_max_init*beta1;

      else

        flag = 0;

      end
    end

    lambda_min_init = 1;
      paramg.lambda = lambda_min_init;
      disp(paramg.lambda );
      [beta_est, ~]=mexFistaFlat(y,Xwlas,W0,paramg);
       lambda_min_init = lambda_min_init*beta2;

     flag  = 1;
     
   while(flag)
      W0 = beta_est;
     paramg.lambda = lambda_min_init;
    disp(paramg.lambda );
     
      [beta_est, ~]=mexFistaFlat(y,Xwlas,W0,paramg);


       if sum(abs(beta_est) < 1e-8) < 1e-9
               lambda_min_init = lambda_min_init*beta2;
       else
          flag = 0;
       end
    end

 C = (log(lambda_max_init) - log(lambda_min_init))/(nlambda - 1);
 

 lambda_seq = lambda_min_init*exp(((1:nlambda) - 1)*C);

  disp("lambda_seq generated")
        end
    else
        error('Invalid package. Supported packages are "slep" and "spam".');
    end
end