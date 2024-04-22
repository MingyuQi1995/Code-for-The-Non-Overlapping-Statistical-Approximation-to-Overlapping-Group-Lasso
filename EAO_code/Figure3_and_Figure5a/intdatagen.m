%Function to generate group, design matrix and response vector
%Input: n (sample size), d (group size), g (number of groups), r (overlapping ratio), variance
%Output: G, X, y,beta_real

function[G,X,y,beta_real] = intdatagen(n, d, g, r, variance)

    overlap_num = floor(d*r);
    p0 = d*g;
    G = zeros(g, p0);

    G(1,1:d) =  1;

    k = d - overlap_num + 1;

    for i = 2:g
        k1 = k + d - 1;
        G(i, k:k1) = 1;
        k = k1 - overlap_num + 1 ;
    end

    G = G(:,1:k1);
  
    disp("G generated")



    toep_cor = variance;

    p = size(G,2);
    g = size(G,1);
    
    mu =   zeros(1,p);
    toep = zeros(p,p);

    G1 = mypar(G);

    k = size(G1,1);
     for i = 1:g
     idx_tmp = find(G(i,:) == 1);
    toep(idx_tmp, idx_tmp) = toep_cor^2;
     end
 
    for i = 1:k
    idx_tmp = find(G1(i,:) == 1);
    toep(idx_tmp, idx_tmp) = toep_cor;
    end
     
    clear G1;
     
    Sigma = toep;  
   
    clear toep;
   
    Sigma = Sigma + (1 - toep_cor ) *diag(repelem(1 , p));
    [V, D] = eig(Sigma);
    D1 = diag(D);
    D1(D1<0.1) = 0.1;
    D = diag(D1);
    Sigma = V*D*V';

    Sigma  = (Sigma  + Sigma') / 2;
   
    disp("Simga generated")
   X = mvnrnd(mu,Sigma,n);

  disp("X generated")


    p = size(G,2);
    g = size(G,1);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
    
    k = floor(0.9*g);
     
    zerogroup = randsample(1:g,k,false);
    
 
     for i = 1:length(zerogroup)
        idx_tmp = find(G(zerogroup(i),:) == 1);
        beta(idx_tmp) = 0;
     end
 
   
    beta_real = beta;
    Expect = X*beta_real';
  
    sig = mean(abs(Expect));
  
    noi =   sig/3;
   disp("beta generated")
% Generate y
   y  = X* beta_real' +  normrnd(0,  noi,[1,n])';

   X = (1/sqrt(n)) * X;
   y = (1/sqrt(n)) * y;
end