%Function to generate group, design matrix and response vector
%Input: n (sample size), d (group size), g (number of groups), r (overlapping ratio), variance
%Output: G,X,y,beta_real

function[G,X,y,beta_real] = treedatagen(n, d, g, r, variance)

%Generate group matrix G
       
% Initialize the matrix G with zeros
    G = zeros(g, d);

    % Set the first row with 1s in the first 'd' positions
    G(1, 1:d) = 1;

    % Generate subsequent rows based on the previous row
    for i = 2:g
        % Set 1s in the current row up to 'r * num_ones'
        G(i, 1:i*d) = 1;
    end


% Generate X

    p = size(G,2);
    g = size(G,1);
    
    mu =   zeros(1,p);
    toep = zeros(p,p);
    G1 = mypar(G);
   
     for i = 1:g
         j = g - i + 1;
      idx_tmp = find(G1(j,:) == 1);
      maxtmp = max(idx_tmp);
      ltmp = length(idx_tmp);
      A = eye(ltmp);
      A(A == 0) = variance;
      toep(idx_tmp, idx_tmp) = A;
     % toep(1:maxtmp, 1:maxtmp) = 0.36^(log(j));
      toep(1:maxtmp, 1:maxtmp) = 0.36;
      toep(idx_tmp, idx_tmp) = A;
     end
 
   Sigma = toep;  
   [V, D] = eig(Sigma);
   D1 = diag(D);
   D1(D1<0.1) = 0.1;
   D = diag(D1);
   Sigma = V*D*V';
  % Make sure A_projected is symmetric
   Sigma  = (Sigma  + Sigma') / 2;

   X = mvnrnd(mu,Sigma,n);


% Generate beta

    p = size(G,2);
    g = size(G,1);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
    
    k = floor(0.9*g);

    idx_tmp = find(G(k,:) == 1);
     
    beta(idx_tmp) = 0;

 
   
    beta_real = beta;
    Expect = X*beta_real';
  
    sig = mean(abs(Expect));
  
    noi =   sig/3;

% Generate y
   y  = X* beta_real' +  normrnd(0,  noi,[1,n])';

   X = (1/sqrt(n)) * X;
   y = (1/sqrt(n)) * y;
end