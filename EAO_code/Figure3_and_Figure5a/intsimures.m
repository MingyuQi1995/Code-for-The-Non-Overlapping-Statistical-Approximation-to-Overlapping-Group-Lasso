%Main function for simulation
%Input: n (sample size), d (group size), g (number of groups), r (overlapping ratio)
%Output: Res (Regularization path computing time and estimation error of each method under the setting).



function[Res] = intsimures(n, d, g, r)

% Change to point six

     [G,X,y,beta_real] = intdatagen(n,d,g,r,0.6);

     disp("data generated")
     
     beta_real = beta_real';
     
     root = pwd;
     
     root1 = strcat(root, "/SLEP-master");
     
     addpath(genpath("./SLEP-master"));
     
     root2 = strcat(root, "/spams-matlab-v2.6 2");
     
     addpath(genpath("./spams-matlab-v2.6 2"));
  
     nlambda = 50;

%----------------------- overlapping group lasso in SLEP-----------------------

opts=[];
opts.init=2; 
opts.method = 1;
opts.tFlag= 0;
opts.maxIter=3e+5;
opts.tol= 1e-5; 
%change : add tol2

opts.rFlag=0;    
opts.nFlag=0;   
[opt_G,opt_ind] = ginf(G,"0");
opts.G = opt_G;
opts.ind = opt_ind;
opts.rFlag=0;  
opts.rsL2=0;


lambda_seq = intlambdagen(X,y,nlambda,opts, "0");

tic;

[lambda_min,lambda_min_pattern] = intcv(X, y, lambda_seq, opts, beta_real, "0");


z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

solve_time_og = toc;

solve_time_og = solve_time_og;
threshold = 1e-20;



RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_ogp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_realpattern);

OG_res = [solve_time_og;ham_og;RSSog];

disp("og done")
%----------------------- wlas -----------------------

opts=[];
opts.init=2; 
opts.mFlag=1;      
opts.lFlag=1;      
opts.tFlag= 0;
opts.maxIter=3e+5;
opts.tol= 1e-5; 
opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

orw = sqrt(length( find(G(1,:) == 1)));
lasweight = sum(G) * orw;

beta_realwlas = beta_real.* lasweight';


Xwlas = X ./repmat(lasweight, size(X, 1), 1);

lambda_seq = intlambdagen(X,y,nlambda,opts, "1");

tic;


[lambda_min,lambda_min_pattern] = intcv(Xwlas, y, lambda_seq, opts, beta_realwlas, "1");



rho =lambda_min;

[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


solve_time_wlas = toc;

solve_time_wlas = solve_time_wlas;

beta_est_wlas = beta_est_wlas./lasweight';
beta_realwlas = beta_realwlas ./ lasweight';

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));



beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_realpattern);



wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];
disp("las done")

%----------------------- wglasso -----------------------


[Group,w] = ginf(G,"2");

vec = Group;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 

opts=[];

opts.init=2;        % starting from a zero point

opts.method = 1;
opts.tFlag= 0;
opts.maxIter=3e+5;


opts.tol= 1e-5; 

opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 
opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w'; % set the weight for the group, a cloumn vector     %
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search



lambda_seqwgl = intlambdagen(X,y,nlambda,opts,"2");

tic;

[lambda_min,lambda_min_pattern] = intcv(X, y, lambda_seqwgl, opts, beta_real, "2");



z = lambda_min;

[beta_est_wgl,~, ~]= glLeastR(X, y, z, opts);

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

solve_time_wgl = toc;

solve_time_wgl = solve_time_wgl;

RSS_wgl = sqrt(sum((beta_est_wgl -  beta_real).^2))./sqrt((sum( beta_real.^2)));


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wgl = hamm/length(beta_realpattern);


wgl_res = [solve_time_wgl;ham_wgl;RSS_wgl];
disp("wgl done")


%----------------------- wglasso1 -----------------------


[Group,w] = ginf(G,"3");

vec = Group;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 

opts=[];

opts.init=2;        % starting from a zero point

opts.method = 1;
opts.tFlag= 0;
opts.maxIter=3e+5;
opts.tol= 1e-5; 

opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 
opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w'; % set the weight for the group, a cloumn vector     %
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search


lambda_seqwgl = intlambdagen(X,y,nlambda,opts,"2");

tic;


[lambda_min,lambda_min_pattern] = intcv(X, y, lambda_seqwgl, opts, beta_real, "2");



z = lambda_min;

[beta_est_wgl,~, ~]= glLeastR(X, y, z, opts);

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

solve_time_wgl = toc;

solve_time_wgl = solve_time_wgl;

RSS_wgl = sqrt(sum((beta_est_wgl -  beta_real).^2))./sqrt((sum( beta_real.^2)));


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wgl = hamm/length(beta_realpattern);


wgl1_res = [solve_time_wgl;ham_wgl;RSS_wgl];
disp("wgl1 done")


%----------------------- wglassod-----------------------


[Group,w] = ginf(G,"4");

vec = Group;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 

opts=[];

opts.init=2;        % starting from a zero point

opts.method = 1;
opts.tFlag= 0;
opts.maxIter=3e+5;
opts.tol= 1e-5; 

opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 
opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w'; % set the weight for the group, a cloumn vector     %
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search


lambda_seqwgl = intlambdagen(X,y,nlambda,opts,"2");

tic;


[lambda_min,lambda_min_pattern] = intcv(X, y, lambda_seqwgl, opts, beta_real, "2");



z = lambda_min;

[beta_est_wgl,~, ~]= glLeastR(X, y, z, opts);

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

solve_time_wgl = toc;

solve_time_wgl = solve_time_wgl;

RSS_wgl = sqrt(sum((beta_est_wgl -  beta_real).^2))./sqrt((sum( beta_real.^2)));


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wgl = hamm/length(beta_realpattern);


wgld_res = [solve_time_wgl;ham_wgl;RSS_wgl];
disp("wgld done")



%% Final

Res = [OG_res wlas_res wgl_res wgl1_res wgld_res];

end