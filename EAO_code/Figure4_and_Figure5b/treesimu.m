%Main function for simulation
%Input: n (sample size)
%Output: Res (Regularization path computing time, estimation error, and support discrepancy of each method under the setting).





function[Res] = treesimu(n)
     n = n;
     d = 4;
     g = 800;
     r = 1;
     var = 0.6;
     
     %%%change point
     
     [G,X,y,beta_real] = treedatagen(n,d,g,r,var);

     beta_real = beta_real';
     root = pwd; 
     root1 = strcat(root, "/SLEP-master");
     addpath(genpath("./SLEP-master"));
     root2 = strcat(root, "/spams-matlab-v2.6");
     addpath(genpath("./spams-matlab-v2.6"));
  
     nlambda = 50; 

     beta_realpattern = abs(beta_real);
     threshold = 1e-10;
     beta_realpattern(abs(beta_realpattern) < threshold) = 0;
     beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
     positions = find(abs(beta_realpattern) > threshold);

    %----------------------- lambda sequence overlapping group lasso -----------------------
    opts=[];

    opts.init=2;        % starting from a zero point

    opts.maxIter=3e+5;  % maximum number of iterations  

    opts.rFlag=0;       % use ratio
    
    opts.nFlag=0;       % without normalization

    %%%%% change point
    [opt_G,opt_ind] = tree_OG_ind1(G);

    opts.G = opt_G;

    opts.ind = opt_ind;

    opts.rsL2=0;


    lambda_seqt = treelambdagen(X,y,1,1,opts,nlambda,"slep","0");


%----------------------- overlapping group lasso in SPAM time-----------------------

B = zeros(size(X));
p =  size(X,2);
% Reverse the order of columns
for i = 1:p
    B(:, i) = X(:, p - i + 1);
end



param.num_threads=-1; % all cores (-1 by default)
param.verbose=false;   % verbosity, false by default
param.it0=50;      % frequency for duality gap computations
param.max_it=3e+5; % maximum number of iterations
param.L0=0.1;
param.tol=1e-5;
param.intercept=false;
param.pos=false;
W0=zeros(size(X,2),size(y,2));
param.loss='square';
param.regul='tree-l2';
param.admm=true;
%%%change point
[tree.own_variables,tree.N_own_variables,tree.eta_g, tree.groups] = treespamcvt_OG_ind1(G);




tic

[fval,lambda_mint,lambda_min_patternt] = treecv(B,y,W0,tree,param,lambda_seqt,beta_real,"spam","0");

solve_time_spam_og = toc;


%----------------------- overlapping group lasso in SLEP time-----------------------

tic;

[~,~,~] = treecv(X,y,W0,fval,opts,lambda_seqt,beta_real,"slep","3");


solve_time_og = toc;


%----------------------- overlapping group lasso  performance-----------------------


%[~, lambda_mint,lambda_min_patternt] = treecv1(X,y,1,1,opts,lambda_seqt,beta_real,"slep","0");


param.lambda = lambda_mint;

W0=zeros(size(X,2),size(y,2));

[beta_spam_og, fv1]=mexFistaTree(y,B,W0,tree,param);

beta_spamog = flip(beta_spam_og);

param.lambda = lambda_min_patternt;

[beta_spam_ogp, fv2]=mexFistaTree(y,B,W0,tree,param);

beta_spamogp = flip(beta_spam_ogp);

f1 = fv1(1);

f2 = fv2(1);

RSSspamog = sqrt(sum((beta_spamog - beta_real ).^2))./sqrt((sum(beta_real.^2)));

beta_estpattern = abs(beta_spamogp);

beta_estpattern (abs(beta_estpattern)  < threshold) = 0;

beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;

hamspam_og = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

OGspam_res = [solve_time_spam_og;hamspam_og;RSSspamog ];

opts.init=2;        % starting from a zero point

opts.tFlag = 2;

opts.tol = f1;

z = lambda_mint;

[beta_slep_og, ~, ~] = tree_LeastR(X, y, z, opts);


opts.tFlag = 2;

opts.tol =  f2;

z = lambda_min_patternt;
[beta_slep_ogp, ~, ~] = tree_LeastR(X, y, z, opts);


RSSog = sqrt(sum((beta_slep_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

beta_estpattern = abs(beta_slep_ogp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
ham_og = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);


OGslep_res = [solve_time_og;ham_og;RSSog];

disp("og done")


%----------------------- lambda sequence lasso -----------------------

   opts.ind = opt_ind;
   weight = opts.ind(3,:);
   Gpar = mypar(G);
   g = size(Gpar,1);
   lasweight = zeros(1,p);
   for  i = 1:g

      tmp = find(Gpar(i,:) == 1);

      lasweight(tmp) = sum(weight);

      weight = weight(2:end);
    end

   Xwlas =  X ./repmat(lasweight, size(X, 1), 1);

   beta_reallas = beta_real .*lasweight';

   opts=[];
   opts.method = 0;
   opts.mFlag=0;      
   opts.lFlag=0;        
   opts.maxIter=3e+5;
   opts.rFlag=0;    
   opts.nFlag=0;   
   opts.rsL2=0;


   lambda_seqlast = treelambdagen(Xwlas,y,1,0,opts,nlambda,"slep","1");

%----------------------- weighted lasso in SPAM time-----------------------

param = [];


param.num_threads=-1; % all cores (-1 by default)
param.verbose=false;   % verbosity, false by default
param.it0=50;      % frequency for duality gap computations
param.max_it=3e+5; % maximum number of iterations
param.L0=0.1;
param.tol=1e-5;
param.intercept=false;
param.pos=false;
param.mode = 2;
param.lambda2 = 0;
W0=zeros(size(X,2),size(y,2));
param.admm=true;
param.loss='square';

param.regul='l1';

tic;

[fvallas,lambda_mint,lambda_min_patternt] = treecv(Xwlas,y,W0,lasweight,param,lambda_seqlast,beta_reallas,"spam","1");

solve_time_spam_wlas = toc;


%----------------------- weighted lasso in SLEP time-----------------------

tic;


[~,~,~] =  treecv(Xwlas,y,W0,fvallas,opts,lambda_seqlast,beta_reallas,"slep","1");

solve_time_slep_wlas = toc;

%----------------------- weighted lasso performance-----------------------

%[~,lambda_mint,lambda_min_patternt] =  treecv1(Xwlas,y,1,1,opts,lambda_seqlast,beta_reallas,"slep","1");


param.lambda = lambda_mint  ;

[beta_spam_las, ~]=mexLasso(y,Xwlas,param);
beta_spam_las = full(beta_spam_las);


param.lambda = lambda_min_patternt;
[beta_spam_lasp, ~]=mexLasso(y,Xwlas,param);
beta_spam_lasp = full(beta_spam_lasp);


r1 = y - Xwlas*beta_spam_las;
f1 = 0.5*(sum(r1.^2)) + lambda_mint*sum(abs(beta_spam_las));
r2= y - Xwlas*beta_spam_lasp;
f2 = 0.5*(sum(r2.^2)) + lambda_min_patternt*sum(abs(beta_spam_lasp));

beta_spam_las = beta_spam_las./lasweight';
beta_reallas = beta_reallas./lasweight';
RSSspamlas = sqrt(sum((beta_spam_las - beta_reallas ).^2))./sqrt((sum(beta_reallas.^2)));

beta_estpattern = abs(beta_spam_lasp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamspam_las = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

lasspam_res = [solve_time_spam_wlas;hamspam_las;RSSspamlas ];

rho =lambda_mint;
opts.tFlag = 2; 
opts.tol = f1 ;
[beta_slep_las, ~, ~]=  LeastR(Xwlas, y, rho, opts);

rhop = lambda_min_patternt;
opts.tFlag = 2;
opts.tol = f2 ;
[beta_slep_lasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_slep_las = beta_slep_las./lasweight';
RSSsleplas = sqrt(sum((beta_slep_las - beta_reallas ).^2))./sqrt((sum(beta_reallas.^2)));

beta_estpattern = abs(beta_slep_lasp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamslep_las = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

lasslep_res = [solve_time_slep_wlas;hamslep_las;RSSsleplas ];


disp("las done")


%----------------------- lambdaseq wgl-----------------------
Gpar = mypar(G);
rownum = size(Gpar,1);
row_indices = 1:rownum ;

A = Gpar.* row_indices';


paramg.groups=int32(sum(A, 1));

vec = paramg.groups;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 
k=length(ind)-1;     % number of groups
w =ones(k,1);

opts=[];

opts.init=2;        % starting from a zero point


opts.method = 1;
opts.mFlag=1;       % treating it as compositive function 
opts.lFlag=1;   
opts.maxIter=3e+5;

opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 
opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w; % set the weight for the group, a cloumn vector  


lambda_seqwglt = treelambdagen(Xwlas,y,1,0,opts,nlambda,"slep","2");

%----------------------- weighted group lasso in SPAM time-----------------------




paramg = [];
paramg.num_threads=-1; % all cores (-1 by default)
paramg.verbose=false;   % verbosity, false by default
paramg.it0=50;      % frequency for duality gap computations
paramg.max_it=2e+5; % maximum number of iterations
paramg.L0= 0.1;
paramg.tol=1e-5;
paramg.intercept=false;
paramg.pos=false;
W0=zeros(size(X,2),size(y,2));
paramg.loss='square';
paramg.regul='group-lasso-l2';
paramg.groups=int32(sum(A, 1));
paramg.admm=true;

beta_realwgl = beta_real .*lasweight';

tic;


[fvalwgl,lambda_mint,lambda_min_patternt] = treecv(Xwlas,y,W0,lasweight,paramg,lambda_seqwglt,beta_realwgl,"spam","2");


solve_time_spam_wgl = toc;


%----------------------- weighted group lasso in SLEP time-----------------------

tic;


[~,~,~] =  treecv(Xwlas,y,W0,fvalwgl,opts,lambda_seqwglt,beta_realwgl,"slep","2");


solve_time_slep_wgl = toc;

%----------------------- weighted group lasso performance-----------------------

%[~,lambda_mint,lambda_min_patternt] =  treecv1(Xwlas,y,1,1,opts,lambda_seqwglt,beta_realwgl,"slep","2");


paramg.lambda = lambda_mint;


[beta_spam_wgl, fv1]=mexFistaFlat(y,Xwlas,W0,paramg);

paramg.lambda = lambda_min_patternt;

[beta_spam_wglp, fv2]=mexFistaFlat(y,Xwlas,W0,paramg);


fg1 = fv1(1);
fg2 = fv2(1);


beta_spam_wgl = beta_spam_wgl./lasweight';

beta_realwgl = beta_realwgl./lasweight';
RSSspamwgl = sqrt(sum((beta_spam_wgl - beta_realwgl ).^2))./sqrt((sum(beta_realwgl.^2)));

beta_estpattern = abs(beta_spam_wglp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamspam_wgl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

wglspam_res = [solve_time_spam_wgl;hamspam_wgl;RSSspamwgl ];



opts.tFlag = 2;
opts.tol = fg1 ;

z = lambda_mint;
[beta_slep_wgl,~, ~]= glLeastR(Xwlas, y, z, opts);


opts.tFlag = 2;
opts.tol =  fg2;


zp = lambda_min_patternt;

[beta_slep_wglp, ~, ~]= glLeastR(Xwlas, y, zp, opts);

beta_slep_wgl = beta_slep_wgl./lasweight';

RSSslepwgl = sqrt(sum((beta_slep_wgl - beta_realwgl ).^2))./sqrt((sum(beta_realwgl .^2)));

beta_estpattern = abs(beta_slep_wglp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamslep_wgl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

wgl_slep_res = [solve_time_slep_wgl;hamslep_wgl;RSSslepwgl ];


disp("wgl slep done")

%----------------------- lamdseq-----------------------

rownum = size(Gpar,1);
row_indices = 1:rownum ;
A = Gpar.* row_indices';
paramg.groups=int32(sum(A, 1));

wadj = 1/sqrt(d);

Xadj = X/wadj;

vec = paramg.groups;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 
k=length(ind)-1;     % number of groups
w =ones(k,1);

opts=[];

opts.init=2;        % starting from a zero point

opts.method = 1;
opts.maxIter=3e+5;


opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 

opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w; % set the weight for the group, a cloumn vector     %
opts.mFlag=1;       % treating it as compositive function 
opts.lFlag=1;   

lambda_seqglt = treelambdagen(Xadj,y,0,0,opts,nlambda,"slep","2");

%----------------------- group lasso in SPAM-----------------------

paramg = [];
paramg.num_threads=-1; % all cores (-1 by default)
paramg.verbose=false;   % verbosity, false by default
paramg.it0=50;      % frequency for duality gap computations
paramg.max_it=3e+5; % maximum number of iterations
paramg.L0=0.1;
paramg.tol=1e-5;
paramg.intercept=false;
paramg.pos=false;
W0=zeros(size(X,2),size(y,2));
paramg.loss='square';
paramg.regul='group-lasso-l2';
paramg.admm=true;

rownum = size(Gpar,1);
row_indices = 1:rownum ;
A = Gpar.* row_indices';
paramg.groups=int32(sum(A, 1));

wadj = 1/sqrt(d);

beta_realgl = beta_real * wadj;

Xadj = X/wadj;


tic;

[fvalgl,lambda_mint,lambda_min_patternt] = treecv(Xadj,y,W0,lasweight,paramg,lambda_seqglt,beta_realgl,"spam","2");


solve_time_spam_gl = toc;

paramg.lambda = lambda_mint  ;



[beta_spam_gl, fv1]=mexFistaFlat(y,Xadj,W0,paramg);

paramg.lambda = lambda_min_patternt;

[beta_spam_glp, fv2]=mexFistaFlat(y,Xadj,W0,paramg);




% disp("gl spam done")

fg1 = fv1(1);
fg2 = fv2(1);


beta_spam_gl = beta_spam_gl/wadj;

beta_realgl = beta_realgl/wadj;
RSSspamgl = sqrt(sum((beta_spam_gl - beta_reallas ).^2))./sqrt((sum(beta_realgl.^2)));

beta_estpattern = abs(beta_spam_glp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamspam_gl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

glspam_res = [solve_time_spam_gl;hamspam_gl;RSSspamgl];


%----------------------- group lasso in SLEP-----------------------
 
beta_realgl = beta_real * wadj;


tic;

[~,~,~] =  treecv(Xadj,y,W0,fvalgl,opts,lambda_seqglt,beta_realgl,"slep","2");



opts.tFlag = 2;
opts.tol = fg1;

z = lambda_mint;



[beta_slep_gl,~, ~]= glLeastR(Xadj, y, z, opts);

solve_time_slep_gl = toc;

zp = lambda_min_patternt;


opts.tFlag = 2;
opts.tol = fg2;

[beta_slep_glp, ~, ~]= glLeastR(Xadj, y, zp, opts);

beta_slep_gl = beta_slep_gl/wadj;

beta_realgl = beta_realgl/wadj;

RSSslepgl = sqrt(sum((beta_slep_gl - beta_reallas ).^2))./sqrt((sum(beta_realgl.^2)));

beta_estpattern = abs(beta_slep_glp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamslep_gl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);



gl_slep_res = [solve_time_slep_gl;hamslep_gl;RSSslepgl ];


%----------------------- lambda seq-----------------------
rownum = size(Gpar,1);
row_indices = 1:rownum ;
A = Gpar.* row_indices';
paramg.groups=int32(sum(A, 1));

wadj = 1;



Xadj = X/wadj;

vec = paramg.groups;
% Initialize variables
unique_numbers = unique(vec);
max_indices = zeros(size(unique_numbers));

% Find the maximum index for each unique number
for i = 1:length(unique_numbers)
    indices = find(vec == unique_numbers(i));
    max_indices(i) = indices(end);
end
ind = [0, max_indices]; 
k=length(ind)-1;     % number of groups
w =ones(k,1);

opts=[];

opts.init=2;        % starting from a zero point

opts.method = 1;
opts.maxIter=3e+5;


opts.rFlag=0;    
opts.nFlag=0;   
opts.rsL2=0;

% Termination 
opts.ind=ind.';       % set the group indices
opts.q=2;           % set the value for q % set the weight for positive and negative s
opts.gWeight=w; % set the weight for the group, a cloumn vector     %
opts.mFlag=1;       % treating it as compositive function 
opts.lFlag=1;       % Nemirovski's line search

lambda_seqglt = treelambdagen(Xadj,y,1,0,opts,nlambda,"slep","2");


%----------------------- group lasso in SPAM weight1-----------------------

paramg = [];
paramg.num_threads=-1; % all cores (-1 by default)
paramg.verbose=false;   % verbosity, false by default
paramg.it0=50;      % frequency for duality gap computations
paramg.max_it=3e+5; % maximum number of iterations
paramg.L0=0.1;
paramg.tol=1e-5;
paramg.intercept=false;
paramg.pos=false;
W0=zeros(size(X,2),size(y,2));
paramg.loss='square';
paramg.regul='group-lasso-l2';
paramg.admm=true;
% Create a row vector with values 1 to n
rownum = size(Gpar,1);
row_indices = 1:rownum ;
A = Gpar.* row_indices';
paramg.groups=int32(sum(A, 1));

wadj = 1;

beta_realgl = beta_real * wadj;

Xadj = X/wadj;




tic;

[fvalgl,lambda_mint,lambda_min_patternt] = treecv(Xadj,y,W0,lasweight,paramg,lambda_seqglt,beta_realgl,"spam","2");


solve_time_spam_gl = toc;


paramg.lambda = lambda_mint  ;

[beta_spam_gl, fv1]=mexFistaFlat(y,Xadj,W0,paramg);

paramg.lambda = lambda_min_patternt;

[beta_spam_glp, fv2]=mexFistaFlat(y,Xadj,W0,paramg);


fg1 = fv1(1);
fg2 = fv2(1);


beta_spam_gl = beta_spam_gl/wadj;

beta_realgl = beta_realgl/wadj;
RSSspamgl = sqrt(sum((beta_spam_gl - beta_reallas ).^2))./sqrt((sum(beta_realgl.^2)));

beta_estpattern = abs(beta_spam_glp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamspam_gl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);

glspam1_res = [solve_time_spam_gl;hamspam_gl;RSSspamgl];


%----------------------- group lasso in SLEP-----------------------


beta_realgl = beta_real;


tic;

[~,~,lambda_min_pattern] =  treecv(Xadj,y,W0,fvalgl,opts,lambda_seqglt,beta_realgl,"slep","2");




opts.tFlag = 2;
opts.tol = fg1;

z = lambda_mint;


[beta_slep_gl,~, ~]= glLeastR(Xadj, y, z, opts);

solve_time_slep_gl = toc;

opts.tFlag = 2;
opts.tol = fg2;


zp = lambda_min_patternt;

[beta_slep_glp, ~, ~]= glLeastR(Xadj, y, zp, opts);

beta_slep_gl = beta_slep_gl/wadj;

beta_realgl = beta_realgl/wadj;

RSSslepgl = sqrt(sum((beta_slep_gl - beta_reallas ).^2))./sqrt((sum(beta_realgl.^2)));

beta_estpattern = abs(beta_slep_glp);
beta_estpattern (abs(beta_estpattern)  < threshold) = 0;
beta_estpattern (abs(beta_estpattern)  >= threshold) = 1;
hamslep_gl = sum(beta_estpattern ~= beta_realpattern)/length(beta_realpattern);


gl1_slep_res = [solve_time_slep_gl;hamslep_gl;RSSslepgl ];




disp("all done")


%% Final
Res = [OGspam_res OGslep_res lasspam_res lasslep_res wglspam_res wgl_slep_res  glspam_res  gl_slep_res   glspam1_res  gl1_slep_res] ;
Res = full(Res);
end
  