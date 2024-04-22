%Main function for simulation
%Input: gene structure 
%Output: Res (Regularization path computing time, estimation error, and support discrepancy of each method under the setting).


function[res] = realsimu(type)

nlambda = 100;
root = pwd;

root1 = strcat(root, "/SLEP-master");

addpath(genpath("./SLEP-master"));

if strcmp(type, '0')

    
load('realdata/simudatabio.mat');




X = data.X;

opt_G = data.opt_G;

opt_ind = data.opt_ind;

   p = size(X,2);

    beta  = normrnd(10,4,[1,p]);

    randomNumbers = randi([0, 1], [1, p]);

    ran = (randomNumbers - 0.5)*2;

    beta = beta .* ran;

     kks = 24;

   indtemp = find(opt_ind(3,:) < 30);

   Smpl = randsample(indtemp, kks);

   optzero = opt_ind;

   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
  beta_real = beta';

  Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;
  
  sgm = sqrt(noi);

  n = size(X,1);
  
  y  = X*beta_real +  normrnd(0,sgm,[1,n])';
  
  
opts=[];

opts.init=2;      

opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;     

opts.rFlag=0;    


opts.nFlag=0;  

opts.G = opt_G;

opts.ind = opt_ind;

opts.maxIter2 = 3e+5;

opts.rFlag=0;     

opts.rsL2=0;


tic;

lambda_seq = reallambdagen(X,y,nlambda,opts,'0');

[lambda_min,lambda_min_pattern] = realcv(X, y, lambda_seq, opts, beta_real,'0');

solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

threshold = 1e-20;

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
beta_estpatternog = abs(beta_est_ogp);

beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_real);

OG_res = [solve_time_og;ham_og;RSSog];
  
  
%% wlasso

load('realdata/simulassodatabio.mat');


opts=[];

opts.init=2;   

opts.tFlag=1;
opts.maxIter=3e+5;  
opts.tol= 1e-5;      
opts.rFlag=0;       
opts.nFlag=0;      
opts.mFlag=1;        
opts.lFlag=1;    

lasweight = data.W;

beta_realwlas = beta_real.* lasweight;

Xwlas = X ./lasweight';



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '1');


[lambda_min, lambda_min_pattern] =  realcv(Xwlas, y, lambda_seq, opts, beta_realwlas, '1');


solve_time_wlas = toc;


opts.init=2;  



rho = lambda_min;


[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);


beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_real);

wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];


%%  wsglasso

load('realdata/simuwsgdatabio.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=optind(3,:)'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg_res = [solve_time_wsg;ham_wsg;RSS_wsg];



%%  wsg1lasso

load('realdata/simuwsgdatabio.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=ones(length(optind(3,:)'),1); 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg1_res = [solve_time_wsg;ham_wsg;RSS_wsg];




%%  wsg1lasso

load('realdata/simuwsgdatabio.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight= sqrt(diff(opts.ind))';
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsgd_res = [solve_time_wsg;ham_wsg;RSS_wsg];

res = [OG_res wlas_res wsg_res wsg1_res wsgd_res];


elseif  strcmp(type, '1')

load('realdata/simudatapid.mat');

X = data.X;
opt_G = data.opt_G;
opt_ind = data.opt_ind;

   p = size(X,2);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
     kks = 13;
   indtemp = find(opt_ind(3,:) < 30);
   Smpl = randsample(indtemp, kks);
   optzero = opt_ind;
   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';

     Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;
  
  sgm = sqrt(noi);

  n = size(X,1);
  
  y  = X*beta_real +  normrnd(0,sgm,[1,n])';
  
  
opts=[];

opts.init=2;      

opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;     

opts.rFlag=0;    


opts.nFlag=0;  

opts.G = opt_G;

opts.ind = opt_ind;

opts.maxIter2 = 3e+5;

opts.rFlag=0;     

opts.rsL2=0;


tic;

lambda_seq = reallambdagen(X,y,nlambda,opts,'0');

[lambda_min,lambda_min_pattern] = realcv(X, y, lambda_seq, opts, beta_real,'0');

solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

threshold = 1e-20;

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
beta_estpatternog = abs(beta_est_ogp);

beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_real);

OG_res = [solve_time_og;ham_og;RSSog];
  
  
%% wlasso

load('realdata/simulassodatapid.mat');


opts=[];

opts.init=2;   

opts.tFlag=1;
opts.maxIter=3e+5;  
opts.tol= 1e-5;      
opts.rFlag=0;       
opts.nFlag=0;      
opts.mFlag=1;        
opts.lFlag=1;    

lasweight = data.W;

beta_realwlas = beta_real.* lasweight;

Xwlas = X ./lasweight';



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '1');


[lambda_min, lambda_min_pattern] =  realcv(Xwlas, y, lambda_seq, opts, beta_realwlas, '1');


solve_time_wlas = toc;


opts.init=2;  



rho = lambda_min;


[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);


beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_real);

wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];


%%  wsglasso

load('realdata/simuwsgdatapid.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=optind(3,:)'; 
opts.q=2;  





tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg_res = [solve_time_wsg;ham_wsg;RSS_wsg];




%%  wsg1lasso

load('realdata/simuwsgdatapid.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
pts.gWeight=ones(length(optind(3,:)'),1);
opts.q=2;  




tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg1_res = [solve_time_wsg;ham_wsg;RSS_wsg];

%%  wsgdlasso

load('realdata/simuwsgdatapid.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight= sqrt(diff(opts.ind))'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsgd_res = [solve_time_wsg;ham_wsg;RSS_wsg];

res = [OG_res wlas_res wsg_res wsg1_res wsgd_res];


 elseif strcmp(type, '2')

load('realdata/simudatakegg.mat');

X = data.X;
opt_G = data.opt_G;
opt_ind = data.opt_ind;

   p = size(X,2);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
     kks = 12;
   indtemp = find(opt_ind(3,:) < 10);
   Smpl = randsample(indtemp, kks);
   optzero = opt_ind;
   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';

     Expect = X*beta_real;
  
  sig = mean(abs(Expect));

    noi =   sig/3;
  
  sgm = sqrt(noi);

  n = size(X,1);
  
  y  = X*beta_real +  normrnd(0,sgm,[1,n])';
  
  
opts=[];

opts.init=2;      

opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;     

opts.rFlag=0;    


opts.nFlag=0;  

opts.G = opt_G;

opts.ind = opt_ind;

opts.maxIter2 = 3e+5;

opts.rFlag=0;     

opts.rsL2=0;


tic;

lambda_seq = reallambdagen(X,y,nlambda,opts,'0');

[lambda_min,lambda_min_pattern] = realcv(X, y, lambda_seq, opts, beta_real,'0');

solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

threshold = 1e-20;

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
beta_estpatternog = abs(beta_est_ogp);

beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_real);

OG_res = [solve_time_og;ham_og;RSSog];

  
%% wlasso

load('realdata/simulassodatakegg.mat');

opts=[];

opts.init=2;   

opts.tFlag=1;
opts.maxIter=3e+5;  
opts.tol= 1e-5;      
opts.rFlag=0;       
opts.nFlag=0;      
opts.mFlag=1;        
opts.lFlag=1;    

lasweight = data.W;

beta_realwlas = beta_real.* lasweight;

Xwlas = X ./lasweight';



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '1');


[lambda_min, lambda_min_pattern] =  realcv(Xwlas, y, lambda_seq, opts, beta_realwlas, '1');


solve_time_wlas = toc;


opts.init=2;  



rho = lambda_min;


[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);


beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_real);

wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];

%%  wsglasso

load('realdata/simuwsgdatakegg.mat');


X = data.X;

X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=optind(3,:)'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg_res = [solve_time_wsg;ham_wsg;RSS_wsg];

%%  wsg1lasso

load('realdata/simuwsgdatakegg.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=ones(length(optind(3,:)'),1); 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg1_res = [solve_time_wsg;ham_wsg;RSS_wsg];


%%  wsg1lasso

load('realdata/simuwsgdatakegg.mat');


X = data.X;

opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight= sqrt(diff(opts.ind))'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsgd_res = [solve_time_wsg;ham_wsg;RSS_wsg];


res = [OG_res wlas_res wsg_res wsg1_res wsgd_res];





 elseif strcmp(type, '3')

load('realdata/simudatawiki.mat');
X = data.X;
opt_G = data.opt_G;
opt_ind = data.opt_ind;

   p = size(X,2);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
     kks = 25;
   indtemp = find(opt_ind(3,:) < 15);
   Smpl = randsample(indtemp, kks);
   optzero = opt_ind;
   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';

     Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;
  
  sgm = sqrt(noi);

  n = size(X,1);
  
  y  = X*beta_real +  normrnd(0,sgm,[1,n])';
  
  
opts=[];

opts.init=2;      

opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;     

opts.rFlag=0;    


opts.nFlag=0;  

opts.G = opt_G;

opts.ind = opt_ind;

opts.maxIter2 = 3e+5;

opts.rFlag=0;     

opts.rsL2=0;


tic;

lambda_seq = reallambdagen(X,y,nlambda,opts,'0');

[lambda_min,lambda_min_pattern] = realcv(X, y, lambda_seq, opts, beta_real,'0');

solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

threshold = 1e-20;

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
beta_estpatternog = abs(beta_est_ogp);

beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_real);

OG_res = [solve_time_og;ham_og;RSSog];
  
%% wlasso

load('realdata/simulassodatawiki.mat');

opts=[];

opts.init=2;   

opts.tFlag=1;
opts.maxIter=3e+5;  
opts.tol= 1e-5;      
opts.rFlag=0;       
opts.nFlag=0;      
opts.mFlag=1;        
opts.lFlag=1;    

lasweight = data.W;

beta_realwlas = beta_real.* lasweight;

Xwlas = X ./lasweight';



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '1');


[lambda_min, lambda_min_pattern] =  realcv(Xwlas, y, lambda_seq, opts, beta_realwlas, '1');


solve_time_wlas = toc;


opts.init=2;  



rho = lambda_min;


[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);


beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_real);

wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];


%%  wsglasso

load('realdata/simuwsgdatawiki.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=optind(3,:)'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg_res = [solve_time_wsg;ham_wsg;RSS_wsg];


%%  wsglasso

load('realdata/simuwsgdatawiki.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=ones(length(optind(3,:)'),1); 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg1_res = [solve_time_wsg;ham_wsg;RSS_wsg];


%%  wsgdlasso

load('realdata/simuwsgdatawiki.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight= sqrt(diff(opts.ind))'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsgd_res = [solve_time_wsg;ham_wsg;RSS_wsg];


res = [OG_res wlas_res wsg_res wsg1_res wsgd_res];



 elseif strcmp(type, '4')
 
load('realdata/simudatareac.mat');
X = data.X;
opt_G = data.opt_G;
opt_ind = data.opt_ind;


p = size(X,2);
beta  = normrnd(10,4,[1,p]);
randomNumbers = randi([0, 1], [1, p]);
ran = (randomNumbers - 0.5)*2;
beta = beta .* ran;
kks = 80;
indtemp = find(opt_ind(3,:) > 8);
Smpl = randsample(indtemp, kks);
optzero = opt_ind;
optzero(:,Smpl) = [];
for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';
   con = 1;

while(con)
   p = size(X,2);
    beta  = normrnd(10,4,[1,p]);
    randomNumbers = randi([0, 1], [1, p]);
    ran = (randomNumbers - 0.5)*2;
    beta = beta .* ran;
     kks = 80;
   indtemp = find(opt_ind(3,:) > 7);
   Smpl = randsample(indtemp, kks);
   optzero = opt_ind;
   optzero(:,Smpl) = [];
   for i = 1:length(optzero)
    idx_tmpa = optzero(1,i);
    idx_tmpb = optzero(2,i);
    idx_tmp = opt_G(idx_tmpa:idx_tmpb);
    beta(idx_tmp) = 0;
   end
   
   beta_real = beta';
   if(nnz(beta_real)< 100)
       con = 1;
   else
      con = 0;
   end
    
end
    
Expect = X*beta_real;
  
  sig = mean(abs(Expect));
  
  noi =   sig/3;

  sgm = sqrt(noi);

  n = size(X,1);
  
  y  = X*beta_real +  normrnd(0,sgm,[1,n])';
  
  
 opts=[];

opts.init=2;      

opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;     

opts.rFlag=0;    


opts.nFlag=0;  

opts.G = opt_G;

opts.ind = opt_ind;

opts.maxIter2 = 3e+5;

opts.rFlag=0;     

opts.rsL2=0;


tic;

lambda_seq = reallambdagen(X,y,nlambda,opts,'0');

[lambda_min,lambda_min_pattern] = realcv(X, y, lambda_seq, opts, beta_real,'0');

solve_time_og = toc;


opts.init=2;  

z = [0, lambda_min];

[beta_est_og, ~, ~]=  overlapping_LeastR(X, y, z, opts);

zp = [0, lambda_min_pattern];

[beta_est_ogp, ~, ~]=  overlapping_LeastR(X, y, zp, opts);

RSSog = sqrt(sum((beta_est_og - beta_real ).^2))./sqrt((sum(beta_real.^2)));

threshold = 1e-20;

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);
beta_estpatternog = abs(beta_est_ogp);

beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_og = hamm/length(beta_real);

OG_res = [solve_time_og;ham_og;RSSog];

  
  
%% wlasso

load('realdata/simulassodatareac.mat');
opts=[];

opts.init=2;   

opts.tFlag=1;
opts.maxIter=3e+5;  
opts.tol= 1e-5;      
opts.rFlag=0;       
opts.nFlag=0;      
opts.mFlag=1;        
opts.lFlag=1;    

lasweight = data.W;

beta_realwlas = beta_real.* lasweight;

Xwlas = X ./lasweight';



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '1');


[lambda_min, lambda_min_pattern] =  realcv(Xwlas, y, lambda_seq, opts, beta_realwlas, '1');


solve_time_wlas = toc;


opts.init=2;  



rho = lambda_min;


[beta_est_wlas, ~, ~]=  LeastR(Xwlas, y, rho, opts);


beta_est_wlas = beta_est_wlas./lasweight;

RSS_wlas = sqrt(sum((beta_est_wlas - beta_realwlas).^2))./sqrt((sum(beta_realwlas.^2)));

rhop = lambda_min_pattern;

[beta_est_wlasp, ~, ~]=  LeastR(Xwlas, y, rhop, opts);


beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wlasp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wlas = hamm/length(beta_real);

wlas_res = [solve_time_wlas;ham_wlas;RSS_wlas];

%%  wsglasso


load('realdata/simuwsgdatareac.mat');


X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=optind(3,:)'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg_res = [solve_time_wsg;ham_wsg;RSS_wsg];



%%  wsg1lasso


load('realdata/simuwsgdatareac.mat');


X = data.X;


opts=[];

opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight=ones(length(optind(3,:)'),1);
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init=2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsg1_res = [solve_time_wsg;ham_wsg;RSS_wsg];

%%  wsgdlasso

load('realdata/simuwsgdatareac.mat');



X = data.X;


opts=[];


  
opts.init=2;        
opts.method = 1;
opts.tFlag=1;
opts.maxIter=3e+5; 
opts.tol= 1e-5;   
opts.rFlag=0;      
opts.nFlag=0;  
optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];
opts.gWeight= sqrt(diff(opts.ind))'; 
opts.q=2;  



tic;

lambda_seq = reallambdagen(X,y,nlambda,opts, '2');

[lambda_min, lambda_min_pattern]  = realcv(X, y, lambda_seq, opts, beta_real, '2');

solve_time_wsg = toc;

opts.init =2;  

z = lambda_min;


[beta_est_wsg, ~, ~]= glLeastR(X, y, z, opts);


RSS_wsg = sqrt(sum((beta_est_wsg -  beta_real).^2))./sqrt((sum( beta_real.^2)));

zp = lambda_min_pattern;

[beta_est_wglp, ~, ~]= glLeastR(X, y, zp, opts);

beta_realpattern = abs(beta_real);
beta_realpattern(abs(beta_realpattern) < threshold) = 0;
beta_realpattern(abs(beta_realpattern) >= threshold) = 1;
positions = find(abs(beta_realpattern) > threshold);

beta_estpatternog = abs(beta_est_wglp);
beta_estpatternog (abs(beta_estpatternog)  < threshold) = 0;
beta_estpatternog (abs(beta_estpatternog) >= threshold) = 1;

hamm = sum(beta_estpatternog ~= beta_realpattern);


ham_wsg = hamm/length(beta_real);
wsgd_res = [solve_time_wsg;ham_wsg;RSS_wsg];

res = [OG_res wlas_res wsg_res wsg1_res wsgd_res];


 elseif strcmp(type, '5')
     disp('no data');
 elseif strcmp(type, '6')
      disp('no data');
 else
     disp('Type is not in the range 0-6');
 end


end