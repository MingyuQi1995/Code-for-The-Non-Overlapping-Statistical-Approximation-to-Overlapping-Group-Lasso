%Function to perform gene analysis based on WIKI pathway
%Input: simudatawiki.mat
%Out: Solving time and predicted AUC


function[biores] = wiki()

%%  Oglasso

root = pwd;
root1 = strcat(root, "/SLEP-master");

addpath(genpath("./SLEP-master"));

%read data

load('realdata/simudatawiki.mat');
%load('realdata/realdata.mat');

X = data.X;
y = data.Y;
y = double(y);
y = 2*(y-0.5);

n = size(X,1);



%randomly assign 200 samples to the training data and 95 samples to the testing data.

samplesize = [1:n];
idxtrain = randsample(samplesize,200,false);
trainX = X(idxtrain,:);
trainy = y(idxtrain);
idxtest = setdiff(samplesize,idxtrain);
testX = X(idxtest,:);
testy = y(idxtest);
X = trainX;
y = trainy;


K = 5;


%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
%opts.tFlag=3; % the relative change is less than opts.tol
opts.tFlag=1;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-6;      % the tolerance parameter
%opts.rStartNum=100;

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


opts.G = double(data.opt_G');
opts.ind = data.opt_ind;

opts.rStartNum=100;

opts.maxIter2 = 3e+5;
opts.tol2= 1e-7;
opts.flag2=2;

fold1 = [1:40];
fold2 = [41:80];
fold3 = [81:120];
fold4 = [121:160];
fold5 = [161:200];
foldinf = [fold1;fold2;fold3;fold4;fold5];
%----------------------- tuning lambda -----------------------


nlambda =100;


%generate regularzation path

tic;

lambda_seq = lambda_gen_og_real(X,y,nlambda,opts);

%tuning parameter

[lambda_LL, ~, lambda_auc, ~,~,~] = cvreal(X, y, lambda_seq, opts, K,foldinf);

solve_time_og = toc;


z = [0, lambda_LL];
    
[beta_est_ogll, c, ~, ~]= overlapping_LogisticR(X, y, z, opts);
nvarll =  length(find(beta_est_ogll > 1e-10));
y_prob = 1./(1+exp(-testX*beta_est_ogll - c));
Ytest1 = 0.5*testy + 0.5;
y_LL = Ytest1.*log(y_prob) + (1- Ytest1).*log(1 - y_prob);

ogll = sum(y_LL);

sg_LL = 0;

t = size(opts.ind,2);
    for k = 1:t
        if any(beta_est_ogll(opts.G(opts.ind(1,k):opts.ind(2,k))) > 1e-10)
         sg_LL = sg_LL + 1;
        end
    end
    
    
z = [0, lambda_auc];

[beta_est_ogauc, c, ~, ~]= overlapping_LogisticR(X, y, z, opts);
nvarauc =  length(find(beta_est_ogauc > 1e-10 ));

y_prob = 1./(1+exp(-testX*beta_est_ogauc - c));
Ytest1 = 0.5*testy + 0.5;

[~,~,~,ogauc] = perfcurve(Ytest1, y_prob, 1);


OG_res = [solve_time_og;ogll;ogauc;nvarauc;nvarll];


%%LOG
%----------------------- Set optional items -----------------------

fold1 = [1:40];
fold2 = [41:80];
fold3 = [81:120];
fold4 = [121:160];
fold5 = [161:200];
foldinf = [fold1;fold2;fold3;fold4;fold5];

opts.G = double(data.opt_G');

inputVector = opts.G;
    % Initialize an empty vector to store the ordered values
    orderedVector = [];

    % Initialize a dictionary to keep track of the order of occurrence
    orderMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
    order = 1;

    % Loop through the input vector
    for i = 1:length(inputVector)
        value = inputVector(i);

        % Check if the value has occurred before
        if isKey(orderMap, value)
            % If yes, get its order and append to the orderedVector
            orderedVector = [orderedVector, orderMap(value)];
        else
            % If not, assign a new order and append to the orderedVector
            orderMap(value) = order;
            orderedVector = [orderedVector, order];
            order = order + 1;
        end
    end


 % Initialize the new matrix with the same number of rows as A
Xdup = zeros(size(X, 1), length(orderedVector));

% Loop through the values in B and duplicate columns from A
for i = 1:length(orderedVector)
    columnIndex = orderedVector(i);  % Get the column index from B
    Xdup(:, i) = X(:, columnIndex);  % Duplicate the column
end

 % Initialize the new matrix with the same number of rows as A
testXdup = zeros(size(testX, 1), length(orderedVector));

% Loop through the values in B and duplicate columns from A
for i = 1:length(orderedVector)
    columnIndex = orderedVector(i);  % Get the column index from B
    testXdup(:, i) = testX(:, columnIndex);  % Duplicate the column
end




% Starting point
opts=[];

opts.init=2;        % starting from a zero point

% Termination 
%opts.tFlag=3; % the relative change is less than opts.tol
opts.tFlag=1;
%opts.maxIter=3e+5;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-6; 
%opts.rStartNum=100;

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];

opts.gWeight=optind(3,:)';    % number of groups
opts.q=2;  

opts.mFlag=1;
% treating it as compositive function 
opts.lFlag=1; 
opts.rStartNum=100;

opts.maxIter2 = 3e+5;
%opts.maxIter2 = 2;
opts.tol2= 1e-7;
opts.flag2=2;

nlambda = 100;

tic;

lambda_seq = lambda_gen_log_real(Xdup,y,nlambda,opts,orderedVector);

[lambda_LL, ~, lambda_auc, ~,~,~] = cvwsgreal(Xdup, y, lambda_seq, opts, K,foldinf);
 

solve_time_log = toc;


z = lambda_LL;
    
[beta_est_logll, c, ~, ~]=  glLogisticR(Xdup, y, z, opts);
nvarll =  length(find(beta_est_logll > 1e-10));

y_prob = 1./(1+exp(-testXdup*beta_est_logll - c));
Ytest1 = 0.5*testy + 0.5;
y_LL = Ytest1.*log(y_prob) + (1- Ytest1).*log(1 - y_prob);

logll = sum(y_LL);

z = lambda_auc;
    
[beta_est_logauc, c, ~, ~]=  glLogisticR(Xdup, y, z, opts);
nvarauc =  length(find(beta_est_logauc > 1e-10));


y_prob = 1./(1+exp(-testXdup*beta_est_logauc - c));
Ytest1 = 0.5*testy + 0.5;

[~,~,~,logauc] = perfcurve(Ytest1, y_prob, 1);



log_res = [solve_time_log;logll; logauc;nvarauc;nvarll];

%%Yu
%----------------------- Set optional items -----------------------

fold1 = [1:40];
fold2 = [41:80];
fold3 = [81:120];
fold4 = [121:160];
fold5 = [161:200];
foldinf = [fold1;fold2;fold3;fold4;fold5];

opts.G = double(data.opt_G');

inputVector = opts.G;
    % Initialize an empty vector to store the ordered values
    orderedVector = [];

    % Initialize a dictionary to keep track of the order of occurrence
    orderMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
    order = 1;

    % Loop through the input vector
    for i = 1:length(inputVector)
        value = inputVector(i);

        % Check if the value has occurred before
        if isKey(orderMap, value)
            % If yes, get its order and append to the orderedVector
            orderedVector = [orderedVector, orderMap(value)];
        else
            % If not, assign a new order and append to the orderedVector
            orderMap(value) = order;
            orderedVector = [orderedVector, order];
            order = order + 1;
        end
    end


 % Initialize the new matrix with the same number of rows as A
Xdup = zeros(size(X, 1), length(orderedVector));

% Loop through the values in B and duplicate columns from A
for i = 1:length(orderedVector)
    columnIndex = orderedVector(i);  % Get the column index from B
    Xdup(:, i) = X(:, columnIndex);  % Duplicate the column
end

 % Initialize the new matrix with the same number of rows as A
testXdup = zeros(size(testX, 1), length(orderedVector));

% Loop through the values in B and duplicate columns from A
for i = 1:length(orderedVector)
    columnIndex = orderedVector(i);  % Get the column index from B
    testXdup(:, i) = testX(:, columnIndex);  % Duplicate the column
end




% Starting point
opts=[];

opts.init=2;        % starting from a zero point

% Termination 
%opts.tFlag=3; % the relative change is less than opts.tol
opts.tFlag=1;
%opts.maxIter=3e+5;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-5; 
%opts.rStartNum=100;

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];

opts.gWeight=optind(3,:)';    % number of groups
opts.q=1e7;  

opts.mFlag=0;
% treating it as compositive function 
opts.lFlag=0; 
opts.rStartNum=100;

opts.maxIter2 = 3e+5;
%opts.maxIter2 = 2;
opts.tol2= 1e-6;
opts.flag2=2;



nlambda = 100;

tic;

lambda_seq = lambda_gen_yu_real(Xdup,y,nlambda,opts,orderedVector);

[lambda_LL, ~, lambda_auc, ~,~,~] = cvwsgreal(Xdup, y, lambda_seq, opts, K,foldinf);
 

solve_time_yu = toc;


z = lambda_LL;
    

[beta_est_yull, c, ~, ~]=  glLogisticR(Xdup, y, z, opts);
nvarll =  length(find(beta_est_yull > 1e-10));

y_prob = 1./(1+exp(-testXdup*beta_est_yull - c));
Ytest1 = 0.5*testy + 0.5;
y_LL = Ytest1.*log(y_prob) + (1- Ytest1).*log(1 - y_prob);

yull = sum(y_LL);

z = lambda_auc;
    
[beta_est_yuauc, c, ~, ~]=  glLogisticR(Xdup, y, z, opts);
nvarauc =  length(find(beta_est_yuauc > 1e-10));


y_prob = 1./(1+exp(-testXdup*beta_est_yuauc - c));
Ytest1 = 0.5*testy + 0.5;

[~,~,~,yuauc] = perfcurve(Ytest1, y_prob, 1);



yu_res = [solve_time_yu;yull; yuauc;nvarauc;nvarll];

%csvwrite("./matlab_res.csv", [solve_timeog;run_Timeog;RSSog;beta_estog])
%% lasso
  
%----------------------- Set optional items -----------------------



load('realdata/simulassodatawiki.mat');



X = data.X;
y = data.Y;
y = double(y);
y = 2*(y-0.5);
trainX = X(idxtrain,:);
trainy = y(idxtrain);

testX = X(idxtest,:);
testy = y(idxtest);

X = trainX;
y = trainy;

K = 5;

n = size(X,1);

%----------------------- Set optional items -----------------------
% Starting point
opts=[];


opts.init=2;        % starting from a zero point

% Termination 
%opts.tFlag=3; % the relative change is less than opts.tol
opts.tFlag=1;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-6;      % the tolerance parameter

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

opts.mFlag=1;
% treating it as compositive function 
opts.lFlag=1; 

opts.rStartNum=100;
%----------------------- tuning lambda -----------------------

nlambda = 100;


tic;
lambda_seq = reallassolambda_gen(X,y,nlambda,opts);


[lambda_LL, ~, lambda_auc, ~,~,~] = cvlassoreal(X, y, lambda_seq, opts, K,foldinf);

solve_time_las = toc;



rho = lambda_LL;
    
[beta_est_lasll, c, ~, ~]= LogisticR(X, y, rho, opts);
nvarll =  length(find(beta_est_lasll > 1e-10 ));
   
y_prob = 1./(1+exp(-testX*beta_est_lasll - c));
Ytest1 = 0.5*testy + 0.5;
y_LL = Ytest1.*log(y_prob) + (1- Ytest1).*log(1 - y_prob);

lasll = sum(y_LL);



    
rho = lambda_auc;

[beta_est_lasauc, c, ~, ~]= LogisticR(X, y, rho, opts);
nvarauc =  length(find(beta_est_lasauc > 1e-10 ));
y_prob = 1./(1+exp(-testX*beta_est_lasauc - c));
Ytest1 = 0.5*testy + 0.5;
[~,~,~,lasauc] = perfcurve(Ytest1, y_prob, 1);


    
las_res = [solve_time_las;lasll; lasauc;nvarauc;nvarll];
    






%% wsglasso

load('realdata/simuwsgdatawiki.mat');

X = data.X;
y = data.Y;
y = double(y);
y = 2*(y-0.5);

trainX = X(idxtrain,:);
trainy = y(idxtrain);

testX = X(idxtest,:);
testy = y(idxtest);

X = trainX;
y = trainy;


K = 5;

n = size(X,1);


%----------------------- Set optional items -----------------------

% Starting point
opts=[];

opts.init=2;        % starting from a zero point

% Termination 
%opts.tFlag=3; % the relative change is less than opts.tol
opts.tFlag=1;
%opts.maxIter=3e+5;
opts.maxIter=3e+5;  % maximum number of iterations
opts.tol= 1e-6; 
%opts.rStartNum=100;

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


optind = data.opt_ind;
ind = optind(2,:);
opts.ind = [0,ind];

opts.gWeight=optind(3,:)';    % number of groups
opts.q=2;  

opts.mFlag=1;
% treating it as compositive function 
opts.lFlag=1; 
opts.rStartNum=100;

opts.maxIter2 = 3e+5;
%opts.maxIter2 = 2;
opts.tol2= 1e-7;
opts.flag2=2;

%----------------------- tuning lambda -----------------------

nlambda = 100;

tic;

lambda_seq = lambda_gen_wsg_real(X,y,nlambda,opts);


[lambda_LL, ~, lambda_auc, ~,~,~] = cvwsgreal(X, y, lambda_seq, opts, K,foldinf);
solve_time_wsg = toc;



z = lambda_LL;
    
[beta_est_wglll, c, ~, ~]=  glLogisticR(X, y, z, opts);
nvarll =  length(find(beta_est_wglll > 1e-10));

y_prob = 1./(1+exp(-testX*beta_est_wglll - c));
Ytest1 = 0.5*testy + 0.5;
y_LL = Ytest1.*log(y_prob) + (1- Ytest1).*log(1 - y_prob);

wglll = sum(y_LL);

z = lambda_auc;
    
[beta_est_wglauc, c, ~, ~]=  glLogisticR(X, y, z, opts);
nvarauc =  length(find(beta_est_wglauc > 1e-10));


y_prob = 1./(1+exp(-testX*beta_est_wglauc - c));
Ytest1 = 0.5*testy + 0.5;

[~,~,~,wglauc] = perfcurve(Ytest1, y_prob, 1);



    

wsg_res = [solve_time_wsg;wglll; wglauc;nvarauc;nvarll];






%% Final
biores = [OG_res las_res wsg_res,log_res,yu_res];


end