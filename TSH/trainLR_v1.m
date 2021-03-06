%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The matlab code for ��A General Framework for Linear Distance Preserving
% Hashing��, in TIP 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Btrain, Btest] = trainLR_v1(Xtrain, Btrain, Utrain, Xdata, Xtest, iterNum, nNeighborNum, numepoch1, numepoch2, nbits)
% Input:
%          Xtrain: n*d, n is the number of training samples
%          Btrain: Pseudo compact codes of training samples, generated by
%          one existing hashing method
%          Utrain: Pseudo binary codes of training samples, generated by
%          one existing hashing method
%          Xdata: database samples
%          Xtest: query samples
%          iterNum: the number of alternative training
%          nNeighborNum: the number of pairs for each training sample
%          numepoch1: iteration number for pretraining networks
%          numepoch2: iteration number for training step
%          nbits---encoding length
% Output:
%          Btrain: the compact codes of database samples
%          Btest: the compact codes of query samples

addpath('./NN/');
[ntrain, dim] = size(Xtrain);

% prepare data
[train_x, mu, sigma] = zscore(Xtrain);
test_xd = normalize(Xdata, mu, sigma);
test_x = normalize(Xtest, mu, sigma);
train_y = Utrain;

% compute the distances between the randomly selected pairs of data points
index = randi(ntrain, ntrain  * nNeighborNum, 2);
index = int32(index - 1);
[distHam, distEuc] = mexHamDisSMP(train_x', Btrain', ntrain * nNeighborNum, index', 1);

% parameters of our LDPSH methods
alpha = 0.1;
lambda = 0.5;
beta = 1e-4;

% initialize the neural network
rand('state',0)
nn = nnsetup([size(Xtrain, 2) size(Xtrain, 2) nbits]);
nn.weightPenaltyL2 = beta;  %  L2 weight decay
nn.dropoutFraction = 0.05;
opts.numepochs =  numepoch1;   %  Number of full sweeps through data
opts.batchsize = 500;  %  Take a mean gradient step over this many samples
[nn, L] = nntrain(nn, train_x, train_y, opts);
opts.numepochs =  numepoch2; 

x = ones(size(distEuc, 1), 2);
x(:, 2) = distEuc;
for i = 1:iterNum
    % linear regression
    b = regress(distHam, x);
    
    % train the neural network
    index = index + 1;
    [nn, L] = nntrain_pairs_v5(nn, train_x, Utrain, opts, index, distEuc, b, alpha, lambda);
    index = int32(index - 1);
    nn.testing = 1;
    nn = nnff(nn, train_x, zeros(size(train_x,1), nn.size(end)));
    nn.testing = 0;
    Utrain = round(nn.a{end});
    Btrain = compactbit(Utrain);
    [distHam, ~] = mexHamDisSMP(train_x', Btrain', ntrain * nNeighborNum, index', 0);
end
clear Xtrain Btrain Utrain Xdata Xtest train_x train_y distEuc distHam index;

nn.testing = 1;
nn = nnff(nn, test_xd, zeros(size(test_xd,1), nn.size(end)));
Utrain = round(nn.a{end});
Btrain = compactbit(Utrain);
clear Utrain;

nn.testing = 1;
nn = nnff(nn, test_x, zeros(size(test_x,1), nn.size(end)));
Utest = round(nn.a{end});
Btest = compactbit(Utest);
clear Utest;
end