function [nn, L]  = nntrain_pairs_v5(nn, train_x, train_y, opts, index, distEuc, b, alpha, lambda)
%NNTRAIN trains a neural net
% [nn, L] = nnff(nn, x, y, opts) trains the neural network nn with input x and
% output y for opts.numepochs epochs, with minibatches of size
% opts.batchsize. Returns a neural network nn with updated activations,
% errors, weights and biases, (nn.a, nn.e, nn.W, nn.b) and L, the sum
% squared error for each training minibatch.

assert(isfloat(train_x), 'train_x must be a float');
assert(nargin == 7 || nargin == 9,'number of input arguments must be 6 or 8')

loss.train.e               = [];
loss.train.e_frac          = [];
loss.val.e                 = [];
loss.val.e_frac            = [];
opts.validation = 0;
if nargin == 9
    opts.validation = 1;
end

m = size(index, 1);
fhandle = [];
if isfield(opts,'plot') && opts.plot == 1
    fhandle = figure();
end

batchsize = opts.batchsize;
numepochs = opts.numepochs;
numbatches = int32(m / batchsize);

assert(rem(numbatches, 1) == 0, 'numbatches must be a integer');

L = zeros(numepochs*numbatches,1);
n = 1;
for i = 1 : numepochs
%     tic;

    for l = 1 : numbatches
        nid = ((l - 1) * batchsize + 1): (l * batchsize);
        iset = reshape(index(nid, :), batchsize * 2, 1);      
        batch_x = train_x(iset, :);
        
        %Add noise to input (for use in denoising autoencoder)
        if(nn.inputZeroMaskedFraction ~= 0)
            batch_x = batch_x.*(rand(size(batch_x))>nn.inputZeroMaskedFraction);
        end
        
        batch_y = train_y(iset, :);
        
        nn = nnff(nn, batch_x, batch_y);
        
        H = nn.a{end};
        distHam = nn.size(end) - sum(round(H(1:batchsize, :)) .* round(H((batchsize + 1):(batchsize * 2), :)), 2) - sum((1 - round(H(1:batchsize, :))) .* (1 - round(H((batchsize + 1):(2 * batchsize), :))), 2);
        distDiff = repmat((distHam - distEuc(nid) * b(2) - b(1)), 2, 1);
        L(n) = sum(sum(distDiff .^ 2)) / batchsize / 4 * alpha + nn.L; % * (1 - alpha)  
        H = round([H((batchsize + 1):(batchsize * 2), :); H(1:batchsize, :)]);
        distDiff = repmat(distDiff, 1, nn.size(end)) .* (1 - 2 * H);
        nn.e = nn.e * alpha - distDiff * lambda / nn.size(end);
        
        nn = nnbp(nn);
        nn = nnapplygrads_orghogonal(nn);
        
        n = n + 1;
    end
    
%     t = toc;
    t = 0;
    if ishandle(fhandle)
        nnupdatefigures(nn, fhandle, loss, opts, i);
    end
    
%     disp(['epoch ' num2str(i) '/' num2str(opts.numepochs) '. Took ' num2str(t) ' seconds' '. Mini-batch mean squared error on training set is ' num2str(mean(L((n-numbatches):(n-1))))]);
    nn.learningRate = nn.learningRate * nn.scaling_learningRate;
end
end