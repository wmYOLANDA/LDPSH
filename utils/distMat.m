%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute Euclidean distance between vectors
% each vector is one row
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distMat(P1, P2)
if nargin == 2
    P1 = double(P1);
    P2 = double(P2);
    X1 = repmat(sum(P1 .^ 2, 2), [1, size(P2, 1)]);
    X2 = repmat(sum(P2 .^ 2, 2), [1, size(P1, 1)]);
    R = P1 * P2';
    D = real(sqrt(X1 + X2' - 2 * R));
else
    P1 = double(P1);
    X1 = repmat(sum(P1 .^ 2, 2), [1, size(P1, 1)]);
    R = P1 * P1';
    D = real(sqrt(X1 + X1' - 2 * R));
end
end