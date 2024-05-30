function [V,D] = qdwheig(A)
% fprintf("My custom qdwh_EIG implementation!\n")
if max(size(A)) == 1
    V = [1];
    D = A;
    return
end
s = medianestimate(A);
dim = max(size(A));
[U,H,it] = qdwh(A-s*eye(dim));
C = (eye(dim)+U)/2;
[V1,V2,W] = subspaceiteration(A,C,10^(-8));
A1 = V1'*A*V1;
A2 = V2'*A*V2;
[W1,D1] = qdwheig(A1);
[W2,D2] = qdwheig(A2);
D = blkdiag(D1, D2);
V = W*blkdiag(W1, W2);
end

function [s] = medianestimate(A)
s = median(diag(A));
end

function [U1,U2,U] = subspaceiteration(A,C,eps)
dim = max(size(C));
r = round(norm(C,"fro")^2);
X = randn(dim,r);
[Q,R] = qr(X);
V1 = Q(1:dim,1:r);
V2 = Q(1:dim,r+1:dim);
E = V2'*A*V1;
while norm(E,"fro")/norm(A,"fro")>eps
    X = C*V1;
    [Q,R] = qr(X);
    V1 = Q(1:dim,1:r);
    V2 = Q(1:dim,r+1:dim);
    E = V2'*A*V1;
end
U1 = V1;
U2 = V2;
U = Q;
end

