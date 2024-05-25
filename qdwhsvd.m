function [U,S,V] = qdwhsvd(A)
% fprintf("My custom qdwh_SVD implementation!\n")
[Up, H,it] = qdwh(A);
[V,S] = qdwheig(H);
U = Up*V;
end