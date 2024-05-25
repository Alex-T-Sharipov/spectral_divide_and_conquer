function [X_k,H,k] = qdwh(A)
    max_iter = 20;
    eps = 2.22e-16;
    tol = (5*eps)^(1/3);

    [m, n] = size(A);

    a_hat = normest(A,0.1);
    X_k = A / a_hat;

    I = eye(n);

    Y = X_k; if m > n, [Q,Y] = qr(X_k,0); end
    smin_est =  norm(Y,1)/condest(Y);  % Actually an upper bound for smin.
    l_k = smin_est/sqrt(n); 
    k = 0;
    while k < max_iter
        a_k = h(l_k);
        b_k = 0.25 * (a_k - 1)^2;
        c_k = a_k + b_k - 1;

        augmented_matrix = [sqrt(c_k) * X_k; I];
        [Q, ~] = qr(augmented_matrix, 0);
        Q1 = Q(1:m, :);
        Q2 = Q(m+1:end, :);

        X_k1 = (b_k / c_k) * X_k + (1 / sqrt(c_k)) * (a_k - b_k / c_k) * Q1 * Q2';

        if norm(X_k1 - X_k, 'fro') < tol
            break;
        end

        l_k1 = l_k * (a_k + b_k * l_k^2) / (1 + c_k * l_k^2);

        X_k = X_k1;
        l_k = l_k1;
        k = k + 1;
    end
    H = X_k'*A; H = (H'+H)/2;
end

function d_val = d(l)
    d_val = nthroot(4 * (1 - l^2) / l^4, 3);
end

function h_val = h(l)
    d_val = d(l);
    inner_sqrt = sqrt(8 - 4 * d_val + (8 * (2 - l^2)) / ((l^2) * sqrt(1 + d_val)));
    h_val = sqrt(1 + d_val) + 0.5 * inner_sqrt;
end