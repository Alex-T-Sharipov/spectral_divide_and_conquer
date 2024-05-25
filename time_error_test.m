function time_error_test(m, n)
    A = rand(m,n) + 1i * rand(m,n);
    tic;
    U_qdwh = qdwh(A);
    qdwh_time = toc;

    tic;
    [U_polar, ~] = my_poldec(A);
    polar_time = toc;

    error_qdwh = norm(U_qdwh' * U_qdwh - eye(size(A, 2)), 'fro');
    error_polar = norm(U_polar' * U_polar - eye(size(A, 2)), 'fro');

    fprintf('QDWH Error: %e, Time: %.2fs\n', error_qdwh, qdwh_time);
    fprintf('Polar Decomposition Error: %e, Time: %.2fs\n', error_polar, polar_time);
end

function [U, P] = my_poldec(A)
    [U, S, V] = svd(A, 'econ');  % Perform economic size SVD
    P = U * S * V';              % This still holds as P = A
    U = U * V';                  % Adjust U to be the unitary part
end