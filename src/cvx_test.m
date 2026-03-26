n = 10;
Q = eye(5);
R = eye(5);

cvx_begin
    % variable x(n);
    % q = x(1:5);
    % u = x(6:10);
    variable q(5);
    variable u(5);
    minimize(quad_form(q, Q) + quad_form(u, R))
    subject to
        0.001 <= u <= 1;
        0 <= q(1) <= 180;
        0 <= q(2) <= 180;
cvx_end
x = [q; u]