function d = LaplacianEW(Dim, n)
d = -4*sin((1:n).'*pi/(2*n+2)).^2;
switch Dim
    case 2
        D = spdiags(d, 0, n, n);
        I = speye(n);
        D = kron(I, D) + kron(D, I);
        d = diag(D);
    case 3
        D = spdiags(d, 0, n, n);
        n2 = n*n;
        I = speye(n);
        I2 = speye(n2);
        D = kron(I2, D) + kron3(I, D, I) + kron(D, I2);
        d = diag(D);
    case 4
        D = spdiags(d, 0, n, n);
        n2 = n*n;
        n3 = n2*n;
        I  = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        D = kron(I3,D) + kron3(I2,D,I) + kron3(I,D,I2) + kron(D,I3);
        d = diag(D);
    case 5
        D = spdiags(d, 0, n, n);
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        I  = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        I4 = speye(n4);
        D = kron(I4,D) + kron3(I3,D,I) + kron3(I2,D,I2) ...
            + kron3(I,D,I3) + kron(D,I4);
        d = diag(D);
    case 6
        D = spdiags(d, 0, n, n);
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        n5 = n4*n;
        I  = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        I4 = speye(n4);
        I5 = speye(n5);
        D = kron(I5, D) + kron3(I4, D, I) + kron3(I3, D, I2) ...
            + kron3(I2, D, I3) + kron3(I, D, I4) + kron(D, I5);
        d = diag(D);
    case 7
        D = spdiags(d, 0, n, n);
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        n5 = n4*n;
        n6 = n5*n;
        I  = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        I4 = speye(n4);
        I5 = speye(n5);
        I6 = speye(n6);
        D = kron(I6, D) + kron3(I5, D, I) + kron3(I4, D, I2) ...
            + kron3(I3, D, I3) + kron3(I2, D, I4) + kron3(I, D, I5) ...
            + kron(D, I6);
        d = diag(D);
end

function M = kron3(M1, M2, M3)
M = kron( kron(M1, M2), M3);
