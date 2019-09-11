function [Lambda, B, x, n] = KolmogorovEW(Dim, dT, x, f, df, h)
dX = x(2)-x(1);
n = size(x, 1);
p = @(x) -f(x);
q = @(x) -( sum(df(x), 2) + 0.5 * ( sum(h(x).^2, 2) ) );
e = ones(n,1);
d = 0.5/dX * spdiags([-e, e], [-1, 1], n, n);
switch Dim
    case 1
        P = spdiags(p(x), 0, n, n);
        Q = spdiags(q(x), 0, n, n);
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1-0.5*dT*Lambda;
        FW = P*d + Q;
        B = speye(n) + dT*FW;
        
    case 2
        [X,Y] = meshgrid(x);
        x = [X(:), Y(:)];
        n2 = n*n;
        I = speye(n);
        D1 = kron(I, d);
        D2 = kron(d, I);
        tmp = p(x);
        P1 = spdiags(tmp(:,1), 0, n2, n2);
        P2 = spdiags(tmp(:,2), 0, n2, n2);
        Q  = spdiags(q(x), 0, n2, n2);
        FW = P1*D1 + P2*D2 + Q;
        B  = speye(n2) + dT*FW;
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1-0.5*dT*Lambda;
        
    case 3
        [Y, Z, X] = meshgrid(x);
        x = [X(:), Y(:), Z(:)];
        n2 = n*n;
        n3 = n2*n;
        I  = speye(n);
        I2 = speye(n2);
        D1 = kron(I2, d);
        D2 = kron3(I, d, I);
        D3 = kron(d, I2);
        tmp = p(x);
        P1 = spdiags(tmp(:,1), 0, n3, n3);
        P2 = spdiags(tmp(:,2), 0, n3, n3);
        P3 = spdiags(tmp(:,3), 0, n3, n3);
        Q  = spdiags(q(x), 0, n3, n3);
        FW = P1*D1 + P2*D2 + P3*D3 + Q;
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1 - 0.5*dT*Lambda;
        B  = speye(n3) + dT*FW;
        
    case 4
        [X4, X3, X2, X1] = ndgrid(x);
        x = [X1(:), X2(:), X3(:), X4(:)];
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        I = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        D1 = kron(I3, d);
        D2 = kron3(I2, d, I);
        D3 = kron3(I, d, I2);
        D4 = kron(d, I3);
        tmp = p(x);
        P1 = spdiags(tmp(:,1), 0, n4, n4);
        P2 = spdiags(tmp(:,2), 0, n4, n4);
        P3 = spdiags(tmp(:,3), 0, n4, n4);
        P4 = spdiags(tmp(:,4), 0, n4, n4);
        Q  = spdiags(    q(x), 0, n4, n4);
        FW = P1*D1 + P2*D2 + P3*D3 + P4*D4 + Q;
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1-0.5*dT*Lambda;
        B  = speye(n4) + dT*FW;
        
    case 5
        [X5, X4, X3, X2, X1] = ndgrid(x);
        x = [X1(:), X2(:), X3(:), X4(:), X5(:)];
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        n5 = n4*n;
        I = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        I4 = speye(n4);
        D1 = kron(I4, d);
        D2 = kron3(I3, d, I);
        D3 = kron3(I2, d, I2);
        D4 = kron3(I, d, I3);
        D5 = kron(d, I4);
        tmp = p(x);
        P1 = spdiags(tmp(:,1), 0, n5, n5);
        P2 = spdiags(tmp(:,2), 0, n5, n5);
        P3 = spdiags(tmp(:,3), 0, n5, n5);
        P4 = spdiags(tmp(:,4), 0, n5, n5);
        P5 = spdiags(tmp(:,5), 0, n5, n5);
        Q  = spdiags(    q(x), 0, n5, n5);
        FW = P1*D1 + P2*D2 + P3*D3 + P4*D4 + P5*D5 + Q;
        B  = speye(n5) + dT*FW;
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1 - 0.5*dT*Lambda;
        
    case 6
        [X6, X5, X4, X3, X2, X1] = ndgrid(x);
        x = [X1(:), X2(:), X3(:), X4(:), X5(:), X6(:)];
        n2 = n*n;
        n3 = n2*n;
        n4 = n3*n;
        n5 = n4*n;
        n6 = n5*n;
        I = speye(n);
        I2 = speye(n2);
        I3 = speye(n3);
        I4 = speye(n4);
        I5 = speye(n5);
        D1 = kron(I5, d);
        D2 = kron3(I4, d, I);
        D3 = kron3(I3, d, I2);
        D4 = kron3(I2, d, I3);
        D5 = kron3(I, d, I4);
        D6 = kron(d, I5);
        tmp = p(x);
        P1 = spdiags(tmp(:,1), 0, n6, n6);
        P2 = spdiags(tmp(:,2), 0, n6, n6);
        P3 = spdiags(tmp(:,3), 0, n6, n6);
        P4 = spdiags(tmp(:,4), 0, n6, n6);
        P5 = spdiags(tmp(:,5), 0, n6, n6);
        P6 = spdiags(tmp(:,6), 0, n6, n6);
        Q  = spdiags(q(x), 0, n6, n6);
        FW = P1*D1 + P2*D2 + P3*D3 + P4*D4 + P5*D5 + P6*D6 + Q;
        B  = speye(n6) + dT*FW;
        Lambda = 1/(dX*dX) * LaplacianEW(Dim, n);
        Lambda = 1 - 0.5*dT*Lambda;
end


function M = kron3(M1, M2, M3)
M = kron( kron(M1, M2), M3);