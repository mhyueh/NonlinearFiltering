function U = DST_Solver(Dim, Lambda, B, U, n)
switch Dim
    case 1
        U = dstn(B*U);
        U = U./Lambda;
        U = dstn(U);
    otherwise
        SIZ = n*ones(1,Dim);
        U = dstn(reshape(B*U, SIZ));
        U = U(:)./Lambda;
        U = dstn(reshape(U, SIZ));
end
U = U(:);
U( U < 0 ) = 0;
U = U / sum(U);