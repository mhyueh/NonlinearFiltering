function [state, obser, s] = SimulateStateObser(T, Tinc, f, h, Dim)
nT     = T/Tinc + 1;
state  = zeros(nT, Dim);
obser  = zeros(nT, Dim);
sqrtdT = sqrt(Tinc);
s = rng('default');
% x0 ~ N(0,1), E[x0] = 0: the initial state
for t = 2:nT
    state(t,:) = state(t-1,:) + ( f( state(t-1,:) ) * Tinc ) + ( sqrtdT * randn(1,Dim) );
end
for t = 2:nT
    obser(t,:) = obser(t-1,:) + ( h( state(t-1,:) ) * Tinc ) + ( sqrtdT * randn(1,Dim) );
end
