clear;clc;close all;
addpath('subprograms');

%% Settings
Dim    = 1;
T      = 20;
dT     = 0.001;
dTau   = 5*dT;
dX     = 0.5;

f  = @(x)   cos(x);
df = @(x) - sin(x);
h  = @(x)   x.^3;

%% Initialize
nTau   = T/dTau;
Tau    = 0:dTau:T;
nT     = dTau/dT;

fprintf('========================================\n');
fprintf([' ' num2str(Dim) '-D Yau-Yau Method using QIEM with DST\n']);
fprintf('========================================\n');

%% Generate States
fprintf('Generating States ...');
tic
[state, obser, s] = SimulateStateObser(T, dT, f, h, Dim);
toc
y = obser(1:nT:end,:);
rX = [min(min(state)), max(max(state))];
fprintf('Range(States) = [%f, %f].\n', rX(1), rX(2));

tstart = tic;
%% Construct Matrix
fprintf('Constructing Matrices ...\n');
x = (rX(1):dX:rX(2)+dX).';
[Lambda, B, x, n] = KolmogorovEW(Dim, dT, x, f, df, h);

%% Solve Kolmogorov Equations
sigma0 = exp( -10 * ( sum(x.^2, 2) ) );
Iu     = zeros((nTau-1)*nT+1, Dim);
Idx    = 1;
U      = sigma0;
U      = U / sum(U);
Iu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);

for jj = 1:nTau
    Idx = Idx + 1;
	fprintf('Computing step %d ... ', Idx);tic
    if jj == 1
        tmp = y(jj,:);
    else
        tmp = y(jj,:) - y(jj-1,:);
    end
    U = NormalizedExp( sum( h(x).*tmp(ones(size(x,1),1), :) , 2) ) .* U;
	U = DST_Solver(Dim, Lambda, B, U, n);
    U = U / sum(U);
	Iu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);
    toc
	for ii = 2:nT
        Idx = Idx + 1;
        fprintf('Computing step %d ... ', Idx);tic
        U = DST_Solver(Dim, Lambda, B, U, n);
        Iu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);
        toc
	end
end
telapsed = toc(tstart);

%% Plot the Result
PlotState(T, dT, state, Iu);
TotalT = 0:dT:T;
for ii = 1:size(state,2)
    figure(ii);
    plot(TotalT, state(:,ii), 'k-'); hold on
    plot(TotalT, Iu(:,ii), 'b-');
    xlabel('time','FontSize',16);
    ylabel('state','FontSize',16);
    legend('States','Estimates');
    hold off
end

Error_RMS = sqrt(mean((sum((state - Iu).^2, 2))/Dim, 1));
Error_M   = mean(sqrt(sum((state - Iu).^2, 2)/Dim), 1);

fprintf('=============================================================\n');
fprintf(['The computation of ' num2str(Dim) '-D Yau-Yau Method is finished.\n']);
fprintf('-------------------------------------------------------------\n');
fprintf(['Terminal Time         : ' num2str(T) ' \n']);
fprintf(['Space Range           : [' num2str(rX(1)) ', ' num2str(rX(2)) '] \n']);
fprintf(['Size of Time Steps    : ' num2str(dT) ' \n']);
fprintf(['Size of Space Steps   : ' num2str(dX) ' \n']);
fprintf(['Root-Mean-Square Error: ' num2str(Error_RMS) ' \n']);
fprintf(['Mean Error            : ' num2str(Error_M) ' \n']);
fprintf(['Time Costs            : ' num2str(telapsed) ' seconds. \n']);
fprintf('=============================================================\n');
