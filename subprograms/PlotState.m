function PlotState(T, dT, state, Iu)
TotalT = 0:dT:T;
for ii = 1:size(state,2)
    figure(ii);
    plot(TotalT, state(:,ii), 'k-'); hold on
    if nargin == 4
        plot(TotalT, Iu(:,ii), 'b-');
    end
    xlabel('time','FontSize',16);
    ylabel('state','FontSize',16);
    if nargin == 4
        legend('States','Estimates');
    else
        legend('States');
    end
    hold off
end