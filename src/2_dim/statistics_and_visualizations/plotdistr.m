function plotdistr(x1,x2,rE)
%---------------------------------------------------------------------------------------------
%this function visualizes the distribution of catalysts
%---------------------------------------------------------------------------------------------

% x1: coordinates of E1 (two rows)
% x2: coordinates of E2 (two rows)

hold off


circle([0,0],1,5000,'k-');

hold on

for i = 1:length(x1(1,:))

filledCircle(x1(:,i),rE,5000,'b');

end

for i = 1:length(x2(1,:))

filledCircle(x2(:,i),rE,5000,'r');

end

hold off