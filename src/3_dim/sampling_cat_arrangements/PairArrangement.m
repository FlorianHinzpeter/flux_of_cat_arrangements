function X = PairArrangement(number,Npair,rE,fd,trials,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of Pair arrangements.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

X = zeros(3*trials,2*Npair);

for j=1:trials

x = zeros(3,2*Npair);

constr = 1;

while constr==1

x(:,1:Npair) = Randomarrangement_dynamic(1,rE+5*rE/100,Npair);

for i = 1:Npair
    
phi = 2*pi*rand;
theta = acos(2*rand-1);
    
x(1,Npair+i) = x(1,i)+fd*cos(phi)*sin(theta);
x(2,Npair+i) = x(2,i)+fd*sin(phi)*sin(theta);
x(3,Npair+i) = x(3,i)+fd*cos(theta);
    
end

constr = constraint(x(1,:),x(2,:),x(3,:),1,rE,2*Npair);


end

X(3*j-2:3*j,:) = x;

end

save(savepath,'X')

end