function Flux = PostProcess(flux_path)

%Flux: alpha  J1 J2 std.dev.

Flux = zeros(7,4);

Alpha = [10^(-1) 10^(-0.5) 10^(0) 10^(0.5) 10^(1) 10^(1.5) 10^(2)];

for i=1:length(Alpha)

alpha = Alpha(i);    
    
load(flux_path)

Avg = sqrt(1/1500*sum((J(1:1500-1,2)-J(1501,2)).^2));

Flux(i,1) = alpha;
Flux(i,2:3) = J(end,:);
Flux(i,4) = Avg;

end