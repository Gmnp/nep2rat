function Z = halfdisksample(n,c,r)
% approximately n perturbed points on a disk 
% with center c and radius r

d0 = sqrt(pi/(2*n)); % est. distance between (unperturbed points)
fun = @(d) length(getpts(d))-n;
d = fzero(fun,d0);
Z = getpts(d);
n = length(Z);
rx = d*(rand(n,1)-0.5); % perturbation
ry = d*(rand(n,1)-0.5); 
Z = Z + rx + 1i*ry;
Z = r*Z + c; % shift and scale

end

function Z = getpts(d)
y = d/2:d:(1);
x = [ -1+d/2:d:1-d/2 ];
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y; Z = Z(:);
Z = Z(abs(Z) <= 1-sqrt(d^2/2));
end
