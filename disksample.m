function Z = disksample(n,c,r)
% approximately n perturbed points on a disk 
% with center c and radius r

d0 = sqrt(pi/n); % est. distance between (unperturbed points)
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
x = d/2:d:(1-d/2);
x = [ -x(end:-1:1) , x ];
[X,Y] = meshgrid(x,x);
Z = X + 1i*Y; Z = Z(:);
Z = Z(abs(Z) <= 1-sqrt(d^2/2));
end
