function ref = channel1D2VReference(cfg,cells,order)
%CHANNEL1D2VREFERENCE Independent characteristic-integral BGK reference.
% For isothermal shear/force flow H=2*vy*phi(x,vx), so
% vx*phi_x + delta*phi = delta*u + a, u=int phi*exp(-vx^2)/sqrt(pi) dvx.
% Integrate transport exactly in each cell with piecewise-constant u,
% then solve the dense moment integral equation (no DG/source iteration).
% Positive speeds use Gauss-Legendre integration on [0,8]; negative speeds
% follow by symmetry. Values returned are exact spatial cell averages for
% this integral discretization, not samples at cell centers.
if nargin<2; cells = 512; end
if nargin<3; order = 160; end
validateattributes(cells,{'numeric'},{'scalar','integer','positive'});
validateattributes(order,{'numeric'},{'scalar','integer','>=',4});
k = (1:order-1)'; b = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(b,1)+diag(b,-1));
[z,idx] = sort(diag(D));
s = 4*(z+1); w = 8*V(1,idx)'.^2.*exp(-s.^2)/sqrt(pi);
dx = cfg.gap/cells; t = cfg.delta*dx./s;
loss = -expm1(-t); attenuation = exp(-t); average = loss./t;
self = 1-average;
small = t<1e-3; ts = t(small);
self(small) = ts/2-ts.^2/6+ts.^3/24-ts.^4/120+ts.^5/720;
kernel = zeros(cells,1); kernel(1) = 2*sum(w.*self);
wall = zeros(cells,1);
for i = 1:cells
    wall(i) = sum(w.*average.*attenuation.^(i-1));
    if i>1
        kernel(i) = sum(w.*average.*loss.*attenuation.^(i-2));
    end
end
K = toeplitz(kernel);
boundary = cfg.wallVelocity(1)*wall+cfg.wallVelocity(2)*flipud(wall);
rhs = boundary+(cfg.acceleration/cfg.delta)*sum(K,2);
u = (eye(cells)-K)\rhs;
ref = struct('x',-cfg.gap/2+((1:cells)'-0.5)*dx, ...
    'velocity',u,'flowRate',dx*sum(u),'cells',cells,'velocityOrder',order, ...
    'method','Characteristic integral equation; piecewise-constant spatial moments', ...
    'positiveSpeeds',s,'positiveWeights',w, ...
    'linearResidual',max(abs(u-K*u-rhs)));
end
