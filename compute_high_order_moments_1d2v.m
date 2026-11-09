function hot = compute_high_order_moments_1d2v(H,velocity)
%COMPUTE_HIGH_ORDER_MOMENTS_1D2V Two-velocity BGK high-order moments.
% H: 2*N-by-Nv, or N-by-2-by-Nv as in the supplied specification.
% velocity.vx, .vy, .weight are vectors; weight ALREADY includes M0*dv.
% Rows may also be kinetic face traces, to form Rhat before differentiating.
vx=velocity.vx(:); vy=velocity.vy(:); w=velocity.weight(:);
assert(numel(vx)==numel(vy) && numel(vx)==numel(w) ...
    && all(isfinite([vx;vy;w])) && all(w>0),'Invalid velocity nodes/weights.');
if ndims(H)==3
    assert(size(H,2)==2,'3D H must be N-by-2-by-Nv.');
    H=reshape(permute(H,[2 1 3]),[],numel(w));
end
assert(isnumeric(H) && isreal(H) && size(H,2)==numel(w) && all(isfinite(H),'all'), ...
    'H must have one column per velocity node.');
v2=vx.^2+vy.^2;
basis=[vx.*(vx.^2-vy.^2-1),vy.*(2*vx.^2-1), ...
    vx.^2.*(v2-2)-(v2-1)];
R=H*(w.*basis);
hot=struct('Rxx',R(:,1),'Rxy',R(:,2),'RT',R(:,3),'R',R);
end
