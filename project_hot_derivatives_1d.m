function derivative = project_hot_derivatives_1d(mesh,R,Rhat,method)
%PROJECT_HOT_DERIVATIVES_1D DG weak derivative including numerical face traces.
% R: 2*N-by-k nodal moments. Rhat: (N+1)-by-k moments of the actual
% kinetic upwind trace, NOT an average of left/right moment values.
% method='weak' (default): M*Rx=-G*R+bR*RhatR-bL*RhatL.
% method='strong': M*Rx=G'*R, the element-interior derivative used by the
% NonOscillatory reference coupling; it does not add numerical face terms.
if nargin<4; method='weak'; end
method=validatestring(method,{'weak','strong'});
edges=mesh.edges(:); N=numel(edges)-1;
assert(all(isfinite(edges)) && all(diff(edges)>0),'Invalid mesh edges.');
assert(isnumeric(R) && isreal(R) && size(R,1)==2*N && all(isfinite(R),'all'), ...
    'R must be a 2*N-by-k nodal array.');
if strcmp(method,'weak')
    assert(isnumeric(Rhat) && isreal(Rhat) ...
        && isequal(size(Rhat),[N+1 size(R,2)]) && all(isfinite(Rhat),'all'), ...
        'Weak derivatives require (N+1)-by-k kinetic face moments.');
end
G=[-1 -1;1 1]/2; derivative=zeros(size(R));
for i=1:N
    nodes=2*i-1:2*i; M=(edges(i+1)-edges(i))*[2 1;1 2]/6;
    if strcmp(method,'weak')
        load=-G*R(nodes,:)+[-Rhat(i,:);Rhat(i+1,:)];
    else
        load=G'*R(nodes,:);
    end
    derivative(nodes,:)=M\load;
end
end
