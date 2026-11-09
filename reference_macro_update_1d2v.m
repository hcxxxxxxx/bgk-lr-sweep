function [value,components,current] = reference_macro_update_1d2v(Q,edges,previous)
%REFERENCE_MACRO_UPDATE_1D2V NonOscillatory's cell macro stopping metric.
% Q: nodal [rho_pert ux uy T_pert], using mixed SOURCE moments in GSIS.
% Cell integrals are [int(1+rho), int(1+T), norm(int(velocity))].
% components: relative L2 updates of [density temperature speed magnitude].
% The reference initializes OLD fields to the stationary reference gas.
% At exactly zero denominator, unchanged zero has update 0; otherwise Inf.
edges=edges(:); h=diff(edges); N=numel(h);
assert(N>=1 && all(isfinite(edges)) && all(h>0),'Invalid cell edges.');
assert(isnumeric(Q) && isreal(Q) && isequal(size(Q),[2*N 4]) ...
    && all(isfinite(Q),'all'),'Q must be finite 2*N-by-4 nodal moments.');
integral=h.*(Q(1:2:end,:)+Q(2:2:end,:))/2;
current=[h+integral(:,1),h+integral(:,4),hypot(integral(:,2),integral(:,3))];
if nargin<3; previous=[h,h,zeros(N,1)]; end
assert(isequal(size(previous),[N 3]) && all(isfinite(previous),'all'), ...
    'Previous cell monitor must be finite N-by-3.');
numerator=sum((current-previous).^2,1); denominator=sum(current.^2,1);
components=zeros(1,3); active=denominator>0;
components(active)=sqrt(numerator(active)./denominator(active));
components(~active & numerator>0)=Inf;
value=max(components);
end
