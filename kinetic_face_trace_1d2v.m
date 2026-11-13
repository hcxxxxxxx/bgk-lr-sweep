function [trace,wall] = kinetic_face_trace_1d2v(H,v,wallVelocity)
%KINETIC_FACE_TRACE_1D2V Current upwind trace, including discrete diffuse walls.
% Recompute reflected density from CURRENT outgoing H, not the previous sweep.
N=size(H,1)/2;
assert(N==floor(N) && size(H,2)==numel(v.weight));
rhoLeft=H(1,v.vx<0)*(v.weight(v.vx<0).*abs(v.vx(v.vx<0)))/v.halfFlux(1);
rhoRight=H(end,v.vx>0)*(v.weight(v.vx>0).*abs(v.vx(v.vx>0)))/v.halfFlux(2);
wall=[rhoLeft;rhoRight]+2*wallVelocity(:)*v.vy';
trace=zeros(N+1,size(H,2));
trace(1,:)=H(1,:); trace(end,:)=H(end,:);
trace(1,v.inLeft)=wall(1,v.inLeft);
trace(end,v.inRight)=wall(2,v.inRight);
trace(2:N,v.inLeft)=H(2:2:end-2,v.inLeft);
trace(2:N,v.inRight)=H(3:2:end,v.inRight);
end
