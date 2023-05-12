function [ind] = out_skin(Xin,params)
%IN_SKIN finds if X is in the fish's skin + gap 
params.fish_length=params.fish_length+params.gap;
params.fish_width=params.fish_width+params.gap/2;
phi=atan2(Xin(:,2),Xin(:,1));

if(Xin(:,2)>0) %head
    ind=find(normrow(Xin)>=params.fish_width);
else       %body
    W=params.fish_width;
    L=params.fish_length;
    a=params.tail_angle;
    p=params.tail_p;
    xv=[-W;-p*W;L*(1-p)*sin(a);W*p;W];
    yv=L*[0;-p;-p-(1-p)*cos(a);-p;0];
    ind=find(~inpolygon(Xin(:,1),Xin(:,2),xv,yv));
end
% [Xs,Ns] = get_skin_polar(phi,params);
% ind=find(normrow(Xs)<(normrow(Xin)));
end

function N=normrow(y)
    N=sqrt(y(:,1).^2+y(:,2).^2);
end
