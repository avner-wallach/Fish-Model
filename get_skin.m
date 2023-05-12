function [X] = get_skin(params)
%GET_SKIN returns position of fish skin 
fish_length=params.fish_length;
fish_width=params.fish_width;
tail_p=params.tail_p;
tail_angle=params.tail_angle;
%%
w=fish_width/2;
L=fish_length*tail_p;

%head
N=1e2;
phi=linspace(0,pi,N)';
X=w*[cos(phi) sin(phi)];

%flanks
t=linspace(0,1,N)';
% x0=w*L/fish_length;
x0=w*(1-tail_p);
% Xrf=[w*(L/fish_length+(1-L/fish_length)*t) L*(t-1)];
% Xlf=[-w*(L/fish_length+(1-L/fish_length)*t) L*(t-1)];
Xrf=[w*((1-tail_p)+tail_p*t) L*(t-1)];
Xlf=[-w*((1-tail_p)+tail_p*t) L*(t-1)];

%tail
t=linspace(0,1,N*(1-tail_p)/tail_p)';
if(tail_angle==0)
    xt=0;
    yt=-fish_length;
else
    yt=-(fish_length-L)/sqrt(1+tan(tail_angle)^2)-L;
    xt=(-yt-L)*tan(tail_angle);
end
Xrt=[x0+(xt-x0)*t -L+(yt+L)*t];
Xlt=[-x0+(xt+x0)*t -L+(yt+L)*t];

%combine
X=[X;flipud(Xlf);Xlt;flipud(Xrt);(Xrf)];

end

