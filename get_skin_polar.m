function [ X,N] = get_skin_polar(phi,params)
%GET_SKIN returns position and norm of fish skin at 'phase' phi
%   phi- phase on skin. 0-tip of head, 0->pi positive side, pi-tip of tail
fish_length=params.fish_length;
fish_width=params.fish_width;
% tail_p=params.tail_p;
% tail_angle=params.tail_angle;

w=fish_width/2;
phi=phi(:);
phi(phi<0)=phi(phi<0)+2*pi;
y_h=w*sin(phi);
x_h=w*cos(phi);
N_h=[cos(phi) sin(phi)];

% alpha=-phi+pi/2;
alpha=phi;
r=fish_length./(-fish_length/w.*cos(alpha)-sin(alpha));
y_p=r.*sin(alpha);
x_p=r.*cos(alpha);
N_p=[-fish_length -w];
N_p=N_p/norm(N_p); %normalize

% alpha=phi-pi/2;
r=fish_length./(fish_length/w.*cos(alpha)-sin(alpha));
y_n=r.*sin(alpha);
x_n=r.*cos(alpha);        
N_n=[fish_length -w];
N_n=N_n/norm(N_n); %normalize        

X=zeros(numel(phi),2);
N=X;
ind=find(phi<=pi);
X(ind,:)=[x_h(ind) y_h(ind)];
N(ind,:)=N_h(ind,:);
ind=find(phi>pi & phi<=3*pi/2);
X(ind,:)=[x_p(ind) y_p(ind)];
N(ind,1)=N_p(1);
N(ind,2)=N_p(2);
ind=find(phi>3*pi/2);
X(ind,:)=[x_n(ind) y_n(ind)];
N(ind,1)=N_n(1);
N(ind,2)=N_n(2);
end

