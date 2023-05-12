function [ X_p,Q_p ] = get_fish_poles( params )
%  GET_FISH_POLES get location of fish charge poles 

%% parameters
fish_length=params.fish_length; %cm
tail_angle=params.tail_angle; %rad
p_density = params.p_density;     %poles/cm
m = params.m;              %number of negative poles
n=p_density*fish_length; %number of poles
tail_p=params.tail_p; %non-bending portion of body (8cm)
y_t=fish_length*tail_p;
%% generate pole positions head is at origin; fish main axis along y-axis
X_p=zeros(n,2);
y=-linspace(0,fish_length,n);
n_tail=find(abs(y)>y_t);
X_p(:,2)=y;
X_p(n_tail,1)=abs((y(n_tail)+y_t))*sin(tail_angle);
X_p(n_tail,2)=-y_t-abs((y(n_tail)+y_t))*cos(tail_angle);
Q_p=ones(n,1)/(n-m);
Q_p((end-m+1):end)=-1/m;
end

