function [ Z,X,Y,F] = get_object_effect(params)
%GET_WALL_EFFECT get the effect of wall in different positions
%% parameters
r_max=params.r_max;
M=params.grid_M;
x_lim=[-r_max r_max];
y_lim=[-r_max r_max];

bgcol=params.bgcol;
circular=params.bgcol;
phi0=params.phi0;
%% generate poles
[ X_p,Q_p ] = get_fish_poles( params );

%% baseline
[X0,N0]=get_skin_polar(phi0,params);
[ v0,e0 ] = get_potential_field(X0,X_p,Q_p);
z0=e0*N0';
% generate skin
[X_s] = get_skin(params);

%% grid
x_v=linspace(x_lim(1),x_lim(2),M);
y_v=linspace(y_lim(1),y_lim(2),M+1);
Z=nan(M,M+1);
X=Z;
Y=Z;
for i=1:numel(x_v)
    x=x_v(i);
    for j=1:numel(y_v)
        y=y_v(j);
        X(i,j)=x;
        Y(i,j)=y;
        params.object_x=[x y];
        params.wall_angle=-atan2(y,x);
        params.wall_dist=params.tank_radius-norm([x y]);
        
        flag=numel(out_skin([x y],params));
%         [ X_o,Q_o] = mirror_object(X_p,Q_p,params);
%         X_all=[X_p;X_o];
%         Q_all=[Q_p;Q_o];
        if(~flag)
            continue;
        end
        
        [ X_w,Q_w,flag] = mirror_wall(X_p,Q_p,params); %wall pole image
        X_all=[X_p;X_w];
        Q_all=[Q_p;Q_w];        
        [ Vo,Eo ] = object_dipole_effect(params,X_all,Q_all,X0);
%         [ X_om,Q_,flag] = mirror_wall([x y],0,params); %wall object
%         image- LET'S WORRY ABOUT THIS LATER
        
        [X_sw,q,flag]= mirror_wall(X_s,zeros(size(X_s,1),1),params); % wall fish image

        if(flag)
            [v1,e1]=get_potential_field(X0,X_all,Q_all);
            Z(i,j)=(e1+Eo)*N0';
        end
    end
end    
Z=Z/z0-1;
% Z = Trans2LFP(Z,params);
if(1)
    F=figure;
    A=axes;
    COL=colormap('lines');
    S=surf(X,Y,zeros(size(Z)),Z,'LineStyle','none','FaceColor','interp');
    view(0,90);

    hold on;
    [x]=get_skin(params);
    plot_skin(F,x,params);
    plot_poles(F,X_p,Q_p,params);

    colormap('parula');
    set(gca,'CLim',[-.5 .5]);
    set(F,'Color',bgcol);
    set(A,'Color',bgcol,'XColor',bgcol,'YColor',bgcol);
end
end

