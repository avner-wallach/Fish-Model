function [ Z,R,TH] = get_wall_effect_polar(params)
%GET_WALL_EFFECT get the effect of wall in different positions
%% parameters
r_max=params.r_max;
M=params.grid_M;
r_lim=[0 r_max];
th_lim=[-pi pi];

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
r_v=linspace(r_lim(1),r_lim(2),M);
th_v=linspace(th_lim(1),th_lim(2),M+1);
Z=nan(M,M+1);
R=Z;
TH=Z;
for i=1:numel(r_v)
    r=r_v(i);
    for j=1:numel(th_v)
        th=th_v(j);
        R(i,j)=r;
        TH(i,j)=th;
        if(params.coordinates)
            params.wall_dist=params.tank_radius-r;
            params.wall_angle=th-pi;            
        else           
            params.wall_dist=r;
            params.wall_angle=th;
        end
        
        [ X_w,Q_w,flag] = mirror_wall(X_p,Q_p,params); %wall pole image
        [X_sw,q,flag]= mirror_wall(X_s,zeros(size(X_s,1),1),params); % wall fish image
        X_all=[X_p;X_w];
        Q_all=[Q_p;Q_w];

        if(flag)
            [v1,e1]=get_potential_field(X0,X_all,Q_all);
            Z(i,j)=e1*N0';
        end
    end
end    
Z=Z/z0-1;
% Z = Trans2LFP(Z,params);
if(params.plotting)
%     F=figure;
%     A=axes;
    cla;
    S=surf(X,Y,zeros(size(Z)),Z,'LineStyle','none','FaceColor','interp');
    view(0,90);

    hold on;
    [x]=get_skin(params);
    plot_skin(x,params);
    plot_poles(X_p,Q_p,params);

    colormap(gca,'parula');
    set(gca,'CLim',[-.65 .25]);
    set(gcf,'Color',bgcol);
    set(gca,'Color',bgcol,'XColor','none','YColor','none','XGrid','off','YGrid','off');
end
end

