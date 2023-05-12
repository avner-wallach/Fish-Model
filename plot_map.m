function plot_map(varargin)
% bgcol=[1 1 1];
params.bgcol=[1 1 1];
nline=0.25;
wline=0.5;
arrowsize=1;

%fish parameters
params.fish_length=17; %cm
params.fish_width=2;%cm
params.tail_angle=-pi/7;%rad
params.p_density=20;     %poles/cm
params.m=1;              %number of negative poles
params.tail_p=0.45; %non-bending portion of body 
params.phi0=.75*pi; %RF location

%wall paramteres
params.tank_radius=23;
params.circular=1;
params.wall_dist=5; %20 for image
params.wall_angle=pi*.5;
% params.wall_angle=pi*1.2; %for image
params.wall_c=1;

%object parameters
% params.object_x=[2.5 -2.5];
params.object_x=-(params.tank_radius-params.wall_dist)*[cos(params.wall_angle) sin(params.wall_angle)];
params.object_R=0;
params.object_c=0;

%local field parameters
params.z0=0.0753;
params.scale=6000;

%ploting parameters
% r_max=20;
params.r_max=20;
params.grid_center=[0 -params.fish_length/2];
% params.grid_center=params.object_x;
params.grid_M=30; %300 for good resolution
params.msize=25;
params.mpos='+';
params.mneg='o';
params.plot_potential=1;
params.plot_field=0;
params.plot_lfield=0;
params.reflection=0;
%% load input args
i=1;
while(i<nargin)
    params.(varargin{i})=varargin{i+1};
    i=i+2;
end    
%% 
% generate poles
[ X_p,Q_p ] = get_fish_poles( params );
X_all=[X_p];
Q_all=[Q_p];

% generate skin
[X_s] = get_skin(params);

if(params.wall_dist>0)
    [ X_w,Q_w,flag] = mirror_wall(X_p,Q_p,params); %wall pole image
    [X_sw,q,flag]= mirror_wall(X_s,zeros(size(X_s,1),1),params); % wall fish image
    X_all=[X_all;X_w];
    Q_all=[Q_all;Q_w];
end

%get potential and field
[ Z_V,Z_E,X,Y ] = get_potential_field_map(X_all,Q_all,params);
[X0,N0]=get_skin_polar(params.phi0,params);
[ v0,e0 ] = get_potential_field(X0,X_all,Q_all);
[vo,eo]=object_dipole_effect(params,X_all,Q_all,X0); %filed and potential due to object
e0=e0+eo;
z0=e0*N0'

%% plot
if(~isfield(params,'axes'))
    F=figure;
    A=axes;
else
    A=params.axes;
    axes(params.axes);
    A.clo;
end

if(params.wall_dist>0)
    plot_wall(params);
end
hold on;

if(params.reflection)
    H=plot_skin(X_sw,params);
    set(H,'FaceAlpha',0.25);
    if(~strcmp(params.mpos,''))
        H=plot_poles(X_w,Q_w,params);
        set(H,'MarkerFaceAlpha',0.25);
    end
end
plot_skin(X_s,params);
hold on;
if(~strcmp(params.mpos,''))
    plot_poles(X_p,Q_p,params);
end
if(params.object_c~=0)
    for k=1:numel(params.object_c)
        R=rectangle('Position',[params.object_x(k,:)-params.object_R(k),[2 2]*params.object_R(k)],'Curvature',[1 1]);
        R.EdgeColor=[1 1 1];
        R.LineWidth=2;
    end
end
%% get potential and field

if(params.plot_field)
    [Xf] = get_skin(params);
    p2=params;
    p2.fish_length=p2.fish_length+5;
    p2.fish_width=p2.fish_width+1.5;
    p2.tail_p=params.fish_length/p2.fish_length*params.tail_p;
    [Xs0] = get_skin(p2);
    if(~params.circular)
        x0=linspace(-params.r_max,params.r_max,5)'+params.grid_center(1);
        y0=-params.r_max*ones(size(x0))+params.grid_center(2);
    else
        phi=linspace(5/4*pi,7/4*pi,10)';    
        p0=-(params.tank_radius-params.wall_dist)*[cos(params.wall_angle) sin(params.wall_angle)]; 

        x0=params.tank_radius*cos(phi)+p0(1);
        y0=params.tank_radius*sin(phi)+p0(2);
    end
    Xs=Xs0(Xs0(:,2)>=min(X_p(:,2)),:);
    Xs=[Xs(1:10:end,:);x0 y0]; 
%     Xs=[Xs(1:5:end,:)];
%     if(~params.reflection)
%         ind=out_wall([X(:),Y(:)],params);
%         z1=Z_E(:,:,1); z2=Z_E(:,:,2); 
%         z1(ind)=nan; z2(ind)=nan; 
%         Z_E(:,:,1)=z1;  Z_E(:,:,2)=z2;         
%     end    
    S=streamline(X',Y',Z_E(:,:,1)',Z_E(:,:,2)',Xs(:,1),Xs(:,2));
%     S=[S;streamline(X',Y',-Z_E(:,:,1)',-Z_E(:,:,2)',Xs(:,1),Xs(:,2))];
    set(S,'LineWidth',nline);
    for i=1:numel(S)  
        ind1=out_wall([S(i).XData' S(i).YData'],params);    
        ind2=find(inpolygon(S(i).XData,S(i).YData,Xs0(:,1),Xs0(:,2)));
        ind=[ind1(:);ind2(:)];
        S(i).XData(ind)=nan;
        S(i).YData(ind)=nan;
    end
    set(S,'Color',.5*[1 1 1]);
end

if(params.plot_potential)
    c=quantile(Z_V(:),linspace(0.02,0.98,22));
%     if(~params.reflection)
    ind=out_wall([X(:),Y(:)],params);
    Z_V(ind)=nan;
%     end
    C=contour(X,Y,Z_V,c,'LineWidth',nline);
    axis([-params.r_max+params.grid_center(1) params.r_max+params.grid_center(1)...
        -params.r_max+params.grid_center(2) params.r_max+params.grid_center(2)]);
    set(A,'Clim',[-0.05 0.05])
end

if(params.plot_lfield)
    z=(z0-params.z0)*params.scale+1;
    ah = annotation('arrow',...
        'headStyle','vback2','HeadLength',arrowsize,'HeadWidth',arrowsize);
    set(ah,'parent',A);
    set(ah,'position',[X0(1) X0(2) ,N0(1)*z,N0(2)*z],'LineWidth',wline,'Color',[236 157 118]/255);    
%     quiver(X0(1),X0(2),N0(1)*z,N0(2)*z,'MaxHeadSize',5,'LineWidth',3,'Color',[236 157 118]/255);
    
end

axis([-params.r_max*.875+params.grid_center(1) params.r_max*.875+params.grid_center(1)...
    -params.r_max+params.grid_center(2) params.r_max+params.grid_center(2)]);
set(A,'Clim',[-0.05 0.05])
% set(F,'Color',params.bgcol);
set(A,'Color',params.bgcol,'XColor','none','YColor','none');
colormap(A,brewermap(12,'BrBG'));
hold off;
% set(gcf,'Position',[380 100 860 860]);
% set(gca,'Position',[.11 .11 .815 .815]);
end

