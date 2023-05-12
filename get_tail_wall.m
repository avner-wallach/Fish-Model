function [F,Z,angle_v,R,TH,N,P1,P2]=get_tail_wall(axes_in,phi0,wall_c)
%% params
if(nargin<1 | numel(axes_in)==0)
    figure;
    A=tight_subplot(1,2);
else
    A(1)=axes_in(1);
    A(2)=axes_in(2);
end

if(nargin<2)
    params.phi0=.75*pi;    
else
    params.phi0=phi0;
end

if(nargin<3)
    params.wall_c=1;
else
    params.wall_c=wall_c;
end

%fish parameters
params.fish_length=19; %cm 23
params.fish_width=1.5;%cm
params.tail_angle=0;%rad
params.p_density=20;     %poles/cm
params.m=1;              %number of negative poles
params.tail_p=0.55; %non-bending portion of body 

K=51; %101
phimax=.4*pi; %.45*pi
angle_v=linspace(-phimax,phimax,K)';

%coordinate system
params.coordinates=0; %1=tank center, 0=tank wall

%wall paramteres
params.tank_radius=23;
params.circular=1;
params.wall_dist=1e4;
params.wall_angle=0;

%object parameters
params.object_x=[5 0];
params.object_R=0;
params.object_c=0.25;

% lfp paramters
params.th=0;
params.sig=.1;

%ploting parameters
r_max=params.tank_radius;
params.r_max=r_max;
params.grid_M=50;%50
params.msize=18;
params.mpos='';
params.mneg='';
bgcol=0*[1 1 1];
params.bgcol=bgcol;
%%
clim=[-.0325 0.0125];
params.plotting=1;
params.minsamps=15;
if(params.plotting)
    params.fig=figure;    
    set(params.fig,'Units','normalized','Color',bgcol);
    set(params.fig,'Position',[0 0 0.55 1]);
end
for i=1:numel(angle_v)
    params.tail_angle=-angle_v(i);
%     [ Z(:,:,i),X,Y ] = get_wall_effect(params);
    [ Z(:,:,i),R,TH ] = get_wall_effect_polar(params);
    if(params.plotting)
        set(gca,'Clim',clim);
        axis('image');
        drawnow;
        F(i)=getframe(gca);
    end

end

Z=(Z-nanmean(Z(:)))/nanstd(Z(:));
% angle_v=zscore(angle_v);

P1=nan(size(Z,1),size(Z,2));
P2=P1;

for i=1:size(Z,1)
    for j=1:size(Z,2)
        z=Z(i,j,:);
        z=z(:);
        idx=find(~isnan(z));
        if(numel(idx)>=params.minsamps);
            ft=fit(angle_v(idx)/phimax,z(idx),'poly1');
            P1(i,j)=ft.p1;
            P2(i,j)=ft.p2;
        else
            P1(i,j)=nan;
            P2(i,j)=nan;
        end
        N(i,j)=numel(idx);    
    end    
end

% plot_map(X,Y,P2,A(1));
% colormap(A(1),'Parula');
% plot_map(X,Y,P1,A(2));
% set(A(2),'CLim',[-.025 .025],'Colormap',flipud(brewermap(64,'RdBu')))

plot_map_polar(R,TH,P2,N,A(1));
colormap(A(1),'Parula');
plot_map_polar(R,TH,P1,N,A(2));
if(bgcol==[1 1 1])
    set(A(2),'CLim',[-.5 .5],'Colormap',flipud(brewermap(64,'RdBu')))
else
    set(A(2),'CLim',[-.5 .5],'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));
end

function plot_map(xx,yy,zz,a)
    axes(a);
    a.clo;
%     COL=colormap('lines');
    S=surf(xx,yy,zeros(size(zz)),zz,'LineStyle','none','FaceColor','interp');
    view(0,90);   
    hold on;    
    params.tail_angle=0;
    [x]=get_skin(params);
    plot_skin(x,params);
%     colormap(modified_jet);
    set(gca,'CLim',[-.1 .1]);
%     set(gcf,'Color',bgcol);
    set(gca,'Color',bgcol,'XColor','none','YColor','none','XGrid','off','YGrid','off');    
    axis('image');
end

    function plot_map_polar(rr,thth,zz,nn,a)
        axes(a);
        a.clo;
        Mzz=[zz zz(:,1)];
        AL=[nn nn(:,1)];

        S=polarplot3d(Mzz,AL,'AngularRange',[pi/2 5*pi/2],'plottype','surfa','RadialRange',[min(rr(:)) max(rr(:))],'AxisLocation','off','MeshScale',.5*[1 1]);       
        view(0,90);   
        hold on;    
        params.tail_angle=0;
        [x]=get_skin(params);
        plot_skin(x,params);
    %     colormap(modified_jet);
        set(gca,'CLim',[-1 1],'Alim',[0 5],'XDir','reverse');
    %     set(gcf,'Color',bgcol);
        set(gca,'Color',bgcol,'XColor','none','YColor','none','XGrid','off','YGrid','off');    
        axis('image');
    end

end