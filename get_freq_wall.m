function F=get_freq_wall()
%% params

%fish parameters
params.fish_length=16; %cm
params.fish_width=2;%cm
params.tail_angle=0;%rad
params.p_density=20;     %poles/cm
params.m=1;              %number of negative poles
params.tail_p=0.6; %non-bending portion of body 
params.phi0=.75*pi;

%wall paramteres
params.tank_radius=23;
params.circular=1;
params.wall_dist=1e4;
params.wall_angle=0;
params.wall_c=1;

%object parameters
params.object_x=[5 0];
params.object_R=0;
params.object_c=0.25;

% lfp paramters
params.th=0;
params.sig=.1;

%ploting parameters
r_max=30;
params.r_max=r_max;
params.grid_M=100;
params.msize=18;
params.mpos='+';
params.mneg='o';
bgcol=[0 0 0];
params.bgcol=bgcol;

%%
K=5;
th_v=linspace(-0.1,-0.05,K)';

plotting=1;
for i=1:numel(th_v)
    params.th=th_v(i);
    [ Z(:,:,i),X,Y,f ] = get_wall_effect(params);
    drawnow;
    F(i)=getframe(f);
    close(f);
end

for i=1:size(Z,1)
    for j=1:size(Z,2)
        z=Z(i,j,:);
        z=z(:);
        idx=find(~isnan(z));
        if(numel(idx)>=3);
            ft=fit(th_v(idx),z(idx),'poly1');
            P1(i,j)=ft.p1;
            P2(i,j)=ft.p2;
        else
            P1(i,j)=nan;
            P2(i,j)=nan;
        end
            
    end    
end

plot_map(X,Y,P1);
set(gca,'CLim',[-.025 .025])
plot_map(X,Y,P2);
colormap('parula');

function plot_map(xx,yy,zz)
    G=figure;
    COL=colormap('lines');
    S=surf(xx,yy,zeros(size(zz)),zz,'LineStyle','none','FaceColor','interp');
    view(0,90);   
    hold on;    
    hold on;
    [x]=get_skin(params);
    plot_skin(G,x,params);
    colormap(modified_jet);
    set(gca,'CLim',[-.1 .1]);
    set(G,'Color',bgcol);
    set(gca,'Color',bgcol,'XColor',bgcol,'YColor',bgcol);
    
end
end