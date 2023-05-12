function [out_strct]=model_tail_wall(dbname)
%% params
% params.phi0=.75*pi;    

%fish parameters
params.fish_length=17; %cm 17
params.fish_width=1.5;%cm
params.tail_angle=0;%rad
params.p_density=20;     %poles/cm
params.m=1;              %number of negative poles
params.tail_p=0.55; %non-bending portion of body 

%coordinate system
params.coordinates=0; %1=tank center, 0=tank wall

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
r_max=params.tank_radius;
params.r_max=r_max;
M=35;
params.grid_M=M;
params.msize=18;
params.mpos='';
params.mneg='';
bgcol=0*[1 1 1];
params.bgcol=bgcol;

% nettype='BP';
% nhidden=100;
nettype='GC';
nhidden=5000;

params.plotting=0;
%% generate data
exp_training=~nargin;
if(exp_training)
    K_tail=25;
    K_phi=25;
%     K_phi=37;
%     nhidden=300;

    phi_v=linspace(0,pi,K_phi);
    angle_v=linspace(-.3*pi,.3*pi,K_tail)';
    K_total=M*(M+1)*K_tail;
    position=nan(K_total,2);
    tail=nan(K_total,1);
    data=nan(K_total,K_phi);
    k=1;
    for i=1:K_tail
        params.tail_angle=-angle_v(i);    
        Z=nan(M*(M+1),K_phi);    
        for j=1:K_phi        
            params.phi0=phi_v(j);          
    %         [z,X,Y] = get_wall_effect(params); 
            [z,R,TH] = get_wall_effect_polar(params); 
            Z(:,j)=z(:);
        end

        position(k:(k+size(Z,1)-1),:)=[R(:),TH(:)];
        tail(k:(k+size(Z,1)-1),:)=ones(M*(M+1),1)*angle_v(i);
        data(k:(k+size(Z,1)-1),:)=Z;
        k=k+size(Z,1);
        if(params.tail_angle==0) %get 'exafference'
            data0=Z;
        end
    end
    
    %standartize
    data=nanzscore(data,1);
    tail=nanzscore(tail,1);
    data0=nanzscore(data0,1);
        
    % train example
    iz=19;
%     iz=28;
    z=data(:,iz);
    zex=repmat(data0(:,iz),K_tail,1);
    zre=data(:,iz)-zex;

    %shuffle data (maintaining context)
    II=tail_shuffle(K_tail,M*(M+1));
    fb_data=data(II,:);
%     fb_data=[fb_data tail(II)];
%     fb_data=repmat(data0,K_tail,1); %filtered fb
    
    ind=find(all(~isnan([position tail fb_data zex]),2));
    indnan=find(any(isnan([position tail fb_data zex]),2));

    inputs={[tail],[tail position],fb_data,[tail fb_data]};
    for k=1:numel(inputs)
        [ net,ex_model(k+1).rsq,ex_model(k+1).sig ] = ni_network_new(inputs{k}(ind,:)',zre(ind)',nettype,[],nhidden);
        z_ni=-net(inputs{k}')';
        z_ni(indnan)=nan;
        [ex_model(k+1).ni_Re,ex_model(k+1).ni_Ex,ex_model(k+1).N]=get_maps(tail,z_ni);
        z_out=z+z_ni;
        [ex_model(k+1).out_Re,ex_model(k+1).out_Ex]=get_maps(tail,z_out);
        ex_model(k+1).total_re_flat=nansum(abs(ex_model(k+1).out_Re(:)));
        ex_model(k+1).var_re_flat=nanvar(ex_model(k+1).out_Re(:));
        ex_model(k+1).total_re_weighted=nansum(abs(ex_model(k+1).out_Re(:).*ex_model(k+1).N(:)))/nansum(ex_model(k+1).N(:));    
    end

    [ex_model(1).out_Re,ex_model(1).out_Ex,ex_model(1).N]=get_maps(tail,z);
    ex_model(1).total_re_flat=nansum(abs(ex_model(1).out_Re(:)));
    ex_model(1).var_re_flat=nanvar(ex_model(1).out_Re(:));
    ex_model(1).total_re_weighted=nansum(abs(ex_model(1).out_Re(:).*ex_model(1).N(:)))/nansum(ex_model(1).N(:));    

    out_strct.ex_model=ex_model;
    out_strct.R=R;
    out_strct.TH=TH;
    save('temp_model.mat');
else
    D=load(dbname);
    out_strct.ex_model=D.strct1.ex_model;
    out_strct.R=D.strct1.R;
    out_strct.TH=D.strct1.TH;    
end
%% train models
db_training=~nargin;
if(db_training)
    for iz=1:K_phi % go over all RFs
        iz
        z=data(:,iz);
        zex=repmat(data0(:,iz),K_tail,1);
        zre=data(:,iz)-zex;

        %shuffle data (maintaining context)
        II=tail_shuffle(K_tail,M*(M+1));
        fb_data=data(II,:);
%         fb_data=[fb_data tail(II)];
        
        ind=find(all(~isnan([position tail fb_data zex]),2));
        indnan=find(any(isnan([position tail fb_data zex]),2));

        inputs={[tail],[tail position],[fb_data],[tail fb_data]};
        names={'mot','motloc','fb','motfb'};
        for k=1:numel(inputs)
            [ net,model(iz,k+1).rsq,model(iz,k+1).sig ] = ni_network_new(inputs{k}(ind,:)',zre(ind)',nettype,[],nhidden);
            z_ni=-net(inputs{k}')';
            z_ni(indnan)=nan;
            [model(iz,k+1).ni_Re,model(iz,k+1).ni_Ex,model(iz,k+1).N]=get_maps(tail,z_ni);
            z_out=z+z_ni;
            [model(iz,k+1).out_Re,model(iz,k+1).out_Ex]=get_maps(tail,z_out);
            model(iz,k+1).total_re_flat=nansum(abs(model(iz,k+1).out_Re(:)));
            model(iz,k+1).var_re_flat=nanvar(model(iz,k+1).out_Re(:));
            model(iz,k+1).total_re_weighted=nansum(abs(model(iz,k+1).out_Re(:).*model(iz,k+1).N(:)))/nansum(model(iz,k+1).N(:));
            model(iz,k+1).name=names{k};
        end

        [model(iz,1).out_Re,model(iz,1).out_Ex,model(iz,1).N]=get_maps(tail,z);
        model(iz,1).total_re_flat=nansum(abs(model(iz,1).out_Re(:)));
        model(iz,1).var_re_flat=nanvar(model(iz,1).out_Re(:));
        model(iz,1).total_re_weighted=nansum(abs(model(iz,1).out_Re(:).*model(iz,1).N(:)))/nansum(model(iz,1).N(:));    
        model(iz,1).name='input';
        save('temp_model.mat');
    end    
    out_strct.models=model;
    out_strct.nettype=nettype;
    out_strct.nhidden=nhidden;
    out_strct.grid_M=M;
    save('temp_model.mat');
else
    D=load(dbname);
    out_strct.models=D.strct1.models;
    out_strct.nettype=D.strct1.nettype;
    out_strct.nhidden=D.strct1.nhidden;
    out_strct.grid_M=D.strct1.grid_M;        
end
%% plot 
out_strct.ex_model=out_strct.ex_model([1 2 4 5]);%for presentation
out_strct.models=out_strct.models(:,[1 2 4 5]);

nline=0.25;
wline=0.5;
bwidth=4;
msize=.5;
fontsize=7;
N=size(out_strct.models,1);
K=size(out_strct.ex_model,2);

Fscat=figure;
set(Fscat,'Position',[564 749 200 164],'Color',[0 0 0]);
ascat=axes;
plot_slope_scatter(out_strct.ex_model(1).out_Re,out_strct.ex_model(1).N,1);
for k=1:K-1
    plot_slope_scatter(out_strct.ex_model(k+1).out_Re,out_strct.ex_model(k+1).N,k+1);
end
set(gca,'Xlim',[.5 K+.5],'Ylim',[-1 1],'Clim',[-.5 .5],'XAxisLocation','origin','XTick',[]);
if(bgcol==[1 1 1])
    set(gca,'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(gca,'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color',bgcol);        

Fni=figure;
set(Fni,'Position',[564 749 164 164],'Color',[0 0 0]);
ani=axes;
xx=out_strct.ex_model(1).out_Re;
yy=out_strct.ex_model(end).ni_Re;
nn=out_strct.ex_model(1).N;
xx=xx(nn>0); yy=yy(nn>0); nn=nn(nn>0);
nn=nn/max(nn)*5;
scatter(xx,yy,nn,xx,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
hold on;
plot([-.5 .5],[.5 -.5],'Color',.5*[1 1 1],'LineStyle','--');
set(gca,'Xlim',[-.5 .5],'Ylim',[-.5 .5],'Clim',[-.5 .5],'XAxisLocation','origin','XTick',[],'YAxisLocation','origin');
if(bgcol==[1 1 1])
    set(gca,'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(gca,'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color',bgcol);        


Fall=figure;
BCOL=brewermap(12,'Paired');
COL=[0.5*[1 1 1];.1216 .4706 .7059;[151 93 166;216 84 39;127 63 152]/255];

set(Fall,'Units','centimeters','Position',[0 1 17 6.8]);
% [h0,pos0]=tight_subplot(2,1,[0 0],[0 0],[0 0],[.67 .33],[]);%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[ha,posa]=tight_subplot(2,K,[0 0],[0.05 0],[0.05 0],[],[]);%,h0(1));
plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(1).out_Re,out_strct.ex_model(1).N,ha(1));
for k=1:K-1
%     plot_map(X,Y,model(k).out_P1,ha(k+1));
%     plot_map(X,Y,model(k).ni_P1,ha(k+numel(inputs)+2));
    plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(k+1).out_Re,out_strct.ex_model(k+1).N,ha(k+1));
    plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(k+1).ni_Re,out_strct.ex_model(k+1).N,ha(k+K+1));
end
set(ha(:),'CLim',[-.5 .5]);
if(bgcol==[1 1 1])
    set(ha(:),'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(ha(:),'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end


ha(K+1).Visible='off';
print('-clipboard','-dbitmap','-r720'); %copy maps
ha(K+1).Visible='on';

%stats
for i=1:N
    re_var_flat(i,:)=cell2mat({out_strct.models(i,:).var_re_flat});
    re_flat(i,:)=cell2mat({out_strct.models(i,:).total_re_flat});
end
re_flat=re_flat./(mean(re_flat(:,1)));
% re_var_flat=re_var_flat./(re_var_flat(:,1)*ones(1,K));

% [hb,posb]=tight_subplot(1,2,[0.05 0.05],[0.05 0.05],[0.05 0.05],[],[],h0(2));%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
X=reshape(re_flat,[],1);
G=reshape(ones(N,1)*[1:K],[],1);
% axes(hb(1));
axes(ha(K+1));
boxplot(X,G,'PlotStyle','compact','Colors',COL,'Symbol','.w');
set(gca,'YLim',[0 1.35],'YTickMode','auto','YTickLabelMode','auto','XTickLabel',[],'Xlim',[0 K]+.5);
set(gca,'TickLength',[.025 .025],'TickDir','both');
set(gca,'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off');
set(gca,'Position',posa{K+1}+[0.015 0.045 -0.06 -0.075]);
% X=reshape(re_var_flat,[],1);
% G=reshape(ones(N,1)*[1:K],[],1);
% axes(hb(2));
% boxplot(X,G,'PlotStyle','compact','Colors',COL,'Symbol','.w');
% set(gca,'YLim',[0 1.05],'YTickMode','auto','YTickLabelMode','auto','XTickLabel',[],'Xlim',[0 K]+.5);
% set(gca,'TickLength',[.025 .025],'TickDir','both');
% set(gca,'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off');
bx=findobj('Tag','Box');
set(bx,'LineWidth',bwidth);
ws=findobj('Tag','Whisker');
set(ws,'LineWidth',wline);
cro=findobj('Tag','MedianInner');
set(cro,'Marker','none');
cri=findobj('Tag','MedianOuter');
set(cri,'MarkerFaceColor',[1 1 1],'MarkerSize',3);
ol=findobj('Tag','Outliers');

srf=findobj('Type','Surface');
img=findobj('Type','Patch');
set(srf,'Visible','off');
set(img,'Visible','off');
% export to metafile
print('-clipboard','-dmeta');

    function [P1,P2,N]=get_maps(tail,zdata)
        I=reshape([1:numel(zdata)],M,M+1,K_tail);
        ZZ=zdata(I);
        TT=tail(I);
        P1=nan(size(ZZ,1),size(ZZ,2));
        P2=P1;
        N=zeros(size(ZZ,1),size(ZZ,2));
        for i=1:size(ZZ,1)
            for j=1:size(ZZ,2)
                zz=ZZ(i,j,:);
                zz=zz(:);
                tt=TT(i,j,:);
                tt=tt(:);                
                idx=find(~isnan(zz));
                if(numel(idx)>=3);
                    ft=fit(tt(idx),zz(idx),'poly1');
                    P1(i,j)=ft.p1;
                    P2(i,j)=ft.p2;
                else
                    P1(i,j)=nan;
                    P2(i,j)=nan;
                end
                N(i,j)=numel(idx);
            end    
        end
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
        set(gca,'CLim',[-.1 .1],'Alim',[0 1.25],'XDir','reverse');
    %     set(gcf,'Color',bgcol);
        set(gca,'Color',bgcol,'XColor','none','YColor','none','XGrid','off','YGrid','off');    
        axis('image');
    end

    function plot_slope_scatter(zz,nn,a)
        m=zz(nn>0);
        n=nn(nn>0);
        n=n/max(n)*5;

        S=scatter(a+randn(size(m))*.05,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
        hold on;
        
    end
end