function [out_strct]=model_tail_conductive_wall(dbname)
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
M=35; %35
params.grid_M=M;
params.msize=18;
params.mpos='';
params.mneg='';
bgcol=1*[1 1 1];
params.bgcol=bgcol;

% nettype='BP';
% nhidden=100;
nettype='GC';
nhidden=1000;

params.plotting=0;
%% generate data
exp_training=~nargin;
if(exp_training)
    K_tail=25;
    K_phi=25;

    phi_v=linspace(0,pi,K_phi);
    angle_v=linspace(-.3*pi,.3*pi,K_tail)';
    K_total=M*(M+1)*K_tail;
    position=nan(K_total,2);
    tail=nan(K_total,1);
    data=nan(K_total,K_phi);
    data_cnd=nan(K_total,K_phi);
    k=1;
    for i=1:K_tail
        params.tail_angle=-angle_v(i);    
        Z=nan(M*(M+1),K_phi);    
        Z_cnd=nan(M*(M+1),K_phi);    
        for j=1:K_phi        
            params.phi0=phi_v(j);              
            
            %generate insulating wall data
            [z,R,TH] = get_wall_effect_polar(params);             
            Z(:,j)=z(:);

            %generate conductive wall data
            params_cond=params;
            params_cond.wall_c=-1;
            [z_cnd,rc,thc] = get_wall_effect_polar(params_cond);             
            Z_cnd(:,j)=z_cnd(:);
            
        end

        position(k:(k+size(Z,1)-1),:)=[R(:),TH(:)];
        tail(k:(k+size(Z,1)-1),:)=ones(M*(M+1),1)*angle_v(i);
        data(k:(k+size(Z,1)-1),:)=Z;
        data_cnd(k:(k+size(Z_cnd,1)-1),:)=Z_cnd;
        k=k+size(Z,1);
        if(params.tail_angle==0) %get 'exafference'
            data0=Z;
        end
    end
    
    %standartize
    dm=nanmean(data,1);
    ds=nanstd(data,1);
    dd=size(data_cnd,1);
    data=nanzscore(data,1);
    data_cnd=(data_cnd-ones(dd,1)*dm)./(ones(dd,1)*ds);
    tail=nanzscore(tail,1);
    data0=nanzscore(data0,1);
        
    % train example
    iz=19;
%     iz=28;
    z=data(:,iz);
    z_cnd=data_cnd(:,iz);
    zex=repmat(data0(:,iz),K_tail,1);
    zre=data(:,iz)-zex;

    %shuffle data (maintaining context)
    II=tail_shuffle(K_tail,M*(M+1));
    fb_data=data(II,:);
    fb_data_cnd=data_cnd(II,:);
%     fb_data=[fb_data tail(II)];
%     fb_data=repmat(data0,K_tail,1); %filtered fb
    
    ind=find(all(~isnan([position tail fb_data zex]),2));
    indnan=find(any(isnan([position tail fb_data zex]),2));
    
    %train (using insulating wall data)
    inputs=[tail fb_data];
    inputs_cnd=[tail fb_data_cnd];
    [ net,ex_model(2).rsq,ex_model(2).sig ] = ni_network_new(inputs(ind,:)',zre(ind)',nettype,[],nhidden);
    
    %test (insulating wall)
    z_ni=-net(inputs')';
    z_ni(indnan)=nan;
    [ex_model(2).ni_Re,ex_model(2).ni_Ex,ex_model(2).N]=get_maps(tail,z_ni);
    z_out=z+z_ni;
    [ex_model(2).out_Re,ex_model(2).out_Ex]=get_maps(tail,z_out);
    ex_model(2).total_re_flat=nansum(abs(ex_model(2).out_Re(:)));
    ex_model(2).var_re_flat=nanvar(ex_model(2).out_Re(:));
    ex_model(2).total_re_weighted=nansum(abs(ex_model(2).out_Re(:).*ex_model(2).N(:)))/nansum(ex_model(2).N(:));    

    
    [ex_model(1).out_Re,ex_model(1).out_Ex,ex_model(1).N]=get_maps(tail,z);
    ex_model(1).total_re_flat=nansum(abs(ex_model(1).out_Re(:)));
    ex_model(1).var_re_flat=nanvar(ex_model(1).out_Re(:));
    ex_model(1).total_re_weighted=nansum(abs(ex_model(1).out_Re(:).*ex_model(1).N(:)))/nansum(ex_model(1).N(:));    

    %test (conductive wall)
    z_ni=-net(inputs_cnd')';
    z_ni(indnan)=nan;
    [ex_model(4).ni_Re,ex_model(4).ni_Ex,ex_model(4).N]=get_maps(tail,z_ni);
    z_out=z_cnd+z_ni;
    [ex_model(4).out_Re,ex_model(4).out_Ex]=get_maps(tail,z_out);
    ex_model(4).total_re_flat=nansum(abs(ex_model(4).out_Re(:)));
    ex_model(4).var_re_flat=nanvar(ex_model(4).out_Re(:));
    ex_model(4).total_re_weighted=nansum(abs(ex_model(4).out_Re(:).*ex_model(4).N(:)))/nansum(ex_model(4).N(:));    

    [ex_model(3).out_Re,ex_model(3).out_Ex,ex_model(3).N]=get_maps(tail,z_cnd);
    ex_model(3).total_re_flat=nansum(abs(ex_model(3).out_Re(:)));
    ex_model(3).var_re_flat=nanvar(ex_model(3).out_Re(:));
    ex_model(3).total_re_weighted=nansum(abs(ex_model(3).out_Re(:).*ex_model(3).N(:)))/nansum(ex_model(3).N(:));    
    
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

%% plot
nline=0.25;
wline=0.5;
bwidth=4;
msize=.5;
fontsize=7;
% K=size(out_strct.ex_model,2);

Fall=figure;
BCOL=brewermap(12,'Paired');
COL=[0.5*[1 1 1];.1216 .4706 .7059;[151 93 166;216 84 39;127 63 152]/255];

set(Fall,'Units','centimeters','Position',[0 1 17 6.8]);
% [h0,pos0]=tight_subplot(2,1,[0 0],[0 0],[0 0],[.67 .33],[]);%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[ha,posa]=tight_subplot(3,2,[0 0],[0.05 0],[0.05 0],[],[]);%,h0(1));
plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(1).out_Re,out_strct.ex_model(1).N,ha(1));
plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(3).out_Re,out_strct.ex_model(3).N,ha(2));
for k=1:2
%     plot_map(X,Y,model(k).out_P1,ha(k+1));
%     plot_map(X,Y,model(k).ni_P1,ha(k+numel(inputs)+2));
    plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(2*k).out_Re,out_strct.ex_model(2*k).N,ha(k+2));
    plot_map_polar(out_strct.R,out_strct.TH,out_strct.ex_model(2*k).ni_Re,out_strct.ex_model(2*k).N,ha(k+4));
end
set(ha(:),'CLim',[-.5 .5]);
if(bgcol==[1 1 1])
    set(ha(:),'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(ha(:),'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end


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
end