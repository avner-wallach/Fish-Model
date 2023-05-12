function [ Xout,Qout,flag ] = mirror_wall(Xin,Qin,params)
% MIRROR_WALL mirror the point Xin due to wall 
%% params
circular=params.circular;
wall_dist=params.wall_dist;
wall_angle=params.wall_angle;
wall_c=params.wall_c;
flag=1;
Xout=zeros(size(Xin));
Qout=zeros(size(Qin));
%% 
if(~circular)   %flat wall
    r=sqrt(Xin(:,1).^2+Xin(:,2).^2);
    alpha=atan2(Xin(:,2),Xin(:,1));
    %check pints do not cross wall
    if(sum(r.*cos(alpha-wall_angle)>=wall_dist)>0) %one or more points cross wall
       flag=0;
    end   
    R=sqrt(sum(Xin.^2,2)); %distance from head of each point
    beta=atan2(Xin(:,2),Xin(:,1));
    gamma=wall_angle-beta-pi;
    D=2*(R.*cos(gamma)+wall_dist);       %distance of points from their image
    Xout=Xin+D*[cos(wall_angle) sin(wall_angle)];
    Qout=wall_c*Qin;
else    %circular wall
    tank_radius=params.tank_radius; %cm        
    p0=-[(tank_radius-wall_dist)*cos(wall_angle) (tank_radius-wall_dist)*sin(wall_angle)]; 
    
    %go over all points
    for i=1:size(Xin,1)
        d=norm(Xin(i,:)-p0); %distance of pole from tank center
        if(d>=tank_radius | wall_dist>=tank_radius)
            flag=0;
        end
        alpha=atan2(Xin(i,2)-p0(2),Xin(i,1)-p0(1)); %azimouth of pole relative to center
        L=tank_radius^2/d;
        Xout(i,:)=[p0(1)+L*cos(alpha),p0(2)+L*sin(alpha)];
        if(wall_c>0)
            Qout(i)=wall_c*Qin(i)*(tank_radius/d);%.^2;
        else
            Qout(i)=wall_c*Qin(i)*tank_radius/d;
        end
    end
end

end

