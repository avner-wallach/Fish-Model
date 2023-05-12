function [ind]=out_wall(X,params)
tank_radius=params.tank_radius;
circular=params.circular;
wall_dist=params.wall_dist;
wall_angle=params.wall_angle;

%% 
    if(~circular)
        R=sqrt(X(:,1).^2+X(:,2).^2);
        phi=atan2(X(:,2),X(:,1));
        Rtan=R.*cos(phi-wall_angle);
        ind=find(Rtan>=wall_dist);
    else        
        %find tank center coordinates
        p0=-[(tank_radius-wall_dist)*cos(wall_angle) (tank_radius-wall_dist)*sin(wall_angle)]; 
        R=sqrt((X(:,1)-p0(1)).^2+(X(:,2)-p0(2)).^2);
        ind=find(R>=tank_radius);

    end
        
end


