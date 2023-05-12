function plot_wall(params)
% colormap('parula');
bgcol=params.bgcol;
tank_radius=params.tank_radius;
circular=params.circular;
wall_dist=params.wall_dist;
wall_angle=params.wall_angle;
% wall_c=params.wall_c;
r_max=params.r_max;
%% plot wall as area
    if(~circular)
        if(wall_angle==0)
            J=area([wall_dist r_max],[r_max r_max],-r_max);
        elseif(wall_angle==pi)
            J=area([-r_max -wall_dist],[r_max r_max],-r_max);
        else
            b=wall_dist/sin(wall_angle);
            a=-1/tan(wall_angle);
            x1=(-r_max-b)/a;
            x2=(r_max-b)/a;
            y1=-a*r_max+b;
            y2=a*r_max+b;
            xx=sortrows([-r_max x1 x2 r_max;y1 -r_max r_max y2]');
            if(wall_angle<pi)
                J=area(xx(:,1),xx(:,2),r_max);
            else
                J=area(xx(:,1),xx(:,2),-r_max);
            end
        end
        J.FaceAlpha=1;
        J.FaceColor=0.25*[1 1 1];

    else
        
%         J=area([-r_max r_max]+params.grid_center(1),[r_max r_max]+params.grid_center(2),-r_max+params.grid_center(2));
%         J.FaceAlpha=1;
%         J.FaceColor=0.25*[1 1 1];
        
        %find tank center coordinates
        p0=-[(tank_radius-wall_dist)*cos(wall_angle) (tank_radius-wall_dist)*sin(wall_angle)]; 
        pos = [p0(1)-tank_radius p0(2)-tank_radius 2*tank_radius 2*tank_radius];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',.8*[1 1 1],'EdgeColor',1-bgcol,'LineWidth',.5);
        hold on;
%         H=plot(p0(1),p0(2),'.k');
        

    end
        
end


