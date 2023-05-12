function [ Z_V,Z_E,X,Y ] = get_potential_field_map(Xp,Qp,params)
%% parameters
r_max=params.r_max;
grid_center=params.grid_center;
M=params.grid_M;

%%
x_v=linspace(-r_max,r_max,M)+grid_center(1);
y_v=linspace(-r_max,r_max,M+1)+grid_center(2);
Z_V=nan(M,M+1);
Z_E=nan(M,M+1,2);
[X,Y]=meshgrid(x_v,y_v);
X=X';
Y=Y';
for i=1:numel(x_v)
    x=x_v(i);
    for j=1:numel(y_v)
        y=y_v(j);
%         r=norm([x y]);
%         phi=atan2(y,x);
%         if(phi<0)
%             phi=phi+2*pi;
%         end
%         [x0,n0]=get_skin_polar(phi);
%         if(r>norm(x0))
%             [Z_V(i,j),E]=get_potential_field([x y],Xp,Qp);
            [V,E]=get_potential_field([x y],Xp,Qp); %filed and potential due to poles
            [Vo,Eo]=object_dipole_effect(params,Xp,Qp,[x y]); %filed and potential due to object
            Z_V(i,j)=V+Vo;
            Z_E(i,j,1)=E(1)+Eo(1);
            Z_E(i,j,2)=E(2)+Eo(2);
%         end
    end
end    

end

