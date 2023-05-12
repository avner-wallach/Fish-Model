function [ Xout,Qout ] = mirror_object(Xin,Qin,params)
% MIRROR_WALL mirror the point Xin due to wall 
%% params
circular=params.circular;
p0=params.object_x;  %coordinates of object center
c=params.object_c;
R=params.object_R;
Xout=zeros(size(Xin));
Qout=zeros(size(Qin));
%% 
    
%go over all points
for i=1:size(Xin,1)
    d=norm(Xin(i,:)-p0); %distance of pole from object center
    d_= R^2/d;  %distance of image from object center
    
    alpha=atan2(Xin(i,2)-p0(2),Xin(i,1)-p0(1)); %azimouth of pole relative to center   
    Xout(i,:)=[p0(1)+d_*cos(alpha),p0(2)+d_*sin(alpha)];
    Qout(i)=c*Qin(i)*R/d;
end

end

