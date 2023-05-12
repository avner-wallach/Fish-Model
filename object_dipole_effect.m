function [ Vo,Eo ] = object_dipole_effect(params,X_p,Q_p,X)
%OBJECT_DIPOLE_EFFECT get the transdermal potential at location X due to dipole induced by object at location params.object_x 
for k=1:numel(params.object_c)
    Chi=params.object_c(k);
    Xo=params.object_x(k,:);
    Ro=params.object_R(k);

    % get potential due to object
    [ v0,e0 ] = get_potential_field(Xo,X_p,Q_p); %potential and field at object
    vo(k)=(Chi*Ro^3/norm(X-Xo)^3)*e0*(X-Xo)';

    % field at skin due to object
    ex=0; ey=0;
    for i=1:numel(Q_p) %go over all poles
        b=Xo(1)-X_p(i,1);
        c=Xo(2)-X_p(i,2);
        K=(Chi*Ro^3/norm(Xo-X_p(i,:))^3);
        ex=ex + K*Q_p(i)*(2*b*(X(1)-Xo(1))^2-b*(X(2)-Xo(2))^2+3*c*(X(1)-Xo(1))*(X(2)-Xo(2)))...
            /((X(1)-Xo(1))^2+(X(2)-Xo(2))^2)^(5/2);
        ey=ey + K*Q_p(i)*(2*c*(X(2)-Xo(2))^2-c*(X(1)-Xo(1))^2+3*b*(X(2)-Xo(2))*(X(1)-Xo(1)))...
            /((X(1)-Xo(1))^2+(X(2)-Xo(2))^2)^(5/2);
    end
    eo(k,:)=[ex ey];
end
Vo=sum(vo);
Eo=sum(eo);

end

