function [ V,E ] = get_potential_field(X,X_p,Q_p)
%GET_POTENTIAL_FIELD get the potential and field at point X due to pole
%collection X_p with charge Q_p

%% potential
XX=ones(size(X_p,1),1)*X-X_p;
D=sqrt(XX(:,1).^2+XX(:,2).^2);
V=sum(1./D.*Q_p);
%% field
D=D.^3;
XX=XX./(D*[1 1]);
E=sum(XX.*(Q_p*[1 1]));

end

