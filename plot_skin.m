function H=plot_skin(X_s,params)
%PLOT_SKIN plot skin
msize=params.msize;
mpos=params.mpos;
mneg=params.mneg;
bgcol=params.bgcol;

[X] = get_skin(params);
% H=plot(X_s(:,1),X_s(:,2));
% set(H,'Color',1-bgcol,'LineWidth',3);
H=patch(X_s(:,1),X_s(:,2),.8-.6*bgcol);
% set(H,'Color',1-bgcol,'LineWidth',3);
set(H,'LineStyle','none');
end

