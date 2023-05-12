function H=plot_poles(X_p,Q_p,params)
%PLOT_POLES plot all poles 
msize=params.msize;
mpos=params.mpos;
mneg=params.mneg;
mskip=20;
hold on;
COL=brewermap(12,'BrBG');

ind=find(Q_p>0);
% H=plot(X_p(ind(1:mskip:end),1),X_p(ind(1:mskip:end),2),mpos);
H1=scatter(X_p(ind(1:mskip:end),1),X_p(ind(1:mskip:end),2),msize,mpos);
H1.MarkerFaceColor=COL(10,:);
H1.MarkerEdgeColor=COL(10,:);
% H.MarkerSize=msize;

ind=find(Q_p<0);
% H=plot(X_p(ind(1:mskip:end),1),X_p(ind(1:mskip:end),2),mneg);
H2=scatter(X_p(ind(1:mskip:end),1),X_p(ind(1:mskip:end),2),msize,mneg);
H2.MarkerFaceColor=COL(3,:);
H2.MarkerEdgeColor=COL(3,:);
% H.MarkerSize=msize;
H=[H1(:);H2(:)]
end

