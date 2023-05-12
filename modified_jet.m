function COL=modified_jet
    a=linspace(0,.5,32)';
    d=linspace(0,1,32)';
    e=[d;flipud(d)];
    f=linspace(0.5,1,32)';
    b=[d;flipud(d)];
    c=[f;flipud(d)];
    COL=[flipud(c) b c];
end