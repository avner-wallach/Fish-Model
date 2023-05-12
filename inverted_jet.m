function COL=inverted_jet
    a=linspace(0,1,32)';
    b=zeros(32,1);
    c=[flipud(a);b];    
    d=zeros(64,1);
    COL=[flipud(c) d c];
end