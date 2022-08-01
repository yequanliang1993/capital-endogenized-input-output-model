%% Data preparations: exiobase data

IOpath = 'D:\\MRIO\\Exiobase\\v3_6\\IOT_%d_pxp.mat';
r=49;
s=200;
ts=21;
S = zeros(8,ts,r*s);
F_hh = zeros(8,ts,r);
for n=1:ts
    IOfile=sprintf(IOpath,n+1994);
    load(IOfile);
    L = inv(eye(r*s)-IO.A);
    fname = sprintf('L%d',n);
    save(fname,'L');
    x(n).x = IO.x;
    for i=1:r
        Yc(n).Yc(:,i) = sum(IO.Y(:,(i-1)*7+1:(i-1)*7+3),2);
        Yk(n).Yk(:,i) = IO.Y(:,(i-1)*7+4);
        Yt(n).Yt(:,i) = sum(IO.Y(:,(i-1)*7+1:i*7),2);
    end
    
    if n==1
        index_eng=find(contains(meta.secLabsF(:,1),'Energy Carrier Use'));
        index_bwc=find(contains(meta.secLabsF(:,1),'Water Consumption Blue'));
        index_co2=find(contains(meta.secLabsF(:,1),'CO2'));
        index_met=find(contains(meta.secLabsF(:,1),'Domestic Extraction Used - Metal Ores'));
        index_land=find(contains(meta.secLabsF(:,3),'km2'));
        index_mat =find(contains(meta.secLabsF(:,1),'Domestic Extraction Used'));
        index_min = find(contains(meta.secLabsF(:,1),'Domestic Extraction Used - Non-Metallic Minerals'));
        index_co2e=find(contains(meta.secLabsF(:,1),'CO2') | contains(meta.secLabsF(:,3),'CO2'));
        index_ch4=find(contains(meta.secLabsF(:,1),'CH4'));
        index_n2o=find(contains(meta.secLabsF(:,1),'N2O'));
    end

    S(1,n,:)=sum(IO.S(index_eng,:));
    S(2,n,:)=sum(IO.S(index_bwc,:));
    S(3,n,:)=sum(IO.S(index_land,:));
    S(4,n,:)=sum(IO.S(index_met,:));
    S(5,n,:)=sum(IO.S(index_min,:));
    S(6,n,:)=sum(IO.S(index_co2e,:))+sum(IO.S(index_ch4,:)).*28+sum(IO.S(index_n2o,:)).*265;
    S(7,n,:)=sum(IO.S(index_co2,:));
    S(8,n,:)=sum(IO.S(index_mat,:));
    
    for i=1:r
        F_hh(1,n,i)=sum(sum(IO.F_hh(index_eng,(i-1)*7+1:i*7)));
        F_hh(2,n,i)=sum(sum(IO.F_hh(index_bwc,(i-1)*7+1:i*7)));
        F_hh(3,n,i)=sum(sum(IO.F_hh(index_land,(i-1)*7+1:i*7)));
        F_hh(4,n,i)=sum(sum(IO.F_hh(index_met,(i-1)*7+1:i*7)));
        F_hh(5,n,i)=sum(sum(IO.F_hh(index_min,(i-1)*7+1:i*7)));
        F_hh(6,n,i)=sum(sum(IO.F_hh(index_co2e,(i-1)*7+1:i*7)))+sum(sum(IO.F_hh(index_ch4,(i-1)*7+1:i*7))).*28+sum(sum(IO.F_hh(index_n2o,(i-1)*7+1:i*7))).*265;
        F_hh(7,n,i)=sum(sum(IO.F_hh(index_co2,(i-1)*7+1:i*7)));
        F_hh(8,n,i)=sum(sum(IO.F_hh(index_mat,(i-1)*7+1:i*7)));
    end
    n
end

save('exioS.mat','S');
save('Yc.mat','Yc');
save('Yk.mat','Yk');
save('Yt.mat','Yt');
save('x.mat','x');

%% capital goods producers

% This is used to inform the making of the concordances - who produce capital goods.

kproducer = zeros(s,1);
for n=1:ts
    for i=1:r
        temp = sum(Yk(n).Yk((i-1)*s+1:i*s,:),2)./sum(sum(Yk(n).Yk((i-1)*s+1:i*s,:),2));
        kproducer = kproducer+(temp>0);
    end
end
kproducer=kproducer';
save('kproducer.txt','kproducer','-ascii');
clear temp
