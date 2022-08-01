clear
r=49;
ts=21;

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
[~,~,country49]=xlsread('country49.xlsx');
country=unique(country49(:,2),'stable');

% cd('D:\WP\CapitalImpacts\Datasets\KLEMS');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
filename = 'all_capital_09I.xlsx';
sheet = 'all_capital_09I';
KLEMS1=importdata(filename,sheet);
data1 = KLEMS1.data.all_capital_09I;
depr1 = KLEMS1.data.depreciationrate;

filename = 'all_capital_17I.xlsx';
sheet = 'all_capital_17I';
KLEMS2=importdata(filename);
data2 = KLEMS2.data.all_capital_17I;
depr2 = KLEMS2.data.depreciationrate;

labs1.country = KLEMS1.textdata.all_capital_09I(2:end,1);
labs1.var = KLEMS1.textdata.all_capital_09I(2:end,2);
labs1.asset = KLEMS1.textdata.all_capital_09I(2:end,3);
labs1.ind = KLEMS1.textdata.all_capital_09I(2:end,4);
labs1.year = KLEMS1.textdata.all_capital_09I(1,5:end);
var1=unique(labs1.var);
year1 = 1970:1:1970+length(labs1.year)-1;
asset1=KLEMS1.textdata.depreciationrate(1,3:end);

labs2.country = KLEMS2.textdata.all_capital_17I(2:end,1);
labs2.var = KLEMS2.textdata.all_capital_17I(2:end,2);
labs2.asset = KLEMS2.textdata.all_capital_17I(2:end,3);
labs2.ind = KLEMS2.textdata.all_capital_17I(2:end,4);
labs2.year = KLEMS2.textdata.all_capital_17I(1,5:end);
var2=unique(labs2.var);
year2 = 1970:1:1970+length(labs2.year)-1;
asset2=KLEMS2.textdata.depreciationrate(1,3:end);

%convert code name from numeric to string forms, e.g. '20'
filename = 'capital_09I_industrycode.xlsx';
[~,~,industry1] = xlsread('capital_09I_industrycode.xlsx');
for i=1:length(industry1(:,1))
    if isstr(industry1{i,1})==0
        industry1{i,1}=num2str(industry1{i,1});
    end
end
% convert code name from numeric to string forms, e.g. '19'
filename = 'capital_17I_industrycode.xlsx';
[~,~,industry2] = xlsread('capital_17I_industrycode.xlsx');
for i=1:length(industry2(:,1))
    if isstr(industry2{i,1})==0
        industry2{i,1}=num2str(industry2{i,1});
    end
end

clear sheet filename

% Indexing all variables
index1.ind = zeros(size(data1,1),1);
index1.country = zeros(size(data1,1),1);
index1.var = zeros(size(data1,1),1);
index1.asset = zeros(size(data1,1),1);
for i=1:size(data1,1)
    if any(strcmp(industry1(:,1),labs1.ind{i}))
        index1.ind(i) = find(strcmp(industry1(:,1),labs1.ind{i}));
    else
        index1.ind(i)=-9999;
    end   
    if any(strcmp(country,labs1.country{i}))
        index1.country(i) = find(strcmp(country,labs1.country{i}));
    else
        index1.country(i)=-9999; 
    end
    if any(strcmp(asset1,labs1.asset{i}))
        index1.asset(i) = find(strcmp(asset1,labs1.asset{i}));
    else
        index1.asset(i)=-9999;  
    end
    if any(strcmp(var1,labs1.var{i}))
        index1.var(i) = find(strcmp(var1,labs1.var{i}));
    else
        index1.var(i)=-9999;
    end
end

index2.ind = zeros(size(data2,1),1);
index2.country = zeros(size(data2,1),1);
index2.var = zeros(size(data2,1),1);
index2.asset = zeros(size(data2,1),1);
for i=1:size(data2,1)
    if any(strcmp(industry2(:,1),labs2.ind{i}))
        index2.ind(i) = find(strcmp(industry2(:,1),labs2.ind{i}));
    else
        index2.ind(i)=-9999;
    end
    if any(strcmp(country,labs2.country{i}))
        index2.country(i) = find(strcmp(country,labs2.country{i}));
    else
        index2.country(i)=-9999;
    end
    if any(strcmp(asset2,labs2.asset{i}))
        index2.asset(i) = find(strcmp(asset2,labs2.asset{i}));
    else
        index2.asset(i)=-9999; 
    end
    if any(strcmp(var2,labs2.var{i}))
        index2.var(i) = find(strcmp(var2,labs2.var{i}));
    else
        index2.var(i)=-9999;
    end
end

data1_0=data1;
index1_0=index1;
temp = [];
for i=1:size(data1,1)
    if any(strcmp(labs1.ind(1,1),labs1.ind{i})) % total industry only
        index1_0.ind(i) = find(strcmp(labs1.ind(1,1),labs1.ind{i}));
    else
        index1_0.ind(i)=-9999;
    end
end

data2_0=data2;
index2_0=index2;
for i=1:size(data2,1)
    if any(strcmp(labs2.ind(1,1),labs2.ind{i})) % total industry only
        index2_0.ind(i) = find(strcmp(labs2.ind(1,1),labs2.ind{i}));
    else
        index2_0.ind(i)=-9999;
    end
end

% redundant industry code in the file
index0 = find(index1.asset<0 | index1.ind<0);
data1(index0,:)=[];
index1.var(index0,:)=[];
index1.ind(index0,:)=[];
index1.asset(index0,:)=[];
index1.country(index0,:)=[];

index0 = find(index2.asset<0 | index2.ind<0);
data2(index0,:)=[];
index2.var(index0,:)=[];
index2.ind(index0,:)=[];
index2.asset(index0,:)=[];
index2.country(index0,:)=[];

index0 = find(index1_0.asset<0 | index1_0.ind<0);
data1_0(index0,:)=[];
index1_0.var(index0,:)=[];
index1_0.ind(index0,:)=[];
index1_0.asset(index0,:)=[];
index1_0.country(index0,:)=[];

index0 = find(index2_0.asset<0 | index2_0.ind<0);
data2_0(index0,:)=[];
index2_0.var(index0,:)=[];
index2_0.ind(index0,:)=[];
index2_0.asset(index0,:)=[];
index2_0.country(index0,:)=[];

clear index0

% Gross fixed capital formation

% Nominal gross fixed capital formation_2009 release
Ki1=zeros(length(asset1),ts,r,length(industry1));
Ki1_0=zeros(length(asset1),ts,r); % total gfcf of assets
for a = 1:length(asset1)
    for i=1:length(country)
        index0 = find(index1_0.var==3 & index1_0.asset==a & index1_0.country==i);
        if isempty(index0)
            Ki1_0(a,:,i)=NaN(1,ts);
        else
            Ki1_0(a,1:length(year1)-25,i)=data1_0(index0,26:end); %26: year1995
        end
        for j=1:length(industry1)
            index=find(index1.var==3 & index1.asset==a & index1.country==i & index1.ind==j);
            if isempty(index)
                Ki1(a,:,i,j)=NaN(1,ts);
            else
                Ki1(a,1:length(year1)-25,i,j)=data1(index,26:end); %26: year1995
            end
        end
    end
end
sum(sum(sum(sum(Ki1,'omitnan'))))

% Nominal gross fixed capital formation_2017 release
Ki2=zeros(length(asset2),ts,r,length(industry2));
Ki2_0=zeros(length(asset2),ts,r);
for a=1:length(asset2)
    for i=1:length(country)
        index0=find(index2_0.var==1 & index2_0.asset==a & index2_0.country==i);
        if isempty(index0)
            Ki2_0(a,:,i)=NaN(1,ts);
        else
            Ki2_0(a,:,i)=data2_0(index0,26:end);
        end
        for j=1:length(industry2)
            index=find(index2.var==1 & index2.asset==a & index2.country==i & index2.ind==j);
            if isempty(index)
                Ki2(a,:,i,j)=NaN(1,ts);
            else
                Ki2(a,:,i,j)=data2(index,26:end);
            end
        end
    end
end
sum(sum(sum(sum(Ki2,'omitnan'))))


% Real gross fixed capital formation, 1995 prices
Ki1_cons=zeros(length(asset1),ts,r,length(industry1));
Ki1_0_cons=zeros(length(asset1),ts,r);
for a=1:length(asset1)
    for i=1:length(country)
        index0=find(index1_0.var==6 & index1_0.asset==a & index1_0.country==i);
        if isempty(index0)
            Ki1_0_cons(a,:,i)=NaN(1,ts);
        else
            Ki1_0_cons(a,1:length(year1)-25,i)=data1_0(index0,26:end); %26: year1995
        end
        for j=1:length(industry1)
            index=find(index1.var==6 & index1.asset==a & index1.country==i & index1.ind==j);
            if isempty(index)
                Ki1_cons(a,:,i,j)=NaN(1,ts);
            else
                Ki1_cons(a,1:length(year1)-25,i,j)=data1(index,26:end); %26: year1995
            end
        end
    end
end
sum(sum(sum(sum(Ki1_cons,'omitnan'))))

% Real gross fixed capital formation, 2010 prices
Ki2_cons=zeros(length(asset2),ts,r,length(industry2));
Ki2_0_cons=zeros(length(asset2),ts,r);
for a=1:length(asset2)
    for i=1:length(country)
        index0=find(index2_0.var==3 & index2_0.asset==a & index2_0.country==i);
        if isempty(index0)
            Ki2_0_cons(a,:,i)=NaN(1,ts);
        else
            Ki2_0_cons(a,:,i)=data2_0(index0,26:end);
        end
        for j=1:length(industry2)
            index=find(index2.var==3 & index2.asset==a & index2.country==i & index2.ind==j);
            if isempty(index)
                Ki2_cons(a,:,i,j)=NaN(1,ts);
            else
                Ki2_cons(a,:,i,j)=data2(index,26:end);
            end
        end
    end
end
sum(sum(sum(sum(Ki2_cons,'omitnan'))))

% Capital formation price

%Gross fixed capital formation price index (1995=1.00)
Kp1=zeros(length(asset1),ts,r,length(industry1));
Kp1_0=zeros(length(asset1),ts,r);
for a=1:length(asset1)
    for i=1:length(country)
        index0=find(index1_0.var==5 & index1_0.asset==a & index1_0.country==i);
        if isempty(index0)
            Kp1_0(a,:,i)=NaN(1,ts);
        else
            Kp1_0(a,1:length(year1)-25,i)=data1_0(index0,26:end); %26: year1995
        end
        for j=1:length(industry1)
            index=find(index1.var==5 & index1.asset==a & index1.country==i & index1.ind==j);
            if isempty(index)
                Kp1(a,:,i,j)=NaN(1,ts);
            else
                Kp1(a,1:length(year1)-25,i,j)=data1(index,26:end); %26: year1995
            end
        end
    end
end
sum(sum(sum(sum(Kp1,'omitnan'))))

%Gross fixed capital formation price index (2010=100)
Kp2=zeros(length(asset2),ts,r,length(industry2));
Kp2_0=zeros(length(asset2),ts,r);
for a=1:length(asset2)
    for i=1:length(country)
        index0=find(index2_0.var==2 & index2_0.asset==a & index2_0.country==i);
        if isempty(index0)
            Kp2_0(a,:,i)=NaN(1,ts);
        else
            Kp2_0(a,:,i)=data2_0(index0,26:end);
        end
        for j=1:length(industry2)
            index=find(index2.var==2 & index2.asset==a & index2.country==i & index2.ind==j);
            if isempty(index)
                Kp2(a,:,i,j)=NaN(1,ts);
            else
                Kp2(a,:,i,j)=data2(index,26:end);
            end
        end
    end
end
sum(sum(sum(sum(Kp2,'omitnan'))))

% Capital stock

% Real fixed capital stock, 1995 prices
Ks1_cons=zeros(length(asset1),ts,r,length(industry1));
for a=1:length(asset1)
    for i=1:length(country)
        for j=1:length(industry1)
            index=find(index1.var==7 & index1.asset==a & index1.country==i & index1.ind==j);
            if isempty(index)
                Ks1_cons(a,:,i,j)=NaN(1,ts);
            else
                Ks1_cons(a,1:length(year1)-25,i,j)=data1(index,26:end); %26: year1995
            end
        end
    end
end
sum(sum(sum(sum(Ks1_cons,'omitnan'))))

% Calculating: Nominal capital stock, in millions of national currency, 2007 release

Ks1 = Ks1_cons.*Kp1; % current = constant * price index
sum(sum(sum(sum(Ks1,'omitnan'))))

% Real fixed capital stock, 2010 prices
Ks2_cons=zeros(length(asset2),ts,r,length(industry2));
for a=1:length(asset2)
    for i=1:length(country)
        for j=1:length(industry2)
            index=find(index2.var==5 & index2.asset==a & index2.country==i & index2.ind==j);
            if isempty(index)
                Ks2_cons(a,:,i,j)=NaN(1,ts);
            else
                Ks2_cons(a,:,i,j)=data2(index,26:end);
            end
        end
    end
end
sum(sum(sum(sum(Ks2_cons,'omitnan'))))


% Nominal capital stock, in millions of national currency (2017 release)
Ks2=zeros(length(asset2),ts,r,length(industry2));
for a=1:length(asset2)
    for i=1:length(country)
        for j=1:length(industry2)
            index=find(index2.var==4 & index2.asset==a & index2.country==i & index2.ind==j);
            if isempty(index)
                Ks2(a,:,i,j)=NaN(1,ts);
            else
                Ks2(a,:,i,j)=data2(index,26:end);
            end
        end
    end
end
sum(sum(sum(sum(Ks2,'omitnan'))))

% Capital depreciation_2007

% Consumption of fixed capital, 1995 prices
% Real fixed capital stock, 1995 prices
Kd1_cons=zeros(length(asset1),ts,r,length(industry1));
for a=1:length(asset1)
    for i=1:length(country)
        for j=1:length(industry1)
            index=find(index1.var==2 & index1.asset==a & index1.country==i & index1.ind==j);
            if isempty(index)
                Kd1_cons(a,:,i,j)=NaN(1,ts);
            else
                Kd1_cons(a,1:length(year1)-25,i,j)=data1(index,26:end); %26: year1995
            end
        end
    end
end
sum(sum(sum(sum(Kd1_cons,'omitnan'))))
% checked: this is calculated from Ks_cons1*depr1

% China

% cd('D:\WP\CapitalImpacts\Datasets\KLEMS');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
filename = 'China_CIP_3.0_(2015).xlsx';
sheet = 'all';
chinadata=importdata(filename,sheet);
data = chinadata.data.all;

deprCHN = chinadata.data.depreciationrate; % equipment, non-residential structure (invested by industrial enterprises), other (invested by agriculture, construction, and services)
deprCHN(:,1:2)=[];
deprCHN(:,1)=[];
deprCHN(:,2)=[];
deprCHN(:,3)=[];
deprCHN = horzcat(deprCHN,NaN(size(deprCHN,1),1)); % for consistent formating with those from the EU KLEMS dataset

data(1,:)=[]; %Col 1. asset (equipment, non-residential structures, other) 
%2. variable (investment in current million yuan, stock in million 1990 yuan, capital input (price) index (1990=1) 
%3. industry (37)

KiCHN = zeros(max(data(:,1)),ts,max(data(:,3))); 
for a = 1:max(data(:,1))
    for j = 1:max(data(:,3))
        index = find(data(:,2)==1 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KiCHN(a,1:16,j)=NaN(1,16);
        else
            KiCHN(a,1:16,j)=data(index,19:end); % years data available: 1980:2010
        end
    end
end
sum(sum(sum(sum(KiCHN,'omitnan'))))

KsCHN_cons=zeros(max(data(:,1)),ts,max(data(:,3))); 
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==2 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KsCHN_cons(a,1:16,j)=NaN(1,16);
        else
            KsCHN_cons(a,1:16,j)=data(index,19:end); % years data available: 1980:2010
        end
    end
end
sum(sum(sum(sum(KsCHN_cons,'omitnan'))))

KpCHN=zeros(max(data(:,3)),ts); 
KpCHN(:,1:16) = data(data(:,2)==3,19:end);
sum(sum(sum(sum(KpCHN,'omitnan'))))

% Korea
% cd('D:\WP\CapitalImpacts\Datasets\KLEMS');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
filename = 'KOR_WK_2015_K.xlsx';
sheet = 'all';
koreadata=importdata(filename,sheet);
data = koreadata.data.all; %Col 1. asset (7 categories) 
%2. variable (real investment in million Korean Won, real net stock in million Korean Won, investment deflator (2000=1) 
%3. industry (37)
deprKOR = koreadata.data.depreciationrate;

KiKOR_cons=zeros(max(data(:,1)),ts,max(data(:,3))); 
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==1 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KiKOR_cons(a,1:18,j)=NaN(1,18);
        else
            KiKOR_cons(a,1:18,j)=data(index,29:end); % years data available: 1970:2012
        end
    end
end
sum(sum(sum(sum(KiKOR_cons,'omitnan'))))

KsKOR_cons=zeros(max(data(:,1)),ts,max(data(:,3))); 
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==2 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KsKOR_cons(a,1:18,j)=NaN(1,18);
        else
            KsKOR_cons(a,1:18,j)=data(index,29:end); % years data available: 1970:2012
        end
    end
end
sum(sum(sum(sum(KsKOR_cons,'omitnan'))))

KpKOR=zeros(max(data(:,1)),ts,max(data(:,3)));  % current = constant * price index
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==3 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KpKOR(a,1:18,j)=NaN(1,18);
        else
            KpKOR(a,1:18,j)=data(index,29:end); % years data available: 1970:2012
        end
    end
end
sum(sum(sum(sum(KpKOR,'omitnan'))))

KiKOR = KiKOR_cons.*KpKOR;
sum(sum(sum(sum(KiKOR,'omitnan'))))

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
c_72x163= xlsread('CapitalProject_Concordances.xlsx','Korea_Industry_Concordance','D3:FJ74');
c_34x163= xlsread('CapitalProject_Concordances.xlsx','Industry_Concordance','D15:FK56');
c_34x163(c_34x163(:,1)==0,:)=[];
c_34x163(:,1)=[];
c_32x34=xlsread('CapitalProject_Concordances.xlsx','Industry_betweenKLEMSreleases','E5:AL36');
c_34x72 = c_34x163*c_72x163';
c_34x72(c_34x72>1)=1;
deprKOR_v1=zeros(size(KiKOR,3),size(KiKOR,1));
for j=1:72
    temp=find(c_34x72(:,j)==1);
    deprKOR_v1(j,:)=mean(deprKOR(temp,:),1);
end

deprKOR=deprKOR_v1;
deprKOR = horzcat(deprKOR,NaN(size(deprKOR,1),1));

% Canada
% cd('D:\WP\CapitalImpacts\Datasets\KLEMS');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
filename = 'can_capital_input_08I.xlsx';
sheet = 'all';
canadadata=importdata(filename,sheet);
data = canadadata.data.all;
%Col 1. asset (4 categories: Non-residential, Residential, ICT, and Non-ICT) 
%2. variable (nominal gfcf, real gfcf (1995 prices), gfcf price index (1995=1), real fixed capital stock (1995 prices) 
%3. industry (31)
deprCAN = canadadata.data.depreciationrate;
deprCAN = horzcat(deprCAN,NaN(size(deprCAN,1),1));

KiCAN=zeros(max(data(:,1)),ts,max(data(:,3))); % in thousand Canadian dollars
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==1 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KiCAN(a,1:10,j)=NaN(1,10);
        else
            KiCAN(a,1:10,j)=data(index,29:end); % years data available: 1970:2004
        end
    end
end
sum(sum(sum(sum(KiCAN,'omitnan'))))

KiCAN_cons=zeros(max(data(:,1)),ts,max(data(:,3))); % in thousand 1995 Canadian dollars
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==2 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KiCAN_cons(a,1:10,j)=NaN(1,10);
        else
            KiCAN_cons(a,1:10,j)=data(index,29:end); % years data available: 1970:2004
        end
    end
end
sum(sum(sum(sum(KiCAN_cons,'omitnan'))))

KpCAN=zeros(max(data(:,1)),ts,max(data(:,3)));  % current = constant * price index
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==3 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KpCAN(a,1:10,j)=NaN(1,10);
        else
            KpCAN(a,1:10,j)=data(index,29:end); % years data available: 1970:2004
        end
    end
end
sum(sum(sum(sum(KpCAN,'omitnan'))))

KsCAN_cons=zeros(max(data(:,1)),ts,max(data(:,3))); %in thousand 1995 Canadian dollars
for a=1:max(data(:,1))
    for j=1:max(data(:,3))
        index=find(data(:,2)==4 & data(:,1)==a & data(:,3)==j);
        if isempty(index)
            KsCAN_cons(a,1:10,j)=NaN(1,10);
        else
            KsCAN_cons(a,1:10,j)=data(index,29:end); % years data available: 1970:2004
        end
    end
end
sum(sum(sum(sum(KsCAN_cons,'omitnan'))))

clear sheet

% Data availability

% cd('D:\WP\CapitalImpacts');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
%Col 1: start year, 2: end year
da1 = data_availability(Ki1,ts,country);
da2 = data_availability(Ki2,ts,country);

% Col 1:4 1) start year; 2) data 1 or 2 for the start year; 3) the first end year; 
% 4) year to move on to another dataset (data2) if necessary 5) the second end year
t0 = zeros(r,5);
for i=1:r
    if da1(i,1)>0
        if da2(i,1)>0
            t0(i,1) = min(da1(i,1),da2(i,1));
        else
            t0(i,1) = da1(i,1);
        end
        if t0(i,1) == da2(i,1)
            t0(i,2) = 2;
            t0(i,3) = da2(i,2);
        else
            t0(i,2)=1;
            t0(i,3)=da1(i,2);
            if da2(i,2)>da1(i,2)
                t0(i,4)=min(max(26,da2(i,1)),da1(i,2)+1);
                t0(i,5)=da2(i,2);
            else
            end
        end
    else
        if da2(i,1)>0
            t0(i,1) = da2(i,1);
            t0(i,2) = 2;
            t0(i,3) = da2(i,2);
        end
    end
end

% turns out only UK needs to use data from two releases;
% industry-specific data in four countries are too sparse to use:  Belgium, Croatia (NOT available in KLEMS at all); Ireland, and Lithuania

c=find(sum(t0,2)==0);
da2_0 = data_availability(Ki2_0,ts,country);
t0(c,1)=da2_0(c,1);
t0(c,2)=3.*(t0(c,1)>0); % using total-industry data
t0(c,3)=da2_0(c,2);
c0=find(t0(:,2)==3);
clear c

% Missing data imputation

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
load('Yk');

Ki1a = Ki1;
Ki2a = Ki2;
Ki2_0a = Ki2_0;
Kp1a = Kp1; 
Kp2a = Kp2;
Kp2_0a = Kp2_0;

% later years are missing
for i=1:r
    if t0(i,2)==1 && t0(i,4)>0 && t0(i,5)<ts % using two releases
        for t=t0(i,5)+1:ts
            temp = max(sum(Ki2(end,t0(i,5),i,:),'omitnan'), sum(sum(Ki2(1:end-1,t0(i,5),i,:),'omitnan')));
            Ki2a(:,t,i,:)= squeeze(Ki2(:,t0(i,5),i,:)).*sum(Yk(t).Yk(:,i))./temp;
            Kp2a(:,t,i,:)= Kp2(:,t0(i,5),i,:);
        end
    end
    
    if t0(i,2)==1 && t0(i,4)==0 && t0(i,3)<ts % using only 2009 release
        for t=t0(i,3)+1:ts
            temp = max(sum(Ki1(end,t0(i,3),i,:),'omitnan'), sum(sum(Ki1(1:end-1,t0(i,3),i,:),'omitnan')));
            Ki1a(:,t,i,:)= squeeze(Ki1(:,t0(i,3),i,:)).*sum(Yk(t).Yk(:,i))./temp;
            Kp1a(:,t,i,:)= Kp1(:,t0(i,3),i,:);
        end
    end
        
    if t0(i,2)>=2 && t0(i,3)<ts % using only 2017 release
        for t=t0(i,3)+1:ts
            if any(c0==i)
                temp = max(sum(Ki2_0(end,t0(i,3),i),'omitnan'), sum(Ki2_0(1:end-1,t0(i,3),i),'omitnan'));
                Ki2_0a(:,t,i)= squeeze(Ki2_0(:,t0(i,3),i)).*sum(Yk(t).Yk(:,i))./temp;
                Kp2_0a(:,t,i)= Kp2_0(:,t0(i,3),i);
            else
                temp = max(sum(Ki2(end,t0(i,3),i,:),'omitnan'), sum(sum(Ki2(1:end-1,t0(i,3),i,:),'omitnan')));
                Ki2a(:,t,i,:)= squeeze(Ki2(:,t0(i,3),i,:)).*sum(Yk(t).Yk(:,i))./temp;
                Kp2a(:,t,i,:)= Kp2(:,t0(i,3),i,:);
            end
        end
    end
end

% earlier years are missing
for i=1:r
    if t0(i,2)==1 && t0(i,4)>0 && t0(i,1)>1 % using two releases
        for t=1:t0(i,1)-1
            temp = max(sum(Ki1(end,t0(i,1),i,:),'omitnan'), sum(sum(Ki1(1:end-1,t0(i,1),i,:),'omitnan')));
            Ki1a(:,t,i,:)= squeeze(Ki1(:,t0(i,1),i,:)).* sum(Yk(t).Yk(:,i))./temp;
            Kp1a(:,t,i,:)= Kp1(:,t0(i,1),i,:);
        end
    end
    
    if t0(i,2)==1 && t0(i,4)==0 && t0(i,1)>1  % using only 2009 release
        for t=1:t0(i,1)-1
            temp = max(sum(Ki1(end,t0(i,1),i,:),'omitnan'), sum(sum(Ki1(1:end-1,t0(i,1),i,:),'omitnan')));
            Ki1a(:,t,i,:)= squeeze(Ki1(:,t0(i,1),i,:)).* sum(Yk(t).Yk(:,i))./temp;
            Kp1a(:,t,i,:)= Kp1(:,t0(i,1),i,:);
        end
    end
    
    if t0(i,2)>=2 && t0(i,1)>1 % using only 2017 release
        for t=1:t0(i,1)-1
            if any(c0==i)
                temp = max(sum(Ki2_0(end,t0(i,1),i),'omitnan'), sum(sum(Ki2_0(1:end-1,t0(i,1),i),'omitnan')));
                Ki2_0a(:,t,i)= squeeze(Ki2_0(:,t0(i,1),i)).* sum(Yk(t).Yk(:,i))./temp;
                Kp2_0a(:,t_adj,i)= Kp2_0(:,t0(i,1),i);
            else
                temp = max(sum(Ki2(end,t0(i,1),i,:),'omitnan'), sum(sum(Ki2(1:end-1,t0(i,1),i,:),'omitnan')));
                Ki2a(:,t,i,:)= squeeze(Ki2(:,t0(i,1),i,:)).* sum(Yk(t).Yk(:,i))./temp;
                Kp2a(:,t,i,:)= Kp2(:,t0(i,1),i,:);
            end
        end
    end
end
clear temp

% Checking data availability again
% cd('D:\WP\CapitalImpacts');
cd('J:\MATLAB\Data\Capital_Endogenization\Codes');
da1a = data_availability(Ki1a,ts,country);
da2a = data_availability(Ki2a,ts,country);

% Col 1:4 1) start year; 2) data 1 or 2 for the start year; 3) the first end year; 
% 4) year to move on to another dataset (data2) if necessary 5) the second end year
t0a = zeros(r,5);
for i=1:r
    if da1a(i,1)>0
        if da2a(i,1)>0
            t0a(i,1) = min(da1a(i,1),da2a(i,1));
        else
            t0a(i,1) = da1a(i,1);
        end
        if t0a(i,1) == da2a(i,1)
            t0a(i,2) = 2;
            t0a(i,3) = da2a(i,2);
        else
            t0a(i,2)=1;
            t0a(i,3)=da1a(i,2);
            if da2a(i,2)>da1a(i,2)
                t0a(i,4)=min(max(26,da2a(i,1)),da1a(i,2)+1);
                t0a(i,5)=da2a(i,2);
            else
            end
        end
    else
        if da2a(i,1)>0
            t0a(i,1) = da2a(i,1);
            t0a(i,2) = 2;
            t0a(i,3) = da2a(i,2);
        end
    end
end

% adding industry-total data: Ireland, and Lithuania

da2_0a = data_availability(Ki2_0a,ts,country);
t0a(c0,1)=da2_0a(c0,1);
t0a(c0,2)=3; % using total-industry data
t0a(c0,3)=da2_0a(c0,2);

% Missing data imputation_China, Canada, and South Korea
% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
load('Yk');

KiCHNa = KiCHN;
KpCHNa = KpCHN;
t1=1;
t2=16;
% later years are missing
i=31;
for t=t2+1:ts
    temp = sum(sum(KiCHN(:,t2,:),'omitnan'));
    KiCHNa(:,t,:)= squeeze(KiCHN(:,t2,:)).*sum(Yk(t).Yk(:,i))./temp;
    KpCHNa(:,t,:)= KpCHN(:,t2,:);
end

KiKORa = KiKOR;
KpKORa = KpKOR;
t1=1;
t2=18;
% later years are missing
i=33;
for t=t2+1:ts
    temp = sum(sum(KiKOR(:,t2,:),'omitnan'));
    KiKORa(:,t,:)= squeeze(KiKOR(:,t2,:)).*sum(Yk(t).Yk(:,i))./temp;
    KpKORa(:,t,:)= KpKOR(:,t2,:);
end

KiCANa = KiCAN;
KpCANa = KpCAN;
t1=1;
t2=10;
% later years are missing
i=32;
for t=t2+1:ts
    temp = sum(sum(KiCAN(:,t2,:),'omitnan'));
    KiCANa(:,t,:)= squeeze(KiCAN(:,t2,:)).*sum(Yk(t).Yk(:,i))./temp;
    KpCANa(:,t,:)= KpCAN(:,t2,:);
end

% Other regions

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder');
filename = 'PWTIc.xlsx';
PWTIcdata=importdata(filename);
data1 = PWTIcdata.data.all;
data2 = PWTIcdata.data.allcountries;
data1(1,:)=[];
data2(:,1)=[];
data2=data2';
deprPWT = PWTIcdata.data.depreciationrate;
deprPWT = horzcat(deprPWT,NaN(size(deprPWT,1),1)); % for consistent formating with those from the EU KLEMS dataset

KiPWT = zeros(4,ts,r); % four assets: structure, machineary, transportation equipment, and other
for i=1:r
    for a=1:size(KiPWT,1)
        if r<45
            index = find(data1(:,1)==a & data1(:,2)==i);
            KiPWT(a,:,i)=data1(index,3:end);
        else
            KiPWT(a,:,i)=data2(a,:);
        end
    end
end
for n=1:ts
    for i=1:r
        KiPWT(:,n,i) = KiPWT(:,n,i).*sum(Yk(n).Yk(:,i))./sum(KiPWT(:,n,i));
    end
end

% EU KLEMS countries with complete data (1995-2015)

da = zeros(length(country),2);
for i=1:length(country)
    temp=[];
    for t=1:ts
        temp1 = squeeze(Ki2(:,t,i,:))';
        if length(find(sum(temp1(:,1:end-1),1,'omitnan')>0))==10 %industry-level data gfcf is available or asset, industry-specfiic data is available
            temp=horzcat(temp,t);
        end
    end
    
    if any(temp)
        da(i,1)= temp(1);
        da(i,2)= temp(end);
    end
end
c = find(sum(da,2)==22);

gfcf_c_proxy = zeros(4,ts,size(Ki2,4));
for i=1:length(c)
    for n=1:ts
        temp = squeeze(Ki2(1:end-1,n,c(i),:))';
        temp(temp<0)=0;
        temp(isnan(temp))=0;
        temp1=zeros(size(Ki2,4),4);
        temp1(:,1)=sum(temp(:,6:7),2);
        temp1(:,2)=sum(temp(:,[1 2 5]),2);
        temp1(:,3)=temp(:,4);
        temp1(:,4)=sum(temp,2)-sum(temp(:,[1 2 4:7]),2);
        
        temp1 = temp1./sum(temp1,'omitnan')./length(c);
        gfcf_c_proxy(:,n,:)=temp1'+squeeze(gfcf_c_proxy(:,n,:));
    end
end

c = setdiff(1:r,horzcat(find(t0(:,1)>0)',31:33));
Ki3 = zeros(4,ts,r,size(Ki2,4)); 
for i=1:length(c)
    for n=1:ts
        for a=1:4
            Ki3(a,n,c(i),:)=KiPWT(a,n,c(i)).*gfcf_c_proxy(a,n,:);
        end
    end
end

save('KLEMS.mat','depr1','depr2','Ki1a','Ki2a','Ki1','Ki2','Kp1a','Kp2a','country','t0a','c0','Ki2_0a','Kp2_0','KiCHNa','KiKORa','KiCANa','KpCHNa','KpKORa','KpCANa','deprCHN','deprCAN','deprKOR','Ki3','deprPWT');
