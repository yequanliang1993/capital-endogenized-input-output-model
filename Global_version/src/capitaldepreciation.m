%% This is code for calculating capital depreciation

% Created by: Quanliang Ye
% Created date: 19/10/2019
% Email add.: yequanliang1993@gmail.com
% Main Contributor: dr. Ranran Wang, CML, Leiden University


%%
clear
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
load('KLEMS.mat')
load('Yk.mat');
load('cfc.mat');
load('x.mat');

c_34x163= xlsread('CapitalProject_Concordances.xlsx','Industry_Concordance','D15:FK56');
c_34x163(c_34x163(:,1)==0,:)=[];
c_34x163(:,1)=[];
c_163x200=xlsread('CapitalProject_Concordances.xlsx','EXIOBASE_industry2product','C3:GT165');
c_32x34=xlsread('CapitalProject_Concordances.xlsx','Industry_betweenKLEMSreleases','E5:AL36');
c_32x200 = c_32x34*c_34x163*c_163x200;
c_32x200(c_32x200>1)=1;
c_34x200 = c_34x163*c_163x200;
c_34x200(c_34x200>1)=1;
c_8x200=xlsread('CapitalProject_Concordances.xlsx','Asset_Concordance','F21:GW28');
c_10x200=xlsread('CapitalProject_Concordances.xlsx','Asset_Concordance','F6:GW15');
c_4x200_pwt=xlsread('CapitalProject_Concordances.xlsx','PWT_Asset_Concordance','C5:GT8');

r=49;
s=200;
ts=21;

t0=t0a;
%%   D_S 1
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
cd_ss = zeros(r,ts,ts,s,s);
c=find(t0(:,2)==1);   % from which dataset, data1 or data2
for i=1:length(c)
    t1 = t0(c(i),1);
    if t0(c(i),4)>0
        t2 = t0(c(i),4)-1;   % the first end year
    else
        t2 = t0(c(i),3);
    end

    for n=1:ts % each exiobase year
        for t=t1:min(n,t2) % each previous investment year
            Ki = squeeze(Ki1a(:,t,c(i),:))';
            
            % check whether there are detailed sector disaggregated data
            if sum(sum(Ki(:,1:end-1),'omitnan'))/sum(Ki(:,end),'omitnan')>=0.95 || length(find(sum(Ki(:,1:end-1),1,'omitnan')>0))==size(Ki1a,1)-1
                Ki(:,end)=[];
            else
                Ki(:,1:end-1)=[];
            end
            Ki(Ki<0)=0;
            cd_ss(c(i),n,t,:,:) = depr_calc(c(i),Ki,depr1,c_8x200,c_32x200,Yk,cfc,t,n);
        end
    end
    country(c(i))
end

sum(sum(sum(sum(sum(cd_ss,'omitnan')))))

%% D_S 2
c = find(t0(:,2)>=2);
for i=1:length(c)
    t1 = t0(c(i),1);
    t2 = t0(c(i),3);
    
    for n=1:ts % each exiobase year
        for t=t1:min(n,t2) % each previous investment year
            if any(c0==c(i))
                Ki = squeeze(Ki2_0a(1:end-1,t,c(i),:))';
            else
                Ki = squeeze(Ki2a(:,t,c(i),:))';
                if sum(sum(Ki(:,1:end-1),'omitnan'))/sum(Ki(:,end),'omitnan')>=0.95 || length(find(sum(Ki(:,1:end-1),1,'omitnan')>0))==size(Ki2a,1)-1
                    Ki(:,end)=[];
                else
                    Ki(:,1:end-1)=[];
                end
            end
            Ki(Ki<0)=0;
            cd_ss(c(i),n,t,:,:)=depr_calc(c(i),Ki,depr2,c_10x200,c_34x200,Yk,cfc,t,n);
        end
    end
    country(c(i))
end
    
sum(sum(sum(sum(sum(cd_ss,'omitnan')))))
        

%% D_S 3 for GBR
c=find(t0(:,4)>0);
for i=1:length(c)
    t1 = t0(c(i),4);
    t2 = t0(c(i),5);    
    for n=1:ts % each exiobase year
        for t=t1:min(n,t2) % each previous investment year
            Ki = squeeze(Ki2a(:,t,c(i),:))';
            if sum(sum(Ki(:,1:end-1),'omitnan'))/sum(Ki(:,end),'omitnan')>=0.95 || length(find(sum(Ki(:,1:end-1),1,'omitnan')>0))==size(Ki2a,1)-1
                Ki(:,end)=[];
            else
                Ki(:,1:end-1)=[];
            end
            Ki(Ki<0)=0;
            cd_ss(c(i),n,t,:,:)=depr_calc(c(i),Ki,depr2,c_10x200,c_34x200,Yk,cfc,t,n);
        end
    end
    country(c(i))
end

sum(sum(sum(sum(sum(cd_ss,'omitnan')))))

%% Plots
c = 31;
ts = 21;
temp=zeros(ts,ts);
for n=1:ts
    for t=1:ts
        temp(t,n) = sum(sum(cd_ss(c,n,t,:,:),'omitnan'));
    end
end

bar(temp','stacked');

%% Plots
test=zeros(ts,r);
for i=1:r %[find(sum(t0a,2)>0)' 31:33]
    for n=1:ts
        test(n,i)=sum(Yk(n).Yk(:,i));
    end
end

test1=zeros(ts,r);
for i=1:r %[find(sum(t0a,2)>0)' 31:33]
    for n=1:ts
        test1(n,i)=sum(sum(sum(cd_ss(i,:,n,:,:),'omitnan')));
    end
end

subplot(1,2,1)
plot(sum(test1,'omitnan')./sum(test,'omitnan'))

subplot(1,2,2)
plot(sum(test1,2,'omitnan')./sum(test,2,'omitnan'))



%% China
% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c_3x200 = xlsread('CapitalProject_Concordances.xlsx','China_Asset_Concordance','B3:GS5');
c_37x200 = xlsread('CapitalProject_Concordances.xlsx','China_Industry_Concordance','G3:GX39');

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c=31;
t1 = 1;
t2 = 21;
for n=1:ts % each exiobase year
    for t=t1:min(n,t2) % each previous investment year
        Ki = squeeze(KiCHNa(:,t,:))';
        Ki(Ki<0)=0;
        cd_ss(c,n,t,:,:)=depr_calc(c,Ki,deprCHN,c_3x200,c_37x200,Yk,cfc,t,n);
    end
end
country(c)

sum(sum(sum(sum(sum(cd_ss,'omitnan')))))

%% Canada

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c_4x200 = xlsread('CapitalProject_Concordances.xlsx','Canada_Asset_Concordance','C5:GT8');
c_31x163 = xlsread('CapitalProject_Concordances.xlsx','Canada_Industry_Concordance','C4:FI34');
c_31x200 = c_31x163*c_163x200;
c_31x200(c_31x200>1)=1;

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c=32;
t1 = 1;
t2 = 21;
for n=1:ts % each exiobase year
    for t=t1:min(n,t2) % each previous investment year
        Ki = squeeze(KiCANa(:,t,:))';
        Ki(Ki<0)=0;
        cd_ss(c,n,t,:,:)=depr_calc(c,Ki,deprCAN,c_4x200,c_31x200,Yk,cfc,t,n);
    end
end
country(c)

sum(sum(sum(sum(sum(cd_ss,'omitnan')))))

%% South Korea

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c_7x200 = xlsread('CapitalProject_Concordances.xlsx','Korea_Asset_Concordance','E5:GV11');
c_72x163 = xlsread('CapitalProject_Concordances.xlsx','Korea_Industry_Concordance','D3:FJ74');
c_72x200 = c_72x163*c_163x200;
c_72x200(c_72x200>1)=1;

% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
c=33;
t1 = 1;
t2 = 21;
for n=1:ts % each exiobase year
    for t=t1:min(n,t2) % each previous investment year
        Ki = squeeze(KiKORa(:,t,:))';
        Ki(Ki<0)=0;
        cd_ss(c,n,t,:,:)=depr_calc(c,Ki,deprKOR,c_7x200,c_72x200,Yk,cfc,t,n);
    end
end
country(c)

sum(sum(sum(sum(sum(cd_ss,'omitnan')))))


%% Other regions
% cd('D:\WP\CapitalImpacts\Datasets');
cd('J:\MATLAB\Data\Capital_Endogenization\New folder')
t1 = 1;
t2 = 21;
c=setdiff(1:r,[find(sum(t0,2)>0)' 31:33]);
for i=1:length(c)
    for n=1:ts % each exiobase year
        for t=t1:min(n,t2) % each previous investment year
            Ki = squeeze(Ki3(:,t,c(i),:))';
            Ki(Ki<0)=0;
            cd_ss(c(i),n,t,:,:)=depr_calc(c(i),Ki,deprPWT,c_4x200_pwt,c_34x200,Yk,cfc,t,n);
        end
    end
    country(c(i))
end
    
sum(sum(sum(sum(sum(cd_ss,'omitnan')))))
    
%%
save('cd_ss.mat','cd_ss','-v7.3');

    
