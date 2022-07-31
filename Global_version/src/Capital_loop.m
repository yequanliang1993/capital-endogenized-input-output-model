%% Capital impacts computation
clear

r=49;
s=200;
ts=21;

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder');
load('cd_ss.mat');
load('Yc.mat');
load('Yk.mat');
load Yt
load('x.mat');
load('exioS.mat');
load('F_hh.mat');

% impactlabs={'Energy use','Blue water consumption','Land use','Metal ores extraction','Non-metallic minerals extraction','GHG emissions','CO2 emissions','Material extractions'};
neeslabs={'Food','Clothing','Shelter','Mobility','Construction','Manufactured products','Services'};

index = 6; % greenhouse gas emission
Fk_n_t=zeros(10,ts,ts,r,r*s); %impacts embedded in depreciated capital goods. 1st dimension: 10 loops
fk_n_t=zeros(10,ts,ts,r,r*s);
for lp = 1:10 % 10 capital loops
    tic
    for t = 1:ts % production year of the depreciated capital goods
        fname = sprintf('L%d.mat',t);
        load(fname);

        for n = t:ts % year of final consumption
            Fk = zeros(r,s*r); % impacts embedded in capital goods depreciated, allocated to year n
            
            if lp == 1
                f0 = squeeze(S(index,t,:))';
            else
                f0 = squeeze(Fk_n_t(lp-1,n,t,:,:))./x(t).x';
                f0(isnan(f0)) = 0;
                f0(isinf(f0)) = 0;
                
                fk_n_t(lp,n,t,:,:) = f0;
            end
            
            for i = 1:r
                gfcf_p = Yk(t).Yk(:,i); %i: investing country
                if sum(sum(cd_ss(i,n,t,:,:),'omitnan'))>0 % capital goods produced in year t and depreciated in year n, invested by country i
                    temp = squeeze(cd_ss(i,n,t,:,:)); % rows (dimension 4): capital goods producing sectors; columns: capital goods consuming sectors
                    temp(isnan(temp)) = 0;
                    temp0 = reshape(gfcf_p,[s,r]); % capital goods purchased in year t by country i; columns: country where the capital goods were purchased from
                    temp0 = temp0./sum(temp0,2);
                    temp0(isnan(temp0)) = 0;
                    temp0(isinf(temp0)) = 0;
                    temp1 = zeros(r*s,s);
                    for m = 1:r % 'finished' capital goods producer
                        for j = 1:s % sectors that produced the 'finished' capital goods
                            temp1((m-1)*s+j,:) = temp0(j,m).*temp(j,:); % for the 'same' capital goods produced in the same year, assuming the cross-country distributions of depreciation and investment are idential
                        end
                    end
                    clear temp temp0
                    
                    if lp == 1
                        temp2 = f0'.*L*temp1;
                        for m = 1:r % m: countries producing for the capital goods; i: countries using the capital goods; columns: capital using sectors in country i
                            Fk(m,(i-1)*s+1:i*s)=sum(temp2((m-1)*s+1:m*s,:));
                        end
                    else
                        temp2 = f0*L*temp1;
                        Fk(:,(i-1)*s+(1:s)) = temp2;
                    end
                    
                    clear temp1 temp2
                end
            end
            
            Fk_n_t(lp,n,t,:,:) = Fk;
            clear Fk f0
        end
        clear L
    end
    clear S
    
    lp
    toc
end
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
save('EF_k_loop','Fk_n_t','fk_n_t','-v7.3')

%% capital-related GHG emissions
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption');
load('EF_k_loop','Fk_n_t')

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder');
load('Yc.mat');
load('x.mat');


r=49;
s=200;
ts=21;

index_p = repmat([1:200]',r,1);

EFk_l_d_c_d_p=zeros(10,ts,ts,r,s,r); 
% 1 dimension: 10 capital loops
% 2 dimension: final demand year
% 3 dimension: capital investment year
% 4 dimension: where final demand consumed
% 5 dimension: 7 human needs
% 6 dimension: where capital environmental pressure occured
% 7 dimension: internal or external 
fk_n_t = zeros(10,ts,ts,r,r*s);
for lp = 1:10
    tic
    for n = 1:ts % year of final consumption
        fname = sprintf('L%d.mat',n);
        load(fname);
        
        y = zeros(r*s,r*s); % final consumption (1:7 human needs categories)
        for i = 1:r
            for j = 1:s
                y(index_p == j,(i-1)*s+j) = Yc(n).Yc(index_p == j,i);
            end
        end
        
        for t = 1:n
            Fk = squeeze(Fk_n_t(lp,n,t,:,:));
            fk = Fk./x(n).x'; % environmental multiplier for capital goods depreciated in year n (rows: countries where ...
            % environemntal impacts occurred directly; columns: country-sector -where capital goods are consumed)
            % Note: depreciated capitals are used to produce both noncapital goods consumed in year n and capital goods to be consumed...
            % in year n, n+1 ...
            
            fk(isnan(fk)) = 0;
            fk(isinf(fk)) = 0;
            fk_n_t(lp,n,t,:,:) = fk;
            
            clear Fk
            
            for m = 1:r
                temp_fk = zeros(1,r*s);
                temp_fk(:,(m-1)*s+1:m*s) = sum(fk(:,(m-1)*s+1:m*s));
                temp = temp_fk*L*y;

                EFk_l_d_c_d_p(lp,n,t,:,:,m) = reshape(temp',s,r)'; %dimension3: country of final consumption; %dimension 5: where capital is consumed; 6: impacts at capital consuming country

                clear temp temp_fk
            end
            
            clear fk
        end
        
        clear y L
    end
    lp
    toc
end
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
save('EF_k_c','EFk_l_d_c_d_p','fk_n_t','-v7.3')

%% assigning to final demand not only final consumption
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption');
load('EF_k_loop','Fk_n_t')

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder');
load('Yc.mat');
load('x.mat');
load Yk


r=49;
s=200;
ts=21;

index_p = repmat([1:200]',r,1);

EFk_d_p_k_fd=zeros(ts,ts,r,s,r,2); 
% 1 dimension: final demand year
% 2 dimension: capital investment year
% 3 dimension: where final demand consumed
% 4 dimension: 200 product
% 5 dimension: where capital environmental pressure occured
% 6 dimension: 1 for final conusmption and 2 for GFCF

for n = 1:ts % year of final consumption
    fname = sprintf('L%d.mat',n);
    load(fname);
    
    yc = zeros(r*s,r*s); % final consumption 
    yk = zeros(r*s,r*s); % GFCF
    for i = 1:r
        for j = 1:s
            yc(index_p == j,(i-1)*s+j) = Yc(n).Yc(index_p == j,i);
            yk(index_p == j,(i-1)*s+j) = Yk(n).Yk(index_p == j,i);
        end
    end
    
    for t = 1:n
        Fk = squeeze(sum(Fk_n_t(:,n,t,:,:),1));
        fk = Fk./x(n).x'; % environmental multiplier for capital goods depreciated in year n (rows: countries where ...
        % environemntal impacts occurred directly; columns: country-sector -where capital goods are consumed)
        % Note: depreciated capitals are used to produce both noncapital goods consumed in year n and capital goods to be consumed...
        % in year n, n+1 ...
        
        fk(isnan(fk)) = 0;
        fk(isinf(fk)) = 0;
        
        clear Fk
        
        for m = 1:r
            temp_fk = zeros(1,r*s);
            temp_fk(:,(m-1)*s+1:m*s) = sum(fk(:,(m-1)*s+1:m*s));
            temp_yc = temp_fk*L*yc;
            temp_yk = temp_fk*L*yk;
            
            EFk_d_p_k_fd(n,t,:,:,m,1) = reshape(temp_yc',s,r)'; %dimension3: country of final consumption; %dimension 5: where capital is consumed; 6: impacts at capital consuming country
            EFk_d_p_k_fd(n,t,:,:,m,2) = reshape(temp_yk',s,r)';
            clear temp_yc temp_fk temp_yk
        end
        
        clear fk
    end
    
    clear yc L yk
    
    n
end

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
save('EFk_fc_gfcf','EFk_d_p_k_fd','-v7.3')


%% EF of final consumption, gfcf and final demand
cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder');
load exioS
load Yk
load Yt
load Yc

EF_fc = zeros(8,ts,r,r,s);
EF_gfcf = zeros(8,ts,r,r,s);
EF_fd = zeros(8,ts,r,r,s);
for n = 1:ts
    fname = sprintf('L%d',n);
    load(fname);
    
    for index = 1:8
        f0 = squeeze(S(index,n,:))';
        f0_l = diag(f0)*L;
        
        temp_ef_fc = f0_l*Yc(n).Yc;
        temp_ef_gfcf = f0_l*Yk(n).Yk;
        temp_ef_fd = f0_l*Yt(n).Yt;
        
        for i = 1:r % regions where ef occured
            EF_fc(index,n,:,i,:) = temp_ef_fc((i-1)*s+(1:s),:)';
            EF_gfcf(index,n,:,i,:) = temp_ef_gfcf((i-1)*s+(1:s),:)';
            EF_fd(index,n,:,i,:) = temp_ef_fd((i-1)*s+(1:s),:)';
        end
        
        clear f0 f0_l temp_ef_fc temp_ef_gfcf temp_ef_fd
    end
    
    clear L
    n
end

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
save('EF_indicator','EF_fc','EF_fd','EF_gfcf','-v7.3')


% by final demand sectors
index_p = repmat([1:200]',r,1);

EF_fc_a = zeros(8,ts,r,s,r);
EF_gfcf_a = zeros(8,ts,r,s,r);
EF_fd_a = zeros(8,ts,r,s,r);
for n = 1:5%ts
    fname = sprintf('L%d',n);
    load(fname);
    
    for index = 6
        f0 = squeeze(S(index,n,:))';
        f0_l = diag(f0)*L;
        
        yc = zeros(r*s,r*s); % final consumption
        yk = zeros(r*s,r*s); % GFCF
        yd = zeros(r*s,r*s); % final demand
        for i = 1:r
            for j = 1:s
                yc(index_p == j,(i-1)*s+j) = Yc(n).Yc(index_p == j,i);
                yk(index_p == j,(i-1)*s+j) = Yk(n).Yk(index_p == j,i);
                yd(index_p == j,(i-1)*s+j) = Yt(n).Yt(index_p == j,i);
            end
        end
        
        temp_ef_fc = f0_l*yc;
        temp_ef_gfcf = f0_l*yk;
        temp_ef_fd = f0_l*yd;
        
        for i = 1:r % regions where ef occured
            EF_fc_a(index,n,:,:,i) = reshape(sum(temp_ef_fc((i-1)*s+(1:s),:))',s,r)';
            EF_gfcf_a(index,n,:,:,i) = reshape(sum(temp_ef_gfcf((i-1)*s+(1:s),:))',s,r)';
            EF_fd_a(index,n,:,:,i) = reshape(sum(temp_ef_fd((i-1)*s+(1:s),:))',s,r)';
        end
        
        clear f0 f0_l temp_ef_fc temp_ef_gfcf temp_ef_fd yc yk yd
    end
    
    clear L
    n
end
        
%% figures
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
load('EF_k_c','EFk_l_d_c_d_p') 
load EF_indicator
load('EF_k_loop','Fk_n_t');

data = importdata('EFk_fc_from_Carl.xlsx');
EFk_fc_carl = data.data.Sheet1(:,2);
p_name = data.textdata.Sheet2;
clear data

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder');
load('EF1a','EFk_d_c_d_a')
load('Yc.mat');
load('Yk.mat');
load Yt
load x
load('F_hh.mat');
load('index_human_needs.mat');
neeslabs={'Clothing','Construction','Food','Manufacturing','Mobility','Services','Shelter'};

r = 49;
s = 200;
ts = 21;
c = 31;

index_hn = index_hn(1:s);
hn_seq = [3 1 7 5 2 4 6];

x_c = zeros(ts,s);
yc_c = zeros(ts,s);
yk_c = zeros(ts,s);
fd_c = zeros(ts,s);
for n = 1:ts
    x_c(n,:) = x(n).x((c-1)*s+(1:s));
    yc_c(n,:) = sum(reshape(Yc(n).Yc(:,c),s,r),2);
    yk_c(n,:) = sum(reshape(Yk(n).Yk(:,c),s,r),2);
    fd_c(n,:) = sum(reshape(Yt(n).Yt(:,c),s,r),2);
end
clear x Yk Yc Yt


fs = 12;
lw = 1;
index = 6;
c = 31;
fgn = [1:30 32:49];
color =[0 1 0; 0 0 1; 1 1 0; 1 0 1; 0.5 0.5 0.5; 0.8 0.6 0.9; 1 0.7 0.7];
figure
temp_efkc = zeros(4*7,10); % product categories that assiged with the largest capital-related GHG emissions
temp_efgfcf = zeros(4*7,1);
temp_efk = squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,:),1),2),3),6));
hn_p_n = {};
for i = 1:7
    hn = hn_seq(i);
    index_hn_p = find(index_hn == hn);
    
    
    [a b] = sort(temp_efk(index_hn_p),'descend');
    
    for j = 1:4
        if j ~= 4
            p = index_hn_p(b(j));
            temp_efkc((i-1)*4+j,:) = squeeze(sum(sum(sum(Fk_n_t(:,:,:,:,(c-1)*s+p),2),3),4));
            temp_efgfcf((i-1)*4+j) = squeeze(sum(sum(EF_gfcf(index,:,c,:,p),2),4));            
            hn_p_n((i-1)*4+j,1) = p_name(p);
        else
            p = index_hn_p(b(j:end));
            temp_efkc((i-1)*4+j,:) = squeeze(sum(sum(sum(sum(Fk_n_t(:,:,:,:,(c-1)*s+p),2),3),4),5));
            temp_efgfcf((i-1)*4+j) = squeeze(sum(sum(sum(EF_gfcf(index,:,c,:,p),2),4),5));
            hn_p_n((i-1)*4+j,1) = {'Others'};
        end
    end
    
    clear index_hn_p a b hn
end

% bar(horzcat(temp_efgfcf,sum(temp_efkc,2))/10^12);
bar([1:28]-0.1,temp_efgfcf/10^12,'BarWidth',0.2);
hold on
bar([1:28]+0.1,temp_efkc'/10^12,'stacked','BarWidth',0.2);

ax = gca;
ax.YLim = [0 8];
ax.XLim = [0.25 4*7+0.75];
ax.XTick = 1:4*7;
ax.XTickLabel = hn_p_n;
ax.XTickLabelRotation = 90;
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
for i = 1:7
    hn = hn_seq(i);
    hold on
    if hn == 4
        text((i-1)*4+0.5,ax.YLim(2)-1,neeslabs(hn));
    elseif hn == 2
        text((i-1)*4+0.8,ax.YLim(2)-1,neeslabs(hn));
    else 
        text((i-1)*4+1.5,ax.YLim(2)-1,neeslabs(hn));
    end
    
    if i < 7
        plot([(i-1)*4+4.5 (i-1)*4+4.5],[ax.YLim(1) ax.YLim(2)],'--','Color', [0.7 0.7 0.7]);
    end
end
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth=lw;
ax.FontSize=fs;
legend('EF^G^F^C^F','F^K_l_p_=_1','F^K_l_p_>_1','','','','','Location','west','Color','none','EdgeColor','none')

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4.5];
print('EFk_EFgfcf','-dpng','-r300')


%% comparing with others
fs = 12;
lw = 1;
index = 6;
c = 31;
fgn = [1:30 32:49];
figure
% subplot(1,2,1)
temp_ghg = vertcat(squeeze(sum(sum(EF_fc(index,:,c,:,:),4),5))+squeeze(F_hh(index,:,c)),...
    squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,:),1),3),5),6)),...
    squeeze(sum(sum(EF_gfcf(index,:,c,:,:),4),5)));
h = area(temp_ghg'/10^12);
h(1).FaceColor = [0.8 0.8 0.8];
h(2).FaceColor = [0.5 0.2 0.6];
h(3).FaceColor = [1 1 1];
hold on
plot(squeeze(sum(sum(sum(EFk_d_c_d_a(index,:,c,:,:,:),4),5),6))/10^12+squeeze(F_hh(index,:,c))/10^12+...
    squeeze(sum(sum(EF_fc(index,:,c,:,:),4),5))/10^12,'Color',[0.47,0.67,0.19],'LineWidth',lw+0.5);
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[1 1 1],'LineWidth',lw+0.5)
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[1 1 1],'LineWidth',lw+0.5)
hold on
plot(EFk_fc_carl,'-','Color','y','LineWidth',lw+0.5);

ax = gca;
ax.YLim = [0 15];
ax.YTick = 0:5:15;
ax.YTickLabel = {'0','5','10','15'};
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[0.7 0.7 0.7],'LineWidth',lw+0.5)
ax.XLim = [1 ts];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth = lw;
ax.FontSize = fs;
legend('EF^C','EF^K^C','EF^G^F^C^F','EF^C^K^C','','','EF^C^K^C',...
    'Orientation','horizontal','NumColumns',3,'Location','northwest','Color','none','EdgeColor','none');
clear temp_ghg


cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 3.5];
print('EFk_fc_trend','-dpng','-r300')
% 


%legend
% subplot(1,2,2)
figure
plot(squeeze(sum(sum(sum(EFk_d_c_d_a(index,:,c,:,:,:),4),5),6))/10^12+squeeze(F_hh(index,:,c))/10^12+...
    squeeze(sum(sum(EF_fc(index,:,c,:,:),4),5))/10^12,'Color',[0.47,0.67,0.19],'LineWidth',lw+0.5);
hold on
plot(EFk_fc_carl,'-','Color','y','LineWidth',lw+0.5);


ax = gca;
ax.YLim = [0 15];
ax.YTick = 0:5:15;
ax.YTickLabel = {'0','5','10','15'};
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[0.7 0.7 0.7],'LineWidth',lw+0.5)
ax.XLim = [1 ts];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth = lw;
ax.FontSize = fs;
legend('EF^C^K^C in Ye et al. (2021)','EF^C^K^C in Sodersten et al. (2018a)',...
    'Location','northwest','Color','none','EdgeColor','none');

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 3.5];
print('EFk_fc_trend_legend','-dpng','-r300')

%% anumal profile by human needs
fs = 12;
lw = 1;
index = 6;
c = 31;
fgn = [1:30 32:49];
% color = [0 0 1; 0.5 0.5 0.5; 0 1 0; 0.8 0.6 0.9; 1 0 1; 1 0.7 0.7; 1 1 0];
color =[0 1 0; 0 0 1; 1 1 0; 1 0 1; 0.5 0.5 0.5; 0.8 0.6 0.9; 1 0.7 0.7];
figure
for i = 1:7
    hn = hn_seq(i);
    hn_p = find(index_hn == hn);
    temp_efk_do(:,i) = squeeze(sum(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,hn_p,c,:),1),3),5),6),7));
    temp_efk_fgn(:,i) = squeeze(sum(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,hn_p,fgn,:),1),3),5),6),7));
    
    clear hn hn_p
end    
bar(temp_efk_do/10^12,'stacked');
hold on
bar(-temp_efk_fgn/10^12,'stacked');
ax = gca;
ax.YLim = [-0.2 1.7];
ax.YTick = -0:0.5:1.5;
ax.YTickLabel = {'0','0.5','1.0','1.5'};
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[0.7 0.7 0.7],'LineWidth',lw+0.5)
ax.XLim = [0.25 ts+0.75];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth=lw;
ax.FontSize=fs;
ax.ColorOrder = color;
legend(neeslabs(hn_seq),'Location','northwest','Color','none','EdgeColor','none');
clear temp_efk_do temp_efk_fgn

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 3.5];
print('Anual_profile_hn','-dpng','-r300')


%% annual profile
figure
colormap = jet(ts);
temp_c = squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,c,:),1),5),6),7));
temp_fgn = squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,fgn,:),1),5),6),7));
temp = temp_c+temp_fgn;
for n = 1:ts
    yr{n,1} = num2str(n+1994);
    age = zeros(21,1);
    age(1:n) = [n-0.5:-1:1-0.5];
    age_av(n) = sum(temp(n,:)./sum(temp(n,:)).*age');
    
    clear age
end
bar(temp_c/10^12,'stacked');
hold on
bar(-temp_fgn/10^12,'stacked');
ax = gca;
ax.YLim = [-0.2 1.7];
ax.YTick = 0:0.5:1.5;
ax.YTickLabel = {'0','0.5','1.0','1.5'};
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[0.7 0.7 0.7],'LineWidth',lw+0.5)
ax.XLim = [0.25 ts+0.75];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth=lw;
ax.FontSize=fs;
% for i = 1:ts
%     if i == 2
%         text(i-0.5,sum(temp_c(i,:))/10^12+0.1,'1.0','FontSize',fs)
%     elseif i == 11
%         text(i-0.5,sum(temp_c(i,:))/10^12+0.1,'4.0','FontSize',fs)
%     else
%         text(i-0.5,sum(temp_c(i,:))/10^12+0.1,num2str(round(age_av(i),1)),'FontSize',fs)
%     end
% end

for i = 1:ts
    if i <= 16
        if i == 2
            text(i-0.45,sum(temp_c(i,:))/10^12+0.1,'1.0','FontSize',fs-2,'Color',[0.7 0.7 0.7]);
        elseif i == 11
            text(i-0.45,sum(temp_c(i,:))/10^12+0.1,'4.0','FontSize',fs-2,'Color',[0.7 0.7 0.7]);
        else
            text(i-0.45,sum(temp_c(i,:))/10^12+0.1,num2str(round(age_av(i),1)),'FontSize',fs-2,'Color',[0.7 0.7 0.7]);
        end
    else
        text(i-0.45,sum(temp_c(i,:))/10^12+0.1,num2str(round(age_av(i),1)),'FontSize',fs);
    end
end



ax.ColorOrder = colormap;
% colorbar('Ticks',[0:0.25:1],'TickLabels',{'1995','2000','2005','2010','2015'},'Direction','reverse')
% legend(yr,'Location','eastoutside','Color','none','EdgeColor','none');
clear temp_c temp_fgn

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 7.5 3.5];
print('Annual profile','-dpng','-r300')



% legend
figure
colormap = parula(ts);
temp_c = squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,c,:),1),5),6),7));
temp_fgn = squeeze(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,fgn,:),1),5),6),7));
temp = temp_c+temp_fgn;
for n = 1:ts
    age = zeros(21,1);
    age(1:n) = [n-0.5:-1:1-0.5];
    age_av(n) = sum(temp(n,:)./sum(temp(n,:)).*age');
    
    clear age
end
bar(temp_c/10^12,'stacked');
hold on
bar(-temp_fgn/10^12,'stacked');
ax = gca;
ax.YLim = [-0.3 1.8];
ax.YTick = -0.3:0.3:1.8;
ax.YTickLabel = {'0.3','0','0.3','0.6','0.9','1.2','1.5','1.8'};
ax.XLim = [0.25 ts+0.75];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth=lw;
ax.FontSize=fs;
for i = 1:ts
    if i == 2
        text(i-0.5,sum(temp_c(i,:))/10^12+0.1,'1.0','FontSize',fs-2)
    elseif i == 11
        text(i-0.5,sum(temp_c(i,:))/10^12+0.1,'4.0','FontSize',fs-2)
    else
        text(i-0.5,sum(temp_c(i,:))/10^12+0.1,num2str(round(age_av(i),1)),'FontSize',fs-2)
    end
end
ax.ColorOrder = colormap;
% colorbar('Ticks',[0:0.25:1],'TickLabels',{'1995','2000','2005','2010','2015'},'Direction','reverse')
legend('','','','','','','','','','','','','','','','','','','','','','Location','eastoutside','Color','none','EdgeColor','none');
clear temp_c temp_fgn

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 5];
print('Fig2_legend','-dpng','-r300')

%% capital loop for EFk
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
load('EF_k_c','EFk_l_d_c_d_p')

r = 49;
s = 200;
ts = 21;
c = 31;


fs = 12;
lw = 1;
figure
temp = squeeze(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,c,:,:),3),5),6));
bar(temp'/10^12,'stacked')
ax = gca;
ax.YLim = [0 1.8];
ax.YTick = 0:0.6:1.8;
ax.YTickLabel = {'0','0.6','1.2','1.8'};
ax.XLim = [0.25 ts+0.75];
ax.XTick = 1:5:ts;
ax.XTickLabel = {'1995','2000','2005','2010','2015'};
ax.XColor = 'black';
ax.YColor = 'black';
ylabel('GHG emissions (Gt CO2eq.)');
% ylabel(sprintf('%s\n%s','GHG emissions','(Gt CO2eq.)'));
ax.LineWidth=lw;
ax.FontSize=fs;
hold on
plot([16.5 16.5],[ax.YLim(1) ax.YLim(2)],'--','Color',[0.7 0.7 0.7],'LineWidth',lw+0.5)
% colorbar('Ticks',[0:0.25:1],'TickLabels',{'1995','2000','2005','2010','2015'},'Direction','reverse')
legend('EF^K^C_l_p_=_1','EF^K^C_l_p_=_2','EF^K^C_l_p_>_2','Location','northwest','Color','none','EdgeColor','none');

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 3.5];
print('EFkc_by_loop','-dpng','-r300')


%% export V.S. import
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
load('EF_k_c','EFk_l_d_c_d_p')

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder')
load country49

reg_n = unique(country49(:,3));
% reg_n{5} = 'OECD';

r = 49;
s = 200;
ts = 21;
c = 31;
temp_exp = zeros(ts,length(reg_n),7);
temp_imp = zeros(ts,length(reg_n),7);
for i = 1:length(reg_n);
    fgn_c{i,1} = setdiff(find(contains(country49(:,3),reg_n(i))),[c]);
end

fs = 12;
lw = 1;
color = [];


%% in different countries
clear
cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
load('EF_k_c','EFk_l_d_c_d_p')
load EF_indicator
load ('EF_k_loop','Fk_n_t');

cd('C:\Users\YeQ\Documents\MATLAB\Capital paper\New folder')
load F_hh

r = 49
s = 200;
ts = 21;

fgn{1,1} = 1:27;
fgn{2,1} = 29;
fgn{3,1} = 34;
fgn{4,1} = 35;

fgn_n = {'EU27','USA','Brazil','India'};

index = 6;
fs = 12;
lw = 1;
color = [];
figure
for c = 1:length(fgn_n)
    subplot(2,4,c)
    temp_gfcf = squeeze(sum(sum(sum(EF_gfcf(index,:,fgn{c},:,:),3),4),5))';
    temp_efk = squeeze(sum(sum(sum(sum(Fk_n_t(:,:,:,:,(fgn{c}'-1)*s+(1:s)),1),3),4),5))';
    plot(horzcat(temp_gfcf,temp_efk)/10^12,'LineWidth',lw+0.5);
    
    ax = gca;
    ax.XLim = [1 ts];
    ax.XTick = 1:5:ts;
    ax.XTickLabel = {'1995','2000','2005','2010','2015'};
    ax.XColor = 'black';
    ax.YColor = 'black';
    ylabel('GHG emissions (Gt CO2eq.)');
    title(fgn_n(c))
    ax.LineWidth=lw;
    ax.FontSize=fs;
    if c == 4
        legend('EF^G^F^C^F','F^K','Location','northwest','Color','none','EdgeColor','none');
    end
    
       
    
    subplot(2,4,c+4)
    temp_ef_c = squeeze(sum(sum(sum(EF_fc(index,:,fgn{c},:,:),3),4),5))'+squeeze(sum(F_hh(index,:,fgn{c}),3))';
    temp_ef_kc_do = squeeze(sum(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,fgn{c},:,fgn{c}),1),3),4),5),6))';
    temp_ef_kc_out = squeeze(sum(sum(sum(sum(sum(EFk_l_d_c_d_p(:,:,:,fgn{c},:,setdiff([1:r],fgn{c})),1),3),4),5),6))';
    h=area(horzcat(temp_ef_c,temp_ef_kc_do,temp_ef_kc_out)/10^12);
    h(1).FaceColor = [0.8 0.8 0.8];
    h(2).FaceColor = [0.5 0.2 0.6];
    h(3).FaceColor = [0.72,0.27,1.00];
    ax = gca;
    ax.XLim = [1 ts];
    ax.XTick = 1:5:ts;
    ax.XTickLabel = {'1995','2000','2005','2010','2015'};
    ax.XColor = 'black';
    ax.YColor = 'black';
    ylabel('GHG emissions (Gt CO2eq.)');
    title(fgn_n(c))
    ax.LineWidth=lw;
    ax.FontSize=fs;
    if c == 4
        legend('EF^C','EF^K^C do','EF^K^C out','Location','northwest','Color','none','EdgeColor','none');
    end
end

cd('C:\Users\YeQ\Documents\MATLAB\Capital loop and current technology assumption')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 6];
print('Comparison other countries','-dpng','-r300')


c = 2;
temp_ef_c_do = squeeze(sum(sum(sum(EF_fc(index,:,fgn{c},fgn{c},:),3),4),5))'+squeeze(sum(F_hh(index,:,fgn{c}),3))';
temp_ef_c_out = squeeze(sum(sum(sum(EF_fc(index,:,fgn{c},setdiff([1:r],fgn{c}),:),3),4),5))'+squeeze(sum(F_hh(index,:,fgn{c}),3))';
sum(temp_ef_c_out)/sum(temp_ef_c_out+temp_ef_c_do)














