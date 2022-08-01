%% codes for global capital endogenization input-output model, based on EXIOBASE v3

% Created by: Quanliang Ye
% Created date: 31/11/2020
% Email add.: yequanliang1993@gamial.com


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
for lp = 1:10 % assume just to run 10 capital loops
    tic
    for t = 1:ts % production year of the depreciated capital goods
        fname = sprintf('L%d.mat',t); % Lenotief inverse matrix, which can be self-created
        load(fname);

        for n = t:ts % year of final consumption
            Fk = zeros(r,s*r); % impacts embedded in capital goods depreciated, allocated to year n
            
            if lp == 1
                f0 = squeeze(S(index,t,:))'; % direct environmental coefficient row vector
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











