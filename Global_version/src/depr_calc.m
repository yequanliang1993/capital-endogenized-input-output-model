%% This is code for the function of calculating capital depreciation

% Created by: Quanliang Ye
% Created date: 19/10/2019
% Email add.: yequanliang1993@gmail.com
% Main Contributor: dr. Ranran Wang, CML, Leiden University

%% 
function cd_ss = depr_calc(c,Ki,depr,c_axp,c_ixs,Yk,cfc,t,n)
s=200;
r=49;

gfcf_a = Ki;
gfcf_a(isnan(gfcf_a))=0;
gfcf_p = sum(reshape(Yk(t).Yk(:,c),[s,r]),2);
cfc_v = cfc(n).cfc((c-1)*s+1:c*s); 
    
% converting KLEMS capital investment data based on exiobase: from assets (8 or 10) to the capital goods producing sectors (200)
c_axp_v1=c_axp; % products producing assets: adjusted concordance - considering the distribution of sectors producing the same assets based on the production structure described in Yk (exiobase)
if size(Ki,1)>1
    if size(Ki,2) == size(depr,2)-1 || any(31:33==c) % asset details are available
        for a=1:size(Ki,2)
            temp = find(c_axp(a,:)==1);
            c_axp_v1(a,temp) = c_axp(a,temp).*gfcf_p(temp)'./sum(gfcf_p(temp));
        end
        clear temp;
        c_axp_v1(isnan(c_axp_v1))=0;
        gfcf_ap=gfcf_a*c_axp_v1;
        gfcf_ap=gfcf_ap.*gfcf_p'./sum(gfcf_ap,'omitnan');       % scale to Yk in exiobase (this adjusts for currency differences, too)
        gfcf_ap(isnan(gfcf_ap))=0;
        
        % adjust when Ki = 0 whereas gfcf of sector is not 0
        for j=1:s
            if gfcf_p(j)>0 && sum(gfcf_ap(:,j))==0 
                gfcf_ap(:,j) = gfcf_p(j).*(sum(gfcf_a,2)./sum(sum(gfcf_a)))';
            end
        end
        
        % check:
        %plot(sum(gfcf_ap,2)
        %hold on
        %plot(sum(gfcf_a,2))
        
        cd=zeros(size(Ki,1),s);
        for j = 1:s
            if sum(gfcf_ap(:,j))>0
                temp = find(c_axp(:,j)>0);
                temp1 = zeros(size(Ki,2),1); % address: one sector producing more than one assets
                if any(temp)
                    temp1(temp)=sum(gfcf_a(:,temp))./sum(sum(gfcf_a(:,temp)));
                    if sum(temp1,'omitnan')==0
                        temp1=repmat(1/size(Ki,2),[size(Ki,2),1]);
                    end
                else
                    temp1= sum(gfcf_a,1)./sum(sum(gfcf_a,1));
                end
                cdr = ((1-depr(:,1:end-1)).^(n-t)).*depr(:,1:end-1);
                cd(:,j) = sum(gfcf_ap(:,j)*c_axp(:,j)'*diag(temp1).*cdr,2);
            end
        end
        clear temp
    else
        gfcf_a = gfcf_a.*sum(gfcf_p)./sum(gfcf_a,'omitnan'); % scale to Yk in exiobase (this adjusts for currency differences, too)
        gfcf_ap = gfcf_a*(gfcf_p./sum(gfcf_p))'; % columns: capital producting sectors; rows: capital consuming (i.e. investing) sectors
        cd=zeros(size(Ki,1),s);
        for j=1:s
            cdr = ((1-depr(:,end)).^(n-t)).*depr(:,end);
            cd(:,j) = gfcf_ap(:,j).*cdr;
        end
    end
    
    c_ixs(isnan(c_ixs))=0;
    c_ixs_v1=c_ixs;
    
    for j = 1:size(Ki,1) % adjusting on the consumption side - going through each capital consuming sector
        temp = find(c_ixs(j,:)==1); % one KLEMS sector is often matched to multiple exiobase sectors
        c_ixs_v1(j,temp) = c_ixs(j,temp).*(cfc_v(temp)./sum(cfc_v(temp))); %assuming sectors with higher capital compensations also consume more capital goods
    end
    
    cd_ss = cd'*c_ixs_v1; % rows: capital producting sectors; columns: capital consuming (i.e. investing) sectors
    cd_ss(isnan(cd_ss))=0;
    cd_ss = cd_ss.*(sum(cd)'./sum(cd_ss,2)); % adjusts for the case where one exiobase sector is matched to multiple KLEMS sectors, by constraining on the amount of capital goods produced

else
    for a=1:size(Ki,2)
        temp = find(c_axp(a,:)==1);
        c_axp_v1(a,temp) = c_axp(a,temp).*gfcf_p(temp)'./sum(gfcf_p(temp));
    end
    clear temp;
    c_axp_v1(isnan(c_axp_v1))=0;
    gfcf_a = gfcf_a.*sum(gfcf_p)./sum(gfcf_a,'omitnan'); % scale to Yk in exiobase (this adjusts for currency differences, too)
    
    gfcf_ap = zeros(size(Ki,2),s);
    for a=1:size(Ki,2)
        gfcf_ap(a,:) = gfcf_a(a)*c_axp_v1(a,:);
    end
    
    gfcf_ap=gfcf_ap.*gfcf_p'./sum(gfcf_ap,1,'omitnan'); 
    gfcf_ap(isnan(gfcf_ap))=0;
    for j=1:s
        if gfcf_p(j)>0 && sum(gfcf_ap(:,j))==0
            gfcf_ap(:,j) = gfcf_p(j).*gfcf_a'./sum(sum(gfcf_a,'omitnan'));
        end
    end
    
    cd=zeros(1,s);
    cdr = ((1-depr(:,1:end-1)).^(n-t)).*depr(:,1:end-1);
    cdr = cdr'*c_ixs;
    cd=sum(gfcf_ap.*cdr);
    cd_ss = cd'*(cfc_v./sum(cfc_v));
end

