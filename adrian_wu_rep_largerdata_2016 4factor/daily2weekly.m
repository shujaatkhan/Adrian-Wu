%--------------------------------------------------------------------------
% Importing daily data and converting it to weekly data by selecting the 
% last observed data for a week.
%
% Daily data (complete)
% Start date: January 4, 1999
% End date:   October 23, 2015
%
% Weekly data (for replication)
% Start date: January 3, 2003 
% End date:   December 31, 2010
%
% NOTE: The first observation should be a Monday
%--------------------------------------------------------------------------
%[data_daily, d_dates] = xlsread('Data_1999.xls','Merged 1999 - Stata');
data_bond_all = xlsread([workdir,'/data/Data_1999_newformat.xlsx']);
data_3mon_all = xlsread([workdir,'/data/FRB_H15_T3.xlsx']);
data_vix_all = xlsread([workdir,'/data/vix_90_17.xlsx']);
data_bond_all(:,1) = x2mdate(data_bond_all(:,1));
data_3mon_all(:,1) = x2mdate(data_3mon_all(:,1));
data_vix_all(:,1) = x2mdate(data_vix_all(:,1));
% Converting the daily data to financial time series
data_bond_fints = fints(data_bond_all);
data_3mon_fints = fints(data_3mon_all);
data_vix_fints = fints(data_vix_all);
% Transforming the daily to weekly data
data_bond_wkly_tmp = toweekly(data_bond_fints);
data_3mon_wkly_tmp = toweekly(data_3mon_fints);
data_vix_wkly_tmp = toweekly(data_vix_fints);
% Converting weekly financial time series to matrix
data_bond_wkly = fts2mat(data_bond_wkly_tmp,1);
data_3mon_wkly = fts2mat(data_3mon_wkly_tmp,1);
data_vix_wkly = fts2mat(data_vix_wkly_tmp,1);
% Limiting weekly data to start-end dates
week_start_datenum = datenum(2003,1,3);
week_end_datenum = datenum(2009,12,31);
week_end_datenum = datenum(2016,12,30);
stind_bond = find(data_bond_wkly(:,1)==week_start_datenum);
endind_bond = find(data_bond_wkly(:,1)==week_end_datenum);
stind_3mon = find(data_3mon_wkly(:,1)==week_start_datenum);
endind_3mon = find(data_3mon_wkly(:,1)==week_end_datenum);
stind_vix = find(data_vix_wkly(:,1)==week_start_datenum);
endind_vix = find(data_vix_wkly(:,1)==week_end_datenum);
data_bond = data_bond_wkly(stind_bond:endind_bond,:);
data_3mon = data_3mon_wkly(stind_3mon:endind_3mon,:);
data_vix = data_vix_wkly(stind_vix:endind_vix,:);
mdates = data_bond(:,1);
T = size(data_bond,1);
clearvars -except data_bond data_3mon data_vix mdates T workdir rundate nx tau K plots_on fignum
%{
% xdates = data_daily(:,1);
% xdates = xdates-xdates(1)+1;
% xweeks = zeros(size(xdates));
% weeknum = 0;
% week_start_datenum = datenum(2003,1,3);
% week_end_datenum = datenum(2009,12,31);
% 
% for i=1:ceil(xdates(end)/7)
%     weeknum = weeknum+1;
%     start_ind = find(( xdates>=(7*(i-1)+1) & xdates<=(7*i) ),1,'first');
%     end_ind = find(( xdates>=(7*(i-1)+1) & xdates<=(7*i) ),1,'last');
%     xweeks(start_ind:end_ind) = weeknum*ones(end_ind-start_ind+1,1);
% end
% 
% for i=1:xweeks(end)
%     yind = find(xweeks==i,1,'last');
%     w_dat_all(i,:) = data_daily(yind,1);
%     y_dat_all(i,:) = data_daily(yind,2:K+1);
%     r_dat_all(i,:) = data_daily(yind,K+2:end);
% end
% 
% for i=1:K
%     yield_cell{i} = [y_dat_all(:,i) r_dat_all(:,i)];
% end
% 
% yields_all = cell2mat(yield_cell);
% T_all = size(y_dat_all,1);
% mdates_all = x2mdate(w_dat_all(:,1));
% 
% week_start_ind = find(mdates_all==week_start_datenum);
% week_end_ind = find(mdates_all==week_end_datenum);
% 
% mdates = mdates_all(week_start_ind:week_end_ind,1);
% yields = yields_all(week_start_ind:week_end_ind,:);
% T = size(yields,1);
% [data_3mtbill, dt_dates] = xlsread([workdir,'/data/Data_1999.xls'],'3mon Nominal weekly 1999-');
%}
for i=1:K
    bond_cell{i} = [data_bond(:,i+1) data_bond(:,i+9)];
end
bonds = cell2mat(bond_cell);
yield_3mon = data_3mon(:,2);
%--------------------------------------------------------------------------
% Running OLS regressions for each of the 8 maturities to get the residuals
% for the BEKK GARCH model.
% The OLS regressions are run with a constant.
%--------------------------------------------------------------------------
Y_res = nan(T-1,K);
R_res = nan(T-1,K);
bekk_cell = cell(1,K);
% %------
% bekk_init=[0.0196
%     0.0336
%     0.0000
%     0.2066
%     0.3659
%     0.9677
%     0.9101];
% %------

for i=1:K
    X = [bonds(1:end-1,2*i-1:2*i)];
    X1 = [ones(length(X),1),X];
    Y = bonds(2:end,2*i-1);
    R = bonds(2:end,2*i);
    Y_hat = X1*(inv(X1'*X1)*X1'*Y);
    R_hat = X1*(inv(X1'*X1)*X1'*R);
    Y_res(:,i) = Y - Y_hat;
    R_res(:,i) = R - R_hat;
    [parm_bekk,ll,ht,vcv,scores] = bekk([Y_res(:,i),R_res(:,i)],[],1,0,1,'Diagonal');
    for j=1:T-1
        bekk_temp = ht(:,:,j);
        bekk_cell{i}(j,:) = bekk_temp(:)';
    end
end

bekk_out = cell2mat(bekk_cell);
data = [mdates(2:end,1), yield_3mon(2:end,1), bonds(2:end,:), bekk_out];
% Need to convert the mdates to excel datenums; otherwise, export problems
data_x = data; data_x(:,1)=m2xdate(data_x(:,1));
data_vix_x = data_vix; data_vix_x(:,1)=m2xdate(data_vix_x(:,1));
% mdates_yields = [mdates, bonds];
% xlswrite('data_all_0310.xls',data,'data');
% xlswrite('data_all_0310.xls',mdates_yields,'mdates_yields');
% Contains: all yield data and bekk output
csvwrite(['./data/data_all_',num2str(rundate(1)),'_',num2str(rundate(2)),'_',num2str(rundate(3)),'.csv'],data_x);
csvwrite(['./data/data_vix_',num2str(rundate(1)),'_',num2str(rundate(2)),'_',num2str(rundate(3)),'.csv'],data_vix_x(2:end,:));
% Contains: matlab dates data and yields (ex. 3mon)
% csvwrite(['./data/mdates_yields_',num2str(rundate(1)),'_',num2str(rundate(2)),'_',num2str(rundate(3)),'.csv'],mdates_yields);
%--------------------------------------------------------------------------
% Summary statistics for variances and covariances of yields (*10^2)
%--------------------------------------------------------------------------
bekk_100 = 100*bekk_out;
for i=1:K
    j = 4*i-3;
    % stats         = [maturity   mean   median   std   max   min]
    stats_nvar(i,:) = [i+2,mean(bekk_100(:,j)),median(bekk_100(:,j)), ...
                         std(bekk_100(:,j)),max(bekk_100(:,j)), ...
                         min(bekk_100(:,j))];
    stats_cov(i,:) = [i+2,mean(bekk_100(:,j+1)),median(bekk_100(:,j+1)), ...
                         std(bekk_100(:,j+1)),max(bekk_100(:,j+1)), ...
                         min(bekk_100(:,j+1))];                 
    stats_rvar(i,:) = [i+2,mean(bekk_100(:,j+3)),median(bekk_100(:,j+3)), ...
                         std(bekk_100(:,j+3)),max(bekk_100(:,j+3)), ...
                         min(bekk_100(:,j+3))];
end
%--------------------------------------------------------------------------
