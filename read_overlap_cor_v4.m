%Read Overlap Corrections for paper
% Calculate a model according to the internal temperature.
% V4: Work for KSE, PAY and SIRTA
% Remove plots not in paper
clear variables;clc;close all;

set(0,'DefaultFigureVisible','on')

%% INPUTS
% Time range for loading Overlap correction dataset
% PAY 01/02/2013 to 21/11/2014
% PAY TUB140016 08/11/2015 to ...
% KSE TUB140005 20/02/2015 to 15/12/2015  
% KSE TUB120011 01/05/2016 to 25/4/2017 (Test Case 17/08/2016)
% KSE TUB120011 25/4/2017 ... Same TUB but different Overlap function after moving instrument (Test Case 17/08/2017)
% KSE TUB170001 25/3/2019 ... 

% 

info.start_day  = 25;
info.start_month= 3;
info.start_year = 2019;

info.end_day  =  1;
info.end_month=  11;
info.end_year =  2020;


% Time range to apply the correction
info_test.start_day  = 1;
info_test.start_month= 4;
info_test.start_year = 2019;

info_test.end_day  =  1;
info_test.end_month=  4;
info_test.end_year =  2019;

station='kse';
suffix = ''; %'' _after_moving_2017

min_nb_good_samples_after_outliers_removal  = 10;

info.plot=1; % Plot all data?
info.save_plots=1;
disp_clouds = true;
disp_pbl = true;
info_reloading=1; %(Reload all data?)

create_netcdf=1;

switch station
    case 'pay'
        info.chm='CHM120106';
%         info.tub='TUB120011';
                info.tub='TUB120011';
        
    case 'kse'
        info.chm='CHM130104';
%         info.tub='TUB140005';
%                 info.tub='TUB120011';
                info.tub='TUB170001';

    case 'lindenberg'
        info.chm='CHM100110';
        info.tub='TUB120001';
    case 'lindenberg2'
        info.chm='CHM140101';
        info.tub='TUB140006';
           case 'SIRTA'
        info.chm='CHM150101';
        info.tub='TUB140013';
    otherwise
        error('Define ceilometer info')
end

corrections_to_analyze = 'good_enough';%{'all','good_enough','well_trusted'};
% folder_output =['../Outputs/' station '/'];
folder_output =['M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites/' station '/'];

folder_ov_ref = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\overlap_functions_Lufft\';
folder_results_algo = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\';
switch station
    case 'pay'
        folder_ncdata = 'M:\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
        if strcmp(info.tub,'TUB120011')
            folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections\';
        else
            folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites\pay\';
        end
    case 'kse'
        folder_ncdata = 'M:\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
        folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites\kse\';
    case 'lindenberg'
        folder_ncdata = 'D:\Ceilometer\CEILINEX\chm15k_100110\';
        folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites\Lindenberg\';
    case 'lindenberg2'
        folder_ncdata = 'D:\Ceilometer\CEILINEX\chm15k_140101\';
        folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites\Lindenberg\';
   
       case 'SIRTA'
        folder_ncdata = 'C:\DATA\Ceilometer\CHM15k\DATA\st\';
        folder_corrections = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections_other_sites\SIRTA\';

    otherwise
        error('Define folder with overlap correction')
end

switch station
    case 'pay'
        if strcmp(info.tub,'TUB120011')
            [results_flag,results_time] = xlsread([folder_results_algo 'results_automatic_algo_paper.xlsx'],1,'A2:B755');
            results_time = datenum(results_time,'dd.mm.yyyy');
            
            ind_dates_nodata = results_flag==-2;
            ind_dates_nodetection = results_flag==-1;
            ind_dates_badquality = results_flag==0;
            ind_dates_acceptablequality = results_flag==1;
            ind_dates_goodquality = results_flag==2;
            ind_dates_verygoodquality = results_flag==3;
            
            list_dates_nodetection = cellstr(datestr(results_time(ind_dates_nodetection),'yyyymmdd'));
            list_dates_badquality = cellstr(datestr(results_time(ind_dates_badquality),'yyyymmdd'));% not better than Lufft
            list_dates_acceptablequality = cellstr(datestr(results_time(ind_dates_acceptablequality),'yyyymmdd'));
            list_dates_goodquality = cellstr(datestr(results_time(ind_dates_goodquality),'yyyymmdd'));
            list_dates_verygoodquality = cellstr(datestr(results_time(ind_dates_verygoodquality),'yyyymmdd'));
            
            % list_dates_outliers_temperature = {'20130406','20130601','20140306','20140607','20140608','20140811','20140925'};
            list_dates_outliers_temperature = {'20130316','20130512','20130523','20130601','20130605','20130627','20130630','20130725',...
                '20130727','20130728','20130906','20131229','20140306','20140313','20140214','20140324','20140325','20140327',...
                '20140417','20140522','20140811','20140823','20140925'};
        else
            list_dates_nodetection = {};
            list_dates_badquality = {};% not better than Lufft
            list_dates_acceptablequality = {};
            list_dates_goodquality = {};
            list_dates_verygoodquality = {};
            list_dates_outliers_temperature = {};
        end
    case 'kse'
        if strcmp(info.tub,'TUB120011')
            list_dates_nodetection = {};
            list_dates_badquality = {'20160526','20160802','20170603','20160924','20170908'};% not better than Lufft
            list_dates_acceptablequality = {};
            list_dates_goodquality = {};
            list_dates_verygoodquality = {};
            list_dates_outliers_temperature = {'20160526','20170603','20170803','20160902','20160707','20170908',...
                '20170827','20170710','20170823','20160814'};
             elseif strcmp(info.tub,'TUB170001')
            list_dates_nodetection = {};
            list_dates_badquality = {'20190401'};% not better than Lufft
            list_dates_acceptablequality = {};
            list_dates_goodquality = {};
            list_dates_verygoodquality = {};
            list_dates_outliers_temperature = {};
        
        else
             list_dates_nodetection = {};
            list_dates_badquality = {};% not better than Lufft
            list_dates_acceptablequality = {};
            list_dates_goodquality = {};
            list_dates_verygoodquality = {};
            list_dates_outliers_temperature = {'20150416','20150705','20150724','20150803','20150908','20150909','20151104'};

        end
    case 'lindenberg'
        list_dates_nodetection = {};
        list_dates_badquality = {};% not better than Lufft
        list_dates_acceptablequality = {};
        list_dates_goodquality = {};
        list_dates_verygoodquality = {};
        list_dates_outliers_temperature = {};
    case 'lindenberg2'
        list_dates_nodetection = {};
        list_dates_badquality = {};% not better than Lufft
        list_dates_acceptablequality = {};
        list_dates_goodquality = {};
        list_dates_verygoodquality = {};
        list_dates_outliers_temperature = {};
    case 'SIRTA'
        list_dates_nodetection = {};
        list_dates_badquality = {};% not better than Lufft
        list_dates_acceptablequality = {};
        list_dates_goodquality = {};
        list_dates_verygoodquality = {};
        list_dates_outliers_temperature = {};
    otherwise
        error('Define flags')
end



%% Loading
if info_reloading == 1 || exist(['all_correction_' station '_' info.tub  suffix '.mat'],'file')==0
    %% Read: all files with a final overlap correction output
    time_vec=datenum(info.start_year,info.start_month,info.start_day):datenum(info.end_year,info.end_month,info.end_day);
    k=1;
    time0 = [];
    quality0 = [];
    nb_candidates0 = [];
    overlap_cor0 = [];
    temp_mean0 = [];
    voltage_PM_mean0 = [];
    for t=1:length(time_vec)
        file=['result_ovp_fc_' info.chm info.tub '_' datestr(time_vec(t),'yyyymmdd') '.mat'];

        folder_day=[folder_corrections datestr(time_vec(t),'yyyy/mm/')];
        
        if exist([folder_day file],'file')>0
            disp([folder_day file]);
            load([folder_day file],'result');
            if isempty(result.index_final) || length(result.index_final)<min_nb_good_samples_after_outliers_removal
                continue;
            end
            overlap_cor0(k,:)=nanmedian(result.ovp_fc(:,result.index_final),2)';
            time0(k)=time_vec(t);
            
            % Get quality
            switch datestr(time_vec(t),'yyyymmdd')
                case list_dates_badquality
                    quality0(k) = 0;
                case list_dates_acceptablequality
                    quality0(k) = 1;
                case list_dates_goodquality
                    quality0(k) = 2;
                case list_dates_verygoodquality
                    quality0(k) = 3;
                otherwise
                    disp('No Quality classification')
                    quality0(k) = 1;
            end
            % Get total number of candidates that passed the fit tests
            nb_candidates0(k) = size(result.ovp_fc,2);
            
            % Get temperature
            
            station_str = station;
            if strcmp(station,'pay')
                if time_vec(t)<=datenum(2013,2,11)
                    station_str = 'NN';
                elseif time_vec(t)<=datenum(2013,6,13)
                    station_str = 'payerne';
                end
            end
            if strcmp(station,'lindenberg2')
                station_str = 'lindenberg';
            end
            
            if strcmp(station,'pay') || strcmp(station,'kse')
                file=[datestr(time_vec(t),'yyyymmdd') '_' station_str '_' info.chm '_000.nc'];
                folder_day=[folder_ncdata datestr(time_vec(t),'yyyy/mm/') ] ;
                
            elseif strcmp(station,'SIRTA')
                folder_day=[folder_ncdata datestr(time_vec(t),'yyyy/mm/dd/') ] ;
                list = dir([folder_day 'chm15k_*.nc']);
                file=list.name;
                
            else
                file=[datestr(time_vec(t),'yyyymmdd') '_' station_str '_' info.chm '.nc'];
                folder_day=[folder_ncdata datestr(time_vec(t),'yyyy/mm/') ] ;
                
            end
            disp([folder_day file])
            if exist([folder_day file],'file')==0
                error('Missing NetCDFfile file')
            end
            
            %convert scale factor
            %             scale_factor = ncreadatt([folder_day file],'temp_int','scale_factor');
            %             if(~isnumeric(scale_factor))
            %                 ncwriteatt([folder_day file],'temp_int','scale_factor',1./str2double(scale_factor));
            %             end
            
            temp_int = ncread([folder_day file],'temp_int');
            nn1 = ncread([folder_day file],'nn1');
                        zenith_angle = ncread([folder_day file],'zenith');

            time_nc = datenum(1904,1,1)+ncread([folder_day file],'time')/3600/24;
            temp_tmp = NaN(length(result.index_final),1);
            voltage_tmp = NaN(length(result.index_final),1);
            range = ncread([folder_day file],'range');
            for i=1:length(result.index_final)
                indt = time_nc>=result.time_start(result.index_final(i)) & time_nc<=result.time_end(result.index_final(i));
                temp_tmp(i) = nanmean(temp_int(indt));
                voltage_tmp(i) = nanmean(nn1(indt));
                
            end
            
            %         temp_mean(k)=mean(ncread([folder_day file],'temp_int'));
            %         voltage_PM_mean(k)=mean(ncread([folder_day file],'nn1'));
            temp_mean0(k)=nanmedian(temp_tmp);
            voltage_PM_mean0(k)=nanmedian(voltage_tmp);
            zenith_angle_mean(k) = nanmedian(zenith_angle);
            
            
            k=k+1;
        else
            disp('No file:')
            disp([folder_day file])
        end
    end
    
    %% Read: overlap correction and scaling in cfg file
    list = dir([folder_ov_ref,info.tub,'_*.cfg']);
    if isempty(list)
        answer=questdlg({'No reference overlap function found','Do you want to calibrate with TUB120011 overlap function'},...
            'No reference overlap function','Yes','No','No');
        if strcmp(answer,'Yes')
            warning('No overlap function found. Using TUB120011 overlap function')
            fid = fopen('TUB120011_20121112_1024.cfg');
            ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
            overlap_ref = cell2mat(ov_cell);
        else
            error('No overlap function found. Ask the corresponding overlap to the manufacturer')
        end
    else
        file=[folder_ov_ref list(1).name];
        fid=fopen(file);
        disp(file);
        data=textscan(fid,'%f','delimiter','\t','headerlines',1);
        overlap_ref=data{1};
    end
    save(['all_correction_' station '_' info.tub '.mat'])
else
    corrections_to_analyze_tmp=corrections_to_analyze;
    info_test_tmp=info_test;
    info_tmp=info;
    
    disp('Load all correction in mat file')
    load(['all_correction_' station '_' info.tub '.mat'])
    corrections_to_analyze=corrections_to_analyze_tmp;
    info=info_tmp;
    info_test=info_test_tmp;
    
    clear info_test_tmp info_tmp
    
end
%% Select: set of overlap functions to analyse
disp('Select functions to analyse')
% all
index_all = 1:length(time0);

% all, without bad quality and without T outliers
index_good_enough=[];
for i=1:length(time0)
    if ~any(strcmp(list_dates_badquality,datestr(time0(i),'yyyymmdd'))) && ~any(strcmp(list_dates_outliers_temperature,datestr(time0(i),'yyyymmdd')))
        index_good_enough(end+1)=i;
    end
end
% only well-trusted (very good quality and not T outliers)
index_well_trusted=[];
for i=1:length(time0)
    if any(strcmp(list_dates_verygoodquality,datestr(time0(i),'yyyymmdd'))) && ~any(strcmp(list_dates_outliers_temperature,datestr(time0(i),'yyyymmdd')))
        index_well_trusted(end+1)=i;
    end
end

if strcmp(corrections_to_analyze,'all')
    index_of_interest = index_all;
elseif strcmp(corrections_to_analyze,'good_enough')
    index_of_interest = index_good_enough;
elseif strcmp(corrections_to_analyze,'well_trusted')
    index_of_interest = index_well_trusted;
else
    error('please specify which corrections to analyse')
end

time = time0(index_of_interest);
quality = quality0(index_of_interest);
nb_candidates = nb_candidates0(index_of_interest);
overlap_cor = overlap_cor0(index_of_interest,:);
temp_mean = temp_mean0(index_of_interest);
voltage_PM_mean = voltage_PM_mean0(index_of_interest);


%% FIG 4 : all overlap functions
disp('Plot fig 4: All overlap functions')
T = temp_mean-273.15;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

figure('Name','fig 4')
set(gcf,'DefaultAxesColorOrder',RGB)
hold on;

patcher(range',nanmin(overlap_cor,[],1),nanmax(overlap_cor,[],1),0.25*ones(1,3),[],'FaceAlpha',0.25);
for i=1:length(time)
    plot(range,overlap_cor(i,:),'displayname',datestr(time(i)));
end
plot(range,overlap_ref,'--k','LineWidth',2)

hc = colorbar;
ylabel(hc,'Median internal Temperature [?C]');
caxis([15,40]);
colormap(jet);
ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200)
grid on;box on;
xlabel('Range [m]')
ylabel('Overlap function')


%% FIG 6a: relative difference between corrected and uncorrected signal
disp('plot fig 6a: Relative difference')
T = temp_mean-273.15;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
% RGB = interp1(1:length(jet),jet,Treduced);
RGB = interp1(1:length(jet),jet,T);

relative_difference=NaN(length(time),length(range));
for i=1:length(time)
    relative_difference(i,:)=(-1./overlap_ref'+1./overlap_cor(i,:))./(1./overlap_ref')*100;
end

figure
set(gcf,'defaultAxesColorOrder',RGB)
colorvec = jet(length(time));
hold on
for i=1:length(time)
    h=scatter(relative_difference(i,:),range,...
        ones(size(overlap_cor(i,:)))*12,...
        repmat(temp_mean(i)-273.15,size(overlap_cor(i,:))),...
        'filled','displayname',datestr(time(i)));
    hp = plot(relative_difference(i,:),range,...
        '.-','displayname',datestr(time(i)));
 end
hc = colorbar;
ylabel(hc,'Median internal Temperature [?C]');
ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-75 75]);
caxis([15 40]);
grid on
box on

%% plot rel diiff changinf with time
disp(' Relative difference with time')
figure
colorvec = jet(length(time));
hold on
for i=1:length(time)
    hp = plot(relative_difference(i,:),range,...
        '.-','displayname',datestr(time(i)),'color',colorvec(i,:));
end
hc = colorbar;
ylabel(hc,'Time');
ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-75 75]);
grid on
box on

% %% plot rel diiff changing with specific time
% disp(' Relative difference with specific time')
% figure
% hold on
% index1 = time < datenum( 2016,10,04);
% index2 = time >= datenum( 2016,10,04) & time <= datenum( 2017,04,25);
% index3 = time > datenum( 2017,04,25);
% 
% 
% % plot(relative_difference(index1,:),range,'k')
% % h1 = plot(relative_difference(~index3,:),range,'r');
% h2= plot(relative_difference(index3,:),range,'b');
% legend([h1(1),h2(1)],'Before 20140425','After 20170425');
% ylim([0,1200])
% set(gca,'ytick',0:200:1200)
% xlabel({'Relative Difference between';'corrected and uncorrected signal [%]' })
% ylabel('Range [m]')
% colormap(jet)
% xlim([-75 75]);
% grid on
% box on


%% FIG 5: Scatter T vs. L2 norm REl. diff
disp('Plot fig 5: Scatter T')
ti = time;
T = temp_mean;
ov_cor = overlap_cor(:,overlap_ref>0.05);
ov = overlap_ref(overlap_ref>0.05);
relative_difference_selected=relative_difference(:,overlap_ref>0.05)/100;
disp_text = false;

A = T-273.15;

% To be coherent with fig 4 in paper Use relative signal
B = sqrt(trapz(abs(relative_difference_selected').^2));


p = polyfit(A,B,1);

SStot = sum((B'-repmat(mean(B),length(B),1)).^2);
SSres = sum((B'-polyval(p,A')).^2);
R2 = 1-SSres/SStot;
RMSE = sqrt(SSres/length(A));
RMSErel = 100*RMSE / sqrt(sum(polyval(p,A').^2)/length(A));

fz = 18;
hf=figure('Units','Normalized','Position',[0.1 0.3 0.75 0.5]);

scatter(A,B,'filled');
if disp_text
    hold on;
    for j=1:length(ti)
        text(A(j),B(j),datestr(ti(j),'dd.mm.yyyy HH:MM'),'HorizontalAlignment','center','FontSize',8);
    end
end
xlim([0 50]);
h_fit = line(xlim,polyval(p,xlim),'Color','k','LineWidth',2,'LineStyle','--');
xlabel('Median internal Temperature [?C]','FontSize',fz);
ylabel({'L^2 Norm of the Relative Difference';'between corrected and uncorrected signal' },'FontSize',fz);
ylim([0 2]);
grid on;
legend(h_fit,['R^2=' num2str(R2) ', RMSE=' num2str(RMSE)]);
set(gca,'FontSize',fz);

%% Fit: Relative difference according to temperature
a = zeros(size(range));
b = zeros(size(range));
r2 = zeros(size(range));
rmse = zeros(size(range));

for j=1:find(range<=1200,1,'last')
    val = overlap_cor(:,j);
    val(val<=0) = eps;
    rel_diff = (-1./overlap_ref(j)+1./val)./(1./overlap_ref(j))*100;
    
    y_fit = rel_diff;
    x_fit = temp_mean(:)-273.15;
    
    x_fit = x_fit(~isnan(y_fit) & ~isinf(y_fit));
    y_fit = y_fit(~isnan(y_fit) & ~isinf(y_fit));
    
    dT = 1;
    [counts,centers] = hist(x_fit,15:dT:40);
    weights = zeros(length(counts),1);
    weights(counts>0) = max(counts)./counts(counts>0);
    weights_fit = zeros(length(x_fit),1);
    for l=1:length(centers)
        indx = x_fit>centers(l)-dT/2 & x_fit<=centers(l)+dT/2;
        weights_fit(indx) = weights(l).^1;
    end
    
    [fo,gof] = fit(x_fit,y_fit,'poly1');
    a(j) = fo.p1(1);
    b(j) = fo.p2(1);
    r2(j) = gof.rsquare;
    rmse(j) = gof.rmse;
end



%% Figure 6b: relative difference model
disp('plot fig 6b: Relative difference using model')

temp_vector = 15:0.1:40;

color_vector=jet(length(temp_vector));
figure
hold on
ov_reconstructed=ones(length(temp_vector),length(range));
for i=1:length(temp_vector)
    rel_diff = NaN(length(range),1);
    for j=1:length(range)
        rel_diff(j) = polyval([a(j),b(j)],temp_vector(i));
    end
    
    %Remove abberant fit results
    if any(abs(a)>10 | abs(b)>100,1)
        rel_diff(range<=range(find(abs(a)>10 | abs(b)>100,1,'last'))) = NaN;
    end
    
    h=scatter(rel_diff,range,...
        ones(length(range),1),repmat(temp_vector(i),length(range),1),...
        'filled','displayname',num2str(temp_vector(i)));
    h=plot(rel_diff,range,'color',color_vector(i,:));
    
end
hc = colorbar;
ylabel(hc,'Internal Temperature [?C]');

ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected (with model) and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-75 75]);
grid on
box on

%% Filtered a and b coefficents
a_filtered=a;
b_filtered=b;
if any(abs(a)>10 | abs(b)>100,1)
    %Remove abberant fit results
    ind_last_as_reference = find(range <=range(find(abs(a)>10 | abs(b)>100,1,'last')),1,'last');
    if ~isempty(ind_last_as_reference)
        a_filtered(range<=range(ind_last_as_reference)) = 0;
        b_filtered(range<=range(ind_last_as_reference)) = 0;
    end
end
% Above 1200m, trust lufft values
a_filtered(range>=1200) = 0;
b_filtered(range>=1200) = 0;

%% write netCDf file of Temp model
if create_netcdf==1
    file=[ folder_output 'Overlap_correction_model_' station '_' info.tub suffix '.nc'];
    if exist(file,'file')
        warning('deleting exisiting NetCDF model file')
        delete(file)
    end
    if exist(folder_output,'dir')==0
        disp(['create directory: ' folder_output])
        mkdir(folder_output)
    end
    
    disp(['Creating ' file])
    nccreate( file,'a','dimensions',{'range', length(range)})
    ncwriteatt( file,'a','long_name','Results of fit (difference = a *Temperature +b )')
    ncwrite( file,'a',a_filtered)
    
    nccreate( file,'b','dimensions',{'range', length(range)})
    ncwriteatt( file,'b','long_name','Results of fit (difference = a * Temperature +b )')
    ncwrite( file,'b',b_filtered)
    
    nccreate( file,'range','dimensions',{'range', length(range)})
    ncwriteatt( file,'range','long_name','Altitude above ground (m)')
    ncwrite( file,'range',range)
    
    nccreate( file,'overlap_ref','dimensions',{'range', length(range)})
    ncwriteatt( file,'overlap_ref','long_name','Reference overlap function')
    ncwrite( file,'overlap_ref',overlap_ref)
    
    nccreate(file,'time','dimensions',{'time',length(time)})
    ncwriteatt( file,'time','long_name','Time with successful calibrations')
    ncwriteatt( file,'time','units','days since 1970-01-01 00:00:00.000 (UTC)')
    ncwrite( file,'time',time-datenum(1970,1,1))
    
    
    ncwriteatt( file,'/','device_name',info.chm)
    ncwriteatt( file,'/','serlom',info.tub)
    ncwriteatt( file,'/','Time_range',['From ' datestr(datenum(info.start_year,info.start_month,info.start_day)) ...
        ' to ' datestr(datenum(info.end_year,info.end_month,info.end_day))])
    
    
    
    
    ncwriteatt(file,'/','description', ['Output for overlap artefact correction. ' ...
        'Use:   Dif(z)= a (z)* T + b(z) and '...
        'Overlap_corrected(z) = 1./ (Dif(z) /100/overlap_ref(z) + 1./overlap_ref(z));'])
end

%% Load: test data
time_vec_test=datenum(info_test.start_year,info_test.start_month,info_test.start_day):...
    datenum(info_test.end_year,info_test.end_month,info_test.end_day);
for t=1:length(time_vec_test)
    date = datestr(time_vec_test(t),'yyyymmdd');
    
    if ~any(time==datenum(date,'yyyymmdd'))
        warning('No Overlap correction for the selected day')
    end
    
    %     [chm,chminfo] = get_chm15k_from_files('pay',[date,'000000'],datestr(datenum(date,'yyyymmdd')+0.99999,'yyyymmddHHMMSS'),folder_ncdata);
    [chm,chminfo]=readcorrectlyncfile3(station_str,date,folder_ncdata);
    if isempty(chm)
        warning('no file')
        continue
    end
    
    RCS_raw = chm.beta_raw;
    if any(time==datenum(date,'yyyymmdd'))
        RCS_dailycorrection = chm.beta_raw./repmat(overlap_cor(time==datenum(date,'yyyymmdd'),:)',1,size(chm.beta_raw,2)).*repmat(overlap_ref,1,size(chm.beta_raw,2));
    else
        RCS_dailycorrection=NaN(size(    RCS_raw));
    end
    %% Calc Grad
    grad_raw = NaN(size(RCS_raw));
    for j=1:length(chm.time)
        grad_raw(:,j) = 1/(2*chm.range_gate)*conv(RCS_raw(:,j),[1 0 -1],'same');
    end
    
    grad_dailycorrection = NaN(size(RCS_dailycorrection));
    if any(time==datenum(date,'yyyymmdd'))
        for j=1:length(chm.time)
            grad_dailycorrection(:,j) = 1/(2*chm.range_gate)*conv(RCS_dailycorrection(:,j),[1 0 -1],'same');
        end
    end
    
    %% calculate RCS and grad(RCS) with the temperature model
    ov_rec_all = ones(size(RCS_raw));
    for i=1:find(chm.range<=1200,1,'last')
        rel_diff = polyval([a_filtered(i),b_filtered(i)],chm.temp_int-273.15);
        ov_rec_all(i,:) = 1./ (rel_diff/100/overlap_ref(i) + 1./overlap_ref(i));
    end
    %     ov_rec_all(find(chm.range<=1200,1,'last')+1:end,:) = repmat(overlap_ref(find(chm.range<=1200,1,'last')+1:end),1,size(RCS_raw,2));
    
    RCS_corr = RCS_raw./ov_rec_all.*repmat(overlap_ref,1,size(RCS_raw,2));
    
    % Calc Grad
    grad_corr = NaN(size(RCS_corr));
    for j=1:length(chm.time)
        grad_corr(:,j) = 1/(2*chm.range_gate)*conv(RCS_corr(:,j),[1 0 -1],'same');
    end
    
    %% Calculate PBL with simple code
    [gradients_raw,~] = simplePBLdetection(chm,overlap_ref,overlap_ref);
    if any(time==datenum(date,'yyyymmdd'))
        [gradients_dailycorrection, ~] = simplePBLdetection(chm,overlap_ref,overlap_cor(time==datenum(date,'yyyymmdd'),:)');
    else
        gradients_dailycorrection=NaN(size(gradients_raw));
    end
    [gradients_corr,xtime] = simplePBLdetection(chm,overlap_ref,ov_rec_all);
    
    %% Write text file with PBL heigt
    disp('Writing text file with PBL height')
    M=[datevec(datenum(date,'yyyymmdd')+xtime) gradients_raw' gradients_dailycorrection' gradients_corr'];
    dlmwrite([folder_output '/PBL_' station '_' date '_' corrections_to_analyze  '.csv' ],M)
    
    if info.plot==1
        %% plot internal temperature
        figure;
        plot(chm.time,chm.temp_int-273.15);
        hold on
        plot(chm.time,chm.temp_ext-273.15);
        datetick;
        box on;grid on;
        xlabel('Time [UT]');ylabel('Temperature [?C]');title(['pay ' date]);
        legend('Internal','External')
        
        if any(time==datenum(date,'yyyymmdd'))
            
            %% load and plot daily overlap correction
            file=['result_ovp_fc_' info.chm info.tub '_' date '.mat'];
            folder_day=[folder_corrections datestr(datenum(date,'yyyymmdd'),'yyyy/mm/')];
            if exist([folder_day file],'file')>0
                disp([folder_day file]);
                load([folder_day file],'result');
                if ~isempty(result.index_final) && length(result.index_final)>=min_nb_good_samples_after_outliers_removal
                    time_start = result.time_start(result.index_final);
                    time_end = result.time_end(result.index_final);
                    range_start = result.range_start(result.index_final);
                    range_end = result.range_end(result.index_final);
                end
            end
            
            
            %% Plot 9: Plot all daily overlaps and manufacturer
            figure
            hold on
            patcher(range,nanmin(result.ovp_fc(:,result.index_final),[],2),nanmax(result.ovp_fc(:,result.index_final),[],2),...
                0.25*ones(1,3),[],...
                'FaceAlpha',0.5,'edgecolor',[0.5 0.5 0.5]);
            h1=plot(range,median(result.ovp_fc(:,result.index_final),2),'color','k','linewidth',2);
            
            h2=plot(range,overlap_ref,'--k','LineWidth',2);
            legend([h1,h2],{'Corrected','Manufacturer'})
            ylim([0 1.2]);
            xlim([0 1200]);
            set(gca,'xtick',0:200:1200)
            grid on;box on;
            xlabel('Range [m]')
            ylabel('Overlap function')
            title(['overlap functions ' date])
            
            %% plot 10: Fit Example
            figure
            subplot(1,2,1)
            hold on
            ylim([0 1000])
            xlim([ 4 6])
            plot(xlim,[min(result.range_start(result.index_final)) min(result.range_start(result.index_final))],'--k')
            plot(xlim,[max(result.range_end(result.index_final)) max(result.range_end(result.index_final))],'--k')
            xlabel('log_{10}(abs(\beta_{raw}))')
            ylabel('Range [m]')
            box on
            grid on
            
            subplot(1,2,2)
            hold on
            ylim([0 1000])
            xlabel('f_c')
%             plot([1 1],ylim,'-k')
            box on
            grid on
            switch station
                case 'pay'
                    station_name='Payerne';
                case'kse'
                    station_name='Kleine Scheidegd';
                otherwise
                    station_name=station;
            end
            
            dif_all=NaN(length(result.time_start(result.index_final)),length(chm.range));
            lrcs_all=NaN(length(result.time_start(result.index_final)),length(chm.range));
            for i=1:length(result.time_start(result.index_final))
                index_time=chm.time>=result.time_start(result.index_final(i)) & chm.time<result.time_end(result.index_final(i));
                index_range =chm.range>=result.range_start(result.index_final(i)) & chm.range<result.range_end(result.index_final(i));
                p = polyfit(chm.range(index_range),mean(log10(abs(chm.beta_raw(index_range,index_time))),2),1);
                dif = mean(log10(abs(chm.beta_raw(:,index_time))),2)-polyval(p,chm.range);
                dif( chm.range>=result.range_end(result.index_final(i)))=0;
                dif_all(i,:)=dif;
                lrcs_all(i,:)=log10(abs(mean(chm.beta_raw(:,index_time),2)));
                overlap_corr_factor = 10.^(dif);
                subplot(1,2,1)
%                 plot(log10(abs(mean(chm.beta_raw(:,index_time),2))),chm.range,'.k')
                h=plot([p(2) p(1)*1000+p(2)],[0 1000],'--r');
                legend(h,'Linear fit')
                
%                 subplot(1,2,2)
%                 plot(1+dif,chm.range,'.k')
            end
            subplot(1,2,1)
            errorshade4(median(lrcs_all),chm.range',...
                nanmin(lrcs_all),nanmax(lrcs_all),...
                zeros(size(chm.range')), zeros(size(chm.range')),'k');
            plot(median(lrcs_all),chm.range,'color','k','linewidth',2);
            
            subplot(1,2,2)
            errorshade4(median(dif_all)+1,chm.range',...
                nanmin(dif_all+1),nanmax(dif_all+1),...
                zeros(size(chm.range')), zeros(size(chm.range')),'k');
            plot(median(dif_all)+1,chm.range,'color','k','linewidth',2);
            
            
            suptitle([station_name ': ' date])  

        else
            warning('No daily correction file')
        end
        
        %% Plot PBL
        figure;
        ymax=2500;
        subplot(3,1,1);
        gradients = gradients_raw;
        correction_type = 'raw';
        list_plots = [];
        list_legends = {};
        hold on;
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        list_legends{end+1} = 'clouds';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'or','markerfacecolor','r','markersize',10);
        list_legends{end+1} = '1. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
        list_legends{end+1} = '2. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
        list_legends{end+1} = '3. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+k');
        list_legends{end+1} = 'min. gradient';
        legend(list_plots,list_legends,'location','eastoutside');
        ylim([0 ymax]);
        xlim([datenum(date,'yyyymmdd') datenum(date,'yyyymmdd')+1]);
        set(gca,'xtick',datenum(date,'yyyymmdd'):4/24:datenum(date,'yyyymmdd')+1);
        datetick('x','HH:MM','keepticks','keeplimits');
        title(['pay ' date ' - ' correction_type]);
        grid on;box on;
        xlabel('Time [UT]');
        ylabel('Range [m]');
        
        subplot(3,1,2);
        gradients = gradients_dailycorrection;
        correction_type = 'daily correction';
        list_plots = [];
        list_legends = {};
        hold on;
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        list_legends{end+1} = 'clouds';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'or');
        list_legends{end+1} = '1. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
        list_legends{end+1} = '2. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
        list_legends{end+1} = '3. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+k');
        list_legends{end+1} = 'min. gradient';
        legend(list_plots,list_legends,'location','eastoutside');
        ylim([0 ymax]);
        xlim([datenum(date,'yyyymmdd') datenum(date,'yyyymmdd')+1]);
        set(gca,'xtick',datenum(date,'yyyymmdd'):4/24:datenum(date,'yyyymmdd')+1);
        datetick('x','HH:MM','keepticks','keeplimits');
        title(['pay ' date ' - ' correction_type]);
        grid on;box on;
        xlabel('Time [UT]');
        ylabel('Range [m]');
        
        subplot(3,1,3);
        gradients = gradients_corr;
        correction_type = 'model';
        list_plots = [];
        list_legends = {};
        hold on;
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        list_legends{end+1} = 'clouds';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'or');
        list_legends{end+1} = '1. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
        list_legends{end+1} = '2. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
        list_legends{end+1} = '3. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+k');
        list_legends{end+1} = 'min. gradient';
        legend(list_plots,list_legends,'location','eastoutside');
        ylim([0 ymax]);
        xlim([datenum(date,'yyyymmdd') datenum(date,'yyyymmdd')+1]);
        set(gca,'xtick',datenum(date,'yyyymmdd'):4/24:datenum(date,'yyyymmdd')+1);
        datetick('x','HH:MM','keepticks','keeplimits');
        title(['pay ' date ' - ' correction_type]);
        grid on;box on;
        xlabel('Time [UT]');
        ylabel('Range [m]');
        
        %% Plot PR2 and gradient witth PBH detection
        
        RCS_list = {RCS_raw,RCS_dailycorrection,RCS_corr};
        grad_list = {grad_raw,grad_dailycorrection,grad_corr};
        gradients_list = {gradients_raw,gradients_dailycorrection,gradients_corr};
        correction_type_list = {'raw','daily correction','model correction'};
        
        for k=1:length(RCS_list)
            %%
            RCS = RCS_list{k};
            grad = grad_list{k};
            gradients = gradients_list{k};
            title_str = [info.chm '/' info.tub '/' station ' - ' date ' - ' correction_type_list{k}];
            
            
            offset = datenum(2000,1,1)-1;
            
            % Choose pcolor parameters
            dxticks = 4;%in hours
            fz = 14;% FontSize
            ylims = [0 2500];
            yticks = ylims(1):250:ylims(2);
            clims = [4.5 5.5];
            clims_grad = [-1000 1000];
            
            % load('ypcmap2','cmap');
            
            figure('Units','normalized','Position',[0.005 0.05 0.8 0.6],'Name',[ 'Pcolor ' correction_type_list{k}])
            handles.axes_pcolor_RCS = subplot(2,1,1);
            handles.surface_pcolor_RCS = pcolor(chm.time-offset,chm.range,log10(abs(RCS)));
            xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset];
            xticks = xlims(1)-mod(mod(xlims(1),1)*24,dxticks)/24:dxticks/24:xlims(2);
            set(gca,'XLim',xlims,'XTick',xticks,'XMinorTick','on');
            set(gca,'YLim',ylims,'YTick',yticks,'YMinorTick','on');
            grid on;
            set(gca,'Layer','top');
            datetick('x','HH:MM','keepticks','keeplimits');
            shading flat
            caxis(clims)
            %     colormap(cmap)
            hc = colorbar;
            xlabel('Time UT [h]','FontSize',fz);
            ylabel(hc,'log10(abs(beta\_raw))','FontSize',fz);
            ylabel('Range (m)','FontSize',fz);
            % daspect(aspect_ratio);
            list_plots = [];
            list_legends = {};
            
            if(disp_pbl)
                hold on;
                list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(4,:),'or','MarkerSize',8,'markerfacecolor','r');
                list_legends{end+1} = '1. strongest gradient';
                list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(3,:),'>g','MarkerSize',7,'markerfacecolor','g');
                list_legends{end+1} = '2. strongest gradient';
                list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(2,:),'<b','MarkerSize',7,'markerfacecolor','b');
                list_legends{end+1} = '3. strongest gradient';
                list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(1,:),'.k','MarkerSize',7,'markerfacecolor','k');
                list_legends{end+1} = 'min. gradient';
            end
            
            if any(time==datenum(date,'yyyymmdd'))
                % add area where the daily overlap function was calculated
                unique_time_start=unique(time_start);
                unique_time_end=unique(time_end);
                
                hold on
                range_start_min=NaN(size(unique_time_start));
                range_start_max=NaN(size(unique_time_start));
                
                % find areas where it the overlap function was calculated
                for i=1:length(unique_time_start)
                    index=time_start==unique_time_start(i);
                    if any(index)
                        range_start_min(i)=min(range_start(index));
                        range_start_max(i)=max(range_end(index));
                    end
                end
                
                % plot these areas without overlaping lines
                range_start_min_all=NaN(size(chm.time));
                range_start_max_all=NaN(size(chm.time));
                for i=1:length(chm.time)
                    index=chm.time(i)>=unique_time_start & chm.time(i)<unique_time_end;
                    if any(index)
                        range_start_min_all(i)=min(range_start_min(index));
                        range_start_max_all(i)=max(range_start_max(index));
                    end
                end
                plot(chm.time-offset,range_start_min_all,'--k')
                list_plots(end+1)=plot(chm.time-offset,range_start_max_all,'--k');
                list_legends{end+1} = 'Ref. zone';
                
                %Vertical bar if there is a new box
                for i=2:length(chm.time)
                    if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i-1))
                        plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
                    end
                end
                
                %Vertical bar if a box is ending
                for i=1:length(chm.time)-1
                    if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i+1))
                        plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
                    end
                end
            end
            % Plot Clouds
            if(disp_clouds)
                for j=1:3
                    cbh = chm.cbh(j,:);
                    cbh(cbh <= chm.cho) = NaN;
                                        hcbh = line(chm.time-offset,cbh,'LineStyle','none','Marker','.','MarkerSize',7,'Color',[0.5 0.5 0.5]);

%                     hcbh = line(chm.time-offset,cbh-chm.cho,'LineStyle','none','Marker','.','MarkerSize',7,'Color',[0.5 0.5 0.5]);
                end
                list_plots(end+1) = hcbh;
                list_legends{end+1} = 'CBH';
            end
            if ~isempty(list_plots)
                legend(list_plots,list_legends,'Orientation','horizontal');
            end
            
            % cosmetics
            set(gca,'FontSize',fz);
            title(title_str,'FontSize',fz);
            handles.axes_pcolor_grad = subplot(2,1,2);
            handles.surface_pcolor_grad = pcolor(chm.time-offset,chm.range,grad);
            xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset];
            xticks = xlims(1)-mod(mod(xlims(1),1)*24,dxticks)/24:dxticks/24:xlims(2);
            set(gca,'XLim',xlims,'XTick',xticks,'XMinorTick','on');
            set(gca,'YLim',ylims,'YTick',yticks,'YMinorTick','on');
            grid on;
            set(gca,'Layer','top');
            datetick('x','HH:MM','keepticks','keeplimits');
            shading flat
            caxis(clims_grad)
            colormap(jet)
            hc = colorbar;
            xlabel('Time UT [h]','FontSize',fz);
            ylabel(hc,'grad(beta\_raw)','FontSize',fz);
            ylabel('Range (m)','FontSize',fz);
            set(gca,'FontSize',fz);
            set(gcf,'PaperPositionMode','auto');
            
            if info.save_plots==1
                disp('Saving figures')
                saveas(gcf,[folder_output  'PR2_grad_' station suffix '_' date '_' corrections_to_analyze '_' correction_type_list{k} '.png' ])
                saveas(gcf,[folder_output  'PR2_grad_' station suffix '_' corrections_to_analyze '_' correction_type_list{k} '.fig' ])
            end
        end
    end
    if  length(time_vec_test) >2
        close all
    end
end