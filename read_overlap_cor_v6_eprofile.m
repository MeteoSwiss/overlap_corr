%Read Overlap Corrections for paper
% Calculate a model according to the internal temperature.
% V4: Work for KSE, PAY and SIRTA
% Remove plots not in paper
clear variables;clc;close all;

set(0,'DefaultFigureVisible','on')

%% Main INPUT
% station='0-20000-0-10393';
% station='0-20000-0-16054';
%station = '0-20008-0-UGR';
% station = '0-20000-0-03808';
% station = '0-20008-0-INO';
%station = '0-20000-0-06348';
%station = '0-20000-0-06610';
% station = '0-20008-0-BRN';
  station = '0-20000-0-06784';
%% INPUTS for each station
switch station
     case '0-20000-0-06784'
        %% Davos    REDO A CONFIG WITH MORE MONTHS IN THE FUTURE
        info.start_day  = 1;
        info.start_month= 7;
        info.start_year = 2021;
        
        info.end_day  =  19;
        info.end_month=  10;
        info.end_year =  2021;
        % Time range to apply the correction
        info_test.start_day  = 9;
        info_test.start_month= 8;
        info_test.start_year = 2021;
        
        info_test.end_day  =  13;
        info_test.end_month=  8;
        info_test.end_year =  2021;
        id = 'A';
        info.tub='TUB120011';
        info.chm = 'TUB120011';
        list_dates_badquality = {};
    case '0-20008-0-BRN'
        %% Bern
        info.start_day  = 14;
        info.start_month= 6;
        info.start_year = 2018;
        
        info.end_day  =  29;
        info.end_month=  7;
        info.end_year =  2021;
        % Time range to apply the correction
        info_test.start_day  = 28;
        info_test.start_month= 5;
        info_test.start_year = 2021;
        
        info_test.end_day  =  28;
        info_test.end_month=  5;
        info_test.end_year =  2021;
        id = 'A';
        info.tub='TUB150046';
        info.chm = 'TUB150046';
        list_dates_badquality = {'20190102','20190221','20190412','20190515','20190716','20190802','20190814',...
                                '20200407','20200528','20200627','20210422'};

    case '0-20000-0-06610'
        %% Payerne
        info.start_day  = 4;
        info.start_month= 3;
        info.start_year = 2021;
        
        info.end_day  =  26;
        info.end_month=  7;
        info.end_year =  2021;
        % Time range to apply the correction
        info_test.start_day  = 3;
        info_test.start_month= 4;
        info_test.start_year = 2021;
        
        info_test.end_day  =  3;
        info_test.end_month=  4;
        info_test.end_year =  2021;
        id = 'A';
        info.tub='TUB200009';
        info.chm = 'TUB200009';
        list_dates_badquality = {'20210506'};% not better than Lufft

switch station
    case '0-20000-0-10393'
        %% lindenberg
        info.start_day  = 1;
        info.start_month= 1;
        info.start_year = 2018;
        
        info.end_day  =  31;
        info.end_month=  12;
        info.end_year =  2020;
        % Time range to apply the correction
        info_test.start_day  = 12;
        info_test.start_month= 5;
        info_test.start_year = 2019;
        
        info_test.end_day  =  12;
        info_test.end_month=  5;
        info_test.end_year =  2019;
        id = '0';
        info.tub='TUB120001';
        info.chm = 'TUB120001';
        list_dates_badquality = {'20140915';'20150315';'20160409';'20180207';...
            '20150704';'20150705';'20150810';'20170928';'20150519'};% not better than Lufft
    case '0-20000-0-16054'
        %% Aosta
        list_dates_badquality = {};% not better than Lufft
        id = '0';
        info.tub=[station '0'];
        info.chm = '';
        
        info.start_day  = 1;
        info.start_month= 1;
        info.start_year = 2018;
        
        info.end_day  =  31;
        info.end_month=  12;
        info.end_year =  2020;
        % Time range to apply the correction
        info_test.start_day  = 12;
        info_test.start_month= 5;
        info_test.start_year = 2019;
        
        info_test.end_day  =  12;
        info_test.end_month=  5;
        info_test.end_year =  2019;
    case '0-20008-0-UGR'
        %% granada
        id = 'A';
        info.tub=[station 'A'];
        info.chm = '';
        
        info.start_day  = 1;
        info.start_month= 1;
        info.start_year = 2018;
        
        info.end_day  =  31;
        info.end_month=  12;
        info.end_year =  2020;
        % Time range to apply the correction
        info_test.start_day  = 12;
        info_test.start_month= 5;
        info_test.start_year = 2019;
        
        info_test.end_day  =  12;
        info_test.end_month=  5;
        info_test.end_year =  2019;
        list_dates_badquality = {'20180425'};% not better than Lufft
        
    case '0-20000-0-03808'
        %% Camborne
        list_dates_badquality = {'20180630','20190423','20190627','20200420'};% not better than Lufft
        id = 'A';
        info.tub='TUB100005';
        info.chm = 'TUB100005';
        info.start_day  = 1;
        info.start_month= 1;
        info.start_year = 2018;
        
        info.end_day  =  31;
        info.end_month=  12;
        info.end_year =  2020;
        % Time range to apply the correction
        info_test.start_day  = 12;
        info_test.start_month= 5;
        info_test.start_year = 2019;
        
        info_test.end_day  =  12;
        info_test.end_month=  5;
        info_test.end_year =  2019;
    case '0-20008-0-INO'
        %% MAGURELE
        id = 'A';
        info.tub='TUB170068';
        info.chm = 'TUB170068';
        info.start_day  = 1;
        info.start_month= 9;
        info.start_year = 2019;
        
        info.end_day  =  31;
        info.end_month=  12;
        info.end_year =  2020;
        % Time range to apply the correction
        info_test.start_day  = 21;
        info_test.start_month= 3;
        info_test.start_year = 2020;
        
        info_test.end_day  =  21;
        info_test.end_month=  3;
        info_test.end_year =  2020;
        
        list_dates_badquality = {};% not better than Lufft
        
         case '0-20000-0-06348'
        %% Cabauw
        id = 'A';
        info.tub='TUB140015';
        info.chm = 'TUB140015';
        info.start_day  = 1;
        info.start_month= 6;
        info.start_year = 2016;
        
        info.end_day  =  31;
        info.end_month=  6;
        info.end_year =  2018;
        
        % Time range to apply the correction
        info_test.start_day  = 4;
        info_test.start_month= 6;
        info_test.start_year = 2018;
        
        info_test.end_day  =  4;
        info_test.end_month=  6;
        info_test.end_year =  2018;
        
        list_dates_badquality = {};% not better than Lufft
        
    otherwise
        error('Define statoon info')
        list_dates_badquality = {};% not better than Lufft
        
end


%% other inputs
%File paths
folder_corrections =['C:/WorkMCH/PAY/stratfinder_external/overlap_correction/' station '/'];
folder_output = folder_corrections;
folder_ncdata = 'C:/WorkMCH/PAY/data/CHM15k_L1_E-PROFILES/';

suffix = ''; %'' _after_moving_2017

min_nb_good_samples_after_outliers_removal  = 10;

info.plot=1; % Plot all data?
info.save_plots=1;
disp_clouds = true;
disp_pbl = 0;
info_reloading=1; %(Reload all data?)

create_netcdf = 1;

%folder_results_algo = 'D:\Projects\2018\201806 Start Finder Overlap\';
%folder_ov_ref = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\overlap_functions_Lufft\';
%folder_ncdata = 'G:\E_PROFILE_ALC\L1_FILES\';

corrections_to_analyze = 'good_enough';%{'all','good_enough','well_trusted'};

list_dates_nodetection = {};
list_dates_acceptablequality = {};
list_dates_goodquality = {};
list_dates_verygoodquality = {};
list_dates_outliers_temperature = {};


%% check availability
time_vec=datenum(info.start_year,info.start_month,info.start_day):datenum(info.end_year,info.end_month,info.end_day);
dtime_vec = datetime(time_vec,'convertfrom','datenum');
availability = false(size(time_vec));
for t = 1:length(time_vec)
    str=[folder_ncdata  station datestr(time_vec(t),'/yyyy/mm/')...
        'L1_' station '_' id datestr(time_vec(t),'yyyymmdd') '*.nc'];
    
    list=dir(str);
    if ~isempty(list)
        availability(t) = true;
    end
end
figure
plot(dtime_vec, availability)
ylim([-0.1 1.1])
title([station ': Availability ' num2str(sum(availability)/length(availability)*100) '%'])
if info.save_plots==1
    filename_figure = [folder_output  '1_availability' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
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
    %     voltage_PM_mean0 = [];
    for t=1:length(time_vec)
        try
            file=['result_ovp_fc_' info.chm info.tub '_' datestr(time_vec(t),'yyyymmddHHMMSS') '.mat'];
            
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
                
                file=['L1_' station '_' id  datestr(time_vec(t),'yyyymmdd') '.nc'];
                folder_day=[folder_ncdata station '/' datestr(time_vec(t),'yyyy/mm/') ] ;
                
                
                disp([folder_day file])
                if exist([folder_day file],'file')==0
                    error('Missing NetCDFfile file')
                end
                
                
                
                temp_int = ncread([folder_day file],'temp_int');
                
                
                time_nc = datenum(1970,1,1)+ncread([folder_day file],'time');
                temp_tmp = NaN(length(result.index_final),1);
                range = ncread([folder_day file],'range');
                for i=1:length(result.index_final)
                    indt = time_nc>=result.time_start(result.index_final(i)) & time_nc<=result.time_end(result.index_final(i));
                    temp_tmp(i) = nanmean(temp_int(indt));
                    %                 voltage_tmp(i) = nanmean(nn1(indt));
                    
                end
                
                %         temp_mean(k)=mean(ncread([folder_day file],'temp_int'));
                %         voltage_PM_mean(k)=mean(ncread([folder_day file],'nn1'));
                temp_mean0(k)=nanmedian(temp_tmp);
                %             voltage_PM_mean0(k)=nanmedian(voltage_tmp);
                % %             zenith_angle_mean(k) = nanmedian(zenith_angle);
                
                
                k=k+1;
            else
                disp('No file:')
                disp([folder_day file])
            end
        catch
            warning('pb')
        end
    end
    
    %% Read: overlap correction and scaling in cfg file
    list = dir([folder_ov_ref,info.tub,'_*.cfg']);
    %     if isempty(list)
    %         answer=questdlg({'No reference overlap function found','Do you want to calibrate with TUB120011 overlap function'},...
    %             'No reference overlap function','Yes','No','No');
    %         if strcmp(answer,'Yes')
    warning('No overlap function found. Using TUB120011 overlap function')
    fid = fopen('TUB120011_20121112_1024.cfg');
    ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
    overlap_ref = cell2mat(ov_cell);
    if length(overlap_ref)~=length(range)
        overlap_ref = interp1(0:14.984999:15344,overlap_ref,range);
        disp('Interpolating')
    end
    %         else
    %             error('No overlap function found. Ask the corresponding overlap to the manufacturer')
    %         end
    %     else
    %         file=[folder_ov_ref list(1).name];
    %         fid=fopen(file);
    %         disp(file);
    %         data=textscan(fid,'%f','delimiter','\t','headerlines',1);
    %         overlap_ref=data{1};
    %     end
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

if k==1
    error('No files found')
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
% voltage_PM_mean = voltage_PM_mean0(index_of_interest);


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
ylabel(hc,'Median internal Temperature [°C]');
caxis([15,40]);
colormap(jet);
ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200)
grid on;box on;
xlabel('Range [m]')
ylabel('Overlap function')
title(station)
if info.save_plots==1
    filename_figure = [folder_output  '2_all overlap functions' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
end

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
ylabel(hc,'Median internal Temperature [°C]');
ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-75 75]);
caxis([15 40]);
title(station)

grid on
box on
if info.save_plots==1
    filename_figure = [folder_output  '3 Relative difference' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
end
%% plot rel diff changing with time
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
title(station)

box on
if info.save_plots==1
    filename_figure = [folder_output  '4_rel diff changing with time_' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
end
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
xlabel('Median internal Temperature [°C]','FontSize',fz);
ylabel({'L^2 Norm of the Relative Difference';'between corrected and uncorrected signal' },'FontSize',fz);
ylim([0 2]);
grid on;
legend(h_fit,['R^2=' num2str(R2) ', RMSE=' num2str(RMSE)]);
set(gca,'FontSize',fz);
title(station)

if info.save_plots==1
    filename_figure = [folder_output  '5_Scatter T_' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
end
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
    
    fo = polyfit(x_fit,y_fit,1);
    a(j) = fo(1);
    b(j) = fo(2);
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
ylabel(hc,'Internal Temperature [°C]');

ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected (with model) and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-75 75]);
grid on
box on
if info.save_plots==1
    filename_figure = [folder_output  '6b_relative difference model_' station suffix  '.png' ];
    disp(filename_figure)
    saveas(gcf,filename_figure)
end

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
    date_str = datestr(time_vec_test(t),'yyyymmdd');
    end_str = datestr(time_vec_test(t)+0.99999999,'yyyymmddHHMMSS');
    
    if ~any(time==datenum(date_str,'yyyymmdd'))
        warning('No Overlap correction for the selected day')
    end
    
    %     [chm,chminfo] = get_chm15k_from_files('pay',[date,'000000'],datestr(datenum(date,'yyyymmdd')+0.99999,'yyyymmddHHMMSS'),folder_ncdata);
    %     [chm,chminfo]=readcorrectlyncfile3(station_str,date,folder_ncdata);
    L1 =read_L1_EPROFILE_v4(station,id,date_str,end_str,folder_ncdata);
    
    if isempty(L1)
        warning('no file')
        continue
    end
    
    RCS_raw = L1.rcs_0';
    if any(time==datenum(date_str,'yyyymmdd'))
        RCS_dailycorrection = L1.rcs_0'./repmat(overlap_cor(time==datenum(date_str,'yyyymmdd'),:)',1,size(L1.rcs_0,1)).*repmat(overlap_ref,1,size(L1.rcs_0,1));
    else
        RCS_dailycorrection=NaN(size(    RCS_raw));
    end
    %% Calc Grad
    grad_raw = NaN(size(RCS_raw));
    for j=1:length(L1.time)
        grad_raw(:,j) = 1./(2*L1.range).*conv(RCS_raw(:,j),[1 0 -1],'same');
    end
    
    grad_dailycorrection = NaN(size(RCS_dailycorrection));
    if any(time==datenum(date_str,'yyyymmdd'))
        for j=1:length(L1.time)
            grad_dailycorrection(:,j) = 1./(2*L1.range).*conv(RCS_dailycorrection(:,j),[1 0 -1],'same');
        end
    end
    
    %% calculate RCS and grad(RCS) with the temperature model
    ov_rec_all = ones(size(RCS_raw));
    for i=1:find(L1.range<=1200,1,'last')
        rel_diff = polyval([a_filtered(i),b_filtered(i)],L1.temp_int-273.15);
        ov_rec_all(i,:) = 1./ (rel_diff/100/overlap_ref(i) + 1./overlap_ref(i));
    end
    %     ov_rec_all(find(chm.range<=1200,1,'last')+1:end,:) = repmat(overlap_ref(find(chm.range<=1200,1,'last')+1:end),1,size(RCS_raw,2));
    
    RCS_corr = RCS_raw./ov_rec_all.*repmat(overlap_ref,1,size(RCS_raw,2));
    
    % Calc Grad
    grad_corr = NaN(size(RCS_corr));
    for j=1:length(L1.time)
        grad_corr(:,j) = 1./(2*L1.range).*conv(RCS_corr(:,j),[1 0 -1],'same');
    end
    
    %% Calculate PBL with simple code
    %     [gradients_raw,~] = simplePBLdetection(L1,overlap_ref,overlap_ref);
    %     if any(time==datenum(date_str,'yyyymmdd'))
    %         [gradients_dailycorrection, ~] = simplePBLdetection(L1,overlap_ref,overlap_cor(time==datenum(date_str,'yyyymmdd'),:)');
    %     else
    %         gradients_dailycorrection=NaN(size(gradients_raw));
    %     end
    %     [gradients_corr,xtime] = simplePBLdetection(L1,overlap_ref,ov_rec_all);
    %
    %% Write text file with PBL heigt
    %     disp('Writing text file with PBL height')
    %     M=[datevec(datenum(date_str,'yyyymmdd')+xtime) gradients_raw' gradients_dailycorrection' gradients_corr'];
    %     dlmwrite([folder_output '/PBL_' station '_' date_str '_' corrections_to_analyze  '.csv' ],M)
    
    if info.plot==1
        %% plot internal temperature
        figure;
        plot(L1.time,L1.temp_int-273.15);
        hold on
        plot(L1.time,L1.temp_ext-273.15);
        datetick;
        box on;grid on;
        xlabel('Time [UT]');ylabel('Temperature [°C]');title(date_str);
        legend('Internal','External')
        
        if any(time==datenum(date_str,'yyyymmdd'))
            
            %% load and plot daily overlap correction
            file=['result_ovp_fc_' info.chm info.tub '_' date_str '.mat'];
            folder_day=[folder_corrections datestr(datenum(date_str,'yyyymmdd'),'yyyy/mm/')];
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
            try
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
                title(['overlap functions ' date_str])
            catch
                warning('pb')
            end
            %% plot 10: Fit Example
            try
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
                
                dif_all=NaN(length(result.time_start(result.index_final)),length(L1.range));
                lrcs_all=NaN(length(result.time_start(result.index_final)),length(L1.range));
                for i=1:length(result.time_start(result.index_final))
                    index_time=L1.time>=result.time_start(result.index_final(i)) & L1.time<result.time_end(result.index_final(i));
                    index_range =L1.range>=result.range_start(result.index_final(i)) & L1.range<result.range_end(result.index_final(i));
                    p = polyfit(L1.range(index_range),mean(log10(abs(L1.beta_raw(index_range,index_time))),2),1);
                    dif = mean(log10(abs(L1.beta_raw(:,index_time))),2)-polyval(p,L1.range);
                    dif( L1.range>=result.range_end(result.index_final(i)))=0;
                    dif_all(i,:)=dif;
                    lrcs_all(i,:)=log10(abs(mean(L1.beta_raw(:,index_time),2)));
                    overlap_corr_factor = 10.^(dif);
                    subplot(1,2,1)
                    %                 plot(log10(abs(mean(chm.beta_raw(:,index_time),2))),chm.range,'.k')
                    h=plot([p(2) p(1)*1000+p(2)],[0 1000],'--r');
                    legend(h,'Linear fit')
                    
                    %                 subplot(1,2,2)
                    %                 plot(1+dif,chm.range,'.k')
                end
                subplot(1,2,1)
                errorshade4(median(lrcs_all),L1.range',...
                    nanmin(lrcs_all),nanmax(lrcs_all),...
                    zeros(size(L1.range')), zeros(size(L1.range')),'k');
                plot(median(lrcs_all),L1.range,'color','k','linewidth',2);
                
                subplot(1,2,2)
                errorshade4(median(dif_all)+1,L1.range',...
                    nanmin(dif_all+1),nanmax(dif_all+1),...
                    zeros(size(L1.range')), zeros(size(L1.range')),'k');
                plot(median(dif_all)+1,L1.range,'color','k','linewidth',2);
                
                
                suptitle([station_name ': ' date_str])
            catch
                warning('pb')
            end
        else
            warning('No daily correction file')
        end
        
        
        %% Plot PBL
        %         figure;
        %         ymax=2500;
        %         subplot(3,1,1);
        %         gradients = gradients_raw;
        %         correction_type = 'raw';
        %         list_plots = [];
        %         list_legends = {};
        %         hold on;
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        %         list_legends{end+1} = 'clouds';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(4,:),'or','markerfacecolor','r','markersize',10);
        %         list_legends{end+1} = '1. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(3,:),'>g');
        %         list_legends{end+1} = '2. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(2,:),'<b');
        %         list_legends{end+1} = '3. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(1,:),'+k');
        %         list_legends{end+1} = 'min. gradient';
        %         legend(list_plots,list_legends,'location','eastoutside');
        %         ylim([0 ymax]);
        %         xlim([datenum(date_str,'yyyymmdd') datenum(date_str,'yyyymmdd')+1]);
        %         set(gca,'xtick',datenum(date_str,'yyyymmdd'):4/24:datenum(date_str,'yyyymmdd')+1);
        %         datetick('x','HH:MM','keepticks','keeplimits');
        %         title(['pay ' date_str ' - ' correction_type]);
        %         grid on;box on;
        %         xlabel('Time [UT]');
        %         ylabel('Range [m]');
        %
        %         subplot(3,1,2);
        %         gradients = gradients_dailycorrection;
        %         correction_type = 'daily correction';
        %         list_plots = [];
        %         list_legends = {};
        %         hold on;
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        %         list_legends{end+1} = 'clouds';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(4,:),'or');
        %         list_legends{end+1} = '1. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(3,:),'>g');
        %         list_legends{end+1} = '2. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(2,:),'<b');
        %         list_legends{end+1} = '3. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(1,:),'+k');
        %         list_legends{end+1} = 'min. gradient';
        %         legend(list_plots,list_legends,'location','eastoutside');
        %         ylim([0 ymax]);
        %         xlim([datenum(date_str,'yyyymmdd') datenum(date_str,'yyyymmdd')+1]);
        %         set(gca,'xtick',datenum(date_str,'yyyymmdd'):4/24:datenum(date_str,'yyyymmdd')+1);
        %         datetick('x','HH:MM','keepticks','keeplimits');
        %         title(['pay ' date_str ' - ' correction_type]);
        %         grid on;box on;
        %         xlabel('Time [UT]');
        %         ylabel('Range [m]');
        %
        %         subplot(3,1,3);
        %         gradients = gradients_corr;
        %         correction_type = 'model';
        %         list_plots = [];
        %         list_legends = {};
        %         hold on;
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
        %         list_legends{end+1} = 'clouds';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(4,:),'or');
        %         list_legends{end+1} = '1. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(3,:),'>g');
        %         list_legends{end+1} = '2. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(2,:),'<b');
        %         list_legends{end+1} = '3. strongest gradient';
        %         list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime,gradients(1,:),'+k');
        %         list_legends{end+1} = 'min. gradient';
        %         legend(list_plots,list_legends,'location','eastoutside');
        %         ylim([0 ymax]);
        %         xlim([datenum(date_str,'yyyymmdd') datenum(date_str,'yyyymmdd')+1]);
        %         set(gca,'xtick',datenum(date_str,'yyyymmdd'):4/24:datenum(date_str,'yyyymmdd')+1);
        %         datetick('x','HH:MM','keepticks','keeplimits');
        %         title(['pay ' date_str ' - ' correction_type]);
        %         grid on;box on;
        %         xlabel('Time [UT]');
        %         ylabel('Range [m]');
        
        %% Plot PR2 and gradient with PBH detection
        
        RCS_list = {RCS_raw,RCS_dailycorrection,RCS_corr};
        grad_list = {grad_raw,grad_dailycorrection,grad_corr};
        %         gradients_list = {gradients_raw,gradients_dailycorrection,gradients_corr};
        correction_type_list = {'raw','daily correction','model correction'};
        
        for k=1:length(RCS_list)
            %%
            RCS = RCS_list{k};
            grad = grad_list{k};
            %             gradients = gradients_list{k};
            title_str = [info.tub '/' station ' - ' date_str ' - ' correction_type_list{k}];
            
            
            offset = datenum(2000,1,1)-1;
            
            % Choose pcolor parameters
            dxticks = 4;%in hours
            fz = 14;% FontSize
            ylims = [0 2500];
            yticks = ylims(1):250:ylims(2);
            clims = [4.5 5.5];
            clims_grad = [-40 40];
            
            % load('ypcmap2','cmap');
            
            figure('Units','normalized','Position',[0.005 0.05 0.8 0.6],'Name',[ 'Pcolor ' correction_type_list{k}])
            handles.axes_pcolor_RCS = subplot(2,1,1);
            handles.surface_pcolor_RCS = pcolor(L1.time-offset,L1.range(L1.range<ylims(2)),log10(abs(RCS(L1.range<ylims(2),:))));
            xlims = [floor(L1.time(1))-offset, ceil(L1.time(end))-offset];
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
            
            %             if(disp_pbl)
            %                 hold on;
            %                 list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime-offset,gradients(4,:),'or','MarkerSize',8,'markerfacecolor','r');
            %                 list_legends{end+1} = '1. strongest gradient';
            %                 list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime-offset,gradients(3,:),'>g','MarkerSize',7,'markerfacecolor','g');
            %                 list_legends{end+1} = '2. strongest gradient';
            %                 list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime-offset,gradients(2,:),'<b','MarkerSize',7,'markerfacecolor','b');
            %                 list_legends{end+1} = '3. strongest gradient';
            %                 list_plots(end+1) = plot(datenum(date_str,'yyyymmdd')+xtime-offset,gradients(1,:),'.k','MarkerSize',7,'markerfacecolor','k');
            %                 list_legends{end+1} = 'min. gradient';
            %             end
            
            if any(time==datenum(date_str,'yyyymmdd'))
                try
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
                    range_start_min_all=NaN(size(L1.time));
                    range_start_max_all=NaN(size(L1.time));
                    for i=1:length(L1.time)
                        index=L1.time(i)>=unique_time_start & L1.time(i)<unique_time_end;
                        if any(index)
                            range_start_min_all(i)=min(range_start_min(index));
                            range_start_max_all(i)=max(range_start_max(index));
                        end
                    end
                    plot(L1.time-offset,range_start_min_all,'--k')
                    list_plots(end+1)=plot(L1.time-offset,range_start_max_all,'--k');
                    list_legends{end+1} = 'Ref. zone';
                    
                    %Vertical bar if there is a new box
                    for i=2:length(L1.time)
                        if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i-1))
                            plot([L1.time(i)  L1.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
                        end
                    end
                    
                    %Vertical bar if a box is ending
                    for i=1:length(L1.time)-1
                        if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i+1))
                            plot([L1.time(i)  L1.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
                        end
                    end
                catch
                    warning('pb')
                end
            end
            % Plot Clouds
            if(disp_clouds)
                for j=1:3
                    cbh = L1.cloud_base_height(:,j);
                    %                     cbh(cbh <= L1.cho) = NaN;
                    hcbh = line(L1.time-offset,cbh,'LineStyle','none','Marker','.','MarkerSize',7,'Color',[0.5 0.5 0.5]);
                    
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
            handles.surface_pcolor_grad = pcolor(L1.time-offset,L1.range(L1.range<ylims(2)),grad(L1.range<ylims(2),:));
            xlims = [floor(L1.time(1))-offset, ceil(L1.time(end))-offset];
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
                saveas(gcf,[folder_output  'PR2_grad_' station suffix '_' date_str '_' corrections_to_analyze '_' correction_type_list{k} '.png' ])
                %                 saveas(gcf,[folder_output  'PR2_grad_' station suffix '_' corrections_to_analyze '_' correction_type_list{k} '.fig' ])
            end
        end
    end
    if  length(time_vec_test) >2
        close all
    end
end