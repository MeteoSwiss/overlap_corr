%Read Overlap Corrections
clear variables;clc;close all;


%% INPUTS
info.start_day  = 1;
info.start_month= 2;
info.start_year = 2013;

info.end_day  =  25;
info.end_month=  11;
info.end_year =  2014;

info.chm='CHM120106';
info.tub='TUB120011';

station='pay';

use_local_data = 0;%{0,1};
corrections_to_analyze = 'all';%{'all','good_enough','well_trusted'};

if use_local_data
    folder_corrections='C:\AllData\SharedData_Maxime\py\';
    folder_ncdata='C:\AllData\';
    folder_ov_ref='M:\pay-home\pay\users\poy\My Documents\workYP\lib_overlap\';
    folder_results_algo = 'C:\AllData\SharedData_Maxime\py\';
else
    folder_corrections = '\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\overlap_correction\corrections\';
    folder_ncdata = '\\meteoswiss.ch\mch\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
    folder_ov_ref = '\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\overlap_correction\overlap_functions_Lufft\';
    folder_results_algo = '\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\overlap_correction\';
end

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

min_nb_good_samples_after_outliers_removal  = 10;


%% read all files with a final overlap correction output
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
        if any(strcmp(list_dates_badquality,datestr(time_vec(t),'yyyymmdd')))
        	quality0(k) = 0;
        end
        if any(strcmp(list_dates_acceptablequality,datestr(time_vec(t),'yyyymmdd')))
            quality0(k) = 1;
        end
        if any(strcmp(list_dates_goodquality,datestr(time_vec(t),'yyyymmdd')))
            quality0(k) = 2;
        end
        if any(strcmp(list_dates_verygoodquality,datestr(time_vec(t),'yyyymmdd')))
            quality0(k) = 3;
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
        
        file=[datestr(time_vec(t),'yyyymmdd') '_' station_str '_' info.chm '_000.nc'];
        folder_day=[folder_ncdata datestr(time_vec(t),'yyyy/mm/') ] ;
        
        disp([folder_day file])
        
        %convert scale factor
        scale_factor = ncreadatt([folder_day file],'temp_int','scale_factor');
        if(~isnumeric(scale_factor))
            ncwriteatt([folder_day file],'temp_int','scale_factor',1./str2num(scale_factor));
        end
        
        temp_int = ncread([folder_day file],'temp_int');
        nn1 = ncread([folder_day file],'nn1');
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
        
        k=k+1;
        
    end
end

%% load overlap correction and scaling in cfg file
list = dir([folder_ov_ref,info.tub,'_*.cfg']);
if isempty(list)
    error('missing overlap function');
else
    file=[folder_ov_ref list(1).name];
    fid=fopen(file);
    disp(file);
    data=textscan(fid,'%f','delimiter','\t','headerlines',1);
    overlap_ref=data{1};
end

save all_correction.mat

%% select set of overlap functions to analyse

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

% % disp which are very good but T outliers
% for i=1:length(time0)
%     if any(strcmp(list_dates_verygoodquality,datestr(time0(i),'yyyymmdd'))) && any(strcmp(list_dates_outliers_temperature,datestr(time0(i),'yyyymmdd')))
%         disp(datestr(time0(i),'yyyymmdd'));
%         disp(temp_mean(i)-273.15);
%     end
% end


%% plot relative difference between corrected and uncorrected signal
figure

T = temp_mean-273.15;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

% set(gcf,'defaultAxesColorOrder',RGB)

hold on

relative_difference=NaN(length(time),length(range));
for i=1:length(time)
    relative_difference(i,:)=(-1./overlap_ref'+1./overlap_cor(i,:))./(1./overlap_ref')*100;
    h=scatter(relative_difference(i,:),range,...
        ones(size(overlap_cor(i,:)))*12,...
        repmat(temp_mean(i)-273.15,size(overlap_cor(i,:))),...
        'filled','displayname',datestr(time(i)));
    hp = plot(relative_difference(i,:),range,'displayname',datestr(time(i)));
    
end
hc = colorbar;
% ylabel(hc,'mean interior T during the day [°C]');
ylabel(hc,'median interior T during calculation of overlap correction [°C]');

ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-80 20]);
caxis([15 40]);
grid on
box on

%% plot all overlap functions

figure;

T = temp_mean-273.15;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

% set(gcf,'DefaultAxesColorOrder',RGB)

hold on;
for i=1:length(time)
%     plot(range,overlap_cor(i,:),'color',0.5*ones(1,3),'displayname',datestr(time(i)));
    plot(range,overlap_cor(i,:),'displayname',datestr(time(i)));
end
hc = colorbar;
% ylabel(hc,'mean interior T during the day [°C]');
ylabel(hc,'median interior T during calculation of overlap correction [°C]');
caxis([15,40]);
colormap(jet);
plot(range,overlap_ref,'--k','LineWidth',2)
patcher(range',nanmin(overlap_cor,[],1),nanmax(overlap_cor,[],1),0.25*ones(1,3),[],'FaceAlpha',0.25);

ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200)
grid on;box on;
xlabel('Range [m]')
title('All corrected overlap functions')

%% Scatter Plot: interior T vs. Diff w.r.t. to Lufft

ti = time;
T = temp_mean;
ov_cor = overlap_cor;
ov = overlap_ref;
disp_text = false;

A = T-273.15;
%l2-norm
B = 100*sqrt(trapz(abs(ov_cor'-repmat(ov,1,length(ti))).^2)) ./ sqrt(trapz(abs(repmat(ov,1,length(ti))).^2));
%linf-norm
% B = 100*max(abs(ov_cor'-repmat(ov,1,length(ti))),[],1) ./ max(abs(repmat(ov,1,length(ti))),[],1);


p = polyfit(A,B,1);

SStot = sum((B'-repmat(mean(B),length(B),1)).^2);
SSres = sum((B'-polyval(p,A')).^2);
R2 = 1-SSres/SStot;
RMSE = sqrt(SSres/length(A));
RMSErel = 100*RMSE / sqrt(sum(polyval(p,A').^2)/length(A));
% disp(['Rsquared of the fit: ' num2str(R2)])
% disp(['RMSE of the fit: ' num2str(RMSE)])

fz = 18;
hf=figure('Units','Normalized','Position',[0 0 1 1],'Visible','on');
scatter(A,B,'filled');% scatter(A,B,36,jet(length(A)),'filled')
if disp_text
    hold on;
    for j=1:length(ti)
        text(A(j),B(j),datestr(ti(j),'dd.mm.yyyy HH:MM'),'HorizontalAlignment','center','FontSize',8);
    end
end
xlim([15 40]);
axis equal;
h_fit = line(xlim,polyval(p,xlim),'Color','k','LineWidth',2,'LineStyle','--');
% h_out = line(xlim,polyval([p(1) p(2)+0.35],xlim),'Color','k','LineStyle','--');
% h_out = line(xlim,polyval([p(1) p(2)-0.35],xlim),'Color','k','LineStyle','--');
% h_out = line(xlim,polyval([p(1) p(2)+3.5],xlim),'Color','k','LineStyle','--');
% h_out = line(xlim,polyval([p(1) p(2)-3.5],xlim),'Color','k','LineStyle','--');

% xlabel('mean interior T during the day [°C]','FontSize',fz);
xlabel('median interior T during calculation of overlap correction [°C]','FontSize',fz);

ylabel('|ov_{corrected} - ov_{manufacturer}|_l^2 / |ov_{manufacturer}|_l^2 (%)','FontSize',fz);
% ylabel('|ov_{corrected} - ov_{manufacturer}|_{l^{\infty}} / |ov_{manufacturer}|_{l^{\infty}} (%)','FontSize',fz);
ylim([0 5]);
% ylim([0 25]);


grid on;

% legend(h_fit,['poly1+-0.35, R^2=' num2str(R2) ', rel. RMSE=' num2str(RMSErel) '%']);
% legend(h_fit,['poly1+-3.5, R^2=' num2str(R2) ', rel. RMSE=' num2str(RMSErel) '%']);
legend(h_fit,['R^2=' num2str(R2) ', RMSE=' num2str(RMSE)]);

set(gca,'FontSize',fz);
set(gcf,'PaperPositionMode','auto');

%% Plots function of time 3D
ov = overlap_ref;
[X,Y] = meshgrid(time-datenum(2000,1,1),double(range));

figure;
surf(X,Y,overlap_cor'-repmat(ov,1,length(time)),repmat(temp_mean-273.15,length(range),1));
% surf(X,Y,overlap_cor',repmat(temp_mean-273.15,length(x),1));
hold on

Ov_350=interp2(X,Y,overlap_cor'-repmat(ov,1,length(time)),time-datenum(2000,1,1),350*ones(size(time)));
OV_350_interped=interp1(time,Ov_350,floor(min(time)):max(time));
OV_350_smoothed=smooth(floor(min(time)):max(time),OV_350_interped,90);
plot3((floor(min(time)):max(time))-datenum(2000,1,1),350*ones(size(floor(min(time)):max(time))),OV_350_smoothed,'r','linewidth',2)

hc = colorbar;
% ylabel(hc,'mean interior T during the day [°C]');
ylabel(hc,'median interior T during calculation of overlap correction [°C]');
colormap(jet);
caxis([15 40])
ylim([0 1200])
set(gca,'ytick',0:200:1200);
box on;grid on;
xlabel('Day')
ylabel('Range [m]')
zlabel('ov\_corrected - ov\_manufacturer');
datetick('x','mm/yy')



%% plot at 4 altitudes
altitudes = [250,350,450,550,650];

figure;
list_plots = [];
list_legends = {};
hold on;
colors = hsv(length(altitudes));
for j=1:length(altitudes)
    [~,ind] = nanmin(abs(range-altitudes(j)));
%     list_plots(end+1) = plot(time,overlap_cor(:,ind),'.-');
    y = (-1./overlap_ref(ind)'+1./overlap_cor(:,ind))./(1./overlap_ref(ind)')*100;
    list_plots(end+1) = plot(time,y,'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:));
    list_legends{end+1} = ['r=' num2str(range(ind),'%2.0f') 'm'];
    hold on;
    plot(time,smooth(y,11),'-','Color',colors(j,:),'linewidth',2);
end
legend(list_plots,list_legends);
set(gca,'xtick',datenum(2013,1:24,1,00,00,00));
set(gca,'xlim',[datenum(2013,4,1),datenum(2014,10,1)]);
box on;grid on;datetick('x','mm/yy','keepticks','keeplimits');
ylim([-40 40]);
ylabel({'Relative difference between';'corrected and uncorrected signal [%]' })
xlabel('month/year');


%% Fit Relative difference according to temperature
a = zeros(size(range));
b = zeros(size(range));
c = zeros(size(range));
r2 = zeros(size(range));
rmse = zeros(size(range));

% figure('nextplot','replacechildren');
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

%         [fo,gof] = fit(x_fit,y_fit,'poly2','Robust','Bisquare');
%         a(j) = fo.p1(1);
%         b(j) = fo.p2(1);
%         c(j) = fo.p3(1);

        [fo,gof] = fit(x_fit,y_fit,'poly1');
%         [fo,gof] = fit(x_fit,y_fit,'poly1','Robust','Bisquare');
%         [fo,gof] = fit(x_fit,y_fit,'poly1','Robust','Bisquare','Weights',weights_fit);
        a(j) = 0;
        b(j) = fo.p1(1);
        c(j) = fo.p2(1);
        r2(j) = gof.rsquare;
        rmse(j) = gof.rmse;
        
%         scatter(x_fit,y_fit,'filled');grid on;box on;line([0 40],feval(fo,[0 40]),'color','k','linewidth',2);
%         xlim([0 40]);xlabel('T [°C]');
% %         ylim([-80 10]);
%         ylabel({'Relative difference between';'corrected and uncorrected signal [%]'})
%         title(['r=' num2str(range(j),'%2.0f') 'm']);
%         pause;
end

temp_vector = 15:0.1:40;
color_vector=jet(length(temp_vector));

figure;

T = temp_vector;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

set(gcf,'DefaultAxesColorOrder',RGB)

hold on
ov_reconstructed=ones(length(temp_vector),length(range));
for i=1:length(temp_vector)
    temp_of_interest=temp_vector(i);
    for j=1:length(range)
        rel_diff = polyval([a(j),b(j),c(j)],temp_of_interest);
        ov_reconstructed(i,j) = 1./ (rel_diff/100/overlap_ref(j) + 1./overlap_ref(j));
    end
    plot(range,ov_reconstructed(i,:),'color',color_vector(i,:));
end

hc = colorbar;
ylabel(hc,'interior T [°C]');
colormap('jet');
caxis([15 40])

h_ref = plot(range,overlap_ref,'--k','LineWidth',2);
ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200)
grid on;box on;
xlabel('Range [m]')
ylabel('Overlap function');
title(['Temperature model, T=' num2str(temp_vector(1)) '...' num2str(temp_vector(end)) '[°C]']);

legend(h_ref,[info.chm info.tub]);

%% plot relative difference model
temp_vector = 15:0.1:40;
color_vector=jet(length(temp_vector));
figure
hold on
ov_reconstructed=ones(length(temp_vector),length(range));
for i=1:length(temp_vector)
    rel_diff = NaN(length(range),1);
    for j=1:length(range)
        rel_diff(j) = polyval([a(j),b(j),c(j)],temp_vector(i));
    end
    
    rel_diff(range<=range(find(b>0 & range <range(find(overlap_ref<0.01,1,'last')),1,'last'))) = NaN;
    
    h=scatter(rel_diff,range,...
        ones(length(range),1)*0,repmat(temp_vector(i),length(range),1),...
        'filled','displayname',num2str(temp_vector(i)));
    h=plot(rel_diff,range,'color',color_vector(i,:));
    
end
hc = colorbar;
ylabel(hc,'interior T [°C]');

ylim([0,1200])
set(gca,'ytick',0:200:1200)
xlabel({'Relative Difference between';'corrected (with model) and uncorrected signal [%]' })
ylabel('Range [m]')
colormap(jet)
xlim([-80 20]);
grid on
box on

%% Test model on long term (overlap functions)
figure;

T = temp_mean-273.15;
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

set(gcf,'DefaultAxesColorOrder',RGB)

subplot(1,2,1)
hold on

subplot(1,2,2)
hold on

subplot(1,2,2)
hold on
ov_reconstructed2=ones(length(time),length(range));

% colorvec=jet(length(time));
for i=1:length(time)
    temp_of_interest=temp_mean(i)-273.15;
    for j=1:length(range)
        rel_diff = polyval([a(j),b(j),c(j)],temp_of_interest);
        ov_reconstructed2(i,j) = 1./ (rel_diff/100/overlap_ref(j) + 1./overlap_ref(j));
    end
    subplot(1,2,1)
    title('corrections')
%     plot(range,overlap_cor(i,:),'color',colorvec(i,:),'displayname',datestr(time(i)));
    plot(range,overlap_cor(i,:),'displayname',datestr(time(i)));
    subplot(1,2,2)
    title('model')
%     plot(range,ov_reconstructed2(i,:),'color',colorvec(i,:));
    plot(range,ov_reconstructed2(i,:));
  
end

subplot(1,2,1)
plot(range,overlap_ref,'--k','LineWidth',2)
ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200);
colormap('jet');
caxis([15,40]);
hc = colorbar;
% ylabel(hc,'mean interior T during the day [°C]');
ylabel(hc,'median interior T during calculation of overlap correction [°C]');
grid on;box on;
xlabel('Range [m]')

subplot(1,2,2)
plot(range,overlap_ref,'--k','LineWidth',2)
ylim([0 1.2]);
xlim([0 1200]);
set(gca,'xtick',0:200:1200);
colormap('jet');
caxis([15,40]);
hc = colorbar;
ylabel(hc,'interior T [°C]');
grid on;box on;
xlabel('Range [m]')

%% Test model on long term (relative difference w.r.t to set of well-trusted functions)

xlims = [-60 60];
ylabel_str = 'Range [m]';
ylims = [0 1200];
index_range_for_stats = overlap_ref>0.005 & overlap_ref<1;

figure;
subplot(2,3,1)
hold on

subplot(2,3,2)
hold on

subplot(2,3,3)
hold on

subplot(2,3,4)
hold on

subplot(2,3,5)
hold on

subplot(2,3,6)
hold on


ov_reconstructed=ones(length(time),length(range));

rel_diff_all1=zeros(length(time),length(range));
rel_diff_all2=zeros(length(time),length(range));
rel_diff_all3=zeros(length(time),length(range));
rel_diff_all4=zeros(length(time),length(range));
rel_diff_all5=zeros(length(time),length(range));
rel_diff_all6=zeros(length(time),length(range));

% consider only well-trusted corrections
index = index_well_trusted;

T = temp_mean(index);
Treduced = (T-min(T))/(max(T)-min(T))*(length(jet)-1)+1;
RGB = interp1(1:length(jet),jet,Treduced);

set(gcf,'DefaultAxesColorOrder',RGB)


% colorvec=jet(length(index));
colorvec=RGB;
for i=1:length(index)
    
    
    subplot(2,3,1)
    
    %calculate relative difference with Lufft reference
    rel_diff_all1(index(i),:) = (-1./overlap_ref'+1./overlap_cor(index(i),:))./(1./overlap_ref')*100;
    plot(rel_diff_all1(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));
    xlabel_rel_diff_all1 = {'Relative difference';'With Lufft reference'};
    
    
    subplot(2,3,2)
    
    %recalculate relative difference with temperature model
    temp_of_interest=temp_mean(index(i))-273.15;
    for j=1:length(range)
        rel_diff = polyval([a(j),b(j),c(j)],temp_of_interest);
        ov_reconstructed(index(i),j) = 1./ (rel_diff/100/overlap_ref(j) + 1./overlap_ref(j));
    end
    rel_diff_all2(index(i),:) = (-1./ov_reconstructed(index(i),:)+1./overlap_cor(index(i),:))./(1./ov_reconstructed(index(i),:))*100;
    plot(rel_diff_all2(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));
    xlabel_rel_diff_all2 = {'Relative difference';'With Temperature model'};
    
    
    subplot(2,3,3)
    
    %recalculate relative difference with last well-trusted overlap
    if i>1
        rel_diff_all3(index(i),:) = (-1./overlap_cor(index(i-1),:)+1./overlap_cor(index(i),:))./(1./overlap_cor(index(i-1),:))*100;
        plot(rel_diff_all3(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));
        xlabel_rel_diff_all3 = {'Relative difference';'With last well-trusted overlap function'};
    end
    
    
    subplot(2,3,4)
    
    %recalculate relative difference with overlap correction of closest temperature
    index_all = 1:length(index);
    index_all_short = index_all(~(i==(1:length(index))));
    [~,index_min]= min(abs(temp_mean(index(index_all_short))-temp_mean(index(i))));
    rel_diff_all4(index(i),:) = (-1./overlap_cor(index(index_all_short(index_min)),:)+1./overlap_cor(index(i),:))./(1./overlap_cor(index(index_all_short(index_min)),:))*100;
    xlabel_rel_diff_all4 = {'Relative difference';'With well-trusted overlap function of closest Temperature'};
    plot(rel_diff_all4(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));  

    
    subplot(2,3,5)
    
    %recalculate relative difference with overlap correction of closest time
    index_all = 1:length(index);
    index_all_short = index_all(~(i==(1:length(index))));
    [~,index_min]= min(abs(time(index(index_all_short))-time(index(i))));
    rel_diff_all5(index(i),:) = (-1./overlap_cor(index(index_all_short(index_min)),:)+1./overlap_cor(index(i),:))./(1./overlap_cor(index(index_all_short(index_min)),:))*100;
    xlabel_rel_diff_all5 = {'Relative difference';'With well-trusted overlap function of closest time'};
    plot(rel_diff_all5(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));  
    
    
    subplot(2,3,6)
    
    %recalculate relative difference with corrected overlap of specific date
    i_dateofinterest = find(time==datenum(2014,06,16));
    if isempty(i_dateofinterest)
        error('date of interest not found')
    end
    rel_diff_all6(index(i),:) = (-1./overlap_cor(i_dateofinterest,:)+1./overlap_cor(index(i),:))./(1./overlap_cor(i_dateofinterest,:))*100;
    xlabel_rel_diff_all6 = {'Relative difference';['With well-trusted overlap of ' datestr(time(i_dateofinterest))]};
    plot(rel_diff_all6(index(i),:),range,'color',colorvec(i,:),'displayname',datestr(time(index(i))));  
    
end


subplot(2,3,1)
hold on;
plot(median(rel_diff_all1(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all1(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all1(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all1);
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all1(index,index_range_for_stats)) *nansum(nansum(rel_diff_all1(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all1(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all1(index,index_range_for_stats),75)-prctile(rel_diff_all1(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

subplot(2,3,2)
plot(median(rel_diff_all2(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all2(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all2(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all2)
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all2(index,index_range_for_stats)) *nansum(nansum(rel_diff_all2(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all2(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all2(index,index_range_for_stats),75)-prctile(rel_diff_all2(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

subplot(2,3,3)
plot(median(rel_diff_all3(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all3(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all3(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all3)
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all3(index,index_range_for_stats)) *nansum(nansum(rel_diff_all3(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all3(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all3(index,index_range_for_stats),75)-prctile(rel_diff_all3(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

subplot(2,3,4)
plot(median(rel_diff_all4(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all4(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all4(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all4)
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all4(index,index_range_for_stats)) *nansum(nansum(rel_diff_all4(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all4(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all4(index,index_range_for_stats),75)-prctile(rel_diff_all4(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

subplot(2,3,5)
plot(median(rel_diff_all5(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all5(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all5(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all5)
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all5(index,index_range_for_stats)) *nansum(nansum(rel_diff_all5(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all5(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all5(index,index_range_for_stats),75)-prctile(rel_diff_all5(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

subplot(2,3,6)
plot(median(rel_diff_all6(index,:)),range,'-k','linewidth',2)
plot(prctile(rel_diff_all6(index,:),25),range,'--k','linewidth',2)
plot(prctile(rel_diff_all6(index,:),75),range,'--k','linewidth',2)
ylim(ylims);
grid on;box on;
ylabel(ylabel_str)
xlim(xlims)
xlabel(xlabel_rel_diff_all6)
hold on;
plot(xlims,min(range(index_range_for_stats))*ones(1,2),'--k');
plot(xlims,max(range(index_range_for_stats))*ones(1,2),'--k');
RMSE=sqrt(1/numel(rel_diff_all6(index,index_range_for_stats)) *nansum(nansum(rel_diff_all6(:,index_range_for_stats).^2,2)));
MEANmedian=mean(median(rel_diff_all6(index,index_range_for_stats)));
MEANiqr=mean(prctile(rel_diff_all6(index,index_range_for_stats),75)-prctile(rel_diff_all6(index,index_range_for_stats),25));
title(['RMSE=' num2str(RMSE) ',MEANmedian=' num2str(MEANmedian) ',MEANiqr=' num2str(MEANiqr)]);
colormap('jet');
caxis([15 40]);

% return


%% Load test data

% date = '20130417';
date = '20140616';

[chm,chminfo] = get_chm15k_from_files('pay',[date,'000000'],datestr(datenum(date,'yyyymmdd')+1,'yyyymmddHHMMSS'),folder_ncdata);
% [chm,chminfo]=readcorrectlyncfile3('pay',date,folder_ncdata);

RCS_raw = chm.beta_raw;
RCS_dailycorrection = chm.beta_raw./repmat(overlap_cor(find(time==datenum(date,'yyyymmdd')),:)',1,size(chm.beta_raw,2)).*repmat(overlap_ref,1,size(chm.beta_raw,2));

% Calc Grad
grad_raw = NaN(size(RCS_raw));
for j=1:length(chm.time)
    grad_raw(:,j) = 1/(2*chm.range_gate)*conv(RCS_raw(:,j),[1 0 -1],'same');
end
grad_dailycorrection = NaN(size(RCS_dailycorrection));
for j=1:length(chm.time)
    grad_dailycorrection(:,j) = 1/(2*chm.range_gate)*conv(RCS_dailycorrection(:,j),[1 0 -1],'same');
end

% plot internal temperature
figure;
plot(chm.time,chm.temp_int-273.15);datetick;ylim([15 40]);box on;grid on;
xlabel('Time [UT]');ylabel('Internal Temperature [°C]');title(['pay ' date]);



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
else warning('No file')
end

%% calculate RCS and grad(RCS) with the model

%correct according to T

ind_last_as_reference = find(b>0 & range <range(find(overlap_ref<0.01,1,'last')),1,'last');
if ~isempty(ind_last_as_reference)
    a(range<=range(ind_last_as_reference)) = 0;
    c(range<=range(ind_last_as_reference)) = 0;
    b(range<=range(ind_last_as_reference)) = 0;
end

ov_rec_all = NaN(size(RCS_raw));
for i=1:find(chm.range<=1200,1,'last')
    rel_diff = polyval([a(i),b(i),c(i)],chm.temp_int-273.15);
    ov_rec_all(i,:) = 1./ (rel_diff/100/overlap_ref(i) + 1./overlap_ref(i));
end
ov_rec_all(find(chm.range<=1200,1,'last')+1:end,:) = repmat(overlap_ref(find(chm.range<=1200,1,'last')+1:end),1,size(RCS_raw,2));

RCS_corr = RCS_raw./ov_rec_all.*repmat(overlap_ref,1,size(RCS_raw,2));

% Calc Grad
grad_corr = NaN(size(RCS_corr));
for j=1:length(chm.time)
    grad_corr(:,j) = 1/(2*chm.range_gate)*conv(RCS_corr(:,j),[1 0 -1],'same');
end

%% Calculate PBL with simple code

[gradients_raw,xtime] = simplePBLdetection(chm,overlap_ref,overlap_ref);
[gradients_dailycorrection,xtime] = simplePBLdetection(chm,overlap_ref,overlap_cor(find(time==datenum(date,'yyyymmdd')),:)');
[gradients_corr,xtime] = simplePBLdetection(chm,overlap_ref,ov_rec_all);

% Plot PBL
figure;

subplot(3,1,1);
gradients = gradients_raw;
correction_type = 'raw';
list_plots = [];
list_legends = {};
hold on;
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(5,:),'.','color',0.5*ones(1,3));
list_legends{end+1} = 'clouds';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'ok');
list_legends{end+1} = '1. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
list_legends{end+1} = '2. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
list_legends{end+1} = '3. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+r');
list_legends{end+1} = 'min. gradient';
legend(list_plots,list_legends,'location','eastoutside');
ylim([0 2500]);
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
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'ok');
list_legends{end+1} = '1. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
list_legends{end+1} = '2. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
list_legends{end+1} = '3. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+r');
list_legends{end+1} = 'min. gradient';
legend(list_plots,list_legends,'location','eastoutside');
ylim([0 2500]);
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
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(4,:),'ok');
list_legends{end+1} = '1. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(3,:),'>g');
list_legends{end+1} = '2. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(2,:),'<b');
list_legends{end+1} = '3. strongest gradient';
list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime,gradients(1,:),'+r');
list_legends{end+1} = 'min. gradient';
legend(list_plots,list_legends,'location','eastoutside');
ylim([0 2500]);
xlim([datenum(date,'yyyymmdd') datenum(date,'yyyymmdd')+1]);
set(gca,'xtick',datenum(date,'yyyymmdd'):4/24:datenum(date,'yyyymmdd')+1);
datetick('x','HH:MM','keepticks','keeplimits');
title(['pay ' date ' - ' correction_type]);
grid on;box on;
xlabel('Time [UT]');
ylabel('Range [m]');

%% Plot

RCS_list = {RCS_raw,RCS_dailycorrection,RCS_corr};
grad_list = {grad_raw,grad_dailycorrection,grad_corr};
gradients_list = {gradients_raw,gradients_dailycorrection,gradients_corr};
correction_type_list = {'raw','daily correction','model correction'};

for k=1:length(RCS_list)
    RCS = RCS_list{k};
    grad = grad_list{k};
    gradients = gradients_list{k};
    title_str = [info.chm '/' info.tub '/' station ' - ' date ' - ' correction_type_list{k}];


    offset = datenum(2000,1,1)-1;
    disp_clouds = true;
    disp_pbl = true;
    % Choose pcolor parameters
    dxticks = 4;%in hours
    fz = 14;% FontSize
    ylims = [0 2500];
    yticks = ylims(1):100:ylims(2);
    clims = [4.5 6];
    clims_grad = [-1000 1000];
    % load('ypcmap2','cmap');
    
    figure('Units','normalized','Position',[0.005 0.05 0.99 0.85],'Name',[ 'Pcolor ' correction_type_list{k}])
    
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
        %     for j=1:3
        %         pblh = chm.pbl(j,:);
        %         pblh(pblh <= chm.cho) = NaN;
        %         hpbl = line(chm.time-offset,pblh-chm.cho,'LineStyle','none','Marker','.','Color','k');
        %     end
        %     list_plots(end+1) = hpbl;
        %     list_legends{end+1} = 'PBL';
        hold on;
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(4,:),'.k','MarkerSize',16);
        list_legends{end+1} = '1. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(3,:),'.m','MarkerSize',16);
        list_legends{end+1} = '2. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(2,:),'.b','MarkerSize',16);
        list_legends{end+1} = '3. strongest gradient';
        list_plots(end+1) = plot(datenum(date,'yyyymmdd')+xtime-offset,gradients(1,:),'.r','MarkerSize',16);
        list_legends{end+1} = 'min. gradient';
    end
    % add area where the overlap function
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
%         plot([unique_time_start(i)  unique_time_end(i)]-offset,repmat(range_start_min(i),2,1),'--k')
%         plot([unique_time_start(i)  unique_time_end(i)]-offset,repmat(range_start_max(i),2,1),'--k')
%         plot([unique_time_start(i)  unique_time_start(i)]-offset,[range_start_min(i) range_start_max(i)],'--k')
%         plot([unique_time_end(i)  unique_time_end(i)]-offset,[range_start_min(i) range_start_max(i)],'--k')
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
    list_legends{end+1} = 'Ref. box';

    %Vertical bar if there is a new box    
    for i=2:length(chm.time)
        if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i-1))
            plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
        end
    end
    
    %Vertical bar if there is a box is ending
    for i=1:length(chm.time)-1
        if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i+1))
            plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
        end
    end
    
    % Plot Clouds
    if(disp_clouds)
        for j=1:3
            cbh = chm.cbh(j,:);
            cbh(cbh <= chm.cho) = NaN;
            hcbh = line(chm.time-offset,cbh-chm.cho,'LineStyle','none','Marker','.','MarkerSize',16,'Color',[0.5 0.5 0.5]);
        end
        list_plots(end+1) = hcbh;
        list_legends{end+1} = 'CBH';
    end
    if ~isempty(list_plots)
        legend(list_plots,list_legends,'Orientation','horizontal');
    end
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
    % daspect(aspect_ratio);
    set(gca,'FontSize',fz);
    
    set(gcf,'PaperPositionMode','auto');

end

