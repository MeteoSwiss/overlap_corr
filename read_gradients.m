clear variables;clc;close all;

%% Inputs
info.start_day  = 1;
info.start_month= 1;
info.start_year = 2013;

info.end_day  =  31;
info.end_month=  12;
info.end_year =  2014;
corrections_to_analyze='all';

folder='C:\DATA\MATLAB\ceilometer\Overlap-function\Outputs\';

%% Load data
time_vec=datenum(info.start_year,info.start_month,info.start_day):datenum(info.end_year,info.end_month,info.end_day);
data=NaN(length(time_vec)*145,21);
k=1;
for t=1:length(time_vec)
    file=[folder 'PBL_' datestr(time_vec(t),'yyyymmdd') '_' corrections_to_analyze  '.csv' ];
    if exist(file,'file')>0
        disp(['Read file: ' file])
        data_daily=dlmread (file);
        data(k:k+length(data_daily)-1,:)=data_daily;
        k=k+length(data_daily);
    else
        disp(['no file  :' datestr(time_vec(t),'yyyymmdd')])
    end
end

data(all(isnan(data),2),:)=[];

%% Parse
disp('Parsing')
time=datenum(data(:,1:6));
grad_raw=data(:,7:11);
gradients_dailycorrection=data(:,12:16);
gradients_model=data(:,17:21);

%% Plot 1 : histogram for strongest gradients
disp('Plotting')
figure('Units','normalized','Position',[0.005 0.05 0.8 0.6],'Name','hist strongest grad')
size_bin=15;
ymax=100;
subplot(1,3,1)
hist(grad_raw(:,2:4),0:size_bin:2500)
title('raw Overlap')
ylim([0 ymax])
xlim([0 2500])

subplot(1,3,2)
hist(gradients_dailycorrection(:,2:4),0:size_bin:2500)
title({'3 strongest gradients';'Daily correction'})
legend('3rd grad.','2nd grad.','Strongest Grad.')
ylim([0 ymax])
xlim([0 2500])

subplot(1,3,3)
hist(gradients_model(:,2:4),0:size_bin:2500)
title('Model correction')
ylim([0 ymax])
xlim([0 2500])

%% Plot1-a : histogram for strongest gradient
figure('Units','normalized','Position',[0.005 0.05 0.8 0.6],'Name','hist strongest grad')
ymax=2500;

subplot(1,3,1)
hist(grad_raw(:,4),0:size_bin:2500);
[counts,centers] = hist(grad_raw(:,4),0:size_bin:2500);
title('raw Overlap')
ylim([0 ymax])
xlim([0 2500])

subplot(1,3,2)
hist(gradients_dailycorrection(:,4),0:size_bin:2500);
[counts2,centers2] = hist(gradients_dailycorrection(:,4),0:size_bin:2500);
title({'Strongest gradient';'Daily correction'})
ylim([0 ymax])
xlim([0 2500])

subplot(1,3,3)
hist(gradients_model(:,4),0:size_bin:2500);
[counts3,centers3] = hist(gradients_model(:,1),0:size_bin:2500);

title('Model correction')
ylim([0 ymax])
xlim([0 2500])
%% Plot1-b : histogram for strongest gradient
figure('Units','normalized','Position',[0.2 0.2 0.8 0.6],'Name','hist strongest grad')
ymax=8000;
size_bin=15;
hist([grad_raw(:,4) gradients_model(:,4)],0:size_bin:2500);
h = findobj(gca,'Type','patch');
h(1).FaceColor ='g';
h(2).FaceColor ='r';
h(1).EdgeColor ='g';
h(2).EdgeColor ='r';
% h(1).FaceAlpha	 =1;

title('Strongest gradient (2013-2014)')
ylim([0 ymax])
xlim([0 2500])
xlabel('Altitude above ground (m)')
ylabel('Number of occurences')
legend('Uncorrected','Corrected by the model')
view(90,-90)
set(gca,'fontsize',16)

xprint('hist strongest gradient','png')

%% Plot 2 : histogram for lowest gradient
figure('Units','normalized','Position',[0.005 0.05 0.8 0.6],'Name','hist Lowest grad')
subplot(1,3,1)
hist(grad_raw(:,1))
title('raw Overlap')

subplot(1,3,2)
hist(gradients_dailycorrection(:,1))
title({'Lowest gradient';'Daily correction'})
subplot(1,3,3)
hist(gradients_model(:,1))
title('Model correction')


%% Plot 3 : Scatter strongest gradients
figure
subplot(1,2,1)
hold on
% scatter(gradients_model(:,2),grad_raw(:,2))
% scatter(gradients_model(:,3),grad_raw(:,3))
scatter(gradients_model(:,4),grad_raw(:,4))
xlabel('Model correction')
ylabel('raw Overlap')


subplot(1,2,2)
hold on
% scatter(gradients_model(:,2),gradients_dailycorrection(:,2))
% scatter(gradients_model(:,3),gradients_dailycorrection(:,3))
scatter(gradients_model(:,4),gradients_dailycorrection(:,4))
xlabel('Model correction')
% legend('3rd grad.','2nd grad.','Strongest Grad.')
ylabel('Daily cor')
% 
% figure
% freq_scatter_v5(gradients_model(:,4),gradients_dailycorrection(:,4),15,15);
% 
% figure
% freq_scatter_v5(gradients_model(:,4),grad_raw(:,4),15,15);
%% Plot 3 : Scatter strongest gradients
% This figure is only plot for case studies shorter than one week
% (optimization)
if length(time_vec)<8
figure
subplot(1,2,1)
hold on
% scatter(gradients_model(:,2),grad_raw(:,2))
% scatter(gradients_model(:,3),grad_raw(:,3))
% scatter(gradients_model(:,1),grad_raw(:,1))
hold on
for i=1:size(gradients_model,1)
    if gradients_model(i,1)~=grad_raw(i,1)
        plot(gradients_model(i,1),grad_raw(i,1),'sb','displayname',datestr(time(i)))
    end
end
plot(gradients_model(gradients_model(:,1)==grad_raw(:,1),1),grad_raw(gradients_model(:,1)==grad_raw(:,1),1),'sb')
xlabel('Model correction')
ylabel('raw Overlap')
ylim([ 0 2500])



subplot(1,2,2)


% scatter(gradients_model(:,2),gradients_dailycorrection(:,2))
% scatter(gradients_model(:,3),gradients_dailycorrection(:,3))
% scatter(gradients_model(:,1),gradients_dailycorrection(:,1))

hold on
for i=1:size(gradients_model,1)
    if gradients_model(i,1)~=gradients_dailycorrection(i,1)
        plot(gradients_model(i,1),gradients_dailycorrection(i,1),'sb','displayname',datestr(time(i)))
    end
end
plot(gradients_model(gradients_model(:,1)==gradients_dailycorrection(:,1),1),...
    gradients_dailycorrection(gradients_model(:,1)==gradients_dailycorrection(:,1),1),'sk')
xlabel('Model correction')
% legend('3rd grad.','2nd grad.','Strongest Grad.')
ylabel('Daily cor')
title('lowest')

% Strongest
figure
% scatter(gradients_model(:,2),gradients_dailycorrection(:,2))
% scatter(gradients_model(:,3),gradients_dailycorrection(:,3))
hold on
for i=1:size(gradients_model,1)
    if gradients_model(i,4)~=gradients_dailycorrection(i,4)
        plot(gradients_model(i,4),gradients_dailycorrection(i,4),'sb','displayname',datestr(time(i)))
    end
end
plot(gradients_model(gradients_model(:,4)==gradients_dailycorrection(:,4),4),...
    gradients_dailycorrection(gradients_model(:,4)==gradients_dailycorrection(:,4),4),'sk')
% scatter(gradients_model(:,1),gradients_dailycorrection(:,1))
xlabel('Model correction')
% legend('3rd grad.','2nd grad.','Strongest Grad.')
ylabel('Daily cor')
title('Strongest')
end

%% daily variability
tmp=datevec(time);
hour_vec=tmp(:,4);
H=0:23;
PBL_uncor=NaN(size(H));
PBL_model=NaN(size(H));
for i=1:length(H)
    PBL_uncor(i)=nanmedian(grad_raw(hour_vec==i,4));
    PBL_model(i)=nanmedian(gradients_model(hour_vec==i,4));
end
figure('Units','normalized','Position',[0.01 0.2 0.8 0.6],'Name','hist Lowest grad')

subplot(1,2,1)
plot(H,PBL_uncor,'r')
hold on
plot(H,PBL_model,'g')
xlim([0 23])
legend('Uncorrected','Corrected with model')
xlabel('Hour UTC')
ylabel('PBL Height (m a.g.l.)')
set(gca,'fontsize',24)

% Monthly variability
month_vec=tmp(:,2);
M=1:12;
PBL_uncor_monthly=NaN(size(M));
PBL_model_monthly=NaN(size(M));
for i=1:length(M)
    PBL_uncor_monthly(i)=nanmedian(grad_raw(month_vec==i,4));
    PBL_model_monthly(i)=nanmedian(gradients_model(month_vec==i,4));
end
subplot(1,2,2)

plot(M,PBL_uncor_monthly,'r')
hold on
plot(M,PBL_model_monthly,'g')
legend('Uncorrected','corrected with model')
xlabel('Month')
ylabel('PBL Height (m a.g.l.)')
xlim([1 12])
set(gca,'fontsize',24)