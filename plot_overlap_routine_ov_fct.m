function plot_overlap_routine_ov_fct(varargin)

nargs = length(varargin);
if nargs==5
    L1 = varargin{1};
    chminfo = varargin{2};
    ov = varargin{3};
    ovp_fc_ok_mat_final = varargin{4};
    ovp_fc_final = varargin{5};
    file_out = [];
elseif nargs==6
    L1 = varargin{1};
    chminfo = varargin{2};
    ov = varargin{3};
    ovp_fc_ok_mat_final = varargin{4};
    ovp_fc_final = varargin{5};
    file_out = varargin{6};
else
    warning('Wrong input arguments');
    return;
end


% for j=1:length(chminfo.Attributes)
%     if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
%         device_name = chminfo.Attributes(j).Value;
%     elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
%         serlom = chminfo.Attributes(j).Value;
%     elseif strcmpi(chminfo.Attributes(j).Name,'location')
%         location = chminfo.Attributes(j).Value;
%     end
% end
if isfield(L1,'instrument_serial_number')
    device_name = L1.instrument_serial_number;
elseif isfield(L1,'device_name')
    device_name = L1.device_name{1};
end
if isfield(L1,'optical_module_id')
    serlom = L1.optical_module_id;
elseif isfield(L1,'serlom')
    serlom = L1.serlom{1};
end
if isfield(L1,'site_location')
location = L1.site_location;
elseif isfield(L1,'location')
location = L1.location{1};
end
stn = [device_name '/' serlom,'/' location];
date_yyyymmdd = datestr(floor(L1.time(1)),'yyyymmdd');


% Choose plot parameters
fz = 10;% FontSize

if isempty(file_out)
    figure('Color','w','Position',[0 0 1920 1024],'Visible','on');
else
    figure('Color','w','Position',[0 0 1920 1024],'Visible','off');
end

hold on;

list_plots = [];
list_legends = {};

if ~isempty(ovp_fc_ok_mat_final)
    colors = hsv(length(ovp_fc_ok_mat_final));
    for k=1:size(ovp_fc_ok_mat_final,2)
        plot(L1.range,ovp_fc_ok_mat_final(:,k),'-','Color',colors(k,:));
    end
end

if ~isempty(ov)
    list_plots(end+1,1) = plot(L1.range,ov,'--k','LineWidth',2);
    list_legends{end+1,1} = [serlom ' Lufft'];
end

if ~isempty(ovp_fc_final)
    list_plots(end+1,1) = plot(L1.range,ovp_fc_final,'-^k','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k');
    list_legends{end+1,1} = [serlom ' corrected'];
end

set(gca,'Ylim',[0 1.5],'YTick',0:0.1:1.5);
set(gca,'Xlim',[0 2500],'XTick',0:50:2500);
grid on;
xlabel('Range (m)','FontSize',fz);
title([stn ', ' date_yyyymmdd],'FontSize',fz);

legend(list_plots,list_legends);

set(gca,'FontSize',fz);

set(gcf,'PaperPositionMode','auto');

if ~isempty(file_out)
    print( '-dpng','-r150','-loose', file_out);
    close;
end

end