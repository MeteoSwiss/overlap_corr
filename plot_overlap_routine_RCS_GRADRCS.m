function handles = plot_overlap_routine_RCS_GRADRCS(varargin)

nargs = length(varargin);
if nargs==4
    chm = varargin{1};
    chminfo = varargin{2};
    RCS = varargin{3};
    grad_raw = varargin{4};
    file_out = [];
elseif nargs==5;
    chm = varargin{1};
    chminfo = varargin{2};
    RCS = varargin{3};
    grad_raw = varargin{4};
    file_out = varargin{5};
elseif nargs==6;
    chm = varargin{1};
    chminfo = varargin{2};
    RCS = varargin{3};
    grad_raw = varargin{4};
    RCS_corr = varargin{5};
    grad_corr = varargin{6};
    file_out = [];
elseif nargs==7
    chm = varargin{1};
    chminfo = varargin{2};
    RCS = varargin{3};
    grad_raw = varargin{4};
    RCS_corr = varargin{5};
    grad_corr = varargin{6};
    file_out = varargin{7};
elseif nargs == 8
    chm = varargin{1};
    chminfo = varargin{2};
    RCS = varargin{3};
    grad_raw = varargin{4};
    RCS_corr = varargin{5};
    grad_corr = varargin{6};
    file_out = varargin{7};
    result = varargin{8};
else
    warning('Wrong input arguments');
    return;
end


for j=1:length(chminfo.Attributes)
    if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
        device_name = chminfo.Attributes(j).Value;
    elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
        serlom = chminfo.Attributes(j).Value;
    elseif strcmpi(chminfo.Attributes(j).Name,'location')
        location = chminfo.Attributes(j).Value;
    end
end

stn = [device_name '/' serlom,'/' location];
date_yyyymmdd = datestr(floor(chm.time(1)),'yyyymmdd');

disp_clouds = true;
disp_pbl = false;
disp_overlap_calculation_areas = true;

offset = datenum(2000,1,1)-1;
            
if nargs==4 || nargs==5

    % Choose pcolor parameters
    aspect_ratio = [1 10000 1];% 'auto'
    dxticks = 1;%in hours
    fz = 10;% FontSize
    ylims = [0 2500];
    yticks = ylims(1):100:ylims(2);
    clims = [5 5.5];
    clims_grad = [-1000 1000];
    load('ypcmap2','cmap');
    

    if isempty(file_out)
        handles.figure_main = figure('Color','w','Position',[0 0 1920 1024],'Visible','on');
    else
        handles.figure_main = figure('Color','w','Position',[0 0 1920 1024],'Visible','off');
    end
        
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
    hc = colorbar;
    xlabel('Time UT [h]','FontSize',fz);
    ylabel(hc,'log10(abs(beta\_raw))','FontSize',fz);
    ylabel('Range (m)','FontSize',fz);
    colormap(cmap)
    daspect(aspect_ratio);
    set(gca,'FontSize',fz);
    
    list_plots = [];
    list_legends = {};
    if(disp_pbl)
        for j=1:3
            pblh = chm.pbl(j,:);
            pblh(pblh <= chm.cho) = NaN;
            hpbl = line(chm.time-offset,pblh-chm.cho,'LineStyle','none','Marker','.','Color','k');
        end
        list_plots(end+1) = hpbl;
        list_legends{end+1} = 'PBL';
    end
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
        legend(list_plots,list_legends);
    end

    title([stn ' ' date_yyyymmdd],'FontSize',fz);

    handles.axes_pcolor_grad = subplot(2,1,2);
    handles.surface_pcolor_grad = pcolor(chm.time-offset,chm.range,grad_raw);
    xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset];
    xticks = xlims(1)-mod(mod(xlims(1),1)*24,dxticks)/24:dxticks/24:xlims(2);
    set(gca,'XLim',xlims,'XTick',xticks,'XMinorTick','on');
    set(gca,'YLim',ylims,'YTick',yticks,'YMinorTick','on');
    grid on;
    set(gca,'Layer','top');
    datetick('x','HH:MM','keepticks','keeplimits');
    shading flat
    caxis(clims_grad)
    hc = colorbar;
    xlabel('Time UT [h]','FontSize',fz);
    ylabel(hc,'grad(beta\_raw)','FontSize',fz);
    ylabel('Range (m)','FontSize',fz);
    colormap(jet)
    daspect(aspect_ratio);
    set(gca,'FontSize',fz);
    
    set(gcf,'PaperPositionMode','auto');
    
    if ~isempty(file_out)
        print( '-dpng','-r150','-loose', file_out);
        close;
    end
  
end

if nargs==6 || nargs==7 || nargs==8

    % Choose pcolor parameters
    aspect_ratio = [1 5000 1];% 'auto'
    dxticks = 4;%in hours
    fz = 14;% FontSize
    ylims = [0 2500];
    yticks = ylims(1):100:ylims(2);
    clims = [4.5 6];
    clims_grad = [-1000 1000];
    load('ypcmap2','cmap');

    
    if isempty(file_out)
        handles.figure_main = figure('Color','w','Position',[0 0 1920 1024],'Visible','on');
    else
        handles.figure_main = figure('Color','w','Position',[0 0 1920 1024],'Visible','off');
    end

    handles.axes_pcolor_RCS = subplot(2,2,1);
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
    colormap(cmap)
    hc = colorbar;
    xlabel('Time UT [h]','FontSize',fz);
    ylabel(hc,'log10(abs(beta\_raw))','FontSize',fz);
    ylabel('Range (m)','FontSize',fz);
    daspect(aspect_ratio);
    list_plots = [];
    list_legends = {};
    if(disp_pbl)
        for j=1:3
            pblh = chm.pbl(j,:);
            pblh(pblh <= chm.cho) = NaN;
            hpbl = line(chm.time-offset,pblh-chm.cho,'LineStyle','none','Marker','.','Color','k');
        end
        list_plots(end+1) = hpbl;
        list_legends{end+1} = 'PBL';
    end
    if(disp_clouds)
        for j=1:3
            cbh = chm.cbh(j,:);
            cbh(cbh <= chm.cho) = NaN;
            hcbh = line(chm.time-offset,cbh-chm.cho,'LineStyle','none','Marker','.','MarkerSize',16,'Color',[0.5 0.5 0.5]);
        end
        list_plots(end+1) = hcbh;
        list_legends{end+1} = 'CBH';
    end
    
    % add area where the final overlap functions where calculated
    if(disp_overlap_calculation_areas)
        time_start = result.time_start(result.index_final);
        time_end = result.time_end(result.index_final);
        range_start = result.range_start(result.index_final);
        range_end = result.range_end(result.index_final);
        
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
            index=chm.time(i)>=unique_time_start & chm.time(i)<=unique_time_end;
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
        
        %Vertical bar if there is a box is ending
        for i=1:length(chm.time)-1
            if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i+1))
                plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
            end
        end

    end
    
    
    
    if ~isempty(list_plots)
%         legend(list_plots,list_legends);
    end
    set(gca,'FontSize',fz);

    title([stn ' ' date_yyyymmdd],'FontSize',fz);

    handles.axes_pcolor_grad = subplot(2,2,3);
    handles.surface_pcolor_grad = pcolor(chm.time-offset,chm.range,grad_raw);
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
    daspect(aspect_ratio);
    set(gca,'FontSize',fz);

    % *************************************************************************

    handles.axes_pcolor_RCS_corr = subplot(2,2,2);
    handles.surface_pcolor_RCS_corr = pcolor(chm.time-offset,chm.range,log10(abs(RCS_corr)));
    xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset];
    xticks = xlims(1)-mod(mod(xlims(1),1)*24,dxticks)/24:dxticks/24:xlims(2);
    set(gca,'XLim',xlims,'XTick',xticks,'XMinorTick','on');
    set(gca,'YLim',ylims,'YTick',yticks,'YMinorTick','on');
    grid on;
    set(gca,'Layer','top');
    datetick('x','HH:MM','keepticks','keeplimits');
    shading flat
    caxis(clims)
    colormap(cmap)
    hc = colorbar;
    xlabel('Time UT [h]','FontSize',fz);
    ylabel(hc,'log10(abs(beta\_raw)) corrected','FontSize',fz);
    ylabel('Range (m)','FontSize',fz);
    daspect(aspect_ratio);
    set(gca,'FontSize',fz);
    
    list_plots = [];
    list_legends = {};
    if(disp_pbl)
        for j=1:3
            pblh = chm.pbl(j,:);
            pblh(pblh <= chm.cho) = NaN;
            hpbl = line(chm.time-offset,pblh-chm.cho,'LineStyle','none','Marker','.','Color','k');
        end
        list_plots(end+1) = hpbl;
        list_legends{end+1} = 'PBL';
    end
    if(disp_clouds)
        for j=1:3
            cbh = chm.cbh(j,:);
            cbh(cbh <= chm.cho) = NaN;
            hcbh = line(chm.time-offset,cbh-chm.cho,'LineStyle','none','Marker','.','MarkerSize',16,'Color',[0.5 0.5 0.5]);
        end
        list_plots(end+1) = hcbh;
        list_legends{end+1} = 'CBH';
    end
    % add area where the final overlap functions where calculated
    if(disp_overlap_calculation_areas)
        time_start = result.time_start(result.index_final);
        time_end = result.time_end(result.index_final);
        range_start = result.range_start(result.index_final);
        range_end = result.range_end(result.index_final);
        
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
            index=chm.time(i)>=unique_time_start & chm.time(i)<=unique_time_end;
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
        
        %Vertical bar if there is a box is ending
        for i=1:length(chm.time)-1
            if ~isnan(range_start_min_all(i)) && isnan(range_start_min_all(i+1))
                plot([chm.time(i)  chm.time(i)]-offset,[range_start_min_all(i) range_start_max_all(i)],'--k')
            end
        end

    end
    
    if ~isempty(list_plots)
%         legend(list_plots,list_legends);
    end

    title([stn ' ' date_yyyymmdd ' corrected'],'FontSize',fz);

    handles.axes_pcolor_grad_corr = subplot(2,2,4);
    handles.surface_pcolor_grad_corr = pcolor(chm.time-offset,chm.range,grad_corr);
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
    ylabel(hc,'grad(beta\_raw) corrected','FontSize',fz);
    ylabel('Range (m)','FontSize',fz);
    daspect(aspect_ratio);
    set(gca,'FontSize',fz);

    set(gcf,'PaperPositionMode','auto');

    if ~isempty(file_out)
        print( '-dpng','-r150','-loose', file_out);
        close;
    end

end