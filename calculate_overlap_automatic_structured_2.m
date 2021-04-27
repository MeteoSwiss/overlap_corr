function [ovp_fc_ok_mat,ovp_fc_ok_mat_time,ovp_fc_ok_mat_range,ind_ovp_fc_ok_good,ind_ovp_fc_ok_final,ovp_fc_final,error_string] = calculate_overlap_automatic_structured_2(chm,chminfo,RCS,ov)

    date_yyyymmdd = datestr(floor(chm.time(1)),'yyyymmdd');
    
    min_time_interval_length = 30;% in minutes
    max_time_interval_length = 30;
    d_time_interval_length = 30;
    min_fit_time = datenum(date_yyyymmdd,'yyyymmdd');
    max_fit_time = datenum(date_yyyymmdd,'yyyymmdd')+1;
    % min_fit_time = datenum([date_yyyymmdd,'093000'],'yyyymmddHHMMSS');max_fit_time = min_fit_time+min_time_interval_length/(24*60);
    d_fit_time = 5;% in minutes {5,15}
    
    %min_fit_range = 600;
    min_fit_range = chm.range(find(ov>0.6,1,'first')); % Default 0.8
    max_fit_range = 1200;
    d_fit_range = 15;%{10}
    min_fit_length = 150;%{200,100}
    max_fit_length = max_fit_range-min_fit_range;
    d_fit_length = 15;%{20}
    min_range_std_over_mean = 225;
    dt_sliding_variance = 10;
    max_std_over_mean = 0.015; % Default 0.01
    max_relgrad = 0.08;%{0.055,0.05} Default 0.05
    max_relgrad_mean = 0.025;%{0.015,0.01} Default 0.015
    min_expected_slope = -2*1/log(10)*10*1e-6;
    max_expected_slope = -2*1/log(10)*0.1*1e-6;
    min_expected_zero_fit_value = 4.75;%{5}
    max_expected_zero_fit_value = 6;
%     thresh_resid = 0.005; % not used
    thresh_resid_rel = 0.0005;%{0.001,0.0005}
%     thresh_resid_whole_zone = 0.01;
    thresh_resid_whole_zone = 0.15;%{0.2,0.15}
    %min_overlap_valid = 750;%{750,1500}
    min_overlap_valid = chm.range(find(ov>=1,1,'first'));
    thresh_overlap_valid_rel_error = 0.01;%{0.015,0.01}
    % first_range_gradY = 500;
    first_range_gradY = min_fit_range;
    min_nb_ok_candidates = 1;%{5,1}
    max_overlap_value = 1.01;
    min_slope = -2.5*1e-4;%{-0.5*1e-4,-2.5*1e-4}
    sgolay_width = 5;%{11,5}
    sgolay_ord = 3;
    whiskers_length = 3;%{3,1.5}
    min_nb_samples = 15;%{10,5,1}
    good_samples_proportion = 1.0;
    min_nb_good_samples = floor(good_samples_proportion*min_nb_samples);%{10,5,1}
    min_nb_good_samples_after_outliers_removal = 10;%{10,5}
    %overestimated = false;
    min_nb_samples_for_skipping_good_test = 100;%{1000,100};
    
    ovp_fc_ok = {};
    ovp_fc_ok_time_interval = {};
    ovp_fc_ok_range_interval = {};
    
    ind_ovp_fc_ok_good = [];
    ind_ovp_fc_ok_final = [];
    
    ovp_fc_final = [];
    ovp_fc_ok_mat_final = [];
    ovp_fc_ok_mat_time_final = [];
    ovp_fc_ok_mat_range_final = [];
    
    ovp_fc_ok_mat = [];
    ovp_fc_ok_mat_time = [];
    ovp_fc_ok_mat_range = [];
    
%     % Each element of %results consists of a cell with 5 elements:
%     % 1) time-interval in datenum format
%     % 2) range-interval if any considered
%     % 3) overlap function if any calculated
%     % 4) error tag
%     %results = {};
 
    
    %% 1. Try succesive time intervals
    for time_interval_length = min_time_interval_length:d_time_interval_length:max_time_interval_length
%         disp(['Time interval length is ' num2str(time_interval_length) ' minutes']);

        error_time = NaN(length(min_fit_time:d_fit_time/(24*60):max_fit_time),1);
        error_flag = zeros(length(min_fit_time:d_fit_time/(24*60):max_fit_time),1);
        error_string = cell(length(min_fit_time:d_fit_time/(24*60):max_fit_time),1);
        ok_candidates = cell(length(min_fit_time:d_fit_time/(24*60):max_fit_time),1);
        it=0;
        for fit_time = min_fit_time:d_fit_time/(24*60):max_fit_time
            it = it+1;
            if fit_time+time_interval_length/(24*60) <= max_fit_time
                index_time = chm.time>=fit_time & chm.time<=fit_time+time_interval_length/(24*60);
                
                disp([datestr(fit_time) '-' datestr(fit_time+time_interval_length/(24*60))]);
                error_string{it} = [datestr(fit_time,'HH:MM:SS') '-' datestr(fit_time+time_interval_length/(24*60),'HH:MM:SS') '\n'];
                error_time(it) = fit_time;
               
                
                % Check if at least one profile in timeseries
                t = chm.time(index_time);
                
                if isempty(t)
                    error_flag(it) = 1;
                    error_string{it} = [error_string{it},'no data in the prescribed time interval'];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning('no data in the prescribed time interval');
                    continue;
                end
                
                error_string{it} = [error_string{it},[datestr(t(1),'HH:MM:SS') '-' datestr(t(end),'HH:MM:SS') ', ' num2str(length(t)) ' timesteps' '\n']];
                avg_t = chm.average_time(index_time)/1000/24/3600;
                tol = 0.1/24/3600;
                
                % if any(diff(t)>=avg_t(2:end)+tol)
                %     warning('no continuity in data for the prescribed time interval');
                %     continue;
                % end
                
                if any(chm.sci(index_time) ~= 0)
                    error_flag(it) = 2;
                    error_string{it} = [error_string{it},'at least one sci~=0 for the prescribed time interval'];
                                      
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning('at least one sci~=0 for the prescribed time interval');
                    continue;
                end
                
                cbh = chm.cbh(1,index_time);
%                 min_cloud_range = min(cbh(cbh>=chm.cho));
                                min_cloud_range = min(cbh(:));

                if isempty(min_cloud_range)
                    min_cloud_range = chm.range(end);
                end
                
                if isfield(chm,'mxd')
                    mxd = chm.mxd(index_time);
                    min_mxd_range = min(mxd(mxd>=chm.cho));
                    if isempty(min_mxd_range)
                        
                        min_mxd_range = chm.range(end);
                    end
                else
                    % No MXD in E-PROFILE
                    min_mxd_range = chm.range(end);
                end
                
                max_available_fit_range = min_cloud_range;
                if max_available_fit_range < min_fit_range+min_fit_length
                    error_flag(it) = 3;
                    error_string{it} = [error_string{it},['low cloud: ' num2str(max_available_fit_range) 'm, should be >=' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low cloud: ' num2str(max_available_fit_range) 'm, should be >=' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                max_available_fit_range = min(max_available_fit_range,min_mxd_range);
                if max_available_fit_range < min_fit_range+min_fit_length
                    error_flag(it) = 4;
                    error_string{it} = [error_string{it},['low mxd: ' num2str(max_available_fit_range) 'm, should be >=' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low mxd: ' num2str(max_available_fit_range) 'm, should be >=' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                max_available_fit_range = min(max_available_fit_range,max_fit_range);
                
                index_range_std_over_mean_check = chm.range >= min_range_std_over_mean & chm.range<=max_available_fit_range;
                
                % std_over_mean = std(log10(abs(RCS(index_range_std_over_mean_check,index_time))),0,2)./mean(log10(abs(RCS(index_range_std_over_mean_check,index_time))),2);
                
                std_over_mean = zeros(length(chm.range(index_range_std_over_mean_check)),1);
                t_sliding_variance = chm.time(index_time);
                for k=1:length(t_sliding_variance)
                    if(t_sliding_variance(k)+dt_sliding_variance /(24*60)<=t_sliding_variance(end))
                        index_time_sliding_variance = chm.time>=t_sliding_variance(k) & chm.time<=t_sliding_variance(k)+dt_sliding_variance/(24*60);
                        %                         std_over_mean_tmp = std(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./mean(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance))),2);
                        
                        std_over_mean_tmp = std(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./...
                            median(log10(abs(RCS(index_range_std_over_mean_check,index_time))),2);
                        std_over_mean = max(std_over_mean,std_over_mean_tmp);
                    end
                end
                
                i_first_over_thresh_std_over_mean = find(std_over_mean>=max_std_over_mean,1,'first');
                if ~isempty(i_first_over_thresh_std_over_mean)
                    rtmp = chm.range(index_range_std_over_mean_check);
                    max_available_fit_range = rtmp(i_first_over_thresh_std_over_mean);
                end
                
                if max_available_fit_range < min_fit_range+min_fit_length
%                         figure;plot(chm.range(index_range_std_over_mean_check),std_over_mean,'.-b');
%                         line(xlim,ones(2,1)*max_std_over_mean,'Color','k','LineStyle','--');
%                         ylim([0 2*max_std_over_mean]);
                    %     error('low std over mean');
                    error_flag(it) = 5;
                    error_string{it} = [error_string{it},['low ' num2str(dt_sliding_variance) '-min sliding std over mean : >= ' num2str(max_std_over_mean) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low ' num2str(dt_sliding_variance) '-min sliding std over mean : >= ' num2str(max_std_over_mean) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                % min_cloud_range = 1200;
                % min_mxd_range = 1200;
                % if any(chm.cbh(1,index_time)-chm.cho >= 0 & chm.cbh(1,index_time)-chm.cho <= min_cloud_range)
                %     error('at least one cbh too low for the prescribed time interval');
                % end
                % if any(chm.mxd(index_time)-chm.cho >= 0 & chm.mxd(index_time)-chm.cho <= min_mxd_range)
                %     error('at least one mxd too low for the prescribed time interval');
                % end
                
                index_range_for_grad = find(chm.range>=min_range_std_over_mean & chm.range <= max_available_fit_range);
                lRCS = log10(abs(RCS(index_range_for_grad,index_time)));
                r = chm.range(index_range_for_grad);
                t = chm.time(index_time);
                
                if length(r) < 3 || length(t) < 3
                    continue;
                end
                
                gradY = conv2(lRCS,[1 2 1;0 0 0;-1 -2 -1],'same');
                relgradY = abs(gradY)./abs(lRCS);
                rconv = r(2:end-1);
                relgradYmax = max(relgradY(2:end-1,2:end-1),[],2);
                i_first_over_thresh_relgradY = find(rconv>=first_range_gradY & relgradYmax>=max_relgrad,1,'first');
                if ~isempty(i_first_over_thresh_relgradY)
                    max_available_fit_range = rconv(i_first_over_thresh_relgradY);
                end
                
                if max_available_fit_range < min_fit_range+min_fit_length
                    %     figure;plot(rconv,relgradYmax,'.-b');
                    %     line(xlim,ones(2,1)*max_relgrad,'Color','k','LineStyle','--');
                    %     ylim([0 2*max_relgrad]);
                    %     error('low grad Y');
                    error_flag(it) = 6;
                    error_string{it} = [error_string{it},['low rel grad Y : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low rel grad Y : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                gradX = conv2(lRCS,[1 0 -1;2 0 -2;1 0 -1],'same');
                relgradX = abs(gradX)./abs(lRCS);
                relgradXmax = max(relgradX(2:end-1,2:end-1),[],2);
                i_first_over_thresh_relgradX = find(rconv>=min_range_std_over_mean & relgradXmax>=max_relgrad,1,'first');
                if ~isempty(i_first_over_thresh_relgradX)
                    max_available_fit_range = min(max_available_fit_range,rconv(i_first_over_thresh_relgradX));
                end
                
                if max_available_fit_range < min_fit_range+min_fit_length
                    %     figure;plot(rconv,relgradXmax,'.-b');
                    %     line(xlim,ones(2,1)*max_relgrad,'Color','k','LineStyle','--');
                    %     ylim([0 2*max_relgrad]);
                    %     error('low grad X');
                    error_flag(it) = 7;
                    error_string{it} = [error_string{it},['low rel grad X : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low rel grad X : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                relgradmagn = sqrt(gradX.^2+gradY.^2)./abs(lRCS);
                relgradmagn_sub = relgradmagn(2:end-1,2:end-1);
                irgradmagn = find(rconv>=first_range_gradY & rconv <= max_available_fit_range);
                for k=length(irgradmagn):-1:1
                    %     if nanmax(nanmax(relgradmagn_sub(irgradmagn(1):irgradmagn(k),:))) <= max_relgrad
                    if nanmax(nanmax(relgradmagn_sub(irgradmagn(1):irgradmagn(k),:))) <= max_relgrad && nanmean(nanmean(relgradmagn_sub(irgradmagn(1):irgradmagn(k),:))) <= max_relgrad_mean
                        max_available_fit_range = rconv(irgradmagn(k));
                        break;
                    end
                end
                
                if max_available_fit_range < min_fit_range+min_fit_length
                    error_flag(it) = 8;
                    error_string{it} = [error_string{it},['low grad magn : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning(['low grad magn : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm']);
                    continue;
                end
                
                
                % lRCStmp = lRCS;
                % rtmp = r;
                % ttmp = t;
                % avg_sz = 2;
                % skip_t = false;
                % while(size(lRCStmp,1)>=3 && size(lRCStmp,2)>=3)
                %     relgradX = abs(conv2(lRCStmp,[1 0 -1;2 0 -2;1 0 -1],'same')./abs(lRCStmp));
                % %     relgradY = abs(conv2(lRCStmp,[1 2 1;0 0 0;-1 -2 -1],'same')./abs(lRCStmp));
                %     rconv = rtmp(2:end-1);
                %
                %     relgradXmax = max(relgradX(2:end-1,2:end-1),[],2);
                %     i_first_over_thresh_relgradX = find(rconv>=min_range_std_over_mean & relgradXmax>=max_relgrad,1,'first');
                %     if ~isempty(i_first_over_thresh_relgradX)
                %         max_available_fit_range = min(max_available_fit_range,rconv(i_first_over_thresh_relgradX));
                %     end
                %
                %     if max_available_fit_range < min_fit_range+min_fit_length
                %         error_flag = 'low grad X';
                % %         error('low grad X');
                %         warning('low grad X');
                %         skip_t = true;
                %         break;
                %     end
                %
                % %     figure;surf(ttmp(2:end-1)-datenum(2000,1,1)+1,rtmp(2:end-1),abs(relgradX(2:end-1,2:end-1)));datetick;colorbar;view(0,90);caxis([0 0.05]);title('relgradX');
                % %     figure;surf(ttmp(2:end-1)-datenum(2000,1,1)+1,rtmp(2:end-1),abs(relgradY(2:end-1,2:end-1)));datetick;colorbar;view(0,90);caxis([0 0.05]);title('relgradY');
                %
                % %     lRCStmp = average_bin(lRCS,avg_sz,false);
                %     lRCStmp = average_bin(lRCS',avg_sz,false);
                %     lRCStmp = lRCStmp';
                % %     rtmp = r(avg_sz:avg_sz:end);
                %     ttmp = t(avg_sz:avg_sz:end);
                %     avg_sz = avg_sz*2;
                % end
                %
                % if(skip_t)
                %     continue;
                % end
                
                
                % lRCSplot = log10(abs(RCS(:,index_time)));
                % relgradX=abs(conv2(lRCSplot,[1 0 -1;2 0 -2;1 0 -1],'same')./abs(lRCSplot));
                % relgradY=abs(conv2(lRCSplot,[1 2 1;0 0 0;-1 -2 -1],'same')./abs(lRCSplot));
                % figure;surf(t(2:end-1)-datenum(2000,1,1)+1,chm.range(2:end-1),abs(relgradX(2:end-1,2:end-1)));datetick;colorbar;view(0,90);caxis([0 0.05]);title('relgradX');ylim([0 1200]);
                % figure;surf(t(2:end-1)-datenum(2000,1,1)+1,chm.range(2:end-1),abs(relgradY(2:end-1,2:end-1)));datetick;colorbar;view(0,90);caxis([0 0.05]);title('relgradY');ylim([0 1200]);
                
                
                error_string{it} = [error_string{it},['the time interval passed the prefit test, and max available fit range is ' num2str(max_available_fit_range) 'm' '\n']];
                
                candidates = {};
                n_tryouts = 0;
                n_polyfit_params_passed = 0;
                
                for fit_length = min_fit_length:d_fit_length:max_fit_length
                    for fit_begin = min_fit_range:d_fit_range:max_available_fit_range
                        if fit_begin+fit_length <= max_available_fit_range
                            
                            n_tryouts = n_tryouts + 1;
                            
                            index_range = chm.range>=fit_begin & chm.range<=fit_begin+fit_length;
                            
                            p = polyfit(chm.range(index_range),mean(log10(abs(RCS(index_range,index_time))),2),1);
                            if p(1)>=min_expected_slope && p(1)<=max_expected_slope && p(2) >= min_expected_zero_fit_value && p(2) <= max_expected_zero_fit_value
                                
                                n_polyfit_params_passed = n_polyfit_params_passed + 1;
                                
                                resid = sqrt(1/nansum(index_range)*nansum((mean(log10(abs(RCS(index_range,index_time))),2)-polyval(p,chm.range(index_range))).^2,1));
                                
                                index_whole_zone = chm.range >= min_fit_range & chm.range <= fit_begin+fit_length;
%                                 p_whole_zone = polyfit(chm.range(index_whole_zone),mean(log10(abs(RCS(index_whole_zone,index_time))),2),1);
%                                 resid_whole_zone = sqrt(1/nansum(index_whole_zone)*nansum((polyval(p,chm.range(index_whole_zone))-polyval(p_whole_zone,chm.range(index_whole_zone))).^2,1));
                                resid_whole_zone = nanmax(abs(mean(log10(abs(RCS(index_whole_zone,index_time))),2)-polyval(p,chm.range(index_whole_zone)))./abs(polyval(p,chm.range(index_whole_zone))));
                                
                                if resid<thresh_resid_rel*nanmean(polyval(p,chm.range(index_range)))
%                                 if resid<thresh_resid
                                    
                                    if resid_whole_zone<thresh_resid_whole_zone
                                        candidates{end+1} = [fit_begin;fit_length;resid;resid_whole_zone;p(1);p(2);nanmean(polyval(p,chm.range(index_range)))];
                                    else
                                        error_string{it} = [error_string{it},['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: whole zone (' num2str(min_fit_range) 'm-' num2str(fit_begin+fit_length) 'm) vs fit zone resid=' num2str(resid) ', should be < ' num2str(thresh_resid_whole_zone) '\n']];
                                        %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[fit_begin fit_begin+fit_length],[],11};
%                                         warning(['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: whole zone (' num2str(min_fit_range) 'm-' num2str(fit_begin+fit_length) 'm) vs fit zone resid=' num2str(resid) ', should be < ' num2str(thresh_resid)]);
                                    end
                                else
                                    error_string{it} = [error_string{it},['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: resid=' num2str(resid) ', should be < ' num2str(thresh_resid_rel*nanmean(polyval(p,chm.range(index_range)))) '\n']];
                                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[fit_begin fit_begin+fit_length],[],10};
%                                     warning(['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: resid=' num2str(resid) ', should be < ' num2str(thresh_resid)]);
                                end
                                
                            else
                                error_string{it} = [error_string{it},['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: p=(' num2str(p(1)/1e-6) '*1e-6,' num2str(p(2)) '), should be ' num2str(min_expected_slope/1e-6) '*1e-6<=p1<=' num2str(max_expected_slope/1e-6) '*1e-6 & ' num2str(min_expected_zero_fit_value) '<=p2<=' num2str(max_expected_zero_fit_value) '\n']];
                                %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[fit_begin fit_begin+fit_length],[],9};
%                                 warning(['fit ' num2str(fit_begin) 'm-' num2str(fit_begin+fit_length) 'm failed: p=(' num2str(p(1)) ',' num2str(p(2)) '), should be ' num2str(min_expected_slope) '<=p1<=' num2str(max_expected_slope) ' & ' num2str(min_expected_zero_fit_value) '<=p2<=' num2str(max_expected_zero_fit_value)]);
                            end
                        end
                    end
                end
                candidates_matrix = cell2mat(candidates);
                ovp_fc_ok_old = ovp_fc_ok;
                ovp_fc_ok_time_interval_old = ovp_fc_ok_time_interval;
                ovp_fc_ok_range_interval_old = ovp_fc_ok_range_interval;
                if(~isempty(candidates_matrix))
                    
                    error_string{it} = [error_string{it},[num2str(size(candidates_matrix,2)) ' range intervals passed the fit test:\n']];
                    for iit = 1:size(candidates_matrix,2)
                        error_string{it} = [error_string{it},[num2str(candidates_matrix(1,iit)) 'm-' num2str(candidates_matrix(1,iit)+candidates_matrix(2,iit)) 'm' ' (' num2str(candidates_matrix(2,iit)) 'm)' ' (alpha=' num2str(candidates_matrix(5,iit)/(-2*1/log(10)*1e-6)) ' Mm^-1, p2=' num2str(candidates_matrix(6,iit)) ', pmean=' num2str(candidates_matrix(7,iit)) ', resid=' num2str(candidates_matrix(3,iit)) ', resid_whole_zone=' num2str(candidates_matrix(4,iit)) ')' '\n']];
                    end

                    ind_ok = [];
%                     ind_results_ok = [];
                    for l=1:size(candidates_matrix,2)
                        range_interval = {candidates_matrix(1,l),candidates_matrix(1,l)+candidates_matrix(2,l)};
                        index_range = chm.range>=range_interval{1} & chm.range<=range_interval{2};
                        result = polyfit(chm.range(index_range),mean(log10(abs(RCS(index_range,index_time))),2),1);
                        dif = mean(log10(abs(RCS(:,index_time))),2)-polyval(result,chm.range);
                        overlap_corr_factor = 10.^(dif);
                        range_stop_correction = range_interval{2};
                        overlap_corr_factor(chm.range>=range_stop_correction)=1;
                        ov_ref = ov;
%                         ov_ref(ov_ref>1) = 1;
                        ovp_fc = ov_ref.*overlap_corr_factor;
                        
                        [max_val_ovp_fc,ind_max_val_ovp_fc] = nanmax(ovp_fc);
                        if(nanmax(ovp_fc))> max_overlap_value*max(ov_ref);
                            error_string{it} = [error_string{it},[num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm' ' (' num2str(candidates_matrix(2,l)) 'm)' ' too high overlap fct value: max = ' num2str(nanmax(ovp_fc)) ' @r=' num2str(chm.range(ind_max_val_ovp_fc)) 'm' ' instead of ' num2str(max_overlap_value) '; Luff ovlap fct is one @r=' num2str(chm.range(find(ov_ref>=1,1,'first'))) 'm' '\n']];
                            %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[candidates_matrix(1,l) candidates_matrix(1,l)+candidates_matrix(2,l)],ovp_fc,13};
%                             warning([num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm: ' 'too high overlap fct value : max = ' num2str(nanmax(ovp_fc)) ' instead of ' num2str(max_overlap_value)]);
                            continue;
                        end
                        
                        
                        index_overlap_valid = chm.range>=min_overlap_valid;
                        r_overlap_valid = chm.range(index_overlap_valid);
                        
                        [valmax,indmax_overlap_valid_rel_error] = nanmax(abs(ov_ref(index_overlap_valid)-ovp_fc(index_overlap_valid))./abs(ov_ref(index_overlap_valid)));
                        
                        
                        if valmax < thresh_overlap_valid_rel_error
                            
                            index_range_for_grad = chm.range>=min_range_std_over_mean & chm.range <= range_interval{2};
                            lRCS_corr = log10(abs(RCS(index_range_for_grad,index_time))./repmat(overlap_corr_factor(index_range_for_grad),1,length(chm.time(index_time))));
                            
                            gradX = conv2(lRCS_corr,[1 0 -1;2 0 -2;1 0 -1],'same');
                            gradY = conv2(lRCS_corr,[1 2 1;0 0 0;-1 -2 -1],'same');
                            relgradmagn = sqrt(gradX.^2+gradY.^2)./abs(lRCS_corr);
                            
                            rtmp = chm.range(index_range_for_grad);rtmp = rtmp(2:end-1);
                            ttmp = chm.time(index_time);ttmp = ttmp(2:end-1);
                            [row,col] = find(relgradmagn(2:end-1,2:end-1)==nanmax(nanmax(relgradmagn(2:end-1,2:end-1))),1,'first');
                            
                            %             if nanmax(nanmax(relgradmagn(2:end-1,2:end-1))) <= max_relgrad
                            %             if round(100*nanmax(nanmax(relgradmagn(2:end-1,2:end-1))))/100 <= max_relgrad && round(100*nanmean(nanmean(relgradmagn(2:end-1,2:end-1))))/100 <= max_relgrad_magn
                            if nanmax(nanmax(relgradmagn(2:end-1,2:end-1))) <= max_relgrad && nanmean(nanmean(relgradmagn(2:end-1,2:end-1))) <= max_relgrad_mean
                                [~,slope,~] = sgolay_smooth(chm.range(1:167),ovp_fc(1:167),sgolay_width,sgolay_ord);
                                r_min_slope = chm.range(1:167);
                                index_range_stop_correction = find(chm.range(1:167)<=range_stop_correction,1,'last');
                                [val_min_slope,index_min_slope] = nanmin(slope);
                                if val_min_slope >= min_slope || index_min_slope>index_range_stop_correction
                                    
%                                     if (nansum(ovp_fc-ov) > 0 && ~overestimated) || (nansum(ovp_fc-ov) < 0 && overestimated)
                                        ind_ok(end+1) = l;
                                        ovp_fc_ok{end+1} = ovp_fc;
                                        ovp_fc_ok_time_interval{end+1} = [fit_time fit_time+time_interval_length/(24*60)];
                                        ovp_fc_ok_range_interval{end+1} = [range_interval{1} range_interval{2}];
                                        %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[candidates_matrix(1,l) candidates_matrix(1,l)+candidates_matrix(2,l)],ovp_fc,0};
%                                         ind_results_ok(end+1) = length(%results);

%                                     else
%                                         error_string{it} = [error_string{it},[num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm' ' (' num2str(candidates_matrix(2,l)) 'm)' ' overlap overestimated/understated failed: sum(ovlap-ov)=' num2str(nansum(ovp_fc-ov)) ', sign should be' num2str(2*double(~overestimated)-1) '\n']];
%                                     end
                                    
                                else
                                    error_string{it} = [error_string{it},[num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm' ' (' num2str(candidates_matrix(2,l)) 'm)' ' slope of overlap fct failed: min=' num2str(nanmin(slope)) ' @r=' num2str(r_min_slope(index_min_slope)) 'm' ' should be >=' num2str(min_slope) '\n']];
                                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[candidates_matrix(1,l) candidates_matrix(1,l)+candidates_matrix(2,l)],ovp_fc,16};
%                                     warning([num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm: ' 'slope failed']);
                                end
                            else
                                error_string{it} = [error_string{it},[num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm' ' (' num2str(candidates_matrix(2,l)) 'm)' ' post grad magn failed: max = ' num2str(nanmax(nanmax(relgradmagn(2:end-1,2:end-1)))) ' @r=' num2str(rtmp(row)) 'm @t=' datestr(ttmp(col),'HH:MM:SS') ' (should be <=' num2str(max_relgrad) ')', ', mean = ' num2str(nanmean(nanmean(relgradmagn(2:end-1,2:end-1)))) ' (should be <=' num2str(max_relgrad_mean) ')' '\n']];
                                %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[candidates_matrix(1,l) candidates_matrix(1,l)+candidates_matrix(2,l)],ovp_fc,15};
%                                 warning([num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm: ' 'post grad magn failed: max = ' num2str(nanmax(nanmax(relgradmagn(2:end-1,2:end-1)))) ' (should be <=' num2str(max_relgrad) ')', ', mean = ' num2str(nanmean(nanmean(relgradmagn(2:end-1,2:end-1)))) ' (should be <=' num2str(max_relgrad_mean) ')']);
                            end
                            
                            %             r = chm.range(index_range_for_grad);
                            %             t = chm.time(index_time);
                            %
                            %             lRCStmp = lRCS_corr;
                            %             rtmp = r;
                            %             ttmp = t;
                            %             avg_sz = 1;
                            %             skip_t = false;
                            %             sz1 = size(lRCS_corr,1);sz2 = size(lRCS_corr,2);
                            % %             sz1 = 3;sz2 = 3;
                            %             while(size(lRCStmp,1)>=sz1 && size(lRCStmp,2)>=sz2)
                            %                 gradX = conv2(lRCStmp,[1 0 -1;2 0 -2;1 0 -1],'same');
                            %                 gradY = conv2(lRCStmp,[1 2 1;0 0 0;-1 -2 -1],'same');
                            %                 relgradmagn = sqrt(gradX.^2+gradY.^2)./abs(lRCStmp);
                            % %                 relgradmagn = sqrt(gradY.^2)./abs(lRCStmp);
                            % %                 relgradmagn = sqrt(gradX.^2+gradY.^2);
                            %                 rconv = rtmp(2:end-1);
                            %
                            % %                 if nanmax(nanmax(relgradmagn(2:end-1,2:end-1))) > max_relgrad
                            %                 if nanmean(nanmean(relgradmagn(2:end-1,2:end-1))) > max_relgrad/5
                            %                     error_flag = 'too high grad magn';
                            %                     warning('too high grad magn');
                            %                     skip_t = true;
                            %                     break;
                            %                 end
                            %
                            %                 lRCStmp = average_bin(lRCS_corr,avg_sz,false);
                            % %                 lRCStmp = average_bin(lRCS_corr',avg_sz,false);
                            %                 lRCStmp = average_bin(lRCStmp',avg_sz,false);
                            %                 lRCStmp = lRCStmp';
                            %                 rtmp = r(avg_sz:avg_sz:end);
                            %                 ttmp = t(avg_sz:avg_sz:end);
                            %                 avg_sz = avg_sz*2;
                            %             end
                            %             if(~skip_t)
                            %                 ind_ok(end+1) = l;
                            %                 ovp_fc_ok{end+1} = ovp_fc;
                            %             end
                            
                        else
                            error_string{it} = [error_string{it},[num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm' ' (' num2str(candidates_matrix(2,l)) 'm)' ' overlap valid test failed: max rel error @r>=' num2str(min_overlap_valid) 'm is ' num2str(valmax) ' @r=' num2str(r_overlap_valid(indmax_overlap_valid_rel_error)) 'm' ' instead of ' num2str(thresh_overlap_valid_rel_error) '\n']];
                            %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[candidates_matrix(1,l) candidates_matrix(1,l)+candidates_matrix(2,l)],ovp_fc,14};
%                             warning([num2str(candidates_matrix(1,l)) 'm-' num2str(candidates_matrix(1,l)+candidates_matrix(2,l)) 'm: ' 'overlap valid test failed: max rel error @r>=' num2str(min_overlap_valid) 'm is ' num2str(valmax) ' instead of ' num2str(thresh_overlap_valid_rel_error)]);
                        end
                        
                    end
                    
                    if length(ovp_fc_ok)>length(ovp_fc_ok_old)
                        tmp = cell2mat(ovp_fc_ok);
                        ok_candidates{it} = tmp(:,length(ovp_fc_ok_old)+1:end);
                    else
                        ok_candidates{it} = [];
                    end
                    
                    
                    if(length(ind_ok)<min_nb_ok_candidates)
                       if(isempty(ind_ok))
                            error_flag(it) = 17;
                            error_string{it} = [error_string{it},['no candidate range intervals passed the post-fit test']];
                            %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
                       else
                            ovp_fc_ok = ovp_fc_ok_old;
                            ovp_fc_ok_time_interval = ovp_fc_ok_time_interval_old;
                            ovp_fc_ok_range_interval = ovp_fc_ok_range_interval_old;
                            error_flag(it) = 18;
                            error_string{it} = [error_string{it},['not enough candidate range intervals passed the post-fit test: ' num2str(length(ind_ok)) ' instead of ' num2str(min_nb_ok_candidates)]];
%                             for j=1:length(ind_results_ok)
%                                res_tmp = results{ind_results_ok(j)};
%                                res_tmp{4} = error_flag(it);
%                                results{ind_results_ok(j)} = res_tmp;
%                             end
%                             warning(['suitable range not found: not enough ok candidates: ' num2str(length(ind_ok)) ' instead of ' num2str(min_nb_ok_candidates)]);
                       end
                    else
%                         [valmin,indmin] = nanmin(candidates_matrix(1,ind_ok),[],2);
%                         ind_find = find(candidates_matrix(1,ind_ok)==valmin);
%                         range_interval = {candidates_matrix(1,ind_ok(ind_find(end))),candidates_matrix(1,ind_ok(ind_find(end)))+candidates_matrix(2,ind_ok(ind_find(end)))};
                        
                        error_string{it} = [error_string{it},[num2str(length(ind_ok)) ' range intervals passed the post-fit test:\n']];
                        for iit = 1:length(ind_ok)-1
                            error_string{it} = [error_string{it},[num2str(candidates_matrix(1,ind_ok(iit))) 'm-' num2str(candidates_matrix(1,ind_ok(iit))+candidates_matrix(2,ind_ok(iit))) 'm' '\n']];
                        end
                        error_string{it} = [error_string{it},[num2str(candidates_matrix(1,ind_ok(end))) 'm-' num2str(candidates_matrix(1,ind_ok(end))+candidates_matrix(2,ind_ok(end))) 'm']];
                    end
                                       
                    
%                     if length(ind_ok) >= min_nb_ok_candidates
%                         [valmin,indmin] = nanmin(candidates_matrix(1,ind_ok),[],2);
%                         ind_find = find(candidates_matrix(1,ind_ok)==valmin);
%                         range_interval = {candidates_matrix(1,ind_ok(ind_find(end))),candidates_matrix(1,ind_ok(ind_find(end)))+candidates_matrix(2,ind_ok(ind_find(end)))};
%                     else
%                         ovp_fc_ok = ovp_fc_ok_old;
%                         warning(['suitable range not found: not enough ok candidates: ' num2str(length(ind_ok)) ' instead of ' num2str(min_nb_ok_candidates)]);
%                     end
                else
                    error_flag(it) = 12;
                    error_string{it} = [error_string{it},['no candidates found because no range intervals passed the fit test']];
                    %results{end+1} = {[fit_time,fit_time+time_interval_length/(24*60)],[],[],error_flag(it)};
%                     warning('no candidates found because no range intervals passed the fit test');
                end
                
                % n_polyfit_params_passed / n_tryouts * 100
                
            end
        end     
        
    end
    
    
    ovp_fc_ok_mat = cell2mat(ovp_fc_ok);
    ovp_fc_ok_mat_time = NaN(2,length(ovp_fc_ok_time_interval));
    ovp_fc_ok_mat_range = NaN(2,length(ovp_fc_ok_time_interval));
    for j=1:length(ovp_fc_ok_time_interval)
        ovp_fc_ok_mat_time(:,j) = ovp_fc_ok_time_interval{j}';
        ovp_fc_ok_mat_range(:,j) = ovp_fc_ok_range_interval{j}';
    end
    
    if(length(ovp_fc_ok)<min_nb_samples)
%           figure
%             flag_list={1,'no data in the prescribed time interval';...
%                 2,'at least one sci~=0 for the prescribed time interval';...
%                 3,['low cloud , should be >=' num2str(min_fit_range+min_fit_length) 'm'];...
%                 4,['low mxd, should be >=' num2str(min_fit_range+min_fit_length) 'm'];...
%                 5,['low ' num2str(dt_sliding_variance) '-min sliding std over mean : >= ' num2str(max_std_over_mean) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm'];...
%                 6,['low rel grad Y : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm'];...
%                 7,['low rel grad X : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(min_range_std_over_mean) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm'];...
%                 8,['low grad magn : >= ' num2str(max_relgrad) ' @ ' num2str(max_available_fit_range) ' m (search started at ' num2str(first_range_gradY) 'm, should reach ' num2str(min_fit_range+min_fit_length) 'm'];...
%                 17,'no candidate range intervals passed the post-fit test';...
%                 18,['not enough candidate range intervals passed the post-fit test (min: ' num2str(min_nb_ok_candidates) ')'];...
%                 };
%             plot(1:length(error_flag),error_flag)
%             ylim([0 20])
%             set(gca,'ytick',cell2mat(flag_list(:,1)),'yticklabel',flag_list(:,2))
%             grid on
        if(isempty(ovp_fc_ok))
            warning('no samples found');
          
            return;
        else
            warning(['unsufficient amount of samples found: ' num2str(length(ovp_fc_ok)) ', should be ' num2str(min_nb_samples)]);
            return;
        end
    else
        disp([num2str(length(ovp_fc_ok)) ' samples found. Looking for good samples.']);
    end
   
    ind_ovp_fc_ok_good = 1:length(ovp_fc_ok);
    if(length(ovp_fc_ok)>=min_nb_samples_for_skipping_good_test)
        ovp_fc_ok_good = ovp_fc_ok;
        ovp_fc_ok_time_interval_good = ovp_fc_ok_time_interval;
        ovp_fc_ok_range_interval_good = ovp_fc_ok_range_interval;
    else
        
        t1_value = zeros(1,length(ovp_fc_ok));
        t2_value = zeros(1,length(ovp_fc_ok));
        rmax_value = zeros(1,length(ovp_fc_ok));
        for l=1:length(ovp_fc_ok)
            time_interval = ovp_fc_ok_time_interval{l};
            range_interval = ovp_fc_ok_range_interval{l};
            rmax = range_interval(2);
            for kk=1:length(ovp_fc_ok)
                time_interval_kk = ovp_fc_ok_time_interval{kk};
                range_interval_kk = ovp_fc_ok_range_interval{kk};
                if isequal(time_interval_kk,time_interval)
                    rmax = max(rmax,range_interval_kk(2));
                end
            end
            rmax_value(l) = rmax;
            t1_value(l) = time_interval(1);
            t2_value(l) = time_interval(2);
        end

        [t1_value_unique,i_t_values,~] = unique(t1_value);
        t2_value_unique =  t2_value(i_t_values);
        rmax_value_unique = rmax_value(i_t_values);


       nb_ok = zeros(1,length(ovp_fc_ok));
       for l=1:length(ovp_fc_ok)
    %        disp([num2str(l) ' of ' num2str(length(ovp_fc_ok))])
           ovp_fc_ok_l = ovp_fc_ok{l};
           for ll=1:length(t1_value_unique)
    %             disp([num2str(ll) ' of ' num2str(length(t1_value_unique))])
                time_interval = [t1_value_unique(ll) t2_value_unique(ll)];

                rmax = rmax_value_unique(ll);

                index_time = chm.time>=time_interval(1) & chm.time<=time_interval(2);
                index_range_for_grad = chm.range>=min_range_std_over_mean & chm.range <= rmax;

                RCS_corr = RCS.*repmat(ov,1,length(chm.time))./repmat(ovp_fc_ok_l,1,length(chm.time));
                lRCS_corr = log10(abs(RCS_corr(index_range_for_grad,index_time)));

                gradX = conv2(lRCS_corr,[1 0 -1;2 0 -2;1 0 -1],'same');
                gradY = conv2(lRCS_corr,[1 2 1;0 0 0;-1 -2 -1],'same');
                relgradmagn = sqrt(gradX.^2+gradY.^2)./abs(lRCS_corr);

                rtmp = chm.range(index_range_for_grad);rtmp = rtmp(2:end-1);
                ttmp = chm.time(index_time);ttmp = ttmp(2:end-1);
    %             [row,col] = find(relgradmagn(2:end-1,2:end-1)==nanmax(nanmax(relgradmagn(2:end-1,2:end-1))),1,'first');

                if nanmax(nanmax(relgradmagn(2:end-1,2:end-1))) <= max_relgrad && nanmean(nanmean(relgradmagn(2:end-1,2:end-1))) <= max_relgrad_mean
                    nb_ok(l) = nb_ok(l)+1;
                else
                    break;
                end
           end
       end

       ind_good = (nb_ok == length(t1_value_unique));

    %    ovp_fc_ok_good = ovp_fc_ok;
       ovp_fc_ok_good = ovp_fc_ok(ind_good);
       ovp_fc_ok_time_interval_good = ovp_fc_ok_time_interval(ind_good);
       ovp_fc_ok_range_interval_good = ovp_fc_ok_range_interval(ind_good);

       disp([num2str(length(ovp_fc_ok_good)) ' good samples found. Looking for std over mean.']);



        t1_value = zeros(1,length(ovp_fc_ok_good));
        t2_value = zeros(1,length(ovp_fc_ok_good));
        rmax_value = zeros(1,length(ovp_fc_ok_good));
        for l=1:length(ovp_fc_ok_good)
            time_interval = ovp_fc_ok_time_interval_good{l};
            range_interval = ovp_fc_ok_range_interval_good{l};
            rmax = range_interval(2);
            for kk=1:length(ovp_fc_ok_good)
                time_interval_kk = ovp_fc_ok_time_interval_good{kk};
                range_interval_kk = ovp_fc_ok_range_interval_good{kk};
                if isequal(time_interval_kk,time_interval)
                    rmax = max(rmax,range_interval_kk(2));
                end
            end
            rmax_value(l) = rmax;
            t1_value(l) = time_interval(1);
            t2_value(l) = time_interval(2);
        end

        [t1_value_unique,i_t_values,~] = unique(t1_value);
        t2_value_unique =  t2_value(i_t_values);
        rmax_value_unique = rmax_value(i_t_values);


        nb_ok_good = zeros(1,length(ovp_fc_ok_good));
        for l=1:length(ovp_fc_ok_good)
    %        disp([num2str(l) ' of ' num2str(length(ovp_fc_ok_good))])
           ovp_fc_ok_l = ovp_fc_ok_good{l};
           for ll=1:length(t1_value_unique)
                time_interval = [t1_value_unique(ll) t2_value_unique(ll)];

                rmax = rmax_value_unique(ll);

                index_time = chm.time>=time_interval(1) & chm.time<=time_interval(2);
                index_range_std_over_mean_check = chm.range>=min_range_std_over_mean & chm.range <= rmax;

                RCS_corr = RCS.*repmat(ov,1,length(chm.time))./repmat(ovp_fc_ok_l,1,length(chm.time));

                std_over_mean = zeros(length(chm.range(index_range_std_over_mean_check)),1);
                t_sliding_variance = chm.time(index_time);
                for k=1:length(t_sliding_variance)
                    if(t_sliding_variance(k)+dt_sliding_variance /(24*60)<=t_sliding_variance(end))
                        index_time_sliding_variance = chm.time>=t_sliding_variance(k) & chm.time<=t_sliding_variance(k)+dt_sliding_variance/(24*60);
        %                         std_over_mean_tmp = std(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./mean(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),2);
                        std_over_mean_tmp = std(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./median(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time))),2);
                        std_over_mean = max(std_over_mean,std_over_mean_tmp);
                        if nanmax(std_over_mean) > max_std_over_mean
                            break;
                        end
                    end
                end
                if nanmax(std_over_mean) <= max_std_over_mean
                    nb_ok_good(l) = nb_ok_good(l)+1;
                else
                    break;
                end

           end
        end
        ind_good_good = (nb_ok_good == length(t1_value_unique));
        ind_ovp_fc_ok_good_tmp = ind_ovp_fc_ok_good(ind_good);
        ind_ovp_fc_ok_good = ind_ovp_fc_ok_good_tmp(ind_good_good);
        ovp_fc_ok_good = ovp_fc_ok_good(ind_good_good);
        ovp_fc_ok_time_interval_good = ovp_fc_ok_time_interval_good(ind_good_good);
        ovp_fc_ok_range_interval_good = ovp_fc_ok_range_interval_good(ind_good_good);

    %     
    %     
    %     
    %    h = waitbar(0,'please wait...');
    %    nb_ok = zeros(1,length(ovp_fc_ok));
    %    for l=1:length(ovp_fc_ok)
    % %        disp(l)
    %        waitbar(l / length(ovp_fc_ok))
    %        ovp_fc_ok_l = ovp_fc_ok{l};
    %        time_interval_l = ovp_fc_ok_time_interval{l};
    %        for ll=1:length(ovp_fc_ok)
    %             time_interval = ovp_fc_ok_time_interval{ll};
    %             range_interval = ovp_fc_ok_range_interval{ll};
    %             
    % %             if ~isequal(time_interval_l,time_interval)
    % 
    %                 rmax = range_interval(2);
    % %                 for kk=1:length(ovp_fc_ok)
    % %                     time_interval_kk = ovp_fc_ok_time_interval{kk};
    % %                     range_interval_kk = ovp_fc_ok_range_interval{kk};
    % %                     if isequal(time_interval_kk,time_interval)
    % %                        rmax = max(rmax,range_interval_kk(2));
    % %                     end
    % %                 end
    % 
    %                 index_time = chm.time>=time_interval(1) & chm.time<=time_interval(2);
    %                 index_range_for_grad = chm.range>=min_range_std_over_mean & chm.range <= rmax;
    %                 
    %                 lRCS_corr = log10(abs(RCS(index_range_for_grad,index_time)).*repmat(ov(index_range_for_grad),1,length(chm.time(index_time)))./repmat(ovp_fc_ok_l(index_range_for_grad),1,length(chm.time(index_time))));
    % 
    %                 gradX = conv2(lRCS_corr,[1 0 -1;2 0 -2;1 0 -1],'same');
    %                 gradY = conv2(lRCS_corr,[1 2 1;0 0 0;-1 -2 -1],'same');
    %                 relgradmagn = sqrt(gradX.^2+gradY.^2)./abs(lRCS_corr);
    % 
    %                 rtmp = chm.range(index_range_for_grad);rtmp = rtmp(2:end-1);
    %                 ttmp = chm.time(index_time);ttmp = ttmp(2:end-1);
    %     %             [row,col] = find(relgradmagn(2:end-1,2:end-1)==nanmax(nanmax(relgradmagn(2:end-1,2:end-1))),1,'first');
    %                 if nanmax(nanmax(relgradmagn(2:end-1,2:end-1))) <= max_relgrad && nanmean(nanmean(relgradmagn(2:end-1,2:end-1))) <= max_relgrad_mean
    % %                     index_range_std_over_mean_check = index_range_for_grad;
    % %                     std_over_mean = zeros(length(chm.range(index_range_std_over_mean_check)),1);
    % %                     t_sliding_variance = chm.time(index_time);
    % %                     for k=1:length(t_sliding_variance)
    % %                         if(t_sliding_variance(k)+dt_sliding_variance /(24*60)<=t_sliding_variance(end))
    % %                             index_time_sliding_variance = chm.time>=t_sliding_variance(k) & chm.time<=t_sliding_variance(k)+dt_sliding_variance/(24*60);
    % %     %                         std_over_mean_tmp = std(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance)).*repmat(ov(index_range_for_grad),1,length(chm.time(index_time)))./repmat(ovp_fc_ok_l(index_range_for_grad),1,length(chm.time(index_time)))),0,2)./mean(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance)).*repmat(ov(index_range_for_grad),1,length(chm.time(index_time)))./repmat(ovp_fc_ok_l(index_range_for_grad),1,length(chm.time(index_time)))),2);
    % %                             std_over_mean_tmp = std(log10(abs(RCS(index_range_std_over_mean_check,index_time_sliding_variance)).*repmat(ov(index_range_for_grad),1,length(chm.time(index_time)))./repmat(ovp_fc_ok_l(index_range_for_grad),1,length(chm.time(index_time)))),0,2)./mean(log10(abs(RCS(index_range_std_over_mean_check,index_time)).*repmat(ov(index_range_for_grad),1,length(chm.time(index_time)))./repmat(ovp_fc_ok_l(index_range_for_grad),1,length(chm.time(index_time)))),2);
    % %                             std_over_mean = max(std_over_mean,std_over_mean_tmp);
    % %                         end
    % %                     end
    % %                     if nanmax(std_over_mean) <= max_std_over_mean
    %                         nb_ok(l) = nb_ok(l)+1;
    % %                     else
    %                         
    % %                     end
    %                 else
    % 
    %                 end
    %             
    % %             end
    %        end
    %    end
    %    close(h)
    % %    good_samples_proportion = 0.0;
    %    ind_good = (nb_ok >= good_samples_proportion*length(ovp_fc_ok));
    % %    ovp_fc_ok_good = ovp_fc_ok;
    %    ovp_fc_ok_good = ovp_fc_ok(ind_good);
    %    ovp_fc_ok_time_interval_good = ovp_fc_ok_time_interval(ind_good);
    %    ovp_fc_ok_range_interval_good = ovp_fc_ok_range_interval(ind_good);
    %    
    %    disp([num2str(length(ovp_fc_ok_good)) ' good samples found. Looking for std over mean.']);
    %    
    % %    nb_ok_good = zeros(1,length(ovp_fc_ok_good));
    % %    for l=1:length(ovp_fc_ok_good)
    % %        disp(l)
    % %        ovp_fc_ok_l = ovp_fc_ok_good{l};
    % %        time_interval_l = ovp_fc_ok_time_interval_good{l};
    % %        for ll=1:length(ovp_fc_ok_good)
    % %             time_interval = ovp_fc_ok_time_interval_good{ll};
    % %             range_interval = ovp_fc_ok_range_interval_good{ll};
    % %             
    % %             rmax = range_interval(2);
    % % 
    % %             index_time = chm.time>=time_interval(1) & chm.time<=time_interval(2);
    % %             index_range_std_over_mean_check = chm.range>=min_range_std_over_mean & chm.range <= rmax;
    % % 
    % %             RCS_corr = RCS.*repmat(ov,1,length(chm.time))./repmat(ovp_fc_ok_l,1,length(chm.time));
    % % 
    % %             std_over_mean = zeros(length(chm.range(index_range_std_over_mean_check)),1);
    % %             t_sliding_variance = chm.time(index_time);
    % %             for k=1:length(t_sliding_variance)
    % %                 if(t_sliding_variance(k)+dt_sliding_variance /(24*60)<=t_sliding_variance(end))
    % %                     index_time_sliding_variance = chm.time>=t_sliding_variance(k) & chm.time<=t_sliding_variance(k)+dt_sliding_variance/(24*60);
    % % %                         std_over_mean_tmp = std(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./mean(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),2);
    % %                     std_over_mean_tmp = std(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time_sliding_variance))),0,2)./mean(log10(abs(RCS_corr(index_range_std_over_mean_check,index_time))),2);
    % %                     std_over_mean = max(std_over_mean,std_over_mean_tmp);
    % %                 end
    % %             end
    % %             if nanmax(std_over_mean) <= max_std_over_mean
    % %                 nb_ok_good(l) = nb_ok_good(l)+1;
    % %             else
    % % 
    % %             end
    % %             
    % %        end
    % %    end
    % %    good_samples_proportion = 1.0;
    % %    ind_good_good = (nb_ok_good >= good_samples_proportion*length(ovp_fc_ok_good));
    % %    ovp_fc_ok_good = ovp_fc_ok_good(ind_good_good);
    %    
    %    

    end
    
    
    if(length(ovp_fc_ok_good)<min_nb_good_samples)
        if(isempty(ovp_fc_ok_good))
            warning('no good samples found');
            return;
        else
            warning(['unsufficient amount of good samples found: ' num2str(length(ovp_fc_ok_good)) ', should be ' num2str(min_nb_good_samples)]);
            return;
        end
    else
        disp([num2str(length(ovp_fc_ok_good)) ' good samples found. Looking for outliers.']);
    end
    
%     whiskers_length = Inf;
    ovp_fc_ok_good_mat = cell2mat(ovp_fc_ok_good);
    pctile_25 = prctile(ovp_fc_ok_good_mat,25,2);
    pctile_75 = prctile(ovp_fc_ok_good_mat,75,2);
    pctile_50 = prctile(ovp_fc_ok_good_mat,50,2);
    outliers_plus = pctile_50 + whiskers_length*(pctile_75-pctile_25);
    outliers_minus = pctile_50 - whiskers_length*(pctile_75-pctile_25);

    ovp_fc_ok_good_mat_time = NaN(2,length(ovp_fc_ok_time_interval_good));
    ovp_fc_ok_good_mat_range = NaN(2,length(ovp_fc_ok_time_interval_good));
    for j=1:length(ovp_fc_ok_time_interval_good)
        ovp_fc_ok_good_mat_time(:,j) = ovp_fc_ok_time_interval_good{j}';
        ovp_fc_ok_good_mat_range(:,j) = ovp_fc_ok_range_interval_good{j}';
    end
    
    for j=1:size(ovp_fc_ok_good_mat,2)
        if ~any(ovp_fc_ok_good_mat(:,j)>outliers_plus | ovp_fc_ok_good_mat(:,j)<outliers_minus)
            ovp_fc_ok_mat_final(:,end+1) = ovp_fc_ok_good_mat(:,j);
            ovp_fc_ok_mat_time_final(:,end+1) = ovp_fc_ok_good_mat_time(:,j);
            ovp_fc_ok_mat_range_final(:,end+1) = ovp_fc_ok_good_mat_range(:,j);
            ind_ovp_fc_ok_final(:,end+1) = ind_ovp_fc_ok_good(:,j);
        end
    end

    
    if size(ovp_fc_ok_mat_final,2)>=min_nb_good_samples_after_outliers_removal
        
        disp([num2str(size(ovp_fc_ok_mat_final,2)) ' good samples found after outliers removal. Taking the median.']);
        
        ovp_fc_final = nanmedian(ovp_fc_ok_mat_final,2);
%         ovp_fc_final(ovp_fc_final>1) = 1;

    else
        %     error('no samples found after outliers removal');
        %     warning('no samples found after outliers removal');continue;
        disp(['unsufficient amount of good samples found after outlier removal: ' num2str(size(ovp_fc_ok_mat_final,2)) ', should be ' num2str(min_nb_samples)]);return;
    end