%% Choose the station and the day

stn_list = {'pay'                    ,'kse'      ,'SIRTA'    ,'Granada'  ,'Lindenberg','Hohenspeissenberg','Hamburg'  ,'Oslo'     };
%          {{'TUB120011','TUB140007'},'TUB140005','TUB140013','TUB120012','TUB120001' ,'TUB070009'        ,'TUB100011','TUB110019'}
stn_list_short = {'py','ks','st','gr','ln','hs','hm','os'};

do_it_manual = false;
% stn = 'pay';
% date_yyyymmdd = '20140616';
% date_yyyymmdd = '20150305';
stn = 'kse';
% date_yyyymmdd = '20150523';
% date_yyyymmdd = '20141123';
% stn = 'SIRTA';
% date_yyyymmdd = '20150420';
% stn = 'Granada';
% date_yyyymmdd = '20130307';
% stn = 'Lindenberg';
% date_yyyymmdd = '20140915';
% data_yyyymmdd = '20150706';
% stn = 'Hamburg';
% date_yyyymmdd = '20150524';


for dn = datenum(2014,08,29):1:datenum(2015,12,15)
    date_yyyymmdd = datestr(dn,'yyyymmdd');
    
    [chm,chminfo,ov] = load_ceilo_and_overlap_data(stn,date_yyyymmdd);
    
    if isempty(chm)
        continue;
    end
    
    if any(isnan(ov))
%         fid = fopen('TUB120001_20120917_1024.cfg');ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
        disp('using default ovp fc');
        fid = fopen('TUB120011_20121112_1024.cfg');ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
        ov_to_use = cell2mat(ov_cell);
        ov = [];
    else
        ov_to_use = ov;
    end
    RCS = chm.beta_raw;
    
    %----------------------------------------------------------------------
    % [ovp_fc_final,results] = calculate_overlap_automatic(chm,chminfo,RCS,ov_to_use);
    
    %----------------------------------------------------------------------
%     [ovp_fc_ok_mat_final,ovp_fc_ok_mat_time_final,ovp_fc_ok_mat_range_final,ind_ovp_fc_ok_good,ind_ovp_fc_ok_final,ovp_fc_final,error_string,results] = calculate_overlap_automatic_structured(chm,chminfo,RCS,ov_to_use);
%     
%     if ~isempty(results)
%         result = struct;
%         result.time = NaN(length(results),2);
%         result.range = NaN(length(results),2);
%         result.tag =  NaN(length(results),1);
%         for k=1:length(results)
%             result.time(k,:) =  results{k}{1};
%             if ~isempty(results{k}{2})
%                 result.range(k,:) =  results{k}{2};
%             else
%                 result.range(k,:) = [NaN NaN];
%             end
%             result.tag(k) = results{k}{4};
%         end
%         
%         istn = strcmpi(stn,stn_list);
% %         out_dir = ['C:\AllData\SharedData_Maxime\' stn_list_short{istn} '\' date_yyyymmdd(1:4) '\' date_yyyymmdd(5:6) '\'];
%         out_dir = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/corrections_other_sites/' stn_list{istn} '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
%         if ~exist(out_dir,'dir')
%             mkdir(out_dir);
%         end
%         
%         for j=1:length(chminfo.Attributes)
%             if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
%                 device_name = chminfo.Attributes(j).Value;
%             elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
%                 serlom = chminfo.Attributes(j).Value;
%             end
%         end
%         
%         out_file = [out_dir 'result_tag_' device_name serlom '_' date_yyyymmdd '.mat'];
%         save(out_file,'result');
%     end
%     
% %     ind_ok_results = zeros(1,length(results));
% %     for k=1:length(results)
% %         if(results{k}{4}==0)
% %             ind_ok_results(k) = 1;
% %         end
% %     end
% %     results_ok = results(ind_ok_results==1);
% %     
% %     if ~isempty(results_ok)
% %         
% %         result = struct;
% %         result.time_start = NaN(length(results_ok),1);
% %         result.time_end = NaN(length(results_ok),1);
% %         result.range_start = NaN(length(results_ok),1);
% %         result.range_end = NaN(length(results_ok),1);
% %         result.ovp_fc = NaN(1024,length(results_ok));
% %         
% %         for k=1:length(results_ok)
% %             time_int = results_ok{k}{1};
% %             result.time_start(k) = time_int(1);
% %             result.time_end(k) = time_int(2);
% %             range_int = results_ok{k}{2};
% %             result.range_start(k) = range_int(1);
% %             result.range_end(k) = range_int(2);
% %             result.ovp_fc(:,k) = results_ok{k}{3};
% %         end
% %         
% %         result.index_good = ind_ovp_fc_ok_good';
% %         result.index_final = ind_ovp_fc_ok_final';
% %         
% %         istn = strcmpi(stn,stn_list);
% % %         out_dir = ['C:\AllData\SharedData_Maxime\' stn_list_short{istn} '\' date_yyyymmdd(1:4) '\' date_yyyymmdd(5:6) '\'];
% %         out_dir = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/corrections_other_sites/' stn_list{istn} '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
% %         if ~exist(out_dir,'dir')
% %             mkdir(out_dir);
% %         end
% %         
% %         for j=1:length(chminfo.Attributes)
% %             if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
% %                 device_name = chminfo.Attributes(j).Value;
% %             elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
% %                 serlom = chminfo.Attributes(j).Value;
% %             end
% %         end
% %         
% %         out_file = [out_dir 'result_ovp_fc_' device_name serlom '_' date_yyyymmdd '.mat'];
% %         save(out_file,'result');
% %         
% %     end
    
    %----------------------------------------------------------------------
    [ovp_fc_ok_mat,ovp_fc_ok_mat_time,ovp_fc_ok_mat_range,ind_ovp_fc_ok_good,ind_ovp_fc_ok_final,ovp_fc_final,error_string] = calculate_overlap_automatic_structured_2(chm,chminfo,RCS,ov_to_use);
    if ~isempty(ovp_fc_ok_mat)
        
        result = struct;
        result.time_start = NaN(size(ovp_fc_ok_mat,2),1);
        result.time_end = NaN(size(ovp_fc_ok_mat,2),1);
        result.range_start = NaN(size(ovp_fc_ok_mat,2),1);
        result.range_end = NaN(size(ovp_fc_ok_mat,2),1);
        result.ovp_fc = NaN(size(ovp_fc_ok_mat,1),size(ovp_fc_ok_mat,2));
        
        for k=1:size(ovp_fc_ok_mat,2)
            result.time_start(k) = ovp_fc_ok_mat_time(1,k);
            result.time_end(k) = ovp_fc_ok_mat_time(2,k);
            result.range_start(k) = ovp_fc_ok_mat_range(1,k);
            result.range_end(k) = ovp_fc_ok_mat_range(2,k);
            result.ovp_fc(:,k) = ovp_fc_ok_mat(:,k);
        end
        
        result.index_good = ind_ovp_fc_ok_good';
        result.index_final = ind_ovp_fc_ok_final';
        
        ovp_fc_ok_mat_final = ovp_fc_ok_mat(:,ind_ovp_fc_ok_final);
        
        istn = strcmpi(stn,stn_list);
%         out_dir = ['C:\AllData\SharedData_Maxime\' stn_list_short{istn} '\' date_yyyymmdd(1:4) '\' date_yyyymmdd(5:6) '\'];
        out_dir = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/corrections_other_sites/' stn_list{istn} '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end
        
        for j=1:length(chminfo.Attributes)
            if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
                device_name = chminfo.Attributes(j).Value;
            elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
                serlom = chminfo.Attributes(j).Value;
            end
        end
        
        out_file = [out_dir 'result_ovp_fc_' device_name serlom '_' date_yyyymmdd '.mat'];
        save(out_file,'result');
        
    end

    %----------------------------------------------------------------------
    
    if ~isempty(ovp_fc_final)

        istn = strcmpi(stn,stn_list);
%         out_dir = ['C:\AllData\SharedData_Maxime\' stn_list_short{istn} '\' date_yyyymmdd(1:4) '\' date_yyyymmdd(5:6) '\'];
        out_dir = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/corrections_other_sites/' stn_list{istn} '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end

        for j=1:length(chminfo.Attributes)
            if strcmpi(chminfo.Attributes(j).Name,'source') || strcmpi(chminfo.Attributes(j).Name,'device_name')
                device_name = chminfo.Attributes(j).Value;
            elseif strcmpi(chminfo.Attributes(j).Name,'serlom')
                serlom = chminfo.Attributes(j).Value;
            end
        end

        grad_raw = NaN(size(RCS));
        for j=1:length(chm.time)
            grad_raw(:,j) = 1/(2*chm.range_gate)*conv(RCS(:,j),[1 0 -1],'same');
        end
        RCS_corr = RCS.*repmat(ov_to_use,1,length(chm.time))./repmat(ovp_fc_final,1,length(chm.time));
        grad_corr = NaN(size(RCS_corr));

        for j=1:length(chm.time)
            grad_corr(:,j) = 1/(2*chm.range_gate)*conv(RCS_corr(:,j),[1 0 -1],'same');
        end

        file_out = [out_dir 'ovp_fc_final_'  stn '_' device_name serlom '_' date_yyyymmdd '.png'];
        plot_overlap_routine_ov_fct(chm,chminfo,ov,ovp_fc_ok_mat_final,ovp_fc_final,file_out);

        file_out = [out_dir stn '_' device_name serlom '_' date_yyyymmdd '.png'];
        plot_overlap_routine_RCS_GRADRCS(chm,chminfo,RCS,grad_raw,RCS_corr,grad_corr,file_out,result);

        out_file = [out_dir 'ovp_fc_final_' device_name serlom '_' date_yyyymmdd '.mat'];
        save(out_file,'ovp_fc_final');
    
    end
 
end