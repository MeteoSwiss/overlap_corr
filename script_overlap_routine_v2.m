% Script To correct overlap functions
% By Y. Poltera and M. Hervo (Meteoswiss) 2016-2019

clear variables; close all; clc

%% Choose the station and the day

stn = 'kse';

start_time = datenum(2019,3,25);
end_time = datenum(2020,11,1);

%% List Of ceilometers previously studied
stn_list = {'pay'                                ,'kse'      ,'SIRTA'    ,'Granada'  ,'Lindenberg','Hohenspeissenberg','Hamburg'  ,'Oslo'     };

%          {{'TUB120011','TUB140007','TUB140016'},'TUB140005','TUB140013','TUB120012','TUB120001' ,'TUB070009'        ,'TUB100011','TUB110019'}
% stn_list_short = {'py','ks','st','gr','ln','hs','hm','os'};


% stn = 'pay';
% date_yyyymmdd = '20140616';
% date_yyyymmdd = '20150305';
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


%% Loop on time
answer=[];
for dn = start_time:1:end_time
    date_yyyymmdd = datestr(dn,'yyyymmdd');
    
    [chm,chminfo,ov] = load_ceilo_and_overlap_data(stn,date_yyyymmdd);
    
    if isempty(chm)
        continue;
    end
    
    if any(isnan(ov))
        if isempty(answer)
            answer=questdlg({'No reference overlap function found','Do you want to calibrate with TUB120011 overlap function'},...
                'No reference overlap function','Yes','No','No');
            if strcmp(answer,'Yes')
                warning('No overlap function found. Using TUB120011 overlap function')
                fid = fopen('TUB120011_20121112_1024.cfg');ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
                ov_to_use = cell2mat(ov_cell);
                ov = [];
            else
                error('No overlap function found. Ask the corresponding overlap to the manufacturer')
            end
        end
    else
        ov_to_use = ov;
    end
    RCS = chm.beta_raw;
    
    
    %----------------------------------------------------------------------
    %% Calculate overlap correction
    [ovp_fc_ok_mat,ovp_fc_ok_mat_time,ovp_fc_ok_mat_range,ind_ovp_fc_ok_good,ind_ovp_fc_ok_final,ovp_fc_final,error_string] = calculate_overlap_automatic_structured_2(chm,chminfo,RCS,ov_to_use);
    
    
    %% Save results
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
    %% Make and Save plots
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