% Script To correct overlap functions
% V3: adapted for E-PROFILE Network
% By Y. Poltera and M. Hervo (Meteoswiss) 2016-2018

clear variables; close all; clc

%% Choose the station and the day

% stn = '0-20008-0-UGR';
% id = 'A';
stn = '0-20000-0-10393';
id = '0';
% stn = '0-20000-0-06610';
% id = 'A';

start_time = datenum(2016,11,2);
end_time = datenum(2016,12,31);


folder_data = 'E:\E_PROFILE_ALC\L1_FILES\';


foulder_out = 'C:\DATA\Projects\201806 Start Finder Overlap/';


%% Loop on time
answer=[];
for dn = start_time:1:end_time
    date_yyyymmdd = datestr(dn,'yyyymmddHHMMSS');
    date_end_yyyymmdd = datestr(dn+0.999999,'yyyymmddHHMMSS');
    
    L1 =read_L1_EPROFILE_v4(stn,id,date_yyyymmdd,date_end_yyyymmdd,folder_data);
    
    if isempty(L1)
        continue;
    end
    
    if isfield(L1,'sky_condition_index')
        L1.sci =L1.sky_condition_index;
    end
    warning('No overlap function found. Using TUB120011 overlap function')
    fid = fopen('TUB120011_20121112_1024.cfg');ov_cell = textscan(fid, '%f','headerLines',1);fclose(fid);
    ov_to_use = cell2mat(ov_cell);
    ov = [];
    
    RCS = L1.rcs_0';
    L1.range_gate = L1.range(2) -L1.range(1);
    L1.cbh = L1.cloud_base_height';
    chminfo = [];
    %----------------------------------------------------------------------
    %% Calculate overlap correction
    [ovp_fc_ok_mat,ovp_fc_ok_mat_time,ovp_fc_ok_mat_range,ind_ovp_fc_ok_good,...
        ind_ovp_fc_ok_final,ovp_fc_final,error_string] ...
        = calculate_overlap_automatic_structured_2(L1,chminfo,RCS,ov_to_use);
    
    
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
        
        out_dir = [foulder_out stn '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end
        
        device_name = L1.instrument_serial_number;
        serlom = L1.optical_module_id;
        
        
        out_file = [out_dir 'result_ovp_fc_' device_name serlom '_' date_yyyymmdd '.mat'];
        save(out_file,'result');
        
    end
    
    %----------------------------------------------------------------------
    %% Make and Save plots
    if ~isempty(ovp_fc_final)
        out_dir = [foulder_out stn '/' date_yyyymmdd(1:4) '/' date_yyyymmdd(5:6) '/'];
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end
        
        grad_raw = NaN(size(RCS));
        for j=1:length(L1.time)
            grad_raw(:,j) = 1/(2*L1.range_gate)*conv(RCS(:,j),[1 0 -1],'same');
        end
        RCS_corr = RCS.*repmat(ov_to_use,1,length(L1.time))./repmat(ovp_fc_final,1,length(L1.time));
        grad_corr = NaN(size(RCS_corr));
        
        for j=1:length(L1.time)
            grad_corr(:,j) = 1/(2*L1.range_gate)*conv(RCS_corr(:,j),[1 0 -1],'same');
        end
        
        file_out = [out_dir 'ovp_fc_final_'  stn '_' device_name serlom '_' date_yyyymmdd '.png'];
        plot_overlap_routine_ov_fct(L1,chminfo,ov,ovp_fc_ok_mat_final,ovp_fc_final,file_out);
        
        file_out = [out_dir stn '_' device_name serlom '_' date_yyyymmdd '.png'];
        plot_overlap_routine_RCS_GRADRCS(L1,chminfo,RCS,grad_raw,RCS_corr,grad_corr,file_out,result);
        
        out_file = [out_dir 'ovp_fc_final_' device_name serlom '_' date_yyyymmdd '.mat'];
        save(out_file,'ovp_fc_final');
        
    end
    
end