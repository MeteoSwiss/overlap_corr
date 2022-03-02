function [L1]=read_L1_EPROFILE_v4(WIGOS_ID,identifier,start_str,end_str,loc_folder,variables,read_all_files)
%==========================================================================
% Function read_L1_EPROFILE_v4(WIGOS_ID,identifier,start_str,[end_str],[loc_folder],[variables],[read_all_files])
% This function read  L1 Netcdf files.
% Preallocation is made using information included in netcdf files.
%
% Input parameters:  start_str,end_str : start and end time as string
% (format : 'yyyymmddHHMMSS' or 'yyyymmdd')
% If end_str is not provided it is equal to start_str+1
%
% Output:            La : structure containing data of interest
%
% Exemple of use :
% read_L1_EPROFILE_v4('16054','0','20160805000000','20160806150000','M:\pay-data\data\pay\REM\ACQ\E_PROFILE_PILOT\DATA\L1/',[],1)
% read_L1_EPROFILE_v4('16054','0','20160805000000','20160806150000')
% read_L1_EPROFILE_v4('16054','0','20160805','20160806')
% read_L1_EPROFILE_v4('16054','0','20160805')
%
% Use of external functions : none
%
% V3: Can select only some variables
% V4: Work with WIGOS IDs
% since May 2020: added input variable "read_all_files" (optional). Default is 0 --> set it to 1 in order to simultaneously read daily and instant files (if available) 

% Author: Maxime Hervo, Simone Bircher-Adrot
% Date: 2016-2018, 2020
%==========================================================================

%% Format Inputs
info.plot=0;

if nargin==0
    %by default values to test the routine
    warning('No Input,Default values')
    WIGOS_ID='0-20000-0-00201';
    identifier='A';
    start_str=['20180528' '000000'];
    end_str  =['20180529' '235909'];
end

if length(start_str)==8
    start_str=[start_str '000000'];
end
start_dn=datenum(start_str,'yyyymmddHHMMSS');

if nargin==3 || isempty(end_str)
    end_str =datestr(start_dn+0.9999,'yyyymmddHHMMSS');
end

if length(end_str)==8
    end_str=[end_str '000000'];
end
end_dn=datenum(end_str,'yyyymmddHHMMSS');

if nargin<5 || isempty(loc_folder)
    %define local paths
    if isunix()
        loc_folder='/data/pay/REM/ACQ/E_PROFILE_ALC/L1_FILES/';
    else
        loc_folder='M:\pay-data\data\pay\REM\ACQ\E_PROFILE_ALC\L1_FILES\';
    end
end

if length(WIGOS_ID)==5
    special_wmo = {'00001','00000','00002','00003'};
    special_wigos = {'0-20008-0-UEX','0-20008-0-UGR','0-20008-0-UVA','0-20008-0-MSA'};
    if any(strcmp(WIGOS_ID,special_wmo))
        WIGOS_ID = special_wigos{strcmp(WIGOS_ID,special_wmo)};
    else
        WIGOS_ID = ['0-20000-0-' WIGOS_ID];
    end
    disp('WMO_nb replaced by WIGOS')
end

if nargin<6
    variables={''};
    %variables={'time','range','attenuated_backscatter','cbh'};
end

if nargin <7 || isempty(read_all_files)
    read_all_files=0;
end

disp(['Reading L1 NetCDF for ' start_str ' to ' end_str])

%% Make file list
disp('Make file list')
time_vec=floor(start_dn):floor(end_dn);
list=[]; folder=[];
for j=1:length(time_vec)
    
    %% Check if Daily file
    str=[loc_folder  WIGOS_ID datestr(time_vec(j),'/yyyy/mm/')...
        'L1_' WIGOS_ID '_' identifier datestr(time_vec(j),'yyyymmdd') '.nc'];
    list1=dir(str);
    
    
    if isempty(list1)
        str=[loc_folder  WIGOS_ID datestr(time_vec(j),'/yyyy/mm/dd/')...
            'L1_' WIGOS_ID '_' identifier datestr(time_vec(j),'yyyymmdd') '*.nc'];
        
        list1=dir(str);
        
        if isempty(list1)
            warning(['No files in ' str ])
        else
            disp([num2str(length(list1)) ' files for ' datestr(time_vec(j))])
        end
        index=false(length(list1),1);
        format='%s%s%s%s.nc';
        
        for i=1:length(list1)
            tmp=list1(i).name;
            data=textscan(tmp,format,'delimiter',{'_','.'});
            
            time_tmp= datenum(data{3}{1}(2:end),'yyyymmddHHMM');
            
            % Hem 23/05/2014: Load data 10 min after to be sure to have entire data set
            if time_tmp<start_dn || time_tmp>end_dn+10/60/24
                index(i)=true;
            end
        end
        disp([num2str(sum(~index)) ' files selected'])
        
        list1(index)=[];
        
        list=[list;list1];
    
    %read daily and 5-min files
    elseif ~isempty(list1) && read_all_files==1
        disp('Using Daily file and 5 min files')
        str_new=[loc_folder  WIGOS_ID datestr(time_vec(j),'/yyyy/mm/dd/')...
            'L1_' WIGOS_ID '_' identifier datestr(time_vec(j),'yyyymmdd') '*.nc'];
        
        list_daily=list1;
        list1=[list1;dir(str_new)];
        
        if isempty(list1)
            warning(['No files in ' str ])
        else
            disp([num2str(length(list_daily)) ' daily files and ' num2str(length(dir(str_new))) ' 5-min files for ' datestr(time_vec(j))])
        end
        index=false(length(list1),1);
        format='%s%s%s%s.nc';
        
        for i=1:length(list1)
            if length(list1(i).name)==35 %5-min files
                tmp=list1(i).name;
                data=textscan(tmp,format,'delimiter',{'_','.'});

                time_tmp= datenum(data{3}{1}(2:end),'yyyymmddHHMM');

                % Hem 23/05/2014: Load data 10 min after to be sure to have entire data set
                if time_tmp<start_dn || time_tmp>end_dn+10/60/24
                    index(i)=true;
                end
            else
                index(i)=false;
            end
        end
        disp([num2str(sum(~index)) ' files selected'])
        
        list1(index)=[];
        
        list=[list;list1];
    else
        disp('Using Daily file')
        list=[list;list1];
    end
end

if isempty(list)
    L1=[];
    warning('No files ' )
    return
end

%% Preallocate assuming files always have the same length
number_of_profiles_total=0;
for u=1:length(list)
    file=list(u).name;
    try
        finfo=ncinfo([list(u).folder '/' list(u).name]);
        for i=1:length(finfo.Variables)
            if strcmpi(finfo.Variables(i).Name,'time')
                number_of_profiles_total=number_of_profiles_total+min(finfo.Variables(i).Size,288);
                %Problem with duplicate data: preallocation is too big

                break
            end
        end
    catch me
        disp(me)
        warning(['Error While Loading L1: ' file])
    end    
end
% time_tmp= ncread([folder(1,:)  list(1).name],'time');
% number_of_profiles_per_file = length(unique(time_tmp));

if isempty(variables) || isempty(variables{1})
    read_all_variables=1;
else
    read_all_variables=0;
    index_deletion = false(length(finfo.Variables),1);
    for i=1:length(finfo.Variables)
        if ~any(strcmp(finfo.Variables(i).Name,variables))
            index_deletion(i) = true;
        end
    end
    finfo.Variables(index_deletion)=[];
end


for i=1:length(finfo.Variables)
    switch length(finfo.Variables(i).Dimensions)
        case 0
            % No preallocation
        case 1
            if strcmp(finfo.Variables(i).Dimensions.Name,'time')
                L1.(finfo.Variables(i).Name)=NaN(number_of_profiles_total,1);
            else % No preallocation
            end
        case 2
            L1.( finfo.Variables(i).Name)=NaN(number_of_profiles_total,finfo.Variables(i).Size(1));
        otherwise
            error('Bad Number of dimensions')
    end
end



%% read (low level functions)
if read_all_variables==1
    disp(['Reading ' num2str(length(list)) ' files (with low level functions)'])
else
    disp(['Reading partially ' num2str(length(list)) ' files (with low level functions)'])
end
N_var=length(finfo.Variables);

kk=1;
for k=1:length(list)
    if floor(k/500)*500-k==0
        disp([num2str(k) ' on ' num2str(length(list))  ])
    end
    file=list(k).name;
    folder = [list(k).folder '/'];
    
    try
        ncid = netcdf.open([folder  file]);
        % read Time to calculate number of profiles in File
        varid=netcdf.inqVarID(ncid,'time');
        time_tmp= netcdf.getVar(ncid,varid);
        number_of_profiles=length(time_tmp);
        L1.time(kk:kk+number_of_profiles-1,1)=time_tmp;
        
        for i=1:N_var
            varname=finfo.Variables(i).Name;
            varid = netcdf.inqVarID(ncid,varname);
            
            switch length(finfo.Variables(i).Dimensions)
                case 0
                    L1. (varname) = double(netcdf.getVar(ncid,varid));
                case 1
                    if strcmp(finfo.Variables(i).Dimensions.Name,'time')
                        if ~strcmp(varname,'time')
                            if ~strcmp(varname,'error_string')
                                L1.(varname)(kk:kk+number_of_profiles-1,1) = netcdf.getVar(ncid,varid);
                            else
                                [~,xtype] = netcdf.inqVar(ncid,varid);
                                if xtype ==12
                                    % Matlab does not suppport to read strings in NEtcdf files
                                    
                                    % need to close files before you can open it again
                                    netcdf.close(ncid)
                                    try
                                        error_str_tmp = h5read([folder  file],'/error_string');
                                        L1.(varname)(kk:kk+number_of_profiles-1,1) = hex2dec(error_str_tmp);
                                    catch
                                        %                                     disp(me)
                                        warning(['Error while reading error string: ' folder  file] )
                                    end
                                    ncid = netcdf.open([folder  file]);
                                else
                                    L1.(varname)(kk:kk+number_of_profiles-1,1) = netcdf.getVar(ncid,varid);
                                end
                            end
                            %
                            
                        else
                            % Already done
                        end
                    else
                        if k==1
                            L1. (varname) = netcdf.getVar(ncid,varid);
                        end
                    end
                case 2
                    L1. (varname)(kk:kk+number_of_profiles-1,:) = netcdf.getVar(ncid,varid)';
                    
                otherwise
                    error('Bad Number of dimensions')
            end
        end
        kk=kk+number_of_profiles;
        
        netcdf.close(ncid)
    catch me
        disp(me)
        warning(['Error While Loading L1: ' file])
    end
end


if end_dn<datenum(2016,09,11)
    % Attribute name changed by Myles
    L1.instrument_type=ncreadatt([folder file],'/','system');
else
    L1.instrument_type=ncreadatt([folder file],'/','instrument_type');
    L1.instrument_firmware_version=ncreadatt([folder file],'/','instrument_firmware_version');
end
try
    L1.instrument_serial_number = ncreadatt([folder file],'/','instrument_serial_number');
catch
    disp('No Info about Serial numbers')
    L1.instrument_serial_number = WIGOS_ID;
end
try
    L1.optical_module_id = ncreadatt([folder file],'/','optical_module_id');
catch
    disp('No Info about optical Serial numbers')
    L1.optical_module_id = identifier;
end
try
    L1.site_location = ncreadatt([folder file],'/','site_location');
catch
    disp('No Info about Location')
    L1.site_location = 'Unknown';
end



%convert time in matlab format
% L1.time_raw=L1.measurement_time;
L1.time_raw=L1.time;
L1.time=datenum(1970,1,1)+L1.time_raw;
L1.nb_file=length(list);

% convert to double
if isfield(L1,'range')
    L1.range = double(L1.range);
end

%% Remove NaN
fill_value = 9.969209968386869e+36;
fields=fieldnames(L1);
for i=1:length(fields)
    L1.(fields{i})(L1.(fields{i})==fill_value) = NaN;
end
%% remove values outside define time_zone & NaNs
index=L1.time<start_dn | L1.time>end_dn | isnan(L1.time);
if any(index)
    fields=fieldnames(L1);
    for i=1:length(fields)
        if size(L1.(fields{i}),1)==size(index,1)
            L1.(fields{i})(index,:)=[];
        end
    end
end

disp('... Reading Done.')

%% Check variables names
if max(time_vec)<datenum(2016,11,2) && max(time_vec)>datenum(2016,9,30)
    L1.cloud_base_height=L1.cbh;
    switch L1.instrument_type
        case {'JENOPTIK NIMBUS'}
            L1.cloud_cover=L1.bcc;
        case {'VAISALA CL51','VAISALA CL31','CAMPBEL SCIENTIFIC CS135'}
            L1.cloud_cover=L1.cf;
    end
end
%% plot
if info.plot==1
    figure
    subplot(1,3,1:2)
    pcolor(L1.time,L1.range,log10(abs(L1.rcs_0))')
    hold on
    
    plot(L1.time,L1.cloud_base_height,'sk','markerfacecolor','k')
    shading flat
    datetick
    colorbar
    %     caxis([-2 1])
    colormap jet
    title([WIGOS_ID '-' identifier ' :' L1.instrument_type ])
    
    subplot(1,3,3)
    
    plot(nanmedian(L1.rcs_0),L1.range)
    %     set(gca,'xscale','log')
    
end



