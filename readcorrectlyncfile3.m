
function [chm,info] = readcorrectlyncfile3(stn,start_str,end_str,root_folder)

do_plot = false;

disp('------- Getting ncfiles -------');

%% make list of daily files
if nargin<2
    [filenames,path] = uigetfile('*.nc', 'Select a single (or consecutive) NetCDF file(s)','MultiSelect','on');

    if(~iscell(filenames))
        if(filenames==0)
            warning('No file selected');
            chm = [];
            info = [];
            return;
        end
        if ischar(filenames)
            filenames = {filenames};
        end
    end
    for k=1:length(filenames)
        filenames{k} = fullfile(path,filenames{k});
        disp(filenames{k});
    end
else
    if length(start_str)==8
        start_str = [start_str '000000'];
    end
    start_dn = datenum(start_str,'yyyymmddHHMMSS');

    if nargin==2
        end_str = datestr(start_dn+0.99999,'yyyymmddHHMMSS');
        if isunix()
            root_folder='/data/pay/REM/ACQ/CEILO_CHM15k/NetCDF/daily/';
        else
            root_folder='M:\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
    %         root_folder='C:\AllData\';
        end
    end
    
    if nargin==3
        if(strcmp(end_str(end),'/') || strcmp(end_str(end),'\'))
            root_folder = end_str;
            end_str = datestr(start_dn+0.99999,'yyyymmddHHMMSS');
        else
            if isunix()
                root_folder='/data/pay/REM/ACQ/CEILO_CHM15k/NetCDF/daily/';
            else
                root_folder='M:\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
%                 root_folder='C:\AllData\';
            end
        end
    end

%     if nargin==3 || nargin==2
%         if isunix()
%             root_folder='/data/pay/REM/ACQ/CEILO_CHM15k/NetCDF/daily/';
%         else
%             root_folder='M:\pay-data\data\pay\REM\ACQ\CEILO_CHM15k\NetCDF\daily\';
%     %         root_folder='C:\AllData\';
%         end
%     end

    if length(end_str)==8 || (length(end_str)==14 && strcmp(end_str(9:14),'000000'))
        end_dn = datenum([end_str '000000'],'yyyymmddHHMMSS')-1+0.99999;
        end_str = datestr(end_dn,'yyyymmddHHMMSS');
    else
        end_dn = datenum(end_str,'yyyymmddHHMMSS');
    end


    % list daily files
    filenames = {};

    for time_day = floor(start_dn):floor(end_dn)
        [year,month,day]=datevec(time_day);
        
        stn_str = {stn};
        if(strcmp(stn,'jenaCHM120106') && time_day<datenum('20121112','yyyymmdd'))
            stn_str = {'einlaufen Laser ME120011'};
        elseif(strcmp(stn,'jenaCHM120106') && time_day==datenum('20121112','yyyymmdd'))
            stn_str = {'einlaufen Laser ME120011','ME120011 24h Test'};
        elseif(strcmp(stn,'jenaCHM120106') && time_day<datenum('20121114','yyyymmdd'))
            stn_str = {'ME120011 24h Test'};
        elseif(strcmp(stn,'pay') && time_day<datenum('20130211','yyyymmdd'))
            stn_str = {'NN'};
        elseif(strcmp(stn,'pay') && time_day==datenum('20130211','yyyymmdd'))
            stn_str = {'NN','payerne'};
        elseif(strcmp(stn,'pay') && time_day<datenum('20130613','yyyymmdd'))
            stn_str = {'payerne'};
        elseif(strcmp(stn,'pay') && time_day==datenum('20130613','yyyymmdd'))
            stn_str = {'payerne','pay'};
        end
        
        if(strcmp(stn,'jenaCHM130104') && time_day<datenum('20140623','yyyymmdd'))
            stn_str = {'NN'};
        elseif(strcmp(stn,'kse') && time_day==datenum('20140623','yyyymmdd'))
            stn_str = {'NN','jun'};
        elseif(strcmp(stn,'kse') && time_day<datenum('20140819','yyyymmdd'))
            stn_str = {'jun'};
        elseif(strcmp(stn,'kse') && time_day==datenum('20140819','yyyymmdd'))
            stn_str = {'jun','kse'};
        end
        
        daystr = num2str(day,'%02.0f');
        monthstr = num2str(month,'%02.0f');
        yearstr = num2str(year,'%02.0f');
        disp(['------- ' yearstr  monthstr  daystr  '  -----'])
        if strcmp(stn_str,'SIRTA')
            folder = [root_folder '' yearstr '/' monthstr '/' daystr '/'];
        elseif strcmp(stn_str,'Granada')
            folder = [root_folder '' yearstr '/' ];

        else
            folder = [root_folder '' yearstr '/' monthstr '/'];
        end
        
        % special cases
        if(strcmp(stn,'jenaCHM120106') && time_day==datenum('20121114','yyyymmdd'))
            filenames_special = {[folder yearstr  monthstr  daystr '_ME120011 24h Test_CHM120106_000.nc'],...
                [folder yearstr  monthstr  daystr '_ME120011 24h Test_CHT120106_001.nc'],...
                [folder yearstr  monthstr  daystr '_ME120011 24h Test_CHM120106_002.nc'],...
                [folder yearstr  monthstr  daystr '_NN_CHM120106_003.nc']};
            for j=1:length(filenames_special)
                list = dir(filenames_special{j});
%                 disp(filenames_special{j});
                disp([filenames_special{j},', ',num2str(list(1).bytes/1024,'%12.4f'),'KB']);
            end
            filenames = {filenames{:},filenames_special{:}};
            continue;
        end
        if(strcmp(stn,'pay') && time_day==datenum('20130508','yyyymmdd'))
            filenames_special = {[folder yearstr  monthstr  daystr '_payerne_CHM120106_000.nc'],...
                [folder yearstr  monthstr  daystr '_NN_CHMstddrd_001.nc'],...
                [folder yearstr  monthstr  daystr '_payerne_CHMstddrd_002.nc'],...
                [folder yearstr  monthstr  daystr '_payerne_CHM120106_003.nc']};
            for j=1:length(filenames_special)
                list = dir(filenames_special{j});
%                 disp(filenames_special{j});
                disp([filenames_special{j},', ',num2str(list(1).bytes/1024,'%12.5f'),'KB']);
            end
            filenames = {filenames{:},filenames_special{:}};
            continue;
        end
        
        for l=1:length(stn_str)
            % check if daily netcdf available
            %             list = dir([folder yearstr  monthstr  daystr '_' lower(stn_str{l}) '_CH*.nc']);
            if strcmp(stn_str{l},'SIRTA')
                str = [folder 'chm15k_*_*_*_' yearstr  monthstr  daystr  '*.nc'];
            else
                str = [folder yearstr  monthstr  daystr '_' stn_str{l} '_CH*.nc'];
            end
            list = dir(str);
            
            if isempty(list)
                %                 disp([folder yearstr  monthstr  daystr '_' lower(stn_str{l}) '_CH*.nc',' : file(s) not found.']);
                disp([str,' : file(s) not found.']);
                continue;
            else
                for j=1:length(list)
                    daily_file = [folder list(j).name];
%                     disp(daily_file);
                    disp([daily_file,', ',num2str(list(j).bytes/1024),'KB']);
                    filenames = {filenames{:},daily_file};
                end
            end
        end
    end
end
    
%% read list of nc files

chm_list = [];

for k=1:length(filenames)
    
    chm = [];
%     ncfilename = fullfile(path,filenames{k});
    ncfilename = filenames{k};
    info = ncinfo(ncfilename);
    
    % read file name
    chm.filename = info.Filename;
    
    for i=1:length(info.Attributes)
        % read global attribute
        chm.(info.Attributes(i).Name) = ncreadatt(ncfilename,'/',info.Attributes(i).Name);
    end
        
    for i=1:length(info.Variables)
        % convert 'scale_factor' attributes to numeric values
        for j=1:length(info.Variables(i).Attributes)
            if(strcmp(info.Variables(i).Attributes(j).Name,'scale_factor'))
                scale_factor = ncreadatt(ncfilename,info.Variables(i).Name,'scale_factor');
                if(~isnumeric(scale_factor))
                    ncwriteatt(ncfilename,info.Variables(i).Name,'scale_factor',1./str2num(scale_factor));
                end
            end
        end
        % read variable
        chm.(info.Variables(i).Name) = double(ncread(ncfilename,info.Variables(i).Name));
        % rewrite 'scale_factor' attributes to their original values
%         for j=1:length(info.Variables(i).Attributes)
%             if(strcmp(info.Variables(i).Attributes(j).Name,'scale_factor'))
%                 ncwriteatt(ncfilename,info.Variables(i).Name,'scale_factor',scale_factor);
%             end
%         end
    end
    
    % convert time to MATLAB time format
    chm.time = datenum(1904,1,1)+chm.time/3600/24;
    if isunix
        root_folder = '/data/pay/PBL4EMPA/overlap_correction/overlap_functions_Lufft/';
    else
        root_folder = '\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\overlap_correction\overlap_functions_Lufft\';    
    end
    % get overlap function
    list = dir([root_folder,chm.serlom,'_*_1024.cfg']);
    if isempty(list)
        warning('missing overlap function');
        ov = NaN(1024,1);
        scaling = NaN;
    else
        if(strcmp(chm.serlom,'TUB120001'))
            if time_day<datenum(2012,09,17)
                fname_overlapfc = list(1).name;
            else
                fname_overlapfc = list(2).name;
            end
        else
            fname_overlapfc = list(1).name;
        end
        disp(['reading ' fname_overlapfc]);
        fid = fopen([root_folder fname_overlapfc]);
        ov_cell = textscan(fid, '%f','headerLines',1);
        frewind(fid);
        scaling_cell = textscan(fid, 'scaling: %f',1,'headerLines',4);
        fclose(fid);
        ov = cell2mat(ov_cell);
        scaling = cell2mat(scaling_cell);
    end
    
    
    if(~isempty(chm.time))
        chm.time_ncfiles = chm.time(1);
    else
%         chm.time_ncfiles = [];
%         chm.time_ncfiles = NaN;
        chm.time_ncfiles = datenum([sprintf('%02.0f',chm.year),sprintf('%02.0f',chm.month),sprintf('%02.0f',chm.day)],'yyyymmdd');
    end

    % to adapt software version '12.03.1 2.13 0.559' and earlier
    if(~isfield(chm,'range_gate_hr'))
        chm.range_gate_hr = 4.994999885559082;
    end
    if(~isfield(chm,'range_hr'))
%         chm.range_hr = [0.00100000004749700;4.99499988555908;9.98999977111816;14.9849996566772;19.9799995422363;24.9750003814697;29.9699993133545;34.9650001525879;39.9599990844727;44.9550018310547;49.9500007629395;54.9449996948242;59.9399986267090;64.9349975585938;69.9300003051758;74.9250030517578;79.9199981689453;84.9150009155273;89.9100036621094;94.9049987792969;99.9000015258789;104.894996643066;109.889999389648;114.885002136230;119.879997253418;124.875000000000;129.869995117188;134.865005493164;139.860000610352;144.854995727539;149.850006103516;154.845001220703];        
        chm.range_hr = [0.0010000000474974513;4.994999885559082;9.9899997711181641;14.984999656677246;19.979999542236328;24.975000381469727;29.969999313354492;34.965000152587891;39.959999084472656;44.955001831054688;49.950000762939453;54.944999694824219;59.939998626708984;64.93499755859375;69.930000305175781;74.925003051757813;79.919998168945313;84.915000915527344;89.910003662109375;94.904998779296875;99.900001525878906;104.89499664306641;109.88999938964844;114.88500213623047;119.87999725341797;124.875;129.8699951171875;134.86500549316406;139.86000061035156;144.85499572753906;149.85000610351562;154.84500122070312];
    end
    if(~isfield(chm,'beta_raw_hr'))
        chm.beta_raw_hr = NaN(32,size(chm.beta_raw,2));
    end
    if(~isfield(chm,'scaling'))
        switch chm.serlom
            case 'TUB120011'
                chm.scaling = 0.35976201295852661;
            case 'TUB140005'
                chm.scaling = 0.29424801468849182;
            case 'TUB140007'
                chm.scaling = 0.129040002822876;
            otherwise
%                 disp([chm.filename,' (',sprintf('%02.0f',chm.year),'-',sprintf('%02.0f',chm.month),'-',sprintf('%02.0f',chm.day),') : serlom ',chm.serlom,' unknown : using scaling=1']);
%                 chm.scaling = 1;
                if isnan(scaling)
                    disp([chm.filename,' (',sprintf('%02.0f',chm.year),'-',sprintf('%02.0f',chm.month),'-',sprintf('%02.0f',chm.day),') : serlom ',chm.serlom,' unknown : using scaling=1']);
                    chm.scaling = 1;
                else
                    chm.scaling = scaling;
                end
        end
    end
    if(~isfield(chm,'p_calc'))
        chm.p_calc = chm.nn2*1e-5;
        if strcmp(chm.serlom,'TUB120011')
            ind_change = chm.time < datenum(2013,09,04,08,46,50);
        else
            ind_change = 1:length(chm.time);
        end
%         fact = 2.9833272964699868;
        fact = 3;% personnal note from Kornelia Pönitz from Lufft
        chm.p_calc(ind_change) = chm.p_calc(ind_change) * fact;
        
        chm.nn3 = NaN(size(chm.nn2));
        if(~isempty(chm.time))
            % background substracted and normalized by laser shot number and stddev
            chm.beta_raw_old = chm.beta_raw;
%             chm.SNR = applyoverlapfc(chm.beta_raw,'none','ref');
            chm.SNR = chm.beta_raw./repmat(ov,1,size(chm.beta_raw,2));
            % normalized range corrected signal
            chm.beta_raw = chm.SNR .* ((ones(size(chm.beta_raw,1),1)*(chm.stddev')) .* ((chm.range.^2)*ones(1,size(chm.beta_raw,2)))) ./ (chm.scaling*ones(size(chm.beta_raw,1),1)*chm.p_calc');
        else
            chm.beta_raw_old = [];
            chm.SNR = [];
            chm.beta_raw = [];
        end
    else
        if(~isempty(chm.time))
            % background substracted and normalized by laser shot number and stddev
%             chm.beta_raw_old = applyoverlapfc(chm.beta_raw,'ref','none') .* (chm.scaling*ones(size(chm.beta_raw,1),1)*chm.p_calc') ./ ((ones(size(chm.beta_raw,1),1)*(chm.stddev')) .* ((chm.range.^2)*ones(1,size(chm.beta_raw,2))));
            chm.beta_raw_old = (chm.beta_raw.*repmat(ov,1,size(chm.beta_raw,2))) .* (chm.scaling*ones(size(chm.beta_raw,1),1)*3*chm.p_calc') ./ ((ones(size(chm.beta_raw,1),1)*(chm.stddev')) .* ((chm.range.^2)*ones(1,size(chm.beta_raw,2))));            
            chm.SNR = chm.beta_raw .* (chm.scaling*ones(size(chm.beta_raw,1),1)*chm.p_calc') ./ ((ones(size(chm.beta_raw,1),1)*(chm.stddev')) .* ((chm.range.^2)*ones(1,size(chm.beta_raw,2))));
        else
            chm.beta_raw_old = [];
            chm.SNR = [];
        end
    end
    
    if(~isempty(chm.time))
        % photon count per laser pulse
%         chm.signal_raw_per_pulse = ...
%         (applyoverlapfc(chm.SNR,'ref','none').*(ones(size(chm.SNR,1),1)*(chm.stddev'))) + ones(size(chm.SNR,1),1)*chm.base';
        chm.signal_raw_per_pulse = ...
        ((chm.beta_raw.*repmat(ov,1,size(chm.beta_raw,2))).*(ones(size(chm.SNR,1),1)*(chm.stddev'))) + ones(size(chm.SNR,1),1)*chm.base';
        % photon count during averaging period
        chm.signal_raw = chm.signal_raw_per_pulse .* (ones(size(chm.SNR,1),1)*chm.laser_pulses');  
    else
        chm.signal_raw_per_pulse = [];
        chm.signal_raw = [];
    end

    % to be consistent with sofware version '12.12.1 2.13 0.723' and later
    if(~isfield(chm,'device_name'))
        chm.device_name = chm.source;
        chm = rmfield(chm,'device');
    end

    chm_list.(['nc',num2str(k)]) = chm;

end

% check_changes(chm_list);

%% append files

chm = [];
% time-dependent variables
chm.time = [];
chm.time_raw = [];
chm.average_time = [];
chm.life_time = [];
chm.laser_pulses = [];
chm.error_ext = [];
chm.nn1 = [];
chm.nn2 = [];
chm.nn3 = [];
chm.base = [];
chm.stddev = [];
chm.p_calc = [];
chm.state_detector = [];
chm.state_laser = [];
chm.state_optics = [];
chm.temp_det = [];
chm.temp_ext = [];
chm.temp_int = [];
chm.temp_lom = [];
chm.pbl = [];
chm.pbs = [];
chm.cbh = [];
chm.cbe = [];
chm.cdp = [];
chm.cde = [];
chm.sci = [];
chm.mxd = [];
chm.vor = [];
chm.voe = [];
chm.tcc = [];
chm.bcc = [];

% time- and range-dependent variables
chm.beta_raw_old = [];
chm.beta_raw = [];
chm.beta_raw_hr = [];
chm.SNR = [];
chm.signal_raw = [];
chm.signal_raw_per_pulse = [];


% same value for all CHM15k-Nimbus files
chm.range = [];
chm.range_gate = [];
chm.range_hr  = [];
chm.range_gate_hr = [];
chm.wavelength = [];

% file-dependent variables
% depends on serlom
chm.scaling = [];

% depends on location and orientation of ceilometer
chm.latitude = [];
chm.longitude = [];
chm.altitude = [];
chm.azimuth = [];
chm.zenith = [];
% depends on location and orientation of ceilometer and of user input
chm.cho = [];
% depends on user input
chm.layer = [];

% name of the ncdf-file
chm.filename = {};
% global attributes of the ncdf-file
chm.day = {};
chm.month = {};
chm.year = {};
chm.location = {};
chm.title = {};
chm.source = {};
chm.device_name = {};
chm.institution = {};
chm.software_version = {};
chm.comment = {};
chm.serlom = {};

chm.time_ncfiles = [];

chm.scaling_ncfiles = {};
chm.latitude_ncfiles = {};
chm.longitude_ncfiles = {};
chm.altitude_ncfiles = {};
chm.azimuth_ncfiles = {};
chm.zenith_ncfiles = {};
chm.cho_ncfiles = {};
chm.layer_ncfiles = {};



for k=1:length(filenames)
    nc = chm_list.(['nc',num2str(k)]);
    
    chm.time = [chm.time;nc.time];
    chm.average_time = [chm.average_time;nc.average_time];
    chm.life_time = [chm.life_time;nc.life_time];
    chm.laser_pulses = [chm.laser_pulses;nc.laser_pulses];
    chm.error_ext = [chm.error_ext;nc.error_ext];
    chm.nn1 = [chm.nn1;nc.nn1];
    chm.nn2 = [chm.nn2;nc.nn2];
    chm.nn3 = [chm.nn3;nc.nn3];
    chm.base = [chm.base;nc.base];
    chm.stddev = [chm.stddev;nc.stddev];
    chm.p_calc = [chm.p_calc;nc.p_calc];
    chm.state_detector = [chm.state_detector;nc.state_detector];
    chm.state_laser = [chm.state_laser;nc.state_laser];
    chm.state_optics = [chm.state_optics;nc.state_optics];
    chm.temp_det = [chm.temp_det;nc.temp_det];
    chm.temp_ext = [chm.temp_ext;nc.temp_ext];
    chm.temp_int = [chm.temp_int;nc.temp_int];
    chm.temp_lom = [chm.temp_lom;nc.temp_lom];
    chm.pbl = [chm.pbl,nc.pbl];
    chm.pbs = [chm.pbs,nc.pbs];
    chm.cbh = [chm.cbh,nc.cbh];
    chm.cbe = [chm.cbe,nc.cbe];
    chm.cdp = [chm.cdp,nc.cdp];
    chm.cde = [chm.cde,nc.cde];
    chm.sci = [chm.sci;nc.sci];
    chm.mxd = [chm.mxd;nc.mxd];
    chm.vor = [chm.vor;nc.vor];
    chm.voe = [chm.voe;nc.voe];
    chm.tcc = [chm.tcc;nc.tcc];
    chm.bcc = [chm.bcc;nc.bcc];

    chm.beta_raw_old = [chm.beta_raw_old,nc.beta_raw_old];
    chm.beta_raw = [chm.beta_raw,nc.beta_raw];
    chm.beta_raw_hr = [chm.beta_raw_hr,nc.beta_raw_hr];
    chm.SNR = [chm.SNR,nc.SNR];
    chm.signal_raw = [chm.signal_raw,nc.signal_raw];
    chm.signal_raw_per_pulse = [chm.signal_raw_per_pulse,nc.signal_raw_per_pulse];

    chm.range = nc.range;
    chm.range_gate = nc.range_gate;
    chm.range_hr = nc.range_hr;
    chm.range_gate_hr = nc.range_gate_hr;
    chm.wavelength = nc.wavelength;


    chm.time_ncfiles = [chm.time_ncfiles;nc.time_ncfiles];
    chm.scaling_ncfiles = {chm.scaling_ncfiles{:},nc.scaling};
    chm.latitude_ncfiles = {chm.latitude_ncfiles{:},nc.latitude};
    chm.longitude_ncfiles = {chm.longitude_ncfiles{:},nc.longitude};
    chm.altitude_ncfiles = {chm.altitude_ncfiles{:},nc.altitude};
    chm.azimuth_ncfiles = {chm.azimuth_ncfiles{:},nc.azimuth};
    chm.zenith_ncfiles = {chm.zenith_ncfiles{:},nc.zenith};
    chm.cho_ncfiles = {chm.cho_ncfiles{:},nc.cho};
    chm.layer_ncfiles = {chm.layer_ncfiles{:},nc.layer};
    chm.filename = {chm.filename{:},nc.filename};
    chm.day = {chm.day{:},nc.day};
    chm.month = {chm.month{:},nc.month};
    chm.year = {chm.year{:},nc.year};
    chm.location = {chm.location{:},nc.location};
    chm.title = {chm.title{:},nc.title};
    chm.source = {chm.source{:},nc.source};
    chm.device_name = {chm.device_name{:},nc.device_name};
    chm.institution = {chm.institution{:},nc.institution};
    chm.software_version = {chm.software_version{:},nc.software_version};
    chm.comment = {chm.comment{:},nc.comment};
    chm.serlom = {chm.serlom{:},nc.serlom};
end

chm.time_raw = (chm.time-datenum(1904,1,1))*(24*3600);


if(isempty(chm.time))
    chm = [];
    info = [];
    return;
else
    % to be consistent with scripts that were adapted with ancient version
    % of get_chm15k_from_files.m
    chm.scaling = chm.scaling_ncfiles{1};
    chm.latitude = chm.latitude_ncfiles{1};
    chm.longitude = chm.longitude_ncfiles{1};
    chm.altitude = chm.altitude_ncfiles{1};
    chm.azimuth = chm.azimuth_ncfiles{1};
    chm.zenith = chm.zenith_ncfiles{1};
    chm.cho = chm.cho_ncfiles{1};
    chm.layer = chm.layer_ncfiles{1};
    
    remove_values = false;
    index = false(size(chm.time));
    if(exist('start_str','var') && exist('end_str','var'))
        %remove values outside defined time interval
        if any(or(chm.time<start_dn, chm.time>end_dn))
            index=or(chm.time<start_dn,chm.time>end_dn);
            remove_values = true;
        end
    else
        if any(chm.time<0)
            index = chm.time<0;
            remove_values = true;
        end
    end

    if(remove_values)
        % time-dependent variables
        chm.time(index) = [];
        chm.time_raw(index) = [];
        chm.average_time(index) = [];
        chm.life_time(index) = [];
        chm.laser_pulses(index) = [];
        chm.error_ext(index) = [];
        chm.nn1(index) = [];
        chm.nn2(index) = [];
        chm.nn3(index) = [];
        chm.base(index) = [];
        chm.stddev(index) = [];
        chm.p_calc(index) = [];
        chm.state_detector(index) = [];
        chm.state_laser(index) = [];
        chm.state_optics(index) = [];
        chm.temp_det(index) = [];
        chm.temp_ext(index) = [];
        chm.temp_int(index) = [];
        chm.temp_lom(index) = [];
        chm.pbl(:,index) = [];
        chm.pbs(:,index) = [];
        chm.cbh(:,index)= [];
        chm.cbe(:,index) = [];
        chm.cdp(:,index) = [];
        chm.cde(:,index) = [];
        chm.sci(index) = [];
        chm.mxd(index) = [];
        chm.vor(index) = [];
        chm.voe(index) = [];
        chm.tcc(index) = [];
        chm.bcc(index) = [];

        % time- and range-dependent variables
        chm.beta_raw_old(:,index) = [];
        chm.beta_raw(:,index) = [];
        chm.beta_raw_hr(:,index) = [];
        chm.SNR(:,index) = [];
        chm.signal_raw(:,index) = [];
        chm.signal_raw_per_pulse(:,index) = [];
    end
    
    [~,iunique,~] = unique(chm.time);
    if(length(iunique)<length(chm.time))
        % time-dependent variables
        chm.time = chm.time(iunique);
        chm.time_raw = chm.time_raw(iunique);
        chm.average_time = chm.average_time(iunique);
        chm.life_time = chm.life_time(iunique);
        chm.laser_pulses = chm.laser_pulses(iunique);
        chm.error_ext = chm.error_ext(iunique);
        chm.nn1 = chm.nn1(iunique);
        chm.nn2 = chm.nn2(iunique);
        chm.nn3 = chm.nn3(iunique);
        chm.base = chm.base(iunique);
        chm.stddev = chm.stddev(iunique);
        chm.p_calc = chm.p_calc(iunique);
        chm.state_detector = chm.state_detector(iunique);
        chm.state_laser = chm.state_laser(iunique);
        chm.state_optics = chm.state_optics(iunique);
        chm.temp_det = chm.temp_det(iunique);
        chm.temp_ext = chm.temp_ext(iunique);
        chm.temp_int = chm.temp_int(iunique);
        chm.temp_lom = chm.temp_lom(iunique);
        chm.pbl = chm.pbl(:,iunique);
        chm.pbs = chm.pbs(:,iunique);
        chm.cbh = chm.cbh(:,iunique);
        chm.cbe = chm.cbe(:,iunique);
        chm.cdp = chm.cdp(:,iunique);
        chm.cde = chm.cde(:,iunique);
        chm.sci = chm.sci(iunique);
        chm.mxd = chm.mxd(iunique);
        chm.vor = chm.vor(iunique);
        chm.voe = chm.voe(iunique);
        chm.tcc = chm.tcc(iunique);
        chm.bcc = chm.bcc(iunique);

        % time- and range-dependent variables
        chm.beta_raw_old = chm.beta_raw_old(:,iunique);
        chm.beta_raw = chm.beta_raw(:,iunique);
        chm.beta_raw_hr = chm.beta_raw_hr(:,iunique);
        chm.SNR = chm.SNR(:,iunique);
        chm.signal_raw = chm.signal_raw(:,iunique);
        chm.signal_raw_per_pulse = chm.signal_raw_per_pulse(:,iunique);        
    end
    
    
end

%% Check changes in ncfiles
disp('------- Checking changes in ncfiles -------');
for j=2:length(chm.time_ncfiles)
    if chm.scaling_ncfiles{j} ~= chm.scaling_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : scaling changed from ',num2str(chm.scaling_ncfiles{j-1}),' to ',num2str(chm.scaling_ncfiles{j})]);
    end
    if chm.latitude_ncfiles{j} ~= chm.latitude_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : latitude changed from ',num2str(chm.latitude_ncfiles{j-1}),' to ',num2str(chm.latitude_ncfiles{j})]);
    end
    if chm.longitude_ncfiles{j} ~= chm.longitude_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : longitude changed from ',num2str(chm.longitude_ncfiles{j-1}),' to ',num2str(chm.longitude_ncfiles{j})]);    
    end
    if chm.altitude_ncfiles{j} ~= chm.altitude_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : altitude changed from ',num2str(chm.altitude_ncfiles{j-1}),' to ',num2str(chm.altitude_ncfiles{j})]);        
    end
    if chm.azimuth_ncfiles{j} ~= chm.azimuth_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : azimuth changed from ',num2str(chm.azimuth_ncfiles{j-1}),' to ',num2str(chm.azimuth_ncfiles{j})]);        
    end
    if chm.zenith_ncfiles{j} ~= chm.zenith_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : zenith changed from ',num2str(chm.zenith_ncfiles{j-1}),' to ',num2str(chm.zenith_ncfiles{j})]);        
    end
    if chm.cho_ncfiles{j} ~= chm.cho_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : cho changed from ',num2str(chm.cho_ncfiles{j-1}),' to ',num2str(chm.cho_ncfiles{j})]);        
    end
    if chm.layer_ncfiles{j} ~= chm.layer_ncfiles{j-1};
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : layer changed from ',num2str(chm.layer_ncfiles{j-1}'),' to ',num2str(chm.layer_ncfiles{j}')]);        
    end
    if ~strcmp(chm.location{j},chm.location{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : location changed from ',chm.location{j-1},' to ',chm.location{j}]);        
    end
    if ~strcmp(chm.title{j},chm.title{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : title changed from ',chm.title{j-1},' to ',chm.title{j}]);  
    end
    if ~strcmp(chm.source{j},chm.source{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : source changed from ',chm.source{j-1},' to ',chm.source{j}]);          
    end
    if ~strcmp(chm.device_name{j},chm.device_name{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : device_name changed from ',chm.device_name{j-1},' to ',chm.device_name{j}]);         
    end
    if ~strcmp(chm.institution{j},chm.institution{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : institution changed from ',chm.institution{j-1},' to ',chm.institution{j}]);          
    end
    if ~strcmp(chm.software_version{j},chm.software_version{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : software_version changed from ',chm.software_version{j-1},' to ',chm.software_version{j}]);         
    end
    if ~strcmp(chm.comment{j},chm.comment{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : comment changed from ',chm.comment{j-1},' to ',chm.comment{j}]);          
    end
    if ~strcmp(chm.serlom{j},chm.serlom{j-1});
        disp([chm.filename{j},' (',sprintf('%02.0f',chm.year{j}),'-',sprintf('%02.0f',chm.month{j}),'-',sprintf('%02.0f',chm.day{j}),'), ',datestr(chm.time_ncfiles(j)),' : serlom changed from ',chm.serlom{j-1},' to ',chm.serlom{j}]);          
    end
end
for j=2:length(chm.average_time)
    if chm.average_time(j) ~= chm.average_time(j-1)
        index = find(chm.time_ncfiles<=chm.time(j),1,'last');
        if(~isempty(index))
            disp([chm.filename{index},' (',sprintf('%02.0f',chm.year{index}),'-',sprintf('%02.0f',chm.month{index}),'-',sprintf('%02.0f',chm.day{index}),'), ',datestr(chm.time(j)),' : average_time changed from ',num2str(chm.average_time(j-1)),' to ',num2str(chm.average_time(j))]);          
        else
            % should not enter here normally...
            disp([datestr(chm.time(j)),' : average_time changed from ',num2str(chm.average_time(j-1)),' to ',num2str(chm.average_time(j))]);          
        end
    end
end

time_plot = [chm.time(1)-chm.average_time(1)/(24*3600*1000);chm.time];
average_time_plot = [0;chm.average_time]/(24*3600*1000);
index = find(diff(time_plot)>average_time_plot(2:end)+2.5/(24*3600));
for i=1:length(index)
       empty_time = [time_plot(index(i))+average_time_plot(index(i)+1):average_time_plot(index(i)+1):time_plot(index(i)+1)-average_time_plot(index(i)+1)]';
       if(~isempty(empty_time))
            ind = find(chm.time_ncfiles<=empty_time(1),1,'last');
            if(~isempty(ind))
                disp([chm.filename{ind},' (',sprintf('%02.0f',chm.year{ind}),'-',sprintf('%02.0f',chm.month{ind}),'-',sprintf('%02.0f',chm.day{ind}),'), ',datestr(empty_time(1)),' : no data from ',datestr(empty_time(1)),' to ',datestr(empty_time(end))]);          
            else
                % should not enter here normally...
                disp([datestr(empty_time(1)),' : no data from ',datestr(empty_time(1)),' to ',datestr(empty_time(end))]);
            end
       end
end



%%

if(~exist('stn','var'))
    stn = chm.location{end};
    if(strcmp(stn,'payerne') || (strcmp(stn,'NN') && chm.time_ncfiles(1)>=datenum('20130208','yyyymmdd')))
        stn = 'pay';
    end
    if(strcmp(stn,'jun'))
        stn = 'kse';
    end
end
if(~exist('start_str','var'))
    start_str = datestr(chm.time(1),'yyyymmddHHMMSS');
end
if(~exist('end_str','var'))
    end_str = datestr(chm.time(end),'yyyymmddHHMMSS');
end



%%
if(do_plot)
    
%%
if isunix
    load('/data/pay/PBL4EMPA/pbl_analysis/bay/pblh_instr_20102014.mat')
    load('/data/pay/PBL4EMPA/pbl_analysis/bay/pblh_sound_20102014.mat')
else
    load('\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\pbl_analysis\bay\pblh_instr_20102014.mat');
    load('\\meteoswiss.ch\mch\pay-data\data\pay\PBL4EMPA\pbl_analysis\bay\pblh_sound_20102014.mat');
end
bR_theta_int = [];
bR_thetav_int = [];
pm_theta = [];
pm_thetav = [];
xtime = [];
bR_theta_woint_sound = [];
bR_thetav_woint_sound = [];
pm_theta_sound = [];
pm_thetav_sound = [];
xtime_sound = [];
for doin = floor(datenum(start_str,'yyyymmddHHMMSS')):ceil(datenum(end_str,'yyyymmddHHMMSS'))
    doi = datestr(doin,'yyyymmdd');
    str_doi = ['d_',doi];
    
    if(isfield(pblh_instr,(str_doi)))
        bR_theta_int = [bR_theta_int,pblh_instr.(str_doi).pay.bR_theta_int];
        bR_thetav_int = [bR_thetav_int,pblh_instr.(str_doi).pay.bR_thetav_woint];
        pm_theta = [pm_theta,pblh_instr.(str_doi).pay.pm_theta];
        pm_thetav = [pm_thetav,pblh_instr.(str_doi).pay.pm_thetav];
        xtime = [xtime,linspace(datenum(doi,'yyyymmdd'),datenum(doi,'yyyymmdd')+1,145)];
    end
    if(isfield(pblh_sound,(str_doi)))
        bR_theta_woint_sound = [bR_theta_woint_sound,pblh_sound.(str_doi).pay.bR_theta_woint];
        bR_thetav_woint_sound = [bR_thetav_woint_sound,pblh_sound.(str_doi).pay.bR_thetav_woint];
        pm_theta_sound = [pm_theta_sound,pblh_sound.(str_doi).pay.pm_theta];
        pm_thetav_sound = [pm_thetav_sound,pblh_sound.(str_doi).pay.pm_thetav];
        xtime_sound = [xtime_sound,datenum(doi,'yyyymmdd')+0.5];
    end
end
%%

offset = datenum(2000,1,1)-1;

if strcmp(chm.serlom,'TUB120011') || strcmp(chm.serlom,'TUB140005') || strcmp(chm.serlom,'TUB140007')
    ZE = 90-chm.zenith*100;
else
    ZE = 90-chm.zenith;
end

min_alt = 0;
% max_alt = 5000;
max_alt = 15000;
xlabel_str = 'Time UT [h]';
ylabel_str = 'Elevation [m asl]';
xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset]; %cf.displayParam.left_right_margin*[-1 1];
xticks = xlims(1)-mod(mod(xlims(1),1)*24,3)/24:2/24:xlims(2);
% yticks = [0:500:5000];
yticks = [0:1000:15000];
% ytickslabels = {'0','500','1000','1500','2000','2500','3000','3500','4000','4500','5000'};
ytickslabels = {'0','1000','2000','3000','4000','5000','6000','7000','8000','9000','10000','11000','12000','13000','14000','15000'};
% cmap = jethighblackcmap;
% clims = [4.5 7.5];
load('ypcmap2');
% cmap(end,:) = [0 0 0];
clims = [3.5 6];


figure('Color','w','Position',[0 0 1680 1050]);%,'Visible','off');
hold on;grid on;box on;

ah1 = subplot(5,1,1:2);


var_to_plot = chm.beta_raw;

% time_plot = [floor(chm.time(1));chm.time;ceil(chm.time(end))];
% time_plot_final = time_plot;
% average_time_plot = [0;chm.average_time;0]/(24*3600*1000);
% var_to_plot_final = [NaN(size(var_to_plot,1),1),var_to_plot,NaN(size(var_to_plot,1),1)];

% time_plot = [floor(chm.time(1));chm.time];
time_plot = [chm.time(1)-chm.average_time(1)/(24*3600*1000);chm.time];
time_plot_final = time_plot;
average_time_plot = [0;chm.average_time]/(24*3600*1000);
var_to_plot_final = [NaN(size(var_to_plot,1),1),var_to_plot];

index = find(diff(time_plot)>average_time_plot(2:end)+2.5/(24*3600));
for i=1:length(index)
       empty_time = [time_plot(index(i))+average_time_plot(index(i)+1):average_time_plot(index(i)+1):time_plot(index(i)+1)-average_time_plot(index(i)+1)]';
       if(~isempty(empty_time))
           ind1 = find(time_plot_final<empty_time(1),1,'last');
           ind2 = find(time_plot_final>empty_time(end),1,'first');
           time_plot_final = [time_plot_final(1:ind1);empty_time;time_plot_final(ind2:end)];
           var_to_plot_final = [var_to_plot_final(:,1:ind1),NaN(size(var_to_plot_final,1),length(empty_time)),var_to_plot_final(:,ind2:end)];
       end
end

time_plot_final = time_plot_final-offset;
var_to_plot_final(var_to_plot_final<=0)= 1000;
var_to_plot_final = log10(var_to_plot_final);

X = time_plot_final;
Y = [0*sind(ZE)+chm.altitude;chm.range*sind(ZE)+chm.altitude];
C = [var_to_plot_final(:,2:end) NaN(size(var_to_plot,1),1);NaN(1,size(var_to_plot_final,2))];
pcolor(X,Y,C);
% pcolor(time_plot_final,chm.range+chm.altitude(end),var_to_plot_final);
colormap(cmap);





% var_to_plot(var_to_plot<=0) = 10^3;
% pcolor(chm.time-offset,chm.range+chm.altitude(end),log10(var_to_plot));

shading flat;caxis(clims);
hcb = colorbar;
ylabel(hcb,'log10(Normalized range corrected signal)','FontWeight','demi');

% plot current time mark
hold on
plot((now-offset)*[1 1],[min_alt max_alt],'r--')
text((now-offset),-0.055*(max_alt-min_alt),'now','Rotation',90,'BackgroundColor',[1 1 1],'FontWeight','demi');


% % plot clouds
% cho = NaN(size(chm.time))';
% alt = NaN(size(chm.time))';
% for k=1:length(chm.time_ncfiles)
%     ind = find(chm.time>=chm.time_ncfiles(k));
%     if(~isempty(ind))
%         cho(ind) = chm.cho_ncfiles{k};
%         alt(ind) = chm.altitude_ncfiles{k};
%     end
% end
% for i=1:size(chm.cbh,1)
%    cb = chm.cbh(i,:);cb(cb<cho) = NaN;cb = (cb-cho)*sind(ZE) + alt;
% %    cb = chm.cbh(i,:);cb(cb<chm.cho) = NaN;cb = (cb-chm.cho)*sind(ZE) + chm.altitude;
%    hold on;
% %    plot(chm.time-offset,cb,'line','none','marker','o','markersize',5,'markeredgecolor','k','markerfacecolor',0.5*ones(1,3));
%    plot(chm.time-offset,cb,'.','color',0.5*ones(1,3));
% end


    
% xlabel(xlabel_str,'FontWeight','demi');
ylabel(ylabel_str,'FontWeight','demi');
set(gca,'XLim',xlims);
set(gca,'XTick',xticks);
set(gca,'YLim',[min_alt,max_alt]);
set(gca,'YTick',yticks,'YTickLabel',ytickslabels);
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'XMinorTick','on');
box on;grid on;
% Cadre            
set(gca,'layer','top');
% Title
title_str = sprintf('PBL, %s (%d m)',chm.location{end},chm.altitude);  
title(title_str);
set(get(gca,'Title'),'FontWeight','bold')

% Copyright
text(xlims(2)-0.01*diff(xlims),max_alt-(max_alt-min_alt)*0.04,'Issued by MeteoSwiss','BackgroundColor',[1 1 1],'FontWeight','demi','HorizontalAlignment','right'); 



% shows synoptic data
disp('------- Getting synoptic conditions in the Alps -------');
startDate_str  = datestr(xticks(find(xticks>=xlims(1),1))+offset,'yyyymmdd');
finishDate_str = datestr(xticks(find(xticks<=xlims(end),1,'last'))+offset,'yyyymmdd');  
synop = get_synop_from_dwh(startDate_str,finishDate_str);
if(isfield(synop,'t'))
    synop.t = synop.t-offset;
for i=1:length(synop.t);
    if synop.t(i) < xlims(2) 
%         text(synop.t(i)+(xlims(2)-xlims(1))*0.01,max_alt-(max_alt-min_alt)*0.04,synop.synop.desc{i},'BackgroundColor',[0.25 0.75 0],'FontWeight','demi','HorizontalAlignment','left'); 
        text(synop.t(i)+(xlims(2)-xlims(1))*0.01,max_alt-(max_alt-min_alt)*0.04,synop.synop.desc{i},'BackgroundColor',[1 1 1],'FontWeight','demi','HorizontalAlignment','left');    
    end
end
end

ah11 = subplot(5,1,3:4);


min_alt = 0;
max_alt = 5000;
% max_alt = 16000;
xlabel_str = 'Time UT [h]';
ylabel_str = 'Elevation [m asl]';
xlims = [floor(chm.time(1))-offset, ceil(chm.time(end))-offset]; %cf.displayParam.left_right_margin*[-1 1];
xticks = xlims(1)-mod(mod(xlims(1),1)*24,3)/24:2/24:xlims(2);
yticks = [0:500:5000];
% yticks = [0:1000:16000];
ytickslabels = {'0','500','1000','1500','2000','2500','3000','3500','4000','4500','5000'};
% ytickslabels = {'0','1000','2000','3000','4000','5000','6000','7000','8000','9000','10000','11000','12000','13000','14000','15000','16000'};


pcolor(X,Y,C);
colormap(cmap);

shading flat;caxis(clims);
% hcb = colorbar;
% ylabel(hcb,'log10(Normalized range corrected signal)','FontWeight','demi');

% plot current time mark
hold on
plot((now-offset)*[1 1],[min_alt max_alt],'r--')
text((now-offset),-0.055*(max_alt-min_alt),'now','Rotation',90,'BackgroundColor',[1 1 1],'FontWeight','demi');


% plot clouds
cho = NaN(size(chm.time))';
alt = NaN(size(chm.time))';
for k=1:length(chm.time_ncfiles)
    ind = find(chm.time>=chm.time_ncfiles(k));
    if(~isempty(ind))
        cho(ind) = chm.cho_ncfiles{k};
        alt(ind) = chm.altitude_ncfiles{k};
    end
end
for i=1:size(chm.cbh,1)
   cb = chm.cbh(i,:);cb(cb<cho) = NaN;cb = (cb-cho)*sind(ZE) + alt;
%    cb = chm.cbh(i,:);cb(cb<chm.cho) = NaN;cb = (cb-chm.cho)*sind(ZE) + chm.altitude;
   hold on;
%    plot(chm.time-offset,cb,'line','none','marker','o','markersize',5,'markeredgecolor','k','markerfacecolor',0.5*ones(1,3));
   hl_cbh = plot(chm.time-offset,cb,'.','color',0.5*ones(1,3));
end


% line(xtime-offset,bR_theta_int+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime-offset,bR_thetav_int+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime-offset,pm_theta+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime-offset,pm_thetav+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime_sound-offset,bR_theta_woint_sound+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime_sound-offset,bR_thetav_woint_sound+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime_sound-offset,pm_theta_sound+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);
% line(xtime_sound-offset,pm_thetav_sound+alt(1),'LineStyle','none','Marker','d','Color','w','MarkerFaceColor','w','MarkerSize',12);

% hold on;
% hp = plot(...
% xtime-offset,bR_theta_int+alt(1),'*',...
% xtime-offset,bR_thetav_int+alt(1),'*',...
% xtime-offset,pm_theta+alt(1),'.',...
% xtime-offset,pm_thetav+alt(1),'.',...
% xtime_sound-offset,bR_theta_woint_sound+alt(1),'^','MarkerFaceColor','m',...
% xtime_sound-offset,bR_thetav_woint_sound+alt(1),'^',...
% xtime_sound-offset,pm_theta_sound+alt(1),'o',...
% xtime_sound-offset,pm_thetav_sound+alt(1),'o');
% legend(hp,{'bR\_theta\_int';'bR\_thetav\_int';...
%     'pm\_theta';'pm\_thetav';'bR\_theta\_woint\_sound';'bR\_thetav\_woint\_sound';'pm\_theta\_sound';'pm\_thetav\_sound'});
% 


% hl_pm_theta_sound = line(xtime_sound-offset,pm_theta_sound+alt(1),'linestyle','none','marker','o','MarkerFaceColor','r','markersize',8);
% hl_bR_theta_woint_sound = line(xtime_sound-offset,bR_theta_woint_sound+alt(1),'linestyle','none','marker','s','MarkerFaceColor','g','markersize',8);
% hl_pm_thetav_sound = line(xtime_sound-offset,pm_thetav_sound+alt(1),'linestyle','none','marker','o','MarkerFaceColor',[255 204 229]/255,'markersize',8);
% hl_bR_thetav_woint_sound = line(xtime_sound-offset,bR_thetav_woint_sound+alt(1),'linestyle','none','marker','s','MarkerFaceColor','y','markersize',8);
% hl_pm_theta = line(xtime-offset,pm_theta+alt(1),'linestyle','none','marker','o','MarkerFaceColor','r','markersize',4);
% hl_bR_theta_int = line(xtime-offset,bR_theta_int+alt(1),'linestyle','none','marker','s','MarkerFaceColor','g','markersize',4);
% hl_pm_thetav = line(xtime-offset,pm_thetav+alt(1),'linestyle','none','marker','o','MarkerFaceColor',[255 204 229]/255,'markersize',4);
% hl_bR_thetav_int = line(xtime-offset,bR_thetav_int+alt(1),'linestyle','none','marker','s','MarkerFaceColor','y','markersize',4);
% legend([hl_pm_theta_sound;hl_bR_theta_woint_sound;hl_pm_thetav_sound;hl_bR_thetav_woint_sound;hl_pm_theta;hl_bR_theta_int;hl_pm_thetav;hl_bR_thetav_int;hl_cbh],...
%     {'pm\_theta\_snd';'bR\_theta\_snd';'pm\_thetav\_snd';'bR\_thetav\_snd';'pm\_theta';'bR\_theta';'pm\_thetav';'bR\_thetav';'cloud base'},'Location','EastOutside');

% ,'Location','NorthOutside','Orientation','horizontal'

legend([hl_cbh],{'cloud base'},'Location','EastOutside');


    % xlabel(xlabel_str,'FontWeight','demi');
ylabel(ylabel_str,'FontWeight','demi');
set(gca,'XLim',xlims);
set(gca,'XTick',xticks);
set(gca,'YLim',[min_alt,max_alt]);
set(gca,'YTick',yticks,'YTickLabel',ytickslabels);
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'XMinorTick','on');
box on;grid on;
% Cadre            
set(gca,'layer','top');

ah2 = subplot(5,1,5);

hold on;            
plot(xlims, [0 0], 'k-');

disp('------- Getting station data -------');
if(strcmp(stn,'kse'))
    stn_data = get_bekse_from_dwh(start_str,end_str);
else
    stn_data = get_surf_from_dwh(stn,start_str,end_str);
end

if(~isempty(stn_data))

% precipitation logical vector
if(isfield(stn_data,'precip'))
    is_rainy = zeros(length(stn_data.precip),1);
    is_rainy(stn_data.precip>0) = 5;
    is_rainy(is_rainy==0) = NaN;
    for i=2:length(is_rainy)-1
        if isnan(is_rainy(i)) && (is_rainy(i+1)>0 || is_rainy(i-1)>0)
            is_rainy(i) = 0;
        end 
    end
end

if(isfield(stn_data,'T'))
    ground_T_grad=[NaN;diff(stn_data.T)*12]; % x6 to have [K/2h]
end


% Monin-Obukhov (param \in [-0.01;0.01] in data.station)
if strcmp(chm.location{end},'pay')
    z=30; % altitude of the device 
    M_O_length=z./stn_data.turb.M_O_stability_param;
    
    M_O_length=M_O_length/300;

    M_O_length_top_plateau = M_O_length;
    M_O_length_top_plateau(M_O_length<10) = NaN;
    M_O_length_top_plateau(M_O_length>10) = 10;
    
    M_O_length_bot_plateau = M_O_length;
    M_O_length_bot_plateau(M_O_length>-10) = NaN;
    M_O_length_bot_plateau(M_O_length<-10) = -10;
    
    M_O_length_mid = M_O_length;
    for i=1:length(M_O_length)-1
        if M_O_length_mid(i)>10 && M_O_length_mid(i+1)<10
            M_O_length_mid(i)=10;
        elseif M_O_length_mid(i)<-10 && M_O_length_mid(i+1)>-10
            M_O_length_mid(i)=-10;
        elseif M_O_length_mid(i)<10 && M_O_length_mid(i+1)>10
            M_O_length_mid(i+1)=10;
        elseif M_O_length_mid(i)>-10 && M_O_length_mid(i+1)<-10
            M_O_length_mid(i+1)=-10;
        end
    end
    M_O_length_mid(M_O_length_mid>10 | M_O_length_mid<-10) = NaN;
    
    
    % vertical heat flux
    v_heat_flux = stn_data.turb.vert_heat_flux;
    
    v_heat_flux_top_plateau = v_heat_flux;
    v_heat_flux_top_plateau(v_heat_flux<10) = NaN;
    v_heat_flux_top_plateau(v_heat_flux>10) = 10;
    
    v_heat_flux_bot_plateau = v_heat_flux;
    v_heat_flux_bot_plateau(v_heat_flux>-10) = NaN;
    v_heat_flux_bot_plateau(v_heat_flux<-10) = -10;
    
    v_heat_flux_mid = v_heat_flux;
    for i=1:length(v_heat_flux)-1
        if v_heat_flux_mid(i)>10 && v_heat_flux_mid(i+1)<10
            v_heat_flux_mid(i)=10;
        elseif v_heat_flux_mid(i)<-10 && v_heat_flux_mid(i+1)>-10
            v_heat_flux_mid(i)=-10;
        elseif v_heat_flux_mid(i)<10 && v_heat_flux_mid(i+1)>10
            v_heat_flux_mid(i+1)=10;
        elseif v_heat_flux_mid(i)>-10 && v_heat_flux_mid(i+1)<-10
            v_heat_flux_mid(i+1)=-10;
        end
    end
    v_heat_flux_mid(v_heat_flux_mid>10 | v_heat_flux_mid<-10) = NaN;

end


list_plots=[];
list_legends={};
if(isfield(stn_data,'sunshine'))
list_plots(end+1) = bar(stn_data.t-offset,stn_data.sunshine,'FaceColor','y','EdgeColor','k');
list_legends{end+1} = 'sunshine [min/10min]';
end

if(isfield(stn_data,'T'))
% list_plots(end+1) = plot(stn_data.t-offset,ground_T_grad,'color','b','line','-', 'linewidth',2);
% list_legends{end+1} = '\Delta T [K/2h]';
end

if(isfield(stn_data,'precip') && nansum(stn_data.precip)>0)
list_plots(end+1) = plot(stn_data.t-offset,is_rainy,'color','k','line','-','linewidth',2);
list_legends{end+1} = 'precip [y/n]';
end

if strcmp(chm.location{end},'pay')
% plot(stn_data.t-offset,M_O_length/300,'color','m','line','-','linewidth',2);
% list_plots(end+1) = plot(stn_data.t-offset,M_O_length_mid,'color','m','line','-','linewidth',2);
% plot(stn_data.t-offset,M_O_length_bot_plateau,'color','m','line','-.','linewidth',2);
% plot(stn_data.t-offset,M_O_length_top_plateau,'color','m','line','-.','linewidth',2);
% list_legends{end+1} = 'M-O-length';
list_plots(end+1) = plot(stn_data.t-offset,v_heat_flux_mid,'color','r','line','-','linewidth',2);
plot(stn_data.t-offset,v_heat_flux_bot_plateau,'color','r','line','-.','linewidth',2);
plot(stn_data.t-offset,v_heat_flux_top_plateau,'color','r','line','-.','linewidth',2); 
list_legends{end+1} = 'vert. heat flux [W/m^2]';
end

legend(list_plots,list_legends,'Location','EastOutside');

end



ylims = [-11,11];


% plot sunrise & sunset
days = unique(floor(chm.time));
for j=1:length(days)
rs = suncycle(chm.latitude,chm.longitude,days(j));
sunrise = days(j)+rs(1)/24;
sunset = days(j)+rs(2)/24;

hold on
    plot((sunrise-offset)*[1 1],[ylims(1) ylims(2)],'k--','linewidth',2)
    text((sunrise-offset),-0.925*(ylims(2)-ylims(1)),'sunrise','Rotation',90,'BackgroundColor',[1 1 1],'FontWeight','demi');

    hold on
    plot((sunset-offset)*[1 1],[ylims(1) ylims(2)],'k--','linewidth',2)
    text((sunset-offset),-0.925*(ylims(2)-ylims(1)),' sunset','Rotation',90,'BackgroundColor',[1 1 1],'FontWeight','demi');
end


axis([xlims,ylims]);
set(gca,'XTick',xticks);
set(gca,'YTick',[-10:10:10]);
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'XMinorTick','on');
xlabel('','FontWeight','demi');
ylabel('Surface Obs.','FontWeight','demi');
title('');
set(get(gca,'Title'),'FontWeight','bold')

box on;grid on;

% xlabel
xlabel(xlabel_str,'FontWeight','demi');

% % dates   
% startDate_str  = datestr(xticks(find(xticks>=xlims(1),1))+offset,'yyyymmdd');
% finishDate_str = datestr(xticks(find(xticks<=xlims(end),1,'last'))+offset,'yyyymmdd');  
% 
% % time axis start tag
% th = text(xlims(1),ylims(1)-(ylims(2)-ylims(1))*0.4,startDate_str,'FontWeight','demi');
% % time axis end tag
% th = text(xlims(2),ylims(1)-(ylims(2)-ylims(1))*0.4,finishDate_str,'FontWeight','demi','HorizontalAlignment','right');

% dates
days = unique([floor(chm.time);ceil(chm.time)]);
for j=1:length(days)
    date_str = datestr(days(j),'yyyymmdd');
    % time axis  tag
    th = text(xticks(find(xticks>=days(j)-offset,1)),ylims(1)-(ylims(2)-ylims(1))*0.4,date_str,'FontWeight','demi','HorizontalAlignment','center');
end

title_str = sprintf('PBL, %s (%d m)',chm.location{end},chm.altitude);  
title(title_str);


%# find current position [x,y,width,height]
pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');
pos11 = get(ah11,'Position');
%# set width of first axes equal to second
pos1(3) = pos2(3);
pos11(3) = pos2(3);
set(ah1,'Position',pos1);
set(ah11,'Position',pos11);

% set(gcf,'PaperOrientation','landscape','PaperType','A4')
set(gcf,'PaperPositionMode','auto')
outfilestr = ['daily_plots/',start_str(1:8),'_',stn,'.png'];
% print('-dpng',outfilestr);

% close all;
end
end