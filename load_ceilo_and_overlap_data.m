function [chm,chminfo,ov] = load_ceilo_and_overlap_data(stn,date_yyyymmdd)

    stations_EPROFILE ={'lindenberg','hohenspeissenberg','hamburg','oslo','flesland'};
    wmoind = {'10393','10962','10140','01492','01311'};
    institute = {'DWD','DWD','DWD','METNO','METNO'};
    root_stations_EPROFILE = [prefix_data_pay,'/data/pay/REM/ACQ/E_PROFILE_PILOT/DATA/'];
    
    root_stations_MCH = [prefix_data_pay,'/data/pay/REM/ACQ/CEILO_CHM15k/NetCDF/daily/'];
    
    root_stations_others = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/test_sites/'];
    
%     root_overlap_cfg = 'M:\pay-home\pay\users\poy\My Documents\workYP\lib_overlap\';
    root_overlap_cfg = [prefix_data_pay,'/data/pay/PBL4EMPA/overlap_correction/overlap_functions_Lufft/'];
    
    chm = [];
    chminfo = [];
    ov = [];

    %% load ceilometer data
    if strcmp(stn,'pay') || strcmp(stn,'kse')
%         [chm,chminfo] = readcorrectlyncfile3(stn,date_yyyymmdd,'C:\AllData\');
        [chm,chminfo] = readcorrectlyncfile3(stn,date_yyyymmdd,root_stations_MCH);
        if isempty(chm)
            warning('no data found');return;
        end
    elseif strcmp(stn,'SIRTA')
%         [chm,chminfo] = readcorrectlyncfile3(stn,[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'C:\AllData\SharedData_Maxime\st\');
        [chm,chminfo] = readcorrectlyncfile3(stn,[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],[root_stations_others,'CHM150101/']);
        if isempty(chm)
            warning('no data found');return;
        end
    elseif strcmp(stn,'Granada')
%         [chm,chminfo] = readcorrectlyncfile3(stn,[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'C:\AllData\SharedData_Maxime\gr\');
        [chm,chminfo] = readcorrectlyncfile3(stn,[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],[root_stations_others,'CHM120107/']);
        if isempty(chm)
            warning('no data found');return;
        end
    elseif strcmp(stn,'Lindenberg')
%         [chm,chminfo] = readcorrectlyncfile3('lindenberg',[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'M:\pay-data\data\pay\REM\ACQ\E_PROFILE_PILOT\DATA\CEILINEX\chm15k_100110\');
%         [chm,chminfo] = readcorrectlyncfile3('lindenberg',[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'M:\pay-data\data\pay\REM\ACQ\E_PROFILE_PILOT\DATA\CEILINEX\chm15k_140101\');
%         [chm,chminfo] = readcorrectlyncfile3('lindenberg',[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'M:\pay-data\data\pay\REM\ACQ\E_PROFILE_PILOT\DATA\CEILINEX\chx15k_080082\');
%         [chm,chminfo] = readcorrectlyncfile3('lindenberg',[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'M:\pay-data\data\pay\REM\ACQ\E_PROFILE_PILOT\DATA\CEILINEX\chx15k_lmu\');

        [chm,chminfo] = readcorrectlyncfile3('10393',[date_yyyymmdd,'000000'],[date_yyyymmdd,'235959'],'C:\AllData\SharedData_Maxime\ln\');
        if isempty(chm)
            warning('no data found');return;
        end
    elseif any(strcmp(stn,stations_EPROFILE))
        istn = strcmpi(stn,stations_EPROFILE);
        
        [chm,chminfo] = read_chm15k(wmoind{istn},institute{istn},[date_yyyymmdd '000000'],[datestr(datenum(date_yyyymmdd,'yyyymmdd')+1,'yyyymmdd') '000000'],root_stations_EPROFILE);
        if isempty(chm)
            warning('no data found');return;
        end
        chm.beta_raw = chm.beta_raw';
        chm.beta_raw_hr = chm.beta_raw_hr';
        chm.pbl = chm.pbl';
        chm.pbs = chm.pbs';
        chm.cbh = chm.cbh';
        chm.cbe = chm.cbe';
        chm.cdp = chm.cdp';
        chm.cde = chm.cde';
    else
        warning(['station ' stn ' unknown']);
        return;
    end
    
    for j=1:length(chminfo.Attributes)
        if strcmpi(chminfo.Attributes(j).Name,'serlom')
            serlom = chminfo.Attributes(j).Value;
        end
    end
    
    %% load overlap correction and scaling in cfg file
    list = dir([root_overlap_cfg,serlom,'_*_1024.cfg']);
    if isempty(list)
        warning('missing overlap function');
        ov = NaN(1024,1);
    else
        % handle change of overlap cfg file for TUB120001 (Lindenberg) on the 20120917
        if(strcmp(serlom,'TUB120001') && datenum(date_yyyymmdd,'yyyymmdd')<datenum(2012,09,17))
            fname_overlapfc = list(1).name;
        else
            fname_overlapfc = list(end).name;
        end
        disp(['reading ' [root_overlap_cfg fname_overlapfc]]);
        fid = fopen([root_overlap_cfg fname_overlapfc]);
        ov_cell = textscan(fid, '%f','headerLines',1);
        frewind(fid);
        scaling_cell = textscan(fid, 'scaling: %f',1,'headerLines',4);
        fclose(fid);
        ov = cell2mat(ov_cell);
        scaling = cell2mat(scaling_cell);
    end
    
    %% handle old software version in Granada
    if strcmp(stn,'Granada')
        chm.beta_raw = chm.beta_raw_old .* repmat(chm.stddev',length(chm.range),1) / scaling ./ repmat(ov,1,length(chm.time)) ./ repmat(chm.p_calc',length(chm.range),1) .* repmat(chm.range,1,length(chm.time)).^2;
    end

end