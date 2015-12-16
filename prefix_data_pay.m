function prefix = prefix_data_pay
% PREFIX_DATA_PAY returns prefix to be used on either Windows or Unix machines for paths on the network beginning with /data/pay/
% EXAMPLE(S):
%       fpath = [prefix_data_pay,'/data/pay/REM/ACQ/CEILO_CHM15k/NetCDF/daily/2015/10/20151025_pay_CHM120106_000.nc'];
%       ncinfo(fpath)
% COPYRIGHT: MeteoSwiss, 2015
% VERSION(S): 10.12.2015, Yann.Poltera@meteoswiss.ch

if ispc
    prefix = '//meteoswiss.ch/mch/pay-data';
else 
    prefix = '';
end

end