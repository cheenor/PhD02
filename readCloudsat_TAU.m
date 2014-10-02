% this code handle the CloudSAT dataset base on the example code
%%% 
% This example code illustrates how to access and visualize CDPC CloudSAT
% Swath file in MATLAB. 
%
%  If you have any questions, suggestions, comments on this example, please 
% use the HDF-EOS Forum (http://hdfeos.org/forums). 
%
%  If you would like to see an  example of any other NASA HDF/HDF-EOS data 
% product that is not listed in the HDF-EOS Comprehensive Examples page 
% (http://hdfeos.org/zoo), feel free to contact us at eoshelp@hdfgroup.org or 
% post it at the HDF-EOS Forum (http://hdfeos.org/forums).

% Tested under: MATLAB R2011b
% Last updated: 2011-11-16

clear
clc 
% list all the fitted dataset
foldpath='X:\Data\Cloudsat\TP\TAU\';
kind='2B-TAU';
namestr=strcat(strcat('*',kind),'_GRANULE_P*.hdf');
file=dir(strcat(foldpath,namestr));
nf=length(file)
%%%%
%     read the envents date 
input=importdata('D:\MyPaper\PhD02\Data\WTP_EventsDate_cloudsat_2010.txt');
evdate=input.data;
ldr=length(evdate(:,1));
ldc=length(evdate(1,:));
%  convert the date to Julian day
ffnm={};
for i0=1:ldr
 y=julia(evdate(i0,1),evdate(i0,2),evdate(i0,3));
 ysrp = num2str(y(1),'%4.4i');
 jsysp = num2str(y(2),'%3.3i');
 ffnm{i0,2}=strcat(ysrp,jsysp);
end
input=importdata('D:\MyPaper\PhD02\Data\ETP_EventsDate_cloudsat_2010.txt');
evdate=input.data;
ldr=length(evdate(:,1));
ldc=length(evdate(1,:));
for i0=1:ldr
 y=julia(evdate(i0,1),evdate(i0,2),evdate(i0,3));
 ysrp = num2str(y(1),'%4.4i');
 jsysp = num2str(y(2),'%3.3i');
 ffnm{i0,1}=strcat(ysrp,jsysp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1=1:125 
  ijx(i1,1)=0.0;
  rf_mean(i1,1)=0.0;
  ijx(i1,2)=0.0;
  rf_mean(i1,2)=0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varname= 'layer_optical_depth'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nf
  for nr=1:2
    ss=0;
  for ix=1:ldr
      filename=file(n).name;
    ss=strfind(filename,ffnm{ix,nr});
     if ss > 0
      break;
     end 
   end
  ixx=n
  if ss>0 
% Open the HDF-EOS2 Swath File.
  FILE_NAME = strcat(foldpath,file(n).name)
  file_id = hdfsw('open', FILE_NAME, 'rdonly')

% Read data.
  SWATH_NAME = kind;
  swath_id = hdfsw('attach', file_id, SWATH_NAME);

% Read lat/lon/height/time data.
  [lon, status] = hdfsw('readfield', swath_id, 'Longitude', [], [], []);
  [lat, status] = hdfsw('readfield', swath_id, 'Latitude', [], [], []);
  [height, status] = hdfsw('readfield', swath_id, 'Height', [], [], []);
  [time, status] = hdfsw('readfield', swath_id, 'Profile_time', [], [], []);

% Make type double for plotting.
  lat=double(lat);
  lon=double(lon);
  time=double(time);


% Read data & attributes.
        DATAFIELD_NAME =varname ;       % 'CloudFraction';
        [data, status] = hdfsw('readfield', swath_id, DATAFIELD_NAME, [],[],[]);
  		
        [long_name, status] = hdfsw('readattr', swath_id, ...
                         strcat(varname,'.long_name'));
        [units, status] = hdfsw('readattr', swath_id, ...
                       strcat(varname,'.units'));
        [scale_factor, status] = hdfsw('readattr', swath_id, ...
                               strcat(varname,'.factor'));
        scale_factor = double(scale_factor);

        [valid_range, status] = hdfsw('readattr', swath_id, ...
                             strcat(varname,'.valid_range'));
%
  data=double(data);
%
  [units_h, status] = hdfsw('readattr', swath_id, ...
                       'Height.units');

  [units_t, status] = hdfsw('readattr', swath_id, ...
                       'Profile_time.units');
  [long_name_t, status] = hdfsw('readattr', swath_id, ...
                        'Profile_time.long_name');

  hdfsw('detach', swath_id);
  hdfsw('close', file_id);

% Process valid_range. Fill value and missing value will be handled by this
% since they are outside of range values.

%  data_rf((data_rf < valid_range_rf(1)) | (data_rf > valid_range_rf(2))) = NaN;
%  data_cpr((data_cpr < valid_range_cpr(1)) | (data_cpr > valid_range_cpr(2))) = NaN;

% Apply scale factor according to [1].
  data = data / scale_factor;
% Apply scale factor according to [1].
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  selected for the region 
%%%  ;;  ETP  lon 90-100   lat  30 37.5
%%%%     WTP  lon 80-90   lat  30 37.5
  lons(1)=90.0;
  lons(2)=80.0;
  lats(1)=30.0;
  lats(2)=30.0;
  lone(1)=100.0;
  lone(2)=90.0;
  late(1)=37.5;
  late(2)=37.5;
    nray=length(lat);
    nbin=length(height(:,1));  
    for nt=1:nray
        if lon(nt)<lone(nr) & lon(nt)>lons(nr)  
          if lat(nt)<late(nr) & lat(nt)>lats(nr)
           for nl=1:nbin
           if data(nl,nt)> valid_range(1) & data(nl,nt)< valid_range(2)
  %            if data_cpr(nl,nt)>5
                rf_mean(nl,nr)=rf_mean(nl,nr)+data(nl,nt);
                ijx(nl,nr)=ijx(nl,nr)+1.0;
 %            end % data_cpr
           end %data_rf
          end % nl
        end % if lat
      end % if lon
    end % nt
   end %% nr
  end %%%% end ss
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:2
 for ix=1:nbin
  if ijx(ix,k) ~= 0
   rf_mean(ix,k)=rf_mean(ix,k)/ijx(ix,k);
  end
 end 
end 
%%%% save the data to txt file
fileOut=strcat(strcat('D:\MyPaper\PhD02\Data\EventsProfile_cloudsat_',kind),'.txt');
outxt=fopen(fileOut,'w');
fprintf(outxt,'%s  ','Height');
fprintf(outxt,'%s  ','ETP');
fprintf(outxt,'%s\n','WTP');
for nl=1:nbin
   fprintf(outxt,'%f  ',height(nl,1));
   fprintf(outxt,'%f  ',rf_mean(nl,1));
   fprintf(outxt,'%f\n',rf_mean(nl,2));
end
sta = fclose(outxt);

%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
