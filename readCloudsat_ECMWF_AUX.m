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
foldpath='X:\Data\Cloudsat\TP\ECMWF_AUX\';
kind='ECMWF-AUX';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varnm={};
varnm{1,1}='Pressure';
varnm{1,2}='Temperature';
varnm{1,3}='Specific_humidity';
vnp=1;
vnv=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nl=1:125  %nl,iv,ip,nr
   for ip=1,vnp
	 for iv=1:vnv
	  for nr=1:2
         ijx(nl,iv,ip,nr)=0.0;
         rf_mean(nl,iv,ip,nr)=0.0;
      end
     end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  ixx=n;
  if ss>0 
% Open the HDF-EOS2 Swath File.
  FILE_NAME = strcat(foldpath,file(n).name)
  file_id = hdfsw('open', FILE_NAME, 'rdonly')

  SWATH_NAME = kind;
  swath_id = hdfsw('attach', file_id, SWATH_NAME);

% Read lat/lon/height/time data.
  [lon, status] = hdfsw('readfield', swath_id, 'Longitude', [], [], []);
  [lat, status] = hdfsw('readfield', swath_id, 'Latitude', [], [], []);
  [height, status] = hdfsw('readfield', swath_id, 'EC_height', [], [], []);
  [time, status] = hdfsw('readfield', swath_id, 'Profile_time', [], [], []);

% Make type double for plotting.
  lat=double(lat);
  lon=double(lon);
  time=double(time);

% Read attributes.
  for ip=1:vnp
  	for iv=1:vnv  		
  		varname=varnm{ip,iv};
        % Read data.
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
%
        [units_h, status] = hdfsw('readattr', swath_id, ...
                       'Height.units');

        [units_t, status] = hdfsw('readattr', swath_id, ...
                       'Profile_time.units');
        [long_name_t, status] = hdfsw('readattr', swath_id, ...
                        'Profile_time.long_name');
        data=double(data);
% Process valid_range. Fill value and missing value will be handled by this
% since they are outside of range values.

%  data_rf((data_rf < valid_range_rf(1)) | (data_rf > valid_range_rf(2))) = NaN;
%  data_cpr((data_cpr < valid_range_cpr(1)) | (data_cpr > valid_range_cpr(2))) = NaN;

% Apply scale factor according to [1].
   if isempty(scale_factor)
       scale_factor=1.
   end
  data= data / scale_factor;
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
              if isempty(valid_range)
                rf_mean(nl,iv,ip,nr)=rf_mean(nl,iv,ip,nr)+data(nl,nt);
                ijx(nl,iv,ip,nr)=ijx(nl,iv,ip,nr)+1.0;
              else
                if data(nl,nt)> valid_range(1) & data(nl,nt)< valid_range(2)
  %            if data_cpr(nl,nt)>5
                rf_mean(nl,iv,ip,nr)=rf_mean(nl,iv,ip,nr)+data(nl,nt);
                ijx(nl,iv,ip,nr)=ijx(nl,iv,ip,nr)+1.0;
 %             end % data_cpr
                end %data_rf
             end
            end % nl
          end % if lat
        end % if lon
      end % nt
     end % iv
    end % ip
    clear data
   hdfsw('detach', swath_id);
   hdfsw('close', file_id);
   end %% end ss 
  end %%%% nr
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ip= 1:vnp
 for iv=1:vnv
  for k=1:2
    for ix=1:nbin
     if ijx(ix,iv,ip,k) ~= 0
       rf_mean(ix,iv,ip,k)=rf_mean(ix,iv,ip,k)/ijx(ix,iv,ip,k);
     end
    end 
   end 
 end
end
%%%% save the data to txt file
fileOut=strcat(strcat('D:\MyPaper\PhD02\Data\EventsProfile_cloudsat_',kind),'.txt');
outxt=fopen(fileOut,'w');
fprintf(outxt,'%s\n','ETP');
fprintf(outxt,'%s  ','Height');
for k=1:vnp
	for iv=1:vnv   %vnp=1;vnv=3;
       fprintf(outxt,'%s  ',varnm{k,iv});
   end
end
fprintf(outxt,'%s\n','');
for k=1:2
	if k==2
		fprintf(outxt,'%s\n','WTP');
	end
   for nl=1:nbin
    fprintf(outxt,'%f  ',height(nl,1));
      for ip=1:vnp
        for iv=1:vnv
          fprintf(outxt,'%f  ',rf_mean(nl,iv,ip,k));
        end
      end
    fprintf(outxt,'%s\n','');
   end
end
sta = fclose(outxt);

%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
