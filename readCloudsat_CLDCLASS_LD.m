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
foldpath='X:\Data\Cloudsat\TP\CLASS_LD\';
kind='2B-CLDCLASS-LIDAR';
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
ncloud=10
for nl=1:ncloud  %nl,iv,ip,nr
	 for iv=1:13
	  for nr=1:2
         ijx(nl,iv,nr)=0.0;
         rf_mean(nl,iv,nr)=0.0;
      end
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vrname={}
vrname{1}= 'Cloudlayer'; % nray 
vrname{2}='CloudLayerBase';
vrname{3}='LayerBaseFlag';
vrname{4}='CloudLayerTop';
vrname{5}='LayerTopFlag';
vrname{6}='CloudFraction';
vrname{7}='CloudPhase';
vrname{8}='CloudPhaseConfidenceLevel';
vrname{9}='CloudLayerType';
vrname{10}='CloudTypeQuality';
vrname{11}='PrecipitationFlag';
vrname{12}='Phase_log';
vrname{13}='Water_layer_top'
vrnm=13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icl=1
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
    for iv=1:vrnm
        varname=vrname{iv}
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
    nbin=length(data(:,1));  
    for nt=1:nray
        if lon(nt)<lone(nr) & lon(nt)>lons(nr)  
          if lat(nt)<late(nr) & lat(nt)>lats(nr)
            if iv >1 
             for nl=1:nbin
                if isempty(valid_range)
                rf_mean(nl,iv,nr)=rf_mean(nl,iv,nr)+data(nl,nt);
                ijx(nl,iv,nr)=ijx(nl,iv,nr)+1.0;
              else
                if data(nl,nt)> valid_range(1) & data(nl,nt)< valid_range(2)
  %            if data_cpr(nl,nt)>5
                rf_mean(nl,iv,nr)=rf_mean(nl,iv,nr)+data(nl,nt);
                ijx(nl,iv,nr)=ijx(nl,iv,nr)+1.0;
 %             end % data_cpr
                end %data_rf
             end
             end % nl
            end %%% if iv
           if iv==1 
            cloudlayer(icl,nr,:)=data;
            icl=icl+1
           end
        end % if lat
      end % if lon
    end % nt
   end %iv
    hdfsw('detach', swath_id);
    hdfsw('close', file_id);
   end %% end ss
  end %%%% nr
end %%% nf
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
      for iv=2:vrnm
          fprintf(outxt,'%f  ',rf_mean(nl,iv,k));
       end
   fprintf(outxt,'%f\n');
end
sta = fclose(outxt);
fileOut=strcat(strcat('D:\MyPaper\PhD02\Data\EventsProfile_cloudsat_',kind),'ETP_cldLayer.txt');
outxt=fopen(fileOut,'w');
fprintf(outxt,'%f  ',cloudlayer(:,1,:));
sta = fclose(outxt);
fileOut=strcat(strcat('D:\MyPaper\PhD02\Data\EventsProfile_cloudsat_',kind),'WTP_cldLayer.txt');
outxt=fopen(fileOut,'w');
fprintf(outxt,'%f  ',cloudlayer(:,2,:));
sta = fclose(outxt);
%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
