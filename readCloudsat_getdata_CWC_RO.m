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

clear;
clc ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varnmall={};
varnmall{1}='RO_liq_effective_radius';
varnmall{2}='RO_liq_number_conc';
varnmall{3}='RO_liq_water_content';
varnmall{4}='RO_ice_effective_radius';
varnmall{5}='RO_ice_number_conc';
varnmall{6}='RO_ice_effective_radius';
%%%%%%%  LO
varnmall{7}='LO_RO_AP_geo_mean_radius';
varnmall{8}='LO_RO_AP_number_conc';
varnmall{9}='LO_RO_effective_radius';
varnmall{10}='LO_RO_number_conc';
varnmall{11}='LO_RO_liquid_water_content';
%%%%% IO
varnmall{12}='IO_RO_AP_log_geo_mean_diameter';
varnmall{13}='IO_RO_AP_log_number_conc';
varnmall{14}='IO_RO_effective_radius';
varnmall{15}='IO_RO_log_number_conc';
varnmall{16}='IO_RO_ice_water_content';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:12
  days(i)=30;
end
days(2)=28;
days(1)=31;
days(3)=31;
days(5)=31;
days(7)=31;
days(8)=31;
days(12)=31;
lons(1)=90.0;
lons(2)=80.0;
lats(1)=30.0;
lats(2)=30.0;
lone(1)=100.0;
lone(2)=90.0;
late(1)=37.5;
late(2)=37.5;
pathout='K:\DATA\CloudSat\TXT\2B-CWC-RVOD\';
rgns{1}='ETP';
rgns{2}='WTP';
for nl=1:125  %nl,iv,ip,nr
   for ip=1,2 
	 for iv=1:8
	  for nr=1:2
         ijx(nl,iv,ip,nr)=0.0;
         rf_mean(nl,iv,ip,nr)=0.0;
      end
     end
   end
end
for iyyy=2015:2015
  yearstr = num2str(iyyy,'%4.4i')
  % list all the fitted dataset
  foldpath=strcat('K:\DATA\CloudSat\',yearstr,'05-09\'); %'X:\Data\Cloudsat\TP\CWC_RO\'
  file=dir(strcat(foldpath,'*2B-CWC-RO_GRANULE_P*.hdf'));
  nf=length(file);
%%%%
  ipp=1;
  for ipp=1:16
    for ir=1:2
      varnm=varnmall{ipp};
      SWATHNAME = '2B-CWC-RO';   
      fileOut=strcat(pathout,rgns{ir},'_EventsCloudsat_',SWATHNAME,'_',varnm,'_',yearstr,'.txt');
      if ipp==1
        fileOut1=strcat(pathout,yearstr,rgns{ir},'EventsCloudsat_',SWATHNAME,'_date.txt');
        outxt1=fopen(fileOut1,'w');
        fprintf(outxt1,'%s  ','year');
        fprintf(outxt1,'%s  ','month');
        fprintf(outxt1,'%s  ','day');
        fprintf(outxt1,'%s  ','hour');
        fprintf(outxt1,'%s  ','minute');                              
        fprintf(outxt1,'%s  ','lon');
        fprintf(outxt1,'%s  ','lat'); 
        fprintf(outxt1,'%s  ','height');
        fprintf(outxt1,'%s\n','');
      end
      outxt=fopen(fileOut,'w');
      fprintf(outxt,'%s  ',varnm);
      fprintf(outxt,'%s\n','');
      for n=1:nf      
        ss=0;
        %for ix=1:ldr
        %  filename=file(n).name;
        %  ss=strfind(filename,ffnm{ix,nr});
        %  if ss > 0
        %    break;
        %  end    
        %end
        ss=1;
        ixx=n;
        if ss>0 
          % Open the HDF-EOS2 Swath File.
          FILE_NAME = strcat(foldpath,file(n).name);
          file_id = hdfsw('open', FILE_NAME, 'rdonly');
          icc=length(foldpath)+1;
          cldsat_yy=str2num(FILE_NAME(icc:icc+3)) ;%YYYYDDDHHMMSS
          cldsat_dd=str2num(FILE_NAME(icc+4:icc+6)); %(5:7)) ;%YYYYDDDHHMMSS
          cldsat_hh=str2num(FILE_NAME(icc+7:icc+8));  %(8:9)) ;%YYYYDDDHHMMSS
          cldsat_mm=str2num(FILE_NAME(icc+9:icc+10));  %(10:11)) ;%YYYYDDDHHMMSS
          cldsat_ss=str2num(FILE_NAME(icc+11:icc+12)); %(12:13)) ;%YYYYDDDHHMMSS
          SWATH_NAME = '2B-CWC-RO';
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
          % Read attributes. 		  		
  		    varname=varnmall{ipp};
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
          data_var= data / scale_factor;
% Apply scale factor according to [1].
          hdfsw('detach', swath_id);
          hdfsw('close', file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  selected for the region 
%%%  ;;  ETP  lon 90-100   lat  30 37.5
%%%%     WTP  lon 80-90   lat  30 37.5
          nray=length(lat);
          nbin=length(height(:,1));
          nray=length(data_var(1,:));
          nbin=length(data_var(:,1)); 
          if nray==1 & nbin>124
            ntmp=nray;
            nray=nbin;
            nbin=ntmp;
          end    
          for nt=1:nray
            if lon(nt)<lone(nr) & lon(nt)>lons(nr)  
              if lat(nt)<late(nr) & lat(nt)>lats(nr)
                dyear=365;
                days(2)=28;
                if mod(cldsat_yy,4)==0 & mod(cldsat_yy,100) ~=0
                  dyear=366;
                  days(2)=29;
                elseif mod(cldsat_yy, 400)==0
                  dyaer=366;
                  days(2)=29;
                end                                        
                sss=cldsat_ss+time(nt);
                amin=sss/60.;
                amin=fix(amin);
                smin= cldsat_mm+amin ;
                cldsat_mmd=mod(smin,60);
                ahh = smin/60.;
                ahh = fix(ahh);
                shh =ahh + cldsat_hh;
                cldsat_hhd=mod(shh,24);
                add = shh/24.;
                add = fix(ahh);
                sdd =add + cldsat_dd;
                cldsat_ddd=mod(sdd,dyear); %366
                ayy=sdd/(dyear*1.0);
                ayy=fix(ayy);
                cldsat_yyd=cldsat_yy+ayy;
%%%%   the following coed convert  Julian day to date MM DD
                julday=cldsat_ddd;
                tempday1=0;
                tempday2=0;
                for im=1:11                                    
                  tempday1=tempday1+days(im); 
                  tempday2=tempday1+days(im+1);
                  if julday-tempday1 > 0 & julday-tempday2<0
                    mm=im+1;
                    dd=julday-tempday1;
                    break;
                  elseif julday-tempday1<1 & im==1
                    mm=im;
                    dd= julday;
                    break;
                  end
                end 
                cldsat_mnd=mm;
                cldsat_ddd=dd; 
                if ipp==1
                  fprintf(outxt1,'%d  ',cldsat_yyd);
                  fprintf(outxt1,'%d  ',cldsat_mnd);
                  fprintf(outxt1,'%d  ',cldsat_ddd);
                  fprintf(outxt1,'%d  ',cldsat_hhd);
                  fprintf(outxt1,'%d  ',cldsat_mmd);                             
                  fprintf(outxt1,'%f  ',lon(nt));
                  fprintf(outxt1,'%f  ',lat(nt));
                end
              %  fprintf(outxt,'%f  ',rain{ir,irain});
                nxxx=length(data_var(:,1));
                for nl=1:nbin
                  if nbin==1
                    fprintf(outxt,'%f  ',nl);
                    fprintf(outxt,'%f  ',data_var(nt,nl));  
                  else
                    if nxxx<nbin
                      fprintf(outxt,'%f  ',nl);
                      fprintf(outxt,'%f  ',data_var(nl,nt));
                    else
                      if ipp==1  
                        fprintf(outxt1,'%f  ',height(nl,nt));
                      end
                      fprintf(outxt,'%f  ',data_var(nl,nt));
%                             if data_var(nl,nt)>0
%                                 a=data_var(nl,nt)
%                             pause(5)
%                              end
                    end
                  end
                  fprintf(outxt,'%s\n','');
                  if ipp==1  
                    fprintf(outxt1,'%s\n','');
                  end
                end %%NL
              end % if lat
            end % if lon
          end % nt
        end % SS   
      end %% end nf
      sta = fclose(outxt);
    end %%%% nr
    clear data_var
    clear height
    if ipp==1
      sta1 = fclose(outxt1);
    end 
  end %ipp
end % iyyy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
