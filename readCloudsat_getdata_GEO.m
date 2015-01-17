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

%%%%
varnm1={};
varnm1{1}='Radar_Reflectivity';
varnm1{3}='CPR_Echo_Top';   %No Determination 0 Clear Profile 1 No hydrometeor layers found in CPR profile
                            %High Cloud 2 Echo top pressure less than 500 mb.;
                            %Mid Level cloud 3 Echo top pressure greater than 500 mb and echo top temperature colder than 273.
                            %Low Level cloud 4 Echo top pressure greater than 500 mb and echo top temperature warmer than 273.
                            %Multi Layer 5 Distinct combinations of above types 
varnm1{2}='CPR_Cloud_mask';
%%%%
%     read the envents date 
input=importdata('D:\MyPaper\PhD02\Data\WTP_EventsDate_cloudsat_2010.txt');
evdate=input.data;
ldr=length(evdate(:,1));
ldc=length(evdate(1,:));
%  convert the date to Julian day
ffnm={};
rain={};
for i0=1:ldr
 y=julia(evdate(i0,1),evdate(i0,2),evdate(i0,3));
 ysrp = num2str(y(1),'%4.4i');
 jsysp = num2str(y(2),'%3.3i');
 ffnm{i0,2}=strcat(ysrp,jsysp);
 rain{1,i0}=evdate(i0,4);
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
 rain{2,i0}=evdate(i0,4);
end
days=[];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1=1:125 
  ijx(i1,1)=0.0;
  rf_mean(i1,1)=0.0;
  ijx(i1,2)=0.0;
  rf_mean(i1,2)=0.0;
end
files={};
foldpath={};
yearstr='2010';
SWATHNAME={};
SWATHNAME{1}='2B-GEOPROF-LIDAR';
SWATHNAME{2}='2B-CWC-RVOD';
SWATHNAME{3}='2B-GEOPROF';
for iys=2006:2010
yearstr= num2str(iys,'%4.4i');  %%%X:\Data\Cloudsat\TP_May2Sep\GEO\2007
foldpath{1}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\');
foldpath{2}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\');
foldpath{3}=strcat('X:\Data\Cloudsat\TP_May2Sep\GEO\',yearstr,'\');
files{1}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\*2B-GEOPROF-LIDAR_GRANULE_P*.hdf');
files{2}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\*2B-CWC-RVOD_GRANULE_P*.hdf');
files{3}=strcat('X:\Data\Cloudsat\TP_May2Sep\GEO\',yearstr,'\*2B-GEOPROF_GRANULE_P*.hdf');
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
pathout='D:\MyPaper\PhD02\Data\';
rgns{1}='ETP';
rgns{2}='WTP';
varnm={};
for ipp=1:3
for ir=1:2
      for ifl=3:3
%      if ifl==1
%        varnm=varnm1
%        np=1;
%        nv=length(varnm);
%      else
%        varnm=varnm2
%        np=length(varnm(:,1));
%        nv=length(varnm(1,:));
%     end
        fileOut=strcat(pathout,rgns{ir},'_EventsCloudsat_',SWATHNAME{ifl},'_',varnm1{ipp},'_',yearstr,'.txt');
        outxt=fopen(fileOut,'w');
     fprintf(outxt,'%s  ','year');
     fprintf(outxt,'%s  ','month');
     fprintf(outxt,'%s  ','day');
     fprintf(outxt,'%s  ','hour');
     fprintf(outxt,'%s  ','minute');                            	
     fprintf(outxt,'%s  ','lon');
     fprintf(outxt,'%s  ','lat'); 
      fprintf(outxt,'%s  ','height');
%      fprintf(outxt,'%s  ','lon');
%      fprintf(outxt,'%s  ','lat');
%      fprintf(outxt,'%s  ','rain');
%       for ip=1:np
%          for iv=1:nv 
            fprintf(outxt,'%s  ',varnm1{ipp});
%          end
%       end
      fprintf(outxt,'%s\n','');
        file=dir(files{ifl});
        nf=length(file);
        for n=1:nf
            ss=0;
 %           for ix=1:ldr
 %                 filename=file(n).name;
 %                 ss=strfind(filename,ffnm{ix,ir});
 %                 if ss > 0
 %               irain=ix
 %                      break;
 %                 end
 %           end %%%%  ix  
 %           ixx=n;
            ss=1;
            if ss>0 
% Open the HDF-EOS2 Swath File.
                FILE_NAME = strcat(foldpath{ifl},file(n).name);
                file_id = hdfsw('open', FILE_NAME, 'rdonly');
                icc=length(foldpath{1})+1;
                cldsat_yy=str2num(FILE_NAME(icc:icc+3)) ;%YYYYDDDHHMMSS
                cldsat_dd=str2num(FILE_NAME(icc+4:icc+6)) ;%(5:7)) ;%YYYYDDDHHMMSS
                cldsat_hh=str2num(FILE_NAME(icc+7:icc+8));  %(8:9)) ;%YYYYDDDHHMMSS
                cldsat_mm=str2num(FILE_NAME(icc+9:icc+10));  %(10:11)) ;%YYYYDDDHHMMSS
                cldsat_ss=str2num(FILE_NAME(icc+11:icc+12)); %(12:13)) ;%YYYYDDDHHMMSS
% Read data.
                SWATH_NAME =SWATHNAME{ifl} ;
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
                [units_h, status] = hdfsw('readattr', swath_id, ...
                       'Height.units');

                [units_t, status] = hdfsw('readattr', swath_id, ...
                       'Profile_time.units');
                [long_name_t, status] = hdfsw('readattr', swath_id, ...
                        'Profile_time.long_name');
% Read attributes.
%                rawdata={};
%                valid_range={};
%                for ip=1:np
%                      for iv=1:nv
%                        if ifl==1
%                            varname=varnm{iv};
%                        else
%                            varname=varnm{ip,iv};
%                        end 
                        varname=varnm1{ipp};
                        DATAFIELD_NAME = varname;
                        [data_var, status] = hdfsw('readfield', swath_id, DATAFIELD_NAME, [],[],[]);
                        [long_name_var, status] = hdfsw('readattr', swath_id, ...
                             strcat(varname,'.long_name'));
                        [units_fc, status] = hdfsw('readattr', swath_id, ...
                               strcat(varname,'.units'));
                        [scale_factor_var, status] = hdfsw('readattr', swath_id, ...
                            strcat(varname,'.factor'));
                        scale_factor_var = double(scale_factor_var);

                        [valid_range_var, status] = hdfsw('readattr', swath_id, ...
                            strcat(varname,'.valid_range'));
%
%
                        data_var=double(data_var);
% Process valid_range. Fill value and missing value will be handled by this
% since they are outside of range values.

%  data_rf((data_rf < valid_range_rf(1)) | (data_rf > valid_range_rf(2))) = NaN;
%  data_cpr((data_cpr < valid_range_cpr(1)) | (data_cpr > valid_range_cpr(2))) = NaN;

% Apply scale factor according to [1].
                        data_var = data_var / scale_factor_var;
%                        nvx=length(data_var(:,1));
%                        nvy=length(data_var(1,:));
%                        for invx=1:nvx
%                            for invy=1:nvy
%                                rawdata{iv,invx,invy}=data_var(invx,invy);
%                            end
%                        end
%                        nvrang=length(valid_range_var(:));
%                        for invrang=1:nvrang
%                            valid_range{iv,invrang}=valid_range_var(invrang);
%                        end
% Apply scale factor according to [1].
%
%                      end % iv
%                end % ip
                hdfsw('detach', swath_id);
                hdfsw('close', file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nray=length(lat);
                nbin=length(height(:,1));
                nray=length(data_var(1,:))
                nbin=length(data_var(:,1));
                if nray==1 & nbin>124
                    ntmp=nray
                    nray=nbin;
                    nbin=ntmp
                end  
                for nt=1:nray
                      if lon(nt)<lone(ir) & lon(nt)>lons(ir)  
                        if lat(nt)<late(ir) & lat(nt)>lats(ir)
%                               for ip=1:np
%                                  for iv=1:nv 
                                    for nl=1:nbin
                                        if nbin==1
                                             if data_var(nt,nl)< valid_range_var(1) & ...
                      							data_var(nt,nl)> valid_range_var(2)
  %            if data_cpr(nl,nt)>5
                        					  data_var(nt,nl)=-9999.0;
                                             end
                                        else    
                      					   if data_var(nl,nt)< valid_range_var(1) & ...
                      							data_var(nl,nt)> valid_range_var(2)
  %            if data_cpr(nl,nt)>5
                        					data_var(nl,nt)=-9999.0;
                      					   end %%%%  if
                                        end
                                    end % nl
%                                  end % iv
%                            end %ip 
%%%%%%%%%%%%%%%%%%%%%%%  verify the time 
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
                                     tempday1=tempday1+days(im) ;
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
                         fprintf(outxt,'%d  ',cldsat_yyd);
                         fprintf(outxt,'%d  ',cldsat_mnd);
                         fprintf(outxt,'%d  ',cldsat_ddd);
                         fprintf(outxt,'%d  ',cldsat_hhd);
                         fprintf(outxt,'%d  ',cldsat_mmd);                            	
                         fprintf(outxt,'%f  ',lon(nt));
%                         fprintf(outxt,'%f  ',lat(nt));
                         nxxx=length(data_var(:,1))
                          for nl=1:nbin
                              if nbin==1
                                  fprintf(outxt,'%f  ',nl);
                                   fprintf(outxt,'%f  ',data_var(nt,nl));
                              else    
                                    fprintf(outxt,'%f  ',height(nl,nt));
                                    fprintf(outxt,'%f  ',data_var(nl,nt));
                              end
%                                    fprintf(outxt,'%f  ',rain{ir,irain});
%                                  for ip=1:np
%                                    for ivc=1:nv 
%                                        fprintf(outxt,'%f  ',rawdata{ivc,nl,nt});
%                                    end
%                                    end
                                  fprintf(outxt,'%s\n','');
                            end %%%%% nl
                        end % if lat
                      end % if lon
                end % nt
              end %%%% end ss
         end %%%% nt
        sta = fclose(outxt);
     end   %%% ifl
end %%%  ir
clear data_var
end %%% ipp
end %%%%% iyyy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
