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
varnm1{1}='CloudFraction';
varnm1{2}='LayerBase';
varnm1{3}='LayerTop';
%%%%
varnm2={};
varnm2{1}='RVOD_liq_effective_radius';
varnm2{2}='RVOD_liq_number_conc';
varnm2{3}='RVOD_liq_water_content';
varnm2{4}='RVOD_liq_water_path'; 
varnm2{5}='RVOD_ice_effective_radius';
varnm2{6}='RVOD_ice_number_conc';
varnm2{7}='RVOD_ice_water_content'; 
varnm2{8}='RVOD_ice_water_path'; 
%%%%%%%  LO
varnm2{9}='LO_RVOD_AP_geo_mean_radius';
varnm2{10}='LO_RVOD_AP_number_conc';
varnm2{11}='LO_RVOD_effective_radius';
varnm2{12}='LO_RVOD_number_conc';
varnm2{13}='LO_RVOD_liquid_water_content';
%%%%% IO
varnm2{14}='IO_RVOD_AP_log_geo_mean_diameter';
varnm2{15}='IO_RVOD_AP_log_number_conc';
varnm2{16}='IO_RVOD_effective_radius';
varnm2{17}='IO_RVOD_log_number_conc';
varnm2{18}='IO_RVOD_ice_water_content';
varnm2{19}='RVOD_CWC_status';
varnm2{20}='N_cloudy_bins';  
varnm2{21}='RVOD_ice_phase_fraction'; 
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
for iyyy=2010:2010
yearstr='2010';
yearstr = num2str(iyyy,'%4.4i')
foldpath{1}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\');
foldpath{2}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\');
files{1}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\*2B-GEOPROF-LIDAR_GRANULE_P*.hdf');
files{2}=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\*2B-CWC-RVOD_GRANULE_P*.hdf');
SWATHNAME={};
SWATHNAME{1}='2B-GEOPROF-LIDAR';
SWATHNAME{2}='2B-CWC-RVOD';
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
pathout='X:/Data/Cloudsat/TP_May2Sep/txt/';
rgns{1}='ETP';
rgns{2}='WTP';
ipp=1;
for ipp=1:21
for ir=1:2
  	for ifl=2:2
      varnm=varnm2{ipp};
      fileOut=strcat(pathout,rgns{ir},'_EventsCloudsat_',SWATHNAME{ifl},'_',varnm,'_',yearstr,'.txt');
      if ipp==1
        fileOut1=strcat(pathout,yearstr,rgns{ir},'EventsCloudsat_',SWATHNAME{ifl},'_date.txt');
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
    	file=dir(files{ifl});
    	nf=length(file);
    	for n=1:nf
        	ss=0;
%        	for ix=1:ldr
%          		filename=file(n).name;
%          		ss=strfind(filename,ffnm{ix,ir});
%          		if ss > 0
%                    irain = ix
%         			break;
%          		end
%        	end %%%%  ix  
           ss =1;
           ixx=n;
        	if ss>0 
% Open the HDF-EOS2 Swath File.
        		FILE_NAME = strcat(foldpath{ifl},file(n).name);
        		file_id = hdfsw('open', FILE_NAME, 'rdonly');
            icc=length(foldpath{ifl})+1;
            cldsat_yy=str2num(FILE_NAME(icc:icc+3)) ;%YYYYDDDHHMMSS
            cldsat_dd=str2num(FILE_NAME(icc+4:icc+6)); %(5:7)) ;%YYYYDDDHHMMSS
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
        		rawdata={};
        		valid_range={};
%        		for ip=1:np
%          			for iv=1:nv
%            			if ifl==1
%                            varname=varnm{iv};
%                        else
%                            varname=varnm{ip,iv};
%                        end                        
                  varname=varnm;
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
%                                rawdata{ip,iv,invx,invy}=data_var(invx,invy);
%                            end
%                        end
%                        nvrang=length(valid_range_var(:));
%                        for invrang=1:nvrang
%                            valid_range{ip,iv,invrang}=valid_range_var(invrang);
%                        end
% Apply scale factor according to [1].
%          			end % iv
%        		end % ip
        		hdfsw('detach', swath_id);
        		hdfsw('close', file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
          			if lon(nt)<lone(ir) & lon(nt)>lons(ir)  
            			if lat(nt)<late(ir) & lat(nt)>lats(ir)
%               				for ip=1:np
%                  				for iv=1:nv 
%                    				for nl=1:nbin
%                                        if nbin==1
%                                             if data_var(nt,nl)< valid_range_var(1) & ...
%                      							data_var(nt,nl)> valid_range_var(2)
%  %            if data_cpr(nl,nt)>5
%                        					  data_var(nt,nl)=-9999.0;
%                                             end
%                                        else 
%                      					   if data_var(nl,nt)< valid_range_var(1) & ...
%                      							data_var(nl,nt)> valid_range_var(2)
%  %            if data_cpr(nl,nt)>5
%                        					data_var(nl,nt)=-9999.0;
%                      					   end %%%%  if
%                                        end
%                    				end % nl
%                  				end % iv
%                			end %ip 
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
%                          fprintf(outxt,'%f  ',rain{ir,irain});
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
%                              if data_var(nl,nt)>0
%                                  a=data_var(nl,nt)
%                              pause(5)
%                               end
                            end
                          end
%                  				for ip=1:np
%                    				for iv=1:nv 
%                    					fprintf(outxt,'%f  ',);
%                					end
%            	  				end
            	  				fprintf(outxt,'%s\n','');
                                if ipp==1  
                                  fprintf(outxt1,'%s\n','');
                                end
                           end %%NL
            			end % if lat
          			end % if lon
        		end % nt
      		end %%%% end ss
     	end %%%% nt
    	sta = fclose(outxt);
 	end   %%% ifl
    end %%%  ir
   clear data_var
  clear height
  if ipp==1
  sta1 = fclose(outxt1);
  end 
end %%%% ipp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
