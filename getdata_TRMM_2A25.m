% this code handles the TRMM 2A25 dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   TRMM   2A25  DATASET
%nscan    var        Number of scans in the granule.
%nray     49         Number of angle bins in each scan.
%ncell1   80         Number of radar range cells at which the rain rate is estimated.
%                    The cells range from 0 to 79. Each cell is 250m apart, with cell 79 at the earth ellipsiod.
%ncell2    5         Number of radar range cells at which the Z-R parameters are output.
%nmeth     2         Number of methods used.
%nestmeth  6         Number of estimation methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　　　　　　　　　
clear
clc 
% 
ncell=80
nscan=49
ncell2=5
nmeth=2
nestmeth=6
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
for iys=2006:2010
for ig=1:2    
    yearstr= num2str(iys,'%4.4i');  %%%X:\Data\Cloudsat\TP_May2Sep\GEO\2007
    foldpath=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\');
    files=strcat('X:\Data\Cloudsat\TP_May2Sep\',yearstr,'\*2B-GEOPROF-LIDAR_GRANULE_P*.hdf');
    fileOut=strcat(pathout,rgns{ir},yearstr,'_TRMM_2A25_TimeSequence.txt');
    outxt=fopen(fileOut,'w');
    fprintf(outxt,'%s  ','year');
    fprintf(outxt,'%s  ','month');
    fprintf(outxt,'%s  ','day');
    fprintf(outxt,'%s  ','hour');
    fprintf(outxt,'%s  ','minute');                            	
    fprintf(outxt,'%s  ','lon');
    fprintf(outxt,'%s  ','lat'); 
    fprintf(outxt,'%s  ','height');
    fprintf(outxt,'%s\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%head=hdfinfo(TRMM2A25file)
% every variable have the branches
%%%%%%   ScanTimeIndex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%Type: 'Scientific Data Set' 
%           Name: 'Minute'
%           Rank: 1
%       DataType: 'int8'
%     Attributes: [1x1 struct]
%           Dims: [1x1 struct]
%          Label: {}
%    Description: {}
%          Index: 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ScanTimeIndex
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,1).Name Year            nray
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,2).Name Month
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,3).Name DayOfMonth
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,4).Name Hour
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,5).Name Minute
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,6).Name Second
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,7).Name MilliSecond
%head.Vgroup(1,2).Vgroup(1,1).SDS(1,8).Name DayOfYear
%%%%%%%   Wanted variables Index   in Vgroup(1,2).SDS(1,1:45)
%           Type: 'Scientific Data Set'
%           Name: 'reliab'
%           Rank: 3
%       DataType: 'int8'
%     Attributes: []
%           Dims: [3x1 struct]
%          Label: {}
%    Description: {}
%          Index: 41
%head.Vgroup(1,2).SDS(1,2).Name Latitude                 nray x nscan
%head.Vgroup(1,2).SDS(1,3).Name Longitude                nray x nscan
%head.Vgroup(1,2).SDS(1,5).Name rain                     ncell x nray x nscan
%head.Vgroup(1,2).SDS(1,6).Name reliab                   ncell x nray x nscan
%head.Vgroup(1,2).SDS(1,16).Name rainFlag                nray x nscan
%head.Vgroup(1,2).SDS(1,18).Name rainAVE                 2 x nray x nscan
%head.Vgroup(1,2).SDS(1,19).Name PrecipWaterSum          2 x nray x nscan
%head.Vgroup(1,2).SDS(1,29).Name freezH                  nray x nscan
%head.Vgroup(1,2).SDS(1,36).Name nearSurfRain            nray x nscan
%head.Vgroup(1,2).SDS(1,38).Name e_SurfRain              nray x nscan
%head.Vgroup(1,2).SDS(1,45).Name rainType                nray x nscan
    WantedIndex=[];
    WantedIndex(1)=2;
    WantedIndex(2)=3;
    WantedIndex(3)=5;
    WantedIndex(4)=6;
    WantedIndex(5)=18;
    WantedIndex(6)=19;
    WantedIndex(7)=16;
    WantedIndex(8)=29;
    WantedIndex(9)=36;
    WantedIndex(10)=38;
    WantedIndex(11)=45;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file=dir(files);
    nf=length(file); 
    FILE_NAME = strcat(foldpath,file(1).name);
    head=hdfinfo(FILE_NAME);
    varout={};
    for i=3:11
        iv=WantedIndex(i);
        varname=head.Vgroup(1,2).SDS(1,j).Name;    
        fileOut=strcat(pathout,rgns{ir},'-',yearstr,'-',varname,'_TRMM_2A25.txt');
        varout{i}=fopen(fileOut,'w');
    end
    clear head
    for n=1:nf
        FILE_NAME = strcat(foldpath,file(n).name);
        head=hdfinfo(FILE_NAME);
%%%%%%%   ScanTime
        varname=head.Vgroup(1,2).Vgroup(1,1).SDS(1,1).Name;
        year=hdfread(FILE_NAME,varname);
        varname=head.Vgroup(1,2).Vgroup(1,1).SDS(1,2).Name;
        month=hdfread(FILE_NAME,varname);
        varname=head.Vgroup(1,2).Vgroup(1,1).SDS(1,3).Name;
        daymon=hdfread(FILE_NAME,varname);
        varname=head.Vgroup(1,2).Vgroup(1,1).SDS(1,4).Name;
        hour=hdfread(FILE_NAME,varname);
        varname=head.Vgroup(1,2).Vgroup(1,1).SDS(1,5).Name;
        minute=hdfread(FILE_NAME,varname);
%%%%%%%%%% latitude  longitude
        varname=head.Vgroup(1,2).SDS(1,2).Name; %Latitude                 nray x nscan
        lat=hdfread(FILE_NAME,varname);
        varname=head.Vgroup(1,2).SDS(1,3).Name; %Longitude                nray x nscan
        lon=hdfread(FILE_NAME,varname);
        nray=length(lon(:,1));  % or nray=head.Vgroup(1,2).SDS(1,2).Dim(1,1).Size
        for i=3:4    %    ncell x nray x nscan
            outvar=varout{i};
            iv=WantedIndex(i);
            varname=head.Vgroup(1,2).SDS(1,j).Name;
            tmp=hdfread(FILE_NAME,varname);



            %% *scale_factor + add_offset  neeed!!!!!!
            for nt=1:nray
                for ns=1:nscan
                    if lon(nt,ns)<lone(ig) & lon(nt,ns)>lons(ig)  
                        if lat(nt,ns)<late(ig) & lat(nt,ns)>lats(ig)
                            for nc=1:ncell
                                fprintf(outvar,'%f  ',tmp(nc,nt,ns));
                            end  % nc 
                            fprintf(outvar,'%s\n','');
                        end %%% lat
                    end % lon
                end % ns 
            end % nt
            clear tmp;
        end %i
        for i=5:6    %    2 x nray x nscan
            outvar=varout{i};
            iv=WantedIndex(i);
            varname=head.Vgroup(1,2).SDS(1,j).Name;
            tmp=hdfread(FILE_NAME,varname);
            for nt=1:nray
                for ns=1:nscan
                    if lon(nt,ns)<lone(ig) & lon(nt,ns)>lons(ig)  
                        if lat(nt,ns)<late(ig) & lat(nt,ns)>lats(ig)
                            for nc=1:2
                                fprintf(outvar,'%f  ',tmp(nc,nt,ns));
                            end  % nc 
                            fprintf(outvar,'%s\n','');
                        end %%% lat
                    end % lon
                end % ns 
            end % nt
            clear tmp;
        end %i
        for i=7:11    %  nray x nscan
            outvar=varout{i};
            iv=WantedIndex(i);
            varname=head.Vgroup(1,2).SDS(1,j).Name;
            tmp=hdfread(FILE_NAME,varname);
            for nt=1:nray
                for ns=1:nscan
                    if lon(nt,ns)<lone(ig) & lon(nt,ns)>lons(ig)  
                        if lat(nt,ns)<late(ig) & lat(nt,ns)>lats(ig)
                            fprintf(outvar,'%f  ',tmp(nt,ns));
                            fprintf(outvar,'%s\n','');
                            if i==11
                                fprintf(outxt,'%d  ',year(nt));
                                fprintf(outxt,'%d  ',month(nt));
                                fprintf(outxt,'%d  ',daymon(nt));
                                fprintf(outxt,'%d  ',hour(nt));
                                fprintf(outxt,'%d  ',minute(nt));
                                fprintf(outxt,'%f  ',lon(nt,ns));
                                fprintf(outxt,'%f  ',lat(nt,ns));
                                fprintf(outxt,'%s\n','');
                            end %%%  i==11
                        end %%% lat
                    end % lon
                end % ns 
            end % nt
            clear tmp;
        end %i

    end %%% for nf
    sta1 = fclose(outxt);
    for i=3:11    %  nray x nscan
        outvar=varout{i};
        sta1 = fclose(outvar);
    end    
end  %%%% ig
end  %%% year        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  References
%
% [1] http://www.cloudsat.cira.colostate.edu/dataSpecs.php
