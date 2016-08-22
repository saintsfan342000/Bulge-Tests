function SyncLVDIC();

% First export a file from aramis that has a column of stage and a
% column of time.  Interpolate pressure and lvdt from the labview file based off time stamp 

%%%% Absolute path of the DIC Stage-Time file and the Labview File
    %(1) Stage (2)Time
    DICfile='F:\BulgeTest_SD-1\StgTime.txt';
    %(1)Time (2)LVDT[in] (3)Pressure(psi)
    LVfile='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-1\LV_BT-1.dat';

%%%% Savepath
    savepath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-1';
    
%Open up DIC File
fid=fopen(DICfile,'r');
ST=cell2mat(textscan(fid,'%f %f','Headerlines',1,'Delimiter',' '));
fclose(fid);
clear fid;

%Open up the labview dat
fid=fopen(LVfile,'r');
TLP=cell2mat(textscan(fid,'%f %f %f','Headerlines',1,'Delimiter',' '));
fclose(fid);
clear fid;

%INterpolate lvdt and pressure at DIC stage time points
lvdt_pres_interped=interp1(TLP(:,1),TLP(:,[2 3]),ST(:,2));

%Export a file 
%(1) Stage (2) Time (3) LVDT (4) Pressure
fid=fopen(sprintf('%s\\STLP.dat',savepath),'w');
fprintf(fid,'(1) Stage (2) Time (3) LVDT (4) Pressure\n');
fprintf(fid,'%.0f %.4f %.8f %.8f\n',[ST lvdt_pres_interped]');
fclose all;