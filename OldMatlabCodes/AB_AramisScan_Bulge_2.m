function AramisScan_Bulge_2();

clear all;
fclose all;

frompath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Recalc\AramisExport';     %Folder where the files you're reading from are
fromprefix='BT2-Recalc-Stage-0-';     %Prefix of the inidividual file names

topath='F:\Martin_Experiments\Sandia_Anis\Bulge Tests\BT-2_Recalc\AramisExport_MissingRemoved';    %Where you want to create the folder of the processed files
toprefix='BT2-Recalc-Stage-0-';

web(topath,'-browser')

%Number of final stage
last=218;

k=1;
for i=0:1:last
%for i=0:1:last;
    tic
    opentext=sprintf('%s\\%s%d.txt',frompath,fromprefix,i);
    fidfrom=fopen(opentext);  
    data=textscan(fidfrom,'%f %f %f %f %f %f %f %f %f %f %f %f %f',...
        'Delimiter',',','EndOfLine','\n','CommentStyle','#','MultipleDelimsAsOne',1,'Headerlines',2);
    L=length(data{13});
    for z=1:13;
        if length(data{z})~=L;
            data{z}(L)=[];
        end;
    end;
    data=cell2mat(data);    %Convert the cell to matrix
    fclose(fidfrom);        %Close the file we read from
    data(any(isnan(data),2),:)=[]; 
    savestring=sprintf('%s\\%s%d.dat',topath,toprefix,i);    %File name with stage number
    fidto=fopen(savestring,'a+');   %a+ simply means we append data to that already in the file, which in this case is none.
    fprintf(fidto,'%d %d %.8f %.8f %.8f %.8f %.8f %.10f %.10f %.10f %.8f %.8f %.8f\n',data');
    fclose(fidto);
    clear data savestring
    k=k+1;
    toc
end;