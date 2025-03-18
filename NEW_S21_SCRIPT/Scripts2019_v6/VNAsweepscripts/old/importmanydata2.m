function data=importmanydata2
%% Used to read in the data in this dir.
directory='C:\Documents and Settings\Jochem\Desktop\NbTa Groningen\Ta APEX_minimal wirebonds\';
files=dir([directory '*.d1']);
start=5470;stop=5772;
step=50/1601;%frequencie ramge all files, please get manually from header
DELIMITER=',';
HEADERLINES=9;
for i=1:length(files)
    bla{i}=files(i).name
end

%%

for i=1:length(files)%NEED TO ADD FREE CHOICE OF FILES -1 to nod reed in last one
    fid=fopen([directory files(i).name]);[directory files(i).name]
    textscan(fid,'%*[^\n]',6);%headerlines does not work....
    bla=textscan(fid,'%*s %f %f %f',1);
    fclose(fid);
    start=bla{1};stop=bla{2};
    npts=bla{3};
    data1=importdata([directory files(i).name], DELIMITER, HEADERLINES);
    [numrows,numcols]=size(data1.data);
    temp=data1.data(1:numrows,1:numcols);
    data((i-1)*(numrows)+1:i*(numrows),2:1+numcols)=temp;
    data((i-1)*(numrows)+1:i*(numrows),1)=linspace(start,stop,numrows);
    [npts numrows]
end

plot(data(:,1),10*log10(data(:,2).^2+data(:,3).^2),'.-')
end
