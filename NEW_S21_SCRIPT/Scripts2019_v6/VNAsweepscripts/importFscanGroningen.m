function [data,power]=importFscanGroningen
% reads in the normal dataformat from the F scans from Groningen, F Re Im. Plots all
% the curves with different color  for easy Popt determination. 
%
% INPUT
% give filename and dir in first cell
%
% OUTPUT
%
% auto saved figure file
%
% data is cell array with in each cell the full F scan
%
% power is a 1D aray with all powers

%% read in and process 
dir='C:\Documents and Settings\Jochem\My Documents\iFolder\KID\APEX camera\Experimental work\F spacing Apex\K97\K97 F7__Chip_1600__31_01_12';
file='k97 big sample power sweep -60-40 dBm  Tbase.dat';

close all;

format('long','e');
fid=fopen([dir '\' file])

textscan(fid,'%*[^\r\n]',2);%skip 2 headerlines
ni=1;
while ni<1000000
    temp=textscan(fid,'%f%f%f',1);
    if isempty(temp{1})
        break;
    end
    ni=ni+1;
end
fseek(fid, 0, 'bof');%reset pointer

textscan(fid,'%*[^\r\n]',1);%skip 1 headerlines
n=1;
while n<20 %(powers)
    temp=textscan(fid,'%*s%f',1);
    if isempty(temp{1})
        break;
    end
    power(n)=temp{1};
    data{n}=cell2mat(textscan(fid,'%f%f%f',ni));%normally F Re Im
    n=n+1;
end
npowers=n-1
fclose(fid);

colortabel1=cellstr([' -r';' -g';' -b';' -c';' -m';' -k';' :r';' :g';' :b';' :c';' :m';' :k';'--r';'--g';'--b';'--c';'--m';'--k';'-.r';'-.g';'-.b';'-.c';'-.m';'-.k';' :r';' :g';' :b';' :c';' :m';' :k';'--r';'--g';'--b';'--c';'--m';'--k']);

%% plot raw data
figure(10)
for n=1:npowers
    plot(data{n}(:,1),10*log10(data{n}(:,2).^2+data{n}(:,3).^2),colortabel1{n});hold on;
end
legend(num2str(power'),'Location','Best');
saveas(figure(10),[dir  '\Figure10' file(1:end-4) '.fig']);

end
