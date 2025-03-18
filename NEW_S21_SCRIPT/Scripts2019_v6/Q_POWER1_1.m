function Q_POWER1_1
% very old script, but usefull for fast assesment of MKID Q factors
%finds optimum power from S21 data. uses assymmetric fitting with
%FitS21_1.m
% modified for power optimum finding and some details Jochem 20-2-2007
% modified again Jochem 1-3-08 to get Q  better Popt finding, it only uses opt Qi method if Qi ai 1.05 times Qi by Fchange method ; 
% also BW range increasd to BW/7 in stead of BW/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear variables
dpath=[cd '/../S21/Power'];%Directory of data, no FSP @ end
refpowerindex=3;%1=lowest P, 2 is single lowest P etc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qrt matrix info: 
% 1 T                                       
% 2 Q (from fitS21)
% 3 Qc
% 4 Qi
% 5 responsivity old method [rad/qp]
% 6 response [degr]
% 7 Re point response (kidcall. circle)
% 8 Im point response
% 9 Fres [Ghz] (from fitS21)
% 10 S21,min (from fitS21)
% 11 phase @ resonance
% 12 qp number
% 13 power
% 14 fres/fres(T=lowest)-1
% 15 responsivity {stephe's method}
% 16 dxfit
% 17 normalised freq. responsivity (f-f0)/f0
% 18 internal power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files=dir([dpath '/KID*.dat']);
file=files(1).name;
info=textscan(file,'%*3s %*f %*1s %*f %s');%bit of filename after power number=info{1}
ChipID=char(info{1}(1)); 
ChipID=ChipID(5:cell2mat(strfind(info{1},'.dat'))-1);%automatically finds Chip ID in Filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading in all the filenames in array of cells
m=1;
if length(files)~=0
    for n=1:length(files)%for loop used to display KID numbers and not all KID files
        temp=cell2mat(textscan(files(n).name,'%*3s %f %*s'));
        if m==1
            kids(m)=temp; 
            m=m+1;
        elseif kids(m-1)~=temp 
            kids(m)=temp; 
            m=m+1;
        end
    end
end
numberofscans=length(files);
disp(['Kids available: ' num2str(kids) ', total number of files ' num2str(numberofscans)]);
for n=1:length(kids)%finds powers used for each KID
    filetjes=dir([dpath filesep 'KID' num2str(kids(n)) '_*.dat']);
    for m=1:length(filetjes)
        KIDpowers{n}(m)=cell2mat(textscan(filetjes(m).name,'%*3c%*f%*1c%f%*s'));
    end
%disp(['KID ' num2str(kids(n)) ' has powers']);
%disp(KIDpowers{n}(:));
end

%%%%%%%%%Reading in the data KID by KID%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:length(kids)%everyKID
    disp(['KID ' num2str(kids(n))]);
    for m=1:length(KIDpowers{n}(:))%Every power. Results of all powers will be combined
        disp(['@' num2str(KIDpowers{n}(m))]);
        file=[dpath filesep 'KID' num2str(kids(n)) '_' num2str(KIDpowers{n}(m)) 'dBm_' ChipID '.dat' ];
        fid=fopen(file);
        %disp(file);
        temp=textscan(fid,'%*s%*s%*c%*c%*c%*c%f%*s',1,'headerlines',1);
        Power=-1*temp{1};
        temp=textscan(fid,'%*2c%f%*s%*s',1,'headerlines',2);
        Q=temp{1};%Q factor
        i=1;bla=[1,1,1];
        temp=textscan(fid,'%*s%*s%*2c%f%*s',1,'headerlines',1);
        T=temp{1};
        temp=textscan(fid,'%*s%*s%*s',1,'headerlines',1);
        while ~isempty(bla)%bla is not empty
            bla=cell2mat(textscan(fid, '%f%f%f', 1));
            if ~isempty(bla)
                data{m}(i,:)=bla;
            end
            i=i+1;
        end
        fclose(fid);
        %%%%Data (GHz dB)is read in for 1 KID now we obtain parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data{m}(:,1) = data{m}(:,1)*1e9;%to Hz
        filtered=smooth(data{m}(:,2),3); %smoothing
        maximum=max(filtered); %find max value
        data{m}(:,2)=data{m}(:,2)-maximum;  %norm S21
        [S21resdB,index] = min(filtered); %find minimum
        S21resdB=S21resdB-maximum ; %normalise S21min
        S21res = 10^(S21resdB/20);
        fres=data{m}(index,1); %get resonance frequency NOW in GHZ!!!   
        
        %fit
        qrt(m,10)=S21res;%used as seeding value
        qrt(m,9)=fres;%used as seeding value
        qrt(m,23)=S21res;%stored
        qrt(m,22)=fres;%stored
        Qi=qrt(m,2)./S21res;%Get Qi=Q/S21min
        Qc=(Qi*Q)/(Qi-Q);%Get Qc=(Qi-Q/Qi*Q)
        addpath('VNAsweepscripts') ;
        [result,fitted] = FitS21simple(data{m}(:,1:2),fres,Q,S21res,0.5);
        fitresult{m}(:,1) = data{m}(:,1);
        fitresult{m}(:,2) = fitted;
%         figure(123)

%         plot(data{m}(:,1),data{m}(:,2),'.');
%         hold on;plot(fitresult{m}(:,1),fitresult{m}(:,2),'-b')
%         
        %Keep smoothed values for S21min and fres, since assymetric fit shifts f. 
        qrt(m,2)=result.Q; %fitted Q factor
        qrt(m,9)=result.Fr; % fitted fres
        qrt(m,10)=result.Smin; % fitted S21min
        qrt(m,4)=qrt(m,2)./qrt(m,10);%Get Qi=Q/S21min
        qrt(m,3)=(qrt(m,4).*qrt(m,2))/(qrt(m,4)-qrt(m,2));%Get Qc=(Qi-Q/QiQ)
        qrt(m,1)=T;
        qrt(m,13)=Power;
        qrt(m,18)=10*log10((2/pi)*(10.^(qrt(m,13)./10).*(qrt(m,2).^2./qrt(m,3))));%internal power
        
    end
    [qrt,sortindex]=sortrows(qrt,13);%now power is increasing with index!!!!!
    %%%%all powers and parameters of 1 KID done
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now we find optimum power%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,z]=sort(qrt(:,13));%sorting vs power!
    y=z(refpowerindex);%same as the one commented out below but allowing for offset use
    %[x,y]=min(qrt(:,13));%minimum power @ KID
    BW=qrt(y,9)/qrt(y,2);%Q=omega/BW so BW = omega/Q
    minF=qrt(y,9)-BW/7;%this is the minimum freq. alllowed and will be used to find power opt.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j=refpowerindex;%starting at refernce power
    warning=1;
    if qrt(1,13) < qrt(2,13) %power is increasing, is allways the case, so no else here (data sorted to power)! 
        exit=0;
        while exit==0
            if  qrt(j,9) < minF;%if res. freq is decreased below fmin
                warning=0;
                exit=1;
            end
            if j==length(qrt(:,13));
                 if warning==1
                    disp(['possibly no good value, max power is used for KID' num2str(kids(n))]);
                    [x,j]=max(qrt(:,13));j=j+1;
                    exit=1;
                else
                    exit=1;
                end
            end
                j=j+1;
        end
        if j>2
            j=j-2;%finds 1 (not2) points lower in power as best; the extra factor due to j=j+1 that is added after exit!!
        elseif j>1
            j=j-1;
        else
            j=1;
            disp(['possibly no good value, min power is used for KID' num2str(kids(n))]);
        end
   end
    
    %15%Higher Qi than the one at minimum power is opt. power
    %Get in jj 5 points around j 

    if  j>=4 && qrt(1,13) < qrt(2,13) && max(qrt(j-3:j,4))>1.05*qrt(j,4);
        %MODIFIED 2-3-08 do only if there is room to play around optimum power as found before AND power increasing with filenumber
        %AND if the optimum Qi is at least 5% higher than the one found previously
        [bla,jj]=max(qrt(j-3:j,4)); %finds maximum Qi at or below in power to prev. measurement
        jj=jj+j-4;
%    elseif  j>=4 && qrt(1,13) > qrt(2,13)%do only if there is room to play around optimum power as found before AND power decreasing with filenumber
 %       [bla,jj]=max(qrt(j:j+4,4)); %finds maximum Qi at or below in power to prev. measurement
 %       jj=jj+j+5;
 %       disp('check Q_power3_2 for correct optimum finding line 102')
    else 
        jj=j;
    end     
    optpower(n,1)=kids(n);   
    optpower(n,2)=qrt(jj,9); %resonance frequency at optimum power, from Fit
    optpower(n,3)=qrt(jj,13); %Popt Best Qi method 
    optpower(n,4)=qrt(jj,4);%Qi
    optpower(n,5)=qrt(jj,3);%Qi
    output(:,2:5)=qrt(:,1:4);
    output(:,1)=kids(n);
    output(:,6)=qrt(:,9)*1E9;%to Hz!
    output(:,8)=qrt(:,13);
    output(:,7)=20*log10(qrt(m,10));%S21min
    output(:,9)=qrt(:,18);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting the results%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
    figure(kids(n));hold off; clf;
    %Set plotranges
    BWfactor=2;
    range=BWfactor*max(qrt(:,9)./qrt(:,2));
    xmin=max(qrt(:,9))-range/2;
    xmax=max(qrt(:,9))+range/2;
    ymax=1;ymin=0;
    
    %PLOT DATA
    subplot(2,2,1);
    colortabel=['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k'];
    for plotbla=1:length(data)
        plot(data{sortindex(plotbla)}(:,1),20*log10(data{sortindex(plotbla)}(:,2)),colortabel(plotbla));hold on;%volt S21
        xlabel('F [GHz]');ylabel('S21 ');title([ChipID ' KID' num2str(kids(n))]);axis([xmin xmax ymin ymax]);axis 'auto y';%S21 vs freq
        legendarray(plotbla)=qrt(plotbla,13);%qrt is sorted!, so sortindex not needed here!
    end
    legend(num2str(legendarray'),-1);
    %plot_lines(minF);plot_lines(qrt(y,9));
    
    subplot(2,2,2)
    semilogy(qrt(:,13),qrt(:,2),'-ok','MarkerSize',2);
    hold on;semilogy(qrt(:,13),qrt(:,3),'-or','MarkerSize',2);
    hold on;semilogy(qrt(:,13),qrt(:,4),'-og','MarkerSize',2);legend('Q','Qc','Qi',2);
    hold on; plot(qrt(j,13),[qrt(j,2),qrt(j,3),qrt(j,4)],'-ob','MarkerSize',5,'MarkerFaceColor','b');
    hold on; plot(qrt(jj,13),[qrt(jj,2),qrt(jj,3),qrt(jj,4)],'-or','MarkerSize',7);
    hold on; plot(qrt(y,13),[qrt(y,2),qrt(y,3),qrt(y,4)],'-og','MarkerSize',8);
    xlabel('P at KID [dBm]');ylabel('Q ');title([ ChipID ' KID' num2str(kids(n)) 'Q factors']);%S
        
    subplot(2,2,3)
    %         plot(data{m}(:,1),data{m}(:,2),'.');
%         hold on;plot(fitresult{m}(:,1),fitresult{m}(:,2),'-b')
    
    
    plot(data{sortindex(jj)}(:,1),(data{sortindex(jj)}(:,2)),'o');hold on;
    legend(['opt power=' num2str(qrt(jj,13))]);
    plot(fitresult{sortindex(jj)}(:,1),fitresult{sortindex(jj)}(:,2),'r')
    xlabel('F [GHz]');ylabel('S21 ');title([ChipID ' KID' num2str(kids(n))]);axis([xmin xmax ymin ymax]);axis 'auto y';
    hold off;
    
    subplot(2,2,4)
    plot(qrt(j,13),qrt(j,9),'-ob','MarkerSize',5,'MarkerFaceColor','b');%Freq method
    hold on; plot(qrt(jj,13),qrt(jj,9),'-or','MarkerSize',7);%opt Qi
    hold on; plot(qrt(y,13),qrt(y,9),'-og','MarkerSize',8);%ref
    xlabel('P at KID [dBm]');ylabel('Fres ');title([ ChipID ' KID' num2str(kids(n)) 'Res. Freq.']);%S
    legend('Freq. meth.','Opt Qi','ref.','Location','West');
    plot(qrt(:,13),qrt(:,9),'-ok','MarkerSize',3);hold on;
    plot(qrt(:,13),qrt(:,22),'xk','MarkerSize',3);hold on;
    minff=ones(length(qrt(:,13)),1);minff(:)=minF;
    hold on;plot(qrt(:,13),minff);
    
    %Write result data
    saveas(gcf,[dpath filesep 'powerdep_KID' '_' num2str(kids(n)) '.fig'],'fig')

    colortabel1=cellstr(['or';'og';'ob';'oc';'om';'ok';'dr';'dg';'db';'dc';'dm';'dk';'xr';'xg';'xb';'xc';'xm';'xk';'or';'og';'ob';'oc';'om';'ok';'dr';'dg';'db';'dc';'dm';'dk';'xr';'xg';'xb';'xc';'xm';'xk']);
    figure(10000)
    subplot(1,2,1)
    semilogy(output(:,9),output(:,3),colortabel1{n});hold on;
    xlabel('Pi [dBm]');ylabel('Q ');
    subplot(1,2,2)
    semilogy(output(:,9),output(:,5),colortabel1{n});hold on;
    xlabel('Pi [dBm]');ylabel('Qi ');
    legend(num2str(kids'));
    saveas(gcf,[dpath filesep 'powerdep_allKIDs' '.fig'],'fig')
    
    clear data;
end;


format('short','g')
display('results from S21 finding optimum Qi using 5 pts around max F-shift: ID Fres Qi Qc Popt')
display(optpower(:,:));
end

% 
% function  [Fresres,Qfres,S21minres,l,stretch,result] = FitS21_1(data,Fres,Qf,S21min)
% %perfoms skewd Lorenztain fit for Q, parabolic fit for Fres and S21min
% %Data should contain F[GHz], S21 magn., (phase optional). Fres is res freq
% %in GHz, Qf is initial Q factor and S21min is original S21min. These values
% %are used as guess values for a fit to Bens eqn. 
% %07-03-07: modified for bad ini guesses bandwidth
% format('long','e');
% 
% no=2; %width in resonator bandwidths to fit over
% datasize=size(data(:,2),1);
% %limit width of fit
% bandwidth=Fres/Qf; %resonator bandwidth
% minfitindex=find(data(:,1)>(Fres-no*bandwidth/2),1);
% maxfitindex=find(data(:,1)>(Fres+no*bandwidth/2),1);
% if minfitindex+4<maxfitindex
%     x=data(minfitindex:maxfitindex,1);y=data(minfitindex:maxfitindex,2).^2;
% else
%     x=data(:,1);y=data(:,2).^2;
% end
% % estimate linear background
% l=(data(1,2)^2-data(datasize,2)^2)/(data(1,1)-data(datasize,1)); %%linear background drift
% 
% %skewterm
% 
% % catch error in finding guess parameter
% low3dbpoint=find(data(:,1)>(Fres-bandwidth/2),1);
% if (isempty(low3dbpoint)==1)
%     low3dbpoint=1;
% end
% high3dbpoint=find(data(:,1)>(Fres+bandwidth/2),1);
% if (isempty(high3dbpoint)==1)
%     high3dbpoint=size(data(:,1),1);
% end
% 
% stretch=(data(low3dbpoint,2)^2-data(high3dbpoint,2)^2)/(bandwidth); % estimate skew parameter
% 
% s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min l stretch],'Lower',[Fres*0.9 Qf/2 S21min/4 -25*abs(l) -25*abs(stretch)],...
%     'Upper',[Fres*1.1 Qf*2 S21min*2 abs(l)*25 +25*abs(stretch)],'MaxFunEvals',1000);
% 
% %Fit Lorentzian like in power space with linear background and stretch
% ftype = fittype('(1-(1-Smin^2)*(1+stretch*(x-Fr)/Fr)/(1+(2*Q*(x-Fr)/Fr)^2)+l*(x-Fr)/Fr)','options', s);  
% %[result,gof,output]=fit(x,y,ftype);
% result=fit(x,y,ftype);
% %figure(1);plot(result);hold on;plot(x,y,'.')
% Qfres=result.Q;
% l=result.l;
% stretch=result.stretch;
% %%new: a simple lorentzian in log space fit to find resonance
% %%frequency/S21min. Is done in power space
% so=0.3;
% 
% minfitindex=find(data(:,1)>(Fres-so*bandwidth/2),1);
% maxfitindex=find(data(:,1)>(Fres+so*bandwidth/2),1);
% if maxfitindex-minfitindex>10
%     x=(data(minfitindex:maxfitindex,1));
%     y=20*log10(data(minfitindex:maxfitindex,2));
% else
%     x=(data(:,1));
%     y=20*log10(data(:,2));
%     disp('bandwidthguess for fit far off, whole S21 curve used for Fo and S21');
% end
% 
% s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint',[Fres Qf S21min], 'MaxFunEvals',100);
% 
% ftype = fittype('10*log10(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2))','options', s);
% [resultFres]=fit(x,y,ftype);
% Fresres=resultFres.Fr;
% S21minres=resultFres.Smin;
