function LIFETIME7
    % Uses test2 gui as a basis to select lifetime files, and then fit
    % them with a single exponential over a specified time range. 
    %22/2/10 v7: fit dRdtheta, export in drdtheta compatible, put amp life
    %into life2 array, some cosmetic changes 
    %modified butoon to select KID and filename extension
    %amplitude fitting
    %modified for new software from 6_3 (only pulse power reading line 281)
    %seeding of fits added to front panel
    %export readout power
    colourtable=['r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k','r','g','b','c','m','y','k'];
    %Initialise usefull variables
    %pathy='J:\HOME\Kid\NIKA development\APEXNIKA5_NIKAdark__29_6_09 12_19';
    pathy=[cd filesep '..'];
%    pathy='J:\HOME\Kid\experiments\BoxinBoxsetup\C2_newsetup__30_11_09 12_36\LifetimehighQ'
    %pathy=pwd;
    KIDnum=0;%0 for all KIDs
    
    cd(pathy);
    % declare variables for latter use

    klist{1}=[]; kidlist{1}=[]; powerlist{1}=[]; VNAlist{1}=[]; Tlist{1}=[];Qlist{1}=[];f0list{1}=[];mintrange=0;maxtrange=3500;
    kmax=2; lengtht=0;lifefit{1}=[]; deltat=0; lifelist{1}=[];lowtlist{1}=[];life2list{1}=[];hightlist{1}=[];   
    KIDdata{1}=[];pulsedata{1}=[]; lowradlist{1}=[];highradlist{1}=[];minradrange=0;maxradrange=0.5;exportselect=0;resultfile=0;
    alist1{1}=[];alist2{1}=[];drdthlist{1}=[];p{1}=[];thdata=[];
    hplot=[];
    radtimepopup={'rad','time','both'};
    singledoublepopup={'single exp','double exp'};
    exportpopup={'Mean','All','dRdTheta'};
    fitpopup={'phase','amplitude','dRdtheta'};
% seed fit:
    seedlife1=1000
    seedlife2=1000
    
    
    %  Initialize and hide the GUI as it is being constructed.
    hf = figure('OuterPosition',[360,600,600,600]);% % different window for figure- ease of saving
   % haxes = axes('parent',hf,'units','pixels','position',[200 60 300 300]) ;       
    haxes = axes('parent',hf) ;       
    % Construct the components.
    hfigmain = figure('Visible','off','units','pixels','Position',[150,400,700,500]); %%
    htextpathy  = uicontrol(hfigmain,'Style','text','String',pathy,...
           'Position',[10,400,350,15]);
    htext  = uicontrol(hfigmain,'Style','text','String','v7, can save out in dRtheta.csv compatible. To create file - calc phase, then amp (stored in life2), then drdtheta',...
           'Position',[10,440,400,30]);
       
    hradtimepopup =uicontrol(hfigmain,'Style','popupmenu','String',radtimepopup,...
            'Position',[300 200 100 20])
    htextradtime = uicontrol(hfigmain,'Style','text','String','Fit in angle or time',...
           'Position',[300,230,100,20]);
    hsingledoublepopup =uicontrol(hfigmain,'Style','popupmenu','String',singledoublepopup,...
            'Position',[300 260 100 20])
    htextsingledoublepopup = uicontrol(hfigmain,'Style','text','String','Fit dbl or sgl',...
           'Position',[300,290,100,20]);
    hexportpopup =uicontrol(hfigmain,'Style','popupmenu','String',exportpopup,...
            'Position',[300 320 100 20])
    htextexportpopup = uicontrol(hfigmain,'Style','text','String','Export',...
           'Position',[300,350,100,20]);
    hfitpopup =uicontrol(hfigmain,'Style','popupmenu','String',fitpopup,...
            'Position',[400 320 100 20])
    htextfitpopup = uicontrol(hfigmain,'Style','text','String','Fit',...
           'Position',[400,350,100,20]);   
       
    % KID list box  
    hlist = uicontrol(hfigmain,'Style','listbox',...
                'String',klist,...
                'Value',1,...
                'Max',1,'Min',0,...
                'Callback',{@update_callback},...
                'Position',[10 10 140 120]);
    htextkidlist = uicontrol(hfigmain,'Style','text','String','KID',...
           'Position',[150,130,30,10]);
    hkidlist = uicontrol(hfigmain,'Style','listbox',...
                'String',kidlist,...
                'Value',[],...
                'Max',2,'Min',0,...
                'Position',[150 10 40 120]);
    htextpowerlist = uicontrol(hfigmain,'Style','text','String','power',...
           'Position',[190,130,30,10]);
    hpowerlist = uicontrol(hfigmain,'Style','listbox',...
                'String',powerlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[190 10 40 120]);
    htextVNAlist = uicontrol(hfigmain,'Style','text','String','VNA',...
               'Position',[230,130,30,10]);
    hVNAlist = uicontrol(hfigmain,'Style','listbox',...
                'String',VNAlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[230 10 40 120]);
    htextTlist = uicontrol(hfigmain,'Style','text','String','T',...
                'Position',[270,130,30,10]);
    hTlist = uicontrol(hfigmain,'Style','listbox',...
                'String',Tlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[270 10 40 120]);
    htextlifelist = uicontrol(hfigmain,'Style','text','String','Life',...
                'Position',[310,130,40,10]);
    hlifelist = uicontrol(hfigmain,'Style','listbox',...
                'String',lifelist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[310 10 50 120]);
    htextlife2list = uicontrol(hfigmain,'Style','text','String','Life2',...
                'Position',[360,130,30,10]);
    hlife2list = uicontrol(hfigmain,'Style','listbox',...
                'String',life2list,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[360 10 50 120]);
    htextdRdtheta = uicontrol(hfigmain,'Style','text','String','dRdth',...
                'Position',[410,130,40,10]);
    hdrdthlist = uicontrol(hfigmain,'Style','listbox',...
                'String',drdthlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[410 10 50 120]);
    htextlowtlist = uicontrol(hfigmain,'Style','text','String','Low t',...
                'Position',[460,130,30,10]); 
    hlowtlist = uicontrol(hfigmain,'Style','listbox',...
                'String',lowtlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[460 10 40 120]);
    htexthightlist = uicontrol(hfigmain,'Style','text','String','High t',...
                'Position',[500,130,30,10]);
    hhightlist = uicontrol(hfigmain,'Style','listbox',...
                'String',hightlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[500 10 40 120]);
    htextlowradlist = uicontrol(hfigmain,'Style','text','String','Low rad',...
                'Position',[540,130,30,10]);
    hlowradlist = uicontrol(hfigmain,'Style','listbox',...
                'String',lowradlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[540 10 40 120]);
    htexthighradlist = uicontrol(hfigmain,'Style','text','String','High rad',...
                'Position',[580,130,30,10]);
    hhighradlist = uicontrol(hfigmain,'Style','listbox',...
                'String',highradlist,...
                'Value',[],...
                'Max',kmax,'Min',0,...
                'Position',[580 10 40 120]);
     
     
            
    % add pushbutton       
    hupdate = uicontrol(hfigmain,'Style','pushbutton',...
                'String','update plot',...
                'Callback',{@update_callback},...
                'Position',[25 150 60 20]);

    haddfile = uicontrol(hfigmain,'Style','pushbutton','String','add file',...
                'Callback',{@add_file_callback},...
                'Position',[90 150 60 20]);
    hremfile = uicontrol(hfigmain,'Style','pushbutton','String','remove file',...
                'Callback',{@rem_file_callback},...
                'Position',[160 150 60 20]);       
    hremfileall = uicontrol(hfigmain,'Style','pushbutton','String','remove all files',...
                'Callback',{@rem_file_callback_all},...
                'Position',[230 150 60 20]);       
    hfit = uicontrol(hfigmain,'Style','pushbutton','String','fit',...
                'Callback',{@fit_callback},...
                'Position',[160 200 60 20]);
    hfitall = uicontrol(hfigmain,'Style','pushbutton','String','fit all',...
                'Callback',{@fit_callback_all},...
                'Position',[230 200 60 20]);
            
   
            
    %Editable text
    hmintrangetext = uicontrol(hfigmain,'Style','text','String','tmin',...
                'Position',[25 230 60 20]);
    hmintrange = uicontrol(hfigmain,'Style','edit','String',num2str(mintrange),...
                'Position',[25 200 60 20]);
    hmaxtrangetext = uicontrol(hfigmain,'Style','text','String','t max',...
                'Position',[90 230 60 20]);
    hmaxtrange = uicontrol(hfigmain,'Style','edit','String',num2str(maxtrange),...
                'Position',[90 200 60 20]);
 
    hminradrangetext = uicontrol(hfigmain,'Style','text','String','angle min',...
                'Position',[25 290 60 20]);
    hminradrange = uicontrol(hfigmain,'Style','edit','String',num2str(minradrange),...
                'Position',[25 260 60 20]);
    hmaxradrangetext = uicontrol(hfigmain,'Style','text','String','angle max',...
                'Position',[90 290 60 20]);
    hmaxradrange = uicontrol(hfigmain,'Style','edit','String',num2str(maxradrange),...
                'Position',[90 260 60 20]);
    hseedlife1 = uicontrol(hfigmain,'Style','edit','String',num2str(seedlife1),...
                'Position',[25 320 60 20]);
    hseedlife1text = uicontrol(hfigmain,'Style','text','String','seed lifetime',...
                'Position',[25 350 60 20]);
    hseedlife2 = uicontrol(hfigmain,'Style','edit','String',num2str(seedlife2),...
                'Position',[90 320 60 20]);
    hseedlife2text = uicontrol(hfigmain,'Style','text','String','seed timeconstant',...
                'Position',[90 360 60 30]);     
            
    
            
    hKIDnumtext = uicontrol(hfigmain,'Style','text','String','KID #',...
                'Position',[155 290 60 20]);
    hKIDnum = uicontrol(hfigmain,'Style','edit','String',num2str(KIDnum),...
                'Position',[155 260 60 20]);
            
    %Menu bar
    hmfile = uimenu(hfigmain,'Label','file');
    hexportinitialise = uimenu(hmfile,'Label','init file',...
                        'Callback',{@exportinitialise_callback});
    hexportparam = uimenu(hmfile,'Label','export param',...
                        'Callback',{@exportparam_callback});

    set(hfigmain,'MenuBar','none');    % Hide standard menu bar menus.
    set(hf,'Toolbar','figure');  % Display the standard toolbar
    % Assign the GUI a name to appear in the window title.
    set(hf,'Name','Lifetime fit');
    set(hfigmain,'Name','Lifetime');
    % Move the GUI to the center of the screen.
    movegui(hf,'center');
    movegui(hfigmain,'center');
    %gui visible
    set(hf,'Visible','on');
    set(hfigmain,'Visible','on');
    
    function rem_file_callback_all(source,eventdata)
        for i=1:size(klist,2);
            set(hlist,'Value',1);
            rem_file_callback;
        end
    end

    function rem_file_callback(source,eventdata)
        %%% function removes data from memory, clears lists
        remfileindex=get(hlist,'Value') ;       %%% get selected plot indexs
        for loopy=1:size(remfileindex,2);
            fucker{loopy}=num2str(kidlist{remfileindex(loopy)})
        end
        %button = questdlg(['Remove ?' fucker],'Remove','Yes','No','Yes')
        button='Yes';
        if button=='Yes' %%%clear data from memory
            VNAlist=remfun(VNAlist,remfileindex);
            kidlist=remfun(kidlist,remfileindex);
            klist=remfun(klist,remfileindex);
            powerlist=remfun(powerlist,remfileindex);
            Tlist=remfun(Tlist,remfileindex);
            pulsedata=remfun(pulsedata,remfileindex);
            KIDdata=remfun(KIDdata,remfileindex);
            lifelist=remfun(lifelist,remfileindex);
            life2list=remfun(lifelist,remfileindex);
            drdthlist=remfun(drdthlist,remfileindex);
            hightlist=remfun(hightlist,remfileindex);
            lowtlist=remfun(lowtlist,remfileindex);
            highradlist=remfun(highradlist,remfileindex);
            lowradlist=remfun(lowradlist,remfileindex);
            highradlist=remfun(alist1,remfileindex);
            lowradlist=remfun(alist2,remfileindex);
            % set lists, and position in lists to nul
            set(hkidlist,'Value',[]);
            set(hkidlist,'String',[kidlist]);
            set(hTlist,'Value',[]);
            set(hTlist,'String',[Tlist]);
            set(hlist,'Value',1);
            set(hlist,'String',[klist]);
            set(hpowerlist,'Value',[]);
            set(hpowerlist,'String',[powerlist]);
            set(hVNAlist,'Value',[]);
            set(hVNAlist,'String',VNAlist);
            set(hlifelist,'Value',[]);
            set(hlifelist,'String',lifelist);
            set(hlife2list,'Value',[]);
            set(hlife2list,'String',life2list);
            set(hlowtlist,'Value',[]);
            set(hlowtlist,'String',lowtlist);
            set(hhightlist,'Value',[]);
            set(hhightlist,'String',hightlist);
            set(hlowradlist,'Value',[]);
            set(hlowradlist,'String',hightlist);
            set(hhighradlist,'Value',[]);
            set(hhighradlist,'String',hightlist);
            set(hdrdthlist,'Value',[]);
            set(hdrdthlist,'String',drdthlist);           
        end
        
        
    end
    function add_file_callback(source,eventdata)
    KIDnum=str2num(get(hKIDnum,'String'))  ;
    %new load
            bob=size(klist,2);
            if isempty(klist{bob})==1;
                bobby=0;
            else
                bobby=1;
            end
            if KIDnum==0
                [klist{bob+bobby},pathy]=uigetfile('KID*mK.dat','Select the file');
            else
                [klist{bob+bobby},pathy]=uigetfile(['KID' num2str(KIDnum) '*mK.dat'],'Select the file');
            end
            cd(pathy);
            fid=fopen(klist{bob+bobby},'r');
            throw2=cell2mat(textscan(fid, '%*s %*s %*4s%f ',1, 'headerLines' ,1));
            throw=cell2mat(textscan(fid, '%*s %*s %*2s%f ',1, 'headerLines' ,5)); %% get T

            C=cell2mat(textscan(fid,'',1, 'headerLines' ,9));
            nopowers=size(C,2)/4; % calculate no. of columns
            fseek(fid,0,'bof'); %reset pointer
            %format read strings
            readstringy=['%*s ' '%*s ' '%*s' '%*s' '%*s']; % remove "LED time: 100, LED currents:,"
            readstringy2=[];
            for i=1:nopowers
                readstringy=[readstringy '%f '];  %read powers 
                readstringy=[readstringy '%*s ']; % ignore , between powers
                readstringy2=[readstringy2 '%f %f %f %f'];  % read string for file, 4 columns per power
            end
            % get data
            
            pulsepowers=cell2mat(textscan(fid,readstringy,1, 'headerLines' ,8));
            loadpulsedata=cell2mat(textscan(fid,readstringy2,'headerlines',1)) ;
            fclose(fid);
            % sort data and lists
            lengtht=size(loadpulsedata,1); %% size of time data 

            for i=0:(nopowers-1)
                powerlist{bob+bobby+i}=throw2;
                Tlist{bob+bobby+i}=throw; %% T
                pulsedata{bob+bobby+i}=loadpulsedata(:,(4*i+1):4*(i+1)) ;   %% sort pulse data
    
                VNAlist{bob+bobby+i}=pulsepowers(i+1);  % update list
                pulsedata{bob+bobby+i}(:,1)=(pulsedata{bob+bobby+i}(:,1)-pulsedata{bob+bobby+i}(int32(lengtht/10),1))*1e6; %change to relative time from pulse
                pulsedata{bob+bobby+i}(:,4)=(pulsedata{bob+bobby+i}(:,4)-pulsedata{bob+bobby+i}(lengtht,4)); %subtract offset
                pulsedata{bob+bobby+i}(:,5)=sqrt(pulsedata{bob+bobby+i}(:,2).^2+pulsedata{bob+bobby+i}(:,3).^2); %amplitude
 
%                pulsedata{bob+bobby+i}(:,5)=pulsedata{bob+bobby+i}(size(pulsedata{bob+bobby+i}(:,5),1),5)-pulsedata{bob+bobby+i}(:,5); %normalise around zero
                pulsedata{bob+bobby+i}(:,5)=pulsedata{bob+bobby+i}(:,5)/mean(pulsedata{bob+bobby+i}((lengtht-lengtht/10):lengtht,5)); %
                klist{bob+bobby+i}=klist{bob+bobby};
            end
            % update window
            set(hlist,'String',klist);
            set(hpowerlist,'String',powerlist);
            kidlist=kidnumber(klist);
            set(hkidlist,'String',[kidlist]);
            set(hVNAlist,'String',VNAlist);
            set(hTlist,'String',Tlist);
                
       
    end
    

    function update_callback(source,eventdata) 
    %update plots

        figure(hf);
        plotlistindex=get(hlist,'Value')  ;      %% get selected plot indexs
        plotlist{:}=klist{plotlistindex};
        fitselect=get(hfitpopup,'Value');  % range in time or radians
        set(hf,'NextPlot','add');
        axes(haxes);
        hold off
        cla
        
            if plotlistindex(1)>size(hightlist,2)
                maxtrange=max(pulsedata{plotlistindex(1)}(:,1));
                mintrange=min(pulsedata{plotlistindex(1)}(:,1));
                pooeyla=1;
            else
                maxtrange=hightlist{plotlistindex(1)};   %convert range to index
                mintrange=lowtlist{plotlistindex(1)};
                pooeyla=isempty(lifefit{plotlistindex(1)});
            end
            fitselect=get(hfitpopup,'Value');  % range in time or radians
            
                     
            if fitselect==1
                hplot=plot(pulsedata{plotlistindex(1)}(:,1),pulsedata{plotlistindex(1)}(:,4));
            elseif fitselect==2
                hplot=plot(pulsedata{plotlistindex(1)}(:,1),pulsedata{plotlistindex(1)}(:,5));
            else
                hplot=plot(pulsedata{plotlistindex(1)}(:,4),pulsedata{plotlistindex(1)}(:,5)); %dRdtheat
            end
            hold on
            if pooeyla==0; %% for just viewing data before plot
                if fitselect==1
                    maxresp(plotlistindex(1))=max(pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),4)); %max resp in fit
                    minresp(plotlistindex(1))=min(pulsedata{plotlistindex(1)}(lengtht/10 +int32(mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),4));
                    axis([mintrange maxtrange minresp(plotlistindex(1)) maxresp(plotlistindex(1))]); %% range to plot fit
                    plot(lifefit{plotlistindex(1)});
                    set(findobj(gca,'Type','line'),'LineWidth',2);
                    xlabel('time (us)');
                    ylabel('Phase Response (rad)');
                elseif fitselect==2
                    maxresp(plotlistindex(1))=max(pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),5)); %max resp in fit
                    minresp(plotlistindex(1))=min(pulsedata{plotlistindex(1)}(lengtht/10 +int32(mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),5));
                    axis([mintrange maxtrange minresp(plotlistindex(1)) maxresp(plotlistindex(1))]); %% range to plot fit
                    plot(lifefit{plotlistindex(1)});
                    set(findobj(gca,'Type','line'),'LineWidth',2);
                    xlabel('time (us)');
                    ylabel('Radius Response (R/<R>)');
                else
                    minbla=int32(lengtht/10 +mintrange/deltat)
                    maxbla=int32(lengtht/10 +maxtrange/deltat)
                  %  axis([pulsedata{plotlistindex(1)}(maxbla,4) pulsedata{plotlistindex(1)}(minbla,4) pulsedata{plotlistindex(1)}(minbla,5) pulsedata{plotlistindex(1)}(maxbla,5)])
                    %axis([pulsedata{plotlistindex(1)}(int32(lengtht/10+maxtrange/deltat),4)-0.1 pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat),4)+0.1 pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat),5)-0.1 pulsedata{plotlistindex(1)}(int32(lengtht/10 +maxtrange/deltat),5)+0.1])
                    minxy=pulsedata{plotlistindex(1)}(int32(lengtht/10+maxtrange/deltat),4);
                    maxxy=pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat),4);
                    dxy=(maxxy-minxy)/4;
                    minyy=pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat),5);
                    maxyy=pulsedata{plotlistindex(1)}(int32(lengtht/10 +maxtrange/deltat),5);
                    dyy=(maxyy-minyy)/4;
                    axis([minxy-2*dxy maxxy+2*dxy minyy-dyy maxyy+dyy])
                    plot(pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),4),p{plotlistindex(1)}(2)+pulsedata{plotlistindex(1)}(int32(lengtht/10 +mintrange/deltat):int32(lengtht/10 +maxtrange/deltat),4).*p{plotlistindex(1)}(1),'r')
                    xlabel('phase')
                    ylabel('Radius')
                end

            end 
                legend(['T=' num2str(Tlist{plotlistindex(1)}) ',Power=' num2str(powerlist{plotlistindex(1)})])    
                figure(hfigmain);
    end

   
    function fit_callback_all(source,eventdata)
        for i=1:size(klist,2);
            set(hlist,'Value',i);
            fit_callback;
        end
    end
        
    function fit_callback(source,eventdata)
    %%fit lifetime
        fitselect=get(hfitpopup,'Value');  % range in time or radians
        radtimeselect=get(hradtimepopup,'Value');  % range in time or radians
        plotlistindex=get(hlist,'Value');        %% get selected plot index
        deltat=pulsedata{plotlistindex(1)}(2,1)-pulsedata{plotlistindex(1)}(1,1); 
        lengtht=size(pulsedata{plotlistindex(1)}(:,1),1); %% size of time data
        if radtimeselect==1   %% radians from angle range text boxes
            maxradrange=str2num(get(hmaxradrange,'String'))  ;  %convert range to index
            minradrange=str2num(get(hminradrange,'String'));
            offset=str2num(get(hmintrange,'String'))/deltat;     % offset by tmin to cut KID ringing
            bobfuck=find(pulsedata{plotlistindex(1)}(int32(lengtht/10+offset):lengtht,4)<maxradrange) ; %array shows where condition is meet in data
            mintrange=pulsedata{plotlistindex(1)}(int32(lengtht/10+offset)+bobfuck(1),1)  %pulse at lengtht/10
            bobfucker=find(pulsedata{plotlistindex(1)}((int32(lengtht/10+offset)+bobfuck(1)):lengtht,4)<minradrange); %shows where minimium fit angle is in data
            maxtrange=pulsedata{plotlistindex(1)}((int32(lengtht/10+offset)+bobfuck(1))+bobfucker(1),1)
        elseif radtimeselect==2   %time from t range text boxes 
            maxtrange=str2num(get(hmaxtrange,'String')) ;   %convert range to index
            mintrange=str2num(get(hmintrange,'String'));
            maxradrange=pulsedata{plotlistindex(1)}(int32(mintrange/deltat)+1,4);
            minradrange=pulsedata{plotlistindex(1)}(int32(maxtrange/deltat)+1,4);
        elseif radtimeselect==3
            maxtrange=str2num(get(hmaxtrange,'String')) ;   %sets time range, creates index
            mintrange=str2num(get(hmintrange,'String'));
            offset=mintrange/deltat;
            maxradrange=str2num(get(hmaxradrange,'String'))  ; %only using max radrange
            bobfuck=find(pulsedata{plotlistindex(1)}(int32(lengtht/10+offset):lengtht,4)<maxradrange) ; %array shows where condition is meet in data
            mintrange=pulsedata{plotlistindex(1)}(int32(lengtht/10+offset)+bobfuck(1),1)  %pulse at lengtht/10
            maxradrange=pulsedata{plotlistindex(1)}(int32(mintrange/deltat)+1,4);
            minradrange=pulsedata{plotlistindex(1)}(int32(maxtrange/deltat)+1,4);
        end
        tdata=pulsedata{plotlistindex(1)}(int32(lengtht/10+mintrange/deltat):int32(lengtht/10+maxtrange/deltat),1);
        if fitselect==1 %phase fit
            rdata=pulsedata{plotlistindex(1)}(int32(lengtht/10+mintrange/deltat):int32(lengtht/10+maxtrange/deltat),4); 
        elseif fitselect==2 %amplitude
            rdata=pulsedata{plotlistindex(1)}(int32(lengtht/10+mintrange/deltat):int32(lengtht/10+maxtrange/deltat),5);
        else
            thdata=pulsedata{plotlistindex(1)}(int32(lengtht/10+mintrange/deltat):int32(lengtht/10+maxtrange/deltat),4);
            rdata=pulsedata{plotlistindex(1)}(int32(lengtht/10+mintrange/deltat):int32(lengtht/10+maxtrange/deltat),5);
        end
        if (fitselect==1)|(fitselect==2) %phase or radius
            seedlife1=str2num(get(hseedlife1,'String'))  ; %% inputed lifetime seeds  
            seedlife2=str2num(get(hseedlife2,'String')) ;
            if get(hsingledoublepopup,'Value')==1  %% single exponent
                if fitselect==1 %phase
                s = fitoptions('Method','NonlinearLeastSquares',...
                               'Upper',[10,90000,1],...
                               'Startpoint',[1 seedlife1 0.001]);
                else %amp
                    s = fitoptions('Method','NonlinearLeastSquares',...
                               'Upper',[10,90000,1],...
                               'Startpoint',[-1 seedlife1 1]);
                end
                
                f = fittype('a*exp(-(x)/b)+c','options',s);
                lifefit{plotlistindex(1)} = fit(tdata,rdata,f);
                if  fitselect==1 %phase fit
                    lifefit{plotlistindex(1)}
                    life2list{plotlistindex(1)}=[];
                    alist1{plotlistindex(1)}=lifefit{plotlistindex(1)}.a;
                    lifelist{plotlistindex(1)}=lifefit{plotlistindex(1)}.b; %% lifetime in fit
                    set(hlifelist,'String',lifelist);
                else % amp
                    lifefit{plotlistindex(1)}
                    %life2list{plotlistindex(1)}=[];
                    alist2{plotlistindex(1)}=lifefit{plotlistindex(1)}.a;
                    life2list{plotlistindex(1)}=lifefit{plotlistindex(1)}.b; %% amp lifetime in 2nd exp fit
                    set(hlife2list,'String',life2list);
                end
            else                                    %%double exponent
                s = fitoptions('Method','NonlinearLeastSquares',...
                               'Startpoint',[1 seedlife1 0.001 1 seedlife2]);
                f = fittype('a*exp(-(x)/b)+c+d*exp(-(x)/e)','options',s);
                lifefit{plotlistindex(1)} = fit(tdata,rdata,f)
                life2list{plotlistindex(1)}=lifefit{plotlistindex(1)}.e; %% life2time in fit
                alist1{plotlistindex(1)}=lifefit{plotlistindex(1)}.a;
                alist2{plotlistindex(1)}=lifefit{plotlistindex(1)}.d;
                lifelist{plotlistindex(1)}=lifefit{plotlistindex(1)}.b; %% lifetime in fit
                set(hlifelist,'String',lifelist);
                set(hlife2list,'String',life2list);
            end
        else
            p{plotlistindex(1)}=polyfit(thdata,rdata,1);
            drdthlist{plotlistindex(1)}=p{plotlistindex(1)}(1);
            set(hdrdthlist,'String',drdthlist);
        end
        
        %sort out lists
        hightlist{plotlistindex(1)}=maxtrange;
        lowtlist{plotlistindex(1)}=mintrange;
        lowradlist{plotlistindex(1)}=minradrange;
        highradlist{plotlistindex(1)}=maxradrange;
        set(hhightlist,'String',hightlist);
        set(hlowtlist,'String',lowtlist);
        set(hlowradlist,'String',lowradlist);
        set(hhighradlist,'String',highradlist);
        update_callback;
    end

    function exportinitialise_callback(source,eventdata)
        exportselect=get(hexportpopup,'Value');
        
        if exportselect==1 %% export mean
            resultfile=uiputfile('.txt','Save lifetime fits','lifefit.txt')
            casewrite('Kid, T, mean lifetime, std dev, mean life2time, std dev, amplitude1, std dev, amplitude2 ,std dev, Powerreadout',resultfile);
        elseif exportselect==2
            resultfile=uiputfile('.txt','Save lifetime fits','lifefit.txt')
            casewrite('Kid, T, lifetime, life2time, PowerKID, Powerhigh, low t (us), high t (us), low angle, high angle, amplitude1, amplitude2',resultfile);
        elseif exportselect==3 
            resultfile=uiputfile('.csv','Save dRdtheta fits','dRdthetalife.csv')
            casewrite('KID,  dRdtheta, std, T[K], readpower, pinternal, Q, Qc, Qi, f0, phase life, std, amp life, std, Nqp, Psky',resultfile);
        end
    end
    function exportparam_callback(source,eventdata)
 
        if exportselect==1 %% export mean
            dlmwrite(resultfile, [mean(cell2mat(kidlist)) mean(cell2mat(Tlist)) mean(cell2mat(lifelist)) std(cell2mat(lifelist)) mean(cell2mat(life2list)) std(cell2mat(life2list)) mean(cell2mat(alist1)) std(cell2mat(alist1)) mean(cell2mat(alist2)) std(cell2mat(alist2)) -mean(cell2mat(powerlist))], '-append','newline', 'pc')
        elseif exportselect==2
            values=[transpose(kidlist),transpose(Tlist),transpose(lifelist),transpose(life2list),transpose(powerlist),transpose(VNAlist),transpose(lowtlist),transpose(hightlist),transpose(lowradlist),transpose(highradlist),transpose(alist1),transpose(alist2)]
            dlmwrite(resultfile, values, '-append','newline', 'pc')
        elseif exportselect==3
            bla=size(kidlist,2);
            dlmwrite(resultfile, [mean(cell2mat(kidlist))  mean(cell2mat(drdthlist)) std(cell2mat(drdthlist)) mean(cell2mat(Tlist)) -mean(cell2mat(powerlist)) 0 0 0 0 0 mean(cell2mat(lifelist)) std(cell2mat(lifelist)) mean(cell2mat(life2list)) std(cell2mat(life2list)) 0 0],'-append','newline', 'pc', 'precision', '%.6g','delimiter', ',')
        end
                
    end

end %% main end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% subfunctions%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kidn]=kidnumber(klist)
%%% read kid number from file name

    for n=1:length(klist)
        place(n)=regexpi(klist{n},'KID','end','once');%after which number follows
        kidn{n}=[cell2mat(textscan(klist{n},['%*' num2str(place(n)) 'c %f %*s']))];
    end
end

function [red] = remfun(org,remindex)
%%%reindex arrays, removing cells specified by remindex
    red{1}=[];
    
    if isempty(org{1})==0 %%ignore empt arrays
        for loopy=1:size(remindex,2);
            org{remindex(loopy)}=[]; % delete original cells
        end
        n=1;
        for fucky=1:size(org,2)    
            if isempty(org{fucky})==0 %% if orginal not empty
                red{n}=org{fucky};   %reindex
                n=n+1;
            end
        end
   end
end

function [fnoise]=freqnoise(phsnoise,Q,f0)
    for ii=1:size(Q,1)
        fnoise{ii}=phsnoise+20*log10(f0(ii)*1e9/4./Q(ii));
        Q
        f0
        +20*log10(f0(ii)/4./Q(ii))
    end
end


function casewrite(strmat,filename)
%CASEWRITE Writes casenames from a string matrix to a file.
%   CASEWRITE(STRMAT,FILENAME) writes a list of names to a file, one per line. 
%   FILENAME is the complete path to the desired file. If FILENAME does not
%   include directory information, the file will appear in the current directory.
%
%   CASEWRITE with just one input displays the File Open dialog box allowing
%   interactive naming of the file.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.9.2.1 $  $Date: 2003/11/01 04:25:22 $

if (nargin == 0)
   error('stats:casewrite:TooFewInputs',...
         'CASEWRITE requires at least one argument.');
end
if nargin == 1
   [F,P]=uiputfile('*');
   filename = [P,F];
end
fid = fopen(filename,'wt');

if fid == -1
   disp('Unable to open file.');
   return
end

if strcmp(computer,'MAC2')
   lf = setstr(13);
else
   lf = setstr(10);
end

lf = lf(ones(size(strmat,1),1),:);
lines  = [strmat lf]';

fprintf(fid,'%s',lines);
fclose(fid);
end

