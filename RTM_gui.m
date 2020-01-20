function RTM_gui

global analysis_gui
times_pushed = 0;

%-------------------------------------------------------------------------%
startstop=false;
critvall=1e-4;
relerr_TT=Inf;
relerr_IF=Inf;
nitmax=100;
[glnodes,glweights]=gauss_legendre(8);
alpha=[1,1];
%-------------------------------------------------------------------------%
t0=0;
%-------------------------------------------------------------------------%
[Vbr,Vbt]=deal(0);
[t,Ca,Ctr,Ct]=deal(0);
[K1rx,K1x,K2x,K3x,K4x,Cxt]=deal(0);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
h = figure('Visible','off','Position',[200,200,850,650],...
    'MenuBar','none','Toolbar','figure');
%-------------------------------------------------------------------------%

% Volumes
%-------------------------------------------------------------------------%
hVpanel = uipanel('Title','Blood volume fractions',...
    'TitlePosition','centertop',...
    'Position',[550/850,550/650,220/850,80/650]);
%-------------------------------------------------------------------------%
hVbrtxt = uicontrol('Style','text','String','Vbr=',...
    'Position',[0,20,35,20],...
    'Parent',hVpanel);
hVbr = uicontrol('Style','edit','String','0.025',...
    'Position',[35,20,70,25],...
    'Parent',hVpanel);
% Vbr = 0.025
%-------------------------------------------------------------------------%
hVbttxt = uicontrol('Style','text','String','Vbt=',...
    'Position',[105,20,35,20],...
    'Parent',hVpanel);
hVbt = uicontrol('Style','edit','String','0.15',...
    'Position',[140,20,70,25],...
    'Parent',hVpanel);
% Vbt = 0.15;
%-------------------------------------------------------------------------%

% TT
%-------------------------------------------------------------------------%
hTTpanel = uipanel('Title','Target Tissue',...
    'TitlePosition','centertop',...
    'Position',[550/850,130/650,220/850,390/650]);
%-------------------------------------------------------------------------%
hTTig = uicontrol('Style','text','String','Initial guess',...
    'Position',[5,340,110,20],...
    'HorizontalAlignment','Left','Parent',hTTpanel);
%-------------------------------------------------------------------------%
hTTigrand = uicontrol('Style','pushbutton','String','Random',...
    'Position',[110,340,105,25],...
    'HorizontalAlignment','Right','Parent',hTTpanel,...
    'Callback',{@TTigrand_Callback});
%-------------------------------------------------------------------------%
hk1rtxt = uicontrol('Style','text','String','K1r=',...
    'Position',[0,310,35,20],...
    'Parent',hTTpanel);
hk1r = uicontrol('Style','edit','String','0.05',...
    'Position',[35,310,70,25],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk2rtxt = uicontrol('Style','text','String','K2r=',...
    'Position',[110,310,35,20],...
    'Parent',hTTpanel);
hk2r = uicontrol('Style','text','String','',...
    'Position',[145,310,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk1txt = uicontrol('Style','text','String','K1=',...
    'Position',[0,280,35,20],...
    'Parent',hTTpanel);
hk1 = uicontrol('Style','edit','String','0.5',...
    'Position',[35,280,70,25],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk2txt = uicontrol('Style','text','String','K2=',...
    'Position',[0,250,35,20],...
    'Parent',hTTpanel);
hk2 = uicontrol('Style','edit','String','0.5',...
    'Position',[35,250,70,25],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk3txt = uicontrol('Style','text','String','K3=',...
    'Position',[110,280,35,20],...
    'Parent',hTTpanel);
hk3 = uicontrol('Style','edit','String','0.1',...
    'Position',[145,280,70,25],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk4txt = uicontrol('Style','text','String','K4=',...
    'Position',[110,250,35,20],...
    'Parent',hTTpanel);
hk4 = uicontrol('Style','edit','String','0.01',...
    'Position',[145,250,70,25],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
% hregtxttit = uicontrol('Style','text','String','Regularization parameter',...
%     'Position',[5,215,210,20],...
%     'Parent',hTTpanel);
% hregtxt = uicontrol('Style','text','String','r=',...
%     'Position',[70,190,30,20],...
%     'Parent',hTTpanel);
% hreg = uicontrol('Style','edit','String','1e6',...
%     'Position',[100,190,70,25],...
%     'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hTTstartstop = uicontrol('Style','pushbutton','String',...
    'Start/Stop iterations','Position',[10,155,200,25],...
    'Parent',hTTpanel,...
    'Callback',{@TTstartstop_Callback});
%-------------------------------------------------------------------------%
hTTnittxt = uicontrol('Style','text','String','Iteration',...
    'Position',[15,130,105,20],...
    'HorizontalAlignment','Right','Parent',hTTpanel);
hTTnit = uicontrol('Style','text','String','',...
    'Position',[125,130,100,20],...
    'HorizontalAlignment','Left','Parent',hTTpanel);
%-------------------------------------------------------------------------%
hTTrec = uicontrol('Style','text','String','Recovered values',...
    'Position',[10,100,200,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk1rtxtrec = uicontrol('Style','text','String','K1r=',...
    'Position',[0,80,35,20],...
    'Parent',hTTpanel);
hk1rrec = uicontrol('Style','text','String','',...
    'Position',[35,80,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk2rtxtrec = uicontrol('Style','text','String','K2r=',...
    'Position',[110,80,35,20],...
    'Parent',hTTpanel);
hk2rrec = uicontrol('Style','text','String','',...
    'Position',[145,80,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk1txtrec = uicontrol('Style','text','String','K1=',...
    'Position',[0,60,35,20],...
    'Parent',hTTpanel);
hk1rec = uicontrol('Style','text','String','',...
    'Position',[35,60,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk2txtrec = uicontrol('Style','text','String','K2=',...
    'Position',[0,40,35,20],...
    'Parent',hTTpanel);
hk2rec = uicontrol('Style','text','String','',...
    'Position',[35,40,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk3txtrec = uicontrol('Style','text','String','K3=',...
    'Position',[110,60,35,20],...
    'Parent',hTTpanel);
hk3rec = uicontrol('Style','text','String','',...
    'Position',[145,60,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hk4txtrec = uicontrol('Style','text','String','K4=',...
    'Position',[110,40,35,20],...
    'Parent',hTTpanel);
hk4rec = uicontrol('Style','text','String','',...
    'Position',[145,40,70,20],...
    'Parent',hTTpanel);
%-------------------------------------------------------------------------%
hrelerrTTtxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,20,210,20],...
    'Parent',hTTpanel);
hrelerrTT = uicontrol('Style','text','String','',...
    'Position',[5,5,210,20],...
    'Parent',hTTpanel);

%-------------------------------------------------------------------------%

% IF
%-------------------------------------------------------------------------%
hIFpanel = uipanel('Title','Input Function',...
    'TitlePosition','centertop',...
    'Position',[550/850,40/650,220/850,60/650]);
%-------------------------------------------------------------------------%
hrelerrIFtxt = uicontrol('Style','text','String','Relative error',...
    'Position',[5,20,210,20],...
    'Parent',hIFpanel);
hrelerrIF = uicontrol('Style','text','String','',...
    'Position',[5,5,210,20],...
    'Parent',hIFpanel);
%-------------------------------------------------------------------------%

%% Callback

%-------------------------------------------------------------------------%
function  TTigrand_Callback(hObject, evendata, handles)
    
    times_pushed = times_pushed + 1;
    
    if times_pushed == 1
        
        Kx=num2cell(rand(1,5));
        [K1rx,K1x,K2x,K3x,K4x]=deal(Kx{:});
        
    else
        
        K1rx = rand(1);
        K1x = rand(1);
        K2x = rand(1);
        K3x = 0.2*rand(1);
        K4x = 0.2*rand(1);
        
    end

    set(hk1r,'String',num2str(K1rx));
    set(hk1,'String',num2str(K1x));
    set(hk2,'String',num2str(K2x));
    set(hk3,'String',num2str(K3x));
    set(hk4,'String',num2str(K4x));
    
end
%-------------------------------------------------------------------------%
function TTstartstop_Callback(hObject, evendata, handles)

    b=-1;
    startstop=~startstop;

    if startstop

        Vbr=str2double(get(hVbr,'String'));
        Vbt=str2double(get(hVbt,'String'));
        
        K1rx=str2double(get(hk1r,'String'));
        K1x=str2double(get(hk1,'String'));
        K2x=str2double(get(hk2,'String'));
        K3x=str2double(get(hk3,'String'));
        K4x=str2double(get(hk4,'String'));
        
        % r=str2double(get(hreg,'String'));
        t = analysis_gui.DATA.ROI.Time_ROI./60; t = t';
            
        Ca = analysis_gui.IF.Averaged_IF;
        Ctr = analysis_gui.RT.Averaged_RT;
        Ct = analysis_gui.DATA.ROI.Averaged_ROI;
        tbool=(t>t0);t=t(tbool);Ca=Ca(tbool);Ctr=Ctr(tbool);Ct=Ct(tbool);

        % Ca and Ctr data (RT) function handle
        Ca=@(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));
        Ctr_smooth = csaps(t,Ctr,1e-2,t); Ctr_smooth(1)=0;  
        Ctr = @(tt)(interp1([0 t],[0 Ctr_smooth],tt,'linear',0));

        % asympthotic Logan plot on Ca, Ctr
        ind_t0_logan_Ctr = compute_logan_t0(t,Ca,Ctr);
        t0_logan_Ctr = t(ind_t0_logan_Ctr);
        [slope_logan_Ctr] = logan_plot_Ctr(Ca,Ctr,t0_logan_Ctr,t(ind_t0_logan_Ctr+1:end),glnodes,glweights);

        % lambda = (1-Vbr)/(slope_Ctr - Vbr)
        lambda_logan = (1-Vbr)/(slope_logan_Ctr - Vbr);
        % sigma = (1-Vbr)/Vbr + lambda
        sigma_logan = (1-Vbr)/Vbr + lambda_logan;

        set(hk2r,'String',num2str(lambda_logan*K1rx));
        
        % Solve DIRECT PROBLEM (with initial parameters)
        Ax = [-(K2x+K3x) K4x; K3x -K4x];
        Cax = concentration_Ca_Ctr(K1rx,Vbr,sigma_logan,Ctr,0,0,t,glnodes,glweights);
        Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
        Cx = concentration_TT(K1x,Ax,Cax,0,[0;0],t,glnodes,glweights);
        Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

        nit=0;
        crit=Inf(6,6);

    end

    while startstop&&any(crit(:)>critvall)&&(nit<nitmax)

        [K1rx,K1x,K2x,K3x,K4x,Cax,Cxt,relerr_TT,relerr_IF,nit,crit,relerrp]=...
            iterate_RTM_data(glnodes,glweights,t,alpha,Vbt,Vbr,sigma_logan,Ct,Ctr,Ca,...
            K1rx,K1x,K2x,K3x,K4x,Cax,Cxt,relerr_TT,nit,crit);
        
        set(hk1rrec,'String',num2str(K1rx));
        set(hk2rrec,'String',num2str(lambda_logan*K1rx));
        set(hk1rec,'String',num2str(K1x));
        set(hk2rec,'String',num2str(K2x));
        set(hk3rec,'String',num2str(K3x));
        set(hk4rec,'String',num2str(K4x));
        set(hrelerrTT,'String',num2str(relerr_TT));
        set(hrelerrIF,'String',num2str(relerr_IF));
        set(hTTnit,'String',num2str(nit));

        set(0,'DefaultAxesFontSize',15)
        set(gcf,'CurrentAxes',hTT);
        plot(t,Ct,'b','Linewidth',2);
        hold on;
        plot(t,Cxt,'g--','Linewidth',1.5);
        xlabel('time [min]'); ylabel('concentration [kBq/cc]');
        title({'Target Tissue'},'Fontsize',15);
        legend('Real data','Reconstructed data','Location','NorthEast');
        hold off;
        drawnow;
        
        set(0,'DefaultAxesFontSize',15)
        set(gcf,'CurrentAxes',hIF);
        plot(t(1:ind_t0_logan_Ctr),Ca(t(1:ind_t0_logan_Ctr)),'r--','Linewidth',2);
        hold on;
        plot(t(ind_t0_logan_Ctr:end),Ca(t(ind_t0_logan_Ctr:end)),'r','Linewidth',2);
        plot(t,Cax(t),'g--','Linewidth',1.5);
        %plot([t(ind_t0_logan_Ctr),t(ind_t0_logan_Ctr)],get(hIF,'YLim'),'k','Linewidth',2)
        xlabel('time [min]'); ylabel('concentration [kBq/cc]');
        title({'Input Function'},'Fontsize',15);
        legend('Real data','(known) asymptotic data','Reconstructed data','Location','NorthEast');
        hold off;
        drawnow;

        b=0;

        if (nit>30)&&(relerrp<relerr_TT)&&(relerrp<0.1)
            startstop =~startstop;
        end

        if (nit>30)&&(relerrp<relerr_TT)&&(relerrp>0.1)
            startstop =~startstop;
            msgbox('Wrong initial guess. Choose other parameters');
            b=1;
        end

    end

    if (b==0)&&(relerrp<0.1)

        msgbox('RTM converged');

        % rename and save variables
        K_parameters = {'K1rx','K1x','K2x','K3x','K4x'};
        K = [K1rx,K1x,K2x,K3x,K4x];

        Ax = [-(K2x+K3x) K4x; K3x -K4x];
        Cax = concentration_Ca_Ctr(K1rx,Vbr,sigma_logan,Ctr,0,0,t,glnodes,glweights);
        Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
        Cx = concentration_TT(K1x,Ax,Cax,0,[0;0],t,glnodes,glweights);
        Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

        TT = analysis_gui.DATA;
        RT = analysis_gui.RT;
        IF = analysis_gui.IF;
        KINETICS.parameters = K_parameters;
        KINETICS.values = K;
        FIT_TT = struct('Fitting_curve',Cxt);
        FIT_TT.Relative_error = relerr_TT;
        FIT_IF = struct('Fitting_curve',Cax(t));
        FIT_IF.Relative_error = relerr_IF;

        format short;
        aux_clock = fix(clock);
        date_time = strcat(num2str(aux_clock(1)),'-',num2str(aux_clock(2)),'-',num2str(aux_clock(3)),'_',num2str(aux_clock(4)),'_',num2str(aux_clock(5)),'_',num2str(aux_clock(6)));
        
        % save .mat file
        save([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_' analysis_gui.DATA.model '_' date_time '.mat'],'TT','RT','IF','KINETICS','FIT_TT','FIT_IF');

        % write .xlsx file
        C = [K_parameters' num2cell(K)'];   
        fid = fopen([analysis_gui.OUTPUTfolder analysis_gui.slash char(analysis_gui.DATA.ROI.Study_Name) '_' analysis_gui.DATA.model '_' date_time '.xls'],'w');
            formatSpec = '%s\t%s\r\n';
            fprintf(fid,formatSpec,'Names','Values');
            formatSpec2 = '%s\t%d\r\n';
            for row = 1:length(K)
                fprintf(fid,formatSpec2,C{row,:});
            end
            fclose(fid);
        
    end

     if (b==0)&&(relerrp>0.1)
       msgbox('Wrong initial guess. Choose other parameters');
     end

    % If stop has not been pressed, that is, the algorithm stopped
    % because of the criterion, we change the startstopl here
    if startstop
        startstop=~startstop;
    end

end

%-------------------------------------------------------------------------%
hTT=axes('Units','Pixels','Position',[73.8 370 434 250]);
hIF=axes('Units','Pixels','Position',[73.8 60 434 250]);
%-------------------------------------------------------------------------%
set([h,...  
    hVpanel,hVbrtxt,hVbr,hVbttxt,hVbt,...
    hTT,hTTpanel,hTTig,...
    hk1rtxt,hk1r,hk2rtxt,hk2r,... 
    hk1rtxtrec,hk1rrec,hk2rtxtrec,hk2rrec,...
    hk1txt,hk1,hk2txt,hk2,hk3txt,hk3,hk4txt,hk4,hTTigrand,...
    hTTstartstop,hTTrec,hk1txtrec,hk1rec,hk2txtrec,hk2rec,...
    hk3txtrec,hk3rec,hk4txtrec,hk4rec,hTTnittxt,hTTnit,...
    hrelerrTTtxt,hrelerrTT,... %hregtxttit,hregtxt,hreg,...
    hIF,hIFpanel,hrelerrIFtxt,hrelerrIF],...
    'Units','normalized');
%-------------------------------------------------------------------------%
set(h,'Name','RTM GUI');
set(h,'NumberTitle','off');
movegui(h,'center');
%-------------------------------------------------------------------------%
set(h,'Visible','on');
%-------------------------------------------------------------------------%

end