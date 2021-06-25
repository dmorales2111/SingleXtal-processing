function [] = XTALFIT2_scale_180(datax,datay,dataz)
    %Specalized Multiple spectral fitting GUI - by P. Stallworth 10/16.  
    %Input for n data sets, the multidata m X 2n matrix 
    %(data(1,n)= angle, data(i:n)= freq array, data(i:n+1)= intensity) to be fit (SPECTset).  
    %Edit 'fitfunc' (at the end) with the appropriate fitting function.
    %datax, datay, dataz = x-axis, y-axis and z-axis data sets

    global Minit Range2 dimtem
    global theta ymin ymax realdatax realdatay realdataz %Orienitx Orienity Orienitz

    %Edit initial slider, max and min values = [height width Txx Tyy Tzz Txy Txz Tyz Offset];
    Minit=[1 0 0;0 1 0;0 0 1]; %predefine lab-frame tensor
    initvalues = [25 100 0 0 0 0 0 0 0];
    maxvalues = [100 500 5000 5000 5000 2000 2000 2000 180];
    minvalues = [5 10 -5000 -5000 -5000 -3000 -3000 -3000 -180];
    res = 200; %resolution of reduced data set
    fitres = 100; %resolution of fittings must be greater than 1
    minRange=6000;
    maxRange=10000;
    del=(maxRange-minRange)/(fitres-1);
    Range2 = (minRange:del:maxRange)'; %defines the range and ppm resolution of the fit
    theta=datax(1,:);
    dimtem=[res 26]; %dimtemp(1) = res rows, dimtemp(2) = 26 cols
    ymin = -1;
    ymax = theta(length(theta))+16;
    datheight = 15; %data spectral intensity adjust
    n=9; %to pick-out baseline points for renormalization

    %Normalize the data and add DC
    tempdatax = zeros(res,26); %preallocate truncated data sets
    tempdatay = zeros(res,26); %preallocate truncated data sets
    tempdataz = zeros(res,26); %preallocate truncated data sets
    for i = 1:13 %renormalizes intensity data to fit between 0 and 360 
        tempdatax(:,2*i-1:2*i) = SECTRES(datax(2:end,2*i-1:2*i),minRange,maxRange,res); 
        basdat = sum(tempdatax(res-n:res,2*i))/10;
        tempdatax(:,2*i) = tempdatax(:,2*i) - basdat;
        maxdat = max(tempdatax(:,2*i));
        if maxdat > 0 %conditional for cases where maxdat = 0
            tempdatax(:,2*i) = datheight*tempdatax(:,2*i)/maxdat + theta(2*i); %normalize and add DC
        else
            tempdatax(:,2*i) = datheight*tempdatax(:,2*i) + theta(2*i); %just add DC
        end
    end
    for i = 1:13
        tempdatay(:,2*i-1:2*i) = SECTRES(datay(2:end,2*i-1:2*i),minRange,maxRange,res);
        basdat = sum(tempdatay(res-n:res,2*i))/10;
        tempdatay(:,2*i) = tempdatay(:,2*i) - basdat;
        maxdat = max(tempdatay(:,2*i));
        if maxdat > 0 %conditional for cases where maxdat = 0
            tempdatay(:,2*i) = datheight*tempdatay(:,2*i)/maxdat + theta(2*i); %normalize and add DC
        else
            tempdatay(:,2*i) = datheight*tempdatay(:,2*i) + theta(2*i); %just add DC
        end
    end
    for i = 1:13
        tempdataz(:,2*i-1:2*i) = SECTRES(dataz(2:end,2*i-1:2*i),minRange,maxRange,res);
        basdat = sum(tempdataz(res-n:res,2*i))/10;
        tempdataz(:,2*i) = tempdataz(:,2*i) - basdat;
        maxdat = max(tempdataz(:,2*i));
        if maxdat > 0 %conditional for cases where maxdat = 0
            tempdataz(:,2*i) = datheight*tempdataz(:,2*i)/maxdat + theta(2*i); %normalize and add DC
        else
            tempdataz(:,2*i) = datheight*tempdataz(:,2*i) + theta(2*i); %just add DC
        end
    end

    realdatax=tempdatax;
    realdatay=tempdatay;
    realdataz=tempdataz;

    %initiate the figure
    scrsz = get(0,'ScreenSize'); %initiate figure 
    hui.fh = figure('units','pixels','position',[0.05*scrsz(3) 0.05*scrsz(4) 0.85*scrsz(3) 0.85*scrsz(4)],...
           'name','XTALFIT2_scale',...  %'menubar','none',...
           'numbertitle','off','ToolBar','none',...
           'resize','on');
       
    fhsize = get(hui.fh,'Position');
    fhheight = fhsize(4);
    fhwidth = fhsize(3);
    fhleft = fhsize(1);
    fhbottom = fhsize(2);

    subplot('position',[(fhleft+25)/fhwidth (fhbottom+50)/fhheight (fhwidth*0.2)/fhwidth fhheight*0.70/fhheight])
    %plot(realdatax(:,1),realdatax(:,2),'.k')
    set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180])
    axis([min(Range2),max(Range2),ymin,ymax])
    hold on %note that "hold on" and "hold off" do not produce messages in the command window
   
    for i=1:dimtem(2)/2
        j=2*i-1;
        plot(realdatax(:,j),realdatax(:,j+1),'.k')
    end
    
    hold off
    xlabel('ppm') 
    ylabel('rotation angle (THETA)')
    title('Rotation about Xtal X-Axis')

    subplot('position',[(fhleft+25+fhwidth*0.2+65)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
    %plot(realdatay(:,1),realdatay(:,2),'.k')
    set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180])
    axis([min(Range2),max(Range2),ymin,ymax])
    hold on 
    
    for i=1:dimtem(2)/2
        j=2*i-1;
        plot(realdatay(:,j),realdatay(:,j+1),'.k')
    end
    
    hold off
    xlabel('ppm') 
    ylabel('rotation angle (THETA)')
    title('Rotation about Xtal Y-Axis')

    subplot('position',[(fhleft+25+2*fhwidth*0.2+130)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
    %plot(realdataz(:,1),realdataz(:,2),'.k')
    set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180])
    axis([min(Range2),max(Range2),ymin,ymax])
    hold on 
    
    for i=1:dimtem(2)/2
        j=2*i-1;
        plot(realdataz(:,j),realdataz(:,j+1),'.k')
    end
    
    hold off
    xlabel('ppm') 
    ylabel('rotation angle (THETA)')
    title('Rotation about Xtal Z-Axis')

    %For example, hui looks like this:
    %hui = 
    %
    %     sl: [209.0018 2.0020 5.0020 8.0020 11.0020 14.0020 17.0020]
    %    txt: [0.0020 3.0020 6.0020 9.0020 12.0020 15.0020 210.0018]
    %     ed: [1.0020 4.0020 7.0020 10.0020 13.0020 16.0020 211.0018]

    %slider1 = height
    hui.txt(1) = uicontrol('style','text','string','fit height','Units','normalized',...
        'Position',[(fhleft+25+fhwidth*0.2+65)/fhwidth (fhbottom-15)/fhheight 0.04 0.02]);
    hui.ed(1) = uicontrol('style','edit','Units','normalized',...
        'string',num2str(initvalues(1)),...
        'position',[(fhleft+25+fhwidth*0.2+65+fhwidth*0.04)/fhwidth (fhbottom-15)/fhheight 0.04 0.02]);
    hui.sl(1)=uicontrol('style','slider','units','normalized',...
        'Position',[(fhleft+25+fhwidth*0.2+65+fhwidth*0.04+fhwidth*0.04)/fhwidth (fhbottom-15)/fhheight 0.1 0.02],...
        'sliderstep',[0.015 1],...
        'value',initvalues(1),...%value has to be within min/max range
        'min',minvalues(1),'max',maxvalues(1)); 

    %slider2 = width
    hui.txt(2) = uicontrol('style','text','string','fit width','units','normalized',...
        'position',[(fhleft+25+2*fhwidth*0.2+130)/fhwidth (fhbottom-15)/fhheight 0.04 0.02]);
    hui.ed(2) = uicontrol('style','edit','units','normalized',...
        'string',num2str(initvalues(2)),...
        'position',[(fhleft+25+2*fhwidth*0.2+130+fhwidth*0.04)/fhwidth (fhbottom-15)/fhheight 0.04 0.02]);
    hui.sl(2)=uicontrol('style','slider','units','normalized',...
    'position',[(fhleft+25+2*fhwidth*0.2+130+fhwidth*0.08)/fhwidth (fhbottom-15)/fhheight 0.1 0.02],...
    'sliderstep',[0.01 1],...
    'value',initvalues(2),...%value has to be within min/max range
    'min',minvalues(2),'max',maxvalues(2)); 

 %%% Sliders 3 through 11, for rotation matrix elements

    % First obtain position vector for all the hui.ax elements and use this data to
    % build the sliders, pushbuttons, and text boxes

%     hui1pos = get(ax1,'Position');
    ax1h = fhheight*0.70;
    ax1w = fhwidth*0.2;
    ax1l = fhleft+25;
    ax1b = fhbottom+50;

%     hui2pos = get(ax2,'Position');
    ax2h = fhheight*0.70;
    ax2w = fhwidth*0.2;
    ax2l = fhleft+25+fhwidth*0.2+65;
    ax2b = fhbottom+50;

%     hui3pos = get(ax3,'Position');
    ax3h = fhheight*0.70;
    ax3w = fhwidth*0.2;
    ax3l = fhleft+25+2*fhwidth*0.2+130;
    ax3b = fhbottom+50;

    %slider3 = Txx
    hui.txt(3) = uicontrol('style','text','string','Txx','units','normalized',...
        'position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(3) = uicontrol('style','edit','string',num2str(0),'units','normalized',...
        'position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(3) = uicontrol('style','slider','units','normalized',...
        'position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(3),...%value has to be within min/max range
    'min',minvalues(3),'max',maxvalues(3));

    %slider4 = Tyy
    hui.txt(4) = uicontrol('style','text','string','Tyy','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(4) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(4) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(4),...%value has to be within min/max range
    'min',minvalues(4),'max',maxvalues(4));

    %slider5 = Tzz
    hui.txt(5) = uicontrol('style','text','string','Tzz','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(5) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(5) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(5),...%value has to be within min/max range
    'min',minvalues(5),'max',maxvalues(5));

    %slider6 = Txy
    hui.txt(6) = uicontrol('style','text','string','Txy','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(6) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(6) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(6),...%value has to be within min/max range
    'min',minvalues(6),'max',maxvalues(6));

    %slider7 = Txz
    hui.txt(7) = uicontrol('style','text','string','Txz','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(7) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(7) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(7),...%value has to be within min/max range
    'min',minvalues(7),'max',maxvalues(7));

    %slider8 = Tyz
    hui.txt(8) = uicontrol('style','text','string','Tyz','units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(8) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(8) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(8),...%value has to be within min/max range
    'min',minvalues(8),'max',maxvalues(8));

    %slider9 = 2nd tensor alpha
    hui.txt(9) = uicontrol('style','text','string','ALF','units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(9) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(9) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+55)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(9),...%value has to be within min/max range
    'min',minvalues(9),'max',maxvalues(9));

    %slider10 = 2nd tensor beta
    hui.txt(10) = uicontrol('style','text','string','BET','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.ed(10) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(10) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+80+ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(9),...%value has to be within min/max range
    'min',minvalues(9),'max',maxvalues(9));

    %slider11 = 2nd tensor gamma
    hui.txt(11) = uicontrol('style','text','string','GAM','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %'fontsize',10);
    hui.ed(11) = uicontrol('style','edit','string',num2str(0),'units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05)/fhheight ax3w*0.15/fhwidth ax3h*0.05/fhheight]); %,'fontsize',10);
    hui.sl(11) = uicontrol('style','slider','units','normalized','position',[(ax3l+ax3w+105+2*ax3w*0.15)/fhwidth (ax3b+ax3h-ax3h*0.25-ax3h*0.05-45-ax3h*0.25-45-2*ax3h*0.05-ax3h*0.25)/fhheight ax3w*0.15/fhwidth ax3h*0.25/fhheight],'sliderstep',[0.005 0.1],...
    'value',initvalues(9),...%value has to be within min/max range
    'min',minvalues(9),'max',maxvalues(9));


    %Euler angles for each orientation
    %x-angles
    hui.txta(1) = uicontrol('style','text','string','alpha','units','normalized','position',[(ax1l+ax1w-0.50*ax1w)/fhwidth (ax1h+ax1b+115)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(1) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax1l+ax1w-0.25*ax1w)/fhwidth (ax1h+ax1b+115)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(2) = uicontrol('style','text','string','beta','units','normalized','position',[(ax1l+ax1w-0.50*ax1w)/fhwidth (ax1h+ax1b+90)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(2) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax1l+ax1w-0.25*ax1w)/fhwidth (ax1h+ax1b+90)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(3) = uicontrol('style','text','string','gamma','units','normalized','position',[(ax1l+ax1w-0.50*ax1w)/fhwidth (ax1h+ax1b+65)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(3) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax1l+ax1w-0.25*ax1w)/fhwidth (ax1h+ax1b+65)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(10) = uicontrol('style','text','string','theta','units','normalized','position',[(ax1l+ax1w-0.50*ax1w)/fhwidth (ax1h+ax1b+40)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(10) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax1l+ax1w-0.25*ax1w)/fhwidth (ax1h+ax1b+40)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    %y-angles
    hui.txta(4) = uicontrol('style','text','string','alpha','units','normalized','position',[(ax2l+ax2w-0.5*ax2w)/fhwidth (ax2h+ax2b+115)/fhheight 0.25*ax2w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(4) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax2l+ax2w-0.25*ax2w)/fhwidth (ax2h+ax2b+115)/fhheight 0.25*ax2w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(5) = uicontrol('style','text','string','beta','units','normalized','position',[(ax2l+ax2w-0.5*ax2w)/fhwidth (ax1h+ax1b+90)/fhheight 0.25*ax2w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(5) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax2l+ax2w-0.25*ax2w)/fhwidth (ax1h+ax1b+90)/fhheight 0.25*ax2w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(6) = uicontrol('style','text','string','gamma','units','normalized','position',[(ax2l+ax2w-0.5*ax2w)/fhwidth (ax1h+ax1b+65)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(6) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax2l+ax2w-0.25*ax2w)/fhwidth (ax1h+ax1b+65)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(11) = uicontrol('style','text','string','theta','units','normalized','position',[(ax2l+ax2w-0.5*ax2w)/fhwidth (ax1h+ax1b+40)/fhheight 0.25*ax1w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(11) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax2l+ax2w-0.25*ax2w)/fhwidth (ax1h+ax1b+40)/fhheight 0.25*ax1w/fhwidth 25/fhheight]); %,'fontsize',10);

    % z-angles
    hui.txta(7) = uicontrol('style','text','string','alpha','units','normalized','position',[(ax3l+ax3w-0.5*ax3w)/fhwidth (ax3h+ax3b+115)/fhheight 0.25*ax3w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(7) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax3l+ax3w-0.25*ax3w)/fhwidth (ax3h+ax3b+115)/fhheight 0.25*ax3w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(8) = uicontrol('style','text','string','beta','units','normalized','position',[(ax3l+ax3w-0.5*ax3w)/fhwidth (ax3h+ax3b+90)/fhheight 0.25*ax3w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(8) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax3l+ax3w-0.25*ax3w)/fhwidth (ax3h+ax3b+90)/fhheight 0.25*ax3w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(9) = uicontrol('style','text','string','gamma','units','normalized','position',[(ax3l+ax3w-0.5*ax3w)/fhwidth (ax3h+ax3b+65)/fhheight 0.25*ax3w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(9) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax3l+ax3w-0.25*ax3w)/fhwidth (ax3h+ax3b+65)/fhheight 0.25*ax3w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.txta(12) = uicontrol('style','text','string','theta','units','normalized','position',[(ax3l+ax3w-0.5*ax3w)/fhwidth (ax3h+ax3b+40)/fhheight 0.25*ax3w/fhwidth 25/fhheight],'fontsize',10);
    hui.txtb(12) = uicontrol('style','text','string',num2str(0),'units','normalized','position',[(ax3l+ax3w-0.25*ax3w)/fhwidth (ax3h+ax3b+40)/fhheight 0.25*ax3w/fhwidth 25/fhheight]); %,'fontsize',10);

    hui.pb(1) = uicontrol('Style','push',... 
                     'Units','pixels',...
                     'units','normalized','position',[ax1l/fhwidth (ax1b+ax1h+40+75/2)/fhheight ax1w*0.3/fhwidth 35/fhheight],...
                     'string', 'X Fit In',...
                     'fontsize',8,'BackgroundColor',[1 0.5 0],'callback',{@pb_load_xdata,hui}); %

    hui.pb(2) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax1l/fhwidth (ax1b+ax1h+40)/fhheight ax1w*0.3/fhwidth 35/fhheight],...
                     'string', 'X Fit Out',...
                     'fontsize',8,'BackgroundColor',[0 0.5 0],'callback',{@pb_dump_xdata,hui}); %

    hui.pb(3) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax2l/fhwidth (ax2b+ax2h+40+75/2)/fhheight ax2w*0.3/fhwidth 35/fhheight],...
                     'string', 'Y Fit In',...
                     'fontsize',8,'BackgroundColor',[1 0.5 0],'callback',{@pb_load_ydata,hui}); %

    hui.pb(4) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax2l/fhwidth (ax2b+ax2h+40)/fhheight ax2w*0.3/fhwidth 35/fhheight],...
                     'string', 'Y Fit Out',...
                     'fontsize',8,'BackgroundColor',[0 0.5 0],'callback',{@pb_dump_ydata,hui}); %

    hui.pb(5) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax3l/fhwidth (ax3b+ax3h+40+75/2)/fhheight ax3w*0.3/fhwidth 35/fhheight],...
                     'string', 'Z Fit In',...
                     'fontsize',8,'BackgroundColor',[1 0.5 0],'callback',{@pb_load_zdata,hui}); %

    hui.pb(6) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax3l/fhwidth (ax3b+ax3h+40)/fhheight ax3w*0.3/fhwidth 35/fhheight],...
                     'string', 'Z Fit Out',...
                     'fontsize',8,'BackgroundColor',[0 0.5 0],'callback',{@pb_dump_zdata,hui}); % 

    hui.pb(7) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[ax1l/fhwidth (ax1b-80)/fhheight ax1w*0.3/fhwidth 35/fhheight],...
                     'string', 'All Fits In',...
                     'fontsize',8,'BackgroundColor','c','callback',{@pb_load_alldata,hui}); %                  

    hui.pb(8) = uicontrol('style','push',... 
                     'units','pix',...
                     'units','normalized','posit',[(ax1l+ax1w*0.3+10)/fhwidth (ax1b-80)/fhheight ax1w*0.3/fhwidth 35/fhheight],...
                     'string', 'All Fits Out',...
                     'fontsize',8,'BackgroundColor',[1 0.75 0],'callback',{@pb_dump_alldata,hui});

    set([hui.sl;hui.txt;hui.ed],'callback',{@update_fitb,hui}); %this line does it all!

end %XTALFIT
%**************************************************************************
function update_fitb(varargin)
global Minit Range2 dimtem R2x R2y R2z
global theta ymin ymax realdatax realdatay realdataz Orienitx Orienity Orienitz
[h,STR] = varargin{[1,3]};  % Get the [address,structure].

switch h  % Who called, editbox or slider? Gives same value to slider and editbox
    case STR.ed(1) %if 1st edit box called
        L = get(STR.sl(1),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(1),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(1) %if slider called
        set(STR.ed(1),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(2) %if 2nd edit box called
        L = get(STR.sl(2),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(2),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(2) %if slider called
        set(STR.ed(2),'string',get(h,'value')) % Set edit box value to current slider value.
  
    case STR.ed(3) %if 3rd edit box called
        L = get(STR.sl(3),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(3),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(3) %if slider called
        set(STR.ed(3),'string',get(h,'value')) % Set edit box value to current slider value.
 
    case STR.ed(4) %if 4th edit box called
        L = get(STR.sl(4),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(4),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(4) %if slider called
        set(STR.ed(4),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(5) %if 5th edit box called
        L = get(STR.sl(5),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(5),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(5) %if slider called
        set(STR.ed(5),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(6) %if 6th edit box called
        L = get(STR.sl(6),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(6),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(6) %if slider called
        set(STR.ed(6),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(7) %if 7th edit box called
        L = get(STR.sl(7),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(7),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(7) %if slider called
        set(STR.ed(7),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(8) %if 8th edit box called
        L = get(STR.sl(8),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(8),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(8) %if slider called
        set(STR.ed(8),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(9) %if 9th edit box called
        L = get(STR.sl(9),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(9),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(9) %if slider called
        set(STR.ed(9),'string',get(h,'value')) % Set edit box value to current slider value.

      case STR.ed(10) %if 10th edit box called
        L = get(STR.sl(10),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(10),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(10) %if slider called
        set(STR.ed(10),'string',get(h,'value')) % Set edit box value to current slider value.

    case STR.ed(11) %if 11th edit box called
        L = get(STR.sl(11),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2} % Set slider value to current edit box value.
            set(STR.sl(11),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set edit box to value out of slider range, 
        end                      % so edit box value remains unchanged. 
    case STR.sl(11) %if slider called
        set(STR.ed(11),'string',get(h,'value')) % Set edit box value to current slider value.
              
    otherwise
        % Do nothing, or whatever.
end

Tenval=Minit;
div1 = get(STR.sl(1),'value'); % get the 1st    "     " = height
div2 = get(STR.sl(2),'value'); %   "  "  2nd    "     " = width

Tenval(1,1) = get(STR.sl(3),'value'); %   "  "  3rd    "     " = Txx
Tenval(2,2) = get(STR.sl(4),'value'); %   "  "  4th    "     " = Tyy
Tenval(3,3) = get(STR.sl(5),'value'); %   "  "  5th    "     " = Tzz
Tenval(1,2) = get(STR.sl(6),'value'); %   "  "  6th    "     " = Txy
Tenval(2,1) = Tenval(1,2);
Tenval(1,3) = get(STR.sl(7),'value'); %   "  "  7th    "     " = Txz
Tenval(3,1) = Tenval(1,3);
Tenval(2,3) = get(STR.sl(8),'value'); %   "  "  8th    "     " = Tyz
Tenval(3,2) = Tenval(2,3);

alf = get(STR.sl(9),'value'); %  = alpha (new tensor alpha)
bet = get(STR.sl(10),'value'); % = beta (new tensor beta)
gam = get(STR.sl(11),'value'); % = gamma (new tensor gamma)
offsetx = str2double(get(STR.txtb(10),'str')); % = x-offset
offsety = str2double(get(STR.txtb(11),'str')); % = y-offset
offsetz = str2double(get(STR.txtb(12),'str')); % = z-offset

%call the main fitting function and plot 
%create the temporary variable arrays: tempx tempy and tempz 
%use inputs = height,width,alpha,beta,gamma,offset (theta),Xorient

scrsz = get(0,'ScreenSize'); %initiate figure
fhheight = 0.85*scrsz(4);
fhwidth = 0.85*scrsz(3);
fhleft = 0.05*scrsz(3);
fhbottom = 0.05*scrsz(4);

tempx=updatefig_2(Tenval,Orienitx,R2x,alf,bet,gam,offsetx,0,div1,div2); 
subplot('position',[(fhleft+25)/fhwidth (fhbottom+50)/fhheight (fhwidth*0.2)/fhwidth fhheight*0.70/fhheight])
%subplot('position',[0.05,0.1,0.25,0.75])
plot(tempx(:,1),tempx(:,2)+theta(1),'r',realdatax(:,1),realdatax(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on
for i=2:dimtem(2)/2
    j=2*i-1;
    Thet=theta(j);
    tempx=updatefig_2(Tenval,Orienitx,R2x,alf,bet,gam,offsetx,Thet,div1,div2);
    plot(tempx(:,1),tempx(:,2)+theta(j),'r',realdatax(:,j),realdatax(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal X-Axis')

tempy=updatefig_2(Tenval,Orienity,R2y,alf,bet,gam,offsety,0,div1,div2); 
subplot('position',[(fhleft+25+fhwidth*0.2+65)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
%subplot('position',[0.35,0.1,0.25,0.75])
plot(tempy(:,1),tempy(:,2)+theta(1),'r',realdatay(:,1),realdatay(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on 
for i=2:dimtem(2)/2
    j=2*i-1;
    Thet=theta(j);
    tempy=updatefig_2(Tenval,Orienity,R2y,alf,bet,gam,offsety,Thet,div1,div2);
    plot(tempy(:,1),tempy(:,2)+theta(j),'r',realdatay(:,j),realdatay(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Y-Axis')

tempz=updatefig_2(Tenval,Orienitz,R2z,alf,bet,gam,offsetz,0,div1,div2); 
subplot('position',[(fhleft+25+2*fhwidth*0.2+130)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
%subplot('position',[0.65,0.1,0.25,0.75])
plot(tempz(:,1),tempz(:,2)+theta(1),'r',realdataz(:,1),realdataz(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on
for i=2:dimtem(2)/2
    j=2*i-1;
    Thet=theta(j);
    tempz=updatefig_2(Tenval,Orienitz,R2z,alf,bet,gam,offsetz,Thet,div1,div2);
    plot(tempz(:,1),tempz(:,2)+theta(j),'r',realdataz(:,j),realdataz(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Z-Axis')
end %function update_fitb

%***************************************************
function z = updatefig_2(Min,Orient,Rtweek,af,bt,gm,offs,tht,h,w) %tht is the offset angle for the fixed tensor
%NOTE difference from updatefig_1 in XTALFIT and XTALFIT1:
%tht is the offset angle for the fixed tensor
%Rtweek is the transformation of 'tweeked' Euler angles from XTALFIT1

I0=[1 0 0;0 1 0;0 0 1];
offs=offs*pi/180;
tht=tht*pi/180;

Roff = [1 0 0;0 cos(offs) sin(offs);0 -sin(offs) cos(offs)]; %offset angle for tensor about lab x-axis
Rot = [1 0 0;0 cos(tht) sin(tht);0 -sin(tht) cos(tht)]; %rotations of THETA about lab x-axis

ROT=Roff*Rot; %lab rotation with offset. Roff and Rot commute

Rotwk= Rtweek*Orient; %Rotwk = 1st-step orientation then 2nd-step tweeked 
R1 = ROT*Rotwk; %complete transformation for 1st tensor (lab rotations about x-axis as 3rd-step)
Min1 = R1*Min*(R1\I0); %1st tensor, transformed to fit data

Reuler = EULer_calc(af,bt,gm); %Euler transformation of 2nd tensor
R2 = ROT*Reuler*Rotwk; %complete transformation for 2nd tensor
Min2 = R2*Min*(R2\I0); %2nd tensor (adujsted by Euler angles), transformed to fit data

z = SPECTfit(Min1,Min2,h,w); %set-up for two tensors.  Can do more, just edit SPECTfit to accept more tensors....
end %updatefig_2

%**************************************************************************
function z = EULer_calc(a,b,g)
%computes composit Euler transformation matrix 
a=a*pi/180;
b=b*pi/180;
g=g*pi/180;

RA=[cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];
RB=[cos(b) 0 -sin(b);0 1 0;sin(b) 0 cos(b)];
RG=[cos(g) sin(g) 0;-sin(g) cos(g) 0;0 0 1];
z = RG*RB*RA; %note order of transformation, these don't commute
end

%**************************************************************************
function z = SPECTfit(Ten1,Ten2,h,w) 
%Default set-up for two tensors.  Can do more, just edit to accept more tensors....

global Range2
F=zeros(size(Range2));

f1 = Ten1(3,3); %frequency shift (zz component) of 1st tensor
f2 = Ten2(3,3); %frequency shift (zz component) of 2nd tensor

for i=1:length(Range2)
   F(i) = 0.3*h*exp(-(Range2(i)-f1).^2/(2*w^2))+... 
          0.3*h*exp(-(Range2(i)-f2).^2/(2*w^2)); %additional Gaussian broadened line   
end

z = [Range2,F];
end %function SPECTfit

%**************************************************************************
function [] = pb_load_xdata(varargin)
%Callback for pushbutton x-data in. Gets values for alpha, beta, gamma and offset (THETA), 
%then updates the slider/editboxes with these angles (in deg) from a data structure IN_data.
global Range2 dimtem theta ymin ymax realdatax R2x Orienitx
S = varargin{3}; % Get structure

load OUT_data %get lab-frame x-tal tensor, Euler angles and lab-frame offset angle from OUT_data
Orienitx = OUT_data.orn; %initial x-tal orientation in lab frame
R2x=EULer_calc(OUT_data.ang(1,1),OUT_data.ang(2,1),OUT_data.ang(3,1)); %tweek transformation from XTALFIT1

%update sliders and editboxes to the New Euler angles (from XTALFIT)
    set(S.sl(3),'value',OUT_data.ten(1,1));  % get the value for Txx.
    set(S.ed(3),'str',num2str(OUT_data.ten(1,1)));
    set(S.sl(4),'value',OUT_data.ten(2,2));  % get the value for Tyy.
    set(S.ed(4),'str',num2str(OUT_data.ten(2,2)));
    set(S.sl(5),'value',OUT_data.ten(3,3));  % get the value for Tzz.
    set(S.ed(5),'str',num2str(OUT_data.ten(3,3)));
    set(S.sl(6),'value',OUT_data.ten(1,2));  % get the value for Txy.
    set(S.ed(6),'str',num2str(OUT_data.ten(1,2)));
    set(S.sl(7),'value',OUT_data.ten(1,3));  % get the value for Txz.
    set(S.ed(7),'str',num2str(OUT_data.ten(1,3)));
    set(S.sl(8),'value',OUT_data.ten(2,3));  % get the value for Tyz.
    set(S.ed(8),'str',num2str(OUT_data.ten(2,3)));
       
    set(S.txtb(1),'str',num2str(OUT_data.ang(1,1)));  %alpha
    set(S.txtb(2),'str',num2str(OUT_data.ang(2,1)));  %beta
    set(S.txtb(3),'str',num2str(OUT_data.ang(3,1)));  %gamma   
    set(S.txtb(10),'str',num2str(OUT_data.ang(4,1))); %theta
    
    alf = get(S.sl(9),'value'); %  = alpha (new tensor alpha)
    bet = get(S.sl(10),'value'); % = beta (new tensor beta)
    gam = get(S.sl(11),'value'); % = gamma (new tensor gamma)
    
h=get(S.sl(1),'value'); % get the 1st slider value = height
w=get(S.sl(2),'value'); %   "  "  2nd    "     " = width

%update the figure (calculate new tensor, x,y,z plots and fits)
Fitx = updatefig_2(OUT_data.ten,Orienitx,R2x,alf,bet,gam,OUT_data.ang(4,1),0,h,w);

scrsz = get(0,'ScreenSize'); %initiate figure
fhheight = 0.85*scrsz(4);
fhwidth = 0.85*scrsz(3);
fhleft = 0.05*scrsz(3);
fhbottom = 0.05*scrsz(4);

%subplot('position',[0.05,0.1,0.25,0.75])
subplot('position',[(fhleft+25)/fhwidth (fhbottom+50)/fhheight (fhwidth*0.2)/fhwidth fhheight*0.70/fhheight])
plot(Fitx(:,1),Fitx(:,2)+theta(1),'r',realdatax(:,1),realdatax(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    Fitx=updatefig_2(OUT_data.ten,Orienitx,R2x,alf,bet,gam,OUT_data.ang(4,1),THET,h,w); %call the main fitting function for each theta
    plot(Fitx(:,1),Fitx(:,2)+theta(j),'r',realdatax(:,j),realdatax(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal X-Axis')
end %pb_load_xdata

%***************************************************
function [] = pb_load_ydata(varargin)
%Callback for pushbutton y-data in. Gets values for alpha, beta, gamma and offset (THETA), 
%then updates the slider/editboxes with these angles (in deg) from a data structure IN_data.
global Range2 dimtem theta ymin ymax realdatay R2y Orienity
S = varargin{3}; % Get structure

load OUT_data %get lab-frame x-tal tensor, Euler angles and lab-frame offset angle from OUT_data
Orienity = OUT_data.orn; %initial x-tal orientation in lab frame
R2y=EULer_calc(OUT_data.ang(1,1),OUT_data.ang(2,1),OUT_data.ang(3,1)); %tweek transformation from XTALFIT1

%update sliders and editboxes to the New Euler angles (from XTALFIT)
    set(S.sl(3),'value',OUT_data.ten(1,1));  % get the value for Txx.
    set(S.ed(3),'str',num2str(OUT_data.ten(1,1)));
    set(S.sl(4),'value',OUT_data.ten(2,2));  % get the value for Tyy.
    set(S.ed(4),'str',num2str(OUT_data.ten(2,2)));
    set(S.sl(5),'value',OUT_data.ten(3,3));  % get the value for Tzz.
    set(S.ed(5),'str',num2str(OUT_data.ten(3,3)));
    set(S.sl(6),'value',OUT_data.ten(1,2));  % get the value for Txy.
    set(S.ed(6),'str',num2str(OUT_data.ten(1,2)));
    set(S.sl(7),'value',OUT_data.ten(1,3));  % get the value for Txz.
    set(S.ed(7),'str',num2str(OUT_data.ten(1,3)));
    set(S.sl(8),'value',OUT_data.ten(2,3));  % get the value for Tyz.
    set(S.ed(8),'str',num2str(OUT_data.ten(2,3)));     
     
    set(S.txtb(4),'str',num2str(OUT_data.ang(1,1)));  %alpha
    set(S.txtb(5),'str',num2str(OUT_data.ang(2,1)));  %beta
    set(S.txtb(6),'str',num2str(OUT_data.ang(3,1)));  %gamma    
    set(S.txtb(11),'str',num2str(OUT_data.ang(4,1)));   %theta

    alf = get(S.sl(9),'value'); %  = alpha (new tensor alpha)
    bet = get(S.sl(10),'value'); % = beta (new tensor beta)
    gam = get(S.sl(11),'value'); % = gamma (new tensor gamma)
    
h=get(S.sl(1),'value'); % get the 1st slider value = height
w=get(S.sl(2),'value'); %   "  "  2nd    "     " = width

%update the figure (calculate new tensor, x,y,z plots and fits)
Fity = updatefig_2(OUT_data.ten,Orienity,R2y,alf,bet,gam,OUT_data.ang(4,1),0,h,w);

scrsz = get(0,'ScreenSize'); %initiate figure
fhheight = 0.85*scrsz(4);
fhwidth = 0.85*scrsz(3);
fhleft = 0.05*scrsz(3);
fhbottom = 0.05*scrsz(4);

%subplot('position',[0.35,0.1,0.25,0.75])
subplot('position',[(fhleft+25+fhwidth*0.2+65)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
plot(Fity(:,1),Fity(:,2)+theta(1),'r',realdatay(:,1),realdatay(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on 
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    Fity=updatefig_2(OUT_data.ten,Orienity,R2y,alf,bet,gam,OUT_data.ang(4,1),THET,h,w); %call the main fitting function for each theta
    plot(Fity(:,1),Fity(:,2)+theta(j),'r',realdatay(:,j),realdatay(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Y-Axis')
end %pb_load_ydata

%***************************************************
function [] = pb_load_zdata(varargin)
%Callback for pushbutton z-data in. Gets values for alpha, beta, gamma and offset (THETA), 
%then updates the slider/editboxes with these angles (in deg) from a data structure IN_data.
global Range2 dimtem theta ymin ymax realdataz R2z Orienitz
S = varargin{3}; % Get structure

load OUT_data %get lab-frame x-tal tensor, Euler angles and lab-frame offset angle from OUT_data (structure format)
Orienitz = OUT_data.orn; %initial x-tal orientation in lab frame
R2z=EULer_calc(OUT_data.ang(1,1),OUT_data.ang(2,1),OUT_data.ang(3,1)); %tweek transformation from XTALFIT1

%update sliders and editboxes to the New Euler angles (from XTALFIT)
    set(S.sl(3),'value',OUT_data.ten(1,1));  % get the value for Txx.
    set(S.ed(3),'str',num2str(OUT_data.ten(1,1)));
    set(S.sl(4),'value',OUT_data.ten(2,2));  % get the value for Tyy.
    set(S.ed(4),'str',num2str(OUT_data.ten(2,2)));
    set(S.sl(5),'value',OUT_data.ten(3,3));  % get the value for Tzz.
    set(S.ed(5),'str',num2str(OUT_data.ten(3,3)));
    set(S.sl(6),'value',OUT_data.ten(1,2));  % get the value for Txy.
    set(S.ed(6),'str',num2str(OUT_data.ten(1,2)));
    set(S.sl(7),'value',OUT_data.ten(1,3));  % get the value for Txz.
    set(S.ed(7),'str',num2str(OUT_data.ten(1,3)));
    set(S.sl(8),'value',OUT_data.ten(2,3));  % get the value for Tyz.
    set(S.ed(8),'str',num2str(OUT_data.ten(2,3)));        
    
    set(S.txtb(7),'str',num2str(OUT_data.ang(1,1)));  %alpha
    set(S.txtb(8),'str',num2str(OUT_data.ang(2,1)));  %beta
    set(S.txtb(9),'str',num2str(OUT_data.ang(3,1)));  %gamma
    set(S.txtb(12),'str',num2str(OUT_data.ang(4,1))); %theta

    alf = get(S.sl(9),'value'); %  = alpha (new tensor alpha)
    bet = get(S.sl(10),'value'); % = beta (new tensor beta)
    gam = get(S.sl(11),'value'); % = gamma (new tensor gamma)    
    
h=get(S.sl(1),'value'); % get the 1st slider value = height
w=get(S.sl(2),'value'); %   "  "  2nd    "     " = width

%update the figure (calculate new tensor, x,y,z plots and fits)
Fitz = updatefig_2(OUT_data.ten,Orienitz,R2z,alf,bet,gam,OUT_data.ang(4,1),0,h,w);

scrsz = get(0,'ScreenSize'); %initiate figure
fhheight = 0.85*scrsz(4);
fhwidth = 0.85*scrsz(3);
fhleft = 0.05*scrsz(3);
fhbottom = 0.05*scrsz(4);

%subplot('position',[0.65,0.1,0.25,0.75])
subplot('position',[(fhleft+25+2*fhwidth*0.2+130)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
plot(Fitz(:,1),Fitz(:,2)+theta(1),'r',realdataz(:,1),realdataz(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on 
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    Fitz = updatefig_2(OUT_data.ten,Orienitz,R2z,alf,bet,gam,OUT_data.ang(4,1),THET,h,w);
    plot(Fitz(:,1),Fitz(:,2)+theta(j),'r',realdataz(:,j),realdataz(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Z-Axis')
end %pb_load_zdata

%**************************************************************************
function [] = pb_load_alldata(varargin)
%Get Alldata from previous XTALFIT session. Input the tensor, all orientations Euler angles, and offsets.
global Range2 dimtem theta ymin ymax 
global realdatax realdatay realdataz Orienitx Orienity Orienitz R2x R2y R2z

S = varargin{3}; % Get structure

load Alldata

%angles in degrees - convert to radians
Orienitx = Alldata.orn(:,1:3); %initial x-tal orientation in lab frame
Orienity = Alldata.orn(:,4:6); %initial y-tal orientation in lab frame
Orienitz = Alldata.orn(:,7:9); %initial z-tal orientation in lab frame

offx = Alldata.ang(4,1);
offy = Alldata.ang(4,2);
offz = Alldata.ang(4,3);

alf = Alldata.ang(1,4);
bet = Alldata.ang(2,4);
gam = Alldata.ang(3,4);

%tweek transformations from XTALFIT1
R2x=EULer_calc(Alldata.ang(1,1),Alldata.ang(2,1),Alldata.ang(3,1));
R2y=EULer_calc(Alldata.ang(1,2),Alldata.ang(2,2),Alldata.ang(3,2));
R2z=EULer_calc(Alldata.ang(1,3),Alldata.ang(2,3),Alldata.ang(3,3));

%update sliders and editboxes to the New Euler angles (from XTALFIT)
    set(S.sl(3),'value',Alldata.ten(1,1));  % get the value for Txx.
    set(S.ed(3),'str',num2str(Alldata.ten(1,1)));
    set(S.sl(4),'value',Alldata.ten(2,2));  % get the value for Tyy.
    set(S.ed(4),'str',num2str(Alldata.ten(2,2)));
    set(S.sl(5),'value',Alldata.ten(3,3));  % get the value for Tzz.
    set(S.ed(5),'str',num2str(Alldata.ten(3,3)));
    set(S.sl(6),'value',Alldata.ten(1,2));  % get the value for Txy.
    set(S.ed(6),'str',num2str(Alldata.ten(1,2)));
    set(S.sl(7),'value',Alldata.ten(1,3));  % get the value for Txz.
    set(S.ed(7),'str',num2str(Alldata.ten(1,3)));
    set(S.sl(8),'value',Alldata.ten(2,3));  % get the value for Tyz.
    set(S.ed(8),'str',num2str(Alldata.ten(2,3)));
    
    set(S.sl(9),'value',alf);  % new tensor alpha
    set(S.ed(9),'str',num2str(alf));
    set(S.sl(10),'value',bet);  % new tensor beta
    set(S.ed(10),'str',num2str(bet));
    set(S.sl(11),'value',gam);  % new tensor gamma
    set(S.ed(11),'str',num2str(gam));
    
    h = get(S.sl(1),'value'); % get the 1st    "     " = height
    w = get(S.sl(2),'value'); %   "  "  2nd    "     " = width

    scrsz = get(0,'ScreenSize'); %initiate figure
    fhheight = 0.85*scrsz(4);
    fhwidth = 0.85*scrsz(3);
    fhleft = 0.05*scrsz(3);
    fhbottom = 0.05*scrsz(4);
    
    %x-axis
    set(S.txtb(1),'str',num2str(Alldata.ang(1,1)));  %alpha
    set(S.txtb(2),'str',num2str(Alldata.ang(2,1)));  %beta
    set(S.txtb(3),'str',num2str(Alldata.ang(3,1)));  %gamma  
    set(S.txtb(10),'str',num2str(offx));  %value for lab x-offset about x-axis (THETA in degrees).
    
    %y-axis
    set(S.txtb(4),'str',num2str(Alldata.ang(1,2)));  %alpha
    set(S.txtb(5),'str',num2str(Alldata.ang(2,2)));  %beta
    set(S.txtb(6),'str',num2str(Alldata.ang(3,2)));  %gamma 
    set(S.txtb(11),'str',num2str(offy));  %value for lab x-offset about y-axis (THETA in degrees).
    
    %z-axis
    set(S.txtb(7),'str',num2str(Alldata.ang(1,3)));  %alpha
    set(S.txtb(8),'str',num2str(Alldata.ang(2,3)));  %beta
    set(S.txtb(9),'str',num2str(Alldata.ang(3,3)));  %gamma        
    set(S.txtb(12),'str',num2str(offz));  %value for lab x-offset about z-axis (THETA in degrees).
    
%update the figure (calculate 1st and 2nd tensor, x,y,z plots and fits.)
%Need the original lab frame tensor, 1st crystal orientation, 1st crystal Euler angle
%adjustments, 2nd crystal Euler angles and offset. 
tempx = updatefig_2(Alldata.ten,Orienitx,R2x,alf,bet,gam,offx,0,h,w);
%subplot('position',[0.05,0.1,0.25,0.75])
subplot('position',[(fhleft+25)/fhwidth (fhbottom+50)/fhheight (fhwidth*0.2)/fhwidth fhheight*0.70/fhheight])
plot(tempx(:,1),tempx(:,2)+theta(1),'r',realdatax(:,1),realdatax(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    tempx=updatefig_2(Alldata.ten,Orienitx,R2x,alf,bet,gam,offx,THET,h,w); %call the main fitting function for each theta
    plot(tempx(:,1),tempx(:,2)+theta(j),'r',realdatax(:,j),realdatax(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal X-Axis')    

tempy = updatefig_2(Alldata.ten,Orienity,R2y,alf,bet,gam,offy,0,h,w);
%subplot('position',[0.35,0.1,0.25,0.75])
subplot('position',[(fhleft+25+fhwidth*0.2+65)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
plot(tempy(:,1),tempy(:,2)+theta(1),'r',realdatay(:,1),realdatay(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on 
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    tempy=updatefig_2(Alldata.ten,Orienity,R2y,alf,bet,gam,offy,THET,h,w); %call the main fitting function for each theta
    plot(tempy(:,1),tempy(:,2)+theta(j),'r',realdatay(:,j),realdatay(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Y-Axis')

%tempz=updatefig_2(Tenval,Orienitz,R2z,alf,bet,gam,offsetz,0,div1,div2); 
tempz = updatefig_2(Alldata.ten,Orienitz,R2z,alf,bet,gam,offz,0,h,w);
%subplot('position',[0.65,0.1,0.25,0.75])
subplot('position',[(fhleft+25+2*fhwidth*0.2+130)/fhwidth (fhbottom+50)/fhheight fhwidth*0.2/fhwidth fhheight*0.70/fhheight])
plot(tempz(:,1),tempz(:,2)+theta(1),'r',realdataz(:,1),realdataz(:,2),'.k')
set(gca,'YTick',[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360])
axis([min(Range2),max(Range2),ymin,ymax])
hold on 
for i=2:dimtem(2)/2
    j=2*i-1;
    THET=theta(j);
    tempz = updatefig_2(Alldata.ten,Orienitz,R2z,alf,bet,gam,offz,THET,h,w);
    plot(tempz(:,1),tempz(:,2)+theta(j),'r',realdataz(:,j),realdataz(:,j+1),'.k')
end
hold off
xlabel('ppm') 
ylabel('rotation angle (THETA)')
title('Rotation about Xtal Z-Axis')
    
end %pb_load_alldata
%**************************************************************************

function [] = pb_dump_xdata(varargin)
%update of Euler angles alpha, beta and gamma is not appropriate since they are different for
%for rotations about xtal-x orientation

global Orienitx

S = varargin{3}; % Get structure
delete('OUT_data.mat')

OUT_data.orn = Orienitx; %initial x-tal orientation in lab frame
OUT_data.ang(1,1) = str2double(get(S.txtb(1),'str')); %alpha
OUT_data.ang(2,1) = str2double(get(S.txtb(2),'str')); %beta
OUT_data.ang(3,1) = str2double(get(S.txtb(3),'str')); %gamma
OUT_data.ang(4,1) = str2double(get(S.txtb(10),'str')); % = x-offset (THETA about x-axis)

%OUT_data.ang(1,2) = get(S.sl(9),'value'); %  = (new tensor alpha)
%OUT_data.ang(2,2) = get(S.sl(10),'value'); % = (new tensor beta)
%OUT_data.ang(3,2) = get(S.sl(11),'value'); % = (new tensor gamma)
%OUT_data.ang(4,2) = 0; % = z-offset (THETA about x-axis)

OUT_data.ten(1,1) = get(S.sl(3),'value');
OUT_data.ten(2,2) = get(S.sl(4),'value');
OUT_data.ten(3,3) = get(S.sl(5),'value');
OUT_data.ten(1,2) = get(S.sl(6),'value');
OUT_data.ten(2,1) = get(S.sl(6),'value');
OUT_data.ten(1,3) = get(S.sl(7),'value');
OUT_data.ten(3,1) = get(S.sl(7),'value');
OUT_data.ten(2,3) = get(S.sl(8),'value');
OUT_data.ten(3,2) = get(S.sl(8),'value');

save OUT_data OUT_data %(structure format)
end %pb_dump_xdata

%**************************************************************************
function [] = pb_dump_ydata(varargin)
%update of Euler angles alpha, beta and gamma is not appropriate since they are different for
%for rotations about xtal-y orientation

global Orienity

S = varargin{3}; % Get structure
delete('OUT_data.mat')

OUT_data.orn = Orienity; %initial x-tal orientation in lab frame
OUT_data.ang(1,1) = str2double(get(S.txtb(4),'str')); %alpha
OUT_data.ang(2,1) = str2double(get(S.txtb(5),'str')); %beta
OUT_data.ang(3,1) = str2double(get(S.txtb(6),'str')); %gamma
OUT_data.ang(4,1) = str2double(get(S.txtb(11),'str')); % = y-offset (THETA about x-axis)

%OUT_data.ang(1,2) = get(S.sl(9),'value'); %  = (new tensor alpha)
%OUT_data.ang(2,2) = get(S.sl(10),'value'); % = (new tensor beta)
%OUT_data.ang(3,2) = get(S.sl(11),'value'); % = (new tensor gamma)
%OUT_data.ang(4,2) = 0; % = z-offset (THETA about x-axis)

OUT_data.ten(1,1) = get(S.sl(3),'value');
OUT_data.ten(2,2) = get(S.sl(4),'value');
OUT_data.ten(3,3) = get(S.sl(5),'value');
OUT_data.ten(1,2) = get(S.sl(6),'value');
OUT_data.ten(2,1) = get(S.sl(6),'value');
OUT_data.ten(1,3) = get(S.sl(7),'value');
OUT_data.ten(3,1) = get(S.sl(7),'value');
OUT_data.ten(2,3) = get(S.sl(8),'value');
OUT_data.ten(3,2) = get(S.sl(8),'value');

save OUT_data OUT_data %(structure format)
end

%**************************************************************************
function [] = pb_dump_zdata(varargin)
%update of Euler angles alpha, beta and gamma is not appropriate since they are different for
%for rotations about xtal-z orientation

global Orienitz

S = varargin{3}; % Get structure
delete('OUT_data.mat')

OUT_data.orn = Orienitz; %initial x-tal orientation in lab frame
OUT_data.ang(1,1) = str2double(get(S.txtb(7),'str')); %alpha
OUT_data.ang(2,1) = str2double(get(S.txtb(8),'str')); %beta
OUT_data.ang(3,1) = str2double(get(S.txtb(9),'str')); %gamma
OUT_data.ang(4,1) = str2double(get(S.txtb(12),'str')); % = z-offset (THETA about x-axis)

%OUT_data.ang(1,2) = get(S.sl(9),'value'); %  = (new tensor alpha)
%OUT_data.ang(2,2) = get(S.sl(10),'value'); % = (new tensor beta)
%OUT_data.ang(3,2) = get(S.sl(11),'value'); % = (new tensor gamma)
%OUT_data.ang(4,2) = 0; % = z-offset (THETA about z-axis)

OUT_data.ten(1,1) = get(S.sl(3),'value');
OUT_data.ten(2,2) = get(S.sl(4),'value');
OUT_data.ten(3,3) = get(S.sl(5),'value');
OUT_data.ten(1,2) = get(S.sl(6),'value');
OUT_data.ten(2,1) = get(S.sl(6),'value');
OUT_data.ten(1,3) = get(S.sl(7),'value');
OUT_data.ten(3,1) = get(S.sl(7),'value');
OUT_data.ten(2,3) = get(S.sl(8),'value');
OUT_data.ten(3,2) = get(S.sl(8),'value');

save OUT_data OUT_data %(structure format)
end

%**************************************************************************
function [] = pb_dump_alldata(varargin)
%Alldataout from XTALFIT session. Send-to-file the tensor, all orientations Euler angles, and offsets.
global Orienitx Orienity Orienitz

S = varargin{3}; % Get structure
%delete Alldata.mat

Alldata.ten=[1 0 0;0 1 0;0 0 1]; %preallocate (is this necessary?)
Alldata.ang = zeros(4,3); %xang = Alldata.ang(:,1), yang = Alldata.ang(:,2), zang = Alldata.ang(:,3)
Alldata.orn = zeros(3,9); %xorien = Alldata.orn(:,1:3), yorien = Alldata.orn(:,4:6), zorien = Alldata.orn(:,7:9),

Alldata.ten(1,1) = get(S.sl(3),'value');
Alldata.ten(2,2) = get(S.sl(4),'value');
Alldata.ten(3,3) = get(S.sl(5),'value');
Alldata.ten(1,2) = get(S.sl(6),'value');
Alldata.ten(2,1) = get(S.sl(6),'value');
Alldata.ten(1,3) = get(S.sl(7),'value');
Alldata.ten(3,1) = get(S.sl(7),'value');
Alldata.ten(2,3) = get(S.sl(8),'value');
Alldata.ten(3,2) = get(S.sl(8),'value');

Alldata.orn = [Orienitx Orienity Orienitz]; %initial x-tal orientation in lab frame

Alldata.ang(1,1) = str2double(get(S.txtb(1),'str')); %alpha
Alldata.ang(2,1) = str2double(get(S.txtb(2),'str')); %beta
Alldata.ang(3,1) = str2double(get(S.txtb(3),'str')); %gamma
Alldata.ang(4,1) = str2double(get(S.txtb(10),'str')); % = x-offset (THETA about x-axis)

Alldata.ang(1,2) = str2double(get(S.txtb(4),'str')); %alpha
Alldata.ang(2,2) = str2double(get(S.txtb(5),'str')); %beta
Alldata.ang(3,2) = str2double(get(S.txtb(6),'str')); %gamma
Alldata.ang(4,2) = str2double(get(S.txtb(11),'str')); % = y-offset (THETA about y-axis)

Alldata.ang(1,3) = str2double(get(S.txtb(7),'str')); %alpha
Alldata.ang(2,3) = str2double(get(S.txtb(8),'str')); %beta
Alldata.ang(3,3) = str2double(get(S.txtb(9),'str')); %gamma
Alldata.ang(4,3) = str2double(get(S.txtb(12),'str')); % = z-offset (THETA about z-axis)

Alldata.ang(1,4) = get(S.sl(9),'value'); %  = alpha (new tensor alpha)
Alldata.ang(2,4) = get(S.sl(10),'value'); % = beta (new tensor beta)
Alldata.ang(3,4) = get(S.sl(11),'value'); % = gamma (new tensor gamma)
Alldata.ang(4,4) = 0; % = z-offset (THETA about z-axis)

save Alldata Alldata %(structure format, you need this syntax with save or else you get an error message)
end

