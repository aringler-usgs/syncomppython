clear all;
clc;
warning off;

datadire='deep';

system(['ls ' datadire '/*LXZ* > tempfile']);
[sys,len]=system('wc -l tempfile');
len=str2double(strrep(len,'tempfile',''));
fid=fopen('tempfile');
system(['mkdir ' datadire 'plots']);
for index=1:len
    display([num2str(index) ' of ' num2str(len)]);
    try
        tempseis=fgetl(fid);
        sta=strrep(tempseis,[datadire '/'],'');
        sta = strrep(sta,'.XX.LXZ.modes.sac.proc','');
        net = ls([datadire '/' sta '*SAC.proc']);
        net = net(length([datadire '/' sta '.'])+1:length([datadire '/' ...
            sta '.'])+2);
        synZ=[datadire '/'  sta '.XX.LXZ.modes.sac.proc'];
        data00Z=[datadire '/'  sta '.' net '.00.LHZ.SAC.proc'];
        data10Z=[datadire '/' sta '.' net '.10.LHZ.SAC.proc'];
        synN=[datadire '/'  sta '.XX.LXN.modes.sac.proc'];
        data00N=[datadire '/'  sta '.' net '.00.LHN.SAC.proc'];
        data10N=[datadire '/'  sta '.' net '.10.LHN.SAC.proc'];
        synE=[datadire '/'  sta '.XX.LXE.modes.sac.proc'];
        data00E=[datadire '/'  sta '.' net '.00.LHE.SAC.proc'];
        data10E=[datadire '/' sta '.' net '.10.LHE.SAC.proc'];

        stsZ=read_sac(synZ,1);
        stsN=read_sac(synN,1);
        stsE=read_sac(synE,1);
        azi=stsZ(12,3);
        jday=num2str(stsZ(15,3));
        dist=stsZ(11,3)/(2*pi*6378/360);
        hour = num2str(stsZ(16,3));
        min = num2str(stsZ(17,3));
        sec = num2str(stsZ(18,3));
        year = num2str(stsZ(14,3));

        dts00Z=read_sac(data00Z,1);
        dts00N=read_sac(data00N,1);
        dts00E=read_sac(data00E,1);

        try
            dts10Z=read_sac(data10Z,1);
            dts10N=read_sac(data10N,1);
            dts10E=read_sac(data10E,1);
            dts10='true';
        catch
            dts10='false';
        end

        if(isempty(dts10E))
            dts10='false';
        end

        [b,a]=butter(2,1/500,'low');
        [d,c]=butter(2,1/100,'high');

        dts00Z(:,2)=filtfilt(b,a,dts00Z(:,2));
        dts00N(:,2)=filtfilt(b,a,dts00N(:,2));
        dts00E(:,2)=filtfilt(b,a,dts00E(:,2));
        dts00Z(:,2)=filtfilt(d,c,dts00Z(:,2));
        dts00N(:,2)=filtfilt(d,c,dts00N(:,2));
        dts00E(:,2)=filtfilt(d,c,dts00E(:,2));

        stsZ(:,2)=filtfilt(b,a,stsZ(:,2));
        stsN(:,2)=filtfilt(b,a,stsN(:,2));
        stsE(:,2)=filtfilt(b,a,stsE(:,2));
        stsZ(:,2)=filtfilt(d,c,stsZ(:,2));
        stsN(:,2)=filtfilt(d,c,stsN(:,2));
        stsE(:,2)=filtfilt(d,c,stsE(:,2));

        if(strcmp(dts10,'true'))
            dts10Z(:,2)=filtfilt(b,a,dts10Z(:,2));
            dts10N(:,2)=filtfilt(b,a,dts10N(:,2));
            dts10E(:,2)=filtfilt(b,a,dts10E(:,2));
            dts10Z(:,2)=filtfilt(d,c,dts10Z(:,2));
            dts10N(:,2)=filtfilt(d,c,dts10N(:,2));
            dts10E(:,2)=filtfilt(d,c,dts10E(:,2));
        end


        fig1=figure('visible','off');
        clf
        subplot(3,1,1);
        p1=plot(stsZ(:,1),stsZ(:,2),'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        p2=plot(dts00Z(:,1),dts00Z(:,2),'LineWidth',1,'Color','g');
        if(strcmp(dts10,'true'))
            p3=plot(dts10Z(:,1),dts10Z(:,2),'LineWidth',1,'Color','r');
            legend([p1 p2 p3],'Synthetic','00 ','10 ','FontSize',5);
        else
            legend([p1 p2],'Synthetic','00 ','FontSize',5);
        end
        xlabel('Time (s)','FontSize',5);
        ylabel('Displacement (m)','FontSize',5);
        xlim([0 4000]);
        title([net ' ' sta ' Z ' year ' ' jday ' ' hour ':' min ':' sec ... 
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);
        subplot(3,1,2);
        p1=plot(stsN(:,1),stsN(:,2),'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        p2=plot(dts00N(:,1),dts00N(:,2),'LineWidth',1,'Color','g');
        if(strcmp(dts10,'true'))
            p3=plot(dts10N(:,1),dts10N(:,2),'LineWidth',1,'Color','r');
            legend([p1 p2 p3],'Synthetic','00 ','10 ','FontSize',5);
        else
            legend([p1 p2],'Synthetic','00 ','FontSize',5);
        end
        xlabel('Time (s)','FontSize',5);
        ylabel('Displacement (m)','FontSize',5);
        xlim([0 4000]);
        title([net ' ' sta ' N ' year ' ' jday ' ' hour ':' min ':' sec ...
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);

        subplot(3,1,3);
        xlim([0 4000])
        p1=plot(stsE(:,1),stsE(:,2),'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        p2=plot(dts00E(:,1),dts00E(:,2),'LineWidth',1,'Color','g');
        if(strcmp(dts10,'true'))
            p3=plot(dts10E(:,1),dts10E(:,2),'LineWidth',1,'Color','r');
            legend([p1 p2 p3],'Synthetic','00 ','10 ','FontSize',5);
        else
            legend([p1 p2],'Synthetic','00 ','FontSize',5);
        end
        xlabel('Time (s)','FontSize',5);
        ylabel('Displacement (m)','FontSize',5);
        xlim([0 4000]);
        title([net ' ' sta ' E ' year ' ' jday ' ' hour ':' min ':' sec ...
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);
        %orient Landscape
        print('-djpeg',[datadire 'plots/' net '.' sta '.' year '.' ...
            jday '.' hour 'Tseries.jpg']);


        [ms00Z,fre]=mscohere(stsZ(:,2),dts00Z(:,2),1000,500,500,1);
        [ms00N,fre]=mscohere(stsN(:,2),dts00N(:,2),1000,500,500,1);
        [ms00E,fre]=mscohere(stsE(:,2),dts00E(:,2),1000,500,500,1);

        if(strcmp(dts10,'true'))
            [ms10Z,fre]=mscohere(stsZ(:,2),dts10Z(:,2),1000,500,500,1);
            [ms10N,fre]=mscohere(stsN(:,2),dts10N(:,2),1000,500,500,1);
            [ms10E,fre]=mscohere(stsE(:,2),dts10E(:,2),1000,500,500,1);
        end

        fig2=figure('visible','off');
        clf
        subplot(3,1,1);
        p1=plot(1./fre,ms00Z,'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        if(strcmp(dts10,'true'))
            p2=plot(1./fre,ms10Z,'LineWidth',1,'Color','r');
            legend([p1 p2],'00 ','10 ','FontSize',5);
        else
            legend(p1,'00 ','FontSize',5);
        end
        xlim([0 250]);
        ylim([.4 1]);
        xlabel('Period (seconds)','FontSize',5);
        ylabel('Coherence (\gamma^2)','FontSize',5);
        title([net ' ' sta ' Z ' year ' ' jday ' ' hour ':' min ':' sec ...
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);
        subplot(3,1,2);
        p1=plot(1./fre,ms00N,'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        if(strcmp(dts10,'true'))
            p2=plot(1./fre,ms10N,'LineWidth',1,'Color','r');
            legend([p1 p2],'00 ','10 ','FontSize',5);
        else
            legend(p1,'00 ','FontSize',5);
        end
        xlabel('Period (seconds)','FontSize',5);
        ylabel('Coherence (\gamma^2)','FontSize',5);
        title([net ' ' sta ' N ' year ' ' jday ' ' hour ':' min ':' sec ...
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);
        xlim([0 250]);
        ylim([.4 1]);
        subplot(3,1,3);
        p1=plot(1./fre,ms00E,'LineWidth',1,'Color','k');
        set(gca,'FontSize',5);
        hold on 
        if(strcmp(dts10,'true'))
            p2=plot(1./fre,ms10E,'LineWidth',1,'Color','r');
            legend([p1 p2],'00 ','10 ','FontSize',5);
        else
            legend(p1,'00 ','FontSize',5);
        end
        xlabel('Period (seconds)','FontSize',5);
        ylabel('Coherence (\gamma^2)','FontSize',5);
        title([net ' ' sta ' E ' year ' ' jday ' ' hour ':' min ':' sec ...
            ' Distance: ' num2str(dist) ' Azimuth: ' num2str(azi)],'FontSize',5);
        xlim([0 250]);
        ylim([.4 1]);
        print('-djpeg',[datadire 'plots/' net '.' sta '.' year '.' ...
            jday '.' hour 'Cohere.jpg']);
    catch  
    end
end
system('rm tempfile');
