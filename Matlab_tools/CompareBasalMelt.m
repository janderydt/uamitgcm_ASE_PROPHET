function CompareBasalMelt

addpath(getenv("froot_tools"));

%% compare melt between 2004 and 2014 for 
%% 1. AS_PROPHET standalone paris2C with default melt
%% 2. AS_PROPHET coupled paris2C with default melt
%% 3. AS_PROPHET coupled paris2C with Dan's modification to the melt (5 cases)
%% 4. PTDC coupled Kimura with default melt

section = ["AS","PIG"];

froot = [getenv("froot_uamitgcm")+"cases/AS_PROPHET_005/output/",...
    getenv("froot_uamitgcm")+"cases/AS_PROPHET_002/output/",...
    getenv("froot_uamitgcm")+"cases/AS_PROPHET_008/output/",...
    getenv("froot_uamitgcm")+"cases/ASE_varmelt/output/"];

titlelabels = {{'standalone paris2c';'default melt'},{'coupled paris2c';'default melt'},...
    {'coupled paris2c';'modified melt'},...
    {'coupled Kimura';'default melt'}};

legendlabels = ["standalone paris2c, default melt","coupled paris2c, default melt",...
    "coupled paris2c, modified melt",...
    "coupled Kimura, default melt"];

for ii=1:numel(froot)

    subd=dir(froot(ii));
    isub = [subd(:).isdir]; %# returns logical vector
    nameFolds = {subd(isub).name}';
    nameFolds(ismember(nameFolds,[".",".."])) = [];

    i1=1; i2=1; i3=1; i4=1; i5=1; starttime=0; nn=0; data(ii).melt=[];
    starttime = 0;
    
    while starttime<=datenum("20150101","yyyymmdd") && nn<numel(nameFolds)

        nn=nn+1;

        MITpath = froot(ii)+nameFolds{nn}+"/MITgcm";
        MITfile = MITpath+"/output.nc";
        nstr = strlength(MITpath);
        startdate = extractBetween(MITpath,nstr-12,nstr-7);    
        starttime = datenum(startdate+"01","yyyymmdd");
        T = double(ncread(MITfile,"time"));
        % read time epoch
        attvalue=ncreadatt(MITfile,"time","units");
        if strfind(attvalue,"seconds")
            epoch = erase(attvalue,"seconds since ");
            epochnum = datenum(epoch);
            MITTime = epochnum + T/(24*60*60);
        elseif strfind(attvalue,"days")
            epoch = erase(attvalue,"days since ");
            epochnum = datenum(epoch);
            MITTime = epochnum + T;
        else
            error("I do not recognise the time format in output.nc");
        end
        
        if any(MITTime>datenum("01012004","ddmmyyyy"))       
        %if starttime>datenum("31122096","ddmmyyyy") 

            Tind = find(MITTime>datenum("01012004","ddmmyyyy"));

            for tt = 1:numel(Tind)

                   [LON,LAT,Melt,~,~,~,~,~,~,Draft] = PlotMeltRates(MITpath,Tind(tt),[1 0 0 0 0 0 0 1]);
                   dx = LON(2,1) - LON(1,1); dy = LAT(1,2) - LAT(1,1);
                   [data(ii).nx,data(ii).ny] = size(Melt);
                   
                   data(ii).melt(i1,:) = Melt(:)*365.25*24*60*60/1e3; %kg/s/m2 to m/yr
                   
                   Melt = Melt*365.25*24*60*60*dx*dy/1e12; %kg/s/m2 to Gt/yr
                   Imelt = find(Melt~=0);

                   for ss=1:numel(section)
                        switch section{ss}
                            case "AS"
                               data(ii).section(ss).melt_int(i1) = sum(Melt(Imelt),"all","omitnan");
                               data(ii).section(ss).melt_mean(i1) = mean(data(ii).melt(i1,Imelt),2);
                               data(ii).section(ss).time(i1) = MITTime(Tind(tt));
                               data(ii).section(ss).bins(i1,:) = [-1800:20:0];
                               Y = discretize(Draft(Imelt),data(ii).section(ss).bins(i1,:));
                               for bb=1:size(data(ii).section(ss).bins(i1,:),2)-1
                                   Itmp = find(Y==bb);
                                   data(ii).section(ss).melt_mean_per_bin(i1,bb) = mean(data(ii).melt(i1,Imelt(Itmp)),"all","omitnan");
                               end
                               i1 = i1+1;
                            case "PIG"
                                xmin = -1699e3; xmax = -1530e3; ymin = -380e3; ymax = -220e3;
                                IPIG = find(LON(Imelt)>xmin & LON(Imelt)<xmax & LAT(Imelt)>ymin & LAT(Imelt)<ymax);
                                data(ii).section(ss).melt_int(i2) = sum(Melt(Imelt(IPIG)),"all","omitnan");
                                data(ii).section(ss).melt_mean(i2) = mean(data(ii).melt(i2,Imelt(IPIG)),2);
                                data(ii).section(ss).time(i2) = MITTime(Tind(tt));
                                data(ii).section(ss).bins(i2,:) = [-1800:20:0];
                                Y = discretize(Draft(Imelt(IPIG)),data(ii).section(ss).bins(i2,:));
                                for bb=1:size(data(ii).section(ss).bins(i2,:),2)-1
                                    Itmp = find(Y==bb);
                                    data(ii).section(ss).melt_mean_per_bin(i2,bb) = mean(data(ii).melt(i2,Imelt(IPIG(Itmp))),"all","omitnan");
                                end
                                i2 = i2+1;
                            case "TW"
                                xmin = -1620e3; xmax = -1500e3; ymin = -520e3; ymax = -380e3;
                                ITW = find(LON(Imelt)>xmin & LON(Imelt)<xmax & LAT(Imelt)>ymin & LAT(Imelt)<ymax);
                                
                            case "CR"
                                xpoly = [-1610 -1485 -1450 -1450 -1610]*1e3;
                                ypoly = [-580 -657 -657 -520 -520]*1e3; 
                                ICR = find(inpoly([LON(Imelt(:)),LAT(Imelt(:))],[xpoly(:) ypoly(:)]));
                                
                            case "DT"
                                xpoly = [-1610 -1485 -1450 -1450 -1625]*1e3;
                                ypoly = [-580 -657 -657 -700 -700]*1e3; 
                                IDT = find(inpoly([LON(Imelt(:)),LAT(Imelt(:))],[xpoly(:) ypoly(:)]));
                                
                        end
                   end
   
            end
            
            fprintf("%s: Done %i \n",froot(ii),nn);
        end
        
    end

    for nn=1:size(data(ii).melt,2)

        data(ii).meanmelt(nn) = mean(data(ii).melt(:,nn),"all","omitnan");

    end

    data(ii).meanmelt = reshape(data(ii).meanmelt,data(ii).nx,data(ii).ny);

end

save dump.mat

CM = [0.4941    0.1843    0.5569;...
    0.5496    0.1936    0.5173;...
    0.6051    0.2029    0.4777;...
    0.6606    0.2123    0.4381;...
    0.7161    0.2216    0.3985;...
    0.7618    0.2315    0.3680;...
    0.7919    0.2424    0.3522;...
    0.8201    0.2534    0.3381;...
    0.8120    0.2503    0.3421;...
    0.8344    0.2684    0.3393;...
    0.9414    0.3656    0.3319;...
    0.9710    0.4627    0.3735;...
    0.9867    0.5547    0.4087;...
    0.9946    0.6400    0.4428;...
    0.9983    0.7197    0.4835;...
    0.9996    0.7900    0.5290;...
    0.9999    0.8634    0.5882;...
    0.9989    0.9432    0.6700;...
    0.9925    0.9683    0.7240;...
    0.9797    0.9854    0.7838;...
    0.9590    0.9954    0.8464;...
    0.9267    0.9982    0.9073;...
    0.8848    0.9970    0.9506;...
    0.8119    0.9864    0.9796;...
    0.7008    0.9600    0.9952;...
    0.5297    0.8851    0.9990;...
    0.3731    0.7412    1.0000;...
    0.2433    0.5450    1.0000];
CM = [interp1(1:size(CM,1),CM(:,1),linspace(1,size(CM,1),21))',...
    interp1(1:size(CM,1),CM(:,2),linspace(1,size(CM,1),21))',...
    interp1(1:size(CM,1),CM(:,3),linspace(1,size(CM,1),21))'];
CM = flipdim(CM,1);

%% pattern
H=fig("units","inches","width",90*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(1,numel(data),"TileSpacing","compact"); n=1;
for i = 1:numel(data)
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for ii=1:numel(data)
    M = data(ii).meanmelt;
    M(M==0)=NaN;
    M(M>0)=0;
    M(M<-200)=-200;
    contourf(ax_fig(ii),LON/1e3,LAT/1e3,-M,[0:5:200],"LineStyle","none"); 
    shading(ax_fig(ii),"flat");
    colormap(ax_fig(ii),CM);
    axis(ax_fig(ii),"equal");
    grid(ax_fig(ii),"on"); 
    box(ax_fig(ii),"on");
    title(ax_fig(ii),titlelabels{ii});
    if ii>1
        yticklabels(ax_fig(ii),'');
    end
end
xlabel(tlo_fig,"psx [km]");
ylabel(tlo_fig,"psy [km]");
cb=colorbar(ax_fig(numel(data)),"Location","eastoutside","Ticks",[0:50:200]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.025;
cb.FontSize=16;
cb.Label.String = "$\mbox{Average basal melt 2004-2014 [m yr}^{-1}\mbox{]}$";
cb.Label.Interpreter = "latex";

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/BasalMeltPatterns_2004-2014";
print(H,fname,'-dpng','-r400');


%% total
H=fig("units","inches","width",90*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");
subplot('position',[0.1 0.1 0.85 0.85]); hold on;

for ii=1:numel(data)
    g(ii)=plot(data(ii).section(1).time,-data(ii).section(1).melt_int,'LineWidth',1);
end
grid on; box on;
xlim([datenum('01012004','ddmmyyyy') datenum('31122014','ddmmyyyy')]);
ylabel('Basal melt [Gt/yr]');
datetick('x','keeplimits');
legend(g(:),legendlabels,'location','northeast');

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/IntegratedMelt_2004-2014";
print(H,fname,'-dpng','-r400');

%% total PIG
H=fig("units","inches","width",90*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");
subplot('position',[0.1 0.1 0.85 0.85]); hold on;

for ii=1:numel(data)
    g(ii)=plot(data(ii).section(2).time,-data(ii).section(2).melt_int,'LineWidth',1);
end

% PIG data:
% Dutrieux 2014 + Heywood 2016 estimates of PIG melting (Gt/y)
% Density of freshwater (kg/m^3)
rho_fw = 1e3;
% Density of ice (kg/m^3)
rho_ice = 917;
dutrieux_sw = [51.3, 79.7, 75.2, 37.3]*rho_ice/rho_fw;
dutrieux_mw = [49.1, 79.4, 69.2, 34.7]*rho_ice/rho_fw;
dutrieux_melt = 0.5*(dutrieux_sw+dutrieux_mw);
dutrieux_err_fac = 0.1;
heywood_melt = [40]*rho_ice/rho_fw;
heywood_err_fac = 0.4;
pig_melt_years = [1994, 2009, 2010, 2012, 2014];
pig_melt_years = datenum(num2str(pig_melt_years'),'yyyy');
pig_melt = [dutrieux_melt, heywood_melt];
pig_err = [dutrieux_melt.*dutrieux_err_fac, heywood_melt.*heywood_err_fac];

g(ii+1)=errorbar(pig_melt_years,pig_melt,pig_err,'.r','LineWidth',2,'markersize',0.1);% yobs,fluxobs,"dk","markersize",10,"markerfacecolor","m"

grid on; box on;
xlim([datenum('01012004','ddmmyyyy') datenum('31122014','ddmmyyyy')]);
ylabel('Basal melt Pine Island [Gt/yr]');
datetick('x','keeplimits');
legend(g(:),[legendlabels,"Observations"],'location','northeast');

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/IntegratedMelt_PIG_2004-2014";
print(H,fname,'-dpng','-r400');


%% melt-vs-depth for PIG
H=fig("units","inches","width",60*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");
%subplot('position',[0.1 0.1 0.85 0.85]); hold on;

tlo_fig = tiledlayout(1,numel(data),"TileSpacing","compact"); n=1;
for i = 1:numel(data)
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

colors = ["Greys","Greys","Greys","Greys"];%"Oranges","Greens","Greens","Greens","Greens","Greens","Purples"];

load("cryos_melt_draft.mat");
bins = [-1800:20:0]; Itmp = find(ml~=0 & ~isnan(ml));
Y = discretize(draft(Itmp),bins);
for bb=1:numel(bins)-1
   Ibb = find(Y==bb);
   MeanMelt_Noel(bb) = mean(ml(Itmp(Ibb)),"all","omitnan");
end

for ii=1:numel(data)
    nt = size(data(ii).section(2).bins,1);
    CM = brewermap(2*nt,colors(ii));
    plot(ax_fig(ii),-ml,draft,'.','color',[0.5 0.5 1],'markersize',6);
    for tt=1:nt
        plot(ax_fig(ii),data(ii).section(2).melt_mean_per_bin(tt,:),(data(ii).section(2).bins(tt,1:end-1)+data(ii).section(2).bins(tt,2:end))/2,'color',CM(floor(nt/2)+tt,:));
    end
    g(ii)=plot(ax_fig(ii),mean(data(ii).section(2).melt_mean_per_bin,1,"omitnan"),(data(ii).section(2).bins(tt,1:end-1)+data(ii).section(2).bins(tt,2:end))/2,'-k','linewidth',2);
    
    h=plot(ax_fig(ii),-MeanMelt_Noel,(bins(2:end)+bins(1:end-1))/2,'--b','linewidth',2);
    grid(ax_fig(ii),"on"); 
    box(ax_fig(ii),"on");
    title(ax_fig(ii),titlelabels{ii},"fontsize",8);
    if ii>1
        yticklabels(ax_fig(ii),'');
    end
    ylim(ax_fig(ii),[-1000 0]);
    xlim(ax_fig(ii),[-150 0]);

end
legend([g(1) h],["Ua-MITgcm","Cryosat2 Noel G."],'location','northwest')

xlabel(tlo_fig,"Melt [m/yr]");
ylabel(tlo_fig,"Depth [m]");

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/MeltvsDepth_2004-2014";
print(H,fname,'-dpng','-r400');



