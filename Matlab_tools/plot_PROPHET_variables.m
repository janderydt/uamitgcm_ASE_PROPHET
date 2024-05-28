function plot_PROPHET_variables

groupID = ["AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET","AS_PROPHET",...
    "PTDC","PTDC","PTDC"];
runID = ["004","005","006","001","002","003","007","008","009",...
    "001","003","002"];
runLabel = ["Static cavity REF","Static cavity Paris2C","Static cavity RCP8.5","REF - default melt","Paris2C - default melt","RCP8.5 - default melt",...
    "REF - modified melt","Paris2C - modified melt","RCP8.5 - modified melt",...
    "Coupled varmelt","Coupled avmelt","Coupled himelt"];

runLineStyle = ["-","-","-","-","-","-","-","-","-",...
    "-","-","-"];
runLineWidth = [1 1 1 1 1 1 1 1 1 1 1 1];
CM1 = brewermap(6,"Greys"); % for stand alone ocean
CM2 = brewermap(6,"Oranges"); % for the default melt parameterization
CM3 = brewermap(6,"Purples"); % for Dan's modified melt parameteriztion
CM4 = brewermap(7,"Greens"); % for De Rydt and Naughten 2023
runLineColor = [CM1([3 4 6],:); CM2([2:2:6],:); CM3([3 4 6],:); CM4([3 5 7],:)];
couplingtimestep = [48 48 48 1 1 1 1 1 1 1 1 1];

basins = [2 3]; nb = numel(basins);
runs_to_plot = [1:9];%[1:9];

groupID = groupID(runs_to_plot);
runID = runID(runs_to_plot);
runLabel = runLabel(runs_to_plot);
runLineStyle = runLineStyle(runs_to_plot);
runLineWidth = runLineWidth(runs_to_plot);
runLineColor = runLineColor(runs_to_plot,:);
couplingtimestep = couplingtimestep(runs_to_plot);

froot_tools = getenv("froot_tools");
froot_uamitgcm = getenv("froot_uamitgcm");

% dVAF
H1=fig('units','inches','width',nb*40*12/72.27,'height',40*12/72.27,'fontsize',16,'font','Helvetica');
tlo1 = tiledlayout(H1,1,numel(basins),"TileSpacing","compact");
for tt=1:numel(basins)
    h1(tt) = nexttile(tlo1); hold(h1(tt),"on");
end

% Basal Melt
H2=fig('units','inches','width',nb*40*12/72.27,'height',40*12/72.27,'fontsize',16,'font','Helvetica');
tlo2 = tiledlayout(H2,1,numel(basins),"TileSpacing","compact");
for tt=1:numel(basins)
    h2(tt) = nexttile(tlo2); hold(h2(tt),"on");
end

% dA
H3=fig('units','inches','width',nb*40*12/72.27,'height',40*12/72.27,'fontsize',16,'font','Helvetica');
tlo3 = tiledlayout(H3,1,numel(basins),"TileSpacing","compact");
for tt=1:numel(basins)
    h3(tt) = nexttile(tlo3); hold(h3(tt),"on");
end

% dVAFdt
H4=fig('units','inches','width',nb*40*12/72.27,'height',40*12/72.27,'fontsize',16,'font','Helvetica');
tlo4 = tiledlayout(H4,1,numel(basins),"TileSpacing","compact");
for tt=1:numel(basins)
    h4(tt) = nexttile(tlo4); hold(h4(tt),"on");
end

for ii = 1:numel(runID)

    Exp = groupID(ii)+"_"+runID(ii);

    if strfind(Exp,"PROPHET")

        data = load_PROPHET_data(froot_tools,froot_uamitgcm,Exp,couplingtimestep(ii));

    elseif strfind(Exp,"PTDC")
        
        addpath(getenv("froot_tools"));
        data = VAF_GroundedArea_GLFlux_Coupled(froot_tools,froot_uamitgcm,Exp,couplingtimestep(ii),0);

    end

    for bb=1:numel(basins)

        % basin panels for BMB
        g2(ii) = plot(h2(bb),data.Basins(basins(bb)).Time,data.Basins(basins(bb)).BasalMassBalanceMITgcm/1e12*365.25*24*60*60,'LineStyle',runLineStyle(ii),'color',runLineColor(ii,:),'LineWidth',runLineWidth(ii)); %kg/s to Gt/yr
        grid(h2(bb),"on"); box(h2(bb),"on"); title(h2(bb),data.Basins(basins(bb)).Name);
        xlim(h2(bb),[datenum('01012004','ddmmyyyy'),datenum('01012100','ddmmyyyy')]);
        datetick(h2(bb),"x","yyyy","keeplimits");

        nonnan_ind = find(~isnan(data.Basins(basins(bb)).VAF));
        if ~isempty(nonnan_ind)
            % basin panels for change in VAF                       
            %plot(h1(bb),data.Basins(basins(bb)).Time,(data.Basins(basins(bb)).VAF-data.Basins(basins(bb)).VAF(nonnan_ind(1)))/3.625e14*1e2,'LineStyle',runLineStyle(ii),'color',runLineColor(ii,:),'LineWidth',runLineWidth(ii)); % change in cm SLE
            yyaxis(h1(bb),"left");
            g1(ii) = plot(h1(bb),data.Basins(basins(bb)).Time,(data.Basins(basins(bb)).VAF-data.Basins(basins(bb)).VAF(nonnan_ind(1)))/1e9,'LineStyle',runLineStyle(ii),'color',runLineColor(ii,:),'LineWidth',runLineWidth(ii),'Marker','none'); % from m3 to Gt
            grid(h1(bb),"on"); box(h1(bb),"on"); title(h1(bb),data.Basins(basins(bb)).Name);
            xlim(h1(bb),[datenum('01012004','ddmmyyyy'),datenum('01012100','ddmmyyyy')]);
            datetick(h1(bb),"x","yyyy","keeplimits");
            ylimtmp = ylim(h1(bb));
            yyaxis(h1(bb),"right");
            ylim(h1(bb),ylimtmp*1e9/3.625e14*1e2);

            % basin panels for change in grounded area
            g3(ii) = plot(h3(bb),data.Basins(basins(bb)).Time,(data.Basins(basins(bb)).GroundedArea-data.Basins(basins(bb)).GroundedArea(nonnan_ind(1)))/1e6,'LineStyle',runLineStyle(ii),'color',runLineColor(ii,:),'LineWidth',runLineWidth(ii)); %in km2
            grid(h3(bb),"on"); box(h3(bb),"on"); title(h3(bb),data.Basins(basins(bb)).Name);
            xlim(h3(bb),[datenum('01012004','ddmmyyyy'),datenum('01012100','ddmmyyyy')]);
            datetick(h3(bb),"x","yyyy","keeplimits");

            % basin panels for VAF rate
            dVAF = (data.Basins(basins(bb)).VAF(2:end)-data.Basins(basins(bb)).VAF(1:end-1))/1e9;
            dt = (data.Basins(basins(bb)).Time(2:end)-data.Basins(basins(bb)).Time(1:end-1))/365.25; % in years
            yyaxis(h4(bb),"left");
            g4(ii) = plot(h4(bb),data.Basins(basins(bb)).Time(1:end-1),dVAF./dt,'LineStyle',runLineStyle(ii),'color',runLineColor(ii,:),'LineWidth',runLineWidth(ii),'Marker','none');
            grid(h4(bb),"on"); box(h4(bb),"on"); title(h4(bb),data.Basins(basins(bb)).Name);
            xlim(h4(bb),[datenum('01012004','ddmmyyyy'),datenum('01012100','ddmmyyyy')]);
            datetick(h4(bb),"x","yyyy","keeplimits");  
            ylimtmp = ylim(h4(bb));
            yyaxis(h4(bb),"right");
            ylim(h4(bb),ylimtmp*1e9/3.625e14*1e3); % in mm SLE/yr

        else
            g1(ii)= plot(-9999,-9999);
            g3(ii)= plot(-9999,-9999);
            g4(ii)= plot(-9999,-9999);
        end

    end
    
end

ind_dynamic = find(~contains(runLabel,"Static"));

for bb=1:numel(basins)
    h1(bb).YAxis(1).Color = 'k';
    h1(bb).YAxis(2).Color = 'k';
end
yyaxis(h1(numel(basins)),"right");
ylabel(h1(numel(basins)),"Change in VAF [cm SLE]");
yyaxis(h1(1),"left");
legend(h1(1),g1(ind_dynamic),runLabel(ind_dynamic),'Location','southwest','fontsize',12);
ylabel(h1(1),"Change in VAF [Gt]");
xlabel(tlo1,"Time");
pos = get(H1,'Position');
set(H1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/dVAF";
for ii=1:numel(runID)
    fname = fname+"_"+runID(ii);
end
print(H1,fname,'-dpng','-r400');

legend(h2(1),g2(:),runLabel,'Location','northeast','fontsize',12);
ylabel(h2(1),"Basal melt [Gt/yr]");
xlabel(tlo2,"Time");
pos = get(H2,'Position');
set(H2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/BasalMelt";
for ii=1:numel(runID)
    fname = fname+"_"+runID(ii);
end
print(H2,fname,'-dpng','-r400');

legend(h3(1),g3(ind_dynamic),runLabel(ind_dynamic),'Location','northeast','fontsize',12);
ylabel(h3(1),"Change in grounded area [km^2]");
xlabel(tlo3,"Time");
pos = get(H3,'Position');
set(H3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/dA";
for ii=1:numel(runID)
    fname = fname+"_"+runID(ii);
end
print(H3,fname,'-dpng','-r400');

for bb=1:numel(basins)
    h4(bb).YAxis(1).Color = 'k';
    h4(bb).YAxis(2).Color = 'k';
end
yyaxis(h4(numel(basins)),"right");
ylabel(h4(numel(basins)),"VAF rate [mm SLE/yr]");
yyaxis(h4(1),"left");
legend(h4(2),g4(ind_dynamic),runLabel(ind_dynamic),'Location','northeast','fontsize',12);
ylabel(h4(1),"VAF rates [Gt/yr]");
xlabel(tlo4,"Time");
pos = get(H4,'Position');
set(H4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/dVAFdt";
for ii=1:numel(runID)
    fname = fname+"_"+runID(ii);
end
print(H4,fname,'-dpng','-r400');






