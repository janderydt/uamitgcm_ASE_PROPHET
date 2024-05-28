function plot_PROPHET_changeinVAF

plotoption='icethickness_speed';
runID='AS_PROPHET_009';
Experiment="RCP85 - modified melt";
basin = 'AS';

switch basin
    case 'PIG'
        %% PIG
        xmin = -1700; xmax = -1530; 
        ymin = -370; ymax = -210;
    case 'TW'
        %% Thwaites
        xmin = -1620; xmax = -1450;
        ymin = -520; ymax = -380;
    case 'DC'
        %% Dotson Crosson
        xmin = -1620; xmax = -1440;
        ymin = -700; ymax = -520;
    case 'AS'
        %% Amundsen
        xmin = -1700; xmax = -1350;
        ymin = -800; ymax = 0;
    otherwise
        error('case unknown.')
end

froot_tools = getenv("froot_tools");
froot_uamitgcm = getenv("froot_uamitgcm");

frootm = [froot_uamitgcm,'/cases/',runID];
%frootm = ['/media/janryd69/mainJDeRydt/UaMITgcm_v2/cases/',runID];
timestep1 = 2; % after 4-year relaxation
timestep2 = 96*12+2;

if timestep1>=timestep2
   error('check timesteps'); 
end

subd=dir([frootm,'/output/']);
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = {subd(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% restartfile = dir([frootm,'/ua_custom/*-RestartFile.mat']);
% files{1}.folder = restartfile(1).folder;
% files{1}.name = restartfile(1).name;

nfolders = length(nameFolds); kk=2;
for jj=1:nfolders
    %disp([frootm,'/output/',nameFolds{jj}]);
    output(jj).filelist=dir([frootm,'/output/',nameFolds{jj},'/Ua/UaDefaultRun_',char(nameFolds{jj}),'*.mat']);
    for ii=1:length(output(jj).filelist)
        files{kk}.folder = output(jj).filelist(ii).folder;
        files{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
end

[files{timestep1}.folder,'/',files{timestep1}.name]
load([files{timestep1}.folder,'/',files{timestep1}.name]);
year0 = files{timestep1+1}.name(end-14:end-11);
month0 = files{timestep1+1}.name(end-10:end-9);
if timestep1==1
    CtrlVar = CtrlVarInRestartFile;
    UserVar = UserVarInRestartFile;
    rho = F.rho; s = F.s; h = F.h; ub = F.ub; vb = F.vb; ab=F.ab;
    b = F.b; B=F.B;
    load([files{timestep1+1}.folder,'/',files{timestep1+1}.name],'time');
end
time0 = datestr(datenum(['01/',month0,'/',year0],'dd/mm/yyyy'),'dd/mm/yyyy');
    
MUA0 = MUA; GF0=GF;
CtrlVar.PlotXYscale = 1e3;
Fh0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),h);
Fs0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),s);
Fwct0 =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),b-B);
speed0 = sqrt(ub.^2+vb.^2);
Fspeed0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed0);
Fub0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ub);
Fvb0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),vb);
Fab0 =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ab);
hfPos = (S>B).*rhow.*(S-B)./rho ; % (positive) flotation thickness
hAF = (h>hfPos).*(h-hfPos) ;
FhAF0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hAF);

if timestep2 > length(files)
    timestep2 = length(files);
end
[files{timestep2}.folder,'/',files{timestep2}.name]
year1 = files{timestep2}.name(end-14:end-11);
month1 = files{timestep2}.name(end-10:end-9);
days1 = str2double(files{timestep2}.name(end-7:end-4));
datenum(['01/',month1,'/',year1]);
time1 = datestr(datenum(['01/',month1,'/',year1],'dd/mm/yyyy')+days1,'dd/mm/yyyy');
load([files{timestep2}.folder,'/',files{timestep2}.name],'h','s','ub','vb','MUA','GF','ab','b','B','S','rhow','rho');
deltah = h-Fh0(MUA.coordinates(:,1),MUA.coordinates(:,2));
hfPos = (S>B).*rhow.*(S-B)./rho ; % (positive) flotation thickness
hAF = (h>hfPos).*(h-hfPos) ;
deltahAF = hAF - FhAF0(MUA.coordinates(:,1),MUA.coordinates(:,2));                              
deltaspeed = sqrt(ub.^2+vb.^2)-Fspeed0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltaab = ab-Fab0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltawct = b-B-Fwct0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltawct(GF.node==1)=NaN;

Fspeed =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),sqrt(ub.^2+vb.^2));
Bcont = B; Bcont(GF.node<1)=NaN;
FBcont =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Bcont);
wctcont = b-B; wctcont(GF.node==1)=NaN;
Fwctcont = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),wctcont);
Fdeltaspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),deltaspeed);
ub(h<=1.5)=NaN; vb(h<=1.5)=NaN;
Fub = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ub);
Fvb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),vb);
Fdeltah = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),deltah);

xmin2 = min(MUA.coordinates(:,1)); xmax2 = max(MUA.coordinates(:,1));
ymin2 = min(MUA.coordinates(:,2)); ymax2 = max(MUA.coordinates(:,2));

[Xm,Ym] = ndgrid([xmin2:500:xmax2],[ymin2:500:ymax2]);
deltaspeedm = Fdeltaspeed(Xm,Ym);
deltahm = Fdeltah(Xm,Ym); 
speedm = Fspeed(Xm,Ym);
Bm = FBcont(Xm,Ym);
wctm = Fwctcont(Xm,Ym);
J = find(isnan(MUA.Boundary.x));
if isempty(J)
    J=length(MUA.Boundary.x)+1;
end
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(1:J(1)-1) MUA.Boundary.y(1:J(1)-1)]));
K = find(abs(deltaspeedm) < 20);
deltaspeedm(I) = NaN;
deltahm(I)=NaN;
du = Fub(Xm,Ym)-Fub0(Xm,Ym); du([I;K])= NaN;
dv = Fvb(Xm,Ym)-Fvb0(Xm,Ym); dv([I;K])= NaN;
speedm(I) = NaN;
Bm(I) = NaN;
wctm(I) = NaN;

    %% plotting 
    Hfig = fig('units','inches','height',60*12/72.27,'width',60*12/72.27,'fontsize',16,'font','Helvetica');
    
    %deltah(GF.node==1)=NaN;  
    %deltah(h<=1.5)=NaN;
    
    switch plotoption
        case 'icethickness_speed'

            tlo = tiledlayout(Hfig,1,2,"TileSpacing","compact");
            for tt=1:2
                hax(tt) = nexttile(tlo); hold(hax(tt),"on");
            end

            PlotNodalBasedQuantities_JDR(hax(1),MUA.connectivity,MUA.coordinates/1e3,-deltahAF); hold on;
            
            %CM = othercolor('Reds7'); %
            %CM=jet(64); CM(1,:)=[1 1 1];

            CM = othercolor('RdBu11');
            CM1 = othercolor('PuRd9',8); CM1 = flipdim(CM1,1);
            CM2 = othercolor('Blues9',8); %
            CM = [CM1 ; CM2]; CM = flipdim(CM,1);
            %CM=flipdim(CM,1);
            colormap(hax(1),CM); 
            
            cb1=colorbar(hax(1),'Position',[0.3 0.90 0.15 0.015],'Location','northoutside','AxisLocation','in','Ticks',[-150:50:150]);
            cb1.XColor='k';
            cb1.YColor='k';
            cb1.TickLength=0.04;
            cb1.FontSize=16;
            cb1.Label.String = '\DeltahAF [m]';
            caxis(hax(1),[-150 150]); 
            
            dspeed = Fdeltaspeed(MUA.coordinates(:,1),MUA.coordinates(:,2)); 
            
            PlotNodalBasedQuantities_JDR(hax(2),MUA.connectivity,MUA.coordinates/1e3,dspeed);
        
            colormap(hax(2),CM);
            cb2=colorbar(hax(2),'Position',[0.73 0.90 0.15 0.015],'Location','northoutside','AxisLocation','in','Ticks',[-1000:500:1000]);
            cb2.XColor='k';
            cb2.YColor='k';
            cb2.TickLength=0.04;
            cb2.FontSize=14;
            cb2.Label.String = '\Deltaspeed [m/yr]';
            caxis(hax(2),[-1000 1000]);

            for ii=1:2
                plot(hax(ii),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');

                CtrlVar.PlotIndividualGLs=0;
                [xGL,yGL] = PlotGroundingLines(CtrlVar,MUA0,GF0);
                plot(hax(ii),xGL/1e3,yGL/1e3,'k','linewidth',1.5);
                [xGL,yGL] = PlotGroundingLines(CtrlVar,MUA,GF);
                plot(hax(ii),xGL/1e3,yGL/1e3,'m','linewidth',1.5);
    
                grid(hax(ii),'on'); box(hax(ii),'on'); axis(hax(ii),'equal');
                xlim(hax(ii),[xmin xmax]); ylim(hax(ii),[ymin ymax]); %xticklabels(gca,{''});
            end

            yticklabels(hax(2),{}); 
            title(tlo,Experiment+" ("+string(time1)+" - "+string(time0)+")");


                
        case 'wcthickness_bathy'
            
            %subplot('position',[0.1 0.1 0.85 0.85]); hold on;
            
            ax1 = axes; hold on;

            %deltah(GF.node<1)=NaN;
            PlotNodalBasedQuantities_JDR(ax1,MUA.connectivity,MUA.coordinates/1e3,deltawct); hold on;

            [M,c]=contour(ax1,Xm/1e3,Ym/1e3,wctm,[0:200:1000]); c.Color=[0.7 0.7 0.7];
            %clabel(M,c,[0:200:1000],'color',[0.5 0.5 0.5],'labelspacing',500);

            Bcont(Bcont<-1600)=-1600; Bcont(Bcont>-200)=-200;
            scale = (-Bcont+max(Bcont,[],'omitnan'))/(max(Bcont,[],'omitnan')-min(Bcont,[],'omitnan'));
            Mask = GF.node;
	        alphaMask = GF.node; I = find(alphaMask<1);
	        alphaMask(I)=0;

            patch('faces',MUA.connectivity,'vertices',MUA.coordinates/1e3,...
            'FaceVertexCData',Mask,'CDataMapping','scaled','EdgeColor','none','FaceColor',[0.5 0.5 0.8],...
            'FaceAlpha','flat','FaceVertexAlphaData',alphaMask.*scale);

            [M2,c2]=contour(ax1,Xm/1e3,Ym/1e3,Bm,[-1600:100:-600]); c2.Color=[0.5 0.5 1];
            %clabel(M2,c2,[-1600:200:-200],'color',[0.3 0.3 1],'labelspacing',500); 
     
            plot(ax1,MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
            
            CtrlVar.PlotGLs=0; CtrlVar.PlotXYscale = 1e3;
            [xGL,yGL,~]=PlotGroundingLines(CtrlVar,MUA,GF);
            plot(ax1,xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'m','linewidth',1.5) ;
            [xGL0,yGL0,~]=PlotGroundingLines(CtrlVar,MUA0,GF0);
            plot(ax1,xGL0/CtrlVar.PlotXYscale,yGL0/CtrlVar.PlotXYscale,'k','linewidth',1.5) ;

            axis(ax1,'equal'); box(ax1,'on');
            xlim(ax1,[xmin xmax]); ylim(ax1,[ymin ymax]);

            set([ax1],'Position',[0.1 0.1 0.85 0.85]);

            %% Give each one its own colormap
            CM = othercolor('RdYlBu11',21); CM(11,:)=[1 1 1];
            colormap(ax1,CM); caxis(ax1,[-750 750]);

            if contains(basin,'AS')
                cb1=colorbar(ax1,'Position',[0.25 0.88 0.15 0.015],'Location','northoutside','AxisLocation','in','Ticks',[-750:250:750]);
                cb1.XColor='k';
                cb1.YColor='k';
                cb1.TickLength=0.04;
                cb1.FontSize=14;
                cb1.Label.String = '\Deltawct [m]';
            end

            %CM2 = othercolor('Bu_10',6); 
            %colormap(ax2,CM2); caxis(ax2,[-1600 -600]);
            
%            title('UaMITgcm');
            ylabel(ax1,'psy [km]','Fontsize',16);
            xlabel(ax1,'psx [km]','Fontsize',16);
            

        case 'melt_thickness'
            
            subplot('position',[0.1 0.1 0.43 0.85]); hold on;
            
            PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates/1e3,deltah); hold on;
            CM = othercolor('RdBu8',25); %CM=flipdim(CM,1); 
            CM(13,:)=[1 1 1];
            colormap(gca,CM); 
            cb1 = colorbar(gca,'fontsize',16);
            cb1.Label.String = '\Deltaice thickness [m]';
            caxis(gca,[-500 500]);
            
            plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
            PlotGroundingLines(CtrlVar,MUA0,GF0,[],[],[],'k','linewidth',1.5);
            PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'b','linewidth',1.5);

            grid on; box on; axis equal;
            xlim([xmin xmax]); ylim([ymin ymax]);

            xlabel('psx [km]','Fontsize',16); ylabel('psy [km]','Fontsize',16);
            
            subplot('position',[0.55 0.1 0.43 0.85]); hold on;
            
            PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates/1e3,-deltaab); hold on;
            
            %caxis([0 1]);

            %     CM = othercolor('RdBu11');
            CM1 = othercolor('PuRd9',32); CM1 = flipdim(CM1,1);
            CM2 = othercolor('Blues9',32); %
            CM = [CM1 ; CM2]; CM = flipdim(CM,1);
            colormap(gca,CM);
            cb2 = colorbar(gca,'fontsize',16);
            cb2.Label.String = '\Deltamelt [m/yr]';
            caxis(gca,[-60 60]);
            
            plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
            PlotGroundingLines(CtrlVar,MUA0,GF0,[],[],[],'k','linewidth',1.5);
            PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'b','linewidth',1.5);

            grid on; box on; axis equal;
            xlim([xmin xmax]); ylim([ymin ymax]);

            xlabel('psx [km]','Fontsize',16); yticklabels(gca,{});
            
            
            %contour(Xm/1e3,Ym/1e3,speedm,[0:500:4500],'color',[0.7 0.7 0.7],'linewidth',1);
            %contour(Xm/1e3,Ym/1e3,deltahm,[-400:100:400],'color',[0.7 0.7 0.7],'linewidth',1);
            
        otherwise 
            error('case unknown')
 
    end
            
    
    
    %title(['Change in ',plotoption,', ',time0,'-',time1],'fontsize',16,'interpreter','none');
    
pos = get(Hfig,'Position');
set(Hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);

folder = ['../Figures/',runID];
if ~exist(folder)
    mkdir(folder);
end
fname = [folder,'/Change_in_',plotoption,'_',basin];
print(Hfig,fname,'-dpng','-r400');