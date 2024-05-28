function plot_BoundaryConditions

addpath(getenv("froot_tools"));

froot_OB = getenv("froot_uamitgcm")+"cases/AS_PROPHET_999/mitgcm_run/input/";
froot_mask = getenv("froot_uamitgcm")+"cases/AS_PROPHET_001/output/200101/MITgcm/output.nc";
ny =  360; nz = 90;
suffix = ["rcp85","paris2c","baseline"];

% load geometry mask
hFacW = ncread(froot_mask,"hFacW"); hFacW(hFacW==0)=nan;
hFacW = squeeze(hFacW(2,:,:)); % only retain western boundary: ny*nz matrix
Zc = ncread(froot_mask,"Z");

colors = ["Oranges","Greens","Greys"];

H=fig('units','inches','width',50*12/72.27,'height',40*12/72.27,'fontsize',14,'font','Helvetica');
subplot('position',[0.1 0.1 0.85 0.85]); hold on;

for ii=1:numel(suffix)
    file=froot_OB+"OBWt_"+suffix(ii)+".bin";
    fid=fopen(file,'r');
    OBWt=fread(fid,'real*4','b');
    n = numel(OBWt); nt = n/(ny*nz);
    OBWt=reshape(OBWt,ny,nz,nt);
    CM = brewermap(2*nt,colors(ii));

    T1 = OBWt(:,:,1).*hFacW; 
    T1profile = mean(T1,1,"omitmissing");
    T1profile(isnan(T1profile))=0;
    Tend =  OBWt(:,:,end).*hFacW;
    Tendprofile = mean(Tend,1,"omitmissing");
    Tendprofile(isnan(Tendprofile))=0;

    patch("XData",[T1profile(:); flipdim(Tendprofile(:),1)],"YData",[Zc(:); flipdim(Zc(:),1)],...
        "EdgeColor",CM(nt,:),"FaceColor",CM(nt,:),"FaceAlpha",0.3);

    dt = 1;
    for tt=1:dt:nt % plot profile every year
        T = OBWt(:,:,tt);
        T = T.*hFacW; % apply mask
        T_profile = mean(T,1,"omitmissing");
        plot(T_profile,Zc,'color',CM(floor(nt/4)+tt,:),'linewidth',0.5);
    end
    %T_mean = mean(mean(OBWt.*repmat(hFacW,1,1,nt),3,"omitmissing"),1,"omitmissing");
    %g(ii)=plot(T_mean,Zc,'color',CM(end,:),'linewidth',2);

end

%% Add De Rydt avmelt profile
froot_OB = getenv("froot_uamitgcm")+"cases/PTDC_003/mitgcm_run/input/";
froot_mask = getenv("froot_uamitgcm")+"cases/PTDC_003/output/202001/MITgcm/output.nc";
ny =  360; nz = 90;

% load geometry mask
hFacW = ncread(froot_mask,"hFacW"); hFacW(hFacW==0)=nan;
hFacW = squeeze(hFacW(2,:,:)); % only retain western boundary: ny*nz matrix
Zc = ncread(froot_mask,"Z");

file=froot_OB+"OBWt.bin";
fid=fopen(file,'r');
OBWt=fread(fid,'real*4','b');
n = numel(OBWt); nt = n/(ny*nz);
OBWt=reshape(OBWt,ny,nz,nt);
    
T = OBWt(:,:,end).*hFacW; 
Tprofile = mean(T,1,"omitmissing");
g(4)=plot(Tprofile,Zc,':b','linewidth',2);

%% Add De Rydt himelt profile
froot_OB = getenv("froot_uamitgcm")+"cases/PTDC_002/mitgcm_run/input/";
froot_mask = getenv("froot_uamitgcm")+"cases/PTDC_002/output/280001/MITgcm/output.nc";
ny = 432; nz = 90;

% load geometry mask
hFacW = ncread(froot_mask,"hFacW"); hFacW(hFacW==0)=nan;
hFacW = squeeze(hFacW(2,:,:)); % only retain western boundary: ny*nz matrix
Zc = ncread(froot_mask,"Z");

file=froot_OB+"OBWt.bin";
fid=fopen(file,'r');
OBWt=fread(fid,'real*4','b');
n = numel(OBWt); nt = n/(ny*nz);
OBWt=reshape(OBWt,ny,nz,nt);
    
T = OBWt(:,:,end).*hFacW; 
Tprofile = mean(T,1,"omitmissing");
g(5)=plot(Tprofile,Zc,'--b','linewidth',2);


for ii=1:numel(suffix)
    uistack(g(ii),'top');
end

grid on; box on;

legend(g(:),[suffix,"avmelt","himelt"],'Location','northeast');

ylim([-800 0]); 
ylabel('Depth [m]'); xlabel('Temperature [C]');

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = "../Figures/BoundaryTemperature_vs_DeRydt2024";
print(H,fname,'-dpng','-r400');
