function data_runID = load_PROPHET_data(froot_tools,froot_UaMITgcm,runID,couplingtimestep)

addpath(froot_tools);

frootm = froot_UaMITgcm+"/cases/"+runID;

if exist("AS_PROPHET_BasinDiagnostics.mat","file")
   load("AS_PROPHET_BasinDiagnostics.mat");
else
   data.runID = runID;  
end

IrunID = find(strcmp([data(:).runID],runID));

if isempty(IrunID)
   IrunID = length(data)+1;
   data(IrunID).runID = runID;
   data(IrunID).Basins = DefineBasins;
   Iend = 0;
else
   if isfield(data(IrunID),"Basins")
       Iend = length(data(IrunID).Basins(1).Time);
   else
       data(IrunID).Basins = DefineBasins;
       Iend=0;
   end
end

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = {subd(isub).name}';
nameFolds(ismember(nameFolds,[".",".."])) = [];

nfolders = length(nameFolds);

if nfolders>0

    output = []; files=[]; kk=1; oceanonly = [];

    for jj=1:nfolders
        %disp([frootm,"/output/",nameFolds{jj}]);
        output(jj).filelist=dir(frootm+"/output/"+nameFolds{jj}+"/Ua/UaDefaultRun_"+char(nameFolds{jj})+"*.mat");
        if ~isempty(output(jj).filelist)
            oceanonly(kk) = 0;
            for ii=1:length(output(jj).filelist)
                files{kk}.folder = output(jj).filelist(ii).folder;
                files{kk}.name = output(jj).filelist(ii).name;
                kk = kk+1;
            end
        else
            % standalone ocean
            oceanonly(kk) = 1;
            output(jj).filelist=dir(frootm+"/output/"+nameFolds{jj}+"/MITgcm/output.nc");
            for ii=1:length(output(jj).filelist)
                files{kk}.folder = output(jj).filelist.folder;
                files{kk}.name = output(jj).filelist(ii).name;
                kk = kk+1;
            end
        end
    end

    nfiles = length(files);

    for ii=Iend/couplingtimestep+1:nfiles

        if oceanonly(ii)
            MITpath = files{ii}.folder;
            Time = ncread(MITpath+"/"+files{ii}.name,'time');
            attvalue=ncreadatt(MITpath+"/"+files{ii}.name,'time','units');
            if strfind(attvalue,'seconds')
                epoch = erase(attvalue,'seconds since ');
                epochnum = datenum(epoch);
                MITTime = double(epochnum+Time/(24*60*60));
            elseif strfind(attvalue,'days')
                epoch = erase(attvalue,'days since ');
                epochnum = datenum(epoch);
                MITTime = double(epochnum+Time);
            else
                error('I do not recognise the time format in output.nc');
            end

            for tt=1:numel(MITTime)

                [LON,LAT,MeltMIT,~,~,~,~,~,~,~] = PlotMeltRates(MITpath,tt,[1 0 0 0 0 0 0 0]);
                dxMIT = LON(2,1)-LON(1,1);
                dyMIT = LAT(1,2)-LAT(1,1);
            
                for bb=1:length(data(IrunID).Basins)
                    NodeMIT = find(inpoly([LON(:),LAT(:)],[data(IrunID).Basins(bb).X(:) data(IrunID).Basins(bb).Y(:)]));
            %             figure; PlotMuaMesh(CtrlVar,MUA); hold on;
            %             plot(Basins(bb).X(:),Basins(bb).Y(:),"-b");
            %             plot(GLgeo(NodeI,7),GLgeo(NodeI,8),"or");
                    data(IrunID).Basins(bb).EleI = [];
                    data(IrunID).Basins(bb).NodeI = [];
                    if isfield(data(IrunID).Basins(bb),"VAF")
                        ind = numel(data(IrunID).Basins(bb).VAF);
                    else
                        ind = 0;
                    end
                    data(IrunID).Basins(bb).VAF(ind+1) = nan;
                    data(IrunID).Basins(bb).GroundedArea(ind+1) = nan;
                    data(IrunID).Basins(bb).BasalMassBalanceUa(ind+1) = nan;
                    data(IrunID).Basins(bb).BasalMassBalanceMITgcm(ind+1) = -sum(MeltMIT(NodeMIT)*dxMIT*dyMIT,"omitnan"); % in kg/s
                    data(IrunID).Basins(bb).SurfaceMassBalance(ind+1) = nan;
                    data(IrunID).Basins(bb).SurfaceMassBalanceGrounded(ind+1) = nan; %m3/yr
                    data(IrunID).Basins(bb).GLFlux(ind+1) = nan;
                    data(IrunID).Basins(bb).BoundaryFlux(ind+1) = nan;
                    data(IrunID).Basins(bb).Time(ind+1) = MITTime(tt);
                end

            end

        else

            MITpath = files{ii}.folder(1:end-3)+"/MITgcm";
            %[LON,LAT,TBottom] = PlotBottomTemperature(MITpath);
            T = mod(ii,couplingtimestep); if T==0; T=couplingtimestep; end
            [LON,LAT,MeltMIT,~,~,~,~,~,~,~] = PlotMeltRates(MITpath,T,[1 0 0 0 0 0 0 0]);
            dxMIT = LON(2,1)-LON(1,1);
            dyMIT = LAT(1,2)-LAT(1,1);
            load(files{ii}.folder+"/"+files{ii}.name,"MUA","B","h","S","rho","rhow","GF","CtrlVar","ub","vb","ab","as","s","b","dhdt");
            xEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
            yEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
            VAF = CalcVAF(CtrlVar,MUA,S,B,h,rho,rhow);  
            GroundedArea.ele = FEintegrate2D([],MUA,GF.node);
            BasalMassBalance=FEintegrate2D([],MUA,ab);
            SurfaceMassBalance=FEintegrate2D([],MUA,as);
            SurfaceMassBalanceGrounded = FEintegrate2D([],MUA,as.*GF.node);
            [LakeNodes,~]=LakeOrOcean2(CtrlVar,MUA,GF);
            [Flux,~,~,~,~,~,~,~,GLgeo,~,~]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,ub,vb,0*ub,0*vb,h,rho,[],[],[],[],LakeNodes);
            [xB,yB,FluxBound]=FluxAcrossModelBoundary(CtrlVar,MUA,ub,vb,h,b,B,rho);
            xBound = mean(xB,2,"omitnan"); yBound = mean(yB,2,"omitnan");
            L=strlength(files{ii}.name);
            yyyymm=extractBetween(files{ii}.name,L-14,L-9);
            dd=extractBetween(files{ii}.name,L-7,L-4);
            startdatenum = datenum(yyyymm+"01","yyyymmdd");

            for bb=1:length(data(IrunID).Basins)
                EleI = find(inpoly([xEle(:) yEle(:)],[data(IrunID).Basins(bb).X(:) data(IrunID).Basins(bb).Y(:)]));
                NodeI = find(inpoly([GLgeo(:,7) GLgeo(:,8)],[data(IrunID).Basins(bb).X(:) data(IrunID).Basins(bb).Y(:)]));
                BoundNodeI = find(inpoly([xBound(:) yBound(:)],[data(IrunID).Basins(bb).X(:) data(IrunID).Basins(bb).Y(:)]));
                NodeMIT = find(inpoly([LON(:),LAT(:)],[data(IrunID).Basins(bb).X(:) data(IrunID).Basins(bb).Y(:)]));
        %             figure; PlotMuaMesh(CtrlVar,MUA); hold on;
        %             plot(Basins(bb).X(:),Basins(bb).Y(:),"-b");
        %             plot(GLgeo(NodeI,7),GLgeo(NodeI,8),"or");
                data(IrunID).Basins(bb).EleI = EleI;
                data(IrunID).Basins(bb).NodeI = NodeI;
                if isfield(data(IrunID).Basins(bb),"VAF")
                    ind = numel(data(IrunID).Basins(bb).VAF);
                else
                    ind = 0;
                end
                data(IrunID).Basins(bb).VAF(ind+1) = sum(VAF.ele(EleI));
                data(IrunID).Basins(bb).GroundedArea(ind+1) = sum(GroundedArea.ele(EleI));
                data(IrunID).Basins(bb).BasalMassBalanceUa(ind+1) = sum(BasalMassBalance(EleI)); % m^3/yr
                data(IrunID).Basins(bb).BasalMassBalanceMITgcm(ind+1) = -sum(MeltMIT(NodeMIT)*dxMIT*dyMIT,"omitnan"); % in kg/s
                data(IrunID).Basins(bb).SurfaceMassBalance(ind+1) = sum(SurfaceMassBalance(EleI));
                data(IrunID).Basins(bb).SurfaceMassBalanceGrounded(ind+1) = sum(SurfaceMassBalanceGrounded(EleI)); %m3/yr
                data(IrunID).Basins(bb).GLFlux(ind+1) = sum(Flux(NodeI));
                data(IrunID).Basins(bb).BoundaryFlux(ind+1) = sum(FluxBound(BoundNodeI));
                data(IrunID).Basins(bb).Time(ind+1) = startdatenum + str2double(dd);
            end

        end

        disp(runID+": "+files{ii}.folder+"- done "+num2str(ii)+" out of "+num2str(nfiles));
    end

    % calculate dVAFdt
    for bb=1:length(data(IrunID).Basins)   
        dt = (data(IrunID).Basins(bb).Time(2:end)-data(IrunID).Basins(bb).Time(1:end-1))/365.25; % convert days to years
        dVAFdt= (data(IrunID).Basins(bb).VAF(2:end)-data(IrunID).Basins(bb).VAF(1:end-1))./dt; % in m3/yr
        data(IrunID).Basins(bb).dVAFdt = [NaN; dVAFdt(:)];
        data(IrunID).Basins(bb).dVAFdt = data(IrunID).Basins(bb).dVAFdt(1:length(data(IrunID).Basins(bb).Time));
    end

    data_runID = data(IrunID);
    save("AS_PROPHET_BasinDiagnostics.mat","data");

else
    disp("No data for "+runID+" - skipping");
    data_runID = [];
end
