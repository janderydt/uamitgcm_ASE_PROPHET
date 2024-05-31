function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,F)

persistent Fsmb_RACMO_climatology DTbmelt
         
%
% When calculating dabdh from ab(b) for floating ice shelves:
% b=S-rho h /rhow
% h=rhow (S-b)/rho
% ab(b)=ab(S-rho h/rhow)
% dab/dh= -(rho/rhow) dab/db
% or:
% dab/dh = dab/db  db/dh = dab/db (-rho/rhow)= -(rho/rhow) dab/db 
          
x = MUA.coordinates(:,1); y=MUA.coordinates(:,2);

%% Surface mass balance: RACMO 2000-2018 climatology
    
if isempty(Fsmb_RACMO_climatology)

    load(UserVar.RACMO_SMB,"Fsmb_RACMO");
    
    % RACMO climatology between 2000 and 2018
    Istart = find(contains(string(Fsmb_RACMO.years),"2000"));
    Iend = find(contains(string(Fsmb_RACMO.years),"2018"));
    dn = Iend-Istart+1;
    RACMO_smb = 0;
    for ii=1:dn
        yrstr = "yr"+string(1999+ii);
        RACMO_smb = RACMO_smb + Fsmb_RACMO.(yrstr).Values;
    end
    Fsmb_RACMO_climatology = Fsmb_RACMO.climatology; Fsmb_RACMO_climatology.Values = RACMO_smb/dn;
end

as = Fsmb_RACMO_climatology(x,y);

if any(isnan(as))
	as(isnan(as)) = nanmean(as);
	fprint('Corrected for NaNs in surface mass balance');
end

%% Basal melt rates using input from MITgcm

%ab is defined on a grid (x,y)  I need to check in the future if coordinates are defined on the same grid and if not then do an
%interpolation

% if isempty(DTbmelt)
%  	DTbmelt = scatteredInterpolant(UserVar.UaMITgcm.MITgcmCGridX(:),...
%          UserVar.UaMITgcm.MITgcmCGridY(:),UserVar.UaMITgcm.MITgcmMelt(:),'linear');
% end
% 
% ab=DTbmelt(x,y);
ab=0*as;

% make sure that ab is zero unless 1) node is afloat and 2) node belongs to ocean
% [GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,F.GF);
% I = [find(GF.NodesCrossingGroundingLines(:)); find(GF.NodesUpstreamOfGroundingLines(:))];
% ab(I)=0;
% 
% [~,LakeNodes,~,~,~,~] = LakeOrOcean(CtrlVar,MUA,GF);
% ab(LakeNodes)=0;

%% Or simply:
[MeltNodes,~]=SpecifyMeltNodes(CtrlVar,MUA,F.GF); % assumes strict treatment of melt nodes, i.e. only nodes of fully floating elements are melting
ab(~MeltNodes) = 0;

h=F.s-F.b;
I=(h<CtrlVar.ThickMin+0.5 & ab<0); ab(I)=0;  dabdh(I)=0;% do not melt ice shelf away where ice thickness is less than 1m.

end
