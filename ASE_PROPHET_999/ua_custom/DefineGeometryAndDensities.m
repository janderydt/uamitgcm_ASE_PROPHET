function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

persistent Fs FB Fb Frho

s=[]; b=[]; S=[]; B=[];
alpha=0 ;

if nargin<5
    FieldsToBeDefined="-s-b-S-B-rho-";
end

MUAnew = MUA;

fprintf("Loading "+FieldsToBeDefined+"\n");

if isempty(Fs)
    load(UserVar.GeometryInterpolants);
end

%% Step I. Bed
if contains(FieldsToBeDefined,'B')
    
    B = FB(MUA.coordinates(:,1),MUA.coordinates(:,2));
    fprintf('Done B \n');
    
    if any(isnan(B))
        error('NaN values in B');
    end
    
end

%% Step II. sea surface
if contains(FieldsToBeDefined,'S')

    S = zeros(MUA.Nnodes,1);
    fprintf('Done S \n');
    
end

%% Step III. ice surface and draft
if contains(FieldsToBeDefined,'s') 

    s = Fs(MUA.coordinates(:,1),MUA.coordinates(:,2));
    fprintf('Done s \n');

    if any(isnan(s))
        error('NaN values in s');
    end

end

if contains(FieldsToBeDefined,'b') 

    b = Fb(MUA.coordinates(:,1),MUA.coordinates(:,2));
    %h_init = Fh_init(MUA.coordinates(:,1),MUA.coordinates(:,2));
    fprintf('Done b \n');

    if any(isnan(b))
        error('NaN values in b');
    end

end

h = s-b;
I = find(h<=CtrlVar.ThickMin);
s(I) = b(I)+CtrlVar.ThickMin;


% Defines densities
if contains(FieldsToBeDefined,{'rho','rhow','g'})

    fprintf('Start loading densities \n');

    if isempty(Frho)
        load(UserVar.DensityInterpolant,'Frho');
    end
       
    rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
    rho(rho<100)=100;
    rho(rho>917)=917;
    
    rhow = 1027;
    
    g=9.81/1000;
    
    fprintf('Done loading densities \n');

end


end



