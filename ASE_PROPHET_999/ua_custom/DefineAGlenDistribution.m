function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,F)

persistent FAGlen

if exist(CtrlVar.NameOfFileForReadingAGlenEstimate,"file")~=2
    error("AGlen file "+ UserVar.NameOfFileForReadingAGlenEstimate+" does not exist");
else
    if isempty(FAGlen)
        load(CtrlVar.NameOfFileForReadingAGlenEstimate,'xA','yA','AGlen');
        FAGlen = scatteredInterpolant(xA,yA,AGlen,'linear');
        fprintf("\n Read rate factor from file "+CtrlVar.NameOfFileForReadingAGlenEstimate+" \n");
    end
    
    load(CtrlVar.NameOfFileForReadingAGlenEstimate,'n');
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    AGlen=FAGlen(x,y);
    n = n(1);
end
end
