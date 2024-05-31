function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F)

persistent FC

if exist(CtrlVar.NameOfFileForReadingSlipperinessEstimate,"file")~=2
    error("C file "+UserVar.NameOfFileForReadingSlipperinessEstimate+" does not exist");
else
    if isempty(FC)
        load(CtrlVar.NameOfFileForReadingSlipperinessEstimate,'xC','yC','C');
        FC = scatteredInterpolant(xC,yC,C,'linear');
        fprintf("\n Read slipperiness from file "+CtrlVar.NameOfFileForReadingSlipperinessEstimate+" \n");
    end

    load(CtrlVar.NameOfFileForReadingSlipperinessEstimate,'m');
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    C = FC(x,y);
    m = m(1);

end

q = 1;      % only needed for Budd
muk = 0.5;  % required for Coulomb friction type sliding law as well as Budd, Tsai, Umbi and Cornford

end
