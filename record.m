function [recordedBesty,recordedBestConV,bestx,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st] = record(offs,offsf,offsVio,recordedBesty,recordedBestConV,bestx,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st)
m = length(offsf);
for i = 1:m
    if isnan(FeasibleFes1st)&&offsVio(i)==0
        FeasibleFes1st = Fes;
        FeasibleObj1st = offsf(i);
        FeasibleX1st = offs(i,:);
    end
    if (offsf(i)<recordedBesty(Fes)&&offsVio(i)==recordedBestConV(Fes))||(offsVio(i)<recordedBestConV(Fes))
        bestx = offs(i,:);
        recordedBesty(Fes) = offsf(i);
        recordedBestConV(Fes) = offsVio(i);
    end
    Fes = Fes + 1;
    recordedBesty(Fes) = recordedBesty(Fes-1);
    recordedBestConV(Fes) = recordedBestConV(Fes-1);
end
end