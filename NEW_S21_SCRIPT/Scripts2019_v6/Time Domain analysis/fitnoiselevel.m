function [afit,sigmalevel] = fitnoiselevel(x,y,a,b,c,inputlog)

if inputlog==1    
    y=10.^(y/10);
    a=10.^(a/10);
end
%fit in real space, not in log space
%s = fitoptions('Method','LinearLeastSquares', 'Startpoint',[a],'Lower',[a/3],...
%    'Upper',[a*3],'MaxFunEvals',10000,'TolFun',1e-21,'TolX',1e-21);
s = fitoptions('Method','NonlinearLeastSquares','Startpoint',a);

ftype = fittype(['a./((1+(x.*2*3.14159*' num2str(b) ').^2))+' num2str(c)],'options', s);
[result bla]=fit(x,y,ftype)

afit=result.a;
result2=confint(result);%confint extracts the 95% confidence intervals from result
sigmalevel=(afit-result2(1,1))/1.96;%convert 95% interval to standard deviation

if inputlog==1
    afit=10*log10(result.a);
    sigmalevel=10*log10(sigmalevel);
end
    




end%function