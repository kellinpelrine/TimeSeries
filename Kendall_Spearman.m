close all;
clearvars;

e = [1.8175974111737432
0.3789054761159727
-1.3678437224172832
-1.8717190432604327
0.9065018557588502
2.265461001660006
-2.278318513482644
0.9754240488741482
-0.8259662269558103
-2.2428200258446815
1.0884016147054296
0.5697456749356591
0.19560644239122174
1.9410490224698607
-0.6003812336663493
-2.797890325341412
0.47897586475743725
-2.0470964924221886
-0.4569214659838777
0.0784757316899054
-0.13310341337868142
-0.20302100744598894
-0.981779180310813
0.07580039272518366
-0.994840971165401
-0.4984761067728956
0.09559920463066253
0.4143652056176483
1.0446308909164423
0.34544393593620953];

n = size(e,1);
t = linspace(1,n,n)';
options = optimoptions(@patternsearch,'Display','off','InitialMeshSize',.5);

aKdefault = patternsearch(@(a)a,1,[],[],[],[],-100,100,@(a)kendalltest(e,n,t,a,.05),options);   
aSdefault = patternsearch(@(a)a,1,[],[],[],[],-100,100,@(a)spearmantest(e,n,t,a,.05),options); 

stepsize = .0001;
i = 1;
for alpha = .01:stepsize:.1
    kendallamin(i) = patternsearch(@(a)a,1,[],[],[],[],-100,100,@(a)kendalltest(e,n,t,a,alpha),options);    
    spearmanamin(i) = patternsearch(@(a)a,1,[],[],[],[],-100,100,@(a)spearmantest(e,n,t,a,alpha),options);  
    i = i + 1;
end

xplot = linspace(0.01,0.1,i - 1);
plot(xplot,kendallamin,'b');
hold on
plot(xplot,spearmanamin,'r');
xlabel('alpha')
ylabel('a')
title('Kendall blue, Spearman red.')


function [c, ceq] = kendalltest(e,n,t,a,alpha)
    X = a*t + e;
    K = 0;
    for i = 1:n
        for j = i + 1:n
            if X(i) < X(j)
                K = K + 1;
            end
            if X(i) == X(j)
                K = K + 1/2;
            end
        end
    end

    T = 4*K/(n*(n - 1)) - 1;
    sigmaT = sqrt(2*(2*n + 5)/((9*n)*(n - 1)));
    tf = abs(T) > norminv((1 - alpha/2),0,1)*sigmaT;
    if tf == 1 && T > 0
        c = -1;
    else
        c = 1;
    end
    ceq = [];
end


function [c, ceq] = spearmantest(e,n,t,a,alpha)
    X = a*t + e;
    
    r = tiedrank(X);
    S = 1 - 6/(n*(n^2 - 1)) * sum((r - linspace(1,n,n)').^2);
    tf = abs(S) > norminv((1 - alpha/2),0,1)*sqrt(1/n);    
    if tf == 1 && S > 0
        c = -1;
    else
        c = 1;
    end
    ceq = [];
end