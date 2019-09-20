X = [-4.023796981061571
-0.9578464348659403
1.3145349536425714
2.3858025089950576
0.9497573705098933
-1.0934690723916378
-0.7305354723701374
1.7495472953141844
4.030276529234811
3.3688691464108764
1.6771515208713361
0.45489039740330844
1.0345439167623791
2.6627861687398098
2.8621401730217677
2.262894273291246
0.32692055679545107
-2.0106303174375517
-3.111872743074652
-2.5910855937934816
-2.3994133164815983
-2.7985669201540966
-1.566368675131364
-0.6812507266091545
-0.37919446969811926
-0.9526062813276055
-0.8544172813413926
0.2641992411485671
2.199737011722542
1.8411542076092697
-0.7834969042483988
-0.006506750260680372
4.076941459852892
4.628900176487166
2.357291983624029
2.9009761080836776
3.2262021476957843
0.12162982332908001
-2.4984201526047163
-3.5644490719548894
-5.109133843649704
-3.316077704571341
1.069302685869149
3.377604851313562
3.2124831646171907
0.9388052217319038
0.5024906457884648
-1.269586504404644
-2.170809763967145
-2.3029781052706886
0.7159836920613467
0.6865634863182508
-1.0074890744290086
-2.772814678850997
-2.2862645935904204
-1.1395505280621385
0.1413857609850674
-0.9862066619777826
-2.6407536780391467
-2.230776322004304
0.15483416696721775
3.2653258107307344
5.627349029331004
3.82895841997205
0.04701097621040673
1.8691896417948282
4.127191905766023
2.387493668565524
-0.10866292425911711
1.4277748520992382
2.0199366909655865
1.9876156208475129
1.988982672996143
2.697774301146138
2.9153082290827474
-0.020652352149987063
-2.4188460406275114
-4.103619415716151
-2.5079884088889077
-0.2859506536111971
1.9971873755537195
1.5688587420810616
2.917796591512333
4.66225509820341
5.274193846380819
3.9449844126940663
3.2812677667523893
1.1027755241656814
-0.40331548612038504
-0.4809942198359875
0.7211387971012782
2.5761334089116668
4.273612636520243
4.78155928787603
4.479372470204414
3.2693997261220034
1.6934283971945394
-0.2445305479009642
0.6279501235052223
1.8752680484740207];

n = size(X,1);
mu = mean(X);
v = var(X,1);




for i = 1:25
    sum = 0;
    for j = 1:n - i
        sum = sum + (X(j) - mu)*(X(j + i) - mu);
    end
    acf(i) = 1/n * sum / v;
end








pacf1 = acf(1);
pacf2 = det([1 acf(1); acf(1) acf(2)]) / det([1 acf(1); acf(1) 1]);
pacf3 = det([1 acf(1) acf(1); acf(1) 1 acf(2); acf(2) acf(1) acf(3)]) / det([1 acf(1) acf(2); acf(1) 1 acf(1); acf(2) acf(1) 1]);





x = size(acf,2);
mirroracf = [fliplr(acf) 1 acf];
for i = 1:x
    A = [];
    B = [];
    for k = 0:i - 1
        rowA = [mirroracf((x + 1) - k : (x + 1) - k + (i - 2))];        
        A = [A; rowA];
        rowB = [mirroracf((x + 1) - k : (x + 1) - k + (i - 1))];        
        B = [B; rowB];
    end
    A = [A mirroracf((x + 2) : (x + 2) + (i - 1))'];
    pacf(i) = det(A) / det(B);
end

    




CB(1) = 2/sqrt(n);
for i = 2:25
    sum = 1;
    for j = 1:i - 1
        sum = sum + 2*acf(j)^2;
    end
    CB(i) = 2/sqrt(n) * sqrt(sum);
end

figure('name','ACF')
subplot(1,2,1);
bar(linspace(1,25,25),acf)
hold on
plot(linspace(1,25,25),CB)
plot(linspace(1,25,25),-CB)
title('ACF');

subplot(1,2,2);
bar([1 2 3],[pacf1 pacf2 pacf3])
hold on
plot([1 2 3],[3/sqrt(n) 3/sqrt(n) 3/sqrt(n)])
plot([1 2 3],-[3/sqrt(n) 3/sqrt(n) 3/sqrt(n)])
title('PACF');



omegaplot = linspace(0,pi,100);
M = 25;
N = size(X,1);
cor = autocorr(X,25);
R = var(X,1).*cor;

triangular = 1/(2*pi)*(R(1));
for k = 1:M
    triangular = triangular + (1/pi) * (1 - k/M) * R(k + 1) * cos(k * omegaplot);
end

figure('name','Spectrum');
plot(omegaplot,triangular);
title('M = 25, triangular window');