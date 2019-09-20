clearvars
clear figures
close all

data = csvread('SVAR_Data.csv'); %note: I've removed the first column (dates) and row (labels) from this file for easier processing
[trend , data] = hpfilter(data,1600);
GDP87 = data(:,1);
USAPGDP = data(:,2);
PCOM = data(:,3);
FF = data(:,4);
NBR = data(:,5);
TOTR = data(:,6);
M1 = data(:,7);
M2 = data(:,8);
M3 = data(:,9);

beta1 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], GDP87(3:146));
lagtest2 = mvregress([ones(143,1), GDP87(3:145), USAPGDP(3:145), PCOM(3:145), FF(3:145), NBR(3:145), TOTR(3:145), M1(3:145), ...
    GDP87(2:144), USAPGDP(2:144), PCOM(2:144), FF(2:144), NBR(2:144), TOTR(2:144), M1(2:144), ...
    GDP87(1:143), USAPGDP(1:143), PCOM(1:143), FF(1:143), NBR(1:143), TOTR(1:143), M1(1:143)], GDP87(4:146));

error = mean((GDP87(3:146) - [ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)] * beta1 ).^2 );
error2 = mean((GDP87(4:146) - [ones(143,1), GDP87(3:145), USAPGDP(3:145), PCOM(3:145), FF(3:145), NBR(3:145), TOTR(3:145), M1(3:145), ...
    GDP87(2:144), USAPGDP(2:144), PCOM(2:144), FF(2:144), NBR(2:144), TOTR(2:144), M1(2:144), ...
    GDP87(1:143), USAPGDP(1:143), PCOM(1:143), FF(1:143), NBR(1:143), TOTR(1:143), M1(1:143)] * lagtest2 ) .^2 );

disp('Error in fit with 2 lags - with 3 lags:')
error - error2

beta2 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], USAPGDP(3:146));
beta3 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], PCOM(3:146));
beta4 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], FF(3:146));
beta5 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], NBR(3:146));
beta6 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], TOTR(3:146));
beta7 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], M1(3:146));

beta = [beta1 beta2 beta3 beta4 beta5 beta6 beta7];
datamat = [data(2:145,1:7) data(1:144,1:7)];
residual = data(3:146,1:7) - ( beta(2:15,:)' * [data(2:145,1:7) data(1:144,1:7)]' )';
sigmahat = zeros(7,7);
for i = 1:144
    sigmahat = sigmahat + residual(i,:)' * residual(i,:);
end
sigmahat = sigmahat / 144;
A = chol(sigmahat,'lower');
Ainvert = inv(A);

G = [beta(2:15,:)' ; eye(7), zeros(7)];
eta = [0 0 0 0 0 0 1]';
u = [Ainvert * eta; zeros(7,1)];
M = 50; %period over which to computing IRF
shock = zeros(14,M);
for i = 1:M
    shock(:,i) = G^i * u;
end

figure(1)
title('IRFs')
hold on
for i = 1:7
    plot(shock(i,1:M)) 
end


%bootstrap
oldshock = shock;
shockbootstrap = zeros(14,M,1000);
for i = 1:1000
    newresidual = zeros(144,7);
    for j = 1:144
        newresidual(j,:) = residual(randi([1 144]),:);
    end
    newseries = zeros(144,7);
    newseries(1:2,:) = data(1:2,1:7);
    for j = 3:146
        newseries(j,:) = [newseries(j - 1,:) newseries(j - 2,:)] * beta(2:15,:) + newresidual(j - 2,:);
    end
    
    GDP87 = newseries(:,1);
    USAPGDP = newseries(:,2);
    PCOM = newseries(:,3);
    FF = newseries(:,4);
    NBR = newseries(:,5);
    TOTR = newseries(:,6);
    M1 = newseries(:,7);

    beta1 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], GDP87(3:146));
    beta2 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], USAPGDP(3:146));
    beta3 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], PCOM(3:146));
    beta4 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], FF(3:146));
    beta5 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], NBR(3:146));
    beta6 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], TOTR(3:146));
    beta7 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
        GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], M1(3:146));

    newbeta = [beta1 beta2 beta3 beta4 beta5 beta6 beta7];
    datamat = [newseries(2:145,:) newseries(1:144,:)];
    newresidual = newseries(3:146,:) - ( newbeta(2:15,:)' * datamat' )';
    newsigmahat = zeros(7,7);
    for j = 1:144
        newsigmahat = newsigmahat + newresidual(j,:)' * newresidual(j,:);
    end
    newsigmahat = newsigmahat / 144;
    newA = chol(newsigmahat,'lower');
    newAinvert = inv(newA);

    newG = [newbeta(2:15,:)' ; eye(7), zeros(7)];
    eta = [0 0 0 0 0 0 1]';
    newu = [newAinvert * eta; zeros(7,1)];
    M = 50; %period over which to computing IRF
    shock = zeros(14,M);
    for j = 1:M
        shock(:,j) = newG^j * newu;
    end
    
    shockbootstrap(:,:,i) = shock;
end

quantileup = zeros(7,50);
quantiledown = zeros(7,50);
for i = 1:7
    quantileup(i,:) = quantile(shockbootstrap(i,:,:),.9,3);
    quantiledown(i,:) = quantile(shockbootstrap(i,:,:),.1,3);  
end

figure(2)
hold on
title('GDP87 with bootstrap confidence interval')
plot(oldshock(1,1:M),'blue') 
plot(quantileup(1,:),'red')
plot(quantiledown(1,:),'red')

figure(3)
hold on
title('USAPGDP with bootstrap confidence interval')
plot(oldshock(2,1:M),'blue') 
plot(quantileup(2,:),'red')
plot(quantiledown(2,:),'red')

figure(4)
hold on
title('PCOM with bootstrap confidence interval')
plot(oldshock(3,1:M),'blue') 
plot(quantileup(3,:),'red')
plot(quantiledown(3,:),'red')

figure(5)
hold on
title('FF with bootstrap confidence interval')
plot(oldshock(4,1:M),'blue') 
plot(quantileup(4,:),'red')
plot(quantiledown(4,:),'red')

figure(6)
hold on
title('NBR with bootstrap confidence interval')
plot(oldshock(5,1:M),'blue') 
plot(quantileup(5,:),'red')
plot(quantiledown(5,:),'red')
        
figure(7)
hold on
title('TOTR with bootstrap confidence interval')
plot(oldshock(6,1:M),'blue') 
plot(quantileup(6,:),'red')
plot(quantiledown(6,:),'red')

disp('General comments: the shocks seems to persist around 20-30 periods, though the biggest effect is in the 1-10 range.')
disp('Most of the IRFs seem to exhibit a rebound around 15-25 (a response with the opposite sign of the earlier effect, for reasons unclear.')
disp('Perhaps this is a mathematical artifact.')


%test a big shock
G = [beta(2:15,:)' ; eye(7), zeros(7)];
bigeta = [0 0 0 0 0 0 5]';
bigu = [Ainvert * bigeta; zeros(7,1)];
M = 50; %period over which to computing IRF
bigshock = zeros(14,M);
for i = 1:M
    bigshock(:,i) = G^i * bigu;
end

figure(8)
title('GDP87, small shock in blue, big in red')
hold on
plot(oldshock(1,1:M),'blue')
plot(bigshock(1,1:M),'red')
disp('Big shock: as noted, the effect of a 5 times larger shock is more than 5 times larger.')
disp('The rebound mentioned above exhibits this behavior as well.')


%experiment with a different order
GDP87 = data(:,1);
USAPGDP = data(:,2);
PCOM = data(:,3);
FF = data(:,4);
NBR = data(:,5);
TOTR = data(:,6);
M1 = data(:,7);
M2 = data(:,8);
M3 = data(:,9);

beta1 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], GDP87(3:146));
beta2 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], FF(3:146));
beta3 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], PCOM(3:146));
beta4 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], USAPGDP(3:146));
beta5 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], NBR(3:146));
beta6 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], TOTR(3:146));
beta7 = mvregress([ones(144,1), GDP87(2:145), USAPGDP(2:145), PCOM(2:145), FF(2:145), NBR(2:145), TOTR(2:145), M1(2:145), ...
    GDP87(1:144), USAPGDP(1:144), PCOM(1:144), FF(1:144), NBR(1:144), TOTR(1:144), M1(1:144)], M1(3:146));

orderbeta = [beta1 beta2 beta3 beta4 beta5 beta6 beta7];
orderdatamat = [data(2:145,1) data(2:145,4) data(2:145,3) data(2:145,2) data(2:145,5:7) ... 
    data(1:144,1) data(1:144,4) data(1:144,3) data(1:144,2) data(1:144,5:7)];
orderresidual = [data(3:146,1) data(3:146,4) data(3:146,3) data(3:146,2) data(3:146,5:7)] - ( orderbeta(2:15,:)' * orderdatamat' )';
ordersigmahat = zeros(7,7);
for i = 1:144
    ordersigmahat = ordersigmahat + residual(i,:)' * residual(i,:);
end
ordersigmahat = ordersigmahat / 144;
orderA = chol(ordersigmahat,'lower');
orderAinvert = inv(orderA);

orderG = [beta(2:15,:)' ; eye(7), zeros(7)];
eta = [0 0 0 0 0 0 1]';
orderu = [orderAinvert * eta; zeros(7,1)];
M = 50; %period over which to computing IRF
ordershock = zeros(14,M);
for i = 1:M
    ordershock(:,i) = G^i * orderu;
end

figure(9)
title('switching order of USAPGDP and FF (two lines, old in blue, new order in red)')
hold on
plot(oldshock(1,1:M),'blue')
plot(ordershock(1,1:M),'red')
disp('Switching order: as expected, changing the order here makes no difference.')

