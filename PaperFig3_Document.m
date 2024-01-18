%% PEA model fit block designed dataset

clear;
clc;
close all;

data = readtable("Exp2.csv");

Deg = unique(data.perturbSize);
[group,perturb,sub,type,Num,hand] = deal(data.group,data.perturbSize,data.sub,...
    data.trialType,data.cycleNum,data.handangle);
for g = 1:length(Deg)
    alldata{g} = reshape(hand(group==g),120,[])';
    Base{g} = nanmean(alldata{g}(:,21:30),2);
    Data(g,:) = nanmean(alldata{g}-Base{g});
    Pertur(g,:) = [zeros(1,30),ones(1,80)*Deg(g),NaN(1,10)];
end


FitRange = [30,120];

func = @(x)BCCTbT_fitting(x,Pertur,Data,FitRange);
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500, ...
    'FiniteDifferenceStepSize',1e-5);

% ini = [0.8 0.2 60 30 0.309 1.853];%initiate point
% LB = [0 0 0 0 0.309 1.853];%lowerbnd
% UB = [1 1 200 200 0.309 1.853];%upperbnd

ini = [0.8 0.2 60 30 0.1 4.2];%initiate point
LB = [0 0 0 0 0.1 0.1 ];%lowerbnd
UB = [1 1 200 200 10 10];%upperbnd

FitParam = fmincon(func,ini,[],[],[],[],LB,UB,[],options);

for g = 1:length(Deg)
    Hand_real(g,:) = BCCTbT(FitParam,Pertur(g,:),size(Pertur,2));
end

SSres = sum(sum((Hand_real(:,FitRange(1):FitRange(2))-Data(:,FitRange(1):FitRange(2))).^2));
SStot = sum(sum((Data-nanmean(Data(:))).^2));
Rsquared = 1-SSres/SStot;

RMS = sqrt(nanmean(nanmean((Hand_real(:,FitRange(1):FitRange(2))-Data(:,FitRange(1):FitRange(2))).^2)));
TrialNum = 1:1:size(Hand_real,2);

deg = 0:1:100;
for g = 1:length(deg)
    angle = deg(g);
    pertur(g,:) = ones(1,80)*angle;
    hand_real = BCCTbT(FitParam,pertur(g,:),80);
    HandSimu(g) = nanmean(hand_real(end-10:end));
   
end


Cols = turbo(length(Deg)+10);
h = figure('color','w','position',[100,100,1200,500]);

subplot(2,4,8)
plot(deg,HandSimu,'linewidth',3,'color',[190, 3, 252,150]./255);hold on;

for g = 1:length(Deg)

    subplot(2,4,g);
    text(8,23,sprintf('%.2f ^o',Deg(g)),'fontsize',12);hold on;
    plot([0,150],[0,0],'color',[0,0,0],'linewidth',1,'linestyle','--');hold on;

    patch([30.5,110.5,110.5,30.5],[-10,-10,25,25],[0.7,0.7,0.7],'facealpha',0.2,'edgecolor','none');

    PatchPlot(TrialNum,alldata{g}-Base{g},Cols(g,:));hold on;
    scatter(TrialNum,Data(g,:),10,Cols(g,:),'filled','markeredgecolor','none', ...
        'markerfacealpha',0.3);hold on;
    plot(TrialNum,Data(g,:),'color',[Cols(g,:),0.3]);hold on;
    plot(TrialNum(FitRange(1)+1:FitRange(2)),Hand_real(g,FitRange(1)+1:FitRange(2)),'color',Cols(g,:)*0.7,'linewidth',2);hold on;
    plot(TrialNum(FitRange(2)+1:120),Hand_real(g,FitRange(2)+1:120),'color',Cols(g,:)*0.7,'linewidth',2,'linestyle','--');hold on;
    set(gca,'fontsize',14,'xtick',[30,110],'ytick',0:10:25);hold on;
    box off;
    axis([0,120,-5,25]);


    subplot(2,4,8)
    EndData = nanmean(alldata{g}(:,101:110)-Base{g},2);
    scatter(Deg(g),EndData ,10,Cols(g,:),'filled','markeredgecolor','none', ...
        'markerfacealpha',0.5);hold on;
    scatter(Deg(g),nanmean(Data(g,101:110)),30,Cols(g,:),'filled','markeredgecolor','none', ...
        'markerfacealpha',0.8,'Marker','d');hold on;
    errorbar(Deg(g),nanmean(EndData),SE(EndData),'capsize',0, ...
        'color',Cols(g,:)*0.7,'linewidth',2);hold on;


end

subplot(2,4,8)
title('BCC model');
xlabel('Perturbation size (deg)');
ylabel('Adapatation extent (deg)');
set(gca,'fontsize',12,'xtick',0:20:100,'ytick',0:5:30,'XTickLabelRotation',0);hold on;
box off;
% grid on;
axis([-4,100,-0.75,35]);

    axes('Position',[0.85 0.3 0.05 0.14]);
    plot(deg,HandSimu,'linewidth',3,'color',[190, 3, 252,150]./255);hold on;
    for g = 1:3
     EndData = nanmean(alldata{g}(:,101:110)-Base{g},2);
       scatter(Deg(g),EndData ,10,Cols(g,:),'filled','markeredgecolor','none', ...
        'markerfacealpha',0.5);hold on;
    scatter(Deg(g),nanmean(Data(g,101:110)),30,Cols(g,:),'filled','markeredgecolor','none', ...
        'markerfacealpha',0.8,'Marker','d');hold on;
    errorbar(Deg(g),nanmean(EndData),SE(EndData),'capsize',0, ...
        'color',Cols(g,:)*0.7,'linewidth',2);hold on;
    end
    box on;
    set(gca,'XTick',[2,4,8],'YTick',[10,20,30]);
    axis([0,10,0,35]);hold on;
%%
function RMSE = BCCTbT_fitting(param,Pertur,Data,FitRange)

sError = [];
for Group = 1:size(Data,1)
    pertur = Pertur(Group,:);
    data = Data(Group,:);
    Index = FitRange(1);
    Index2 = FitRange(2);
    Hand_real = BCCTbT(param,pertur,length(data));
    sError = [sError,(Hand_real(Index:Index2)-data(Index:Index2)).^2];
end

RMSE = sqrt(nanmean(sError));

end
function Hand_real = BCCTbT(param,pertur,trialNum)
a = param(1);
b = param(2);
sig_h = param(3);
sig_g = param(4);
slope = param(5);
intercept = param(6);

s_g = 0;
Hand_real = zeros(1,trialNum);
Hand_est = zeros(1,trialNum);

for t = 1:trialNum
    s_v = -pertur(t);
    sig_v = (slope*pertur(t)+intercept).^2;
    s_h = Hand_real(t);
    if ~isnan(s_v)
        [c_mean,~] = CueCombine(s_h,s_v,s_g,sig_h,sig_v,sig_g);
    else
        [c_mean,~] = CueCombine(s_h,s_g,0,sig_h,sig_g,0,1);
    end
    Hand_est(t) = c_mean;

    if t<trialNum
        Hand_real(t+1) = a*Hand_real(t)+b*(0-Hand_est(t));
    end
end

end
function [c_mean,c_sig] = CueCombine(s_1,s_2,s_3,sig_1,sig_2,sig_3,twocues)

s_12 = s_1+1/(sig_2/sig_1+1)*(s_2-s_1);
sig_12 = sig_2*sig_1/(sig_2+sig_1);
s_123 = s_12+1/(1+sig_3/sig_12)*(s_3-s_12);
sig_123 = sig_12*sig_3/(sig_12+sig_3);

c_mean = s_123;
c_sig = sig_123;
if nargin>6
    if twocues==1
        c_mean = s_12;
        c_sig = sig_12;
    end
end

end
function p = PatchPlot(x,y,Color,dim)
    if nargin == 2
        dim = 1;
        Color = 'r';
    elseif nargin == 3
        dim = 1;
    end
    m = nanmean(y,dim);
    se = SE(y,dim);

    X = [x,fliplr(x)];
    Y = [m+se,fliplr(m-se)];
    

    patch(X,Y,Color,'edgecolor','none','facealpha',0.2);hold on;
end