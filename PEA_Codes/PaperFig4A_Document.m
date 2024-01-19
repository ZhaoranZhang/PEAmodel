
%% load data
clear;


data = readtable("Tsay_report.csv");

Angle = abs(data.cond(1));
sub = data.SN;
hand = data.hand_theta;
prop = data.hand_report;
Cond = data.cond;

Hand = [];
Report = [];
for s = 1:max(sub)
    if mean(Cond(sub==s))>0
        d = -hand(sub==s)';
        r = -prop(sub==s)';
    else
        d = hand(sub==s)';
        r = prop(sub==s)';
    end
    IsOutlier = OutlierMAD (d,2.5);
    d(IsOutlier) = NaN;
    % bin the data
    Hand(s,:) = nanmean(reshape(d,4,[]));
    Report(s,:) = nanmean(reshape(r,4,[]));
end


handreal = nanmean(Hand(:,6:85))';
handreport = nanmean(Report(:,6:85))';

%% fit with PEA model
Range = 21:70;
trialNum = length(handreal);
pertur = [zeros(1,20),ones(1,50)*Angle,NaN(1,10)];


func = @(x)ReportHand(x,pertur,handreal,handreport,Range);


%parameters: a b R1 R2

% u = 1/sig_u v = 1/sig_v p = 1/sig_h
% R1 = v/(u+v+p);
% R2 = p/(u+v+p);
% hand_est(t) = R1*s_v+R2*x(t);
% x(t+1) = -b*R1*s_v+(-b*R2+a)*x(t)
% sig_sq = 1/(u+v+p);
% Hand_report(t) = R1*s_v/(1+R2)+2*R2*x(t)/(1+R2);

ini = [0.5  0.2  0.5  0.5];
LB =  [0.01   0.01 0.05 0.05];
UB =  [0.999 0.999    0.95 0.95];
FitParam = fmincon(func,ini,[],[],[],[],LB,UB);
 [handreal_sim,handreport_sim,handest_sim] = ReportHand_Simu(FitParam,pertur,trialNum);

Residue = [handreal_sim(Range)-handreal(Range)' handreport_sim(Range)-handreport(Range)'];
Residue = Residue(:);
Residue = Residue(~isnan(Residue));
NumofParam = size(FitParam,2);

pd = fitdist(Residue,'normal');
[aic,bic] = aicbic(-pd.NLogL,NumofParam,1);


SSres = sum(Residue.^2);
Data = [handreal(Range)' handreport(Range)'];
SStot = sum((Data-nanmean(Data)).^2);
Rsquared = 1-SSres/SStot;
RMS = sqrt(nanmean(Residue.^2));

cols = turbo(20);
figure('Position',[100,100,450,350],'color','w');
set(gca,'fontsize',14);
hold on;
plot([0 trialNum+1],[0 0],'color',[0.8,0.8,0.8]);
plot([20.5 20.5],[-10 30],'color',[0.8,0.8,0.8]);
plot([70.5 70.5],[-10 30],'color',[0.8,0.8,0.8]);
tsayCol1 = [56,114,70]./255;
tsayCol2 = [156,89,160]./255;

scatter(1:trialNum,handreal,30,'o','MarkerFaceColor',cols(1,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(1:trialNum,handreport,30,'o','MarkerFaceColor',cols(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
hhandPerc = plot(1:70,handest_sim(1:70),'color','r','LineWidth',2);
hhand = plot(1:70,handreal_sim(1:70),'color',cols(1,:),'LineWidth',2);
hVerbal = plot(1:70,handreport_sim(1:70),'color',cols(3,:),'LineWidth',2);

plot(70:80,handest_sim(70:80),'color','r','LineWidth',2,'LineStyle',':');
plot(70:80,handreal_sim(70:80),'color',cols(1,:),'LineWidth',2,'LineStyle',':');
plot(70:80,handreport_sim(70:80),'color',cols(3,:),'LineWidth',2,'LineStyle',':');

title('BCC model fitting');
xlabel('Movement Cycle');
ylabel('Hand Angle (deg)');
axis([-1,trialNum+1,-7,22]);
legend([hhand hVerbal hhandPerc],{'Actual','Report','Bias'},...
    'Location','northwest','box','off','fontsize',12);



%%
function [Hand_real,Hand_report,Hand_est] = ReportHand_Simu(param,pertur,trialNum)

Hand_real = zeros(1,trialNum);
Hand_report = zeros(1,trialNum);
Hand_est = zeros(1,trialNum);
a = param(1);
b = param(2);
R1 = param(3);
R2 = param(4);

for t = 1:trialNum
    s_v = -pertur(t);
    s_h = Hand_real(t);

    if ~isnan(s_v)
        Hand_est(t) = R1*s_v+R2*s_h;
        Hand_report(t) = R1*s_v/(1+R2)+2*R2*s_h/(1+R2);
    else
        Hand_est(t) = R2/(1-R1)*s_h;
        Hand_report(t) = (1-R1)/(1-R1+R2)*Hand_est(t)+R2/(1-R1+R2)*s_h;
    end

    if t<trialNum
        Hand_real(t+1) = a*Hand_real(t)+b*(0-Hand_est(t));
    end

end

end

function RMSE = ReportHand(param,pertur,ydata1,ydata2,Range)


a = param(1);
b = param(2);
R1 = param(3);
R2 = param(4);

trialNum = length(ydata1);
Hand_real = zeros(1,trialNum);
Hand_report = zeros(1,trialNum);
Hand_est = zeros(1,trialNum);

for t = 1:trialNum
    s_v = -pertur(t);
    s_h = Hand_real(t);
    if ~isnan(s_v)
        Hand_est(t) = R1*s_v+R2*s_h;
        Hand_report(t) = R1*s_v/(1+R2)+2*R2*s_h/(1+R2);
    else
        Hand_est(t) = R2/(1-R1)*s_h;
        Hand_report(t) = (1-R1)/(1-R1+R2)*Hand_est(t)+R2/(1-R1+R2)*s_h;
    end
    Hand_real(t+1) = a*Hand_real(t)+b*(0-Hand_est(t));
end
D1 = Hand_real(:);
D2 = Hand_report(:);
ydata1 = ydata1(:);
ydata2 = ydata2(:);

RMSE = sqrt(nanmean([(D1(Range)-ydata1(Range)).^2;(D2(Range)-ydata2(Range)).^2]));

end

function IsOutlier = OutlierMAD (data,crit)
MAD = 1.4862*nanmedian(abs(data-nanmedian(data)));
IsOutlier = abs(data-nanmedian(data))>crit*MAD;
end

