%% load data
clear;
data = readtable("Tsay_STL.csv");


[sub,perturb,Stl] = deal(data.sub,data.perturb,data.stl);
stl = reshape(Stl,9,[])';
meanEnd = nanmean(stl);
Deg = unique(perturb)';


%% PEA model
% x(1): sig_h; x(2):slope; x(3): intercept
% func = @(x)ReportHand(x,pertur,VU,handreal,handreport);

func = @(x)PEASTL(x,Deg,meanEnd);
ini = [10 0.1 0.35 2.5];%%% A, B, prop uncertainty, planning uncertainty/prop uncertainty
LB = [0  0 0.05 0.01];
UB = [50 1 1 20];

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'FiniteDifferenceStepSize',1e-5);
% options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'FiniteDifferenceStepSize',1e-5);

% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
%     'MaxFunctionEvaluations',1500,'FiniteDifferenceStepSize',1e-5);
FitParam = fmincon(func,ini,[],[],[],[],LB,UB,[],options);


s_v = linspace(0,100,100);

sig_i = FitParam(1);
B = FitParam(2);
slope = FitParam(3);
intercept = FitParam(4);
% x(1): sig_i;  integrated sigma g and h 
% x(2): learning rate
% x(3): slope; 0.51
% x(4): intercept 5.20
sig_v = slope*s_v+intercept;

for i  = 1:length(s_v)
    BCC_model_STL(i) = s_v(i)*(sig_i^2/(sig_v(i)^2+sig_i^2))*B;
    BCC_model_hand_est(i) = s_v(i)*(sig_i^2/(sig_v(i)^2+sig_i^2));
%     hand_final(i) = sig_h^2/sig_v(i)^2*s_v(i);
end


for j = 1:length(Deg)
    hand_est(j) = Deg(j)*(sig_i^2/((slope*Deg(j)+intercept)^2+sig_i^2));
    simuData(j) = hand_est(j)*B;
%     simuData(j) = sig_h^2/(FitParam(2)*Deg(j)+FitParam(3))^2*Deg(j);
end
Res = simuData-meanEnd;
SSres = nansum(Res.^2);
SStot = nansum((meanEnd-mean(meanEnd)).^2);
Rsquared_bcc = 1-SSres/SStot;
% ResVR = tempVerbal_fit-tempVerbal';
Res = Res(:);
Res = Res(~isnan(Res));
NumofParam = size(FitParam,2);

pd = fitdist(Res,'normal');
[aic_BCC,~] = aicbic(-pd.NLogL,NumofParam,1);

%% PReMo model

% Tsay model
func = @(x)TsaySTL(x,Deg,meanEnd);
ini = [0.05 0.4];%%% A, B, prop uncertainty, planning uncertainty/prop uncertainty
LB = [0 0];
% UB = [0.5 1];
    UB = [10 10];

FitParam = fmincon(func,ini,[],[],[],[],LB,UB,[],options);

slope = FitParam(1);
sat = FitParam(2);
s_v = linspace(0,100,100);
for i  = 1:length(s_v)
    if s_v(i)*slope<sat
        Tsay_model(i) = s_v(i)*slope;
    else
        Tsay_model(i) = sat;
    end
end
simuData = [];
for j = 1:length(Deg)
    if Deg(j)*slope<sat
        simuData(j) = Deg(j)*slope;
    else
        simuData(j) = sat;
    end
end

% SSres = nansum(nansum((Hand_real(:,FitRange(1):FitRange(2))-Data(:,FitRange(1):FitRange(2))).^2));
% SStot = nansum(nansum((Data-nanmean(Data(:))).^2));
% Rsquared = 1-SSres/SStot;
Res = simuData-meanEnd;

SSres = nansum(Res.^2);
SStot = nansum((meanEnd-mean(meanEnd)).^2);
Rsquared_tsay = 1-SSres/SStot;
% ResVR = tempVerbal_fit-tempVerbal';
Res = Res(:);
Res = Res(~isnan(Res));
NumofParam = size(FitParam,2);

pd = fitdist(Res,'normal');
[aic_Tsay,~] = aicbic(-pd.NLogL,NumofParam,1);


%% Causal inference model


func = @(x)RelSTL(x,Deg,meanEnd);

ini=[0.7674   80.6253    0.0036];%initiate point
LB=[0 0 0 ];%lowerbnd
UB=[inf inf inf];%upperbnd

FitParam = fmincon(func,ini,[],[],[],[],LB,UB,[],options);

s_v = linspace(0,100,100);
e = s_v;
x(1) = FitParam(1);
x(2) = FitParam(2);
x(3) = FitParam(3);
p = x(1).*(normpdf(e,0,x(2))./(normpdf(e,0,x(2))+x(3)));

Rel_model = e.*p;

simuData = [];
for j = 1:length(Deg)
p = x(1).*(normpdf(Deg(j),0,x(2))./(normpdf(Deg(j),0,x(2))+x(3)));

simuData(j) = Deg(j).*p;
end




Res = simuData-meanEnd;
SSres = nansum(Res.^2);
SStot = nansum((meanEnd-mean(meanEnd)).^2);
Rsquared_rel = 1-SSres/SStot;
% ResVR = tempVerbal_fit-tempVerbal';
Res = Res(:);
Res = Res(~isnan(Res));
NumofParam = size(FitParam,2);

pd = fitdist(Res,'normal');
[aic_Rel,~] = aicbic(-pd.NLogL,NumofParam,1);

figure('color','w','Position',[100,100,1100,300]);
% sgtitle('Single Trial Learning');
Axis = [-2,100,-0.01,0.8];
subplot(1,3,1);
box off;
set(gca,'fontsize',12,'xtick',0:20:100);hold on;
axis(Axis);
b = plot(s_v,BCC_model_STL,'color',[247, 54, 109]./255,'LineWidth',3);hold on;
 errorbar(Deg,meanEnd,SE(stl),'o','color',[139, 192, 214]./255, ...
    'markerfacecolor',[139, 192, 214]./255,...
        'MarkerEdgeColor','none','linewidth',2,'capsize',0);hold on;
   
%     xticks(tempPertur(1+(uncer_i-1)*10:uncer_i*10));
    xlabel('Perturbation (deg)');
    ylabel('\Delta Hand angle (deg/trial)');
subtitle('BCC model');
% legend(b,{sprintf('AIC = %.2f',aic_BCC)},'Location','northeast',...
%     'box','off','FontSize',10);
text(70,0.7,sprintf('\n{\\itR}^2 = %.2f\nAIC = %.2f',Rsquared_bcc,aic_BCC),'FontSize',10);hold on;

subplot(1,3,2);
box off;
set(gca,'fontsize',12,'xtick',0:20:100);hold on;
axis(Axis);


 b = plot(s_v,Tsay_model,'color',[247, 54, 109]./255,'LineWidth',3);hold on;  
 errorbar(Deg,meanEnd,SE(stl),'o','color',[139, 192, 214]./255, ...
    'markerfacecolor',[139, 192, 214]./255,...
        'MarkerEdgeColor','none','linewidth',2,'capsize',0);hold on;
%     xticks(tempPertur(1+(uncer_i-1)*10:uncer_i*10));
    xlabel('Perturbation (deg)');
    % ylabel('Final angle (deg)');
subtitle('PReMo model');
% legend(b,{sprintf('\nR^2 = %.2f\nAIC = %.2f',Rsquared_tsay,aic_Tsay)},'Location','northeast',...
%     'box','off','FontSize',10);
text(70,0.7,sprintf('\n{\\itR}^2 = %.2f\nAIC = %.2f',Rsquared_tsay,aic_Tsay),'FontSize',10);hold on;

subplot(1,3,3);
box off;
set(gca,'fontsize',12,'xtick',0:20:100);hold on;
axis(Axis);

b = plot(s_v,Rel_model,'color',[247, 54, 109]./255,'LineWidth',3);hold on;
 errorbar(Deg,meanEnd,SE(stl),'o','color',[139, 192, 214]./255, ...
    'markerfacecolor',[139, 192, 214]./255,...
        'MarkerEdgeColor','none','linewidth',2,'capsize',0);hold on;
%     xticks(tempPertur(1+(uncer_i-1)*10:uncer_i*10));
    xlabel('Perturbation (deg)');
    % ylabel('Final angle (deg)');
subtitle('Causal inference model');
% legend(b,{sprintf('AIC = %.2f',aic_Rel)},'Location','northeast',...
%     'box','off','FontSize',10);
text(70,0.7,sprintf('\n{\\itR}^2 = %.2f\nAIC = %.2f',Rsquared_rel,aic_Rel),'FontSize',10);hold on;
%% functions

function RMSE = PEASTL(x,s_v,data)


% x(1): sig_i;  integrated sigma g and h 
% x(2): learning rate
% x(3): slope;
% x(4): intercept

sig_i = x(1);
B = x(2);
slope = x(3);
intercept = x(4);

sig_v = slope*s_v+intercept;

for i  = 1:length(s_v)
    hand_est(i) = s_v(i)*(sig_i^2/(sig_v(i)^2+sig_i^2));
    STL(i) = hand_est(i)*B;
%     hand_final(i) = sig_h^2/sig_v(i)^2*s_v(i);
end

RMSE = sqrt(mean((data-STL).^2));

end

function RMSE = TsaySTL(x,s_v,data)


% x(1):slope; x(2): saturate


slope = x(1);
sat = x(2);


for i  = 1:length(s_v)
    if s_v(i)*slope<sat
        hand_final(i) = s_v(i)*slope;
    else
        hand_final(i) = sat;
    end
end

RMSE = sqrt(mean((data-hand_final).^2));

end

function RMSE = RelSTL(x,s_v,data)


e = s_v;
p = x(1).*(normpdf(e,0,x(2))./(normpdf(e,0,x(2))+x(3)));

hand_final = e./(1./p-1);


RMSE = sqrt(mean((data-hand_final).^2));

end








