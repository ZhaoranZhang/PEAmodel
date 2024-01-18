
%% load data
clear;

data = readtable("Exp3.csv");
Deg = [10,20,40,80];

[perturb,sub,Num,propBias] = deal(data.perturbSize,data.sub,...
    data.trialNum,data.propBias);

Data1 = reshape(propBias(Num==1),4,[])';
Data2 = reshape(propBias(Num==2),4,[])';
Data3 = reshape(propBias(Num==3),4,[])';

%% simulation with PEA model

fitParam(1)= 0.2;% b
fitParam(2) = 0.97;% a
fitParam(4) = 5.048^2;% sig_g free para
fitParam(3) = 11.19^2;%FitParam(1);% sig_h
fitParam(5) = 0.309;%FitParam(2);% slope free para
fitParam(6) = 1.853;%FitParam(3);% intercept
fitParam(7) = 0.4;

deg = 1:1:100;
Tn = 100;

HandEstSimu = [];
for g = 1:length(deg)
    angle = deg(g);
     pertur(g,:) = ones(1,Tn)*angle;
    [hand_real,hand_est] = BCCTbTProp(fitParam,pertur(g,:),Tn);
   
   HandSimu(g) = nanmean(hand_real(end));
   HandEstSimu(g) = nanmean(hand_est(6));
   HandEstSimuS(g) = nanmean(hand_est(6));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HandEstSimuAll = [];
Ratios = 0.1:0.1:0.8;
Ratios = [0.05,Ratios];
for r = 1:length(Ratios)
    
    for g = 1:length(deg)
        angle = deg(g);
        pertur(g,:) = ones(1,Tn)*angle;
        [~,hand_est] = BCCTbTProp(fitParam,pertur(g,:),Tn);
        HandEstSimuAll(r,g) = hand_est(6);
    end
    HandEstSimuAll(r,:) = HandEstSimuAll(r,:).*Ratios(r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color','w','Position',[300,100,1050,275]);
Cols_ = hot(40)*0.9;
Cols = Cols_(1:40,:);
L1 = 35;
L2 = 15;

Axis = [-2,100,-4,2];
subplot(1,3,3);
box off;
set(gca,'fontsize',12,'xtick',0:20:100,'XTickLabelRotation',0,'ytick',-4:1:1);hold on;
axis(Axis);
subtitle('Experiment');hold on;
% scatter(Deg,absBlur,20,'r','filled','MarkerFacealpha',0.2,'markeredgecolor','none');hold on;
% scatter(Deg,absClear,20,'b','filled','MarkerFacealpha',0.2,'markeredgecolor','none');hold on;
% plot(Deg,AllData,'color',[Cols(1,:),0.1],'LineWidth',1.5);hold on;

plot([0,100],[0,0],'color',[0,0,0,0.3],'linestyle','--','LineWidth',1);hold on;
Leg(2) = errorbar(Deg-2,nanmean(Data2),SE(Data2),'CapSize',0,'color',Cols(25,:),...
    'LineWidth',2,'marker','diamond','markerfacecolor',Cols(25,:),...
    'MarkerEdgeColor','none','markersize',7);hold on;


Leg(3) = errorbar(Deg+2,nanmean(Data3),SE(Data3),'CapSize',0,'color',Cols(L1,:),...
    'LineWidth',2,'marker','diamond','markerfacecolor',Cols(L1,:),...
    'MarkerEdgeColor','none','markersize',7);hold on;

Leg(1) = errorbar(Deg,nanmean(Data1),SE(Data1),'CapSize',0,'color',Cols(L2,:),...
    'LineWidth',2,'marker','diamond','markerfacecolor',Cols(L2,:),...
    'MarkerEdgeColor','none','markersize',7);hold on;
legend(Leg,{'1^{st} Trial','2^{nd} Trial','3^{rd} Trial'},'FontSize',10,'box','off',...
    'location','southeast');


subplot(1,3,1);
box off;
set(gca,'fontsize',12,'xtick',0:20:100,'XTickLabelRotation',0,'ytick',-4:1:1);hold on;
axis(Axis);
subtitle('BCC Model');hold on;
plot([0,100],[0,0],'color',[0,0,0,0.3],'linestyle','--','LineWidth',1);hold on;

itv = (L1-L2)/(0.5-0.05);
Step = L1-round((Ratios-0.05)*itv);
ccols = Cols(Step,:);

for i = 1:size(HandEstSimuAll,1)
plot(deg,HandEstSimuAll(i,:),'color',[ccols(i,:),0.3],'LineWidth',1.5,'LineStyle','-');hold on;
end
plot(deg,HandEstSimuAll(6,:),'color',[Cols(L2,:),0.9],'LineWidth',3);hold on;
plot(deg,HandEstSimuAll(1,:),'color',[Cols(L1,:),0.9],'LineWidth',3);hold on;


%% Simulation with Tasy's model

% a = 0.4361;% a 0-1 (sigma_u^2/(sigma_u^2+sigma_v^2))
% b = 0.1424; % b 0-1 (sigma_u^2/(sigma_u^2+sigma_p^2))
% eta = 0.7238;% eta 0-1
% sat = 2.0756;% sat 0-100
% B = 0.3774;
% A = 0.9987;

fitParam(1)= 0.4361;% a
fitParam(2) = 0.1424;% b
fitParam(3) = 0.7238;% eta 0-1
fitParam(4) = 2.0756;% sat 0-100
fitParam(5) = 0.3774;% B
fitParam(6) = 0.9987;% A


deg = 0:1:100;
trialNum = 80;
Ratios = 0.1:0.1:0.8;
Ratios = [0.05,Ratios];
for j = 1:length(Ratios)

    R = Ratios(j);
    for i = 1:length(deg)
        s_v = ones(1,trialNum).*deg(i);
        [Hand_real, Hand_est] = TsayTbTProp(fitParam,s_v,trialNum);
        Prop(j,i) = R*Hand_est(6);
    end
end



subplot(1,3,2);
box off;
set(gca,'fontsize',12,'xtick',0:20:100,'XTickLabelRotation',0,'ytick',-4:1:1);hold on;
axis(Axis);
subtitle('PReMo Model');hold on;
plot([0,100],[0,0],'color',[0,0,0,0.3],'linestyle','--','LineWidth',1);hold on;

for i = 1:size(HandEstSimuAll,1)
plot(deg,Prop(i,:),'color',[ccols(i,:),0.3],'LineWidth',1.5,'LineStyle','-');hold on;
end

plot(deg,Prop(6,:),'color',[Cols(L2,:),0.9],'LineWidth',3);hold on;
plot(deg,Prop(1,:),'color',[Cols(L1,:),0.9],'LineWidth',3);hold on;
%%
function [Hand_real,Hand_est] = BCCTbTProp(param,pertur,trialNum)

Hand_real = zeros(1,trialNum);
Hand_est = zeros(1,trialNum);
Hand_est2 = zeros(1,trialNum);
a = param(1);
b = param(2);

s_g = 0;
sig_h = param(3);
sig_g = param(4);
slope = param(5);
intercept = param(6);


for t = 1:trialNum
    s_v = -pertur(t);
    sig_v = (slope*pertur(t)+intercept).^2;
    s_h = Hand_real(t);
    if ~isnan(s_v)
        [c_mean,c_sig_sq] = CueCombine(s_h,s_v,s_g,sig_h,sig_v,sig_g);
%         [c_mean,c_sig_sq] = CueCombine(s_h,s_v,s_g,sig_h,sig_v,0,1);
    else
        [c_mean,c_sig_sq] = CueCombine(s_h,s_g,0,sig_h,sig_g,0,1);
    end
   
    
[c_mean2,c_sig_sq2] = CueCombine(c_mean,0,0,c_sig_sq,sig_h,0,1);
Hand_est2(t) = c_mean2;
 Hand_est(t) = c_mean;
    if t<trialNum
        Hand_real(t+1) = a*Hand_real(t)+b*(0-Hand_est(t));
    end
end


end
function [Hand_real, Hand_est] = TsayTbTProp(x,s_v,trialNum)

% trialNum = 40;
% s_v  = Pertur_(g,:);
% x = FitParam;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = x(1);% a 0-1 (sigma_u^2/(sigma_u^2+sigma_v^2))
b = x(2); % b 0-1 (sigma_u^2/(sigma_u^2+sigma_p^2))
eta = x(3);% eta 0-1
sat = x(4);% sat 0-100
B = x(5);
A = x(6);
Hand_real = zeros(1,trialNum);
xper = zeros(1,trialNum);
xv = -s_v;
for i = 1:trialNum
    xp = Hand_real(i);
    if abs(xv(i))==1000
%         beta = min([sat,abs(eta*(a*xp-b*xp))]);
%         xper(i) = beta+xp*b;
        xper(i) = xp;
    elseif ~isnan(xv(i)) 
        beta = -min([sat,abs(eta*(a*xv(i)-b*xp))]);
        xper(i) = beta+xp*b;
    else
        xper(i) = b*xp;
    end
    if i<trialNum
        Hand_real(i+1) = A*Hand_real(i)+B*(0-xper(i));
    end
end

if nargout>1
    Hand_est = xper;
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