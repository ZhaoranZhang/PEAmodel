%% experiment data
clear;
data = readtable("Exp4data.csv");
Deg = [-64, -16, -4,4,16,64];

[sub,day,perturb,blur,stl] = deal(data.sub,data.day,data.perturb, ...
    data.blur,data.stl);
subnum = unique(sub);
perturb = round(perturb);
% blur: 1 is clear 2 is blur
for pert = 1:3
    for s = 1:length(subnum)
        % all 3 days
        temp1 = stl(sub==subnum(s)&perturb==Deg(4-pert)&blur==1);
        temp2 = stl(sub==subnum(s)&perturb==Deg(3+pert)&blur==1);
        absClear(s,pert) = (nanmean(temp1)-nanmean(temp2))./2;
        temp1 = stl(sub==subnum(s)&perturb==Deg(4-pert)&blur==2);
        temp2 = stl(sub==subnum(s)&perturb==Deg(3+pert)&blur==2);
        absBlur(s,pert) = (nanmean(temp1)-nanmean(temp2))./2;
        % for d = 1:3
        %     temp = stl(sub==subnum(s)&abs(perturb)==Deg(pert)&day==d);
        % end
    end
end

%
Deg = [4,16,64];
figure('color','w','Position',[300,100,380,300]);
colr = [205, 24, 24]./255;
colb = [25, 167, 206]./255;

Axis = [0,70,0,2];

box off;
set(gca,'fontsize',16,'xtick',0:20:100);hold on;
axis(Axis);
subtitle('Experiment');hold on;
% scatter(Deg,absBlur,20,'r','filled','MarkerFacealpha',0.2,'markeredgecolor','none');hold on;
% scatter(Deg,absClear,20,'b','filled','MarkerFacealpha',0.2,'markeredgecolor','none');hold on;
plot(Deg,absBlur,'color',[colr,0.1],'LineWidth',2);hold on;
plot(Deg,absClear,'color',[colb,0.1],'LineWidth',2);hold on;
errorbar(Deg,nanmean(absBlur),SE(absBlur),'-d','linewidth',2,...
    'CapSize',0,'markerfacecolor',colr,'color',colr);hold on;
    errorbar(Deg,nanmean(absClear),SE(absClear),'-d','linewidth',2,...
    'CapSize',0,'markerfacecolor',colb,'color',colb);hold on;
%% simulation with PEA model


Axis = [-2,70,-0.1,1.5];
colr = [205, 24, 24]./255;
colb = [25, 167, 206]./255;

Step = 12;
colmap = [linspace(1,colr(1),Step)' linspace(1,colr(2),Step)' linspace(1,colr(3),Step)';
    linspace(colr(1),0,Step)' linspace(colr(2),0,Step)' linspace(colr(3),0,Step)'];
colmap(10,:) = [];
% colmap = flip(colmap,2);

Deg = [4,16,64];

s_v = linspace(0,100,100);

ratios = linspace(1.1,3,10);
Ratio = ratios(4);
% fitParam(1)= 0.2;%Constant(1);% a
% fitParam(2) = 0.97;%Constant(2);% b
% fitParam(4) = 5.048^2;%Constant(3);% sig_g free para
% fitParam(3) = 11.19^2;%FitParam(1);% sig_h
% fitParam(5) = 0.309;%FitParam(2);% slope free para
% fitParam(6) = 1.853;%FitParam(3);% intercept
% fitParam(7) = 0.4;
% sig_hg = 5.048^2*11.19^2/(11.19^2+5.048^2); 21.17

sig_i = sqrt(21.17);
slope = 0.309;
intercept = 1.853;
B = 0.21;
sig_v = slope*s_v+intercept;
% sig_v_Blur = (slope*s_v+intercept).*Ratio;
% 
for i  = 1:length(s_v)
    BCC_model_STL(i) = B*s_v(i)*(sig_i^2/(sig_v(i)^2+sig_i^2));
    % BCC_model_STL_Blur(i) = B*s_v(i)*(sig_i^2/(sig_v_Blur(i)^2+sig_i^2));
end

BCC_Blur = [];
for rs = 1:length(ratios)
    Ratio = ratios(rs);
    sig_v_Blur = (slope*s_v+intercept).*Ratio;
    for i  = 1:length(s_v)
        BCC_Blur(rs,i) = B*s_v(i)*(sig_i^2/(sig_v_Blur(i)^2+sig_i^2));
    end
end


figure('color','w','Position',[300,100,350,300]);

% subplot(1,3,1);
box off;
set(gca,'fontsize',16,'xtick',0:20:100);hold on;
axis(Axis);
subtitle('BCC Model');hold on;
for j = 1:size(BCC_Blur,1)
plot(s_v,BCC_Blur(j,:),'color',[colmap(2+2*(j-1),:),0.9],'LineWidth',2);hold on;
end

b = plot(s_v,BCC_model_STL,'color',[colb,0.9],'LineWidth',3);hold on;
% plot(s_v,BCC_model_STL_Blur,'color',[colr,0.9],'LineWidth',3);hold on;


%% simulation Tsay's model
colr = [205, 24, 24]./255;
colb = [25, 167, 206]./255;

% a = 0.4361;% a 0-1 (sigma_u^2/(sigma_u^2+sigma_v^2))
% b = 0.1424; % b 0-1 (sigma_u^2/(sigma_u^2+sigma_p^2))
% eta = 0.7238;% eta 0-1
% sat = 2.0756;% sat 0-100
% B =0.3774;
% A =0.9987;
Ratio = ratios(4);
a_blur = 1/(1+(Ratio)^2*(1/0.4361-1));

s_v = linspace(0,100,100);
FitParams = [0.4361,0.1424,0.7238,2.0756,0.3774];
FitParams_blur = [a_blur ,0.1424,0.7238,2.0756,0.3774];

for g = 1:length(s_v)
    Tsay_STL(g) = TsaySTL_VU(FitParams,s_v(g));
    Tsay_STL_Blur(g) = TsaySTL_VU(FitParams_blur,s_v(g));
end



Tsay_Blur = [];
for rs = 1:length(ratios)
    Ratio = ratios(rs);
    a_blur = 1/(1+(Ratio)^2*(1/0.4361-1));
    FitParams_blur = [a_blur ,0.1424,0.7238,2.0756,0.3774];

    for i  = 1:length(s_v)
        Tsay_Blur(rs,i) = TsaySTL_VU(FitParams_blur,s_v(i));
    end
end



figure('color','w','Position',[300,100,350,300]);

box off;
set(gca,'fontsize',16,'xtick',0:20:100);hold on;
axis(Axis);
subtitle('PReMo Model');hold on;

for j = 1:size(BCC_Blur,1)
plot(s_v,Tsay_Blur(j,:),'color',[colmap(1+2*(j-1),:),0.9],'LineWidth',2);hold on;
end

b(2) = plot(s_v,Tsay_Blur(6,:),'color',[colmap(1+2*(6-1),:),0.9],'LineWidth',2);hold on;
b(1) = plot(s_v,Tsay_STL,'color',[colb,0.9],'LineWidth',3);hold on;


legend(b,{'Clear target','Blurred target'},'fontsize',10,'location','northwest','box','off');

% plot(Deg,simdata_Tsay_STL_Blur,'-d','color',colr,'linewidth',2,'markerfacecolor',colr);hold on;
% plot(Deg,simdata_Tsay_STL,'-d','color',colb,'linewidth',2,'markerfacecolor',colb);hold on;
% xlabel('Perturbation (deg)');
% ylabel('Learning rate (deg/trial)');
%% simulation with causal inference model

colr = [205, 24, 24]./255;
colb = [25, 167, 206]./255;

% S = 201.7668 ;
% Sigma = 17.6013;
% C = 31.4169;
% A = 0.9486;
% B  = 0.8544;


s_v = linspace(0,100,100);

FitParams = [201.7668 17.6013 31.4169 0.9486 0.8544];

sp = 40;% assume sigma_prop = 40 sigma_prop is 2 times of sigma_v;
si = 17.6013;
vp = sqrt(sp^2*si^2/(sp^2-si^2));
Ratio = 2;
% Sigma_blur = sqrt(sp^2/(1+(1/Ratio^2)*(sp^2-si^2)/si^2));
Sigma_blur = sqrt(5*Ratio^2/(Ratio^2+4))*si;
FitParams_blur = [201.7668 Sigma_blur 31.4169 0.9486 0.8544];


for g = 1:length(s_v)
    Rel_STL(g) = RelSTL_VU(FitParams,s_v(g));
    Rel_STL_Blur(g) = RelSTL_VU(FitParams_blur,s_v(g));
end


Rel_Blur = [];



for rs = 1:length(ratios)
    Ratio = ratios(rs);
   % Sigma_blur = sqrt(sp^2/(1+(1/Ratio^2)*(sp^2-si^2)/si^2));
   Sigma_blur = sqrt(5*Ratio^2/(Ratio^2+4))*si;
    FitParams_blur = [201.7668 Sigma_blur 31.4169 0.9486 0.8544];

    for i  = 1:length(s_v)
        Rel_Blur(rs,i) = RelSTL_VU(FitParams_blur,s_v(i));
    end
end


figure('color','w','Position',[300,100,350,300]);
box off;
set(gca,'fontsize',16,'xtick',0:20:100);hold on;
axis(Axis);
subtitle('Causal Inference Model');hold on;

for j = 1:size(Rel_Blur,1)
plot(s_v,Rel_Blur(j,:),'color',[colmap(1+2*(j-1),:),0.9],'LineWidth',2);hold on;
end

b(2) = plot(s_v,Rel_STL,'color',[colb,0.9],'LineWidth',3);hold on;
% b(1) = plot(s_v,Rel_STL_Blur,'color',[colr,0.9],'LineWidth',2);hold on;
% b(2) = plot(Deg,simdata_Rel_STL_Blur,'-d','color',colr,'linewidth',2,'markerfacecolor',colr);hold on;
% b(1) = plot(Deg,simdata_Rel_STL,'-d','color',colb,'linewidth',2,'markerfacecolor',colb);hold on;
% xlabel('Perturbation (deg)');
% ylabel('Hand angle (deg/trial)');


%%
function [Hand_stl] = TsaySTL_VU(x,s_v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = x(1);% a 0-1 (sigma_u^2/(sigma_u^2+sigma_v^2))
b = x(2); % b 0-1 (sigma_u^2/(sigma_u^2+sigma_p^2))
eta = x(3);% eta 0-1
sat = x(4);% sat 0-100
B = x(5);



xv = -s_v;

xp = 0;

beta = min([sat,abs(eta*(a*xv-b*xp))]);
xper = -beta+xp*b;

Hand_stl = B*(0-xper);




end
function Hand_stl = RelSTL_VU(x,s_v)

S = x(1);
Sigma = x(2);
C = x(3);
A = x(4);
B  = x(5);

xv = -s_v;
e = xv;
p = S.*(normpdf(e,0,Sigma)./(normpdf(e,0,Sigma)+C));
xper = e*p;

Hand_stl = B*(0-xper);


end





































