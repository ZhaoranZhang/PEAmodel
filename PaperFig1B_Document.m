
%% load data
clear;
Fig1B = readtable('Fig1B.csv');

%% fitting with PEA mdoel
T = {'Kim,2018 exp1','Kim,2018 exp2','Morehead,2017'};

% R1 = x(1);% sig_h/inter
% R2 = x(2);% slope/inter 

ini = [10/1.8 0.3/1.8];
LB =  [1    0];
UB =  [5    1];
PLB = [1.5  0.02];
PUB = [3    0.4];
s_v = linspace(0,100,100);


Cols = turbo(4);
figure('color','w','Position',[100,100,400,300]);

for i = 1:3
    Deg = Fig1B.perturbSize(Fig1B.study==i)';
    meanEnd = Fig1B.extent(Fig1B.study==i)';
    Title = T{i};

    func = @(x)BCCEnd(x,Deg,meanEnd);
    FitParam = fmincon(func,ini,[],[],[],[],LB,UB);

   
    R1 = FitParam(1);
    R2 = FitParam(2);


    for j  = 1:length(s_v)

        BCC_model(j) = (R1/(1+R2*s_v(j))).^2*s_v(j);
    end


    simuData = [];
    for j = 1:length(Deg)
        simuData(j) = (R1/(1+R2*Deg(j))).^2*Deg(j);
    end

    SSres = sum((simuData-meanEnd).^2);
    SStot = sum((meanEnd-nanmean(meanEnd)).^2);
    Rsquared(i) = 1-SSres/SStot;
    RMS(i) = sqrt(nanmean((simuData-meanEnd).^2));



    plot(s_v,BCC_model,'color',[Cols(i,:)*0.8,0.5],'LineWidth',3);hold on;
    plot(Deg,meanEnd,'-o','markerfacecolor',Cols(i,:),...
        'MarkerEdgeColor','k','linewidth',1,'color',Cols(i,:));hold on;


end

box off;
set(gca,'fontsize',16,'xtick',0:20:100);hold on;

axis([-5,100,0,45]);
xlabel('Perturbation (deg)');
ylabel('Hand angle (deg)');
%%
function RMSE = BCCEnd(x,s_v,data)

R1 = x(1);% sig_h/inter
R2 = x(2);% slope/inter 


for i  = 1:length(s_v)
    hand_final(i) = (R1/(1+R2*s_v(i))).^2*s_v(i);
end

RMSE = sqrt(mean((data-hand_final).^2));

end



