function F = Box_16(X11,X12,X13,X14,X21,X22,X23,X24,X31,X32,X33,X34,...
    X41,X42,X43,X44)

% Handling values more than 1
X11 = beyond1(X11); X12 = beyond1(X12);X13 = beyond1(X13);X14 = beyond1(X14);
X21 = beyond1(X21); X22 = beyond1(X22);X23 = beyond1(X23);X24 = beyond1(X24);
X31 = beyond1(X31); X32 = beyond1(X32);X33 = beyond1(X33);X34 = beyond1(X34);
X41 = beyond1(X41); X42 = beyond1(X42);X43 = beyond1(X43);X44 = beyond1(X44);

% X1:4 - data; Name1 - name of the distribution

group10 = [zeros(length(X11),1)',ones(length(X12),1)',zeros(length(X13),1)'+2,...
    zeros(length(X14),1)'+3];
group25 = [zeros(length(X21),1)'+4,zeros(length(X22),1)'+5,zeros(length(X23),1)'+6,...
    zeros(length(X24),1)'+7];
group50 = [zeros(length(X31),1)'+8,zeros(length(X32),1)'+9,zeros(length(X33),1)'+10,...
    zeros(length(X34),1)'+11];
group100 = [zeros(length(X41),1)'+12,zeros(length(X42),1)'+13,zeros(length(X43),1)'+14,...
    zeros(length(X44),1)'+15];


group =[group10 group25 group50 group100];


positions = [1 1.25 1.5 1.75 2.25 2.5 2.75 3 3.5 3.75 4 4.25 4.75 5 5.25 5.5];

X11 = reshape(X11,[],1);X12 = reshape(X12,[],1);X13 = reshape(X13,[],1);X14 = reshape(X14,[],1);
X21 = reshape(X21,[],1);X22 = reshape(X22,[],1);X23 = reshape(X23,[],1);X24 = reshape(X24,[],1);
X31 = reshape(X31,[],1);X32 = reshape(X32,[],1);X33 = reshape(X33,[],1);X34 = reshape(X34,[],1);
X41 = reshape(X41,[],1);X42 = reshape(X42,[],1);X43 = reshape(X43,[],1);X44 = reshape(X44,[],1);

X1 = [X11; X12; X13; X14];
X2 = [X21; X22; X23; X24];
X3 = [X31; X32; X33; X34];
X4 = [X41; X42; X43; X44];



X = [X1' X2' X3' X4'];

% difference
figure;

hB = boxplot(X,group,'Notch','marker', 'positions', positions,'Color','kkkkkkkkkkkkkkkk','symbol','');
h = findobj(gca,'Tag','Box');
for j=[1,5,9,13] % it ll take revesre order
    hP = patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5,'EdgeColor','none');
    
    % Get patch objects
    hPatch2 = findobj(hP, 'Type', 'patch');
    
    % Apply Hatch Fill
    
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.95,0.87,0.73],'FaceAlpha',1,'EdgeColor','none');
    hold on
    hh2 = hatchfill(hPatch2,'single',90,3,'none',[0.4,0.4,0.4]);
    
end

for j=[2,6,10,14]
    hP = patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5,'EdgeColor','none');
    
    % Get patch objects
    hPatch2 = findobj(hP,'Type','patch');
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.95,0.87,0.73],'FaceAlpha',1,'EdgeColor','none');
    hold on
    hh2 = hatchfill(hPatch2,'speckle',5,2,'none',[0.4,0.4,0.4]); %'single',135,2,'none','w'
end

for j=[3,7,11,15] % it ll take revesre order
    hP = patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5,'EdgeColor','none');
    
    % Get patch objects
    hPatch2 = findobj(hP,'Type','patch');
    
    % Apply Hatch Fill
    
    patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',1,'EdgeColor','none');
    hold on
    hh2 = hatchfill(hPatch2, 'single', 60, 2,'none','w');
    
end

for j=[4,8,12,16] % it ll take revesre order
    hP = patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',1,'EdgeColor','none');
end

hB = boxplot(X,group,'Notch','marker', 'positions', positions,'Color','kkkkkkkkkkkk','symbol','');
lw = 1.5;
set(gca,'fontsize',12,'FontWeight','bold','FontName', 'Times');
set(hB,'linew',lw);

set(gca,'xTick',[1.125 1.625 2.375 2.875 3.625 4.175 4.875 5.375])
set(gca,'xTickLabel',{'26','25','51','50','76','75','101','100'})

ylabel('Med(\pi)','fontsize',12,'FontWeight',...
    'bold','FontName', 'Times');

% title(strcat(Name1));
set(gca,'fontsize',12,'FontWeight','bold','FontName', 'Times')
xlabel('Sample size','fontsize',12,'FontWeight','bold','FontName', 'Times')
% xtickangle(45)

[hLegend,objh,~] = legend(findall(gca,'Tag','Box'),...
    {'C-Moment: with extreme','L-Moment: with extreme',...
    'C-Moment','L-Moment'});

set(objh,'linewidth',1.5);
% hChildren = findall(get(hLegend,'Children'), 'Type','Line');
set(gca,'fontsize',12,'FontWeight','bold','FontName', 'Times')
for j = 1:16
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines(j), 'Color', 'w','LineStyle','-','linewidth',1.5);
end
ylim([-0.2 1.5])

line([2,2],ylim,'Color','k','LineWidth',0.5,'LineStyle','--')
line([3.25,3.25],ylim,'Color','k','LineWidth',0.5,'LineStyle','--')
line([4.5,4.5],ylim,'Color','k','LineWidth',0.5,'LineStyle','--')
%line(xlim,[1,1],'Color','k','LineWidth',1.5,'LineStyle','-.')
end
function XO = beyond1(X)
XO = X;
XO(X>1)=[];
XO(XO<0)=[];
XO(isnan(XO))=[];
XO = round(XO,4);
end



