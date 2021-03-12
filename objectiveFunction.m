function scarto = objectiveFunction(reactionParameters,Data)

reactionParameters=abs(reactionParameters);

vTime=Data(:,1);
vB=Data(:,2);
vD=Data(:,3);
% vB=Data(:,4);
% vC=Data(:,5);
% vD=Data(:,6);


y0 = [vB(1) vD(1) 1 1] ;
[Time,sol] = ode45(@(t,y) odeSystem(t,y,reactionParameters), vTime, y0); % da inserire la function (modello da calibrare)
solB=sol(:,1);
solD=sol(:,2);
% solB=sol(:,3);
% solC=sol(:,4);
% solD=sol(:,5);

% size_vTime = size(vTime)
% size_sold = size(solD)
% size_vd = size(vD)

figure(6);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
subplot(2,1,1)
plot(Time,solD,'k',vTime,vD,'ok' , 'LineWidth', 2);
%plot(vTime,solD,'k', 'LineWidth', 2);
subplot(2,1,2)
plot(Time,solB,'r', vTime,vB,'or', 'LineWidth', 2);

shg

%check length of vTime and Time

if length(vTime)  == length(Time)

method=3;
switch method
    case 1
        scarto=(1/length(vTime))*sum((solD-vD).^2+(solD-vD).^2+(solB-vB).^2)/2;
    case 2
        scarto=(1/length(vTime))*sum((vTAC.*(vTAC-solD)).^2+(vB.*(vB-solB)).^2)/2;
    case 3
        x1D=vD(2:length(vD));
        x2D=[x1D;x1D(end)];
        dx2D=x2D-vD;
        x1B=vB(2:length(vB));
        x2B=[x1B;x1B(end)];
        dx2B=x2B-vB;
        
        scarto=(1/length(vTime))*sum((dx2D.*(vD-solD)).^2 + (dx2B.*(vB-solB)).^2);
end

else
    
    scarto =1;
            
end

end




