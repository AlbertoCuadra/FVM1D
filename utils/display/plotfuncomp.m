function plotfuncomp(x,t,wn,w,cfl,nx,model,framerate,hObject,handles)
% Plot configuration
% go to axis a
fig = handles.axesPlot;
axes(fig)
hold on
% if model(1)==1
%     tit0 = 'Rieman Problem';
% elseif model(3)==1
%     tit0 = 'Cauchy Problem';
% end
color = colours;
fig = gca;
sub1 = subplot(2,1,1,fig);
cla %Limpiar el gráfico
grid on; box on; hold on;
set(sub1,'LineWidth',1.2,'FontSize',12,'BoxStyle','full')
xlabel('$x$','FontSize',16,'interpreter','latex');
ylabel('$h(x,t)\ [m]$','FontSize',16,'interpreter','latex');
minwn1 = min(min(wn(:,:,1)));
maxwn1 = max(max(wn(:,:,1)));
minwn2 = min(min(wn(:,:,2)));
maxwn2 = max(max(wn(:,:,2)));
if (minwn1 ~= 0 || maxwn1 ~= 0)
ylim(sub1,[minwn1-0.1*abs(minwn1) maxwn1+0.1*abs(maxwn1)]);
end
sub2 = subplot(2,1,2);
cla %Limpiar el gráfico
grid on; box on; hold on;
set(sub2,'LineWidth',1.2,'FontSize',12,'BoxStyle','full')
xlabel('$x$','FontSize',16,'interpreter','latex');
ylabel('$u(x,t)\ [m/s]$','FontSize',16,'interpreter','latex');
if (minwn2 ~= 0 || maxwn2 ~= 0)
ylim(sub2,[minwn2-0.1*abs(minwn2) maxwn2+0.1*abs(maxwn2)]);
end
% mov=VideoWriter(strcat(fpath,filename),'MPEG-4');
% set(mov,'FrameRate',10);
% open(mov);
% Plot
leg{1} = 'h: Godunov';
leg{2} = 'h: Exact';
leg{3} = 'u: Godunov';
leg{4} = 'u: Exact';
for i=1:length(t)
    if get(handles.pushbuttonStop,'userdata') % stop condition
        break;
    end
	if exist('dl1')
        delete(dl1);
        delete(dl12);
        delete(dl2);
        delete(dl22);
    end
    tit1 = ['Time: ' num2str(sprintf('%.3f',round(t(i),3))) ' s',...
        ', CFL: ' num2str(cfl) '$, N_x = $' num2str(nx)];
    title(sub1,{tit1},'Interpreter','latex','FontSize',14);
	dl1 = plot(x,wn(:,i,1),'-v','Parent',sub1,'color',color(2,:),...
      'LineWidth',1.2,'MarkerSize',2);
	dl12 = plot(x,w(:,i,1),'Parent',sub1,'color',color(3,:),'LineWidth',1.2);
    dl2 = plot(x,wn(:,i,2),'-v','Parent',sub2,'color',color(1,:),...
      'LineWidth',1.2,'MarkerSize',2); 
    dl22 = plot(x,w(:,i,2),'Parent',sub2,'color',color(3,:),'LineWidth',1.2);
    legend(sub1,leg(1:2),'FontSize',12,'Location','northeastoutside','interpreter','latex')
    legend(sub2,leg(3:4),'FontSize',12,'Location','northeastoutside','interpreter','latex')
    pause(framerate)
    %     M(i)=getframe(fig);
%     filename2 = strcat(fpath,filename,'-',num2str(i));
%     saveas(fig,filename2,'epsc');
end
% writeVideo(mov,M);
% movie(M);
% close(mov);
 % enable the start button
 set(hObject,'Enable','on');
 % disable the halt button
 set(handles.pushbuttonStop,'Enable','off');