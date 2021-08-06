function plotfun(x,t,w,framerate,hObject,handles)
% Plot configuration
% go to axis a
fig = handles.axesPlot;
axes(fig)
hold on
color = colours;
fig = gca;
sub1 = subplot(2,1,1,fig);
cla %Limpiar el gráfico
grid on; box on; hold on;
set(sub1,'LineWidth',1.2,'FontSize',12,'BoxStyle','full')
xlabel('$x$','FontSize',16,'interpreter','latex');
ylabel('$h(x,t)\ [m]$','FontSize',16,'interpreter','latex');
ylim([min(min(w(:,:,1)))-0.2 max(max(w(:,:,1)))+0.2]);
sub2 = subplot(2,1,2);
cla %Limpiar el gráfico
grid on; box on; hold on;
set(sub2,'LineWidth',1.2,'FontSize',12,'BoxStyle','full')
xlabel('$x$','FontSize',16,'interpreter','latex');
ylabel('$u(x,t)\ [m/s]$','FontSize',16,'interpreter','latex');
ylim(sub2,[min(min(w(:,:,2)))-0.2 max(max(w(:,:,2)))+0.2]);
% mov=VideoWriter(strcat(fpath,filename),'MPEG-4');
% set(mov,'FrameRate',framerate);
% open(mov);
% Plot
for i=1:length(t)
    if get(handles.pushbuttonStop,'userdata') % stop condition
        break;
    end
    if exist('dl1')
        delete(dl1);
        delete(dl2);
    end 
    dl1 = plot(x,w(:,i,1),'Parent',sub1,'color',color(2,:),'LineWidth',1.4);
    tit = ['Time: ',num2str(sprintf('%.3f',round(t(i),3))),' s'];
    title({tit},'Parent',sub1,'Interpreter','latex','FontSize',18)
    dl2 = plot(x,w(:,i,2),'Parent',sub2,'color',color(1,:),'LineWidth',1.4); 
    pause(framerate)
%   M(i)=getframe(fig);
%   filename2 = strcat(fpath,filename,'-',num2str(i));
%   saveas(fig,filename2,'epsc');
end
% writeVideo(mov,M);
% movie(M);
% close(mov);
% enable the start button
set(hObject,'Enable','on');
% disable the halt button
set(handles.pushbuttonStop,'Enable','off');