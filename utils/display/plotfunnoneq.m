function plotfunnoneq(x,t,framerate,hObject,handles,...
    model,wl,wr,xb,tb,s,w0,dtdx,nx)
% Dimensionalization
w_exact = zeros(length(x),length(t));
wnq = zeros(length(x),length(t));
wnnc = zeros(length(x),length(t));
wncen = zeros(length(x),length(t));
wnqr = zeros(length(x),length(t));
wngod = zeros(length(x),length(t));
% Inizialization
% w0 is implicit in the first result of w_exact
w_exact(:,1) = exact(x,t(1),xb,tb,wl,wr,s); 
waq(1,:)= w_exact(:,1);% Q-scheme (Roe=van leer for Burgers equation) 
wanc(1,:)= w_exact(:,1);% Non conservative scheme
wacen(1,:)= w_exact(:,1);% Centred scheme
waqr(1,:)= w_exact(:,1);% Q-scheme Local Lax Friedrich Regularization
wagod(1,:)= w_exact(:,1);% Godunov method
% Plot configuration
% go to axis a
fig = handles.axesPlotA2;
axes(fig)
fig = gca;
cla %Limpiar el gráfico
hold on
color = colours;
set(fig,'LineWidth',1.2,'FontSize',12,'BoxStyle','full')

grid on; box on; hold on;

descr = {'Description of the problem:';
    'Burgers equation:';
    ' ';
    '\hspace{0.3cm}$\frac{\partial w}{\partial t} + \frac{\partial}{\partial x} \left( \frac{w^2}{2}\right) =0.$';
    ' '};
% text(3.2,0.3,descr,'interpreter','latex','FontSize',18)
wmin=min(wl,wr)-0.2; %lower limit for the y axi
wmax=max(wl,wr)+0.2; %upper limit for the y axi
ylim([wmin,wmax]);
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$w(x,t)$','FontSize',20,'interpreter','latex');
% mov=VideoWriter(strcat(fpath,filename),'MPEG-4');
% set(mov,'FrameRate',framerate);
% open(mov);
% Plot
cont=0;
for i=1:length(model)
    if model(i)==1
        cont=cont+1;
    end
end
i=1;
while i<=cont
    if model(1)==1
        leg{i} = ['Exact'];
        i=i+1;
    end
    if model(2)==1
        leg{i} = ['Q-Scheme'];
        i=i+1;
    end
    if model(3)==1
        leg{i} = ['Non-conservative']; 
        i=i+1;
    end
    if model(4)==1
        leg{i} = ['Centred'];
        i=i+1;
    end 
    if model(5)==1
        leg{i} = ['Q-scheme LLF'];
        i=i+1;
    end
    if model(6)==1
        leg{i} = ['Godunov'];
        i=i+1;
    end
end
for i=1:length(t)
    w_exact(:,i) = exact(x,t(i),xb,tb,wl,wr,s);
    wnq(:,i)=qscheme_btbc(waq,dtdx,nx);
    wnnc(:,i)=ncon_btbc(wanc,dtdx,nx);
    wncen(:,i)=centred(wacen,dtdx,nx);
    wnqr(:,i)=qscheme_llfr_btbc(waqr,dtdx,nx);
    wngod(:,i)=god_btbc(wagod,dtdx,nx);
    waq(1,:)=wnq(:,i);
    wanc(1,:)=wnnc(:,i);
    wacen(1,:)=wncen(:,i);
    waqr(1,:)=wnqr(:,i);
    wagod(1,:)=wngod(:,i);
end
for i=1:length(t)
    if get(handles.pushbuttonStopA2,'userdata') % stop condition
        break;
    end
    if model(1)==1
        if exist('dl1')
            delete(dl1);
        end
%     w_exact(:,i) = exact(x,t(i),xb,tb,wl,wr,s);
    dl1 = plot(x,w_exact(:,i),'color',color(2,:),'LineWidth',3,...
         'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    % Q-scheme 
    if model(2)==1
        if exist('dl2')
            delete(dl2);
        end
%     wnq=qscheme_btbc(waq,dtdx,nx);
    dl2 = plot(x,wnq(:,i),'--s','color',color(1,:),'LineWidth',1.4,...
         'MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    % Non conservative scheme 
    if model(3)==1
        if exist('dl3')
            delete(dl3);
        end
%     wnnc=ncon_btbc(wanc,dtdx,nx);
    dl3 = plot(x,wnnc(:,i),'-.o','color',color(3,:),'LineWidth',1.4,...
         'MarkerFaceColor',color(3,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    % Centred scheme 
    if model(4)==1
        if exist('dl4')
            delete(dl4);
        end
%     wncen=centred(wacen,dtdx,nx);
    dl4 = plot(x,wncen(:,i),'-v','color',color(4,:),'LineWidth',1.4,...
         'MarkerFaceColor',color(4,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    % Q-scheme with Roe linearization
    if model(5)==1
        if exist('dl5')
            delete(dl5);
        end
%     wnq=qscheme_btbc(waq,dtdx,nx);
    dl5 = plot(x,wnqr(:,i),'-d','color',color(5,:),'LineWidth',1.4,...
         'MarkerFaceColor',color(5,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    if model(6)==1
        if exist('dl6')
            delete(dl6);
        end
%     wnq=qscheme_btbc(waq,dtdx,nx);
    dl6 = plot(x,wngod(:,i),'-.<','color',color(6,:),'LineWidth',1.4,...
         'MarkerFaceColor',color(6,:),'MarkerEdgeColor','black',...
         'MarkerSize',4);
    end
    tit = ['Time: ',num2str(sprintf('%.3f',round(t(i),3))),' s'];
    title([{tit}],'Interpreter','latex','FontSize',18)
    legend(leg,'FontSize',18,'Location','northeastoutside','interpreter','latex') 
    pause(framerate)
%   M(i)=getframe(fig);
%   filename2 = strcat(fpath,filename,'-',num2str(i));
%   saveas(fig,filename2,'epsc');
%     if model(2)==1
%         waq=wnq;
%     elseif model(3)==1
%         wanc=wnnc;
%     elseif model(4)==1
%         wacen=wncen;
%     end
end
% writeVideo(mov,M);
% movie(M);
% close(mov);
 % enable the start button
 set(hObject,'Enable','on');
 % disable the halt button
 set(handles.pushbuttonStopA2,'Enable','off');