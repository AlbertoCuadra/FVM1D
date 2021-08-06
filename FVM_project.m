function varargout = FVM_project(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FVM_project_OpeningFcn, ...
                   'gui_OutputFcn',  @FVM_project_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% --- Executes just before FVM_project is made visible.
function FVM_project_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
logo1 = imread('logo_m2i.jpg');
imshow(logo1,'Parent', handles.axesLogo1); 
imshow(logo1,'Parent', handles.axesLogoA2); 
%%%%%%%%%%%%%%%%%%%% INICIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise tabs
handles.tabManager = TabManager(hObject);
% Set-up a selection changed function on the create tab groups
tabGroups = handles.tabManager.TabGroups;
for tgi=1:length(tabGroups)
    set(tabGroups(tgi),'SelectionChangedFcn',@tabChangedCB)
end
% TEXT annotations need an axes as parent so create an invisible axes which
% is as big as the figure
% Find all static text UICONTROLS whose 'Tag' starts with latex_
lbls = findobj(hObject,'-regexp','tag','latex_*');
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      f = get(l,'FontSize');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex');
      if i==1
        text(handles.Description,p(1),p(2),s,'interpreter','latex',...
            'FontSize',f)
      elseif i==2
        text(handles.Description2,p(1),p(2),s,'interpreter','latex',...
            'FontSize',f)
      elseif i==3
        text(handles.Description3,p(1),p(2),s,'interpreter','latex',...
            'FontSize',f)
      end
end
axes(handles.axesPlotA2);
cla reset;
%%%%%%%%%%% Clear axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(handles.axesPlot);
% cla reset
% axes(handles.axesSeccion);
% cla reset
% axes(handles.axesTemperatura);
% cla reset
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes FVM_project wait for user response (see UIRESUME)
% uiwait(handles.figureMain);

% Called when a user clicks on a tab
function tabChangedCB(src, eventdata)
% disp(['Changing tab from ' eventdata.OldValue.Title ' to ' eventdata.NewValue.Title ] );

% --- Outputs from this function are returned to the command line.
function varargout = FVM_project_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
close;
% --- Executes on button press in pushbuttonRepresentar.
function pushbuttonPlot_Callback(hObject, eventdata, handles)
% disable the start button
set(hObject,'Enable','off');
set(handles.pushbuttonStop,'Enable','on');
set(handles.pushbuttonStop,'userdata',0);
% set(handles.pushbuttonPlot,'userdata',1);
guidata(hObject,handles);
%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver
model1 = get(handles.checkbox1,'Value');
model2 = get(handles.checkbox2,'Value');
model3 = get(handles.checkbox3,'Value');
model4 = get(handles.checkbox4,'Value');
model = [model1, model2, model3, model4];
% Mesh
a = str2num(get(handles.editxa,'string'));
b = str2num(get(handles.editxb,'string'));
nx = str2num(get(handles.editnx,'string'));
t0 = str2num(get(handles.editt0,'string'));
tmax = str2num(get(handles.edittf,'string'));
cfl = str2num(get(handles.editcfl,'string'));
% Matrix A
a11 = str2num(get(handles.edita11,'string'));
a12 = str2num(get(handles.edita12,'string'));
a21 = str2num(get(handles.edita21,'string'));
a22 = str2num(get(handles.edita22,'string'));
% Initial condition
wl1 = str2num(get(handles.editwl1,'string'));
wl2 = str2num(get(handles.editwl2,'string'));
wr1 = str2num(get(handles.editwr1,'string'));
wr2 = str2num(get(handles.editwr2,'string'));
% Pause
framerate = str2num(get(handles.editPause,'string'));
% Construction data
wl = [wl1;wl2];
wr = [wr1;wr2];
deltax = abs((b-a))/(nx-1);
x = a:deltax:b; 
A = [a11 a12;a21 a22];
[P_aux,D] = eig(A); 
lambda= diag(D);
deltat=cfl*deltax/max(abs(lambda)); % Deltat from Courant number
dtdx= deltat/deltax;
nt=tmax/deltat;
t = t0:deltat:tmax;
% Text editing
set(handles.textdx,'String',round(deltax,3));
set(handles.textdt,'String',round(deltat,3));
set(handles.textnt,'String',round(nt));
%%%%%%%%%%%% Solver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(model)>0
    % SOLVING THE SYSTEM - RIEMANN PROBLEM
    if model(1)==1 || model(2)==1
        % Matrix of the hyperbolic system
        A = [a11 a12;a21 a22];
        % Eigenvectors and eigenvalues 
        [P_aux,D] = eig(A); 
        % Ordering eigenvalues and eigenvectors 
        % lambda1 < lambda2 < ... < lambdap
        % x-lambda1*t> ... > 0 > ... > x-lambdap*t
        [lambda,ind] = sort(diag(D));
        P=P_aux(:,ind);
        alphal = P\wl; %inv(P)*wl; pinv(P)*wl; P\wl;
        alphar = P\wr; %inv(P)*wr; pinv(P)*wr; P\wr;
        % Summoning w_riemann
        w = zeros(length(x),length(t),length(lambda));
        for i=1:length(t)
            w(:,i,:) = funw(x,t(i),lambda,wl,wr,alphal,alphar,P);
        end
    end
    if model(2)==1
        % GODUNOV METHOD FOR LINEAR SYSTEMS
        % Riemman Problem
        wa = zeros(length(lambda),length(x));
        for i=1:length(x)
            if x(i)<=0 % Case 1:
                wa(:,i) = wl;
            else  % Case 2:
                wa(:,i) = wr;
            end
        end
        % Summoning Godunov Method
        wn = zeros(length(lambda),length(x),length(t));
        wn (:,:,1) = wa;
        for i=2:length(t)
            wn(:,:,i) = god(x,A,dtdx,wa);
            wa = wn(:,:,i);
        end
        wn = permute(wn,[2 3 1]);
    end
    if model(1)==1 && model(2)==1 
        plotfuncomp(x,t,wn,w,cfl,nx,model,framerate,hObject,handles);
    elseif model(1)==1
        plotfun(x,t,w,framerate,hObject,handles);
    elseif model(2)==1
        plotfun(x,t,wn,framerate,hObject,handles);
    end
    % SOLVING THE SYSTEM - CAUCHY PROBLEM
    if model(3)==1 || model(4)==1
        % Matrix of the hyperbolic system
        A = [a11 a12;a21 a22];
        % Eigenvectors and eigenvalues 
        [P,D] = eig(A); 
        lambda = diag(D);
        % Summoning w_cauchy
        for i=1:length(t)
            for j=1:length(x)
                w(j,i,:) = w_cauchy(x(j),t(i),P,lambda);
            end
        end
    end
    if model(4)==1
        % GODUNOV METHOD FOR LINEAR SYSTEMS
        % Cauchy Problem
            wa=w0_cauchy(x);
        % Summoning Godunov Method
        wn = zeros(length(lambda),length(x),length(t));
        wn (:,:,1) = wa;
        for i=2:length(t)
            wn(:,:,i) = god(x,A,dtdx,wa);
            wa = wn(:,:,i);
        end
        wn = permute(wn,[2 3 1]);
    end
    if model(3)==1 && model(4)==1 
        plotfuncomp(x,t,wn,w,cfl,nx,model,framerate,hObject,handles);
    elseif model(3)==1
        plotfun(x,t,w,framerate,hObject,handles);
    elseif model(4)==1
        plotfun(x,t,wn,framerate,hObject,handles);
    end
end
set(hObject,'Enable','on');
guidata(hObject, handles);

function edit4_CreateFcn(hObject, eventdata, handles)
function ayuda_menu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function version_Callback(hObject, eventdata, handles)
% hObject    handle to version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(hObject,'Value') == 1
    set(handles.checkbox3,'Enable','off');
    set(handles.checkbox4,'Enable','off');
else
    set(handles.checkbox3,'Enable','on');
    set(handles.checkbox4,'Enable','on');
end
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editRech4_Callback(hObject, eventdata, handles)
% hObject    handle to editRech4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRech4 as text
%        str2double(get(hObject,'String')) returns contents of editRech4 as a double


% --- Executes during object creation, after setting all properties.
function editRech4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRech4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editReair_Callback(hObject, eventdata, handles)
% hObject    handle to editReair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editReair as text
%        str2double(get(hObject,'String')) returns contents of editReair as a double


% --- Executes during object creation, after setting all properties.
function editReair_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editReair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to editdosado (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdosado as text
%        str2double(get(hObject,'String')) returns contents of editdosado as a double


% --- Executes during object creation, after setting all properties.
function editdosado_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdosado (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEje_Callback(hObject, eventdata, handles)
% hObject    handle to editEje (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEje as text
%        str2double(get(hObject,'String')) returns contents of editEje as a double


% --- Executes during object creation, after setting all properties.
function editEje_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEje (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edith2o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edith2o as text
%        str2double(get(hObject,'String')) returns contents of edith2o as a double


% --- Executes during object creation, after setting all properties.
function edith2o_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edith2o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to editco2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editco2 as text
%        str2double(get(hObject,'String')) returns contents of editco2 as a double


% --- Executes during object creation, after setting all properties.
function editco2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editco2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to editch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editch4 as text
%        str2double(get(hObject,'String')) returns contents of editch4 as a double


% --- Executes during object creation, after setting all properties.
function editch4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edito2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edito2 as text
%        str2double(get(hObject,'String')) returns contents of edito2 as a double


% --- Executes during object creation, after setting all properties.
function edito2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edito2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to editn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editn2 as text
%        str2double(get(hObject,'String')) returns contents of editn2 as a double


% --- Executes during object creation, after setting all properties.
function editn2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to editTsalida (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTsalida as text
%        str2double(get(hObject,'String')) returns contents of editTsalida as a double


% --- Executes during object creation, after setting all properties.
function editTsalida_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTsalida (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to editTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTmax as text
%        str2double(get(hObject,'String')) returns contents of editTmax as a double


% --- Executes during object creation, after setting all properties.
function editTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to editTmed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTmed as text
%        str2double(get(hObject,'String')) returns contents of editTmed as a double


% --- Executes during object creation, after setting all properties.
function editTmed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to editTAF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTAF as text
%        str2double(get(hObject,'String')) returns contents of editTAF as a double


% --- Executes during object creation, after setting all properties.
function editTAF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTAF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edito (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edito as text
%        str2double(get(hObject,'String')) returns contents of edito as a double


% --- Executes during object creation, after setting all properties.
function edito_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edito (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to editoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editoh as text
%        str2double(get(hObject,'String')) returns contents of editoh as a double


% --- Executes during object creation, after setting all properties.
function editoh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to editco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editco as text
%        str2double(get(hObject,'String')) returns contents of editco as a double


% --- Executes during object creation, after setting all properties.
function editco_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function Info_Callback(hObject, eventdata, handles)
% hObject    handle to Info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox1.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(hObject,'Value') == 1
    set(handles.checkbox3,'Enable','off');
    set(handles.checkbox4,'Enable','off');
else
    set(handles.checkbox3,'Enable','on');
    set(handles.checkbox4,'Enable','on');
end

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
if get(hObject,'Value') == 1
    set(handles.checkbox1,'Enable','off');
    set(handles.checkbox2,'Enable','off');
else
    set(handles.checkbox1,'Enable','on');
    set(handles.checkbox2,'Enable','on');
end

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
if get(hObject,'Value') == 1
    set(handles.checkbox1,'Enable','off');
    set(handles.checkbox2,'Enable','off');
else
    set(handles.checkbox1,'Enable','on');
    set(handles.checkbox2,'Enable','on');
end

% --- Executes on button press in axis51.
function axis51_Callback(hObject, eventdata, handles)
% hObject    handle to axis51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis51


% --- Executes on button press in axis61.
function axis61_Callback(hObject, eventdata, handles)
% hObject    handle to axis61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis61


% --- Executes on button press in axis12.
function axis12_Callback(hObject, eventdata, handles)
% hObject    handle to axis12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis12


% --- Executes on button press in axis22.
function axis22_Callback(hObject, eventdata, handles)
% hObject    handle to axis22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis22


% --- Executes on button press in axis32.
function axis32_Callback(hObject, eventdata, handles)
% hObject    handle to axis32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis32


% --- Executes on button press in axis42.
function axis42_Callback(hObject, eventdata, handles)
% hObject    handle to axis42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis42


% --- Executes on button press in axis52.
function axis52_Callback(hObject, eventdata, handles)
% hObject    handle to axis52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis52


% --- Executes on button press in axis62.
function axis62_Callback(hObject, eventdata, handles)
% hObject    handle to axis62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axis62

% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Autor_Callback(hObject, eventdata, handles)
web('https://www.linkedin.com/in/albertocuadralara/');
pause(1);
web('mailto::albertocuadralara@gmail.com');

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when figureMain is resized.
function figureMain_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figureMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonClose.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonPlot.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu14.
function popupmenu14_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu14


% --- Executes during object creation, after setting all properties.
function popupmenu14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editxa_Callback(hObject, eventdata, handles)
% hObject    handle to editxa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editxa as text
%        str2double(get(hObject,'String')) returns contents of editxa as a double


% --- Executes during object creation, after setting all properties.
function editxa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editxa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editxb_Callback(hObject, eventdata, handles)
% hObject    handle to editxb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editxb as text
%        str2double(get(hObject,'String')) returns contents of editxb as a double


% --- Executes during object creation, after setting all properties.
function editxb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editxb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editt0_Callback(hObject, eventdata, handles)
% hObject    handle to editt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editt0 as text
%        str2double(get(hObject,'String')) returns contents of editt0 as a double


% --- Executes during object creation, after setting all properties.
function editt0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editt0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittf_Callback(hObject, eventdata, handles)
% hObject    handle to edittf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittf as text
%        str2double(get(hObject,'String')) returns contents of edittf as a double


% --- Executes during object creation, after setting all properties.
function edittf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editnx_Callback(hObject, eventdata, handles)
% hObject    handle to editnx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editnx as text
%        str2double(get(hObject,'String')) returns contents of editnx as a double


% --- Executes during object creation, after setting all properties.
function editnx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editnx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editcfl_Callback(hObject, eventdata, handles)
% hObject    handle to editcfl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editcfl as text
%        str2double(get(hObject,'String')) returns contents of editcfl as a double


% --- Executes during object creation, after setting all properties.
function editcfl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcfl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edita11_Callback(hObject, eventdata, handles)
% hObject    handle to edita11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita11 as text
%        str2double(get(hObject,'String')) returns contents of edita11 as a double


% --- Executes during object creation, after setting all properties.
function edita11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edita21_Callback(hObject, eventdata, handles)
% hObject    handle to edita21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita21 as text
%        str2double(get(hObject,'String')) returns contents of edita21 as a double


% --- Executes during object creation, after setting all properties.
function edita21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edita12_Callback(hObject, eventdata, handles)
% hObject    handle to edita12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita12 as text
%        str2double(get(hObject,'String')) returns contents of edita12 as a double


% --- Executes during object creation, after setting all properties.
function edita12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edita22_Callback(hObject, eventdata, handles)
% hObject    handle to edita22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita22 as text
%        str2double(get(hObject,'String')) returns contents of edita22 as a double


% --- Executes during object creation, after setting all properties.
function edita22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editwl1_Callback(hObject, eventdata, handles)
% hObject    handle to editwl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwl1 as text
%        str2double(get(hObject,'String')) returns contents of editwl1 as a double


% --- Executes during object creation, after setting all properties.
function editwl1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editwr1_Callback(hObject, eventdata, handles)
% hObject    handle to editwr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwr1 as text
%        str2double(get(hObject,'String')) returns contents of editwr1 as a double


% --- Executes during object creation, after setting all properties.
function editwr1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editwl2_Callback(hObject, eventdata, handles)
% hObject    handle to editwl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwl2 as text
%        str2double(get(hObject,'String')) returns contents of editwl2 as a double


% --- Executes during object creation, after setting all properties.
function editwl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editwr2_Callback(hObject, eventdata, handles)
% hObject    handle to editwr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwr2 as text
%        str2double(get(hObject,'String')) returns contents of editwr2 as a double


% --- Executes during object creation, after setting all properties.
function editwr2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when TabA01Main is resized.
function TabA01Main_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to TabA01Main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function checkbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editPause_Callback(hObject, eventdata, handles)
% hObject    handle to editPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPause as text
%        str2double(get(hObject,'String')) returns contents of editPause as a double


% --- Executes during object creation, after setting all properties.
function editPause_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonStop.
function pushbuttonStop_Callback(hObject, eventdata, handles)
set(handles.pushbuttonStop,'userdata',1);
set(handles.pushbuttonPlot,'Enable','on');
% Update handles structure
guidata(hObject,handles);
% hObject    handle to pushbuttonStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonStopA2.
function pushbuttonStopA2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStopA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbuttonStopA2,'userdata',1);
set(handles.pushbuttonPlotA2,'Enable','on');
% Update handles structure
guidata(hObject,handles);

function editPauseA2_Callback(hObject, eventdata, handles)
% hObject    handle to editPauseA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPauseA2 as text
%        str2double(get(hObject,'String')) returns contents of editPauseA2 as a double


% --- Executes during object creation, after setting all properties.
function editPauseA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPauseA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPlotA2.
function pushbuttonPlotA2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set(handles.pushbuttonStopA2,'Enable','on');
set(handles.pushbuttonStopA2,'userdata',0);
guidata(hObject,handles);
%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver
model1 = get(handles.checkbox5,'Value');
model2 = get(handles.checkbox6,'Value');
model3 = get(handles.checkbox7,'Value');
model4 = get(handles.checkbox8,'Value');
model5 = get(handles.checkbox9,'Value');
model6 = get(handles.checkbox10,'Value');
model = [model1, model2, model3, model4, model5, model6];
% Mesh
a = str2num(get(handles.editxaA2,'string'));
b = str2num(get(handles.editxbA2,'string'));
nx = str2num(get(handles.editnxA2,'string'));
t0 = str2num(get(handles.editt0A2,'string'));
tmax = str2num(get(handles.edittfA2,'string'));
cfl = str2num(get(handles.editcflA2,'string'));
% Initial condition
wl = str2num(get(handles.editwlA2,'string'));
wr = str2num(get(handles.editwrA2,'string'));
% Discontinuity
xb = str2num(get(handles.editxb2,'string')); % Initial location of the discontinuity
tb = str2num(get(handles.edittbA2,'string'));
% Pause
framerate = str2num(get(handles.editPause,'string'));
% Construction data
deltax = abs((b-a))/(nx-1);
x = a:deltax:b;
s=(wl+wr)/2;
% Deltat from Courant number
deltat=cfl*deltax/max(abs(wl),abs(wr));
dtdx= deltat/deltax;
nt=round(tmax/deltat);
t = t0:deltat:tmax;
% Text editing
set(handles.textdxA2,'String',round(deltax,3));
set(handles.textdtA2,'String',round(deltat,3));
set(handles.textntA2,'String',nt);
% Initial condition
for i=1:length(x)
    if x(i)<=0
        w0(i)=wl;
    elseif x(i)<xb
        w0(i)=wl-x(i);
    else
        w0(i)=wr; 
    end
end

% Initial condition - Riemann Problem
% for i=1:length(x)
%     if x(i)<=xb
%         w0(i)=wl;
%     else
%         w0(i)=wr; 
%     end
% end
%%%%%%%%%%%% Solver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(model)>0
    plotfunnoneq(x,t,framerate,hObject,handles,...
        model,wl,wr,xb,tb,s,w0,dtdx,nx);
end
set(hObject,'Enable','on');
guidata(hObject, handles);

% --- Executes on button press in pushbuttonCloseA2.
function pushbuttonCloseA2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCloseA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --- Executes on button press in checkbox5.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function editwlA2_Callback(hObject, eventdata, handles)
% hObject    handle to editwlA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwlA2 as text
%        str2double(get(hObject,'String')) returns contents of editwlA2 as a double


% --- Executes during object creation, after setting all properties.
function editwlA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwlA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editwrA2_Callback(hObject, eventdata, handles)
% hObject    handle to editwrA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwrA2 as text
%        str2double(get(hObject,'String')) returns contents of editwrA2 as a double


% --- Executes during object creation, after setting all properties.
function editwrA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwrA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit68_Callback(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit68 as text
%        str2double(get(hObject,'String')) returns contents of edit68 as a double


% --- Executes during object creation, after setting all properties.
function edit68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit69_Callback(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit69 as text
%        str2double(get(hObject,'String')) returns contents of edit69 as a double


% --- Executes during object creation, after setting all properties.
function edit69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit70_Callback(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit70 as text
%        str2double(get(hObject,'String')) returns contents of edit70 as a double


% --- Executes during object creation, after setting all properties.
function edit70_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editxaA2_Callback(hObject, eventdata, handles)
% hObject    handle to editxaA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editxaA2 as text
%        str2double(get(hObject,'String')) returns contents of editxaA2 as a double


% --- Executes during object creation, after setting all properties.
function editxaA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editxaA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editxbA2_Callback(hObject, eventdata, handles)
% hObject    handle to editxbA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editxbA2 as text
%        str2double(get(hObject,'String')) returns contents of editxbA2 as a double


% --- Executes during object creation, after setting all properties.
function editxbA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editxbA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editt0A2_Callback(hObject, eventdata, handles)
% hObject    handle to editt0A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editt0A2 as text
%        str2double(get(hObject,'String')) returns contents of editt0A2 as a double


% --- Executes during object creation, after setting all properties.
function editt0A2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editt0A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittfA2_Callback(hObject, eventdata, handles)
% hObject    handle to edittfA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittfA2 as text
%        str2double(get(hObject,'String')) returns contents of edittfA2 as a double


% --- Executes during object creation, after setting all properties.
function edittfA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittfA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editnxA2_Callback(hObject, eventdata, handles)
% hObject    handle to editnxA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editnxA2 as text
%        str2double(get(hObject,'String')) returns contents of editnxA2 as a double


% --- Executes during object creation, after setting all properties.
function editnxA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editnxA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editcflA2_Callback(hObject, eventdata, handles)
% hObject    handle to editcflA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editcflA2 as text
%        str2double(get(hObject,'String')) returns contents of editcflA2 as a double


% --- Executes during object creation, after setting all properties.
function editcflA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcflA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenunoneq.
function popupmenunoneq_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenunoneq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenunoneq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenunoneq


% --- Executes during object creation, after setting all properties.
function popupmenunoneq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenunoneq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editxb2_Callback(hObject, eventdata, handles)
% hObject    handle to editxb2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editxb2 as text
%        str2double(get(hObject,'String')) returns contents of editxb2 as a double


% --- Executes during object creation, after setting all properties.
function editxb2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editxb2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittbA2_Callback(hObject, eventdata, handles)
% hObject    handle to edittbA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittbA2 as text
%        str2double(get(hObject,'String')) returns contents of edittbA2 as a double


% --- Executes during object creation, after setting all properties.
function edittbA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittbA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --------------------------------------------------------------------
function FVM1_Callback(hObject, eventdata, handles)
% hObject    handle to FVM1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('FVM1.pdf')


% --------------------------------------------------------------------
function FVM2_Callback(hObject, eventdata, handles)
% hObject    handle to FVM2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('FVM2.pdf')
