function varargout = untitled1(varargin)
% UNTITLED1 MATLAB code for untitled1.fig
%      UNTITLED1, by itself, creates a new UNTITLED1 or raises the existing
%      singleton*.
%
%      H = UNTITLED1 returns the handle to a new UNTITLED1 or the handle to
%      the existing singleton*.
%
%      UNTITLED1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED1.M with the given input arguments.
%
%      UNTITLED1('Property','Value',...) creates a new UNTITLED1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled1

% Last Modified by GUIDE v2.5 25-Nov-2019 18:45:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled1_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled1_OutputFcn, ...
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


% --- Executes just before untitled1 is made visible.
function untitled1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled1 (see VARARGIN)

% Choose default command line output for untitled1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
evalin('base','clear');
evalin('base','clc');

flagPlot        = 1;                                                        % (Flag) indicating whether you can plot
typeNACA        = str2num(get(handles.airfoil,'String'));

assignin('base','flagPlot',flagPlot);
assignin('base','typeNACA',typeNACA);

set(handles.graficas,'Enable','off');


% --- Executes on button press in graficas.
function graficas_Callback(hObject, eventdata, handles)
    xu           = evalin('base','xu');                                         % Upper surface X values
    xl           = evalin('base','xl');                                         % Lower surface X values
    yu           = evalin('base','yu');                                         % Upper surface Y values
    yl           = evalin('base','yl');  
    AoA          = evalin('base','AoA'); 
    AoARad          = evalin('base','AoARad');
    NumbAoA = evalin('base','NumbAoA');
    Uinf = evalin('base','Uinf');
    miu = evalin('base','miu');
    Rho = evalin('base','Rho');
    MaxZ=max(max(yu))
    assignin('base','MaxZ',MaxZ);
    [m n]=size(xu)
    for j=1:1:m
        if yu(j)==MaxZ
            MaxY=j;
        end
    end
    


    assignin('base','MaxY',MaxY);
    r1=m-MaxY;
    MaxX=xu(MaxY);
    
    assignin('base','MaxX',MaxX);
    AMeanx=[MaxY];
    DMeanx=[r1];
    AMeanz=[MaxY];
    DMeanz=[r1];
    for j=1:(MaxY)
        AMeanx(j)=xu(j);
        AMeanz(j)=yu(j);
    end
    for j=MaxY:m
        DMeanx(j-MaxY+1)=xu(j);
        DMeanz(j-MaxY+1)=yu(j);
    end
    DMeanz(m-MaxY+1)=0;
    AMeanx=rot90(AMeanx,3);
    DMeanx=rot90(DMeanx,3);
    AMeanz=rot90(AMeanz,3);
    DMeanz=rot90(DMeanz,3);
    assignin('base','AMeanx',AMeanx);
    assignin('base','AMeanz',AMeanz);
    assignin('base','DMeanx',DMeanx);
    assignin('base','DMeanz',DMeanz);
    for j=1:(MaxY)
        if isnan(AMeanx(j))
            AMeanx(j)=0;
        end
        if isnan(AMeanz(j))
            AMeanz(j)=0;
        end
    end
    
    c=1
    MaxXA=acosd(1-((2*MaxX)));
    MaxXR=deg2rad(MaxXA);
    assignin('base','MaxXA',MaxXA);
    assignin('base','MaxXR',MaxXR)
    dzdx1= polyfit(AMeanx, AMeanz, 2);
    dzdx2= polyfit(DMeanx, DMeanz, 2);
    Aa=dzdx1(1);
    B1=dzdx1(2);
    C1=dzdx1(3);
    Aaa=dzdx2(1);
    B2=dzdx2(2);
    C2=dzdx2(3);
    assignin('base','Aaa',Aaa);
    fun11=@(theta) -Aa.*cos(theta)+B1+Aa;
    fun12=@(theta) -Aaa.*cos(theta)+B2+Aaa;
    fun21=@(theta) (-Aa.*cos(theta)+B1+Aa).*cos(theta);
    fun22=@(theta) (-Aaa.*cos(theta)+B2+Aaa).*cos(theta);
    fun31=@(theta) (-Aa.*cos(theta)+B1+Aa).*cos(2*theta);
    fun32=@(theta) (-Aaa.*cos(theta)+B2+Aaa).*cos(2*theta);
    fun41=@(theta) (-Aa.*cos(theta)+B1+Aa).*(cos(theta)-1);
    fun42=@(theta) (-Aaa.*cos(theta)+B2+Aaa).*(cos(theta)-1);

    Aos=size(AoA);
    
    assignin('base','dzdx1',dzdx1)
    assignin('base','dzdx2',dzdx2)
    for j=1:Aos(2)

        Cl(j)=(2*pi*AoARad(j))-(2*(integral(fun11,0,MaxXR)+integral(fun12,MaxXR,pi)))-(2*(integral(fun21,0,MaxXR)+integral(fun22,MaxXR,pi)));

    end

    for j=1:Aos(2)
        Alphazerol=AoARad(j)-(Cl(j)/(2*pi));
    end
    
    AlphazerolA=radtodeg(Alphazerol);
    assignin('base','Alpha0A',AlphazerolA);
    assignin('base','Cl',Cl);
    assignin('base','Alphazerol',Alphazerol);
    axes(handles.axes1);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Lift Coefficient');                                                      % Title
    xlabel('angle of attack []');                                                           % X-axis label
    ylabel('Lift Coefficient []');                                                           % Y-axis label
    axis equal;
    plotU = plot(AoA,Cl,'k-');
    Re=(Rho*Uinf)/miu
    Cdf=2*(1.328/(sqrt(Re)));                                               % Plot the upper surface (black-line)
    po=(2*Aaa*xu(m)+B2);
    for j=1:NumbAoA
        Vl(j)=((Uinf)*sin(AoARad(j)-atan(po)));
    end
    assignin('base','po',po);
    fun6=@(x) (Uinf^8.21);
    dsaoda=integral(fun6,0,1, 'ArrayValued',true);
    Teta=sqrt(((0.440*miu)/(Rho*Uinf^9.21))*integral(fun6,0,1, 'ArrayValued',true));
    for j=1:NumbAoA
        Cd1(j)=2*(Teta)*(Vl(j)/Uinf)^(6.45/2);
        Ar=9.5/9.5;
        Cd1ar(j)=Cl(j)^2/(pi*(Ar));
        %Plotear
        Cd(j)=(Cdf+Cd1(j)+Cd1ar(j));
    end
    assignin('base','Cd',Cd);
    assignin('base','dskad',dsaoda);
    assignin('base','Teta',Teta);
    assignin('base','Cd1',Cd1);
    assignin('base','Vl',Vl);
    axes(handles.axes2);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Drag Coefficient');                                                      % Title
    xlabel('angle of attack []');                                                           % X-axis label
    ylabel('Drag Coefficient []');                                                           % Y-axis label
    ylim([-0.2,0.2])
    axis equal;
    plotU = plot(AoA,Cd,'k-');
     axes(handles.axes3);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Drag Polar');                                                      % Title
    xlabel('Drag Coefficient []');                                                           % X-axis label
    ylabel('Lift Coefficient []');                                                           % Y-axis label
    axis equal;
    plotU = plot(Cd,Cl,'k-');
    
% --- Executes on button press in perfil.
function perfil_Callback(hObject, eventdata, handles)
typeNACA = get(handles.airfoil,'String');
[~,sizeInput] = size(typeNACA);
AoA=[-10:1:16];
AoARad=deg2rad(AoA);
Mach=str2num(get(handles.Mach,'string'));
Height=str2num(get(handles.altura,'string')); %metros
[T,a,P,Rho]=atmosisa(Height);
Uinf=a*Mach;
q=(1/2)*Rho*Uinf^2;
miu=(1.458e-6)*T^(3/2)*((1/(T+110.4)));
Re=(Rho*Uinf)/miu
Cdf=2*(1.328/(sqrt(Re)));
assignin('base','AoA',AoA);
assignin('base','AoARad',AoARad);
assignin('base','Mach',Mach);
assignin('base','Height',Height);
assignin('base','Uinf',Uinf);
assignin('base','miu',miu);
assignin('base','Re',Re);
assignin('base','Cdf',Cdf);
assignin('base','Rho',Rho);
assignin('base','T',T);

NumbAoA=size(AoA,2);
assignin('base','NumbAoA',NumbAoA);
if (sizeInput == 4)
    Minit = str2double(typeNACA(1));
    Pinit = str2double(typeNACA(2));
    Tinit = str2double(typeNACA(3:4));
    M = Minit/100;
    P = Pinit/10;
    T = Tinit/100;
    a0 = 0.2969;
    a1 = -0.1260;
    a2 = -0.3516;
    a3 = 0.2843;
    a4= -0.1036;
    gridPts=57;
    beta = linspace(0,pi,gridPts)';
    x = (0.5*(1-cos(beta)));
    yc     = ones(gridPts,1);
    dyc_dx = ones(gridPts,1);
    theta  = ones(gridPts,1);
    for i = 1:1:gridPts
        if (x(i) >= 0 && x(i) < P)
            yc(i)     = (M/P^2)*((2*P*x(i))-x(i)^2);
            dyc_dx(i) = ((2*M)/(P^2))*(P-x(i));
        elseif (x(i) >=P && x(i) <=1)
            yc(i)     = (M/(1-P)^2)*(1-(2*P)+(2*P*x(i))-(x(i)^2));
            dyc_dx(i) = ((2*M)/((1-P)^2))*(P-x(i));
        end
        theta(i) = atan(dyc_dx(i));
    end
    yt = 5*T.*((a0.*sqrt(x)) + (a1.*x) + (a2.*x.^2) + (a3.*x.^3) + (a4.*x.^4));
    xu = x(:)  - yt(:).*sin(theta);
    yu = yc(:) + yt(:).*cos(theta);
    xl = x(:) + yt(:).*sin(theta);
    yl = yc(:) - yt(:).*cos(theta);
    xc = x;
    x_qc = 0.25;
    y_qc = 0;
    axisOfRot = [0 0 1];                                                        % Rotate about Z-axis
    origin    = [x_qc y_qc 0];                                                  % Origin of rotation is at the quarter-chord
    thetaRot  = 0;  
    axes(handles.axes4);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Airfoil Plot');                                                      % Title
    xlabel('x/c []');                                                           % X-axis label
    ylabel('y/c []');                                                           % Y-axis label
    axis equal;
    plotU = plot(xu,yu,'k-');                                                   % Plot the upper surface (black-line)
    plotL = plot(xl,yl,'k-');
    plotC = plot(xc,yc,'r-');
    z_x =polyfit(xc,yc,2)
    set(handles.graficas,'Enable','on');
    assignin('base','xu',xu);
    assignin('base','xl',xl);
    assignin('base','yu',yu);
    assignin('base','yl',yl);
    assignin('base','theta',theta);
    
end
if (sizeInput == 5)
    % Typical Inputs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    name=get(handles.airfoil,'String');
    gridPts=56;
    Espaciado=1;
    Tiene_TE=1;
    % % [[Calculating key parameters-----------------------------------------]]
    cld=str2num(name(1))*(3/2)/10;
    p=5*str2num(name(2))/100;
    r=str2num(name(3));
    t=str2num(name(4:5))/100;
    a0= 0.2969;
    a1=-0.1260;
    a2=-0.3516;
    a3= 0.2843;
    if Tiene_TE ==1
        a4=-0.1015; % For finite thick TE
    else
        a4=-0.1036;  % For zero thick TE
    end
    % % [[Giving x-spacing---------------------------------------------------]]
    if Espaciado==1
        beta=linspace(0,pi,gridPts+1)';
        x=(0.5*(1-cos(beta))); % Half cosine based spacing
        encabezado=['NACA' name ' : [' num2str(2*gridPts) 'panels,Half cosine x-spacing]'];
    else
        x=linspace(0,1,gridPts+1)';
        encabezado=['NACA' name ' : [' num2str(2*gridPts) 'panels,Uniform x-spacing]'];
    end
    yt=(t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);
    if r
    P=[ 0.1 0.15 0.2 0.25 ];
    M=[ 0.13 0.2170 0.318 0.441 ];
    K=[ 51.99 15.793 6.520 3.191 ];
    else
    P=[0.05 0.1 0.15 0.2 0.25];
    M=[0.0580 0.1260 0.2025 0.2900 0.3910];
    K=[361.4 51.64 15.957 6.643 3.230];
    end
    m=spline(P,M,p);
    k1=spline(M,K,m);
    xc1=x(x<m);
    xc2=x(x>=m);
    xc=[xc1 ; xc2];
    if p==0
        xu=x;
        yu=yt;
        xl=x;
        yl=-yt;

        yc=zeros(size(xc));
    else
        yc1=(1/6)*k1*( xc1.^3-3*m*xc1.^2+m^2*(3-m)*xc1 );
        yc2=(1/6)*k1*m^3*(1-xc2);
        yc=(cld/0.3)*[yc1 ; yc2];
        dyc1_dx=(1/6)*k1*( 3*xc1.^2-6*m*xc1+m^2*(3-m) );
        dyc2_dx=repmat(-(1/6)*k1*m^3,size(xc2));
        dyc_dx=[dyc1_dx ; dyc2_dx];
        theta=atan(dyc_dx);
        xu=x-yt.*sin(theta);
        yu=yc+yt.*cos(theta);
        xl=x+yt.*sin(theta);
        yl=yc-yt.*cos(theta);
    end

    name=['NACA ' name];
    x=[flipud(xu) ; xl(2:end)];
    z=[flipud(yu) ; yl(2:end)];
    indx1=1:min( find(x==min(x)) );  % Upper surface indices
    indx2=min( find(x==min(x)) ):length(x); % Lower surface indices
    xu=x(indx1); % Upper Surface x
    yu=z(indx1); % Upper Surface z
    xl=x(indx2); % Lower Surface x
    yl=z(indx2); % Lower Surface z
    axes(handles.axes4);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Airfoil Plot');                                                      % Title
    xlabel('x/c []');                                                           % X-axis label
    ylabel('y/c []');                                                           % Y-axis label
    axis equal;
    plotU = plot(xu,yu,'k-');                                                   % Plot the upper surface (black-line)
    plotL = plot(xl,yl,'k-');
    plotC = plot(xc,yc,'r-');
    assignin('base','xu',xu);
    assignin('base','xl',xl);
    assignin('base','yu',yu);
    assignin('base','yl',yl);
    assignin('base','yc',yc);
    assignin('base','xc',x);
    set(handles.graficas,'Enable','on');
end
if (sizeInput == 6)
    n=str2num(get(handles.airfoil,'String'));
    c=1
    s=56
    beta=linspace(0,pi,s);  % Angle for cosine spacing
    x=(1-cos(beta))/2;
    alpha=0
    t=rem(n,100)/100;   % Maximum thickness as fraction of chord (two last digits)
    sym=0;  % Symetric airfoil variable
    alpha=alpha/180*pi;
    yc=zeros(1,s); % Mean camber vector prelocation
    dyc_dx=zeros(1,s); 
    y_t=t/0.2*(0.2969*sqrt(x)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4); % Thickness y coordinate with closed trailing edge
    ser=floor(n/100000);    % Number of series (1st digit)
    a=rem(floor(n/10000),10)/10;  % Chordwise position of minimum pressure (2nd digit)
    c_li=rem(floor(n/100),10)/10;  % Design lift coefficient (4th digit)
    g=-1/(1-a)*(a^2*(1/2*log(a)-1/4)+1/4);  % G constant calculation
    h=1/(1-a)*(1/2*(1-a)^2*log(1-a)-1/4*(1-a)^2)+g; % H constant calculation
    yc=c_li/(2*pi*(a+1))*(1/(1-a)*(1/2*(a-x).^2.*log(abs(a-x))-1/2*(1-x).^2.*log(1-x)+1/4*(1-x).^2-1/4*(a-x).^2)-x.*log(x)+g-h*x)+(1/2-x)*sin(alpha); % Mean camber y coordinate
    dyc_dx=-(c_li*(h+log(x)-(x/2-a/2+(log(1-x).*(2*x-2))/2+(log(abs(a-x)).*(2*a-2*x))/2+(sign(a-x).*(a-x).^2)./(2*abs(a-x)))/(a-1)+1))/(2*pi*(a+1)*cos(alpha))-tan(alpha);    % Mean camber first derivative
    theta=atan(dyc_dx); % Angle for modifying x coordinate
    x=1/2-(1/2-x)*cos(alpha);   % X coordinate rotation
    %----------------------- COORDINATE ASSIGNATION ---------------------------
    xu=(x-y_t.*sin(theta))*c;
    xu=rot90(xu)% X extrados coordinate
    xl=(x+y_t.*sin(theta))*c;
    xl=rot90(xl)% X intrados coordinate
    yu=(yc+y_t.*cos(theta))*c;
    yu=rot90(yu)% Y extrados coordinate
    yl=(yc-y_t.*cos(theta))*c
    yl=rot90(yl);    % Y intrados coordinate
    axes(handles.axes4);                                                  % Select axis for plotting
    cla; hold on; grid on;
    title('Airfoil Plot');                                                      % Title
    xlabel('x/c []');                                                           % X-axis label
    ylabel('y/c []');                                                           % Y-axis label
    axis equal;
    plotU = plot(xu,yu,'k-');                                                   % Plot the upper surface (black-line)
    plotL = plot(xl,yl,'k-');
    plotC = plot(x,yc,'r-');
    assignin('base','xu',xu);
    assignin('base','xl',xl);
    assignin('base','yu',yu);
    assignin('base','yl',yl);
    set(handles.graficas,'Enable','on');
end
    
function airfoil_Callback(hObject, eventdata, handles)
typeNACA = get(handles.airfoil,'String');
assignin('base','typeNACA',typeNACA);

% Check if it's exactly 4 digits
[~,sizeInput] = size(typeNACA);
if (sizeInput == 4 || sizeInput == 5 || sizeInput == 6)
    set(handles.textStatus,'String','Listo para graficar');
    set(handles.textStatus,'ForegroundColor',[0 0 0]);
    drawnow();
    set(handles.perfil,'Enable','on');
else
    set(handles.textStatus,'String','Introduce un perfil de 4 o 5 dígitos');
    set(handles.textStatus,'ForegroundColor',[1 0 0]);
    drawnow();
    set(handles.perfil,'Enable','off');
end
% Hints: get(hObject,'String') returns contents of airfoil as text
%        str2double(get(hObject,'String')) returns contents of airfoil as a double


% --- Executes during object creation, after setting all properties.
function airfoil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to airfoil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mach_Callback(hObject, eventdata, handles)
% hObject    handle to Mach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mach as text
%        str2double(get(hObject,'String')) returns contents of Mach as a double


% --- Executes during object creation, after setting all properties.
function Mach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function altura_Callback(hObject, eventdata, handles)
% hObject    handle to altura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of altura as text
%        str2double(get(hObject,'String')) returns contents of altura as a double


% --- Executes during object creation, after setting all properties.
function altura_CreateFcn(hObject, eventdata, handles)
% hObject    handle to altura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
