function myGui
% создание графического окна с заголовком myplotgui и без надписи Figure 1
hF = figure('Name', 'Interface', 'NumberTitle', 'on');
set(hF, 'Position', [50 200 1500 700])

% Matrix A
uicontrol('Style', 'frame', ...
 'Position', [40 665 55 35], 'Tag', 'hAframe');
uicontrol('Style', 'text', ...
 'Position', [42 668 52 30],...
 'String', 'A()',...
 'FontSize', 18, 'Tag', 'hAmatrix');
% cells of matrix A
uicontrol('Style', 'edit', ...
 'Position', [40 630 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hAmatrix_00');
uicontrol('Style', 'edit', ...
 'Position', [75 630 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hAmatrix_01');
uicontrol('Style', 'edit', ...
 'Position', [40 590 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hAmatrix_10');
uicontrol('Style', 'edit', ...
 'Position', [75 590 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hAmatrix_11');

% Matrix B
uicontrol('Style', 'frame', ...
 'Position', [145 665 55 35], 'Tag', 'hBframe');
uicontrol('Style', 'text', ...
 'Position', [147 668 52 30],...
 'String', 'B()',...
 'FontSize', 18, 'Tag', 'hBmatrix');
% cells of matrix B
uicontrol('Style', 'edit', ...
 'Position', [145 630 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hBmatrix_00');
uicontrol('Style', 'edit', ...
 'Position', [180 630 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hBmatrix_01');
uicontrol('Style', 'edit', ...
 'Position', [145 590 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hBmatrix_10');
uicontrol('Style', 'edit', ...
 'Position', [180 590 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hBmatrix_11');

% function f(t)
uicontrol('Style', 'frame', ...
 'Position', [250 665 90 35], 'Tag', 'hfframe');
uicontrol('Style', 'text', ...
 'Position', [252 668 85 30],...
 'String', 'f()',...
 'FontSize', 18, 'Tag', 'hftext');
uicontrol('Style', 'edit', ...
 'Position', [250 630 90 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hf_00');
uicontrol('Style', 'edit', ...
 'Position', [250 590 90 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hf_01');

% q1, q2, q3, q4
%q1
uicontrol('Style', 'text', ...
 'Position', [35 530 30 30],...
 'String', 'q1',...
 'FontSize', 16, 'Tag', 'hQ1');
uicontrol('Style', 'edit', ...
 'Position', [40 500 20 30],...
 'String', '4',... 
 'FontSize', 14, 'Tag', 'hQ1_0');
uicontrol('Style', 'edit', ...
 'Position', [40 450 20 30],...
 'String', '4',... 
 'FontSize', 14, 'Tag', 'hQ1_1');

%q2
uicontrol('Style', 'text', ...
 'Position', [70 530 30 30],...
 'String', 'q2',...
 'FontSize', 16, 'Tag', 'hQ2');
uicontrol('Style', 'edit', ...
 'Position', [75 500 20 30],...
 'String', '6',... 
 'FontSize', 14, 'Tag', 'hQ2_0');
uicontrol('Style', 'edit', ...
 'Position', [75 450 20 30],...
 'String', '7',... 
 'FontSize', 14, 'Tag', 'hQ2_1');

%q3
uicontrol('Style', 'text', ...
 'Position', [105 530 30 30],...
 'String', 'q3',...
 'FontSize', 16, 'Tag', 'hQ3');
uicontrol('Style', 'edit', ...
 'Position', [110 500 20 30],...
 'String', '7',... 
 'FontSize', 14, 'Tag', 'hQ3_0');
uicontrol('Style', 'edit', ...
 'Position', [110 450 20 30],...
 'String', '6',... 
 'FontSize', 14, 'Tag', 'hQ3_1');

%q4
uicontrol('Style', 'text', ...
 'Position', [140 530 30 30],...
 'String', 'q4',...
 'FontSize', 16, 'Tag', 'hQ4');
uicontrol('Style', 'edit', ...
 'Position', [145 500 20 30],...
 'String', '8',... 
 'FontSize', 14, 'Tag', 'hQ4_0');
uicontrol('Style', 'edit', ...
 'Position', [145 450 20 30],...
 'String', '4',... 
 'FontSize', 14, 'Tag', 'hQ4_1');

% constants
%alpha
uicontrol('Style', 'text', ...
 'Position', [200 530 30 30],...
 'String', 'al',...
 'FontSize', 16, 'Tag', 'halpha');
uicontrol('Style', 'edit', ...
 'Position', [205 500 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'halpha_0');

%beta
uicontrol('Style', 'text', ...
 'Position', [235 530 30 30],...
 'String', 'bt',...
 'FontSize', 16, 'Tag', 'hbeta');
uicontrol('Style', 'edit', ...
 'Position', [240 500 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hbeta_0');

%a
uicontrol('Style', 'text', ...
 'Position', [270 530 30 30],...
 'String', 'a',...
 'FontSize', 16, 'Tag', 'ha');
uicontrol('Style', 'edit', ...
 'Position', [275 500 20 30],...
 'String', '4',... 
 'FontSize', 14, 'Tag', 'ha_0');

%b
uicontrol('Style', 'text', ...
 'Position', [305 530 30 30],...
 'String', 'b',...
 'FontSize', 16, 'Tag', 'hb');
uicontrol('Style', 'edit', ...
 'Position', [310 500 20 30],...
 'String', '4',... 
 'FontSize', 14, 'Tag', 'hb_0');

%c
uicontrol('Style', 'text', ...
 'Position', [340 530 30 30],...
 'String', 'c',...
 'FontSize', 16, 'Tag', 'hc');
uicontrol('Style', 'edit', ...
 'Position', [345 500 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hc_0');

%r
uicontrol('Style', 'text', ...
 'Position', [375 530 30 30],...
 'String', 'r',...
 'FontSize', 16, 'Tag', 'hr');
uicontrol('Style', 'edit', ...
 'Position', [380 500 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hr_0');

% x0
%{
uicontrol('Style', 'text', ...
 'Position', [375 668 30 30],...
 'String', 'x0',...
 'FontSize', 16, 'Tag', 'hx0');
uicontrol('Style', 'edit', ...
 'Position', [380 630 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hx0_0');
uicontrol('Style', 'edit', ...
 'Position', [380 590 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'hx0_1');
%}
% time
%t0
uicontrol('Style', 'text', ...
 'Position', [35 390 30 30],...
 'String', 't0',...
 'FontSize', 16, 'Tag', 'ht0');
uicontrol('Style', 'edit', ...
 'Position', [40 360 20 30],...
 'String', '0',... 
 'FontSize', 14, 'Tag', 'ht0_0');

%tmax
uicontrol('Style', 'text', ...
 'Position', [70 390 30 30],...
 'String', 'tm',...
 'FontSize', 16, 'Tag', 'htMax');
uicontrol('Style', 'edit', ...
 'Position', [75 360 20 30],...
 'String', '0.7',... 
 'FontSize', 14, 'Tag', 'htMax_0');

%n
uicontrol('Style', 'text', ...
 'Position', [105 390 30 30],...
 'String', 'n',...
 'FontSize', 16, 'Tag', 'hn');
uicontrol('Style', 'edit', ...
 'Position', [110 360 40 30],...
 'String', '15',... 
 'FontSize', 14, 'Tag', 'hn_0');

%%{
%RelTol
uicontrol('Style', 'text', ...
 'Position', [35 330 60 30],...
 'String', 'RelT',...
 'FontSize', 16, 'Tag', 'hRelT');
uicontrol('Style', 'edit', ...
 'Position', [40 300 55 30],...
 'String', '1e-3',... 
 'FontSize', 14, 'Tag', 'hRelT_0');

%AbsTol
uicontrol('Style', 'text', ...
 'Position', [35 270 60 30],...
 'String', 'AbsT',...
 'FontSize', 16, 'Tag', 'hAbsT');
uicontrol('Style', 'edit', ...
 'Position', [40 240 55 30],...
 'String', '1e-6',... 
 'FontSize', 14, 'Tag', 'hAbsT_0');
%}
 
%start culculation
uicontrol('Style', 'pushbutton', ...
 'Position', [205 410 195 70],...
 'String', 'Решить!',...
 'FontSize', 16,...
 'Tag', 'hdes',...
 'Callback', @GetValues);
%записали указали на объекты в структуру
handles = guihandles(hF);
guidata(hF, handles)
set(hF, 'HandleVisibility', 'callback');

%---------------------------------------------------------------------%
function GetValues(src,evt)
handles = guidata(src);

values = containers.Map;

%getting values of matrix A with ckecking and converting to double
Amatrix_00 = get(handles.hAmatrix_00, 'String');
if strcmp(Amatrix_00, '')
    error('No value at Amatrix_00');
end
Amatrix_00 = str2double(Amatrix_00);
if isnan(Amatrix_00) 
    error('Some trash at Amatrix_00');
end   
values('Amatrix_00') = Amatrix_00;

%
Amatrix_01 = get(handles.hAmatrix_01, 'String');
if strcmp(Amatrix_01, '')
    error('No value at Amatrix_01');
end
Amatrix_01 = str2double(Amatrix_01);
if isnan(Amatrix_01) 
    error('Some trash at Amatrix_01');
end  
values('Amatrix_01') = Amatrix_01;

%
Amatrix_10 = get(handles.hAmatrix_10, 'String');
if strcmp(Amatrix_10, '')
    error('No value at Amatrix_10');
end
Amatrix_10 = str2double(Amatrix_10);
if isnan(Amatrix_10) 
    error('Some trash at Amatrix_10');
end  
values('Amatrix_10') = Amatrix_10;

%
Amatrix_11 = get(handles.hAmatrix_11, 'String');
if strcmp(Amatrix_11, '')
    error('No value at Amatrix_11');
end
Amatrix_11 = str2double(Amatrix_11);
if isnan(Amatrix_11) 
    error('Some trash at Amatrix_11');
end  
values('Amatrix_11') = Amatrix_11;

%getting values of matrix B with ckecking and converting to double
Bmatrix_00 = get(handles.hBmatrix_00, 'String');
if strcmp(Bmatrix_00, '')
    error('No value at Bmatrix_00');
end
Bmatrix_00 = str2double(Bmatrix_00);
if isnan(Bmatrix_00) 
    error('Some trash at Bmatrix_00');
end   
values('Bmatrix_00') = Bmatrix_00;

%
Bmatrix_01 = get(handles.hBmatrix_01, 'String');
if strcmp(Bmatrix_01, '')
    error('No value at Bmatrix_01');
end
Bmatrix_01 = str2double(Bmatrix_01);
if isnan(Bmatrix_01) 
    error('Some trash at Bmatrix_01');
end  
values('Bmatrix_01') = Bmatrix_01;
%

Bmatrix_10 = get(handles.hBmatrix_10, 'String');
if strcmp(Bmatrix_10, '')
    error('No value at Bmatrix_10');
end
Bmatrix_10 = str2double(Bmatrix_10);
if isnan(Bmatrix_10) 
    error('Some trash at Bmatrix_10');
end  
values('Bmatrix_10') = Bmatrix_10;

%
Bmatrix_11 = get(handles.hBmatrix_11, 'String');
if strcmp(Bmatrix_11, '')
    error('No value at Bmatrix_11');
end
Bmatrix_11 = str2double(Bmatrix_11);
if isnan(Bmatrix_11) 
    error('Some trash at Bmatrix_11');
end  
values('Bmatrix_11') = Bmatrix_11;

% getting values of vector q1
Q1_0 = get(handles.hQ1_0, 'String');
if strcmp(Q1_0, '')
    error('No value at Q1_0');
end
Q1_0 = str2double(Q1_0);
if isnan(Q1_0) 
    error('Some trash at Q1_0');
end 
values('Q1_0') = Q1_0;

%
Q1_1 = get(handles.hQ1_1, 'String');
if strcmp(Q1_1, '')
    error('No value at Q1_1');
end
Q1_1 = str2double(Q1_1);
if isnan(Q1_1) 
    error('Some trash at Q1_1');
end 
values('Q1_1') = Q1_1;

% getting values of vector q2
Q2_0 = get(handles.hQ2_0, 'String');
if strcmp(Q2_0, '')
    error('No value at Q2_0');
end
Q2_0 = str2double(Q2_0);
if isnan(Q2_0) 
    error('Some trash at Q2_0');
end 
values('Q2_0') = Q2_0;

%
Q2_1 = get(handles.hQ2_1, 'String');
if strcmp(Q2_1, '')
    error('No value at Q2_1');
end
Q2_1 = str2double(Q2_1);
if isnan(Q2_1) 
    error('Some trash at Q2_1');
end 
values('Q2_1') = Q2_1;

% getting values of vector q3
Q3_0 = get(handles.hQ3_0, 'String');
if strcmp(Q3_0, '')
    error('No value at Q3_0');
end
Q3_0 = str2double(Q3_0);
if isnan(Q3_0) 
    error('Some trash at Q3_0');
end 
values('Q3_0') = Q3_0;

%
Q3_1 = get(handles.hQ3_1, 'String');
if strcmp(Q3_1, '')
    error('No value at Q3_1');
end
Q3_1 = str2double(Q3_1);
if isnan(Q3_1) 
    error('Some trash at Q3_1');
end 
values('Q3_1') = Q3_1;

% getting values of vector q4
Q4_0 = get(handles.hQ4_0, 'String');
if strcmp(Q4_0, '')
    error('No value at Q4_0');
end
Q4_0 = str2double(Q4_0);
if isnan(Q4_0) 
    error('Some trash at Q4_0');
end 
values('Q4_0') = Q4_0;

%
Q4_1 = get(handles.hQ4_1, 'String');
if strcmp(Q4_1, '')
    error('No value at Q4_1');
end
Q4_1 = str2double(Q4_1);
if isnan(Q4_1) 
    error('Some trash at Q4_1');
end 
values('Q4_1') = Q4_1;

% getting values of constants 
%alpha
alpha_0 = get(handles.halpha_0, 'String');
if strcmp(alpha_0, '')
    error('No value at alpha_0');
end
alpha_0 = str2double(alpha_0);
if isnan(alpha_0) 
    error('Some trash at alpha_0');
end 
values('alpha_0') = alpha_0;

%beta
beta_0 = get(handles.hbeta_0, 'String');
if strcmp(beta_0, '')
    error('No value at beta_0');
end
beta_0 = str2double(beta_0);
if isnan(beta_0) 
    error('Some trash at beta_0');
end 
values('beta_0') = beta_0;

%a
a_0 = get(handles.ha_0, 'String');
if strcmp(a_0, '')
    error('No value at a_0');
end
a_0 = str2double(a_0);
if isnan(a_0) 
    error('Some trash at a_0');
end 
values('a_0') = a_0;

%b
b_0 = get(handles.hb_0, 'String');
if strcmp(b_0, '')
    error('No value at b_0');
end
b_0 = str2double(b_0);
if isnan(b_0) 
    error('Some trash at b_0');
end 
values('b_0') = b_0;

%c
c_0 = get(handles.hc_0, 'String');
if strcmp(c_0, '')
    error('No value at c_0');
end
c_0 = str2double(c_0);
if isnan(c_0) 
    error('Some trash at c_0');
end 
values('c_0') = c_0;

%r
r_0 = get(handles.hr_0, 'String');
if strcmp(r_0, '')
    error('No value at r_0');
end
r_0 = str2double(r_0);
if isnan(r_0) 
    error('Some trash at r_0');
end 
values('r_0') = r_0;

% getting value of x0
%{
% x0[0]
x0_0 = get(handles.hx0_0, 'String');
if strcmp(x0_0, '')
    error('No value at x0_0');
end
x0_0 = str2double(x0_0);
if isnan(x0_0) 
    error('Some trash at x0_0');
end 
values('x0_0') = x0_0;

%x0[1]
x0_1 = get(handles.hx0_1, 'String');
if strcmp(x0_1, '')
    error('No value at x0_1');
end
x0_1 = str2double(x0_1);
if isnan(x0_1) 
    error('Some trash at x0_1');
end 
values('x0_1') = x0_1;
%}

%getting value of t0
t0_0 = get(handles.ht0_0, 'String');
if strcmp(t0_0, '')
    error('No value at t0_0');
end
t0_0 = str2double(t0_0);
if isnan(t0_0) 
    error('Some trash at t0_0');
end 
values('t0_0') = t0_0;

%getting value of tMax
tMax_0 = get(handles.htMax_0, 'String');
if strcmp(tMax_0, '')
    error('No value at tMax_0');
end
tMax_0 = str2double(tMax_0);
if isnan(tMax_0) 
    error('Some trash at tMax_0');
end 
values('tMax_0') = tMax_0;

%getting value of n
n = get(handles.hn_0, 'String');
if strcmp(n, '')
    error('No value at n');
end
n = str2double(n);
if isnan(n) 
    error('Some trash at n');
end 
values('n') = n;

%getting value of RelTol
RelT = get(handles.hRelT_0, 'String');
if strcmp(RelT, '')
    error('No value at RelT');
end
RelT = str2double(RelT);
if isnan(RelT) 
    error('Some trash at RelT');
end 
values('RelT') = RelT;

%getting value of AbsTol
AbsT = get(handles.hAbsT_0, 'String');
if strcmp(AbsT, '')
    error('No value at AbsT');
end
AbsT = str2double(AbsT);
if isnan(n) 
    error('Some trash at AbsT');
end 
values('AbsT') = AbsT;

%getting value of f(t)
f_00 = get(handles.hf_00, 'String');
if strcmp(f_00, '')
    error('No value at f_00');
end
f_00 = str2double(f_00);
if isnan(f_00) 
    error('Some trash at f_00');
end 
values('f_00') = f_00;

f_01 = get(handles.hf_01, 'String');
if strcmp(f_01, '')
    error('No value at f_01');
end
f_01 = str2double(f_01);
if isnan(f_01) 
    error('Some trash at f_01');
end 
values('f_01') = f_01;
StartCulculation(values)



