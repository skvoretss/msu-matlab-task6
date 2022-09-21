function myGui
% создание графического окна с заголовком myplotgui и без надписи Figure 1
hF = figure('Name', 'Interface', 'NumberTitle', 'on');
set(hF, 'Position', [100 300 400 300])

%variant of the task
%var
uicontrol('Style', 'text', ...
 'Position', [255 165 40 30],...
 'String', 'var',...
 'FontSize', 16, 'Tag', 'hvar');
uicontrol('Style', 'edit', ...
 'Position', [260 130 30 30],...
 'String', '2',... 
 'FontSize', 14, 'Tag', 'hvar_0');

%%{
% explanation
uicontrol('Style', 'text', ...
 'Position', [10 215 400 30],...
 'String', 'var = 1 is u1 >= 0; var = 2 is any u1',...
 'FontSize', 12, 'Tag', 'hvar');
 %}
% constants
%k1
uicontrol('Style', 'text', ...
 'Position', [10 165 30 30],...
 'String', 'k1',...
 'FontSize', 16, 'Tag', 'hk1');
uicontrol('Style', 'edit', ...
 'Position', [15 130 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hk1_0');

%k2
uicontrol('Style', 'text', ...
 'Position', [45 165 30 30],...
 'String', 'k2',...
 'FontSize', 16, 'Tag', 'hk2');
uicontrol('Style', 'edit', ...
 'Position', [50 130 20 30],...
 'String', '2',... 
 'FontSize', 14, 'Tag', 'hk2_0');

%eps
uicontrol('Style', 'text', ...
 'Position', [80 165 30 30],...
 'String', 'ep',...
 'FontSize', 16, 'Tag', 'hep');
uicontrol('Style', 'edit', ...
 'Position', [85 130 20 30],...
 'String', '0.18',... 
 'FontSize', 14, 'Tag', 'hep_0');

%L
uicontrol('Style', 'text', ...
 'Position', [115 165 30 30],...
 'String', 'L',...
 'FontSize', 16, 'Tag', 'hL');
uicontrol('Style', 'edit', ...
 'Position', [120 130 20 30],...
 'String', '3.4',... 
 'FontSize', 14, 'Tag', 'hL_0');

%S
uicontrol('Style', 'text', ...
 'Position', [150 165 30 30],...
 'String', 'S',...
 'FontSize', 16, 'Tag', 'hS');
uicontrol('Style', 'edit', ...
 'Position', [155 130 20 30],...
 'String', '7',... 
 'FontSize', 14, 'Tag', 'hS_0');

%alpha
uicontrol('Style', 'text', ...
 'Position', [185 165 30 30],...
 'String', 'al',...
 'FontSize', 16, 'Tag', 'hal');
uicontrol('Style', 'edit', ...
 'Position', [190 130 20 30],...
 'String', '1',... 
 'FontSize', 14, 'Tag', 'hal_0');

%T
uicontrol('Style', 'text', ...
 'Position', [220 165 30 30],...
 'String', 'T',...
 'FontSize', 16, 'Tag', 'hT');
uicontrol('Style', 'edit', ...
 'Position', [225 130 20 30],...
 'String', '2',... 
 'FontSize', 14, 'Tag', 'hT_0');

%%{
%RelTol
uicontrol('Style', 'text', ...
 'Position', [185 100 62 30],...
 'String', 'RelT',...
 'FontSize', 16, 'Tag', 'hRelT');
uicontrol('Style', 'edit', ...
 'Position', [190 70 57 30],...
 'String', '1e-3',... 
 'FontSize', 14, 'Tag', 'hRelT_0');

%AbsTol
uicontrol('Style', 'text', ...
 'Position', [185 40 62 30],...
 'String', 'AbsT',...
 'FontSize', 16, 'Tag', 'hAbsT');
uicontrol('Style', 'edit', ...
 'Position', [190 10 57 30],... % 190 380 57 30
 'String', '1e-6',... 
 'FontSize', 14, 'Tag', 'hAbsT_0');
%%}
 
%start culculation
uicontrol('Style', 'pushbutton', ...
 'Position', [15 35 160 70],... %15 400 160 70
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

% getting values of constants 
%getting value of alpha
alpha_0 = get(handles.hal_0, 'String');
if strcmp(alpha_0, '')
    error('No value at alpha_0');
end
alpha_0 = str2double(alpha_0);
if isnan(alpha_0) 
    error('Some trash at alpha_0');
end 
values('alpha_0') = alpha_0;

%getting value of T
T_0 = get(handles.hT_0, 'String');
if strcmp(T_0, '')
    error('No value at T_0');
end
T_0 = str2double(T_0);
if isnan(T_0) 
    error('Some trash at T_0');
end 
values('T_0') = T_0;

%getting value of k1
k1_0 = get(handles.hk1_0, 'String');
if strcmp(k1_0, '')
    error('No value at k1_0');
end
k1_0 = str2double(k1_0);
if isnan(k1_0) 
    error('Some trash at k1_0');
end 
values('k1_0') = k1_0;

%getting value of k2
k2_0 = get(handles.hk2_0, 'String');
if strcmp(k2_0, '')
    error('No value at k2_0');
end
k2_0 = str2double(k2_0);
if isnan(k2_0) 
    error('Some trash at k2_0');
end 
values('k2_0') = k2_0;

%getting value of eps
eps_0 = get(handles.hep_0, 'String');
if strcmp(eps_0, '')
    error('No value at eps_0');
end
eps_0 = str2double(eps_0);
if isnan(eps_0) 
    error('Some trash at eps_0');
end 
values('eps_0') = eps_0;

%getting value of L
L_0 = get(handles.hL_0, 'String');
if strcmp(L_0, '')
    error('No value at L_0');
end
L_0 = str2double(L_0);
if isnan(L_0) 
    error('Some trash at L_0');
end 
values('L_0') = L_0;

%getting value of S
S_0 = get(handles.hS_0, 'String');
if strcmp(S_0, '')
    error('No value at S_0');
end
S_0 = str2double(S_0);
if isnan(S_0) 
    error('Some trash at S_0');
end 
values('S_0') = S_0;

%getting value of var
var_0 = get(handles.hvar_0, 'String');
if strcmp(var_0, '')
    error('No value at var_0');
end
var_0 = str2double(var_0);
if isnan(var_0) 
    error('Some trash at var_0');
end 
values('var_0') = var_0;

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
if isnan(AbsT) 
    error('Some trash at AbsT');
end 
values('AbsT') = AbsT;

StartCalculation(values)



