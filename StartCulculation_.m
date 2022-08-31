function StartCulculation(values)
if values('alpha_0') < 0
    error('alpha < 0')
end
if values('T_0') <= 0
    error('T <= 0')
end
if values('L_0') <= 0
    error('L <= 0')
end
if values('eps_0') <= 0
    error('eps <= 0')
end

if values('k1_0') >= values('k2_0')
    error('k2 <= k1')
end

if values('k1_0') <= 0
    error('k1 <= 0')
end

if values('k2_0') <= 0
    error('k2 <= 0')
end

if values('S_0') <= 0
    error('S <= 0')
end

if (values('var_0') ~= 1) || (values('var_0') ~= 2)
    error('var is out of range')
end

alpha = values('alpha_0');
k1 = values('k1_0');
k2 = values('k2_0');
L = values('L_0');
S = values('S_0');
eps = values('eps_0');
T = values('eps_0');
var = values('var_0');

RelT = values('RelT');
AbsT = values('AbsT');

step = 0.01;
if var == 1 %u1 >= 0
    % ������ ��� ������������
    for tSwitch = 0:step:T
        % ������ 1 \
        A1 = [; ];
        B1 = [; S-eps]; 
        % ������ 2
        % ������ 3
        % ������ 4
        % ������ 5
        % ������ 6
        % ������ 7
        % ������ 8
        % ������ 9
    end
end
end