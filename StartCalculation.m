function StartCalculation(values)
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

if (values('var_0') ~= 1) && (values('var_0') ~= 2)
    disp(values('var_0'));
    error('var is out of range')
end

alpha = values('alpha_0');
k1 = values('k1_0');
k2 = values('k2_0');
L = values('L_0');
S = values('S_0');
eps = values('eps_0');
T = values('T_0');
var = values('var_0');

%Реализация алгоритма
step = 0.001;
t = 0:step:T;

x2_func = @(k, t, psi0) psi0(1)*(1/(1+k)^2 - exp((1+k)*t)/(2*(1+k)^2) - ...
        exp(-(1+k)*t)/(2*(1+k)^2)) + ...
        psi0(2)*(exp((1+k)*t)/(2*(1+k)) - ...
        exp(-(1+k)*t)/(2*(1+k))) - ...
        (alpha/(2*(1+k))) * (1-exp(-(1+k)*t));
x1_func = @(k, t, psi0) psi0(1)*(t/(1+k)^2 - exp((1+k)*t)/(2*(1+k)^3) + ...
        exp(-(1+k)*t)/(2*(1+k)^3)) + ...
        psi0(2)*(exp((1+k)*t)/(2*(1+k)^2) + ...
        exp(-(1+k)*t)/(2*(1+k)^2) - 1/((1+k)^2)) - ...
        (alpha*t/(2*(1+k)) + alpha*exp(-(1+k)*t)/(2*(1+k)^2) - ...
        alpha/(2*(1+k)^2));
    
J_func = @(k, p, psi0) (psi0(1)/(1+k) + (psi0(2) - psi0(1)/(1+k))*...
    exp((1+k)*p) - alpha/2).^2 + (psi0(1)/(1+k) + ...
(psi0(2) - psi0(1)/(1+k))*exp((1+k)*p) - alpha/2)*alpha;

psi2_func = @(k, p, psi0) psi0(1)/(1+k) + (psi0(2) - psi0(1)/(1+k))*exp((1+k)*p);
% 1st var
if var == 1 %u1 >= 0
    
    % psi_2(T) > 0 => k = k1, x2_T = S-eps
    psi0 = GetPsi0(S-eps, k1, L, alpha, T);
    
    if (psi0(2) >= alpha/2) && (psi0(2) > psi0(1)/(1+k1)) % no switches
        disp("No switches");
        k = k1;
        J_res = integral(@(p)J_func(k1, p, psi0), 0, T);
        [t, y] = ode45(@(t, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:T, [0, 0]);
        x2 = y(:, 2);
        x1 = y(:, 1);
        u2 = ones(length(t), 1)*k1;
        [t, psi2] = ode45(@(t, y) -psi0(1) + y *(1+k), 0:step:T, psi0(2));
        u1 = psi2 - alpha/2;
    elseif (psi0(2) >= 0) && (psi0(2) < alpha/2) && ...
            (psi0(2) > psi0(1)/(1+k1)) %one switch at u1
        disp('Switch only at u1');
        k = k1;
        tSwitch = log((alpha/2 - psi0(1)/(1+k1))/...
            (psi0(2) - psi0(1)/(1+k1))) / (1+k1);
        if tSwitch < 0
            error('tSwitch < 0; var 1, ch.2')
        elseif tSwitch >= T
            error('tSwitch >= T, u1 = 0, task cannot be solved: x1 = 0, x1(T) > 0 - ?!')
        else
            k = k1;
            [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
            [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+k), tSwitch:step:T, alpha/2);
            psi2 = [psi2_1; psi2_2];
            t = [t_1; t_2];
            pos = find(t == tSwitch, 1) - 1;
            u1_2 = psi2_2 - alpha/2;
            u1 = [zeros(pos, 1); u1_2];
            u2 = ones(length(t), 1)*k1;
            J_res = integral(@(p)J_func(k1, p, psi0), tSwitch, T);
            
            [t_2, y] = ode45(@(t_2, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_2) - ...
                    alpha/2 - y(2)*(1+k)], tSwitch:step:T, [0, 0]);
            x2_2 = y(:, 2);
            x1_2 = y(:, 1);
            x2_1 = zeros(pos, 1);
            x2 = [x2_1; x2_2];
            
            x1_1 = zeros(pos, 1);
            x1 = [x1_1; x1_2];
        end
        
    elseif (psi0(2) > alpha/2) && (psi0(2) < psi0(1)/(1+k1)) % switch both at u1 and u2
        disp('psi2 go down, psi2 > alpha/2');
        tSwitch_u1 = log((alpha/2 - psi0(1)/(1+k1))/...
            (psi0(2) - psi0(1)/(1+k1))) / (1+k1);
        if tSwitch_u1 ~= T
            error('System stoped moving(');
        else
            psi2 = psi2_func(k1, t, psi0);
            u1 = psi2 - alpha/2;
            u2 = ones(length(t), 1)*k1;
            J_res = integral(@(p)J_func(k1, p, psi0), 0, T);
            
            if ~isempty(find(u1<0, 1))
                error('u1 < 0, var 1');
            end

            x2 = x2_func(k1, t, psi0);

            x1 = x1_func(k1, t, psi0);
        end
    elseif (psi0(2) < 0) && (psi0(2) > psi0(1)/(1+k2))
        disp('psi2 go up, psi2(0) < 0');
        tSwitch_u2 = log((-psi0(1)/(1+k2))/...
            (psi0(2) - psi0(1)/(1+k2))) / (1+k2);
        if tSwitch_u2 < 0
            error('tSwitch_u2 < 0')
        elseif tSwitch_u2 >= T
            error('tSwitch_u2 >= T')
        end
        tSwitch_u1 = log((alpha/2 - psi0(1)/(1+k1))/...
            (psi0(2) - psi0(1)/(1+k1))) / (1+k1);
        if tSwitch_u1 < 0
            error('tSwitch_u1 < 0')
        elseif tSwitch_u1 >= T
            error('tSwitch_u1 >= T, system hasnot started moving')
        end
        t1 = 0:step:tSwitch_u2;
        t2 = tSwitch_u2:step:tSwitch_u1;
        t3 = tSwitch_u1:step:T;
        t1 = t1';
        t2 = t2';
        t3 = t3';
        t = [t1; t2; t3];
        t = unique(t);
        pos_u2 = find(t == tSwitch_u2, 1);
        pos_u2 = pos_u2 - 1;
        pos_u1 = find(t == tSwitch_u1, 1);
        pos_u1 = pos_u1 - 1;
        
        psi_2_t3 = psi2_func(k1,t3, psi0);
        
        u1 = [zeros(pos_u2, 1); zeros(pos_u1 - pos_u2, 1); psi_2_t3 - alpha/2]; 
        u2 = [ones(pos_u2,1)* k2; ones(length(t) - pos_2, 1)*k1];
        J_res = integral(@(p)J_func(k1, p, psi0), t(pos_u1), T);
        
        x2_t3 = x2_func(k1, t3, psi0);
        x2 = [zeros(pos_u2, 1); zeros(pos_u1 - pos_u2, 1); x2_t3];
        
        x1_t3 = x1_func(k1, t3, psi0);
        x1 = [zeros(pos_u2, 1); zeros(pos_u1 - pos_u2, 1); x1_t3];
    elseif (psi0(2) <= alpha/2) && (psi0(2) > 0) && (psi0(2) < psi0(1)/(1+k1))
        error('no move has been started, task cannot be solved: x1 = 0, x1(T) > 0 - ?!')
    elseif (psi0(2) <= 0) && (psi0(1) >= 0)
        error('no move has been started, task cannot be solved: x1 = 0, x1(T) > 0 - ?!')
    elseif (psi0(2) <= 0) && (psi0(2) < psi0(1)/(1+k2))
        error('no move has been started, task cannot be solved: x1 = 0, x1(T) > 0 - ?!')
    elseif (psi0(2) < 0) && (psi0(2) == psi0(1)/(1+k2))
        error('no move has been started, task cannot be solved: x1 = 0, x1(T) > 0 - ?!')
    elseif (psi0(2) == 0) && (psi0(1) == 0)
        error('task cannot be solved with PMP')
    elseif (psi0(2) > 0) && (psi0(2) == psi0(1)/(1+k1))
        disp("psi2 == const > 0");
        psi_2 = psi0(2);
        
        u1 = ones(length(t), 1)*(psi_2 - alpha/2);
        u2 = ones(length(t), 1)*k1;
        x2 = (psi_2 - alpha/2)/(1+k1)*exp((1+k1)*t) - ...
            (psi_2 - alpha/2)/(1+k1)*exp(-(1+k1)*t);
        x1 = (psi_2 - alpha/2)/((1+k1)^2)*(exp((1+k1)*t) + exp(-(1+k1)*t)) - ...
            2*(psi_2 - alpha/2)/((1+k1)^2);
        
        J_res = integral(@(p)J_func(k1, p, psi0), 0, T);
    end
elseif var == 2 % u1 любое
    sum = k1 + k2;
    J_res = intmax('int64');
    for i = 1:2
        if i == 1
            psi0 = GetPsi0(S-eps, k1, L, alpha, T);
        else
            psi0 = GetPsi0(S+eps, k2, L, alpha, T);
        end
        if (psi0(2) < alpha/2) && (psi0(2) > 0)
            k = k2;
        else
            k = k1;
        end
        if ((-psi0(1)/(1+k))/(psi0(2) - psi0(1)/(1+k))) <= 0
            tSwitch_u2psi2 = -1;
        else
            tSwitch_u2psi2 = log((-psi0(1)/(1+k))/...
            (psi0(2) - psi0(1)/(1+k))) / (1+k);
        end
        if ((alpha-psi0(2)-psi0(1)/(1+k))/(psi0(2)-psi0(1)/(1+k))) <= 0
            tSwitch_u2x2 = -1;
        else
            tSwitch_u2x2 = log((alpha-psi0(2)-psi0(1)/(1+k))/...
            (psi0(2)-psi0(1)/(1+k))) /(1+k);
        end
        if ((tSwitch_u2psi2 <= 0) || (tSwitch_u2psi2 >= T)) && ((tSwitch_u2x2 <= 0) || (tSwitch_u2x2 >= T))
            disp("No switch at u2");
            J = integral(@(p)J_func(k, p, psi0), 0, T);
            disp("J = ");
            disp(J);
            if (i == 1) || (J < J_res)
                [t, y] = ode45(@(t, y)[y(2); psi0(1)/(1+k) + ...
                (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t) - ...
                alpha/2 - y(2)*(1+(sum-k))], 0:step:T, [0, 0]);
                x2 = y(:, 2);
                x1 = y(:, 1);
                [t, psi2] = ode45(@(t, y) -psi0(1) + y*(1+k), 0:step:T, psi0(2));
                u1 = psi2 - alpha/2;
                u2 = ones(length(t), 1)*k;
            end
            if i == 1
                J_res = J;
            elseif J < J_res
                J_res = J;
            end
        elseif ((tSwitch_u2psi2 <= 0) || (tSwitch_u2psi2 >= T)) && ...
            (tSwitch_u2x2 > 0) && (tSwitch_u2x2 < T)
            disp("Possible switch only by x2");
            tSwitch = tSwitch_u2x2;
            [t_1, y] = ode45(@(t_1, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_1) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:tSwitch, [0, 0]);
            x2_1 = y(:, 2);
            x1_1 = y(:, 1);
            [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+(sum-k)) + ...
                (psi0(2)-psi0(1)/(1+(sum-k)))*exp((1+(sum-k))*t_2) - ...
                alpha/2 - y(2)*(1+(sum-k))], tSwitch:step:T, [x1_1(end), x2_1(end)]);
            x2_2 = y(:, 2);
            x1_2 = y(:, 1);
            if (x2_1(end-1)*x2_2(2) > 0)
                disp('Nope, no switch');
                [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+k) + ...
                (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_2) - ...
                alpha/2 - y(2)*(1+k)], tSwitch:step:T, [x1_1(end), x2_1(end)]);
                x2_2 = y(:, 2);
                x1_2 = y(:, 2);
                if (x2_1(end-1)*x2_2(2) < 0)
                    disp('Task cannot be solved with this parameters, got to the point where x2 must be 0 all the time');
                    J = J_res;
                end
            else
                J_1 = integral(@(p)J_func(k, p, psi0), 0, tSwitch);
                J_2 = integral(@(p)J_func((sum-k), p, psi0), tSwitch, T);
                J = J_1 + J_2;
                disp("J = ");
                disp(J);
            end
            if (i == 1) || (J < J_res)
                t = [t_1; t_2];
                x2 = [x2_1 ; x2_2];
                x1 = [x1_1 ; x1_2];
                
                pos_u2 = find(t == tSwitch, 1) - 1;
                [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                
                [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum-k), tSwitch:step:T, psi2_1(end));
                t = [t_1; t_2];
                psi2 = [psi2_1; psi2_2];
                u1 = psi2 - alpha/2;
                u2_t1 = ones(pos_u2, 1)*k;
                u2_t2 = ones(length(t) - pos_u2, 1)*(sum - k);
                u2 = [u2_t1; u2_t2];
            end
            if i == 1
                J_res = J;
            elseif J < J_res
                J_res = J;
            end

        elseif ((tSwitch_u2psi2 > 0) && (tSwitch_u2psi2 < T)) && ...
            ((tSwitch_u2x2 <= 0) || (tSwitch_u2x2 >= T))
            disp("Switch only by psi2");
            tSwitch = tSwitch_u2psi2;
            [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                
            [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum - k), tSwitch:step:T, psi2_1(end));
            t = [t_1; t_2];
            if (psi2_1(end - 1)*psi2_2(2) > 0)
                [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+k), tSwitch:step:T, psi2_1(end));
                if (psi2_1(end - 1)*psi2_2(2) < 0)
                    disp('Task cannot be solved with this parameters, got to the point where psi2 must be 0 all the time');
                    J = J_res;
                end
            else
                J_1 = integral(@(p)J_func(k, p, psi0), 0, tSwitch);
                J_2 = integral(@(p)J_func((sum-k), p, psi0), tSwitch, T);
                J = J_1 + J_2;
                disp("J = ");
                disp(J);
            end
            if (i == 1) || (J < J_res)
                [t_1, y] = ode45(@(t_1, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_1) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:tSwitch, [0, 0]);
                x2_1 = y(:, 2);
                x1_1 = y(:, 1);
                [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+(sum-k)) + ...
                    (psi0(2)-psi0(1)/(1+(sum-k)))*exp((1+(sum-k))*t_2) - ...
                    alpha/2 - y(2)*(1+(sum-k))], tSwitch:step:T, [x1_1(end), x2_1(end)]);
                x2_2 = y(:, 2);
                x1_2 = y(:, 1);
                t = [t_1; t_2];
                x2 = [x2_1 ; x2_2];
                x1 = [x1_1 ; x1_2];
                [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum - k), tSwitch:step:T, psi2_1(end));
                pos_u2 = find(t == tSwitch, 1) - 1;
                psi2 = [psi2_1; psi2_2];
                t = [t_1; t_2];
                u1 = psi2 - alpha/2;
                u2_t1 = ones(pos_u2, 1)*k;
                u2_t2 = ones(length(t) - pos_u2, 1)*(sum - k);
                u2 = [u2_t1; u2_t2];
            end
            if i == 1
                J_res = J;
            elseif J < J_res
                J_res = J;
            end

        elseif ((tSwitch_u2psi2 > 0) && (tSwitch_u2psi2 < T)) && ...
            ((tSwitch_u2x2 > 0) && (tSwitch_u2x2 < T))
            if (psi0(2) > alpha)
                tSwitch = tSwitch_u2x2;
                [t_1, y] = ode45(@(t_1, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_1) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:tSwitch, [0, 0]);
                x2_1 = y(:, 2);
                x1_1 = y(:, 1);
                [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+(sum-k)) + ...
                    (psi0(2)-psi0(1)/(1+(sum-k)))*exp((1+(sum-k))*t_2) - ...
                    alpha/2 - y(2)*(1+(sum-k))], tSwitch:step:T, [x1_1(end), x2_1(end)]);
                x2_2 = y(:, 2);
                x1_2 = y(:, 1);
                if (x2_1(end-1)*x2_2(2) > 0)
                    disp('No, not even first switch');
                    [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_2) - ...
                    alpha/2 - y(2)*(1+k)], tSwitch:step:T, [x1_1(end), x2_1(end)]);
                    x2_2 = y(:, 2);
                    x1_2 = y(:, 2);
                    if (x2_1(end-1)*x2_2(2) < 0)
                        disp('Task cannot be solved with this parameters, got to the point where x2 must be 0 all the time');
                        tSwitch = -1;
                    end
                end
                if ((-psi0(1)/(1+(sum - k)))/(psi0(2) - psi0(1)/(1+(sum - k)))) <= 0
                    tSwitch_u2psi2 = -1;
                else
                    tSwitch_u2psi2 = log((-psi0(1)/(1+k))/...
                    (psi0(2) - psi0(1)/(1+k))) / (1+k);
                end
                tSwitch_2 = tSwitch_u2psi2;
            else
                tSwitch = tSwitch_u2psi2;
                [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                
                [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum - k), tSwitch:step:T, psi2_1(end));
                t = [t_1; t_2];
                if (psi2_1(end - 1)*psi2_2(2) > 0)
                    [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+k), tSwitch:step:T, psi2_1(end));
                    if (psi2_1(end - 1)*psi2_2(2) < 0)
                        disp('Task cannot be solved with this parameters, got to the point where psi2 must be 0 all the time');
                        tSwitch = -1;
                    end
                end
                if ((alpha-psi0(2)-psi0(1)/(1+(sum - k)))/(psi0(2)-psi0(1)/(1+(sum - k)))) <= 0
                    tSwitch_u2x2 = -1;
                else
                    tSwitch_u2x2 = log((alpha-psi0(2)-psi0(1)/(1+k))/...
                    (psi0(2)-psi0(1)/(1+k))) /(1+k);
                end
                tSwitch_2 = tSwitch_u2x2;
            end
            if tSwitch < 0
                J = intmax('int64');
            elseif (tSwitch_2 <= tSwitch) || (tSwitch_2 <= 0)
                disp("No 2nd switch at u2");
                J_1 = integral(@(p)J_func(k, p, psi0), 0, tSwitch);
                J_2 = integral(@(p)J_func((sum-k), p, psi0), tSwitch, T);
                J = J_1 + J_2;
                if (i == 1) || (J < J_res)
                    [t_1, y] = ode45(@(t_1, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_1) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:tSwitch, [0, 0]);
                    x2_1 = y(:, 2);
                    x1_1 = y(:, 1);
                    [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+(sum-k)) + ...
                        (psi0(2)-psi0(1)/(1+(sum-k)))*exp((1+(sum-k))*t_2) - ...
                        alpha/2 - y(2)*(1+(sum-k))], tSwitch:step:T, [x1_1(end), x2_1(end)]);
                    t = [t_1; t_2];
                    x2_2 = y(:, 2);
                    x1_2 = y(:, 1);
                    x2 = [x2_1 ; x2_2];
                    x1 = [x1_1 ; x1_2];
                    [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                
                    [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum - k), tSwitch:step:T, psi2_1(end));
                    psi2 = [psi2_1; psi2_2];
                    t = [t_1; t_2];
                    u1 = psi2 - alpha/2;
                    pos_u2 = find(t == tSwitch, 1) - 1;
                    u2_t1 = ones(pos_u2, 1)*k;
                    u2_t2 = ones(length(t) - pos_u2, 1)*(sum - k);
                    u2 = [u2_t1; u2_t2];
                end
                if i == 1
                    J_res = J;
                elseif J < J_res
                    J_res = J;
                end
                
            else 
                disp("We have 2 switches at u2");
                J_1 = integral(@(p)J_func(k, p, psi0), 0, tSwitch);
                J_2 = integral(@(p)J_func((sum-k), p, psi0), tSwitch, tSwitch_2);
                J_3 = integral(@(p)J_func(k, p, psi0), tSwitch_2, T);
                J = J_1 + J_2 + J_3;
                if (i == 1) || (J < J_res)
                    [t_1, y] = ode45(@(t_1, y) [y(2); psi0(1)/(1+k) + ...
                    (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_1) - ...
                    alpha/2 - y(2)*(1+k)], 0:step:tSwitch, [0, 0]);
                    x2_1 = y(:, 2);
                    x1_1 = y(:, 1);
                    [t_2, y] = ode45(@(t_2, y)[y(2); psi0(1)/(1+(sum-k)) + ...
                        (psi0(2)-psi0(1)/(1+(sum-k)))*exp((1+(sum-k))*t_2) - ...
                        alpha/2 - y(2)*(1+(sum-k))], tSwitch:step:tSwitch_2, [x1_1(end), x2_1(end)]);
                    x2_2 = y(:, 2);
                    x1_2 = y(:, 1);
                    [t_3, y] = ode45(@(t_3, y)[y(2); psi0(1)/(1+k) + ...
                        (psi0(2)-psi0(1)/(1+k))*exp((1+k)*t_3) - ...
                        alpha/2 - y(2)*(1+k)], tSwitch_2:step:T, [x1_2(end), x2_2(end)]);
                    x2_3 = y(:, 2);
                    x1_3 = y(:, 1);
                    t = [t_1; t_2; t_3];
                    x2 = [x2_1 ; x2_2; x2_3];
                    x1 = [x1_1 ; x1_2; x1_3];
                    if x2_2(end - 1)*x2_3(2) > 0
                        J = intmax('int64');
                    end
                    [t_1, psi2_1] = ode45(@(t_1, y) -psi0(1) + y*(1+k), 0:step:tSwitch, psi0(2));
                
                    [t_2, psi2_2] = ode45(@(t_2, y) -psi0(1) + y*(1+sum-k), tSwitch:step:tSwitch_2, psi2_1(end));
                    
                    [t_3, psi2_3] = ode45(@(t_3, y) -psi0(1) + y*(1+k), psi2_2(end));
                    psi2 = [psi2_1; psi2_2; psi2_3];
                    if psi2_2(end - 1)*psi2_3(2) > 0
                        J = intmax('int64');
                    end
                    t = [t_1; t_2; t_3];
                    u1 = psi2 - alpha/2;
                    pos1_u2 = find(t == tSwitch, 1) - 1;
                    pos2_u2 = find(t == tSwitch_2, 1) - 1;
                    u2_t1 = ones(pos1_u2, 1)*k;
                    u2_t2 = ones(pos2_u2 - pos1_u2, 1)*(sum - k);
                    u2_t3 = ones(length(t) - pos2_u2, 1)*k;
                    u2 = [u2_t1; u2_t2; u2_t3];
                end
                if i == 1
                    J_res = J;
                elseif J < J_res
                    J_res = J;
                end
                
            end
        end
    end
end
if J_res == intmax('int64')
    disp("There is no a decision(");
else
    %u1
    figure;
    plot(t, u1);
    title('1st optimal parameter');
    xlabel('t');
    ylabel('u_1' );

    %u2
    figure;
    plot(t, u2);
    title('2nd optimal parameter');
    xlabel('t');
    ylabel('u_2' );

    %x2
    figure;
    plot(t, x2);
    title('Speed');
    xlabel('t');
    ylabel('x_2');

    %x1
    figure;
    plot(t, x1);
    title('Trajectory');
    xlabel('t');
    ylabel('x_1');

    %psi2

    figure;
    plot(t, psi2);
    hold on;
    plot(t, ones(length(t),1)*alpha/2);
    title('Conjugate variable');
    xlabel('t');
    ylabel('\psi_2');
    legend('\psi_2', '\alpha*0.5')

    disp("J_res = ");
    disp(J_res);
end
end
