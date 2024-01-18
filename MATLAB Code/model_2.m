
%Example for getting a single solution
%[infected_A, total_A, infected_B, total_B] = solve_ODES_params(0.005, 2 * 10^2, 2 * 10^2, sol_interval);

%Left plot in paper Figure 10
%sol_interval = [0, 1200];
%plot_infection_over_m(0, 0.005, sol_interval);
%Middle plot in paper Figure 10
sol_interval = [0, 800];
%plot_infection_boundary_fixed_IC(0, 0.016, 2 * 10^2, 1.6, 12, sol_interval);
%Right plot in paper Figure 10
%plot_infection_boundary_fixed_CCs(0.001, 0.005, 200, 200, 8, 12, sol_interval);

[infected_A, total_A, infected_B, total_B] = solve_ODES_params(0.005, 2 * 10^2, 2 * 10^2, sol_interval);
infected_A
infected_B
total_A - infected_A
total_B - infected_B

function evaluateMosquitoODEs(m, axes, interval)
    
    %parameters from original paper
    b_f = 0.5;
    b_m = 0.5;
    sigma = 1;
    phi_u = 13;
    phi_w = 11;
    v_w = 0.95;
    v_u = 0.05;
    psi = 1/8.75;
    mu_a = 0.02;
    mu_fu = 1/17.5;
    mu_fw = 1/15.8;
    mu_mu = 1/10.5;
    mu_mw = 1/10.5;
    K_a_A = 2 * 10^5;
    K_a_B = 2 * 10^5;
    % K_a_B = 4 * 10^5;

    function val = eta_u(x, y, hab_K_a)
        val = phi_u * (1 - (x + y)/hab_K_a);
    end
    
    function val = eta_w(x, y, hab_K_a)
        val = phi_w * (1 - (x + y)/hab_K_a);
    end

    function val = m_u(x, y)
        if (x + y == 0) 
            val = 0;
        else
            val = x/(x + y);
        end
    end

    function val = m_w(x, y)
        if (x + y == 0) 
            val = 0;
        else
            val = y/(x + y);
        end
    end

    % The aquatic stages will have no crosstalk, but all adult stages will.
    ODESystem = @(t, Y) [
        eta_u(Y(1), Y(2), K_a_A) * Y(5) + v_u * eta_w(Y(1), Y(2), K_a_A) * Y(6) - (mu_a + psi) * Y(1);  %Y(1) is A_u_A
        v_w * eta_w(Y(1), Y(2), K_a_A) * Y(6) - (mu_a + psi) * Y(2);                                    %Y(2) is A_w_A
        b_f * psi * Y(1) - (sigma + mu_fu) * Y(3) + m * (Y(12) - Y(3));                                 %Y(3) is F_u_A
        b_f * psi * Y(2) - (sigma + mu_fw) * Y(4) + m * (Y(13) - Y(4));                                 %Y(4) is F_w_A
        sigma * m_u(Y(7), Y(8)) * Y(3) - mu_fu * Y(5) + m * (Y(14) - Y(5));                             %Y(5) is F_pu_A
        sigma * Y(4) - mu_fw * Y(6) + m * (Y(15) - Y(6));                                               %Y(6) is F_pw_A
        b_m * psi * Y(1) - mu_mu * Y(7) + m * (Y(16) - Y(7));                                           %Y(7) is M_u_A
        b_m * psi * Y(2) - mu_mw * Y(8) + m * (Y(17) - Y(8));                                           %Y(8) is M_w_A
        sigma * m_w(Y(7), Y(8)) * Y(3) - mu_fu * Y(9) + m * (Y(18) - Y(9));                             %Y(9) is F_ps_A
        
        eta_u(Y(10), Y(11), K_a_B) * Y(14) + v_u * eta_w(Y(10), Y(11), K_a_B) * Y(15) - (mu_a + psi) * Y(10);  %Y(10) is A_u_B
        v_w * eta_w(Y(10), Y(11), K_a_B) * Y(15) - (mu_a + psi) * Y(11);                                       %Y(11) is A_w_B
        b_f * psi * Y(10) - (sigma + mu_fu) * Y(12) + m * (Y(3) - Y(12));                                      %Y(12) is F_u_B
        b_f * psi * Y(11) - (sigma + mu_fw) * Y(13) + m * (Y(4) - Y(13));                                      %Y(13) is F_w_B
        sigma * m_u(Y(16), Y(17)) * Y(12) - mu_fu * Y(14) + m * (Y(5) - Y(14));                                %Y(14) is F_pu_B
        sigma * Y(13) - mu_fw * Y(15) + m * (Y(6) - Y(15));                                                    %Y(15) is F_pw_B
        b_m * psi * Y(10) - mu_mu * Y(16) + m * (Y(7) - Y(16));                                                %Y(16) is M_u_B
        b_m * psi * Y(11) - mu_mw * Y(17) + m * (Y(8) - Y(17));                                                %Y(17) is M_w_B
        sigma * m_w(Y(16), Y(17)) * Y(12) - mu_fu * Y(18) + m * (Y(9) - Y(18));                                %Y(18) is F_ps_B
    ];
    
    % starting populations
    initialConditions = [
        K_a_A       ,... % A_u_A
        3 * 10^5    ,... % A_w_A
        3 * 10^5    ,... % F_u_A
        5 * 10^5    ,... % F_w_A
        1 * 10^5    ,... % F_pu_A
        1 * 10^5    ,... % F_pw_A
        4 * 10^5    ,... % M_u_A
        0 * 10^5    ,... % M_w_A
        0 * 10^5    ,... % F_ps_A
        K_a_B       ,... % A_u_B
        0 * 10^5    ,... % A_w_B
        8 * 10^5    ,... % F_u_B
        0 * 10^5    ,... % F_w_B
        1 * 10^5    ,... % F_pu_B
        0 * 10^5    ,... % F_pw_B
        4 * 10^5    ,... % M_u_B
        0 * 10^5    ,... % M_w_B
        0 * 10^5    ,... % F_ps_B
    ];
    
    solution = ode45(ODESystem, interval, initialConditions);

    
    function [props_w_A, props_w_B, props_u_A, props_u_B] = getInfectionProportions(solution, tVals)
        props_u_A = [];
        props_u_B = [];
        props_w_A = [];
        props_w_B = [];
        ws_A = @(t)deval(solution, t, [2, 4, 6, 8]);
        ws_B = @(t)deval(solution, t, [11, 13, 15, 17]);
        all_sol_A = @(t)deval(solution, t, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        all_sol_B = @(t)deval(solution, t, [10, 11, 12, 13, 14, 15, 16, 17, 18]);
        
        for i=1:length(tVals)
            num_w_A = sum(ws_A(tVals(i)));
            num_w_B = sum(ws_B(tVals(i)));
            num_tot_A = sum(all_sol_A(tVals(i)));
            num_tot_B = sum(all_sol_B(tVals(i)));
            prop_w_A = num_w_A/num_tot_A;
            prop_w_B = num_w_B/num_tot_B;
            props_w_A = [props_w_A, prop_w_A];
            props_w_B = [props_w_B, prop_w_B];
            props_u_A = [props_u_A, 1 - prop_w_A];
            props_u_B = [props_u_B, 1 - prop_w_B];
        end
    end
    
    tVals = linspace(interval(1), interval(2), 1000);
    
    [props_w_A, props_w_B, props_u_A, props_u_B] = getInfectionProportions(solution, tVals);

    hold on;
    plot(axes, tVals, props_u_A, tVals, props_w_A, tVals, props_u_B, tVals, props_w_B);
    xlabel(axes, "Time (Days)", "FontName", "Times New Roman");
    ylabel(axes, "Proportion out of Total Mosquitoes", "FontName", "Times New Roman");
    title(axes, "Infection Densities Over Time", "FontName", "Times New Roman");
    legend_keys = {'Uninfected, Habitat A', 'Infected, Habitat A', 'Uninfected, Habitat B', 'Infected, Habitat B'};
    legend(axes, legend_keys, 'Location', 'best');
    hold off

end

function updateSolution(event, axes, interval) 
    evaluateMosquitoODEs(event.Value, axes, interval);
end
















function [infected_A, total_A, infected_B, total_B] = solve_ODES_params(m, CC_A, CC_B, interval, init_infected)
    %constant parameters from original paper
    b_f = 0.5;
    b_m = 0.5;
    sigma = 1;
    phi_u = 13;
    phi_w = 11;
    v_w = 0.95;
    v_u = 0.05;
    psi = 1/8.75;
    mu_a = 0.02;
    mu_fu = 1/17.5;
    mu_fw = 1/15.8;
    mu_mu = 1/10.5;
    mu_mw = 1/10.5;
    K_a_A = CC_A;
    K_a_B = CC_B;

    if(nargin == 4)
        init_infected = 5 * CC_A;
    end

    function val = eta_u(x, y, hab_K_a)
        val = phi_u * (1 - (x + y)/hab_K_a);
    end
    
    function val = eta_w(x, y, hab_K_a)
        val = phi_w * (1 - (x + y)/hab_K_a);
    end

    function val = m_u(x, y)
        if (x + y == 0) 
            val = 0;
        else
            val = x/(x + y);
        end
    end

    function val = m_w(x, y)
        if (x + y == 0) 
            val = 0;
        else
            val = y/(x + y);
        end
    end

    % The aquatic stages will have no crosstalk, but all adult stages will.
    ODESystem = @(t, Y) [
        eta_u(Y(1), Y(2), K_a_A) * Y(5) + v_u * eta_w(Y(1), Y(2), K_a_A) * Y(6) - (mu_a + psi) * Y(1);  %Y(1) is A_u_A
        v_w * eta_w(Y(1), Y(2), K_a_A) * Y(6) - (mu_a + psi) * Y(2);                                    %Y(2) is A_w_A
        b_f * psi * Y(1) - (sigma + mu_fu) * Y(3) + m * (Y(12) - Y(3));                                 %Y(3) is F_u_A
        b_f * psi * Y(2) - (sigma + mu_fw) * Y(4) + m * (Y(13) - Y(4));                                 %Y(4) is F_w_A
        sigma * m_u(Y(7), Y(8)) * Y(3) - mu_fu * Y(5) + m * (Y(14) - Y(5));                             %Y(5) is F_pu_A
        sigma * Y(4) - mu_fw * Y(6) + m * (Y(15) - Y(6));                                               %Y(6) is F_pw_A
        b_m * psi * Y(1) - mu_mu * Y(7) + m * (Y(16) - Y(7));                                           %Y(7) is M_u_A
        b_m * psi * Y(2) - mu_mw * Y(8) + m * (Y(17) - Y(8));                                           %Y(8) is M_w_A
        sigma * m_w(Y(7), Y(8)) * Y(3) - mu_fu * Y(9) + m * (Y(18) - Y(9));                             %Y(9) is F_ps_A
        
        eta_u(Y(10), Y(11), K_a_B) * Y(14) + v_u * eta_w(Y(10), Y(11), K_a_B) * Y(15) - (mu_a + psi) * Y(10);  %Y(10) is A_u_B
        v_w * eta_w(Y(10), Y(11), K_a_B) * Y(15) - (mu_a + psi) * Y(11);                                       %Y(11) is A_w_B
        b_f * psi * Y(10) - (sigma + mu_fu) * Y(12) + m * (Y(3) - Y(12));                                      %Y(12) is F_u_B
        b_f * psi * Y(11) - (sigma + mu_fw) * Y(13) + m * (Y(4) - Y(13));                                      %Y(13) is F_w_B
        sigma * m_u(Y(16), Y(17)) * Y(12) - mu_fu * Y(14) + m * (Y(5) - Y(14));                                %Y(14) is F_pu_B
        sigma * Y(13) - mu_fw * Y(15) + m * (Y(6) - Y(15));                                                    %Y(15) is F_pw_B
        b_m * psi * Y(10) - mu_mu * Y(16) + m * (Y(7) - Y(16));                                                %Y(16) is M_u_B
        b_m * psi * Y(11) - mu_mw * Y(17) + m * (Y(8) - Y(17));                                                %Y(17) is M_w_B
        sigma * m_w(Y(16), Y(17)) * Y(12) - mu_fu * Y(18) + m * (Y(9) - Y(18));                                %Y(18) is F_ps_B
    ];
    
    % starting populations
    initialConditions = [
        K_a_A    ,... % A_u_A
        0 * K_a_A    ,... % A_w_A
        1.5 * K_a_A    ,... % F_u_A
        init_infected    ,... % F_w_A
        0.5 * K_a_A    ,... % F_pu_A
        0 * K_a_A    ,... % F_pw_A
        1.5 * K_a_A    ,... % M_u_A
        0 * K_a_A    ,... % M_w_A
        0 * K_a_A    ,... % F_ps_A
        K_a_B       ,... % A_u_B
        0 * K_a_B    ,... % A_w_B
        1.5 * K_a_B    ,... % F_u_B
        0 * K_a_B    ,... % F_w_B
        0.5 * K_a_B    ,... % F_pu_B
        0 * K_a_B    ,... % F_pw_B
        1.5 * K_a_B    ,... % M_u_B
        0 * K_a_B    ,... % M_w_B
        0 * K_a_B    ,... % F_ps_B
    ];
    
    solution = ode45(ODESystem, interval, initialConditions);
    all_sol_A = deval(solution, interval(2), [1, 2, 3, 4, 5, 6, 7, 8, 9]);
    all_sol_B = deval(solution, interval(2), [10, 11, 12, 13, 14, 15, 16, 17, 18]);
    
    %all_sol_A
    %all_sol_B

    infected_A = all_sol_A(2) + all_sol_A(4) + all_sol_A(6) + all_sol_A(8);
    infected_B = all_sol_B(2) + all_sol_B(4) + all_sol_B(6) + all_sol_B(8);
    total_A = sum(all_sol_A);
    total_B = sum(all_sol_B);

end


function plot_infection_over_m(m_min, m_max, sol_interval)
    fig = uifigure('Name', 'Solution Curves');
    axes = uiaxes(fig);
    
    m_num = 100;
    m_vals = linspace(m_min, m_max, m_num);
    infected_A_vals = zeros(m_num, 1);
    infected_B_vals = zeros(m_num, 1);
    
    tic;
    for i = 1:m_num
        [infected_A, total_A, infected_B, total_B] = solve_ODES_params(m_vals(i), 2 * 10^2, 2 * 10^2, sol_interval);
        infected_A_vals(i) = infected_A / total_A;
        infected_B_vals(i) = infected_B / total_B;
    end
    toc;

    hold on;
    plot(axes, m_vals, infected_A_vals, m_vals, infected_B_vals);
    xlabel(axes, "m");
    ylabel(axes, "Infection Density");
    title(axes, "Equilibrium Infection Densities vs Migration Parameter");
    legend_keys = {'Habitat A', 'Habitat B'};
    legend(axes, legend_keys, 'Location', 'best');
    hold off
end


function plot_infection_boundary_fixed_IC(m_min, m_max, CC_A, max_scale, iterations, sol_interval)
    fig = uifigure('Name', 'Solution Curves');
    axes = uiaxes(fig);
    
    m_num = 200;
    m_vals = linspace(m_min, m_max, m_num);
    boundary_vals = zeros(m_num, 1);
    
    tic;
    for i = 1:m_num
        % ratio is CC_B / CC_A
        ratio = 0.001;
        [~, ~, infected_B, total_B] = solve_ODES_params(m_vals(i), CC_A, ratio * CC_A, sol_interval);
        infected_B_prop = infected_B / total_B;
        if(infected_B_prop < 0.5)
            boundary_vals(i) = ratio;
            continue
        end
        ratio = max_scale;
        [~, ~, infected_B, total_B] = solve_ODES_params(m_vals(i), CC_A, ratio * CC_A, sol_interval);
        infected_B_prop = infected_B / total_B;
        if(infected_B_prop > 0.5)
            disp('No habitat proportion was found for the following value of m:')
            m_vals(i)
            continue
        end
        lower_bound = 0;
        upper_bound = max_scale;
        for j = 1:iterations
            ratio = (upper_bound - lower_bound) / 2 + lower_bound;
            [~, ~, infected_B, total_B] = solve_ODES_params(m_vals(i), CC_A, ratio * CC_A, sol_interval);
            infected_B_prop = infected_B / total_B;
            if(infected_B_prop > 0.5)
                % Successfully invaded Habitat B. Set lower bound to ratio
                lower_bound = ratio;
            else
                % Failed to invade Habitat B. Set upper bound to ratio
                upper_bound = ratio;
            end
        end
        ratio = (upper_bound - lower_bound) / 2 + lower_bound;
        boundary_vals(i) = ratio;
    end
    toc;

    hold on;
    plot(axes, m_vals, boundary_vals);
    xlabel(axes, "m");
    ylabel(axes, "K_B/K_A");
    title(axes, "Successful Infection Boundary as a Function of m and K_B/K_A");
    legend_keys = {'Infection Boundary'};
    legend(axes, legend_keys, 'Location', 'best');
    hold off
end


function plot_infection_boundary_fixed_CCs(m_min, m_max, CC_A, CC_B, max_scale, iterations, sol_interval)
    fig = uifigure('Name', 'Solution Curves');
    axes = uiaxes(fig);
    
    m_num = 100;
    m_vals = linspace(m_min, m_max, m_num);
    boundary_vals = zeros(m_num, 1);

    initial_uninfected = (1 + 1.5 + 0.5 + 1.5) * CC_A + (1 + 1.5 + 0.5 + 1.5) * CC_B;
    
    tic;
    for i = 1:m_num
        % ratio is initial infected / initial uninfected
        ratio = max_scale;
        [~, ~, infected_B, total_B] = solve_ODES_params(m_vals(i), CC_A, CC_B, sol_interval, initial_uninfected * ratio);
        infected_B_prop = infected_B / total_B;
        if(infected_B_prop < 0.5)
            disp('No initial infection was found for the following value of m:')
            m_vals(i)
            continue
        end
        lower_bound = 0;
        upper_bound = max_scale;
        for j = 1:iterations
            ratio = (upper_bound - lower_bound) / 2 + lower_bound;
            [~, ~, infected_B, total_B] = solve_ODES_params(m_vals(i), CC_A, CC_B, sol_interval, initial_uninfected * ratio);
            infected_B_prop = infected_B / total_B;
            if(infected_B_prop > 0.5)
                % Successfully invaded Habitat B. Set upper bound to ratio
                upper_bound = ratio;
            else
                % Failed to invade Habitat B. Set lower bound to ratio
                lower_bound = ratio;
            end
        end
        ratio = (upper_bound - lower_bound) / 2 + lower_bound;
        boundary_vals(i) = ratio;
    end
    toc;

    hold on;
    plot(axes, m_vals, boundary_vals);
    xlabel(axes, "m");
    ylabel(axes, "initial infected / initial uninfected");
    title(axes, "Infection Boundary for Given m and Initial Infection Density");
    legend_keys = {'Infection Boundary'};
    legend(axes, legend_keys, 'Location', 'best');
    hold off
end
