function l_min = comp_l_min(alpha_0, l_0)
% copmute minimum muscle length
global alpha_1
    % adopted from MuscleFixedWidthPennationModel.cpp (OpenSim)
    h = l_0 * sin(alpha_0);
    p = sin(alpha_1);

    if alpha_0 > 0, l_min = h / p;
    else,           l_min = l_0 * 0.01;
    end

end
