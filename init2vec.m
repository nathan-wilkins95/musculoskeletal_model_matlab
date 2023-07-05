 function [l_vec, q_vec] = init2vec(activations)
global m_specs

muscle_names = fieldnames(m_specs);

l_vec = NaN(6, 1);      q_vec = NaN(6, 1);
for i = 1:length(muscle_names)
    l_vec(i) = m_specs.(string(muscle_names{i})).init_fiber_length;
    q_vec(i) = activations.(string(muscle_names{i}))(1); % take initial conditions from OpenSim 
end
end