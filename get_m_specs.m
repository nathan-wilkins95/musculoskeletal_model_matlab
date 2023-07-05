function m_specs = get_m_specs(arm_data)

m_specs = struct;

muscle_data = struct(arm_data{"muscles"});

muscle_names = fieldnames(muscle_data);
for i = 1:length(muscle_names)
   muscle_name = muscle_names{i};
   params = struct(muscle_data.(string(muscle_name)){"params"});
   init_fiber_length = muscle_data.(string(muscle_name)){"init_fiber_length"};
   
   m_specs.(string(muscle_name)).params = params;     
   m_specs.(string(muscle_name)).init_fiber_length = init_fiber_length;
end
end