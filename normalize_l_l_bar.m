function [l_bar_n, l_n] = normalize_l_l_bar(l_bar, l)

global l_bar_min l_bar_max l_min l_max nm

l_bar_n = NaN(nm, 1);       l_n = NaN(nm, 1);

for i = 1:nm
    l_bar_n(i)      = (l_bar(i) - l_bar_min(i)) / ( l_bar_max(i) - l_bar_min(i) ); 
    l_n(i)          = (l(i) - l_min(i)) / (l_max(i) - l_min(i));
end
end
