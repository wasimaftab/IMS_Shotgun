%% Map edge attributes to the closed interval [desired_lower_bound, desired_upper_bound]
function mapped_edge_attr = map_numbers_to_an_interval(edge_attr, current_lower_bound, current_upper_bound, desired_lower_bound, desired_upper_bound)
if current_lower_bound < current_upper_bound
    mapped_edge_attr = ((edge_attr-current_lower_bound)/(current_upper_bound-current_lower_bound))*(desired_upper_bound-desired_lower_bound)+desired_lower_bound;
else
    mapped_edge_attr = edge_attr;
    mapped_edge_attr(:) = desired_upper_bound;
end
