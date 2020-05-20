function calc_ext_input_corr_based(Ncells, Ne, ext_corr_input_series, n_ext, r_corr_poisson, dt, wext_e, wext_i, pattern_id, popmembers, pattern_loop_idx)
    ext_input = zeros(Ncells)
    not_popmembers = setdiff(1:Ne, popmembers[pattern_id, :])
    ext_input[not_popmembers] = sum(rand(n_ext, size(not_popmembers, 1)) .< (r_corr_poisson * dt), 1) * wext_e
    ext_input[popmembers[pattern_id, :]] = ext_corr_input_series[:, pattern_loop_idx] * wext_e
    ext_input[(1+Ne):Ncells] = sum(rand(n_ext, Ncells-Ne) .< (r_corr_poisson * dt), 1) * wext_i
    return ext_input
end
