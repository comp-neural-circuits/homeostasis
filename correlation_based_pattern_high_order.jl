function correlation_based_pattern_high_order(r_corr_pattern, n_ext, dt, n_pattern_loops, pattern_id, popmembers, p1, p2)
    # r_corr_pattern, firing rate of template neuron for correlation based pattern
    # n_ext, number of external neurons
    # dt, time step, 0.0001s
    # n_pattern_loop, number of loops for this pattern
    # n_assembly_cells, number of cells in each assembly
    # p1, p2, copy probabilities
    # ext_corr_input_series, input for assembly cells during the activation of pattern
    pop_idx = popmembers[pattern_id, :]
    n_assembly_cells = length(pop_idx)
    poisson_spike_template_mat = rand(n_pattern_loops, n_assembly_cells) .< (r_corr_pattern * dt)
    joint_poisson_spike_template = rand(n_pattern_loops) .< (r_corr_pattern * dt)
    ext_corr_input_mat = zeros(n_assembly_cells, n_ext, n_pattern_loops)
    ext_corr_input_series = zeros(n_assembly_cells, n_pattern_loops)
    for ii = 1:n_assembly_cells
        poisson_spike_idx_1 = find(poisson_spike_template_mat[:, ii] .== true)
        poisson_spike_idx_2 = find(joint_poisson_spike_template .== true)

        n_copy_1 = floor(Int, p1 * length(poisson_spike_idx_1))
        n_copy_2 = floor(Int, p2 * length(poisson_spike_idx_2))

        for jj = 1:n_ext
            poisson_spike_idx_jj_1 = poisson_spike_idx_1[randperm(length(poisson_spike_idx_1))[1:n_copy_1]]
            poisson_spike_idx_jj_2 = poisson_spike_idx_2[randperm(length(poisson_spike_idx_2))[1:n_copy_2]]
            ext_input_ii_jj = zeros(n_pattern_loops)
            ext_input_ii_jj[poisson_spike_idx_jj_1] = 1
            ext_input_ii_jj[poisson_spike_idx_jj_2] = 1
            ext_corr_input_mat[ii, jj, :] = ext_input_ii_jj
        end
    end
    ext_corr_input_series = squeeze(sum(ext_corr_input_mat, 2), 2)
    return ext_corr_input_series
end
