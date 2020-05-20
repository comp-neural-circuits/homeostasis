function calc_ext_input_poisson_based(Ne, Ni, Ncells, n_ext, r_ex, r_ix, dt, wext_e, wext_i)
    ext_input = zeros(Ncells)
    ext_input[1:Ne] = sum(rand(n_ext, Ne) .< (r_ex * dt), 1) * wext_e
    ext_input[(1+Ne):Ncells] = sum(rand(n_ext, Ni) .< (r_ix * dt), 1) * wext_i
    return ext_input
end
