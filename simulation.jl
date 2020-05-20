function simulation(weights, connections, popmembers)
    println("setting up parameters")
    Ne, Ni, wee, wei, wie, wii, prob = initialize_network()
    Ncells = Ne+Ni
    wext_e_ini = 0.78 # external to e synaptic strength at baseline
    wext_e_1 = 0.78
    wext_e_2 = 0.78 * 0.92 # external to e synaptic strength at md
    wext_i_ini = 0.85 # external to i synaptic strength at baseline
    wext_i_1 = 0.85
    wext_i_2 = 0.85 * 0.85 # external to i synaptic strength at md
    t_total = 200 # duration of early md in s, 200s
    w_min_ee = 0 # minimum synaptic strength of e to e
    w_max_ee = 1.2 # maximum synaptic strength of e to e
    w_min_ie = 0 # minimum synaptic strength of i to e
    w_max_ie = 6 # maximum synaptic strength of i to e
    w_min_ei = 0 # minimum synaptic strength of e to i
    w_max_ei = 1 # maximum synaptic strength of e to i

    # simulation
    dt = 0.0001 # integration timestep
    T = 2000 # simulation time in s
    Nsteps = round(Int,T/dt) # number of time steps
    
    # membrane dynamics
    v_exc = 0.000 # e reversal potential
    v_rest = -0.070 # resting potential
    v_inh = -0.080 # i reversal potential
    tau_m_e = 0.020 # e membrane time constant
    tau_m_i = 0.010 # i membrane time constant

    # receptors
    tau_ampa = 0.005 # time constant for ampa
    tau_gaba = 0.010 # time constant for gaba
    tau_nmda = 0.100 # time constant for nmda
    alpha = 0.5 # weights for ampa and nmda, determine the overall exc conductance

    # firing threshold
    v_thr = ones(Ncells) * (-0.050) # firing threshold
    v_spike_rest = -0.07 # resetting membrane potential after spike
    t_ref = 0.005 # refractory period
    t_allow_spike = zeros(Ncells) # allowed spike timing
    
    # conductance based synapse
    v = ones(Ncells) * (-0.070) # membrane potential
    g_inh = zeros(Ncells) # inhibitory conductance
    g_ampa = zeros(Ncells) # ampa conductance
    g_nmda = zeros(Ncells) # nmda conductance
    g_exc = zeros(Ncells) # exc conductance

    # target firing rate
    rho_0_e = ones(Ne) * 5 # target firing rate of excitatory neurons
    rho_0_i = ones(Ni) * 13 # target firing rate of inhibitory neurons
    
    # triplet stdp parameters for E-E
    b_learning_active_e_e = false
    tau_p = 0.0168 # time constant of r1
    tau_m = 0.0337 # time constant of r2
    tau_x = 0.101 # time constant of o1
    tau_y = 0.114 # time constant of o2
    r1 = zeros(Ncells) # presynaptic spike detector (faster decay)
    r2 = zeros(Ncells) # presynaptic spike detector (slower decay)
    o1 = zeros(Ncells) # postsynaptic spike detector (faster decay)
    o2 = zeros(Ncells) # postsynaptic spike detector (slower decay)
    A2_p = zeros(Ncells) # A2 positive
    A2_m = ones(Ncells) * (0.0071) # A2 negative
    A2_m_initial = 0.0071 # A2 negative initial
    A3_p = ones(Ncells) * (0.0065) # A3 positive
    A3_m = zeros(Ncells) # A3 negative09
    pre_spike_time = -100*ones(Ncells) # last spike timing

    # metaplasticity e-e
    b_metaplasticity_ee = false
 
    # istdp
    b_learning_active_istdp = false
    tau_istdp = 0.02 # width of istdp curve
    eta_istdp = 1 # istdp learning rate
    trace_istdp = zeros(Ncells) # low-pass filtered spike train for istdp

    # heterosynaptic plasticity
    b_active_normalization = true
    eta_norm = ones(Ne) # learning rate
    tau_norm = 10 # time constant

    sum_wee_initial = zeros(Ne) # intial weight sum
    sum_wee_var = zeros(Ne) # variable weight sum
    n_steps_normalization = Int(1/dt) # normalization frequency
    Nee = zeros(Int,Ne) # number of excitatory connections to postsynaptic excitatory neuron
    for cc = 1:Ne 
        for dd = 1:Ne
            sum_wee_var[cc] += weights[dd,cc]
            if connections[dd,cc] > 0 
                sum_wee_initial[cc] += connections[dd,cc] + 0.08
                Nee[cc] += 1
            end
        end
    end

    # synaptic scaling
    b_synaptic_scaling_e_e = false
    tau_r = 20 # time constant for postsynaptic activity estimator
    r_t_e_ss = ones(Ne) * 5 * 20 # postsynaptic excitatory activity estimator
    r_t_i_ss = ones(Ni) * 13 * 20 # postsynaptic inhibitory activity estimator
    eta_ss_e = 1 # e-e learning rate
    eta_ss_i = 1 # e-i learning rate
    tau_ss = t_total # time constant

    # intrinsic plasticity
    b_intrinsic_plasticity = false
    eta_ip_e = 0.00025 # learning rate
    eta_ip_i = 0.00025 # learning rate
    tau_ip = t_total # time constant
    
    # external poisson input
    n_ext = 1000 # number of external neurons
    r_exc_e = 5 # external input rate to e (Hz)
    r_exc_i = 5 # external input rate to i (Hz)

    # external poisson input
    b_only_poisson = true
    b_poisson = true
    b_poisson_1st_time = true
    n_start_poisson_loops = Int(1/dt) # start to impose input pattern after 5s
    poisson_loop_idx = 0
    
    # external pattern input
    Npop = size(popmembers,1) # number of assemblies
    n_patterns = Npop
    pattern_id = 1
    n_pattern_loops = Int(1/dt) # duration of pattern, 1s
    pattern_loop_idx = 1
    n_gap_loops = Int(3/dt) # gap between different patterns, 1s
    
    # externalorthogonal rate based input
    b_active_pattern = false
    
    # external correlated based input
    b_correlated_input = true
    b_new_correlated_pattern = true
    b_shuffle_pattern = false
    r_corr_pattern = r_exc_e # 20 Hz, so that the firing rate of each neuron in one assembly is 20 Hz
    r_corr_poisson = r_exc_e
    ext_corr_input_series = zeros(size(popmembers, 1), n_pattern_loops)
    p1 = 0.4
    p2 = 0.6

    # save information
    b_save_output = true
    loop_counter = 1 # loop counter in every second, from 1 to 10000
    n_save_steps_spikes = Int(1/dt) # save spike info every second
    spike_time_mat = zeros((n_save_steps_spikes, Ncells)) # spike train matrix
    weights_ss = zeros(Ncells,Ncells) # weights change in 1 s induced by synaptic scaling
    weights_ltp = zeros(Ncells,Ncells) # weights change in 1 s induced by LTP
    weights_ltd = zeros(Ncells,Ncells) # weights change in 1 s induced by LTD
    v_me = zeros((n_save_steps_spikes)) # membrane potential of excitatory neuron
    v_mi = zeros((n_save_steps_spikes)) # membrane potential of inhibitory neuron

    # correlated based input pattern
    for tt = 1:Nsteps

        if mod(tt,Nsteps/100) == 1
            @printf("\r%d%%",round(Int,100*tt/Nsteps)) # print out the progress percentage
        end

        if b_active_normalization
            if mod(tt,n_steps_normalization) == 0 #excitatory synaptic normalization
                for cc = 1:Ne
                    sum_wee = 0.
                    for dd = 1:Ne
                        sum_wee += weights[dd,cc]
                    end
                    
                    if sum_wee > sum_wee_initial[cc]
                        for dd = 1:Ne
                            if connections[dd,cc] > 0.
                                weights[dd,cc] -= eta_norm[cc] * (sum_wee-sum_wee_initial[cc])/Nee[cc]
                                if weights[dd,cc] < w_min_ee
                                    weights[dd,cc] = w_min_ee
                                elseif weights[dd,cc] > w_max_ee
                                    weights[dd,cc] = w_max_ee
                                end
                            end
                        end
                    end
                end
            end
        end # end normalization

        if tt > 30000
            b_learning_active_istdp = true
        end

        if 100000 < tt < (t_total * 10000 + 100000)
            b_learning_active_e_e = true
            b_metaplasticity_ee = true
	    b_synaptic_scaling_e_e = true
	    b_intrinsic_plasticity = true
            wext_e_1 = wext_e_ini - (wext_e_ini - wext_e_2) * (tt - 100000)/(t_total * 10000)
            wext_i_1 = wext_i_ini - (wext_i_ini - wext_i_2) * (tt - 100000)/(t_total * 10000)
        end
            
        if tt >= (t_total * 10000 + 100000)
            wext_e_1 = wext_e_2
            wext_i_1 = wext_i_2
        end

        if b_save_output
            if mod(tt,n_save_steps_spikes) == 0
                file = matopen("data/weight_matrix_$(tt).mat", "w")
                write(file, "weight_matrix", weights)
                close(file)
                file = matopen("data/weight_matrix_ss_$(tt).mat", "w")
                write(file, "weight_matrix_ss", weights_ss)
                close(file)
                file = matopen("data/weight_matrix_ltp_$(tt).mat", "w")
                write(file, "weight_matrix_ltp", weights_ltp)
                close(file)
                file = matopen("data/weight_matrix_ltd_$(tt).mat", "w")
                write(file, "weight_matrix_ltd", weights_ltd)
                close(file)
                file = matopen("data/thr_info_$(tt).mat", "w")
                write(file, "thr_info", v_thr)
                close(file)
                file = matopen("data/A2_$(tt).mat", "w")
                write(file, "A2", A2_m)
                close(file)
                weights_ss = zeros(Ncells,Ncells)
                weights_ltp = zeros(Ncells,Ncells)
                weights_ltd = zeros(Ncells,Ncells)
            end

            if loop_counter == (n_save_steps_spikes + 1)
                file = matopen("data/spike_matrix_$(tt-1).mat", "w")
                write(file, "spike_matrix", spike_time_mat)
                close(file)
                loop_counter = 1 # reset loop counter
                v_me = zeros((n_save_steps_spikes))  # reset spike train matrix
                v_mi = zeros((n_save_steps_spikes))  # reset spike train matrix
                spike_time_mat = zeros((n_save_steps_spikes, Ncells))  # reset spike train matrix
            end
        end

        t = dt*tt # real time
        t_prev = dt*(tt-1)
        
        # update external input, for poisson input
        if b_only_poisson
            ext_input = calc_ext_input_poisson_based(Ne, Ni, Ncells, n_ext, r_exc_e, r_exc_i, dt, wext_e_1, wext_i_1)
        else
            if b_poisson
                ext_input = calc_ext_input_poisson_based(Ne, Ni, Ncells, n_ext, r_exc_e, r_exc_i, dt, wext_e_1, wext_i_1)
                poisson_loop_idx += 1
                if b_poisson_1st_time
                    if poisson_loop_idx > n_start_poisson_loops
                        b_active_pattern = true
                        b_poisson_1st_time = false
                        b_poisson = false
                        poisson_loop_idx = 0
                    end
                else
                    if poisson_loop_idx > n_gap_loops
                        b_active_pattern = true
                        b_poisson = false
                        poisson_loop_idx = 0
                    end
                end
            end

            if b_active_pattern

                if b_correlated_input
                    if b_new_correlated_pattern
                        ext_corr_input_series = correlation_based_pattern_high_order(r_corr_pattern, n_ext, dt, n_pattern_loops, pattern_id, popmembers, p1, p2)
                        b_new_correlated_pattern = false
                    end
                    ext_input = calc_ext_input_corr_based(Ncells, Ne, ext_corr_input_series, n_ext, r_corr_poisson, dt, wext_e_1, wext_i_1, pattern_id, popmembers, pattern_loop_idx)
                end

                pattern_loop_idx += 1
                if pattern_loop_idx > n_pattern_loops
                    pattern_loop_idx = 1
                    b_new_correlated_pattern = true
                    b_active_pattern = false
                    b_poisson = true
                    if b_shuffle_pattern
                        pattern_id = sample(1: n_patterns)
                    else
                        pattern_id = mod(pattern_id + 1, n_patterns)
                        if pattern_id == 0
                            pattern_id = copy(n_patterns)
                        end
                    end
                end
            end
        end

        # update recurrent input
        b_last_spike = pre_spike_time .== t_prev # boolean indicate whether neurons spiked at last time step
        exc_input = weights[1:Ne, :]'*b_last_spike[1:Ne]
        inh_input = weights[(1+Ne):Ncells, :]'*b_last_spike[(1+Ne):Ncells]
        exc_input = exc_input + ext_input

        # update conductance
        g_inh = g_inh + (-g_inh + inh_input)/float(tau_gaba) * dt
        g_ampa = g_ampa + ((-g_ampa + exc_input) / float(tau_ampa)) * dt
        g_nmda = g_nmda + ((-g_nmda + g_ampa) / float(tau_nmda)) * dt
        g_exc = alpha * g_ampa + (1 - alpha) * g_nmda

        # update membrane potantial and firing thresholds
        
        # seperate inhibitory and excitatory
        v[1:Ne] = v[1:Ne] + ((v_rest - v[1:Ne]) + g_exc[1:Ne] .* (v_exc - v[1:Ne]) + g_inh[1:Ne] .* (v_inh - v[1:Ne])) / float(tau_m_e) * dt
        v[(1+Ne):Ncells] = v[(1+Ne):Ncells] + ((v_rest - v[(1+Ne):Ncells]) + g_exc[(1+Ne):Ncells] .* (v_exc - v[(1+Ne):Ncells]) + g_inh[(1+Ne):Ncells] .* (v_inh - v[(1+Ne):Ncells])) / float(tau_m_i) * dt
             
        spike_info = (v .> v_thr) .& (t .> t_allow_spike)
        unspike_info = .~spike_info
        
        if b_save_output
            v_me[loop_counter] = v[1]
            v_mi[loop_counter] = v[1+Ne]
            spike_time_mat[loop_counter, :] = spike_info
            loop_counter = loop_counter + 1
        end
        
        spike_neuron_idx = find(spike_info .== true)
        spike_neuron_exc_idx = intersect(1:Ne, spike_neuron_idx)
        spike_neuron_inh_idx = intersect(Ne+1:Ncells, spike_neuron_idx)
        pre_spike_time[spike_neuron_idx] = t
        t_allow_spike[spike_neuron_idx] = t + t_ref
        v[spike_neuron_idx] = v_spike_rest # reset membrane potential

        if tt > 30000

	    if  mod(tt, Int(30/dt)) == 0
            	if b_metaplasticity_ee
		       # update A2_m
			for cc = 1:Ne
		           A2_m[cc] = A2_m[cc] * ((r_t_e_ss[cc]/float(tau_r))/float(rho_0_e[cc]))
		           if A2_m[cc] > A2_m_initial * 1.85
                                A2_m[cc] = A2_m_initial * 1.85
                           elseif A2_m[cc] < A2_m_initial * 0.15
                                A2_m[cc] = A2_m_initial * 0.15
                           end
			end
		end
            end

            # update postsynaptic activity estimator
            for cc = 1:Ne
                if cc in spike_neuron_exc_idx
                    r_t_e_ss[cc] = r_t_e_ss[cc] - r_t_e_ss[cc]/float(tau_r) * dt + 1
                else
                    r_t_e_ss[cc] = r_t_e_ss[cc] - r_t_e_ss[cc]/float(tau_r) * dt
                end
            end
            
            for cc = Ne+1:Ncells
                if cc in spike_neuron_inh_idx
                    r_t_i_ss[cc-Ne] = r_t_i_ss[cc-Ne] - r_t_i_ss[cc-Ne]/float(tau_r) * dt + 1
                else
                    r_t_i_ss[cc-Ne] = r_t_i_ss[cc-Ne] - r_t_i_ss[cc-Ne]/float(tau_r) * dt
                end
            end
        end

        # triplet STDP e-e
        for cc = 1:Ne
            if cc in spike_neuron_exc_idx
                r1[cc] = r1[cc] - r1[cc]/float(tau_p) * dt + 1
                o1[cc] = o1[cc] - o1[cc]/float(tau_m) * dt + 1
                r2[cc] = r2[cc] - r2[cc]/float(tau_x) * dt
                o2[cc] = o2[cc] - o2[cc]/float(tau_y) * dt
            else
                r1[cc] = r1[cc] - r1[cc]/float(tau_p) * dt
                o1[cc] = o1[cc] - o1[cc]/float(tau_m) * dt
                r2[cc] = r2[cc] - r2[cc]/float(tau_x) * dt
                o2[cc] = o2[cc] - o2[cc]/float(tau_y) * dt
            end
        end
        
        if b_learning_active_e_e
            
            # triplet
            for cc in spike_neuron_exc_idx
                
                # E-E
                for dd = 1:Ne # cc is presynaptic neuron
                    if connections[cc, dd] == 0. # check connection
                        continue
                    end
                    weights_ltd_temp = copy(weights[cc, dd])
                    weights[cc, dd] = weights[cc, dd] - o1[dd] * (A2_m[dd] + A3_m[dd] * r2[cc])
                    if weights[cc, dd] > w_max_ee
                        weights[cc, dd] = w_max_ee
                    elseif weights[cc, dd] < w_min_ee
                        weights[cc, dd] = w_min_ee
                    end
                    weights_ltd[cc, dd] = weights_ltd[cc, dd] + weights[cc, dd] - weights_ltd_temp
                end

                for dd = 1:Ne # cc is postsynaptic neuron
                    if connections[dd, cc] == 0. # check connection
                        continue
                    end
                    weights_ltp_temp = copy(weights[dd, cc])
                    weights[dd, cc] = weights[dd, cc] + r1[dd] * (A2_p[cc] + A3_p[cc] * o2[cc])
                    if weights[dd, cc] > w_max_ee
                        weights[dd, cc] = w_max_ee
                    elseif weights[dd, cc] < w_min_ee
                        weights[dd, cc] = w_min_ee
                    end
                    weights_ltp[dd, cc] = weights_ltp[dd, cc] + weights[dd, cc] - weights_ltp_temp
                end
            end

            # synaptic scaling
            if b_synaptic_scaling_e_e
                for cc = 1:Ne
                    for dd = 1:Ne
                    weights_ss_temp = copy(weights[dd, cc])
                        weights[dd, cc] = weights[dd, cc] + eta_ss_e * weights[dd, cc] * (1 - r_t_e_ss[cc]/float(tau_r * rho_0_e[cc])) * dt / tau_ss
                        if weights[dd, cc] > w_max_ee
                            weights[dd, cc] = w_max_ee
                        elseif weights[dd, cc] < w_min_ee
                            weights[dd, cc] = w_min_ee
                        end
                    weights_ss[dd, cc] = weights_ss[dd, cc] + weights[dd, cc] - weights_ss_temp
                    end
                end
            end

            # intrinsic plasticity
            if b_intrinsic_plasticity
                for cc = 1:Ne
                    v_thr[cc] = v_thr[cc] + eta_ip_e * (r_t_e_ss[cc]/float(tau_r) - rho_0_e[cc]) * dt / tau_ip
                end
                
                for cc = Ne+1:Ncells
                    v_thr[cc] = v_thr[cc] + eta_ip_i * (r_t_i_ss[cc-Ne]/float(tau_r) - rho_0_i[cc-Ne]) * dt / tau_ip
                end
            end
        end

        for cc in spike_neuron_idx
            r2[cc] += 1
            o2[cc] += 1
        end

        # iSTDP part, i-e plasticity
        trace_istdp = trace_istdp - trace_istdp * dt / tau_istdp
        trace_istdp[spike_neuron_idx] = trace_istdp[spike_neuron_idx] + 1 # increase 1 if the neuron is spiking
        
        if b_learning_active_istdp
            for cc in spike_neuron_idx
                if cc <= Ne
                    for dd = (Ne+1):Ncells # cc is excitatory, postsynaptic neuron, dd is inhibitory, presynaptic neuron
                        if connections[dd, cc] == 0.
                            continue
                        end
                        weights[dd, cc] += eta_istdp * trace_istdp[dd]
                        if weights[dd, cc] > w_max_ie # weights can only increase in this loop
                            weights[dd, cc] = w_max_ie
                        end
                    end
                else
                    for dd = 1:Ne # cc is inhibitory, presynaptic neuron, dd is excitatory, postsynaptic neuron
                        if connections[cc, dd] == 0.
                            continue
                        end
                        weights[cc, dd] += eta_istdp * (trace_istdp[dd] - 2 * rho_0_e[dd] * tau_istdp)
                        if weights[cc, dd] > w_max_ie
                            weights[cc, dd] = w_max_ie
                        elseif weights[cc, dd] < w_min_ie
                            weights[cc, dd] = w_min_ie
                        end
                    end
                end
            end
        end
    end
end
