function simulation_new()
    println("Initializing Network")
    
    Ne, Ni, wee, wei, wie, wii, prob = initialize_network()
    Ncells = Ne+Ni
    
    # load initial weight matrix
    initial_weight_file = matopen("initial_weight_matrix.mat")
    connections = read(initial_weight_file, "weight_matrix")
    close(initial_weight_file)

    # load popmembers
    popmembers_file = matopen("popmembers.mat")
    popmembers = read(popmembers_file, "popluation_members")
    close(popmembers_file)

    # load pretrained weight matrix
    weight_file = matopen("weight_matrix_8000000.mat")
    weights = read(weight_file, "weight_matrix")
    close(weight_file)
    
    simulation(weights, connections, popmembers)

end
