function initialize_network()
    Ne = 800 # number of excitatory neuron
    Ni = 200 # number of inhibitory neuron
    wee = 0.2 # e-e weight
    wei = 0.2 # e-i weight
    wie = 2.0 # i-e weight
    wii = 2.0 # i-i weight
    prob = 0.2 # connectivity probability
    return Ne, Ni, wee, wei, wie, wii, prob
end
