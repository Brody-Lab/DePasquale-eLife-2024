using DSP, Missings, Statistics, Distributions

"""
"""
pad_xcorr(r,nT,max,cells) = map((r,nT)-> vcat(missings(max - nT), 
        xcorr(r[cells[1]], r[cells[2]]; padmode=:none) ./ 
        (length(r[cells[1]]) .- vcat(length(r[cells[1]])-1:-1:0, 1:length(r[cells[1]])-1)), 
        missings(max-nT)), r, nT)


"""
"""
function do_xcorrs_real(data, sess, μ_λ_PSTH)
    
    dt = data[sess][1].input_data.dt
    ncells = data[sess][1].ncells
    cells = vcat([repeat([i], 2) for i in 1:ncells], collect(combinations(1:ncells,2)))
    
    nT = map(x-> length(x[1]), μ_λ_PSTH[sess])
    mymax = maximum(nT);
    spikes_PSTH = μ_λ_PSTH[sess]  
    
    mspikes_PSTH = mean.(map(i-> mean.(getindex.(spikes_PSTH, i)), 1:ncells))
    
    xcorr_spikes_PSTH = map(cells-> 
        mapslices(x-> mean(skipmissing(x))/mspikes_PSTH[cells[1]] .- mspikes_PSTH[cells[2]], 
            hcat(pad_xcorr(spikes_PSTH, nT, mymax, cells)...); dims=2), 
        cells); 

    nT = map(x-> x.input_data.binned_clicks.nT, data[sess]) .+ 2 * data[sess][1].input_data.pad;
    mymax = maximum(nT);
    spikes = map(x-> x.spikes/dt, data[sess]);
    
    mspikes = mean.(map(i-> mean.(getindex.(spikes, i)), 1:ncells))
    
    xcorr_spikes = map(cells-> 
        mapslices(x-> mean(skipmissing(x))/mspikes[cells[1]] .- mspikes[cells[2]], 
            hcat(pad_xcorr(spikes, nT, mymax, cells)...); dims=2), 
        cells);
        
    result = xcorr_spikes - xcorr_spikes_PSTH;

    return mymax, result, cells, spikes, nT
    
end


"""
"""
function do_xcorrs_synth(data, sess, μ_λ_syn, μ_λ_syn_PSTH)
    
    dt = data[sess][1].input_data.dt
    ncells = data[sess][1].ncells
    cells = vcat([repeat([i], 2) for i in 1:ncells], collect(combinations(1:ncells,2)))

    nT = map(x-> length(x[1]), μ_λ_syn_PSTH[sess])
    mymax = maximum(nT);
    spikes_syn_PSTH = μ_λ_syn_PSTH[sess]
    
    mspikes_syn_PSTH = mean.(map(i-> mean.(getindex.(spikes_syn_PSTH, i)), 1:ncells))
    
    xcorr_spikes_syn_PSTH = map(cells-> 
        mapslices(x-> mean(skipmissing(x))/mspikes_syn_PSTH[cells[1]] .- mspikes_syn_PSTH[cells[2]], 
            hcat(pad_xcorr(spikes_syn_PSTH, nT, mymax, cells)...); dims=2), 
        cells);  
    
    nT = map(x-> length(x[1]), μ_λ_syn[sess])
    mymax = maximum(nT);
    spikes_syn = μ_λ_syn[sess]  
    
    mspikes_syn = mean.(map(i-> mean.(getindex.(spikes_syn, i)), 1:ncells))
    
    xcorr_spikes_syn = map(cells-> 
        mapslices(x-> mean(skipmissing(x))/mspikes_syn[cells[1]] .- mspikes_syn[cells[2]], 
            hcat(pad_xcorr(spikes_syn, nT, mymax, cells)...); dims=2), 
        cells); 
        
    result = xcorr_spikes_syn - xcorr_spikes_syn_PSTH;
    
    return mymax, result, cells, spikes_syn, nT
    
end