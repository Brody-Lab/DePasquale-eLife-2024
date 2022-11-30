"""
"""
function CoM(x, t, thresh, margin)
    
    isevent = false
    go = true
    tprime = 5
    y = []
    
    while go
                
        if ((t - tprime) >= 1) & ((t + tprime) <= length(x))
        
            y = x[t - tprime : t + tprime]

            if (x[t - 1] < thresh) & (x[t + 1] > thresh)
                
                if all(x[t - tprime : t - 1] .< thresh) & all(x[t + 1 : t + tprime] .> thresh)
                    if any(x[t - tprime : t - 1] .< (thresh - margin)) & any(x[t + 1 : t + tprime] .> (thresh + margin))
                        isevent = true
                        go = false
                    else
                        tprime = tprime + 1
                    end
                else
                    isevent = false
                    go = false
                end
                
            elseif (x[t - 1] > thresh) & (x[t + 1] < thresh)
                
                if all(x[t - tprime : t - 1] .> thresh) & all(x[t + 1 : t + tprime] .< thresh)
                    if any(x[t - tprime : t - 1] .> (thresh + margin)) & any(x[t + 1 : t + tprime] .< (thresh - margin))
                        isevent = true
                        go = false
                    else
                        tprime = tprime + 1
                    end
                else
                    isevent = false
                    go = false
                end
            else
                tprime = tprime + 1
            end
            
        else
            
            isevent = false
            go = false
            
        end
        
    end
    
    return y, tprime, go, isevent
end