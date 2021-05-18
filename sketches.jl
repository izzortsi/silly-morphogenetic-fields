function get_modvn_neighborhood(x, y, N, M) #get the four Von Neumann neighbors of a point in a lattice toroidally

    x_left = mod1(x-1, N)
    x_right = mod1(x+1, N)
    y_up = mod1(y+1, M)
    y_down = mod1(y-1, M)

    neighbors = [(x_left, y), (x_right, y), (x, y_up), (x, y_down)]

    return neighbors

end


# %%

function update_states_vn!(mfield::MorphogeneticField)
    
    I = CartesianIndices(mfield.CellsTypes)
    
    types = deepcopy(mfield.CellsTypes)
    #states = deepcopy(mfield.CellsStates)

    for i in I
    
        tup_i = Tuple(i)
        nbhood = get_modvn_neighborhood(tup_i..., mfield.Dimensions...)
        nbs_types = [types[nb...] for nb in nbhood]
        nbstypes_mode = mode(nbs_types)
        mfield.CellsTypes[tup_i...] = nbstypes_mode

    end

end
# %%

function update_types_vn!(mfield::MorphogeneticField)
    
    I = CartesianIndices(mfield.CellsTypes)
    
    types = deepcopy(mfield.CellsTypes)
    states = deepcopy(mfield.CellsStates)

    for i in I
        
        tup_i = Tuple(i)

        if states[tup_i...] == 0
            tup_i = Tuple(i)
            nbhood = get_modvn_neighborhood(tup_i..., mfield.Dimensions...)
            nbs_types = [types[nb...] for nb in nbhood]
            nbstypes_mode = mode(nbs_types)
            mfield.CellsTypes[tup_i...] = nbstypes_mode
            mfield.CellsStates[tup_i...] = 1
        else
            mfield.CellsStates[tup_i...] = 0
        end

    end

end
# %%

function update_types_t!(mfield::MorphogeneticField)
    
    I = CartesianIndices(mfield.CellsTypes)
    
    types = deepcopy(mfield.CellsTypes)
    states = deepcopy(mfield.CellsStates)

    xdim = mfield.Dimensions[1] 
    ydim = mfield.Dimensions[2]

    for i in I
        
        tup_i = Tuple(i)

        if states[tup_i...] == 0

            x_0, y_0 = tup_i
            nbhood = []

            for i in x_0-1:x_0+1
                for j in y_0-1:y_0+1
                    if i != j
                        push!(nbhood, (mod1(i, xdim), mod1(j, ydim)))
                    end
                end
            end

            nbs_types = [types[nb...] for nb in nbhood]
            nbstypes_mode = mode(nbs_types)
            mfield.CellsTypes[tup_i...] = nbstypes_mode
            mfield.CellsStates[tup_i...] = 1
        else
            mfield.CellsStates[tup_i...] = 0
        end

    end

end
# %%

function update_types!(mfield::MorphogeneticField)
    
    I = CartesianIndices(mfield.CellsTypes)
    
    types = deepcopy(mfield.CellsTypes)
    states = deepcopy(mfield.CellsStates)

    for i in I
        
        tup_i = Tuple(i)

        if states[tup_i...] == 0

            x_0, y_0 = tup_i(i)
            nbhood = []

            for i in x_0-1:x_0+1
                for j in y_0-1:y_0+1
                    if i != j && i*j != 0 && i <= mfield.Dimensions[1] && j <= mfield.Dimensions[2]
                        push!(nbhood, (i, j))
                    end
                end
            end

            nbs_types = [types[nb...] for nb in nbhood]
            nbstypes_mode = mode(nbs_types)
            mfield.CellsTypes[tup_i...] = nbstypes_mode
            mfield.CellsStates[tup_i...] = 1
        else
            mfield.CellsStates[tup_i...] = 0
        end

    end

end
# %%

function update_states_t!(mfield::MorphogeneticField)
    
    I = CartesianIndices(mfield.CellsTypes)
    
    types = deepcopy(mfield.CellsTypes)
    #states = deepcopy(mfield.CellsStates)
    
    xdim = mfield.Dimensions[1] 
    ydim = mfield.Dimensions[2]

    for i in I
        tup_i = Tuple(i)
        x_0, y_0 = tup_i
        nbhood = []


        for i in x_0-1:x_0+1
            for j in y_0-1:y_0+1
                if i != j
                    push!(nbhood, (mod1(i, xdim), mod1(j, ydim)))
                end
            end
        end
        
        nbs_types = [types[nb...] for nb in nbhood]
        nbstypes_mode = mode(nbs_types)
        mfield.CellsTypes[tup_i...] = nbstypes_mode

    end

end