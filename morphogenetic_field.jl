# %%
using Plots
using Random
using StatsBase
# %%

mutable struct MorphogeneticField{T <: Integer, S <: Real}
    Dimensions::Tuple{Int64, Int64}
    Types::Array{T}
    States::Array{S}
    CellsTypes::Array{T, 2}
    CellsStates::Array{S, 2}
end
# %%
tarr = Array{Float64}([])
# %%
tarr
typeof(tarr)
# %%

function MorphogeneticField(
    dimensions::Tuple{Int64, Int64},
    types::Array{T},
    states::Array{S};
    cells_types::Array{T} = Array{T}([]),
    cells_states::Array{S} = Array{S}([]),
) where {T <: Integer, S <: Real}

    if length(cells_types) == 0
        cells_types = [
            rand(types) for x in 1:dimensions[1],
            y in 1:dimensions[2]
        ]
    end

    if length(cells_states) == 0
        cells_states = [
            rand(states) for x in 1:dimensions[1],
            y in 1:dimensions[2]
        ]
    end
    MorphogeneticField(
        dimensions,
        types,
        states,
        cells_types,
        cells_states,
    )
end

# %%

function update_types!(mfield::MorphogeneticField)

    I = CartesianIndices(mfield.CellsTypes)

    types = copy(mfield.CellsTypes)
    #states = deepcopy(mfield.CellsStates)

    for i in I
        tup_i = Tuple(i)
        x_0, y_0 = tup_i
        nbhood = []

        for i in (x_0 - 1):(x_0 + 1)
            for j in (y_0 - 1):(y_0 + 1)
                if i != j &&
                   i * j != 0 &&
                   i <= mfield.Dimensions[1] &&
                   j <= mfield.Dimensions[2]
                    push!(nbhood, (i, j))
                end
            end
        end
        nbs_types = [types[nb...] for nb in nbhood]
        nbstypes_mode = mode(nbs_types)
        mfield.CellsTypes[tup_i...] = nbstypes_mode
    end
end

# %%

# %%

mfield = MorphogeneticField((200, 200), [0, 1, 2], [0, 1])

# %%

heatmap(mfield.CellsTypes)

# %%

anim = @animate for i in 1:20
    heatmap!(mfield.CellsTypes)
    update_types!(mfield)
end
# %%

gif(anim, "mfield.gif", fps = 2)


