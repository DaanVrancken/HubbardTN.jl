####################
# State properties #
####################

"""
    dim_state(ψ::InfiniteMPS)

Determine the bond dimensions in an infinite MPS.
"""
function dim_state(ψ::InfiniteMPS)
    dimension = Int64.(zeros(length(ψ)))
    for i in 1:length(ψ)
        dimension[i] = dim(space(ψ.AL[i],1))
    end
    return dimension
end

"""
    density_e(ψ::InfiniteMPS, symm::SymmetryConfig)

Compute the number of electrons per site in the unit cell.
"""
function density_e(ψ::InfiniteMPS, symm::SymmetryConfig)
    n = number_e(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)
    bands = Int(length(ψ)/symm.cell_width)

    Ne = zeros(bands,symm.cell_width)
    for i in 1:bands
        for j in 1:symm.cell_width
            Ne[i,j] = real(expectation_value(ψ, (i+(j-1)*bands) => n))
        end
    end

    check = (sum(Ne)/(length(ψ)) ≈ symm.filling[1]/symm.filling[2])
    println("Filling is conserved: $check")

    return Ne
end


"""
    density_spin(ψ::InfiniteMPS, symm::SymmetryConfig)

Compute the spin density per site in the unit cell.
"""
function density_spin(ψ::InfiniteMPS, symm::SymmetryConfig)
    n_up = number_up(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)
    n_down = number_down(symm.particle_symmetry, symm.spin_symmetry; filling=symm.filling)

    bands = Int(length(ψ)/symm.cell_width)

    Nup = zeros(bands,symm.cell_width);
    Ndown = zeros(bands,symm.cell_width);
    for i in 1:bands
        for j in 1:symm.cell_width
            Nup[i,j] = real(expectation_value(ψ, (i+(j-1)*bands) => n_up))
            Ndown[i,j] = real(expectation_value(ψ, (i+(j-1)*bands) => n_down))
        end
    end

    return Nup, Ndown
end

"""
    calc_ms(ψ::InfiniteMPS, symm::SymmetryConfig)

Compute the staggered magnetization in an InfiniteMPS.
"""
function calc_ms(ψ::InfiniteMPS, symm::SymmetryConfig)
    up, down = density_spin(ψ, symm)
    Mag = up - down
    println(Mag)
    if !all(x -> isapprox(abs(x),abs(Mag[1,1]),rtol=10^(-6)), vec(Mag))
        @warn "Staggerd magnetization varies across unit cell: returning value for first site only."
    end
    return abs(Mag[1,1])
end


##########
# Saving #
##########

"""
    save_computation(d::Dict{String, Any}, path::String, file_name::String)

Save the output dictionary of e.g. `compute_groundstate` as a `.jld2` file at the specified path.
"""
function save_computation(d::Dict{String, Any}, path::String, file_name::String)
    ispath(path) || mkdir(path)
    @save joinpath(path, file_name*".jld2") d
end

"""
    load_computation(path_to_file::String)

Load the output dictionary of e.g. `compute_groundstate` stored as a `.jld2` file.
"""
function load_computation(path_to_file::String)
    @load joinpath(path_to_file) d
    return d
end

"""
    save_state(ψ::InfiniteMPS, path::String, name::String)

Save the tensors of an `InfiniteMPS` object to disk as individual `.jld2` files.

# Arguments
- `ψ::InfiniteMPS`: The infinite matrix product state (MPS) whose tensors will be saved.
- `path::String`: The base directory where the state folder will be created.
- `name::String`: The name of the subdirectory under `path` where the tensors will be stored.

# Description
This function creates a subdirectory `joinpath(path, name)` and saves each tensor 
`ψ.AL[i]` as a `.jld2` file named `state<i>.jld2` inside it. Each tensor is converted 
to a `Dict` before saving for serialization compatibility. The function prints 
a message after each tensor is successfully saved. This approach may be useful for
storing large MPS objects.
"""
function save_state(ψ::InfiniteMPS, path::String, name::String)
    path = joinpath(path,name)
    mkdir(path)
    for i in 1:length(ψ)
        d = convert(Dict,ψ.AL[i])
        @save joinpath(path,"state$i.jld2") d
        println("State $i saved.")
    end
end

"""
    load_state(path::String) -> InfiniteMPS

Load an `InfiniteMPS` object from a directory of saved `.jld2` tensor files.

# Arguments
- `path::String`: Path to the directory containing the saved MPS tensor files 
  (e.g., `state1.jld2`, `state2.jld2`, ...).

# Description
This function reconstructs an `InfiniteMPS` object previously saved with [`save_state`](@ref).  
It reads all `.jld2` files in the specified directory, converts each stored `Dict`
back into a `TensorMap`, and combines them into a periodic array before wrapping 
the result in an `InfiniteMPS` object.

# Returns
- `InfiniteMPS`: The reconstructed infinite matrix product state.
"""
function load_state(path::String)
    entries = readdir(path)
    file_count = count(entry -> isfile(joinpath(path, entry)), entries)

    @load joinpath(path,"state1.jld2") d
    A = [convert(TensorMap, d)]
    for i in 2:file_count
        @load joinpath(path,"state$i.jld2") d
        push!(A, convert(TensorMap, d))
    end

    return InfiniteMPS(PeriodicArray(A))
end