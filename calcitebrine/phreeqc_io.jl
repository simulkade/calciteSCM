# A script to facilitate reading input concentrations, 
# creating phreeqc file, running it, reading the output, and plotting it

# Read the input file
using CSV, DataFrames
include("SCM_functions.jl")

input_data = CSV.read("sample_input.csv")
(N, C) = size(input_data)
phqc_solution = Array{Any}(N)
ph_res        = Array{Any}(N)

MW_calcite = 100.09 # g/mol
#modified_so and wolthers
capacitance = [1.3, 4.5] 
#heberling
#capacitance = [0.45, 999.99] 

# create the SCM function
#modified wolthers constants[-3.58, -2.8, -2.2, 12.85, -24.73, 10.15, 1.55, 0.35]
#modified Sø constants [-3.52108, -2.61239, -2.72276, 14.5456, -24.73, 9.26094, 1.37606, 0.205842, -9.8144, -8.8383]
#heberling constants [0.501, 0.56, 1.68, 1.48, -0.05, 0.04, -7.07, 0.4, 0.56, 1.68, 1.48]

chalk_eq=[-3.52108, -2.61239, -2.72276, 14.5456, -24.73, 9.26094, 1.37606, 0.205842, -9.8144, -8.8383]

surface_reactions=modified_so(log_k=chalk_eq)

# create the phreeqc solutions and solid phases
for i in 1:N
    phqc_solution[i] = """
    SOLUTION $i
        temp  $(input_data[:temp][i])
        pressure  $(input_data[:pressure][i])
        -units mol/kgw
        pH     $(input_data[:pH][i])   charge
        Ba     $(input_data[:Ba][i])
        Ca     $(input_data[:Ca][i])
        Mg     $(input_data[:Mg][i])
        Na     $(input_data[:Na][i])
        K      $(input_data[:K][i])
        Cl     $(input_data[:Cl][i])
        S      $(input_data[:S][i])
        C      $(input_data[:C][i])
    EQUILIBRIUM_PHASES $i
        CO2(g)   $(-input_data[:p_CO2][i])  1.0
        Calcite  0.0  $(input_data[:calcite][i]/MW_calcite)
    # END
    SURFACE $i
		-sites_units density
		-cd_music 
		-equilibrate $i
		Chalk_a      $(input_data[:site_dens][i])      $(input_data[:sp_area][i])    $(input_data[:calcite][i])
		Chalk_c      $(input_data[:site_dens][i])
		-capacitance $(capacitance[1]) $(capacitance[2]) 
		#-donnan # uncomment, and see the effect on Charge Balance...
    """
end

# creae a selected output block
sel_output = """
SELECTED_OUTPUT 1
    -file phqc_out.csv
    -reset false
    USER_PUNCH
			-headings dielec sigma_chalk psi_chalk sigma_clay psi_clay psi_chalk1 psi_chalk2
			10 PUNCH EPS_R
			30 PUNCH EDL("sigma", "Chalk")
			50 PUNCH EDL("psi", "Chalk")
            51 PUNCH EDL("sigma", "Clay")
			52 PUNCH EDL("psi", "Clay")
			100 PUNCH EDL("psi1", "Chalk")
            110 PUNCH EDL("psi2", "Chalk")
    END
END
"""
 
# write everything into a file
zeta_pot = zeros(N)
for i in 1:N
    ph_file=open("phreeqc_file", "w")
    write(ph_file, surface_reactions*phqc_solution[i]*sel_output)
    close(ph_file)
    run(`phreeqc phreeqc_file`)
    ph_res[i]=CSV.read("phqc_out.csv"; delim = "\t  ")
    zeta_pot[i] = ph_res[i][:psi_chalk2][end]
end

# add the calculated zeta potential to the table
input_data[:zeta] = zeta_pot

# write to a file
CSV.write("results.csv", input_data)