# An script to facilitate reading input concentrations, 
# creating phreeqc file, running it, reding the output, and plotting it

# load the required packages
using CSV, DataFrames, LsqFit 
# using LsqFit, Optim
include("SCM_functions.jl")

# read the data file
input_data = CSV.read("sample_input.csv")
capacitance = [1.3, 4.5]
zeta_exp = input_data[:zeta]

# calculate zeta potential for all data

function zeta_calc(chalk_eq, input_data, capacitance)
    (N, C) = size(input_data)
    MW_calcite = 100.09 # g/mol
    phqc_solution = Array{Any}(N)
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
            C      $(input_data[:C][i]) as HCO3
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
    # create the surface reactions
    surface_reactions=modified_so(log_k=chalk_eq)
    # write everything into a file
    zeta_pot = zeros(N)
    for i in 1:N
        ph_file=open("phreeqc_file", "w")
        write(ph_file, surface_reactions*phqc_solution[i]*sel_output)
        close(ph_file)
        run(`phreeqc phreeqc_file`)
        ph_res = CSV.read("phqc_out.csv"; delim = "\t  ")
        zeta_pot[i] = ph_res[:psi_chalk2][end]
    end
    return zeta_pot
end # end of zeta_calc function

complex_chalk_brine = param -> zeta_calc(param, input_data, capacitance)
model = (x, param) -> zeta_calc(param, input_data, capacitance)

function objective_function(param::Vector{Float64}, grad::Vector{Float64})
    zeta_calc = complex_chalk_brine(param)
    zeta_error=(zeta_calc-zeta_exp)
    if length(grad)>0
      dx=1e-5
      for i in 1:length(param)
        x=copy(param)
        x[i]+=dx
        grad_zeta=(complex_chalk_brine(x)-zeta_calc)/dx
        grad[i]=2.0*sum(grad_zeta.*zeta_error)
      end
    end
    return sum(zeta_error.*zeta_error)
end

function myconstraint(results::Vector{Float64}, param::Vector{Float64}, grad::Array{Float64,2})
    # zeta_exp.*zeta_calc > 0.0
    # -zeta_exp*zeta_calc<0
    # grad nxm [j,i] dci/dxj
    zeta_calc = complex_chalk_brine(param)
    zeta_mult = zeta_calc.*zeta_exp
    m=length(zeta_exp)
    n=length(param)
    if length(grad)>0
      dx=1e-5
      for i in 1:n
        x=copy(param)
        x[i]+=dx
        grad_zeta=(complex_chalk_brine(x)-zeta_calc)/dx
        grad[i,:]=-zeta_exp.*grad_zeta
      end
    end
    results[:]=-zeta_mult
end
  
  # LN_COBYLA works fine
  # LD_SLSQP
chalk_logk_init =[-3.58, -2.8, -2.2, 12.85, -24.73, 10.15, 1.55, 0.35, -9.8144, -7.5]

chalk_logk_lb=chalk_logk_init-0.5*abs(chalk_logk_init)
chalk_logk_ub=chalk_logk_init+0.5*abs(chalk_logk_init)
N_data = length(zeta_exp)
w = ones(N_data)
x_dummy = collect(linspace(1,length(zeta_exp)))
fit = curve_fit(model, x_dummy, zeta_exp, w, chalk_logk_init, lower = chalk_logk_lb, upper = chalk_logk_ub)

println(fit.param)