# PKGs
include("GrowthModels.jl")
using ModelingToolkit, DifferentialEquations, Main.GrowthModels, Plots
plotly()

##################
# generate differential equations for base model
base_eqs_dict, base_lam = Main.GrowthModels.base_model_eqs();

# compose the base model from proteome fractions
base_model = Main.GrowthModels.compose_model(base_eqs_dict)

prob = ODEProblem(base_model, Main.GrowthModels._format(Main.GrowthModels.base_model_ss_values), (0, 1e4), Main.GrowthModels._format(Main.GrowthModels.base_model_parameters); jac=true);
sol  = solve(prob, Rodas4());
plot(sol)

# simulate host + heterologous
# generate differential equations for base model extended with a single heterologous protein
ha_eqs, ha_lam = Main.GrowthModels.het_model_eqs(input_eqs_dict = base_eqs_dict, input_lam = base_lam);
ha_het_model = Main.GrowthModels.compose_model(ha_eqs)

# time span
tspan = (0.0, 1e4)
# exposed relevant parameters
ns = 0.5
kurep = exp10(-2.6575)
kbrep = exp10(-1.3335)
wmaxrep = 150
kappa_ini = 0.1
# oganize in a vector
het_p = [ns, wmaxrep, kbrep, kurep, dmrep, dprep, kappa_ini]

# assemble into single parameter vector
append!(p, het_p)