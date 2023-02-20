# PKGs
include("GrowthModels.jl")
using ModelingToolkit, DifferentialEquations, Main.GrowthModels, Plots
plotly()

##################

# generate differential equations for base model
base_eqs_dict, base_lam, base_params, base_ss = Main.GrowthModels.base_model_eqs();

# compose the base model from proteome fractions
base_model = Main.GrowthModels.compose_model(base_eqs_dict)

prob = ODEProblem(base_model, Main.GrowthModels._format(base_ss), (0, 1e4), Main.GrowthModels._format(base_params); jac=true);
sol  = solve(prob, Rodas4());
plot(sol)

# simulate host + heterologous
# generate differential equations for base model extended with a single heterologous protein
ha_eqs, ha_lam, ha_params, ha_ss = Main.GrowthModels.het_model_eqs(input_eqs_dict   = base_eqs_dict, 
                                                                   input_lam        = base_lam, 
                                                                   input_param_vals = base_params,
                                                                   input_ss_vals    = base_ss);
ha_het_model = Main.GrowthModels.compose_model(ha_eqs)

prob_h = ODEProblem(ha_het_model, Main.GrowthModels._format(ha_ss), (0, 1e4), Main.GrowthModels._format(ha_params); jac=true);
sol_h  = solve(prob_h, Rodas4());
plot(sol_h)


# host aware vectors need to get the steady state values of base model
function map_vals(model::ODESystem, sol::ODESolution)
    return [states(model)[idx] => sol[:, end][idx] for idx in range(1, length(sol[:, end]))];
end;

function update_multiple!(vect::Vector, dict::Dict)
    [Main.GrowthModels._update!(dict, key, val) for (elm, val) in vect for key in keys(dict) if isequal(elm, dict[key][1])];
end;


update_multiple!(map_vals(base_model, sol), base_ss)