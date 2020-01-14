include("BaselineModel.jl")
include("EstimationBaseline.jl")
include("SetupBaseline.jl")

# TODO:  next, update choice probability functions to allow for aging of child, and for ageing out
# TODO: update the meaning of NK to make sure it includes 0
# TODO: fix the fact that the budget function (1) does not distinguish between arms with same budget and (2) May have fewer years than we wanted
# option 1 here: go back and remove unnecessary years from the moment file
# option 2: code to just extrapolate for years outside of budget (this works for future stuff)
# TODO: check the moment simulation code and make sure it works as intended
# TODO: check the criterion function and make sure it works
