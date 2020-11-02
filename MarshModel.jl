#New NonA marsh model based on Morris et al 2002 model
using DifferentialEquations
using Plots


function NonAMarsh(dM,M,p,t)

    #params
    #intrinsic growth rate
    SL,r = p

    #set the params (Alizada et al 2016)
    #Plant params, in g/cm^2 (max biomass of 0.1088 g/cm2)
    a = 0.1000
    b = -0.3718
    c = 0.1021

    #dep params (in cm/yr)
    q = 0.0018
    k = 2.5*(10^-5)

    #model EQNS
    #plant biomass EQN
    dM[1] = r*M[1]*(1.0 - (M[1]/((a*M[2]) + (b*M[2]*M[2]) + (c))))
    #marsh platform depth EQN w/ SLR
    dM[2] = SL - (q+(k*M[1]))*M[2]
end

M0 = [0.01, 0.3]
p = (0.001,0.01)
tspan = (0.0,10000.0)
prob = ODEProblem(NonAMarsh,M0,tspan,p)
sol = solve(prob)

p1 = plot(sol,vars=(0,1), label = "Biomass", xlabel = "Time (yr)", ylabel = "Biomass (g/cm^2)",lw = 3, legend = true)
p2 = plot(sol,vars=(0,2), label = "Depth", xlabel = "Time (yr)", ylabel = "Depth (cm)", lw = 3, legend = true)
p3 = plot(sol,vars=(1,2), title= "phase space", xlabel = "Biomass (g/cm^2)", ylabel = "Depth (cm)", legend = false, lw = 3)

plot(p1, p2, p3, layout = (3, 1))