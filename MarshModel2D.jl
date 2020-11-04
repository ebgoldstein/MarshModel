#New NonA marsh model based on Morris et al 2002 model
using DifferentialEquations
using Plots

function NonAMarsh(dM,M,p,t)
    ########################################################
    #params
    #intrinsic growth rate
    SL,r = p

    #set the params (from Morris et al 2002)
    #max biomass is ~1800 g/m^2
    #Plant params, in g/m^2
    a = 15500
    b = -18550
    c = -1364

    #dep params (in m/yr)
    q = 0.0018
    k = 1.5*(10^-5)

    ########################################################
    #carrying capacity based on Morris et al 2002
    κ = (a*M[2]) + (b*M[2]*M[2]) + (c)

    #plant biomass EQN
    dM[1] = r*M[1]*(1.0 - (M[1]/κ))

    ########################################################

    #BIOMASS EQN — — time varying biomass to use for depth eqn 
    #biomass varies between 1 and 0.5
    #Bt = M[1] * ((0.25)*sin(2*π*t)+0.75)
    #biomass varies between 1 and 0
    Bt = M[1] * ((0.5)*sin(2*π*t)+0.5)
    #biomass is constant
    #Bt = M[1]

    ########################################################

    #trapping EQN
    # at 0 biomass, 0 efficiency
    # at 1000 biomass, 100% efficiency
    # at 2000 biomass, 0 efficiency

    pp1 = 500
    pp2 = 1000000

    TE = (Bt/pp1) - ((Bt*Bt)/pp2)
    #no blocking trapping term
    #TE = 1

    #marsh platform depth EQN w/ SLR
    dM[2] = SL - (q+(TE*k*Bt))*M[2]
    ########################################################

end

#parameters (SLR, r)
#SLR is Charleston,SC is 3.32 mm/yr or 0.0032
p = (0.003, 0.1)
#ICs (biomass, depth)
M0 = [100, 0.6]

#time span
tspan = (0.0,1000.0)
prob = ODEProblem(NonAMarsh,M0,tspan,p)
sol = solve(prob)

#plots
p1 = plot(sol,vars=(0,1), label = "Biomass", xlabel = "Time (yr)", ylabel = "Biomass (g/m^2)",lw = 3, legend = true)
p2 = plot(sol,vars=(0,2), label = "Depth", xlabel = "Time (yr)", ylabel = "Depth (m)", lw = 3, legend = true)
p3 = plot(sol,vars=(1,2), title= "Phase space", xlabel = "Biomass (g/m^2)", ylabel = "Depth (m)", legend = false, lw = 3)
plot(p1, p2, p3, layout = (3, 1))