#New NonA marsh model based on Morris et al 2002 model
using DifferentialEquations
using Plots

function NonAMarsh(dM,M,p,t)
    ########################################################
    #params
    #intrinsic growth rate
    SL = p

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
    #Biomass based on Morris et al 2002
    B = (a*M[1]) + (b*M[1]*M[1]) + (c)

    ########################################################

    #BIOMASS EQN — — time varying biomass to use for depth eqn 
    #biomass varies between 1 and 0.5
    Bt = M[1] * ((0.25)*sin(2*π*t)+0.75)
    #biomass varies between 1 and 0
    #Bt = B * ((0.5)*sin(2*π*t)+0.5)
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
    dM[1] = SL - (q+(TE*k*Bt))*M[1]
    ########################################################

end

#parameters (SLR)
#SLR is Charleston,SC is 3.32 mm/yr or 0.0032
p = (0.003)
#ICs (initial depth)
M0 = [0.4]

#time span
tspan = (0.0,1000.0)
prob = ODEProblem(NonAMarsh,M0,tspan,p)
sol = solve(prob)

#plots
plot(sol, xlabel = "Time (yr)", ylabel = "Depth (m)",lw = 3, legend = true)


