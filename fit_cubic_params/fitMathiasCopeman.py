from __future__ import print_function

import os, sys, subprocess
import numpy as np, matplotlib.pyplot as plt, scipy.optimize

try:
    from PureFluid import PureFluid
except ImportError as IE:
    print('building PureFluid module, please be patient...')
    if not os.path.exists('build'): os.mkdir('build')
    subprocess.call('cmake .. -G "Visual Studio 14 2015 Win64"', cwd='build', stdout = sys.stdout)
    subprocess.call('cmake --build . --config Release', cwd='build', stdout = sys.stdout)
    # Try again to import, should work
    from PureFluid import PureFluid

import random
from deap import algorithms, base, creator, tools
from deap.algorithms import eaSimple
from operator import attrgetter

def objective(c, pf, pexp, Texp, pcorr, rhoLcorr):
    if c is not None:
        pf.set_c123(np.array(c).tolist())
    for i,T in enumerate(Texp):
        try:
            pcorr[i] = pf.saturation_pressure(T)
            rhoLcorr[i] = pf.rhomolar();
        except BaseException as BE:
            pcorr[i] = 1000*pexp[i];
            rhoLcorr[i] = np.nan;
    err = (pcorr - pexp)/pexp
    ssq = np.sum(np.power(err, 2))
    return (ssq,)

def minimize_deap(f, p0, kwargs):
    # weight is -1 because we want to minimize the error
    creator.create("FitnessMax", base.Fitness, weights=(-1.0,)) 
    creator.create("Individual", list, fitness=creator.FitnessMax)
    
    # See: http://deap.readthedocs.org/en/master/tutorials/basic/part1.html#a-funky-one    
    toolbox = base.Toolbox()

    toolbox.register("attr", random.uniform, -2, 2)
    toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.attr,), n = len(p0))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", f, **kwargs)
    if 'DeltaPenality' in dir(tools):
        toolbox.decorate("evaluate", tools.DeltaPenality(feasible, 100000))
    # If two individuals mate, interpolate between them, allow for a bit of extrapolation
    toolbox.register("mate", tools.cxBlend, alpha = 0.3)
    sigma = np.array([0.01]*len(p0)).tolist()

    toolbox.register("mutate", tools.mutGaussian, mu = 0, sigma = sigma, indpb=1.0)
    toolbox.register("select", tools.selTournament, tournsize=5) # The greater the tournament size, the greater the selection pressure 
    
    hof = tools.HallOfFame(50)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("min", np.min)
    
    pop = toolbox.population(n=500)
    pop, log = eaSimple(pop, 
                        toolbox, 
                        cxpb=0.5, # Crossover probability
                        mutpb=0.3, # Mutation probability
                        ngen=20, 
                        stats=stats, 
                        halloffame=hof, 
                        verbose=True)
    best = max(pop, key=attrgetter("fitness"))
    print(best, best.fitness)
    return list(best)


def do_fluid(fluid):
    import CoolProp.CoolProp as CP
    Tc,Tt,acentric,pc = [CP.PropsSI(k,fluid) for k in ['Tcrit','T_triple','acentric','p_critical']]
    pf = PureFluid([Tc],[pc],[acentric],8.3144598,True)
    Texp = np.linspace(Tt, Tc-5, 200)
    pexp = CP.PropsSI('P','T',Texp,'Q',0,fluid)
    rhoLexp = CP.PropsSI('Dmolar','T',Texp,'Q',0,fluid)
    pcorr = pexp.copy()
    rhoLcorr = pexp.copy()

    # Actually do the fit and obtain the attractive parameters for EOS
    a = minimize_deap(objective, [0.0]*3, kwargs = dict(pf=pf, 
                                                        pexp=pexp, 
                                                        pcorr=pcorr, 
                                                        rhoLcorr=rhoLcorr, 
                                                        Texp=Texp))

    # Calculate the values based on the obtained coefficients
    objective(a, pf, pexp, Texp, pcorr, rhoLcorr)

    # Pressure deviation plot
    plt.plot(Texp, (pexp/pcorr-1)*100)
    plt.xlabel('T')
    plt.ylabel(r'$(p_{\rm exp}/p_{\rm pred}-1)\times 100$ (%)')
    plt.savefig(fluid+'_pdev.png')
    plt.close()

    # And now we solve for the volume translation term c by evaluating the difference between 
    # the actual saturated liquid density (from EOS) and the saturated liquid density 
    # predicted by cubic
    pf.saturation_pressure(0.7*Tc)
    c = 1/CP.PropsSI('Dmolar','T',0.7*Tc,'Q',0,fluid) - 1/pf.rhomolar();
    rhoLvt = 1/(1/rhoLcorr+c); # corrected liquid density after volume translation

    # Volume deviation
    plt.plot(Texp, (rhoLexp/rhoLcorr-1)*100, label = 'predicted')
    plt.plot(Texp, (rhoLexp/rhoLvt-1)*100, label = 'volume-translated')
    plt.axvline(0.7*Tc, dashes=[2,2])
    plt.xlabel('T')
    plt.ylabel(r'$(\rho_{\rm exp}/\rho_{\rm pred}-1)\times 100$ (%)')
    plt.legend(loc='best')
    plt.savefig(fluid+'_rhodev.png')
    plt.close()

if __name__=='__main__':
    do_fluid('ARGON')
    do_fluid('R125')
    do_fluid('WATER')