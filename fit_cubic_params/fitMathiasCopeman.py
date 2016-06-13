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

def objective(c, pf, pexp, Texp, full_output = False):
    if c is not None:
        pf.set_c123(np.array(c).tolist())
    pcorr = pexp.copy()
    for i,T in enumerate(Texp):
        try:
            pcorr[i] = pf.saturation_pressure(T)
        except:
            pcorr[i] = 1000*pexp[i]
    err = (pcorr - pexp)/pexp
    ssq = np.sum(np.power(err, 2))
    if full_output:
        return dict(ssq = ssq,
                    pcorr = pcorr)
    else:
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
    
    pop = toolbox.population(n=1500)
    pop, log = eaSimple(pop, 
                        toolbox, 
                        cxpb=0.5, # Crossover probability
                        mutpb=0.3, # Mutation probability
                        ngen=30, 
                        stats=stats, 
                        halloffame=hof, 
                        verbose=True)
    best = max(pop, key=attrgetter("fitness"))
    print(best, best.fitness)
    return list(best)

if __name__=='__main__':
    import CoolProp.CoolProp as CP 
    fluid = 'R125'
    Tc,Tt,acentric,pc = [CP.PropsSI(k,fluid) for k in ['Tcrit','T_triple','acentric','p_critical']]
    print(Tc, Tt, pc, acentric)
    pf = PureFluid([Tc],[pc],[acentric],8.3144598,True)
    Texp = np.linspace(Tt, Tc-5, 500)
    pexp = CP.PropsSI('P','T',Texp,'Q',0,fluid)

    c = minimize_deap(objective, [0.0]*3, kwargs = dict(pf=pf, pexp=pexp, Texp=Texp))
    #c = [0.68638,-0.42475,0.72531]

    pp = objective(None, pf, pexp, Texp, full_output = True)['pcorr']
    plt.plot(Texp, pexp)
    plt.plot(Texp, pp)
    plt.yscale('log')
    plt.show()