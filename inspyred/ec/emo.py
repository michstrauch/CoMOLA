"""
    ======================================================
    :mod:`emo` -- Evolutionary multiobjective optimization
    ======================================================
    
    This module provides the framework for making multiobjective evolutionary 
    computations.
    
    .. Copyright 2012 Inspired Intelligence Initiative

    .. This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

    .. This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

    .. You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>.
       
    .. module:: emo
    .. moduleauthor:: Aaron Garrett <aaron.lee.garrett@gmail.com>
"""
from inspyred.ec import ec
import math


class Pareto(object):
    """Represents a Pareto multiobjective solution.
    
    A Pareto solution is a set of multiobjective values that can be 
    compared to other Pareto values using Pareto preference. This means 
    that a solution dominates, or is better than, another solution if it 
    is better than or equal to the other solution in all objectives and
    strictly better in at least one objective.
    
    Since some problems may mix maximization and minimization among
    different objectives, an optional `maximize` parameter may be
    passed upon construction of the Pareto object. This parameter
    may be a list of Booleans of the same length as the set of 
    objective values. If this parameter is used, then the `maximize`
    parameter of the evolutionary computation's ``evolve`` method 
    should be left as the default True value in order to avoid
    confusion. (Setting the `evolve`'s parameter to False would
    essentially invert all of the Booleans in the Pareto `maximize`
    list.) So, if all objectives are of the same type (either
    maximization or minimization), then it is best simply to use
    the `maximize` parameter of the `evolve` method and to leave
    the `maximize` parameter of the Pareto initialization set to
    its default True value. However, if the objectives are mixed
    maximization and minimization, it is best to leave the ``evolve``'s
    `maximize` parameter set to its default True value and specify
    the Pareto's `maximize` list to the appropriate Booleans.
    
    The typical usage is as follows::
    
        @inspyred.ec.evaluators.evaluator
        def my_evaluator(candidate, args):
            obj1 = 1 # Calculate objective 1
            obj2 = 2 # Calculate objective 2
            obj3 = 3 # Calculate objective 3
            return emo.Pareto([obj1, obj2, obj3])
    
    """
    def __init__(self, values=None, maximize=True):
        if values is None:
            values = []
        self.values = values
        try:
            iter(maximize)
        except TypeError:
            maximize = [maximize for v in values]
        self.maximize = maximize
        
    def __len__(self):
        return len(self.values)
    
    def __getitem__(self, key):
        return self.values[key]
        
    def __iter__(self):
        return iter(self.values)
    
    def __lt__(self, other):
        if len(self.values) != len(other.values):
            raise NotImplementedError
        else:
            not_worse = True
            strictly_better = False
            for x, y, m in zip(self.values, other.values, self.maximize):
                if m:
                    if x > y:
                        not_worse = False
                    elif y > x:
                        strictly_better = True
                else:
                    if x < y:
                        not_worse = False
                    elif y < x:
                        strictly_better = True
            return not_worse and strictly_better
            
    def __le__(self, other):
        return self < other or not other < self
        
    def __gt__(self, other):
        return other < self
        
    def __ge__(self, other):
        return other < self or not self < other
        
    def __eq__(self, other):
        return self.values == other.values
        
    def __ne__(self, other):
        return not (self == other)
    
    def __str__(self):
        return str(self.values)
        
    def __repr__(self):
        return str(self.values)


class NSGA2(ec.EvolutionaryComputation):
    """Evolutionary computation representing the nondominated sorting genetic algorithm.
    
    This class represents the nondominated sorting genetic algorithm (NSGA-II)
    of Kalyanmoy Deb et al. It uses nondominated sorting with crowding for 
    replacement, binary tournament selection to produce *population size*
    children, and a Pareto archival strategy. The remaining operators take 
    on the typical default values but they may be specified by the designer.
    
    """
    def __init__(self, random):
        ec.EvolutionaryComputation.__init__(self, random)
        self.archiver = ec.archivers.best_archiver
        self.replacer = ec.replacers.nsga_replacement
        self.selector = ec.selectors.tournament_selection
    
    def evolve(self, generator, evaluator, pop_size=100, seeds=None, maximize=True, bounder=None, is_available = None, custom_individual= [],archive = None, num_generation= 0 , **args):
        args.setdefault('num_selected', pop_size)
        args.setdefault('tournament_size', 2)
        previous_arc = archive
        return ec.EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, maximize, bounder, is_available, custom_individual, previous_arc, num_generation, **args)

    
class PAES(ec.EvolutionaryComputation):
    """Evolutionary computation representing the Pareto Archived Evolution Strategy.
    
    This class represents the Pareto Archived Evolution Strategy of Joshua
    Knowles and David Corne. It is essentially a (1+1)-ES with an adaptive
    grid archive that is used as a part of the replacement process. 
    
    """
    def __init__(self, random):
        ec.EvolutionaryComputation.__init__(self, random)
        self.archiver = ec.archivers.adaptive_grid_archiver
        self.selector = ec.selectors.default_selection
        self.variator = ec.variators.gaussian_mutation
        self.replacer = ec.replacers.paes_replacement  

    def evolve(self, generator, evaluator, pop_size=1, seeds=None, maximize=True, bounder=None, **args):
        final_pop = ec.EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, maximize, bounder, **args)
        try:
            del self.archiver.grid_population
        except AttributeError:
            pass
        try:
            del self.archiver.global_smallest
        except AttributeError:
            pass
        try:
            del self.archiver.global_largest
        except AttributeError:
            pass
        return final_pop
    

