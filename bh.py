"""
basinhopping: The basinhopping global optimization algorithm
"""
from __future__ import division, print_function, absolute_import

import numpy as np
import math
from numpy import cos, sin
import scipy.optimize
import collections
from scipy._lib._util import check_random_state

__all__ = ['basinhopping']


class Storage(object):
    """
    Class used to store the lowest energy structure
    """
    def __init__(self, minres):
        self._add(minres)

    def _add(self, minres):
        self.minres = minres
        self.minres.x = np.copy(minres.x)

    def update(self, minres):
        if minres.fun < self.minres.fun:
            self._add(minres)
            return True
        else:
            return False

    def get_lowest(self):
        return self.minres


class BasinHoppingRunner(object):
    def __init__(self, x0, minimizer, step_taking, accept_tests, disp=False):
        self.x = np.copy(x0)
        self.minimizer = minimizer
        self.step_taking = step_taking
        self.accept_tests = accept_tests
        self.disp = disp

        self.nstep = 0

        # initialize return object
        self.res = scipy.optimize.OptimizeResult()
        self.res.minimization_failures = 0

        # do initial minimization
        minres = minimizer(self.x)
        if not minres.success:
            self.res.minimization_failures += 1
            if self.disp:
                print("warning: basinhopping: local minimization failure")
        self.x = np.copy(minres.x)
        self.energy = minres.fun
        if self.disp:
            print("basinhopping step %d: f %g" % (self.nstep, self.energy))

        # initialize storage class
        self.storage = Storage(minres)

        if hasattr(minres, "nfev"):
            self.res.nfev = minres.nfev
        if hasattr(minres, "njev"):
            self.res.njev = minres.njev
        if hasattr(minres, "nhev"):
            self.res.nhev = minres.nhev

    def _monte_carlo_step(self):
        
        x_after_step = np.copy(self.x)
        x_after_step = self.step_taking(x_after_step)

        # do a local minimization
        minres = self.minimizer(x_after_step)
        x_after_quench = minres.x
        energy_after_quench = minres.fun
        if not minres.success:
            self.res.minimization_failures += 1
            if self.disp:
                print("warning: basinhopping: local minimization failure")

        if hasattr(minres, "nfev"):
            self.res.nfev += minres.nfev
        if hasattr(minres, "njev"):
            self.res.njev += minres.njev
        if hasattr(minres, "nhev"):
            self.res.nhev += minres.nhev

        accept = True
        for test in self.accept_tests:
            testres = test(f_new=energy_after_quench, x_new=x_after_quench,
                           f_old=self.energy, x_old=self.x)
            if testres == 'force accept':
                accept = True
                break
            elif testres is None:
                raise ValueError("accept_tests must return True, False, or "
                                 "'force accept'")
            elif not testres:
                accept = False

        
        if hasattr(self.step_taking, "report"):
            self.step_taking.report(accept, f_new=energy_after_quench,
                                    x_new=x_after_quench, f_old=self.energy,
                                    x_old=self.x)

        return accept, minres

    def one_cycle(self):
        """Do one cycle of the basinhopping algorithm
        """
        self.nstep += 1
        new_global_min = False

        accept, minres = self._monte_carlo_step()

        if accept:
            self.energy = minres.fun
            self.x = np.copy(minres.x)
            new_global_min = self.storage.update(minres)

        # print some information
        if self.disp:
            self.print_report(minres.fun, accept)
            if new_global_min:
                print("found new global minimum on step %d with function"
                      " value %g" % (self.nstep, self.energy))

        # save some variables as BasinHoppingRunner attributes
        self.xtrial = minres.x
        self.energy_trial = minres.fun
        self.accept = accept

        return new_global_min

    def print_report(self, energy_trial, accept):
        """print a status update"""
        minres = self.storage.get_lowest()
        print("basinhopping step %d: f %g trial_f %g accepted %d "
              " lowest_f %g" % (self.nstep, self.energy, energy_trial,
                                accept, minres.fun))


class AdaptiveStepsize(object):
    
    def __init__(self, takestep, accept_rate=0.5, interval=50, factor=0.9,
                 verbose=True):
        self.takestep = takestep
        self.target_accept_rate = accept_rate
        self.interval = interval
        self.factor = factor
        self.verbose = verbose

        self.nstep = 0
        self.nstep_tot = 0
        self.naccept = 0

    def __call__(self, x):
        return self.take_step(x)

    def _adjust_step_size(self):
        old_stepsize = self.takestep.stepsize
        accept_rate = float(self.naccept) / self.nstep
        if accept_rate > self.target_accept_rate:
            # We're accepting too many steps.  This generally means we're
            # trapped in a basin.  Take bigger steps
            self.takestep.stepsize /= self.factor
        else:
            # We're not accepting enough steps.  Take smaller steps
            self.takestep.stepsize *= self.factor
        if self.verbose:
            print("adaptive stepsize: acceptance rate %f target %f new "
                  "stepsize %g old stepsize %g" % (accept_rate,
                  self.target_accept_rate, self.takestep.stepsize,
                  old_stepsize))

    def take_step(self, x):
        self.nstep += 1
        self.nstep_tot += 1
        if self.nstep % self.interval == 0:
            self._adjust_step_size()
        return self.takestep(x)

    def report(self, accept, **kwargs):
        "called by basinhopping to report the result of the step"
        if accept:
            self.naccept += 1


class RandomDisplacement(object):
    """
    Add a random displacement of maximum size `stepsize` to each coordinate
    Calling this updates `x` in-place.
    Parameters
    ----------
    stepsize : float, optional
        Maximum stepsize in any dimension
    random_state : None or `np.random.RandomState` instance, optional
        The random number generator that generates the displacements
    """
    def __init__(self, stepsize=0.5, random_state=None):
        self.stepsize = stepsize
        self.random_state = check_random_state(random_state)

    def __call__(self, x):
        x += self.random_state.uniform(-self.stepsize, self.stepsize,
                                       np.shape(x))
        return x


class MinimizerWrapper(object):
    """
    wrap a minimizer function as a minimizer class
    """
    def __init__(self, minimizer, func=None, **kwargs):
        self.minimizer = minimizer
        self.func = func
        self.kwargs = kwargs

    def __call__(self, x0):
        if self.func is None:
            return self.minimizer(x0, **self.kwargs)
        else:
            return self.minimizer(self.func, x0, **self.kwargs)


class Metropolis(object):
    """
    Metropolis acceptance criterion
    Parameters
    ----------
    T : float
        The "temperature" parameter for the accept or reject criterion.
    random_state : None or `np.random.RandomState` object
        Random number generator used for acceptance test
    """
    def __init__(self, T, random_state=None):
        # Avoid ZeroDivisionError since "MBH can be regarded as a special case
        # of the BH framework with the Metropolis criterion, where temperature
        # T = 0."  (Reject all steps that increase energy.)
        self.beta = 1.0 / T if T != 0 else float('inf')
        self.random_state = check_random_state(random_state)

    def accept_reject(self, energy_new, energy_old):
        """
        If new energy is lower than old, it will always be accepted.
        If new is higher than old, there is a chance it will be accepted,
        less likely for larger differences.
        """
        w = math.exp(min(0, -float(energy_new - energy_old) * self.beta))
        rand = self.random_state.rand()
        return w >= rand

    def __call__(self, **kwargs):
        """
        f_new and f_old are mandatory in kwargs
        """
        return bool(self.accept_reject(kwargs["f_new"],
                    kwargs["f_old"]))


def basinhopping(func, x0, niter=100, T=1.0, stepsize=0.5,
                 minimizer_kwargs=None, take_step=None, accept_test=None,
                 callback=None, interval=50, disp=False, niter_success=None,
                 seed=None):
   
    x0 = np.array(x0)
    rng = check_random_state(seed)

    if minimizer_kwargs is None:
        minimizer_kwargs = dict()
    wrapped_minimizer = MinimizerWrapper(scipy.optimize.minimize, func,
                                         **minimizer_kwargs)

    # set up step-taking algorithm
    if take_step is not None:
        if not isinstance(take_step, collections.Callable):
            raise TypeError("take_step must be callable")
        # if take_step.stepsize exists then use AdaptiveStepsize to control
        # take_step.stepsize
        if hasattr(take_step, "stepsize"):
            take_step_wrapped = AdaptiveStepsize(take_step, interval=interval,
                                                 verbose=disp)
        else:
            take_step_wrapped = take_step
    else:
        # use default
        displace = RandomDisplacement(stepsize=stepsize, random_state=rng)
        take_step_wrapped = AdaptiveStepsize(displace, interval=interval,
                                             verbose=disp)

    # set up accept tests
    if accept_test is not None:
        if not isinstance(accept_test, collections.Callable):
            raise TypeError("accept_test must be callable")
        accept_tests = [accept_test]
    else:
        accept_tests = []
    # use default
    metropolis = Metropolis(T, random_state=rng)
    accept_tests.append(metropolis)

    if niter_success is None:
        niter_success = niter + 2

    bh = BasinHoppingRunner(x0, wrapped_minimizer, take_step_wrapped,
                            accept_tests, disp=disp)

    # start main iteration loop
    count, i = 0, 0
    message = ["requested number of basinhopping iterations completed"
               " successfully"]
    for i in range(niter):
        new_global_min = bh.one_cycle()

        if isinstance(callback, collections.Callable):
            # should we pass a copy of x?
            val = callback(bh.xtrial, bh.energy_trial, bh.accept)
            if val is not None:
                if val:
                    message = ["callback function requested stop early by"
                               "returning True"]
                    break

        count += 1
        if new_global_min:
            count = 0
        elif count > niter_success:
            message = ["success condition satisfied"]
            break

    # prepare return object
    res = bh.res
    res.lowest_optimization_result = bh.storage.get_lowest()
    res.x = np.copy(res.lowest_optimization_result.x)
    res.fun = res.lowest_optimization_result.fun
    res.message = message
    res.nit = i + 1
    return res


def _test_func2d_nograd(x):
    f = (cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0]
         + 1.010876184442655)
    return f


def _test_func2d(x):
    f = (cos(14.5 * x[0] - 0.3) + (x[0] + 0.2) * x[0] + cos(14.5 * x[1] -
         0.3) + (x[1] + 0.2) * x[1] + x[0] * x[1] + 1.963879482144252)
    df = np.zeros(2)
    df[0] = -14.5 * sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2 + x[1]
    df[1] = -14.5 * sin(14.5 * x[1] - 0.3) + 2. * x[1] + 0.2 + x[0]
    return f, df


if __name__ == "__main__":
    print("\n\nminimize a 2d function without gradient")
    # minimum expected at ~[-0.195, -0.1]
    kwargs = {"method": "L-BFGS-B"}
    x0 = np.array([1.0, 1.])
    scipy.optimize.minimize(_test_func2d_nograd, x0, **kwargs)
    ret = basinhopping(_test_func2d_nograd, x0, minimizer_kwargs=kwargs,
                       niter=200, disp=False)
    print("minimum expected at  func([-0.195, -0.1]) = 0.0")
    print(ret)

    print("\n\ntry a harder 2d problem")
    kwargs = {"method": "L-BFGS-B", "jac": True}
    x0 = np.array([1.0, 1.0])
    ret = basinhopping(_test_func2d, x0, minimizer_kwargs=kwargs, niter=200,
                       disp=False)
    print("minimum expected at ~, func([-0.19415263, -0.19415263]) = 0")
    print(ret)