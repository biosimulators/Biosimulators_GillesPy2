""" Data model for GillesPy2 algorithms and their parameters

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-12-09
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

import enum
import gillespy2

__all__ = [
    'Algorithm', 'AlgorithmParameter', 'VodeMethod', 'HybridTauIntegrationMethod',
    'kisao_algorithm_map',
]


class Algorithm(object):
    """ Simulation algorithm

    Attributes:
        name (:obj:`str`): name
        solver (:obj:`type`): solver
        solver_args (:obj:`dict`): solver arguments
        parameters (:obj:`dict`): dictionary that maps KiSAO ids to :obj:`AlgorithmParameter`\\ s
    """

    def __init__(self, name, solver, parameters=None, **solver_args):
        """
        Args:
            name (:obj:`str`): name
            solver (:obj:`type`): solver
            **solver_args (:obj:`dict`): solver arguments
            parameters (:obj:`dict`, optional): dictionary that maps KiSAO ids to :obj:`AlgorithmParameter`\\ s
        """
        self.name = name
        self.solver = solver
        self.solver_args = solver_args
        self.parameters = parameters or {}


class AlgorithmParameter(object):
    """ Simulation algorithm parameter

    Attributes:
        name (:obj:`str`): name
        key (:obj:`str`): key
        data_type (:obj:`type`): data type
        default (:obj:`object`): defualt value
    """

    def __init__(self, name, key, data_type, default):
        """
        Args:
            name (:obj:`str`): name
            key (:obj:`str`): key
            data_type (:obj:`type`): data type
            default (:obj:`float`): defualt value
        """
        self.name = name
        self.key = key
        self.data_type = data_type
        self.default = default

    def set_value(self, solver_args, str_value):
        """ Apply the value of a parameter to a data structure of solver arguments

        Args:
            solver_args (:obj:`dict`): solver arguments
            str_value (:obj:`string`): string representation of parameter value

        Raises:
            :obj:`ValueError`: if :obj:`str_value` is not a valid string representation
                of the data type of the parameter
            :obj:`NotImplementedError`: if :obj:`str_value` is not a valid value of an
                enumerated parameter
        """
        keys = self.key.split('.')
        for key in keys[0:-1]:
            if key in solver_args:
                nested_solver_args = solver_args[key]
            else:
                nested_solver_args = {}
                solver_args[key] = nested_solver_args
            solver_args = nested_solver_args

        if not str_value:
            value = None

        elif self.data_type == bool:
            if str_value.lower() == 'false' or str_value == '0':
                value = False
            elif str_value.lower() == 'true' or str_value == '1':
                value = True
            else:
                raise ValueError("Value '{}' is not a valid Boolean".format(str_value))

        elif self.data_type == int:
            try:
                value = int(str_value)
            except ValueError:
                raise ValueError("Value '{}' is not a valid integer".format(str_value))

        elif self.data_type == float:
            try:
                value = float(str_value)
            except ValueError:
                raise ValueError("Value '{}' is not a valid float".format(str_value))

        elif issubclass(self.data_type, enum.Enum):
            try:
                value = self.data_type(str_value).name
            except ValueError:
                raise NotImplementedError(
                    '{} is not a supported value of {}. The value of {} must be one of the following:\n  - {}'.format(
                        str_value, self.name, self.name,
                        '\n  - '.join('{}: {}'.format(value, name) for name, value in self.data_type.__members__.items())))

        else:
            raise NotImplementedError('Data type {} is not supported'.format(self.data_type.__name__))

        solver_args[keys[-1]] = value


class VodeMethod(str, enum.Enum):
    """ VODE method """
    bdf = 'KISAO_0000288'  # stiff
    adams = 'KISAO_0000289'  # non-stiff


class HybridTauIntegrationMethod(str, enum.Enum):
    """ Hybrid tau integration method """
    BDF = 'KISAO_0000288'
    LSODA = 'KISAO_0000088'
    Radau = 'KISAO_0000304'
    RK23 = 'KISAO_0000537'
    RK45 = 'KISAO_0000032'


kisao_algorithm_map = {
    'KISAO_0000088': Algorithm("LSODA", gillespy2.ODESolver, integrator="lsoda", parameters={
        'KISAO_0000211': AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12),
        'KISAO_0000209': AlgorithmParameter("relative tolerance", 'integrator_options.rtol', float, 1e-6),
        'KISAO_0000480': AlgorithmParameter("lower half bandwith", 'integrator_options.lband', int, None),
        'KISAO_0000479': AlgorithmParameter("upper half bandwith", 'integrator_options.uband', int, None),
        'KISAO_0000415': AlgorithmParameter("maximum number of steps", 'integrator_options.nsteps', int, 500),
        'KISAO_0000559': AlgorithmParameter("initial step size", 'integrator_options.first_step', float, 0.0),
        'KISAO_0000485': AlgorithmParameter("minimum step size", 'integrator_options.min_step', float, 0.0),
        'KISAO_0000467': AlgorithmParameter("maximum step size", 'integrator_options.max_step', float, float("inf")),
        'KISAO_0000219': AlgorithmParameter("maximum non-stiff order (Adams order)", 'integrator_options.max_order_ns', int, 12),
        'KISAO_0000220': AlgorithmParameter("maximum stiff order (BDF order)", 'integrator_options.max_order_s', int, 5),
    }),

    'KISAO_0000087': Algorithm("dopri5", gillespy2.ODESolver, integrator="dopri5", parameters={
        'KISAO_0000211': AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12),
        'KISAO_0000209': AlgorithmParameter("relative tolerance", 'integrator_options.rtol', float, 1e-6),
        'KISAO_0000415': AlgorithmParameter("maximum number of steps", 'integrator_options.nsteps', int, 500),
        'KISAO_0000559': AlgorithmParameter("initial step size", 'integrator_options.first_step', float, 0.0),
        'KISAO_0000467': AlgorithmParameter("maximum step size", 'integrator_options.max_step', float, float("inf")),
        'KISAO_0000538': AlgorithmParameter("safety factor on new step selection", 'integrator_options.safety', float, 0.9),
        'KISAO_0000540': AlgorithmParameter("maximum factor to increase/decrease step size by in one step",
                                            'integrator_options.ifactor', float, 10.),
        'KISAO_0000539': AlgorithmParameter("minimum factor to increase/decrease step size by in one step",
                                            'integrator_options.dfactor', float, 0.2),
        'KISAO_0000541': AlgorithmParameter("Beta parameter for stabilised step size control", 'integrator_options.beta', float, 0.),
    }),

    'KISAO_0000436': Algorithm("dop835", gillespy2.ODESolver, integrator="dop835", parameters={
        'KISAO_0000211': AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12),
        'KISAO_0000209': AlgorithmParameter("relative tolerance", 'integrator_options.rtol', float, 1e-6),
        'KISAO_0000415': AlgorithmParameter("maximum number of steps", 'integrator_options.nsteps', int, 500),
        'KISAO_0000559': AlgorithmParameter("initial step size", 'integrator_options.first_step', float, 0.0),
        'KISAO_0000467': AlgorithmParameter("maximum step size", 'integrator_options.max_step', float, float("inf")),
        'KISAO_0000538': AlgorithmParameter("safety factor on new step selection", 'integrator_options.safety', float, 0.9),
        'KISAO_0000540': AlgorithmParameter("maximum factor to increase/decrease step size by in one step",
                                            'integrator_options.ifactor', float, 6.),
        'KISAO_0000539': AlgorithmParameter("minimum factor to increase/decrease step size by in one step",
                                            'integrator_options.dfactor', float, 0.333),
        'KISAO_0000541': AlgorithmParameter("Beta parameter for stabilised step size control", 'integrator_options.beta', float, 0.),
    }),

    'KISAO_0000535': Algorithm("vode", gillespy2.ODESolver, integrator="vode", parameters={
        'KISAO_0000211': AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12),
        'KISAO_0000209': AlgorithmParameter("relative tolerance", 'integrator_options.rtol', float, 1e-6),
        'KISAO_0000480': AlgorithmParameter("lower half bandwith", 'integrator_options.lband', int, None),
        'KISAO_0000479': AlgorithmParameter("upper half bandwith", 'integrator_options.uband', int, None),
        'KISAO_0000415': AlgorithmParameter("maximum number of steps", 'integrator_options.nsteps', int, 500),
        'KISAO_0000559': AlgorithmParameter("initial step size", 'integrator_options.first_step', float, 0.0),
        'KISAO_0000485': AlgorithmParameter("minimum step size", 'integrator_options.min_step', float, 0.0),
        'KISAO_0000467': AlgorithmParameter("maximum step size", 'integrator_options.max_step', float, float("inf")),
        'KISAO_0000484': AlgorithmParameter("order", 'integrator_options.order', int, 12),
        'KISAO_0000475': AlgorithmParameter("integration method", 'integrator_options.method', VodeMethod, VodeMethod.adams),
        'KISAO_0000542': AlgorithmParameter("with Jacobian", 'integrator_options.with_jacobian', bool, False),
    }),

    'KISAO_0000536': Algorithm("zvode", gillespy2.ODESolver, integrator="zvode", parameters={
        'KISAO_0000211': AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12),
        'KISAO_0000209': AlgorithmParameter("relative tolerance", 'integrator_options.rtol', float, 1e-6),
        'KISAO_0000480': AlgorithmParameter("lower half bandwith", 'integrator_options.lband', int, None),
        'KISAO_0000479': AlgorithmParameter("upper half bandwith", 'integrator_options.uband', int, None),
        'KISAO_0000415': AlgorithmParameter("maximum number of steps", 'integrator_options.nsteps', int, 500),
        'KISAO_0000559': AlgorithmParameter("initial step size", 'integrator_options.first_step', float, 0.0),
        'KISAO_0000485': AlgorithmParameter("minimum step size", 'integrator_options.min_step', float, 0.0),
        'KISAO_0000467': AlgorithmParameter("maximum step size", 'integrator_options.max_step', float, float("inf")),
        'KISAO_0000484': AlgorithmParameter("order", 'integrator_options.order', int, 12),
        'KISAO_0000475': AlgorithmParameter("integration method", 'integrator_options.method', VodeMethod, VodeMethod.adams),
        'KISAO_0000542': AlgorithmParameter("with Jacobian", 'integrator_options.with_jacobian', bool, False),
    }),

    'KISAO_0000029': Algorithm("SSA", gillespy2.SSACSolver, parameters={
        'KISAO_0000488': AlgorithmParameter("seed", 'seed', int, None),
    }),

    'KISAO_0000039': Algorithm("tau-leaping", gillespy2.TauLeapingSolver, parameters={
        'KISAO_0000488': AlgorithmParameter("seed", 'seed', int, None),
        'KISAO_0000228': AlgorithmParameter("epsilon", 'tau_tol', float, 0.03),
    }),

    'KISAO_0000575': Algorithm("hybrid tau solver", gillespy2.TauHybridSolver, parameters={
        'KISAO_0000488': AlgorithmParameter("seed", 'seed', int, None),
        'KISAO_0000228': AlgorithmParameter("epsilon", 'tau_tol', float, 0.03),
        'KISAO_0000475': AlgorithmParameter("integration method", 'integrator',
                                            HybridTauIntegrationMethod, HybridTauIntegrationMethod.LSODA),
    }),
}
