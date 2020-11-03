""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-26
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from Biosimulations_utils.simulation.data_model import Simulation  # noqa: F401
from Biosimulations_utils.simulation.sedml import modify_xml_model_for_simulation
from Biosimulations_utils.simulator.utils import exec_simulations_in_archive
from gillespy2.sbml import SBMLimport
import enum
import gillespy2
import os
import tempfile

__all__ = [
    'Algorithm', 'AlgorithmParameter', 'VodeMethod', 'HybridTauIntegrationMethod',
    'kisao_algorithm_map',
    'InputError'
    'exec_combine_archive', 'exec_simulation',
]


class Algorithm(object):
    """ Simulation algorithm

    Attributes:
        name (:obj:`str`): name
        solver (:obj:`type`): solver
        solver_args (:obj:`dict`): solver arguments
    """

    def __init__(self, name, solver, parameters=None, **solver_args):
        """
        Args:
            name (:obj:`str`): name
            solver (:obj:`type`): solver
            **solver_args (:obj:`dict`): solver arguments
            parameters (:obj:`dict`): dictionary that maps KiSAO ids to :obj:`AlgorithmParameter`s
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

    def set_value(self, solver_args, value):
        """ Apply the value of a parameter to a data structure of solver arguments

        Args:
            solver_args (:obj:`dict`): solver arguments
            value (:obj:`object`): parameter value
        """
        keys = self.key.split('.')
        for key in keys[0:-1]:
            if key in solver_args:
                nested_solver_args = solver_args[key]
            else:
                nested_solver_args = {}
                solver_args[key] = nested_solver_args
            solver_args = nested_solver_args

        if isinstance(self.data_type, enum.Enum):
            try:
                value = self.data_type(value).name
            except ValueError:
                raise InputError(expression=value,
                                 message="{} option '{}' is not supported".format(self.data_type.__name__, value))

        solver_args[key[-1]] = value


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
        'KISAO_0000540': AlgorithmParameter("maximum factor to increase/decrease step size by in one step",                          'integrator_options.ifactor', float, 10.),
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
    'KISAO_0000028': Algorithm("hybrid tau solver", gillespy2.TauHybridSolver, parameters={
        'KISAO_0000488': AlgorithmParameter("seed", 'seed', int, None),
        'KISAO_0000228': AlgorithmParameter("epsilon", 'tau_tol', float, 0.03),
        'KISAO_0000475': AlgorithmParameter("integration method", 'integrator',
                                            HybridTauIntegrationMethod, HybridTauIntegrationMethod.LSODA),
    }),
}


class InputError(Exception):
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


def exec_combine_archive(archive_file, out_dir):
    """ Execute the SED tasks defined in a COMBINE archive and save the outputs

    Args:
        archive_file (:obj:`str`): path to COMBINE archive
        out_dir (:obj:`str`): directory to store the outputs of the tasks
    """
    exec_simulations_in_archive(archive_file, exec_simulation, out_dir)


def exec_simulation(model_filename, model_sed_urn, simulation, working_dir, out_filename, out_format):
    ''' Execute a simulation and save its results

    Args:
       model_filename (:obj:`str`): path to the model
       model_sed_urn (:obj:`str`): SED URN for the format of the model (e.g., `urn:sedml:language:sbml`)
       simulation (:obj:`Simulation`): simulation
       working_dir (:obj:`str`): directory of the SED-ML file
       out_filename (:obj:`str`): path to save the results of the simulation
       out_format (:obj:`str`): format to save the results of the simulation (e.g., `csv`)
    '''

    # Read the model located at `os.path.join(working_dir, model_filename)` in the format
    # with the SED URN `model_sed_urn`.
    if model_sed_urn != "urn:sedml:language:sbml":
        format = model_sed_urn.split("language:")
        raise InputError(
            expression=format, message="Model language with URN '{}' is not supported".format(model_sed_urn))

    # If necessary, apply the model parameter changes specified by `simulation.model_parameter_changes`
    original_model_abs_filename = os.path.join(working_dir, model_filename)
    if simulation.model_parameter_changes:
        file_handle, model_abs_filename = tempfile.mkstemp(suffix='.xml')
        os.close(file_handle)
        modify_xml_model_for_simulation(
            simulation, original_model_abs_filename, model_abs_filename)
    else:
        model_abs_filename = original_model_abs_filename

    # Convert SBML into a GillesPy2 model
    model = SBMLimport.convert(model_abs_filename, simulation.model.name)[0]
    if model is None:
        raise InputError(expression="model", message=model[1])

    # Load the algorithm specified by `simulation.algorithm`
    algorithm_id = simulation.algorithm.kisao_term.id

    # Todo remove this when Biosimulations_utils/pull/24 is merged
    if (len(algorithm_id) == 7):
        algorithm_id = "KISAO_" + algorithm_id

    algorithm = kisao_algorithm_map.get(algorithm_id, None)
    print(algorithm_id)

    if algorithm is None:

        raise InputError(expression=algorithm_id,
                         message="Algorithm with KISAO id '{}' is not supported".format(algorithm_id))

    # Apply the algorithm parameter changes specified by `simulation.algorithm_parameter_changes`
    algorithm_params = {}
    for change in simulation.algorithm_parameter_changes:
        parameter = algorithm.parameters.get(
            change.parameter.kisao_term.id, None)
        if parameter is None:
            raise InputError(
                expression=change.parameter.kisao_term.id,
                message="Algorithm parameter with KiSAO id '{}' is not supported".format(change.parameter.kisao_term.id))
        parameter.set_value(algorithm_params, change.value)

    # Validate that start time is 0 because this is the only option that GillesPy2 supports
    if simulation.start_time > 0:
        raise InputError(expression=simulation.start_time,
                         message='Start time must be at least 0')
    print(simulation.end_time)
    # Simulate the model from `simulation.start_time` to `simulation.end_time` and record `simulation.num_time_points` + 1 time points
    increment = (simulation.end_time - simulation.start_time) / \
        simulation.num_time_points
    results = model.run(algorithm.solver, **algorithm.solver_args,
                        **algorithm_params, t=simulation.end_time)
    print(results)

    # TODO: ignore all time points before `simulation.start_time`

    # TODO: Save a report of the results of the simulation with `simulation.num_time_points` time points
    # beginning at `simulation.output_start_time` to `out_filename` in `out_format` format.
    # This should save all of the variables specified by `simulation.model.variables`.


if __name__ == "__main__":

    exec_simulation("tests/fixtures/BIOMD0000000028.xml",
                    "urn:sedml:language:sbml", "simulation", "", "", "")
