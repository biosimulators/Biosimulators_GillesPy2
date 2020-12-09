""" Methods for using GillesPy2 to execute SED tasks in COMBINE archives and save their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Date: 2020-12-09
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.report.data_model import ReportFormat, DataGeneratorVariableResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  DataGeneratorVariable, DataGeneratorVariableSymbol)
import enum
import gillespy2
import math
import numpy
import os
import re

__all__ = [
    'Algorithm', 'AlgorithmParameter', 'VodeMethod', 'HybridTauIntegrationMethod',
    'kisao_algorithm_map',
    'exec_sedml_docs_in_combine_archive', 'exec_sed_task',
]


class Algorithm(object):
    """ Simulation algorithm

    Attributes:
        name (:obj:`str`): name
        solver (:obj:`type`): solver
        solver_args (:obj:`dict`): solver arguments
        parameters (:obj:`dict`): dictionary that maps KiSAO ids to :obj:`AlgorithmParameter`s
    """

    def __init__(self, name, solver, parameters=None, **solver_args):
        """
        Args:
            name (:obj:`str`): name
            solver (:obj:`type`): solver
            **solver_args (:obj:`dict`): solver arguments
            parameters (:obj:`dict`, optional): dictionary that maps KiSAO ids to :obj:`AlgorithmParameter`s
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

        if self.data_type == bool:
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

    'KISAO_0000028': Algorithm("hybrid tau solver", gillespy2.TauHybridSolver, parameters={
        'KISAO_0000488': AlgorithmParameter("seed", 'seed', int, None),
        'KISAO_0000228': AlgorithmParameter("epsilon", 'tau_tol', float, 0.03),
        'KISAO_0000475': AlgorithmParameter("integration method", 'integrator',
                                            HybridTauIntegrationMethod, HybridTauIntegrationMethod.LSODA),
    }),
}


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, report_formats=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., CSV or HDF5)
    """
    exec_sedml_docs_in_archive(archive_filename, exec_sed_task, out_dir,
                               apply_xml_model_changes=True,
                               report_formats=report_formats)


def exec_sed_task(task, variables):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`DataGeneratorVariable`): variables that should be recorded

    Returns:
        :obj:`DataGeneratorVariableResults`: results of variables

    Raises:
        :obj:`ValueError`: if the task or an aspect of the task is not valid, or the requested output variables
            could not be recorded
        :obj:`NotImplementedError`: if the task is not of a supported type or involves an unsuported feature
    '''
    # check that task is an instance of Task
    if not isinstance(task, Task):
        raise NotImplementedError('Task type {} is not supported'.format(task.__class__.__name__))

    # check that task has model
    if not task.model:
        raise ValueError('Task must have a model')

    # check that model is encoded in SBML
    if not task.model.language or not re.match('^{}($|:)'.format(ModelLanguage.SBML.value), task.model.language):
        raise NotImplementedError("Model language {} is not supported. Model language must be '{}'.".format(
            task.model.language, ModelLanguage.SBML.value))

    # check that model parameter changes have already been applied (because handled by :obj:`exec_sedml_docs_in_archive`)
    if task.model.changes:
        raise NotImplementedError('Model changes are not supported')

    # check that task has model
    simulation = task.simulation
    if not simulation:
        raise ValueError('Task must have a simulation')

    # check that simulation is a time course simulation
    if not isinstance(simulation, UniformTimeCourseSimulation):
        raise NotImplementedError('Simulation type {} is not supported. Simulation must be an instance of {}.'.format(
            simulation.__class__.__name__, UniformTimeCourseSimulation.__name__))

    # Read the model located at `task.model.source` in the format
    # with the SED URN `model_language_urn`.
    # Convert SBML into a GillesPy2 model
    if not task.model.source or not os.path.isfile(task.model.source):
        raise FileNotFoundError("Model source '{}' must be a file".format(task.model.source or ''))

    model, errors = gillespy2.import_SBML(task.model.source)
    if model is None or errors:
        raise ValueError('Model at {} could not be imported:\n  - {}'.format(
            task.model.source, '\n  - '.join(message for message, code in errors)))

    # Load the algorithm specified by `simulation.algorithm`
    if not simulation.algorithm:
        raise ValueError('Simulation must have an algorithm')

    algorithm_kisao_id = simulation.algorithm.kisao_id
    algorithm = kisao_algorithm_map.get(algorithm_kisao_id, None)
    if algorithm is None:
        raise NotImplementedError("".join([
            "Algorithm with KiSAO id '{}' is not supported. ".format(algorithm_kisao_id),
            "Algorithm must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                '{}: {} ({})'.format(kisao_id, algorithm.name, algorithm.solver.__name__)
                for kisao_id, algorithm in kisao_algorithm_map.items())),
        ]))

    solver = algorithm.solver
    if solver == gillespy2.SSACSolver and (model.get_all_events() or model.get_all_assignment_rules()):
        solver = gillespy2.NumPySSASolver

    # Apply the algorithm parameter changes specified by `simulation.algorithm.parameter_changes`
    algorithm_params = {}
    for change in simulation.algorithm.changes:
        parameter = algorithm.parameters.get(change.kisao_id, None)
        if parameter is None:
            raise NotImplementedError("".join([
                "Algorithm parameter with KiSAO id '{}' is not supported. ".format(change.kisao_id),
                "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                    '{}: {}'.format(kisao_id, parameter.name) for kisao_id, parameter in algorithm.parameters.items())),
            ]))
        parameter.set_value(algorithm_params, change.new_value)

    # Validate that start time is 0 because this is the only option that GillesPy2 supports
    if simulation.initial_time != 0:
        raise NotImplementedError('Initial simulation time {} is not supported. Initial time must be 0.'.format(simulation.initial_time))

    if simulation.output_start_time < simulation.initial_time:
        raise ValueError('Output start time {} must be at least the initial time {}.'.format(
            simulation.output_start_time, simulation.initial_time))

    if simulation.output_end_time < simulation.output_start_time:
        raise ValueError('Output end time {} must be at least the output start time {}.'.format(
            simulation.output_end_time, simulation.output_start_time))

    # set the simulation time span
    number_of_points = (simulation.output_end_time - simulation.initial_time) / \
        (simulation.output_end_time - simulation.output_start_time) * simulation.number_of_points
    if number_of_points != math.floor(number_of_points):
        raise NotImplementedError('Time course must specify an integer number of time points')
    number_of_points = int(number_of_points)
    model.timespan(numpy.linspace(simulation.initial_time, simulation.output_end_time, number_of_points + 1))

    # Simulate the model from ``simulation.start_time`` to ``simulation.output_end_time``
    # and record ``simulation.number_of_points`` + 1 time points
    results_dict = model.run(solver, **algorithm.solver_args, **algorithm_params)[0]

    # transform the results to an instance of :obj:`DataGeneratorVariableResults`
    variable_results = DataGeneratorVariableResults()
    unpredicted_symbols = []
    unpredicted_targets = []
    for variable in variables:
        if (variable.symbol and variable.target) or (not variable.symbol and not variable.target):
            raise ValueError('Variable must define a symbol or target')

        if variable.symbol:
            if variable.symbol == DataGeneratorVariableSymbol.time:
                variable_results[variable.id] = results_dict['time'][-(simulation.number_of_points + 1):]
            else:
                unpredicted_symbols.append(variable.symbol)

        elif variable.target:
            match = (
                re.match(r"^/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species\[@id='(.*?)'\]$", variable.target) or
                re.match(r'^/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species\[@id="(.*?)"\]$', variable.target)
            )
            if match and match.group(1) in results_dict:
                variable_results[variable.id] = results_dict[match.group(1)][-(simulation.number_of_points + 1):]
            else:
                unpredicted_targets.append(variable.target)

    if unpredicted_symbols:
        raise NotImplementedError("".join([
            "The following variable symbols are not supported:\n  - {}\n\n".format(
                '\n  - '.join(sorted(unpredicted_symbols)),
            ),
            "Symbols must be one of the following:\n  - {}".format(DataGeneratorVariableSymbol.time),
        ]))

    if unpredicted_targets:
        raise ValueError(''.join([
            'The following variable targets could not be recorded:\n  - {}\n\n'.format(
                '\n  - '.join(sorted(unpredicted_targets)),
            ),
            'Targets must be one of the following:\n  - {}'.format(
                '\n  - '.join("/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='{}']".format(id)
                              for id in sorted(results_dict.keys()) if id != 'time'),
            )
        ]))

    # return results
    return variable_results
