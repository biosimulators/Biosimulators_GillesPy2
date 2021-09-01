""" Methods for using GillesPy2 to execute SED tasks in COMBINE archives and save their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Date: 2020-12-09
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  Variable, Symbol)
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import gillespy2
import math
import numpy
import functools

__all__ = [
    'exec_sedml_docs_in_combine_archive', 'exec_sed_task',
]


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_task(task, variables, log=None, config=None):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task
       config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`ValueError`: if the task or an aspect of the task is not valid, or the requested output variables
            could not be recorded
        :obj:`NotImplementedError`: if the task is not of a supported type or involves an unsuported feature
    '''
    config = config or get_config()

    if config.LOG and not log:
        log = TaskLog()

    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(task.model.changes, ()),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(task.model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(task.simulation, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(task.simulation),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    target_x_paths_ids = validation.validate_variable_xpaths(variables, task.model.source, attr='id')

    # Read the SBML-encoded model located at `task.model.source`
    model, errors = gillespy2.import_SBML(task.model.source)
    if model is None or errors:
        raise ValueError('Model at {} could not be imported:\n  - {}'.format(
            task.model.source, '\n  - '.join(message for message, code in errors)))

    # Load the algorithm specified by `simulation.algorithm`
    simulation = task.simulation
    algorithm_kisao_id = simulation.algorithm.kisao_id
    algorithm_substitution_policy = get_algorithm_substitution_policy(config=config)
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        algorithm_kisao_id, KISAO_ALGORITHM_MAP.keys(),
        substitution_policy=algorithm_substitution_policy)
    algorithm = KISAO_ALGORITHM_MAP[exec_kisao_id]

    solver = algorithm.solver
    if solver == gillespy2.SSACSolver and (model.get_all_events() or model.get_all_assignment_rules()):
        solver = gillespy2.NumPySSASolver

    # Apply the algorithm parameter changes specified by `simulation.algorithm.parameter_changes`
    algorithm_params = {}
    if exec_kisao_id == algorithm_kisao_id:
        for change in simulation.algorithm.changes:
            parameter = algorithm.parameters.get(change.kisao_id, None)
            if parameter:
                try:
                    parameter.set_value(algorithm_params, change.new_value)
                except (NotImplementedError, ValueError) as exception:
                    if (
                        ALGORITHM_SUBSTITUTION_POLICY_LEVELS[algorithm_substitution_policy]
                        <= ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                    ):
                        raise
                    else:
                        warn('Unsuported value `{}` for algorithm parameter `{}` was ignored:\n  {}'.format(
                            change.new_value, change.kisao_id, str(exception).replace('\n', '\n  ')),
                            BioSimulatorsWarning)
            else:
                if (
                    ALGORITHM_SUBSTITUTION_POLICY_LEVELS[algorithm_substitution_policy]
                    <= ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                ):
                    msg = "".join([
                        "Algorithm parameter with KiSAO id '{}' is not supported. ".format(change.kisao_id),
                        "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                            '{}: {}'.format(kisao_id, parameter.name) for kisao_id, parameter in algorithm.parameters.items())),
                    ])
                    raise NotImplementedError(msg)
                else:
                    msg = "".join([
                        "Algorithm parameter with KiSAO id '{}' was ignored because it is not supported. ".format(change.kisao_id),
                        "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                            '{}: {}'.format(kisao_id, parameter.name) for kisao_id, parameter in algorithm.parameters.items())),
                    ])
                    warn(msg, BioSimulatorsWarning)

    # Validate that start time is 0 because this is the only option that GillesPy2 supports
    if simulation.initial_time < 0:
        raise NotImplementedError(
            'Negative initial simulation time {} is not supported. Initial time must be >= 0.'.format(simulation.initial_time))

    # set the simulation time span
    number_of_points = (simulation.output_end_time - simulation.initial_time) / \
        (simulation.output_end_time - simulation.output_start_time) * simulation.number_of_points
    if number_of_points != math.floor(number_of_points):
        raise NotImplementedError('Time course must specify an integer number of time points')
    number_of_points = int(number_of_points)
    timespan = numpy.linspace(simulation.initial_time, simulation.output_end_time, number_of_points + 1)
    model.timespan(timespan)

    # determine allowed variable targets
    predicted_ids = list(model.get_all_species().keys()) + list(model.get_all_parameters().keys())
    unpredicted_symbols = set()
    unpredicted_targets = set()
    for variable in variables:
        if variable.symbol:
            if variable.symbol != Symbol.time:
                unpredicted_symbols.add(variable.symbol)

        else:
            if target_x_paths_ids[variable.target] not in predicted_ids:
                unpredicted_targets.add(variable.target)

    if unpredicted_symbols:
        raise NotImplementedError("".join([
            "The following variable symbols are not supported:\n  - {}\n\n".format(
                '\n  - '.join(sorted(unpredicted_symbols)),
            ),
            "Symbols must be one of the following:\n  - {}".format(Symbol.time),
        ]))

    if unpredicted_targets:
        raise ValueError(''.join([
            'The following variable targets could not be recorded:\n  - {}\n\n'.format(
                '\n  - '.join(sorted(unpredicted_targets)),
            ),
            'Targets must have one of the following ids:\n  - {}'.format(
                '\n  - '.join(sorted(predicted_ids)),
            ),
        ]))

    # Simulate the model from ``simulation.start_time`` to ``simulation.output_end_time``
    # and record ``simulation.number_of_points`` + 1 time points
    results_dict = model.run(solver, **algorithm.solver_args, **algorithm_params)[0]

    # transform the results to an instance of :obj:`VariableResults`
    variable_results = VariableResults()
    parameters = model.get_all_parameters()
    for variable in variables:
        if variable.symbol:
            variable_results[variable.id] = results_dict['time'][-(simulation.number_of_points + 1):]

        elif variable.target:
            sbml_id = target_x_paths_ids[variable.target]
            dynamics = results_dict.get(sbml_id, None)
            if dynamics is None:
                variable_results[variable.id] = numpy.full((simulation.number_of_points + 1,), parameters[sbml_id].value)
            else:
                variable_results[variable.id] = dynamics[-(simulation.number_of_points + 1):]

    # log action
    if config.LOG:
        log.algorithm = exec_kisao_id
        log.simulator_details = {
            'method': solver.__module__ + '.' + solver.__name__,
            'arguments': dict(**algorithm.solver_args, **algorithm_params),
        }

    # return results and log
    return variable_results, log
