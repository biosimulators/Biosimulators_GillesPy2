""" Methods for using GillesPy2 to execute SED tasks in COMBINE archives and save their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Date: 2020-12-09
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import kisao_algorithm_map
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.plot.data_model import PlotFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  Variable, Symbol)
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.utils.core import raise_errors_warnings
import gillespy2
import math
import numpy
import functools

__all__ = [
    'exec_sedml_docs_in_combine_archive', 'exec_sed_task',
]


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                       report_formats=None, plot_formats=None,
                                       bundle_outputs=None, keep_individual_outputs=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., csv or h5)
        plot_formats (:obj:`list` of :obj:`PlotFormat`, optional): report format (e.g., pdf)
        bundle_outputs (:obj:`bool`, optional): if :obj:`True`, bundle outputs into archives for reports and plots
        keep_individual_outputs (:obj:`bool`, optional): if :obj:`True`, keep individual output files

    Returns:
        :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      report_formats=report_formats,
                                      plot_formats=plot_formats,
                                      bundle_outputs=bundle_outputs,
                                      keep_individual_outputs=keep_individual_outputs)


def exec_sed_task(task, variables, log=None):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`ValueError`: if the task or an aspect of the task is not valid, or the requested output variables
            could not be recorded
        :obj:`NotImplementedError`: if the task is not of a supported type or involves an unsuported feature
    '''
    log = log or TaskLog()

    model = task.model
    sim = task.simulation

    raise_errors_warnings(validation.validate_task(task),
                          error_summary='Task `{}` is invalid.'.format(task.id))
    raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                          error_summary='Language for model `{}` is not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_change_types(task.model.changes, ()),
                          error_summary='Changes for model `{}` are not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_changes(task.model),
                          error_summary='Changes for model `{}` are invalid.'.format(model.id))
    raise_errors_warnings(validation.validate_simulation_type(task.simulation, (UniformTimeCourseSimulation, )),
                          error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
    raise_errors_warnings(validation.validate_simulation(task.simulation),
                          error_summary='Simulation `{}` is invalid.'.format(sim.id))
    raise_errors_warnings(validation.validate_data_generator_variables(variables),
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

    # set the simulation time span
    number_of_points = (simulation.output_end_time - simulation.initial_time) / \
        (simulation.output_end_time - simulation.output_start_time) * simulation.number_of_points
    if number_of_points != math.floor(number_of_points):
        raise NotImplementedError('Time course must specify an integer number of time points')
    number_of_points = int(number_of_points)
    model.timespan(numpy.linspace(simulation.initial_time, simulation.output_end_time, number_of_points + 1))

    # determine allowed variable targets
    predicted_ids = list(model.get_all_species().keys())
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
    for variable in variables:
        if variable.symbol:
            variable_results[variable.id] = results_dict['time'][-(simulation.number_of_points + 1):]

        elif variable.target:
            variable_results[variable.id] = results_dict[target_x_paths_ids[variable.target]][-(simulation.number_of_points + 1):]

    # log action
    log.algorithm = algorithm_kisao_id
    log.simulator_details = {
        'method': solver.__module__ + '.' + solver.__name__,
        'arguments': dict(**algorithm.solver_args, **algorithm_params),
    }

    # return results and log
    return variable_results, log
