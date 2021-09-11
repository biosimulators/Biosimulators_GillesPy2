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
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog, StandardOutputErrorCapturerLevel  # noqa: F401
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange, UniformTimeCourseSimulation,  # noqa: F401
                                                  Variable, Symbol)
from biosimulators_utils.sedml.exec import exec_sed_doc as base_exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import gillespy2
import lxml
import math
import numpy

__all__ = [
    'exec_sedml_docs_in_combine_archive', 'exec_sed_doc', 'exec_sed_task', 'preprocess_sed_task',
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
    return exec_sedml_docs_in_archive(exec_sed_doc, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_doc(doc, working_dir, base_out_path, rel_out_path=None,
                 apply_xml_model_changes=True,
                 log=None, indent=0, pretty_print_modified_xml_models=False,
                 log_level=StandardOutputErrorCapturerLevel.c, config=None):
    """ Execute the tasks specified in a SED document and generate the specified outputs

    Args:
        doc (:obj:`SedDocument` or :obj:`str`): SED document or a path to SED-ML file which defines a SED document
        working_dir (:obj:`str`): working directory of the SED document (path relative to which models are located)

        base_out_path (:obj:`str`): path to store the outputs

            * CSV: directory in which to save outputs to files
              ``{base_out_path}/{rel_out_path}/{report.id}.csv``
            * HDF5: directory in which to save a single HDF5 file (``{base_out_path}/reports.h5``),
              with reports at keys ``{rel_out_path}/{report.id}`` within the HDF5 file

        rel_out_path (:obj:`str`, optional): path relative to :obj:`base_out_path` to store the outputs
        apply_xml_model_changes (:obj:`bool`, optional): if :obj:`True`, apply any model changes specified in the SED-ML file before
            calling :obj:`task_executer`.
        log (:obj:`SedDocumentLog`, optional): log of the document
        indent (:obj:`int`, optional): degree to indent status messages
        pretty_print_modified_xml_models (:obj:`bool`, optional): if :obj:`True`, pretty print modified XML models
        log_level (:obj:`StandardOutputErrorCapturerLevel`, optional): level at which to log output
        config (:obj:`Config`, optional): BioSimulators common configuration
        simulator_config (:obj:`SimulatorConfig`, optional): tellurium configuration

    Returns:
        :obj:`tuple`:

            * :obj:`ReportResults`: results of each report
            * :obj:`SedDocumentLog`: log of the document
    """
    return base_exec_sed_doc(exec_sed_task, doc, working_dir, base_out_path,
                             rel_out_path=rel_out_path,
                             apply_xml_model_changes=apply_xml_model_changes,
                             log=log,
                             indent=indent,
                             pretty_print_modified_xml_models=pretty_print_modified_xml_models,
                             log_level=log_level,
                             config=config)


def exec_sed_task(task, variables, preprocessed_task=None, log=None, config=None):
    ''' Execute a task and save its results

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        preprocessed_task (:obj:`dict`, optional): preprocessed information about the task, including possible
            model changes and variables. This can be used to avoid repeatedly executing the same initialization
            for repeated calls to this method.
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

    if preprocessed_task is None:
        preprocessed_task = preprocess_sed_task(task, variables, config=config)

    sim = task.simulation

    # get model
    model = preprocessed_task['model']['model']

    # modify model
    raise_errors_warnings(validation.validate_model_change_types(task.model.changes, (ModelAttributeChange,)),
                          error_summary='Changes for model `{}` are not supported.'.format(task.model.id))
    change_target_model_obj_map = preprocessed_task['model']['change_target_model_obj_map']
    for change in task.model.changes:
        model_obj = change_target_model_obj_map[change.target]
        new_value = float(change.new_value)
        if isinstance(model_obj, gillespy2.core.parameter.Parameter):
            model_obj.value = new_value
        else:
            model_obj.initial_value = new_value

    # Validate that start time is 0 because this is the only option that GillesPy2 supports
    if sim.initial_time < 0:
        raise NotImplementedError(
            'Negative initial simulation time {} is not supported. Initial time must be >= 0.'.format(sim.initial_time))

    # set the simulation time span
    number_of_points = (sim.output_end_time - sim.initial_time) / \
        (sim.output_end_time - sim.output_start_time) * sim.number_of_points
    if number_of_points != math.floor(number_of_points):
        raise NotImplementedError('Time course must specify an integer number of time points')
    number_of_points = int(number_of_points)
    timespan = numpy.linspace(sim.initial_time, sim.output_end_time, number_of_points + 1)
    model.timespan(timespan)

    # Simulate the model from ``sim.start_time`` to ``sim.output_end_time``
    # and record ``sim.number_of_points`` + 1 time points
    solver = preprocessed_task['simulation']['solver']
    solver_args = preprocessed_task['simulation']['solver_args']
    results_dict = model.run(solver, **solver_args)[0]

    # transform the results to an instance of :obj:`VariableResults`
    variable_target_sbml_id_map = preprocessed_task['model']['variable_target_sbml_id_map']
    variable_results = VariableResults()
    parameters = model.get_all_parameters()
    for variable in variables:
        if variable.symbol:
            variable_results[variable.id] = results_dict['time'][-(sim.number_of_points + 1):]

        elif variable.target:
            sbml_id = variable_target_sbml_id_map[variable.target]
            dynamics = results_dict.get(sbml_id, None)
            if dynamics is None:
                variable_results[variable.id] = numpy.full((sim.number_of_points + 1,), parameters[sbml_id].value)
            else:
                variable_results[variable.id] = dynamics[-(sim.number_of_points + 1):]

    # log action
    if config.LOG:
        log.algorithm = preprocessed_task['simulation']['algorithm_kisao_id']
        log.simulator_details = {
            'method': solver.__module__ + '.' + solver.__name__,
            'arguments': solver_args,
        }

    # return results and log
    return variable_results, log


def preprocess_sed_task(task, variables, config=None):
    """ Preprocess a SED task, including its possible model changes and variables. This is useful for avoiding
    repeatedly initializing tasks on repeated calls of :obj:`exec_sed_task`.

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`dict`: preprocessed information about the task
    """
    config = config or get_config()

    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(task.model.id))
        raise_errors_warnings(validation.validate_model_change_types(task.model.changes, (ModelAttributeChange,)),
                              error_summary='Changes for model `{}` are not supported.'.format(task.model.id))
        raise_errors_warnings(*validation.validate_model_changes(task.model),
                              error_summary='Changes for model `{}` are invalid.'.format(task.model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    model_etree = lxml.etree.parse(task.model.source)
    change_target_sbml_id_map = validation.validate_target_xpaths(task.model.changes, model_etree, attr='id')
    variable_target_sbml_id_map = validation.validate_target_xpaths(variables, model_etree, attr='id')

    # Read the SBML-encoded model located at `task.model.source`
    model, errors = gillespy2.import_SBML(task.model.source)
    if model is None or errors:
        raise ValueError('Model at {} could not be imported:\n  - {}'.format(
            task.model.source, '\n  - '.join(message for message, code in errors)))

    # preprocess model changes
    parameters = model.get_all_parameters()
    species = model.get_all_species()
    change_target_model_obj_map = {}
    invalid_changes = []
    for change in task.model.changes:
        sbml_id = change_target_sbml_id_map[change.target]
        model_obj = parameters.get(sbml_id, species.get(sbml_id, None))
        if model_obj is None:
            invalid_changes.append(change.target)
        else:
            change_target_model_obj_map[change.target] = model_obj

    if invalid_changes:
        raise ValueError(''.join([
            'The following model targets cannot be changed:\n  - {}\n\n'.format(
                '\n  - '.join(sorted(invalid_changes)),
            ),
            'Model change targets must have one of the following SBML ids:\n  - {}'.format(
                '\n  - '.join(sorted(list(parameters.keys()) + list(species.keys()))),
            ),
        ]))

    # Load the algorithm specified by `sim.algorithm`
    algorithm_substitution_policy = get_algorithm_substitution_policy(config=config)
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        sim.algorithm.kisao_id, KISAO_ALGORITHM_MAP.keys(),
        substitution_policy=algorithm_substitution_policy)
    algorithm = KISAO_ALGORITHM_MAP[exec_kisao_id]

    solver = algorithm.solver
    if solver == gillespy2.SSACSolver and (model.get_all_events() or model.get_all_assignment_rules()):
        solver = gillespy2.NumPySSASolver

    # Apply the algorithm parameter changes specified by `sim.algorithm.parameter_changes`
    algorithm_params = {}
    if exec_kisao_id == sim.algorithm.kisao_id:
        for change in sim.algorithm.changes:
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

    # determine allowed variable targets
    predicted_ids = list(species.keys()) + list(parameters.keys())
    unpredicted_symbols = set()
    unpredicted_targets = set()
    for variable in variables:
        if variable.symbol:
            if variable.symbol != Symbol.time:
                unpredicted_symbols.add(variable.symbol)

        else:
            if variable_target_sbml_id_map[variable.target] not in predicted_ids:
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
            'Targets must have one of the following SBML ids:\n  - {}'.format(
                '\n  - '.join(sorted(predicted_ids)),
            ),
        ]))

    # return preprocessed information about the task
    return {
        'model': {
            'model': model,
            'change_target_model_obj_map': change_target_model_obj_map,
            'variable_target_sbml_id_map': variable_target_sbml_id_map,
        },
        'simulation': {
            'algorithm_kisao_id': exec_kisao_id,
            'solver': solver,
            'solver_args': dict(**algorithm.solver_args, **algorithm_params),
        }
    }
