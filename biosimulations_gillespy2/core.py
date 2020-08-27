""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Author name <email@organization>
:Date: YYYY-MM-DD
:Copyright: YYYY, Owner
:License: <License, e.g., MIT>
"""

from Biosimulations_utils.simulation.data_model import Simulation  # noqa: F401
from Biosimulations_utils.simulator.utils import exec_simulations_in_archive
from Biosimulations_utils.simulation.sedml import modify_xml_model_for_simulation
import os
from gillespy2.sbml import SBMLimport
from numpy import linspace
# Supported Solvers
from gillespy2 import TauLeapingSolver
from gillespy2 import ODESolver
from gillespy2 import SSACSolver

__all__ = ['exec_combine_archive', 'exec_simulation']


params = {
    "0000211": 1e-12,
    "0000209": 1e-6,
    "0000480": None,
    "0000479": None,
    "0000415": 500,
    "0000483": 0.0,
    "0000485": 0.0,
    "0000467": float('inf'),
    "000219": 12,
    "0000220": 5,
    "0000488": None,
    "0000228": 0.03,
}
ode_params = {
    "atol": "0000211",
    "rtol": "0000209",
    "lband": "0000480",
    "uband": "0000479",
    "nsteps": "0000415",
    "first_step": "0000483",
    "min_step": "0000485",
    "max_step": "0000467",
    "max_order_ns": "000219",
    "max_order_s": "0000220"
}
ssac_params = {
    "seed": "0000488",
}
tau_leap_params = {
    "seed": "0000488",
    "tau_tol": "0000228"
}
algorithm_map = {
    "0000088": (ODESolver, ode_params),
    "0000029": (SSACSolver, ssac_params),
    "0000039": (TauLeapingSolver, tau_leap_params)
}


class InputError(Exception):
    def __init__(self, message):

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
        raise InputError(expression=format, message="Unsupported")

    file_name = os.path.join(working_dir, model_filename)

    new_file_name = file_name + "_modified"

    # Apply the model parameter changes specified by `simulation.model_parameter_changes`
    modify_xml_model_for_simulation(
        simulation, model_filename, new_file_name, "sbml")

    # Create a gilespy2 model
    model_error = SBMLimport.convert(new_file_name, simulation.model.name)
    model = model_error[0]

    if (model is None):
        print(model_error[1])
        raise InputError(message=model_error[1])

    # Load the algorithm specified by `simulation.algorithm`
    algorithm_id = simulation.algorithm.kisao_term.id
    algorithm = algorithm_map.get(algorithm_id)

    if algorithm is None:
        raise InputError("Unsupported Algorithm")
    # Apply the algorithm parameter changes specified by `simulation.algorithm_parameter_changes`
    for algorithm_change in simulation.algorithm_parameter_changes:

        kisao_id = algorithm_change.parameter.kisao_term.id

        if(kisao_id in params.keys()):
            params[kisao_id] = algorithm_change.value
        else:
            raise InputError("parameter with kisao id" +
                             kisao_id + " is not supported")

    solver = algorithm[0]
    mapping = algorithm[1]
    solver_params = {}
    if solver == ODESolver:
        integrator_options = {}
        for param, kisao in mapping.items():
            integrator_options[param] = params[kisao]
        solver_params["integrator_options"] = integrator_options
    else:
        for param, kisao in mapping.items():
            solver_params[param] = params[kisao]
    # Simulate the model from `simulation.start_time` to `simulation.end_time`
    start_time = simulation.start_time
    end_time = simulation.end_time
    num_time_points = simulation.num_time_points
    timespan = linspace(start_time, end_time, num=num_time_points)
    model.timespan(timespan)
    # Save a report of the results of the simulation with `simulation.num_time_points` time points
    # beginning at `simulation.output_start_time` to `out_filename` in `out_format` format.
    # This should save all of the variables specified by `simulation.model.variables`.

    results = model.run(solver, **solver_params)
    print(simulation.model.variables)
    for variable in simulation.model.variables:
        target = variable.target
        print(target)
    print(results)


if __name__ == "__main__":
    class Model:
        name = "testModel"

    class Simulation:
        model = Model()

    exec_simulation("tests/fixtures/BIOMD0000000028.xml",
                    "urn:sedml:language:sbml", Simulation(), "", "", "")
