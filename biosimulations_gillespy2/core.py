""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Author name <email@organization>
:Date: YYYY-MM-DD
:Copyright: YYYY, Owner
:License: <License, e.g., MIT>
"""

from Biosimulations_utils.simulation.data_model import Simulation  # noqa: F401
from Biosimulations_utils.simulator.utils import exec_simulations_in_archive
import os
from gillespy2.sbml import SBMLimport
from gillespy2 import TauHybridSolver

__all__ = ['exec_combine_archive', 'exec_simulation']
algorithm_map={
    
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
        format=model_sed_urn.split("language:")
        raise InputError(expression=format, message="Unsupported")

    
    file_name=os.path.join(working_dir,model_filename)
    
    # Create a gilespy model
    model= SBMLimport.convert(file_name,simulation.model.name)[0]
    if(model is None):
        raise InputError(expression="model", message=model[1])
    
    # Apply the model parameter changes specified by `simulation.model_parameter_changes`
    if(simulation.model_parameter_changes):

        for change in simulation.model_parameter_changes:
            target= change.parameter.target
            target=target.split(".")[1]
            # TODO see if this cast should be handled by utils. The value is set as an expression
            value= str(change.value)
            model.set_parameter(target,value)



    # Load the algorithm specified by `simulation.algorithm`
    algorithm_id= simulation.algorithm.kisao_term.id
    algorithm= algorithm_map.get(algorithm_id)
    

    # Apply the algorithm parameter changes specified by `simulation.algorithm_parameter_changes`

    # Simulate the model from `simulation.start_time` to `simulation.end_time`

    # Save a report of the results of the simulation with `simulation.num_time_points` time points
    # beginning at `simulation.output_start_time` to `out_filename` in `out_format` format.
    # This should save all of the variables specified by `simulation.model.variables`.

    results = model.run(TauHybridSolver)
    print(results)
    
if __name__ == "__main__":
  
    exec_simulation("tests/fixtures/BIOMD0000000028.xml", "urn:sedml:language:sbml", "simulation","", "","")