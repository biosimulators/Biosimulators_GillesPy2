""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_gillespy2 import __main__
from biosimulators_gillespy2 import core
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.simulator.exec import exec_sedml_docs_in_archive_with_containerized_simulator
from biosimulators_utils.sedml import data_model as sedml_data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.sedml.utils import append_all_nested_children_to_doc
from biosimulators_utils.utils.core import are_lists_equal
from unittest import mock
import datetime
import dateutil.tz
import numpy.testing
import os
import shutil
import tempfile
import unittest


class TestCase(unittest.TestCase):
    DOCKER_IMAGE = 'ghcr.io/biosimulators/biosimulators_gillespy2/gillespy2:latest'

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_exec_sed_task(self):
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml'),
                language=sedml_data_model.ModelLanguage.SBML.value,
                changes=[],
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000029',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000488',
                            new_value='10',
                        ),
                    ],
                ),
                initial_time=0.,
                output_start_time=10.,
                output_end_time=20.,
                number_of_points=20,
            ),
        )

        variables = [
            sedml_data_model.DataGeneratorVariable(id='time', symbol=sedml_data_model.DataGeneratorVariableSymbol.time),
            sedml_data_model.DataGeneratorVariable(id='BE', target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']"),
            sedml_data_model.DataGeneratorVariable(id='BUD', target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BUD"]'),
            sedml_data_model.DataGeneratorVariable(id='Cdc20', target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']"),
        ]

        variable_results = core.exec_sed_task(task, variables)

        self.assertTrue(sorted(variable_results.keys()), sorted([var.id for var in variables]))
        self.assertEqual(variable_results[variables[0].id].shape, (task.simulation.number_of_points + 1,))
        numpy.testing.assert_equal(
            variable_results['time'],
            numpy.linspace(task.simulation.output_start_time, task.simulation.output_end_time, task.simulation.number_of_points + 1),
        )

    def test_exec_sed_task_errors(self):
        task = None
        variables = []
        with self.assertRaisesRegex(NotImplementedError, 'is not supported'):
            core.exec_sed_task(task, variables)

        task = sedml_data_model.Task()
        with self.assertRaisesRegex(ValueError, 'must have a model'):
            core.exec_sed_task(task, variables)

        task.model = sedml_data_model.Model()
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Model language must be'):
            core.exec_sed_task(task, variables)

        task.model.language = 'sbml'
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Model language must be'):
            core.exec_sed_task(task, variables)

        task.model.language = sedml_data_model.ModelLanguage.SBML
        task.model.changes = [sedml_data_model.ModelAttributeChange()]
        with self.assertRaisesRegex(NotImplementedError, 'changes are not supported'):
            core.exec_sed_task(task, variables)

        task.model.changes = []
        with self.assertRaisesRegex(ValueError, 'must have a simulation'):
            core.exec_sed_task(task, variables)

        task.simulation = mock.Mock()
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Simulation must be'):
            core.exec_sed_task(task, variables)

        task.simulation = sedml_data_model.UniformTimeCourseSimulation(
            initial_time=10.,
            output_start_time=-10.,
            output_end_time=-20.,
            number_of_points=10,
        )
        with self.assertRaisesRegex(FileNotFoundError, 'must be a file'):
            core.exec_sed_task(task, variables)

        task.model.source = os.path.join(self.dirname, 'invalid-model.xml')
        with open(task.model.source, 'w') as file:
            file.write('!')
        with self.assertRaisesRegex(ValueError, 'could not be imported'):
            core.exec_sed_task(task, variables)

        task.model.source = os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml')
        with self.assertRaisesRegex(ValueError, 'must have an algorithm'):
            core.exec_sed_task(task, variables)

        task.simulation.algorithm = sedml_data_model.Algorithm()
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Algorithm must'):
            core.exec_sed_task(task, variables)

        task.simulation.algorithm.kisao_id = 'KISAO_0000029'
        task.simulation.algorithm.changes = [
            sedml_data_model.AlgorithmParameterChange(kisao_id='unsupported'),
        ]
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Parameter must'):
            core.exec_sed_task(task, variables)

        task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000488'
        task.simulation.algorithm.changes[0].new_value = ''
        with self.assertRaisesRegex(ValueError, 'not a valid integer'):
            core.exec_sed_task(task, variables)

        task.simulation.algorithm.changes[0].new_value = '10'
        with self.assertRaisesRegex(NotImplementedError, 'is not supported. Initial time must be 0'):
            core.exec_sed_task(task, variables)

        task.simulation.initial_time = 0.
        with self.assertRaisesRegex(ValueError, 'must be at least'):
            core.exec_sed_task(task, variables)

        task.simulation.output_start_time = 10.
        with self.assertRaisesRegex(ValueError, 'must be at least'):
            core.exec_sed_task(task, variables)

        task.simulation.output_end_time = 20.1
        with self.assertRaisesRegex(NotImplementedError, 'must specify an integer'):
            core.exec_sed_task(task, variables)

        task.simulation.output_end_time = 20.
        variables = [
            sedml_data_model.DataGeneratorVariable()
        ]
        with self.assertRaisesRegex(ValueError, 'must define a symbol or target'):
            core.exec_sed_task(task, variables)

        variables = [
            sedml_data_model.DataGeneratorVariable(symbol='x', target='y')
        ]
        with self.assertRaisesRegex(ValueError, 'must define a symbol or target'):
            core.exec_sed_task(task, variables)

        variables = [
            sedml_data_model.DataGeneratorVariable(symbol='unsupported')
        ]
        with self.assertRaisesRegex(NotImplementedError, 'Symbols must be'):
            core.exec_sed_task(task, variables)

        variables = [
            sedml_data_model.DataGeneratorVariable(symbol=sedml_data_model.DataGeneratorVariableSymbol.time),
            sedml_data_model.DataGeneratorVariable(target='--undefined--'),
        ]
        with self.assertRaisesRegex(ValueError, 'Targets must be'):
            core.exec_sed_task(task, variables)

        variables = [
            sedml_data_model.DataGeneratorVariable(id='time', symbol=sedml_data_model.DataGeneratorVariableSymbol.time),
            sedml_data_model.DataGeneratorVariable(id='BE', target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']"),
            sedml_data_model.DataGeneratorVariable(id='BUD', target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BUD"]'),
            sedml_data_model.DataGeneratorVariable(id='Cdc20', target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']"),
        ]

        variable_results = core.exec_sed_task(task, variables)

        self.assertTrue(sorted(variable_results.keys()), sorted([var.id for var in variables]))
        self.assertEqual(variable_results[variables[0].id].shape, (task.simulation.number_of_points + 1,))
        numpy.testing.assert_equal(
            variable_results['time'],
            numpy.linspace(task.simulation.output_start_time, task.simulation.output_end_time, task.simulation.number_of_points + 1),
        )

    def test_exec_sedml_docs_in_combine_archive(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = os.path.join(self.dirname, 'out')
        core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir, report_formats=[report_data_model.ReportFormat.HDF5])

        self._assert_combine_archive_outputs(doc, out_dir)

    def _build_combine_archive(self):
        doc = self._build_sed_doc()

        archive_dirname = os.path.join(self.dirname, 'archive')
        os.mkdir(archive_dirname)

        model_filename = os.path.join(archive_dirname, 'model_1.xml')
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml'),
            model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim_1.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        updated = datetime.datetime(2020, 1, 2, 1, 2, 3, tzinfo=dateutil.tz.tzutc())
        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model_1.xml', combine_data_model.CombineArchiveContentFormat.SBML.value, updated=updated),
                combine_data_model.CombineArchiveContent(
                    'sim_1.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value, updated=updated),
            ],
            updated=updated,
        )
        archive_filename = os.path.join(self.dirname, 'archive.omex')
        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        return (doc, archive_filename)

    def _build_sed_doc(self):
        doc = sedml_data_model.SedDocument()
        doc.models.append(sedml_data_model.Model(
            id='model_1',
            source='model_1.xml',
            language=sedml_data_model.ModelLanguage.SBML.value,
            changes=[],
        ))
        doc.simulations.append(sedml_data_model.UniformTimeCourseSimulation(
            id='sim_1_time_course',
            algorithm=sedml_data_model.Algorithm(
                kisao_id='KISAO_0000029',
                changes=[
                    sedml_data_model.AlgorithmParameterChange(
                        kisao_id='KISAO_0000488',
                        new_value='10',
                    ),
                ],
            ),
            initial_time=0.,
            output_start_time=10.,
            output_end_time=20.,
            number_of_points=20,
        ))
        doc.tasks.append(sedml_data_model.Task(
            id='task_1',
            model=doc.models[0],
            simulation=doc.simulations[0],
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_time',
            variables=[
                sedml_data_model.DataGeneratorVariable(
                    id='var_time',
                    symbol=sedml_data_model.DataGeneratorVariableSymbol.time,
                    task=doc.tasks[0],
                    model=doc.models[0],
                ),
            ],
            math='var_time',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_BE',
            variables=[
                sedml_data_model.DataGeneratorVariable(
                    id='var_BE',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']",
                    task=doc.tasks[0],
                    model=doc.models[0],
                ),
            ],
            math='var_BE',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_BUD',
            variables=[
                sedml_data_model.DataGeneratorVariable(
                    id='var_BUD',
                    target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BUD"]',
                    task=doc.tasks[0],
                    model=doc.models[0],
                ),
            ],
            math='var_BUD',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_Cdc20',
            variables=[
                sedml_data_model.DataGeneratorVariable(
                    id='var_Cdc20',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']",
                    task=doc.tasks[0],
                    model=doc.models[0],
                ),
            ],
            math='var_Cdc20',
        ))
        doc.outputs.append(sedml_data_model.Report(
            id='report_1',
            datasets=[
                sedml_data_model.Dataset(id='data_set_time', label='Time', data_generator=doc.data_generators[0]),
                sedml_data_model.Dataset(id='data_set_BE', label='BE', data_generator=doc.data_generators[1]),
                sedml_data_model.Dataset(id='data_set_BUD', label='BUD', data_generator=doc.data_generators[2]),
                sedml_data_model.Dataset(id='data_set_Cdc20', label='Cdc20', data_generator=doc.data_generators[3]),
            ],
        ))

        append_all_nested_children_to_doc(doc)

        return doc

    def _assert_combine_archive_outputs(self, doc, out_dir):
        self.assertEqual(os.listdir(out_dir), ['reports.h5'])

        report = ReportReader().run(out_dir, 'sim_1.sedml/report_1', format=report_data_model.ReportFormat.HDF5)

        self.assertEqual(sorted(report.index), sorted([d.id for d in doc.outputs[0].datasets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(report.shape, (len(doc.outputs[0].datasets), sim.number_of_points + 1))
        numpy.testing.assert_equal(
            report.loc['data_set_time', :].to_numpy(),
            numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1),
        )

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')

    def test_exec_sedml_docs_in_combine_archive_with_cli(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        env = self._get_combine_archive_exec_env()

        with mock.patch.dict(os.environ, env):
            with __main__.App(argv=['-i', archive_filename, '-o', out_dir]) as app:
                app.run()

        self._assert_combine_archive_outputs(doc, out_dir)

    def _get_combine_archive_exec_env(self):
        return {
            'REPORT_FORMATS': 'HDF5'
        }

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = self._get_combine_archive_exec_env()

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        self._assert_combine_archive_outputs(doc, out_dir)

    def assert_outputs_created(self, dirname):
        name = "BIOMD0000000297"
        self.assertEqual(set(os.listdir(dirname)), set(['ex1', 'ex2']))
        self.assertEqual(set(os.listdir(os.path.join(dirname, 'ex1'))), set([name]))
        self.assertEqual(set(os.listdir(os.path.join(dirname, 'ex2'))), set([name]))
