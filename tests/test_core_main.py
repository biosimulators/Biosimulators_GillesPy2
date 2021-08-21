""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_gillespy2 import __main__
from biosimulators_gillespy2 import core
from biosimulators_gillespy2.data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.archive.io import ArchiveReader
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.log.data_model import TaskLog
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.simulator.exec import exec_sedml_docs_in_archive_with_containerized_simulator
from biosimulators_utils.simulator.specs import gen_algorithms_from_specs
from biosimulators_utils.sedml import data_model as sedml_data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.sedml.utils import append_all_nested_children_to_doc
from biosimulators_utils.utils.core import are_lists_equal
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from unittest import mock
import enum
import datetime
import dateutil.tz
import numpy
import numpy.testing
import os
import shutil
import tempfile
import unittest


class TestCase(unittest.TestCase):
    DOCKER_IMAGE = 'ghcr.io/biosimulators/biosimulators_gillespy2/gillespy2:latest'
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level2/version4',
    }

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
                    kisao_id='KISAO_0000088',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000211',
                            new_value='1e-10',
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
            sedml_data_model.Variable(
                id='time',
                symbol=sedml_data_model.Symbol.time,
                task=task,
            ),
            sedml_data_model.Variable(
                id='BE',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='Cdh1',
                target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="Cdh1"]',
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='Cdc20',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='kswe_prime',
                target="/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='kswe_prime']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='kswe',
                target="/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='kswe']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='compartment',
                target="/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment[@id='compartment']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
        ]

        variable_results, _ = core.exec_sed_task(task, variables, TaskLog())

        self.assertTrue(sorted(variable_results.keys()), sorted([var.id for var in variables]))
        self.assertEqual(variable_results[variables[0].id].shape, (task.simulation.number_of_points + 1,))
        numpy.testing.assert_almost_equal(
            variable_results['time'],
            numpy.linspace(task.simulation.output_start_time, task.simulation.output_end_time, task.simulation.number_of_points + 1),
        )
        for variable in variables:
            self.assertFalse(numpy.any(numpy.isnan(variable_results[variable.id])))
        numpy.testing.assert_almost_equal(
            variable_results['kswe_prime'],
            numpy.full((task.simulation.number_of_points + 1,), 2.)
        )
        numpy.testing.assert_almost_equal(
            variable_results['compartment'],
            numpy.full((task.simulation.number_of_points + 1,), 1.)
        )

    def test_exec_sed_task_errors(self):
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            task = sedml_data_model.Task()
            task.model = sedml_data_model.Model(id='model')
            task.model.source = os.path.join(self.dirname, 'valid-model.xml')
            with open(task.model.source, 'w') as file:
                file.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
                file.write('<sbml2 xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">')
                file.write('  <model id="model">')
                file.write('  </model>')
                file.write('</sbml2>')
            task.model.language = sedml_data_model.ModelLanguage.SBML
            task.model.changes = []
            task.simulation = sedml_data_model.UniformTimeCourseSimulation(
                id='simulation',
                algorithm=sedml_data_model.Algorithm(kisao_id='KISAO_0000448'),
                initial_time=10.,
                output_start_time=10.,
                output_end_time=20.1,
                number_of_points=10,
            )

            variables = []

            with self.assertRaisesRegex(ValueError, 'could not be imported'):
                core.exec_sed_task(task, variables, TaskLog())
            task.model.source = os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml')

            with self.assertRaisesRegex(AlgorithmCannotBeSubstitutedException, 'No algorithm can be substituted'):
                core.exec_sed_task(task, variables, TaskLog())
            task.simulation.algorithm.kisao_id = 'KISAO_0000029'
            task.simulation.algorithm.changes = [
                sedml_data_model.AlgorithmParameterChange(kisao_id='KISAO_0000531'),
            ]

            with self.assertRaisesRegex(NotImplementedError, 'is not supported. Parameter must'):
                core.exec_sed_task(task, variables, TaskLog())
            task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000488'
            task.simulation.algorithm.changes[0].new_value = 'abc'

            with self.assertRaisesRegex(ValueError, 'not a valid integer'):
                core.exec_sed_task(task, variables, TaskLog())
            task.simulation.algorithm.changes[0].new_value = '10'

            with self.assertRaisesRegex(NotImplementedError, 'is not supported. Initial time must be 0'):
                core.exec_sed_task(task, variables, TaskLog())
            task.simulation.initial_time = 0.

            with self.assertRaisesRegex(NotImplementedError, 'must specify an integer'):
                core.exec_sed_task(task, variables, TaskLog())
            task.simulation.output_end_time = 20.
            variables = [
                sedml_data_model.Variable(id='var_1', symbol='unsupported', task=task)
            ]

            with self.assertRaisesRegex(NotImplementedError, 'Symbols must be'):
                core.exec_sed_task(task, variables, TaskLog())
            variables = [
                sedml_data_model.Variable(id='var_1', symbol=sedml_data_model.Symbol.time, task=task),
                sedml_data_model.Variable(id='var_2', target='/invalid:target', target_namespaces={'invalid': 'invalid'}, task=task),
            ]

            with self.assertRaisesRegex(ValueError, 'XPaths must reference unique objects.'):
                core.exec_sed_task(task, variables, TaskLog())
            variables = [
                sedml_data_model.Variable(id='var_1', symbol=sedml_data_model.Symbol.time, task=task),
                sedml_data_model.Variable(
                    id='BE',
                    target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R1']",
                    target_namespaces=self.NAMESPACES,
                    task=task
                ),
            ]

            with self.assertRaisesRegex(ValueError, 'Targets must have'):
                core.exec_sed_task(task, variables, TaskLog())
            variables = [
                sedml_data_model.Variable(
                    id='time',
                    symbol=sedml_data_model.Symbol.time,
                    task=task
                ),
                sedml_data_model.Variable(
                    id='BE',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']",
                    target_namespaces=self.NAMESPACES,
                    task=task,
                ),
                sedml_data_model.Variable(
                    id='Cdh1',
                    target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="Cdh1"]',
                    target_namespaces=self.NAMESPACES,
                    task=task,
                ),
                sedml_data_model.Variable(
                    id='Cdc20',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']",
                    target_namespaces=self.NAMESPACES,
                    task=task,
                ),
            ]

            task.simulation.algorithm.kisao_id = 'KISAO_0000088'
            task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000211'
            task.simulation.algorithm.changes[0].new_value = '1e-10'
            variable_results, _ = core.exec_sed_task(task, variables, TaskLog())

            self.assertTrue(sorted(variable_results.keys()), sorted([var.id for var in variables]))
            self.assertEqual(variable_results[variables[0].id].shape, (task.simulation.number_of_points + 1,))
            numpy.testing.assert_almost_equal(
                variable_results['time'],
                numpy.linspace(task.simulation.output_start_time, task.simulation.output_end_time, task.simulation.number_of_points + 1),
            )

        # algorithm substitution
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml'),
                language=sedml_data_model.ModelLanguage.SBML.value,
                changes=[],
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000088',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000211',
                            new_value='1e-10',
                        ),
                    ],
                ),
                initial_time=0.,
                output_start_time=10.,
                output_end_time=20.,
                number_of_points=20,
            ),
        )

        variables = []

        task.simulation.algorithm.changes[0].new_value = 'not a number'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaisesRegex(ValueError, 'is not a valid'):
                core.exec_sed_task(task, variables, TaskLog())

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarnsRegex(BioSimulatorsWarning, 'Unsuported value'):
                core.exec_sed_task(task, variables, TaskLog())

        task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000531'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaisesRegex(NotImplementedError, 'is not supported'):
                core.exec_sed_task(task, variables, TaskLog())

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarnsRegex(BioSimulatorsWarning, 'was ignored because it is not supported'):
                core.exec_sed_task(task, variables, TaskLog())

    def test_exec_sedml_docs_in_combine_archive(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = os.path.join(self.dirname, 'out')
        core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                report_formats=[
                                                    report_data_model.ReportFormat.h5,
                                                    report_data_model.ReportFormat.csv,
                                                ],
                                                bundle_outputs=True,
                                                keep_individual_outputs=True)

        self._assert_combine_archive_outputs(doc, out_dir)

    def test_exec_sedml_docs_in_combine_archive_with_all_algorithms(self):
        for alg in gen_algorithms_from_specs(os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')).values():
            alg_props = KISAO_ALGORITHM_MAP[alg.kisao_id]
            alg.changes = []
            for param_kisao_id, param_props in alg_props.parameters.items():

                if param_kisao_id == 'KISAO_0000488':
                    new_value = '1'
                else:
                    new_value = param_props.default

                if isinstance(new_value, enum.Enum):
                    new_value = new_value.value
                if new_value is None:
                    new_value = ''
                else:
                    new_value = str(new_value)

                alg.changes.append(sedml_data_model.AlgorithmParameterChange(
                    kisao_id=param_kisao_id,
                    new_value=new_value,
                ))
            doc, archive_filename = self._build_combine_archive(algorithm=alg)

            variables = []
            for data_gen in doc.data_generators:
                for var in data_gen.variables:
                    variables.append(var)
            doc.tasks[0].model.source = os.path.join(os.path.dirname(__file__),
                                                     'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml')
            results, _ = core.exec_sed_task(doc.tasks[0], variables, TaskLog())
            self.assertEqual(set(results.keys()), set(var.id for var in variables))

            alg.changes = []
            for param_kisao_id, param_props in alg_props.parameters.items():
                new_value = param_props.default
                if isinstance(new_value, enum.Enum):
                    new_value = new_value.value
                if new_value is not None:
                    new_value = str(new_value)
                    alg.changes.append(sedml_data_model.AlgorithmParameterChange(
                        kisao_id=param_kisao_id,
                        new_value=new_value,
                    ))
            doc, archive_filename = self._build_combine_archive(algorithm=alg)

            out_dir = os.path.join(self.dirname, alg.kisao_id)
            core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                    report_formats=[
                                                        report_data_model.ReportFormat.h5,
                                                        report_data_model.ReportFormat.csv,
                                                    ],
                                                    bundle_outputs=True,
                                                    keep_individual_outputs=True)

            self._assert_combine_archive_outputs(doc, out_dir)

    def _build_combine_archive(self, algorithm=None):
        doc = self._build_sed_doc(algorithm=algorithm)

        archive_dirname = os.path.join(self.dirname, 'archive')
        if not os.path.isdir(archive_dirname):
            os.mkdir(archive_dirname)

        model_filename = os.path.join(archive_dirname, 'model_1.xml')
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited', 'ex1', 'BIOMD0000000297.xml'),
            model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim_1.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model_1.xml', combine_data_model.CombineArchiveContentFormat.SBML.value),
                combine_data_model.CombineArchiveContent(
                    'sim_1.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value),
            ],
        )
        archive_filename = os.path.join(self.dirname,
                                        'archive.omex' if algorithm is None else 'archive-{}.omex'.format(algorithm.kisao_id))
        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        return (doc, archive_filename)

    def _build_sed_doc(self, algorithm=None):
        if algorithm is None:
            algorithm = sedml_data_model.Algorithm(
                kisao_id='KISAO_0000088',
                changes=[
                    sedml_data_model.AlgorithmParameterChange(
                        kisao_id='KISAO_0000211',
                        new_value='1e-10',
                    ),
                ],
            )

        doc = sedml_data_model.SedDocument()
        doc.models.append(sedml_data_model.Model(
            id='model_1',
            source='model_1.xml',
            language=sedml_data_model.ModelLanguage.SBML.value,
            changes=[],
        ))
        doc.simulations.append(sedml_data_model.UniformTimeCourseSimulation(
            id='sim_1_time_course',
            algorithm=algorithm,
            initial_time=0.,
            output_start_time=0.1,
            output_end_time=0.2,
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
                sedml_data_model.Variable(
                    id='var_time',
                    symbol=sedml_data_model.Symbol.time,
                    task=doc.tasks[0],
                ),
            ],
            math='var_time',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_BE',
            variables=[
                sedml_data_model.Variable(
                    id='var_BE',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='BE']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_BE',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_Cdh1',
            variables=[
                sedml_data_model.Variable(
                    id='var_Cdh1',
                    target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="Cdh1"]',
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_Cdh1',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_Cdc20',
            variables=[
                sedml_data_model.Variable(
                    id='var_Cdc20',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='Cdc20']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_Cdc20',
        ))
        doc.outputs.append(sedml_data_model.Report(
            id='report_1',
            data_sets=[
                sedml_data_model.DataSet(id='data_set_time', label='Time', data_generator=doc.data_generators[0]),
                sedml_data_model.DataSet(id='data_set_BE', label='BE', data_generator=doc.data_generators[1]),
                sedml_data_model.DataSet(id='data_set_Cdh1', label='Cdh1', data_generator=doc.data_generators[2]),
                sedml_data_model.DataSet(id='data_set_Cdc20', label='Cdc20', data_generator=doc.data_generators[3]),
            ],
        ))

        append_all_nested_children_to_doc(doc)

        return doc

    def _assert_combine_archive_outputs(self, doc, out_dir):
        self.assertEqual(set(['reports.h5', 'reports.zip', 'sim_1.sedml']).difference(set(os.listdir(out_dir))), set())

        report = doc.outputs[0]

        # check HDF report
        report_results = ReportReader().run(report, out_dir, 'sim_1.sedml/report_1', format=report_data_model.ReportFormat.h5)

        self.assertEqual(sorted(report_results.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(len(report_results[report.data_sets[0].id]), sim.number_of_points + 1)
        numpy.testing.assert_almost_equal(
            report_results[report.data_sets[0].id],
            numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1),
        )

        for data_set_result in report_results.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

        # check CSV report
        report_results = ReportReader().run(report, out_dir, 'sim_1.sedml/report_1', format=report_data_model.ReportFormat.csv)

        self.assertEqual(sorted(report_results.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(len(report_results[report.data_sets[0].id]), sim.number_of_points + 1)
        numpy.testing.assert_almost_equal(
            report_results[report.data_sets[0].id],
            numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1),
        )

        for data_set_result in report_results.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

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
            'REPORT_FORMATS': 'h5,csv'
        }

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = self._get_combine_archive_exec_env()

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        self._assert_combine_archive_outputs(doc, out_dir)

    def test_more_complex_archive(self):
        archive_filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000297.edited.omex')
        core.exec_sedml_docs_in_combine_archive(archive_filename, self.dirname,
                                                report_formats=[
                                                    report_data_model.ReportFormat.csv,
                                                    report_data_model.ReportFormat.h5,
                                                ],
                                                plot_formats=[],
                                                bundle_outputs=True,
                                                keep_individual_outputs=True)

        self.assertEqual(set(['reports.zip', 'reports.h5', 'ex1', 'ex2']).difference(set(os.listdir(self.dirname))), set())
        self.assertEqual(set(os.listdir(os.path.join(self.dirname, 'ex1'))), set(['BIOMD0000000297.sedml']))
        self.assertEqual(set(os.listdir(os.path.join(self.dirname, 'ex2'))), set(['BIOMD0000000297.sedml']))
        self.assertEqual(set(os.listdir(os.path.join(self.dirname, 'ex1', 'BIOMD0000000297.sedml'))),
                         set(['two_species.csv', 'three_species.csv']))
        self.assertEqual(set(os.listdir(os.path.join(self.dirname, 'ex2', 'BIOMD0000000297.sedml'))),
                         set(['one_species.csv', 'four_species.csv']))

        archive = ArchiveReader().run(os.path.join(self.dirname, 'reports.zip'))

        self.assertEqual(
            sorted(file.archive_path for file in archive.files),
            sorted([
                'ex1/BIOMD0000000297.sedml/two_species.csv',
                'ex1/BIOMD0000000297.sedml/three_species.csv',
                'ex2/BIOMD0000000297.sedml/one_species.csv',
                'ex2/BIOMD0000000297.sedml/four_species.csv',
            ]),
        )

        report = sedml_data_model.Report(
            data_sets=[
                sedml_data_model.DataSet(id='data_set_time_two_species', label='time'),
                sedml_data_model.DataSet(id='data_set_Cln4', label='Cln4'),
                sedml_data_model.DataSet(id='data_set_Swe13', label='Swe13'),
            ]
        )

        report_results = ReportReader().run(report, self.dirname,
                                            'ex1/BIOMD0000000297.sedml/two_species',
                                            format=report_data_model.ReportFormat.h5)
        self.assertEqual(sorted(report_results.keys()), sorted(['data_set_time_two_species', 'data_set_Cln4', 'data_set_Swe13']))
        numpy.testing.assert_almost_equal(report_results['data_set_time_two_species'], numpy.linspace(0., 1., 10 + 1))
