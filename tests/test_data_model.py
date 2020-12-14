""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_gillespy2 import data_model
import enum
import gillespy2
import json
import os
import unittest


class DataModelTestCase(unittest.TestCase):
    def test_algorithm(self):
        alg = data_model.Algorithm("LSODA", gillespy2.ODESolver, integrator="lsoda", parameters=None)
        self.assertEqual(alg.name, "LSODA")
        self.assertEqual(alg.solver, gillespy2.ODESolver)
        self.assertEqual(alg.solver_args, {'integrator': 'lsoda'})
        self.assertEqual(alg.parameters, {})

        alg = data_model.Algorithm("LSODA", gillespy2.ODESolver, integrator="lsoda", parameters={'param_1': 'value'})
        self.assertEqual(alg.name, "LSODA")
        self.assertEqual(alg.solver, gillespy2.ODESolver)
        self.assertEqual(alg.solver_args, {'integrator': 'lsoda'})
        self.assertEqual(alg.parameters, {'param_1': 'value'})

    def test_algorithm_parameter(self):
        param = data_model.AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12)
        self.assertEqual(param.name, "absolute tolerance")
        self.assertEqual(param.key, "integrator_options.atol")
        self.assertEqual(param.data_type, float)
        self.assertEqual(param.default, 1e-12)

    def test_algorithm_parameter_boolean(self):
        param = data_model.AlgorithmParameter('name', 'integrator_options', bool, True)

        solver_args = {}
        param.set_value(solver_args, 'false')
        self.assertEqual(solver_args['integrator_options'], False)

        solver_args = {}
        param.set_value(solver_args, '1')
        self.assertEqual(solver_args['integrator_options'], True)

        with self.assertRaises(ValueError):
            param.set_value(solver_args, 'f')

    def test_algorithm_parameter_integer(self):
        param = data_model.AlgorithmParameter('name', 'integrator_options.atol', int, 10)
        solver_args = {}
        param.set_value(solver_args, '11')
        self.assertEqual(solver_args['integrator_options']['atol'], 11)

        with self.assertRaises(ValueError):
            param.set_value(solver_args, '1.1')

        with self.assertRaises(ValueError):
            param.set_value(solver_args, '1.')

        with self.assertRaises(ValueError):
            param.set_value(solver_args, 'a')

    def test_algorithm_parameter_float(self):
        param = data_model.AlgorithmParameter("absolute tolerance", 'integrator_options.atol', float, 1e-12)
        solver_args = {}
        param.set_value(solver_args, '1e-14')
        self.assertEqual(solver_args['integrator_options']['atol'], 1e-14)

        with self.assertRaises(ValueError):
            param.set_value(solver_args, 'a')

    def test_algorithm_parameter_enum(self):
        param = data_model.AlgorithmParameter('name', 'integrator_options.atol', data_model.VodeMethod, data_model.VodeMethod.bdf.value)
        solver_args = {}
        param.set_value(solver_args, data_model.VodeMethod.bdf.value)
        self.assertEqual(solver_args['integrator_options']['atol'], data_model.VodeMethod.bdf.name)

        with self.assertRaises(NotImplementedError):
            param.set_value(solver_args, '--invalid--')

    def test_algorithm_parameter_invalid_type(self):
        param = data_model.AlgorithmParameter('name', 'integrator_options.atol', str, 'default')
        solver_args = {}
        with self.assertRaises(NotImplementedError):
            param.set_value(solver_args, 'value')

    def test_data_model_matches_specifications(self):
        with open(os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json'), 'r') as file:
            specs = json.load(file)

        self.assertEqual(
            set(data_model.kisao_algorithm_map.keys()),
            set(alg_specs['kisaoId']['id'] for alg_specs in specs['algorithms']))

        for alg_specs in specs['algorithms']:
            alg_props = data_model.kisao_algorithm_map[alg_specs['kisaoId']['id']]

            self.assertEqual(set(alg_props.parameters.keys()), set(param_specs['kisaoId']['id'] for param_specs in alg_specs['parameters']))

            for param_specs in alg_specs['parameters']:
                param_props = alg_props.parameters[param_specs['kisaoId']['id']]

                param_props_specs_match = (
                    (param_props.data_type is bool and param_specs['type'] == 'boolean')
                    or (param_props.data_type is int and param_specs['type'] == 'integer')
                    or (param_props.data_type is float and param_specs['type'] == 'float')
                    or (
                        issubclass(param_props.data_type, enum.Enum) and param_specs['type'] == 'kisaoId'
                        and set(param_props.data_type.__members__.values()) == set(param_specs['recommendedRange'])
                    )
                )
                self.assertTrue(param_props_specs_match)
