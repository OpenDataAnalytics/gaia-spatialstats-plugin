#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#  Copyright Kitware Inc. and Epidemico Inc.
#
#  Licensed under the Apache License, Version 2.0 ( the "License" );
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
###############################################################################
import json
import os
import unittest

import pysal

from gaia import formats
from gaia.parser import deserialize
from gaia.core import config

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestGaiaSpatialStatsViaParser(unittest.TestCase):
    """Tests for the Gaia SpatialStats plugin via Parser"""

    def test_process_cluster(self):
        """Test Cluster Process"""
        with open(os.path.join(testfile_path,
                               'cluster_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = json.loads(process.output.read(format=formats.JSON))
            with open(os.path.join(
                    testfile_path,
                    'cluster_process_results.json')) as gj:
                expected_json = json.load(gj)
            self.assertIn('features', output)
            self.assertEquals(len(expected_json['features']),
                              len(output['features']))
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()

    def test_process_autocorrelation(self):
        """Test Autocorrelation Process"""
        with open(os.path.join(testfile_path,
                               'autocorrelation_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = process.output.read(format=formats.JSON)
            with open(os.path.join(
                    testfile_path,
                    'autocorrelation_process_results.json')) as exp:
                expected_json = json.load(exp)
            self.assertIn('I', output)
            self.assertEquals(expected_json['I'],
                              output['I'])
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()

    def test_process_gamma(self):
        """Test GammaProcess Process"""
        with open(os.path.join(testfile_path,
                               'gamma_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = process.output.read(format=formats.JSON)
            with open(os.path.join(
                    testfile_path,
                    'gamma_process_results.json')) as exp:
                expected_json = json.load(exp)
            self.assertIn('g', output)
            self.assertEquals(expected_json['g'],
                              output['g'])
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()

    def test_process_gearyc(self):
        """Test GearyCProcess Process"""
        with open(os.path.join(testfile_path,
                               'gearyc_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = process.output.read(format=formats.JSON)
            with open(os.path.join(
                    testfile_path,
                    'gearyc_process_results.json')) as exp:
                expected_json = json.load(exp)
            self.assertIn('C', output)
            self.assertEquals(expected_json['C'],
                              output['C'])
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()

    def test_process_weight(self):
        """Test Weight Process"""
        with open(os.path.join(testfile_path,
                               'weight_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = process.output.read(format=formats.WEIGHT)
            exp = pysal.open(os.path.join(testfile_path,
                                          'weight_process_result.gal'), 'r')
            expected_w = exp.read()
            exp.close()
            self.assertEquals(expected_w.n,
                              output.n)
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()

    def test_process_classifier(self):
        """Test ClassifierProcess Process"""
        with open(os.path.join(testfile_path,
                               'classifier_process.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = process.output.read(format=formats.JSON)
            with open(os.path.join(
                    testfile_path,
                    'classifier_process_results.json')) as exp:
                expected_json = json.load(exp)
            self.assertIn('classifier', output)
            self.assertEquals(expected_json['classifier'],
                              output['classifier'])
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()
