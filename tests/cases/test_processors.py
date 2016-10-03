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
from gaia.geo.geo_inputs import VectorFileIO
from gaia_spatialstats.processes import ClusterProcess, WeightProcess, \
    AutocorrelationProcess, GearyCProcess, GammaProcess, ClassifierProcess

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestGaiaSpatialStatsProcessors(unittest.TestCase):

    def test_cluster(self):
        """
        Test ClusterProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = ClusterProcess('num_hospitals', inputs=[vector_io])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'cluster_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = json.loads(process.output.read(format=formats.JSON))
            self.assertEquals(len(expected_json['features']),
                              len(actual_json['features']))
        finally:
            if process:
                process.purge()

    def test_autocorrelation(self):
        """
        Test AutocorrelationProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = AutocorrelationProcess('num_hospitals',
                                         inputs=[vector_io])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'autocorrelation_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = process.output.read(format=formats.JSON)
            self.assertEquals(expected_json['I'],
                              actual_json['I'])
        finally:
            if process:
                process.purge()

    def test_weight(self):
        """
        Test WeightProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = WeightProcess('knnW', inputs=[vector_io])
        try:
            process.compute()
            exp = pysal.open(os.path.join(testfile_path,
                                          'weight_process_result.gal'), 'r')
            expected_w = exp.read()
            exp.close()
            actual = process.output.read(format=formats.WEIGHT)
            self.assertEquals(expected_w.n,
                              actual.n)
        finally:
            if process:
                process.purge()

    def test_gearyc(self):
        """
        Test GearyCProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = GearyCProcess('num_hospitals',
                                inputs=[vector_io])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'gearyc_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = process.output.read(format=formats.JSON)
            self.assertEquals(expected_json['C'],
                              actual_json['C'])
        finally:
            if process:
                process.purge()

    def test_gamma(self):
        """
        Test GammaProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = GammaProcess('num_hospitals',
                               inputs=[vector_io])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'gamma_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = process.output.read(format=formats.JSON)
            self.assertEquals(expected_json['g'],
                              actual_json['g'])
        finally:
            if process:
                process.purge()

    def test_classifier(self):
        """
        Test ClassifierProcess for vector inputs
        """
        vector_io = VectorFileIO(
            name='input', uri=os.path.join(testfile_path,
                                           'baghdad_hospitals.geojson'))
        process = ClassifierProcess('num_hospitals', inputs=[vector_io])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'classifier_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = process.output.read(format=formats.JSON)
            self.assertEquals(expected_json['classifier']['classifier_name'],
                              actual_json['classifier']['classifier_name'])
            self.assertEquals(expected_json['classifier']['gadf'],
                              actual_json['classifier']['gadf'])
        finally:
            if process:
                process.purge()
