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
import logging
import numpy as np
import pandas as pd
import json
import pysal

import gaia.formats as formats
from gaia.inputs import JsonFileIO
from gaia.geo import GaiaProcess
from gaia.geo.geo_inputs import VectorFileIO
import gaia_spatialstats.pysal_weights as wt
from gaia_spatialstats.inputs import WeightFileIO

logger = logging.getLogger('gaia')


class ClusterProcess(GaiaProcess):
    """
    Local Moran's I Calculation (Local Indicators of Spatial Association,
    or LISAs) to identifying clusters in data.
    https://pysal.readthedocs.org/en/latest/users/tutorials/autocorrelation.html#local-moran-s-i
    http://pysal.readthedocs.org/en/latest/library/esda/moran.html

    Returns original vector layer with associated Moran's I statistics,
    including:

    lm_Is: float, Moran's I
    lm_q: float, quadrat location
    lm_p_sims: float, p-value based on permutations (low p-value means observed
    Is differ from expected Is significantly)
    lm_sig: boolean, True if p_sims is below 0.05
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('var_col')
    optional_args = ('adjust_by_col')
    default_output = formats.JSON
    adjust_by_col = None

    def __init__(self, var_col, **kwargs):
        self.var_col = var_col
        super(ClusterProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())
        first_df = self.inputs[0].read()
        col = self.var_col
        adjust_by_col = self.adjust_by_col

        # filter out null fields or else weight functions won't work
        keep = first_df[col].notnull()
        filtered_df = first_df[keep].reset_index()

        # get Local Moran's I
        f = np.array(filtered_df[col])
        w = wt.gpd_contiguity(filtered_df)
        if adjust_by_col:
            adjust_by = np.array(filtered_df[adjust_by_col])
            lm = pysal.esda.moran.Moran_Local_Rate(e=f, b=adjust_by, w=w,
                                                   permutations=9999)
        else:
            lm = pysal.Moran_Local(y=f, w=w, permutations=9999)

        sig = lm.p_sim < 0.05
        filtered_df['lm_sig'] = sig
        filtered_df['lm_p_sim'] = lm.p_sim
        filtered_df['lm_q'] = lm.q
        filtered_df['lm_Is'] = lm.Is

        self.output.data = filtered_df
        self.output.write()
        logger.debug(self.output)


class AutocorrelationProcess(GaiaProcess):
    """
    Calculate Moran's I global autocorrelation for the input data.
    Default number of permutations = 999
    Uses contiguity weight (queen) by default.
    https://pysal.readthedocs.org/en/latest/users/tutorials/autocorrelation.html#moran-s-i
    http://pysal.readthedocs.org/en/latest/library/esda/moran.html

    Returns the following Moran's I attributes as json:
    I: float, value of Moran's I
    EI: float, expected value of I under normality assumption
    p_norm: float, p-value of I under normality assumption
    EI_sim: float, average value of I from permutations
    p_sim: array, p-value based on permutations (one-tailed)
    z_sim: float, standardized I based on permutations
    p_z_sim: float, p-value based on standard normal approximation
    from permutations
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('var_col')
    optional_args = ('adjust_by_col', 'permutations')
    default_output = formats.JSON
    adjust_by_col = None
    permutations = None

    def __init__(self, var_col, **kwargs):
        self.var_col = var_col
        super(AutocorrelationProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = JsonFileIO(name='result',
                                     uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())
        for input in self.inputs:
            if input.name == 'input':
                first_df = input.read()
        col = self.var_col
        adjust_by_col = self.adjust_by_col
        permutations = self.permutations
        if not permutations:
            permutations = 999

        # filter out null fields
        keep = first_df[col].notnull()
        filtered_df = first_df[keep]

        # get Global Moran's I
        f = np.array(filtered_df[col])
        w = wt.gpd_contiguity(filtered_df)
        if adjust_by_col:
            adjust_by = np.array(filtered_df[adjust_by_col])
            mi = pysal.esda.moran.Moran_Rate(e=f, b=adjust_by, w=w,
                                             permutations=permutations)
        else:
            mi = pysal.Moran(y=f, w=w, permutations=permutations)

        keep = ['I', 'EI', 'p_norm', 'EI_sim', 'p_sim', 'z_sim', 'p_z_sim']
        mi_dict = {k: getattr(mi, k) for k in keep}

        self.output.data = mi_dict
        self.output.write()
        logger.debug(self.output)


class WeightProcess(GaiaProcess):
    """
    Calculate spatial weight.
    weight_type available includes: contiguity, knnW, distanceBandW, kernel
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('weight_type')
    default_output = formats.WEIGHT

    def __init__(self, weight_type, **kwargs):
        self.weight_type = weight_type
        super(WeightProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = WeightFileIO(name='result',
                                       uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())
        for input in self.inputs:
            if input.name == 'input':
                first_df = input.read()
        weight_type = self.weight_type
        if weight_type == 'contiguity':
            w = wt.gpd_contiguity(first_df)
        elif weight_type == 'knnW':
            w = wt.gpd_knnW(first_df)
        elif weight_type == 'distanceBandW':
            w = wt.gpd_distanceBandW(first_df)
        elif weight_type == 'kernel':
            w = wt.gpd_kernel(first_df)
        # TODO: add params related to dif weight types
        else:
            print(u'weight type {0} not available'.format(weight_type))
        self.output.data = w
        self.output.write()
        logger.debug(self.output)


class GearyCProcess(GaiaProcess):
    """
    Calculate Geary's C statistic for spatial autocorrelation for input data.
    Default number of permutations = 999
    Uses contiguity weight (queen) and binary transformation by default.
    http://pysal.readthedocs.io/en/v1.11.0/library/esda/geary.html

    Returns the following Geary's C attributes as json:
    C: float, value of Geary's C
    EC: float, expected value of C
    p_norm: float, p-value of C under normality assumption
    EC_sim: float, average value of C from permutations (if permutations!=0)
    p_sim: array, p-value based on permutations (one-tailed)
    z_sim: float, standardized C based on permutations
    p_z_sim: float, p-value based on standard normal approximation
    from permutations (one-tailed)
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('var_col')
    optional_args = ('transformation', 'permutations')
    default_output = formats.JSON
    transformation = None
    permutations = None

    def __init__(self, var_col, **kwargs):
        self.var_col = var_col
        super(GearyCProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = JsonFileIO(name='result',
                                     uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())
        for input in self.inputs:
            if input.name == 'input':
                first_df = input.read()
        col = self.var_col
        transformation = self.transformation
        permutations = self.permutations
        if not permutations:
            permutations = 999
        if not transformation:
            transformation = 'B'

        # filter out null fields
        keep = first_df[col].notnull()
        filtered_df = first_df[keep]

        # get Geary's C
        f = np.array(filtered_df[col])
        w = wt.gpd_contiguity(filtered_df)

        c = pysal.esda.geary.Geary(y=f, w=w,
                                   transformation=transformation,
                                   permutations=permutations)

        keep = ['C', 'EC', 'p_norm', 'EC_sim', 'p_sim', 'z_sim', 'p_z_sim']
        c_dict = {k: getattr(c, k) for k in keep}
        self.output.data = c_dict
        self.output.write()
        logger.debug(self.output)


class GammaProcess(GaiaProcess):
    """
    Calculate Gamma index for spatial autocorrelation for the input data.
    Default number of permutations = 999
    Uses contiguity weight (queen), cross product similarity function.
    http://pysal.readthedocs.io/en/v1.11.0/library/esda/gamma.html

    Returns the following Gamma index attributes as json:
    gamma: float, value of Gamma index
    p_sim_g: array, p-value based on permutations (one-sided)
    mean_g: float, average of permuted Gamma values
    min_g: float, minimum of permuted Gamma values
    max_g: float, maximum of permuted Gamma values
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('var_col')
    optional_args = ('operation', 'standardize', 'permutations')
    default_output = formats.JSON
    operation = None
    permutations = None
    standardize = None

    def __init__(self, var_col, **kwargs):
        self.var_col = var_col
        super(GammaProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = JsonFileIO(name='result',
                                     uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = VectorFileIO(name='result',
                                       uri=self.get_outpath())
        for input in self.inputs:
            if input.name == 'input':
                first_df = input.read()
        col = self.var_col
        permutations = self.permutations
        operation = self.operation
        standardize = self.standardize
        if not permutations:
            permutations = 999
        if not operation:
            operation = 'c'
        if not standardize:
            standardize = 'no'

        # filter out null fields
        keep = first_df[col].notnull()
        filtered_df = first_df[keep]

        # get Gamma's index
        f = np.array(filtered_df[col])
        w = wt.gpd_contiguity(filtered_df)

        g = pysal.esda.gamma.Gamma(y=f, w=w,
                                   operation=operation,
                                   standardize=standardize,
                                   permutations=permutations)

        keep = ['g', 'p_sim_g', 'mean_g', 'min_g', 'max_g']
        g_dict = {k: getattr(g, k) for k in keep}
        self.output.data = g_dict
        self.output.write()
        logger.debug(self.output)


class ClassifierProcess(GaiaProcess):
    """
    Classify for choropleth mapping.
    default map_classifier: Natural_Breaks. Other options: Equal_Interval,
    Fisher_Jenks, Fisher_Jenks_Sampled, HeadTail_Breaks, Jenks_Caspall,
    Jenks_Caspall_Forced, Jenks_Caspall_Sampled, Max_P_Classifier,
    Maximum_Breaks, Natural_Breaks, Quantiles, Percentiles, Std_Mean,
    User_Defined

    http://pysal.readthedocs.io/en/v1.11.0/library/esda/mapclassify.html

    Returns GeoJson of input data with "class" property for each feature
    as well as a high level "classifier" attribute that includes:
    classifier_name: string, name of map classifier used
    classes: dict, with class label as key and counts of features in
    class in each value
    gadf: float, Goodness of absolute deviation of fit
    """
    required_inputs = (('input', formats.VECTOR),)
    required_args = ('var_col')
    optional_args = ('map_classifier', 'k', 'pick_optimal',
                     'hinge', 'pct', 'initial', 'mindiff', 'multiples', 'bins')
    default_output = formats.JSON
    map_classifier = None
    k = None
    pick_optimal = None

    # optional args for different classifier
    hinge = None
    pct = None
    initial = None
    mindiff = None
    multiples = None
    bins = None

    def __init__(self, var_col, **kwargs):
        self.var_col = var_col
        super(ClassifierProcess, self).__init__(**kwargs)
        if not self.output:
            self.output = JsonFileIO(name='result',
                                     uri=self.get_outpath())

    def compute(self):
        if not self.output:
            self.output = JsonFileIO(name='result',
                                     uri=self.get_outpath())
        for input in self.inputs:
            if input.name == 'input':
                first_df = input.read()
        col = self.var_col
        map_classifier = self.map_classifier
        pick_optimal = self.pick_optimal

        # optional args for different classifier
        opt_kwargs = dict(
            k=self.k,
            hinge=self.hinge,
            pct=self.pct,
            initial=self.initial,
            mindiff=self.mindiff,
            multiples=self.multiples,
            bins=self.bins
        )

        opt_kwargs = {k: v for k, v in opt_kwargs.items() if v}

        if pick_optimal and not map_classifier:
            # pick optimal k-classifier based on k and gadf
            y = np.array(first_df[col])
            ks = pysal.esda.mapclassify.K_classifiers(y, pct=0.8)
            map_classifier = ks.best.name
            opt_kwargs['k'] = ks.best.k

        elif not pick_optimal and map_classifier:
            map_classifier = self.map_classifier
        else:
            map_classifier = 'Natural_Breaks'

        y = np.array(first_df[col])

        # get classifier
        classifier = getattr(pysal.esda.mapclassify, map_classifier)
        result = classifier(y, **opt_kwargs)

        # create labels and bin data
        result_bins = list(result.bins)
        labels = []
        num_bins = len(result_bins) - 1
        for i in range(0, num_bins):
            if i == (num_bins):
                labels.append('> {}'.format(result_bins[i]))
            else:
                labels.append('{} - {}'.format(result_bins[i],
                                               result_bins[i+1]))

        first_df['class'] = pd.cut(y,
                                   result_bins,
                                   labels=labels,
                                   include_lowest=True)

        # format output
        result_counts = getattr(result, 'counts', None)
        result_classes = {label: {'count': result_counts[idx]}
                          for idx, label in enumerate(labels)}
        output_dict = json.loads(first_df.to_json())
        output_dict['classifier'] = dict(
            classifier_name=map_classifier,
            classes=result_classes
        )

        try:
            output_dict['classifier']['gadf'] = result.get_gadf()
        except NotImplementedError:
            output_dict['classifier']['gadf'] = None

        self.output.data = output_dict
        self.output.write()
        logger.debug(self.output)


PLUGIN_CLASS_EXPORTS = [
    ClusterProcess,
    AutocorrelationProcess,
    WeightProcess,
    GearyCProcess,
    GammaProcess,
    ClassifierProcess
]
