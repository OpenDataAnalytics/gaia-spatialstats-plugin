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
import pysal
from gaia import formats, types
from gaia.inputs import FileIO, UnsupportedFormatException


class WeightFileIO(FileIO):
    """Read vector and write weight file data (such as .gal)"""

    #: Data type (vector or raster)
    type = types.WEIGHT

    #: acceptable data format extensions
    format = formats.WEIGHT

    #: Default output format
    default_output = formats.WEIGHT

    def read(self, format=None):
        if self.ext not in formats.WEIGHT:
            raise UnsupportedFormatException(
                "Only the following weight formats are supported: {}".format(
                    ','.join(formats.WEIGHT)
                )
            )
        if self.data is None:
            weightfile = pysal.open(self.uri, 'r')
            self.data = weightfile.read()
            weightfile.close()
        return self.data

    def write(self, filename=None, as_type='gal'):
        """
        Write data (assumed pysal weight object) to gal binary weight files
        :param filename: Base filename
        :param as_type: gal
        :return: location of file
        """
        if not filename:
            filename = self.uri
        self.create_output_dir(filename)
        if as_type == 'gal':
            gal = pysal.open(filename, 'w')
            gal.write(self.data)
            gal.close()
        else:
            raise NotImplementedError('{} not a valid type'.format(as_type))
        return self.uri

PLUGIN_CLASS_EXPORTS = [
    WeightFileIO
]
