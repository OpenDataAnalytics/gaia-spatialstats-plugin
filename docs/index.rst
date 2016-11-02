Gaia SpatialStats Plugin
================================

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that
provides spatial statistics processes using PySAL, including:

.. Geary's C Process for spatial autocorrelation
.. Gamma Process for spatial autocorrelation
.. Classifier Process for choropleth mapping
.. Local Moran's I Calculation for identifying spatial cluster
.. Global Moran's I Calculation for aptial autocorrelation
.. Spatial weight calculation by contiguity, knnW, distanceBandW, or kernel

Processes
-----------------

* Vector processes
   * Spatial Weight
      * Calculates spatial weight and returns a spatial weight object. Weight type available includes: contiguity, knnW, distanceBandW, kernel
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Weight-Process>`__
   * Geary's C Process
      * Calculates Geary's C statistic for spatial autocorrelation for input data.
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Geary's-C-Process>`__
   * Gamma Process
      * Calculates Gamma Index for spatial autocorrelation for the input data. Uses contiguity weight (queen) and cross product similarity function by default.
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Gamma-Process>`__
   * Classifier Process
      * Classify for choropleth mapping. Default map classifier: Natural_breaks.
Other options: Equal_Interval, Fisher_Jenks, Fisher_Jenks_Sampled, HeadTail_Breaks, Jenks_Caspall, Jenks_Caspall_Forced, Jenks_Caspall_Sampled, Max_P_Classifier, Maximum_Breaks, Natural_Breaks, Quantiles, Percentiles, Std_Mean, User_Defined
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Classifier-Process>`__
   * Cluster Process
      * Calculates the Local Moran's I (Local Indicators of Spatial Association, or LISAs) to identify clusters in data. Returns original vector layer with associated Moran's I statistics
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Cluster-Process>`__
   * Autocorrelation Process
      * Calculates the Moran's I global autcorrelation for the input data. Uses contiguity weight (queen) by default. Returns Moran's I attributes as json.
      * `Example <examples/gaia_spatialstats_plugin.ipynb#Autocorrelation-Process>`__

Installation
-----------------

- git clone https://github.com/OpenDataAnalytics/gaia-spatialstats-plugin.git
- cd gaia-spatialstats-plugin
- pip install -e .
- pip install -r requirements


Testing
-----------------

- python -m unittest discover


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   gaia_spatialstats


.. _Gaia: http://www.github.com/opendataanalytics/gaia
.. _Kitware: http://www.kitware.com
.. _Epidemico: http://epidemico.com
