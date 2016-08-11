from setuptools import setup, find_packages

setup(
  name="gaia-spatialstats-plugin",
  version="0.0.1",
  description="""Gaia Spatial Statistics plugin""",
  author="Matt Bertrand",
  install_requires=["gaia>=0.0.0"],
  packages=find_packages(),
  include_package_data=True,
  entry_points={
    'gaia.plugins': [
            "gaia_spatialstats.inputs = gaia_spatialstats.inputs",
            "gaia_spatialstats.processes = gaia_spatialstats.processes"
        ]
  }
)
