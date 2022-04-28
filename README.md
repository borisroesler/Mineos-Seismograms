# Mineos-Seismograms
A Python wrapper for Mineos to generate synthetic seismograms

Mineos, available from https://geodynamics.org, computes synthetic seismograms in a spherically symmetric non-rotating Earth by summing normal modes. Attenuation, gravity and transversal anisotropy effects may be optionally taken into account. Mineos is executed through a command line using appropriate parameters in a parameter file. This Python wrapper generates the necessary input files for Mineos and provides the input parameters for the generation of synthetic seismograms in acceleration. The event file is any moment tensor catalog in the same format as obtained with the GCMT project. Four different Earth models are provided to generate synthetic seismograms for.

The binary files included in this project are the latest files for MacOS. The wrapper can be used different binary files, which can be obtained from https://geodynamics.org.
