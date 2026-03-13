# Spatiotemporal Data Fusion for Global Surface Ozone

A Bayesian Maximum Entropy (BME) data fusion framework for global ozone mapping, rigorously incorporating multi-source data from ground-based stations, chemical transport models, and satellite observations.

## Overview
This repository contains the core computational codebase for the BME data fusion framework. It includes massive spatiotemporal data pre-processing pipelines written in **R**, and the main data fusion architecture implemented in **MATLAB**.

The framework is designed to integrate heterogeneous environmental datasets:
* **Hard Data:** High-fidelity ground observations and vertical profiles (e.g., Ozonesondes).
* **Soft Data:** Satellite retrievals accompanied by probabilistic uncertainties (treating satellite data as Gaussian probability density functions rather than exact points).
* **Background Models:** Multi-model composites (e.g., M3Fusion) used to extract the Global Offset (GO).
