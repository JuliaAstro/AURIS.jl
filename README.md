# AURIS

Astronomical Utilities for Radio Interferometry and Simulation 

[![Build Status](https://github.com/barrettp/AURIS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/barrettp/AURIS.jl/actions/workflows/CI.yml?query=branch%3Amain)

AURIS.jl is designed for high performance processing of radio interferometry data. Some performance features are:
* using half-precision (16-bit) floats for visibility data, which reduces memory and data throughput with little to no loss of precision,
* memory-mapped files, which allows large files to be processed at one time,
* high-performance optimization algorithms that include built-in auto-differentation for improved accuracy and performance.

The following two images compare the CASA calibration pipeline (left) to the AURIS prototype calibration pipeline (right). The CASA pipeline took 140 minutes to complete. The AURIS pipeline took 1.5 minutes, nearly two orders of magnitude faster.

<img width="344" height="300" alt="XCSJ1040+3957_CASA image" src="https://github.com/user-attachments/assets/36847ed5-dece-4e43-adba-c8f6de8b8361" />

<img width="344" height="300" alt="XCSJ1040+3957_AURIS image" src="https://github.com/user-attachments/assets/575f6d10-b173-486c-b6d1-a61bb703a66f" />
