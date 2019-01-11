====================================================================================================================
MRST's Dual-Porosity Module
====================================================================================================================
This module was developed by Rafael March and Victoria Spooner, of the Carbonates Reservoir Groups of
Heriot-Watt University (https://carbonates.hw.ac.uk/). Any inquiries, sugestions or feedback should be 
sent to rmc1@hw.ac.uk, ves1@hw.ac.uk, s.geiger@hw.ac.uk and/or f.doster@hw.ac.uk.

We strongly suggest to look at the examples/1ph/single_phase_depletion.m to understand how to initialize a dual-porosity
model.

Description of folders and files:

- models:
--- DualPorosityReservoirModel: Base class for all reservoirs that show dual-porosity behaviour
--- ThreePhaseBlackOilDPModel:  Three-phase black-oil model. This serves as a base class for the two phase model. 
							    However, the equations for the 3ph case are still not implemented and will come 
								in a future release. *****Do not instantiate this model directly, as it will give 
								an error because the equations function does not exist.*****
--- TwoPhaseOilWaterDPModel: Two-phase compressible model. This is ready to be used and the basic usage is shown
							 in an example.							
--- SequentialPressureTransportDPModel: A sequential model for the two-phase compressible model. This is a
										dual-porosity implementation of the SequentialPressureTransportModel in the
										blackoil-sequential module.
--- TransportOilWaterDPModel: A model that represents the transport equation in the sequential model.
--- PressureOilWaterWaterDPModel: A model that represents the pressure equation in the sequential model.

- models\equations:
--- equationsOilWaterDP: Equation for the TwoPhaseOilWaterDPModel.
--- transportEquationOilWaterDP: Transport equation for the TransportOilWaterDPModel.
--- pressureEquationOilWaterDP: Pressure equation for the PressureOilWaterWaterDPModel.

- transfer_models:
--- TransferFunction: this is a base class for all the transfer models. All the transfer models should extend this class.
					  The other files ending with "...TransferFunction" are special implementations of transfer functions
					  available in the literature. The most traditional one is KazemiOilWaterGasTransferFunction (see
					  references).
					  
- shape_factors:
--- ShapeFactor: this is a base class for all the shape factors. All the shape factors should extend this class.
					  The other files ending with "...ShapeFactor" are special implementations of shape factors
					  available in the literature. The most traditional one is KazemiShapeFactor (see
					  references).					  
					  
- utils:
--- getSequentialDPModelFromFI: useful function to transform a fully coupled model (TwoPhaseOilWaterDPModel) into
								a sequential one (SequentialPressureTransportDPModel).
								
- examples:
--- single_phase_depletion: here we show that for a single-phase 0D model the dual-porosity module provides
							the exact solution.
---	two_phase_water_injection: a two-phase water injection example.
--- three_phase_water_alternating_gas: a three-phase WAG example.

References:
[1] Kazemi et al. Numerical Simulation of Water-Oil Flow in Naturally Fractured Reservoirs, SPE Journal, 1976
[2] Quandalle and Sabathier. Typical Features of a Multipurpose Reservoir Simulator, SPE Journal, 1989
[3] Lu et al. General Transfer Functions for Multiphase Flow in Fractured Reservoirs, SPE Journal, 2008

