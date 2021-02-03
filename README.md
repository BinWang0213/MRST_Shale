MRST-Shale: A Open-source Shale Gas Simulator
==============================================================================================
Bin Wang (binwang.0213@gmail.com)

Craft & Hawkins Department of Petroleum Engineering, Louisiana State University, US

<p align="center">
  <img src = "https://github.com/BinWang0213/MRST_Shale/blob/master/doc/demo.png" height="400">
</p>

`OpenShale` is a light-weight open-source library for simulating flow in naturally-fractured shale-gas reservoirs based on the Matlab Reservoir Simulation Toolbox (MRST) provided by SINTEF ICT (http://www.sintef.no/projectweb/mrst/). It requires MRST to be added to the matlab path by running "startup.m".

Currently available features include:

* Explicit fracture modeling (EFM)
* 2D Embedded discrete fracture modeling (EDFM)
* Fracture mesh generation with log-LGR
* Random 2D/3D discrete natural fracture generation (ADFNE, Alghalandis, 2017) 
* Arbitrary user defined sorption&transport function
* Fully implicit compositional/black-oil solver with automatic differentiation (Krogstad et al, 2015)
* Verified solver against commercial simulator (CMG)

Features available in our new 3D code (https://github.com/UnconvRS/shale/tree/MRST):
* Improved pEDFM ([SPE-201243-PA](https://doi.org/10.2118/201243-PA))
* 3D verification
* Compositional simulation 

This documentation of this code is presented in an upcoming book chapter [Advanced Modelling with the MATLAB Reservoir Simulation Toolbox](https://www.cambridge.org/core/books/advanced-modelling-with-the-matlab-reservoir-simulation-toolbox/7AC2425C73F6F729DB88DB1A504FA1E7)  

## Install & Usage

The code is tested and developed based on `MRST 2018a`. New version of MRST may cause problmes.

### Install ###
1. Install [MRST2018a](http://www.mrst.no). 
2. Add the MRST_Shale folder to the "modules" folder of MRST.
3. Update the "mrst-2018a/modules/hfm" and "mrst-2018a/modules/compositional"  from  "MRST_Shale/3rdParty"
4. Update the "mrst-2018a/utils/ADI.m"  from  "MRST_Shale/3rdParty"

Once MRST and our module are installed and our module can be used like any other MRST module. 

### Getting start ###

To run examples in MRST_Shale:
1. Before any script that relies on the repository is run, MRST must be started. This is done by running the file startup.m which is loacted inside of your MRST directory.
2. Navigate to  "mrst-2018a/modules/" and add the "MRST_Shale" folder to Path by  "Add to Path" - > "Selected Folders and Subfolderes"
3. Run any example script in the "MRST_Shale/examples"


this will make the contents of the geochemistry directory available in the workspace.

## Example List
| Example | FileName  | Comments |
|---|---|---|
| Case 1a/1b  | Case1_MRST_EDFM/LGR_benchmark.m  | Single vertical HF with EFM and EDFM |
| Case 1c  | Case1_MRST_EDFM/LGR_NF_benchmark.m  | Single vertical HF and 3 horizontal NFs with EFM and EDFM |
| Case 2  | Case2_MRST_LGR_benchmark.m  | 5 vertical HFs verification case (Jiang et al, 2015) |
| Case 3  | Case3_HistoryMatching_Barnett.m  | 28 vertical HFs history matching case of Barnett Shale (Cao et al, 2016) |
| Case 4  | Case4_NaturalFrac_Geomechanics_Barnett.m  | Same with Case 3 but with non-planar HFs, random generated NFs and geomechanics effect |

## License
Please cite this project when you using it in your project

Wang, B. MRST-Shale: An Open-Source Framework for Generic Numerical Modeling of Unconventional Shale and Tight Gas Reservoirs. Preprints 2020, 2020010080 (doi: 10.20944/preprints202001.0080.v1).

## Reference

* Alghalandis, Y.F., 2017. ADFNE: Open source software for discrete fracture network engineering, two and three dimensional applications. Computers & Geosciences, 102, pp.1-11.
* Cao, P., Liu, J. and Leong, Y.K., 2016. A fully coupled multiscale shale deformation-gas transport model for the evaluation of shale gas extraction. Fuel, 178, pp.103-117.
* Chai, Z., 2018. An Efficient Method for Fractured Shale Reservoir Simulation & History Matching: The cEDFM Approach (Doctoral dissertation).
* Jiang, J. and Younis, R.M., 2015. Numerical study of complex fracture geometries for unconventional gas reservoirs using a discrete fracture-matrix model. Journal of Natural Gas Science and Engineering, 26, pp.1174-1186.
* Krogstad, S., Lie, K.A., Møyner, O., Nilsen, H.M., Raynaud, X. and Skaflestad, B., 2015, February. MRST-AD–an open-source framework for rapid prototyping and evaluation of reservoir simulation problems. In SPE reservoir simulation symposium. Society of Petroleum Engineers.
