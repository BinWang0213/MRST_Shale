OpenShale: A Open-source Shale Gas Simulator
==============================================================================================
Bin Wang (binwang.0213@gmail.com), Femi olorode (folorode@lsu.edu)

Craft & Hawkins Department of Petroleum Engineering, Louisiana State University, US

<p align="center">
  <img src = "https://github.com/UnconvRS/OpenShale/blob/master/doc/demo.png" height="400">
</p>

`OpenShale` is a light-weight open source library for simulating flow in naturally-fractured shale-gas reservoirs. 

Currently available features include:

* Explicit fracture modeling (EFM)
* 2D Embedded discrete fracture modeling (EDFM)
* Fracture mesh generation with log-LGR
* Random 2D/3D discrete natural fracture generation (ADFNE, Alghalandis, 2017) 
* Arbitrary user defined sorption&transport function
* Fully implicit compositional/black-oil solver with automatic differentiation (Krogstad et al, 2015)
* Verified solver against commercial simulator (CMG)

TODO List:
* 3D EDFM verification
* Improved EDFM such as, projected EDFM (Tene et al, 2017) or compartmental EDFM (Chai, 2018)
* More verification of compositional simulation 

## Install & Usage

`Anaconda 5.3` (https://www.anaconda.com/download/) and `Shapely 1.5` (https://anaconda.org/scitools/shapely) is required. The library works on both Windows and Linux.

After downloading and unzipping the current <a href="https://github.com/BinWang0213/PyGRDECL/archive/master.zip">repository</a>, navigate to the library directory and run a simple example contained in `Example_GettingStart.ipynb`:

## Example List
| Example | FileName  | Comments |
|---|---|---|
| Case 1a/1b  | Case1_MRST_EDFM/LGR_benchmark.m  | Single vertical HF with EFM and EDFM |
| Case 1c  | Case1_MRST_EDFM/LGR_NF_benchmark.m  | Single vertical HF and 3 horizontal NFs with EFM and EDFM |
| Case 2  | Case2_MRST_LGR_benchmark.m  | 5 vertical HFs verification case (Jiang et al, 2015) |
| Case 3  | Case3_HistoryMatching_Barnett.m  | 28 vertical HFs history matching case of Barnett Shale (Cao et al, 2016) |
| Case 4  | Case4_NaturalFrac_Geomechanics_Barnett.m  | Same with Case 3 but with non-planar HFs, random generated NFs and geomechanics effect |

## Getting Start

```matlab
from GRDECL2VTK import * 

#Read GRDECL File
Model=GeologyModel(filename='./ExampleData/dome.grdecl')

#Convert ECLIPSE grdecl format into VTK
Model.GRDECL2VTK()

#Decompose the model into sub-volumes in terms of faults automatically (this function requires shapely library)
Model.decomposeModel()

#Output to VTK format
Model.Write2VTU()
```

## Reference

* Alghalandis, Y.F., 2017. ADFNE: Open source software for discrete fracture network engineering, two and three dimensional applications. Computers & Geosciences, 102, pp.1-11.
* Cao, P., Liu, J. and Leong, Y.K., 2016. A fully coupled multiscale shale deformation-gas transport model for the evaluation of shale gas extraction. Fuel, 178, pp.103-117.
* Chai, Z., 2018. An Efficient Method for Fractured Shale Reservoir Simulation & History Matching: The cEDFM Approach (Doctoral dissertation).
* Ţene, M., Bosma, S. B., Al Kobaisi, M. S., & Hajibeygi, H. (2017). Projection-based embedded discrete fracture model (pEDFM). Advances in Water Resources, 105, 205-216.
* Jiang, J. and Younis, R.M., 2015. Numerical study of complex fracture geometries for unconventional gas reservoirs using a discrete fracture-matrix model. Journal of Natural Gas Science and Engineering, 26, pp.1174-1186.
* Krogstad, S., Lie, K.A., Møyner, O., Nilsen, H.M., Raynaud, X. and Skaflestad, B., 2015, February. MRST-AD–an open-source framework for rapid prototyping and evaluation of reservoir simulation problems. In SPE reservoir simulation symposium. Society of Petroleum Engineers.

## License

This code is compatible with the Opensource Matlab Reservoir Simulation Toolbox (MRST) provided by SINTEF ICT (http://www.sintef.no/projectweb/mrst/). It requires MRST to be added to the matlab path by running "startup.m".
