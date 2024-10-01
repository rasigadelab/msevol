# msevol
Multi-Scale EVOLution of antibiotic resistance

Developed at the Centre International de Recherche en Infectiologie (CIRI), University of Lyon, Hospices Civils de Lyon.

Concepts and original implementation, 2017-2024, Jean-Philippe Rasigade, jean-philippe.rasigade@chu-lyon.fr. Msevol development supported by the Agence Nationale de la Recherche under grant ANR-20-CE35-0012 and by the University of Lyon under grant LabEx EcoFect.


## Quickstart

### Generating the executable
msevol is tested under Windows with MSVC 2022, although any C++11 compliant compiler should be OK with some work. With MSVC installed, clone the repository, open the MSVC solution and build the "msevolcli" project. The location of the executable, called msevol.exe, depends on MSVC configuration. 

If it is not automatically done, this executbale must then be copy in the "msevolr/bin" folder in order the simulator to be use with R.


### Using msevol with R
A set of companion R scripts was developed to configure an MSEvol model, run simulations on this model, and analyze the results. Please refer to the tutorial 'doc/Rmse1_tutorial.docx' and the example scripts provided to learn how to use MSEvol with R.

The scripts have been tested on Windows using R version 4.1.2, including the following packages: data.table (1.14.2), tidyr (1.1.4), ggplot2 (3.4.4), visNetwork (2.1.0), doParallel (1.0.17), foreach (1.5.2).



## Licence
The msevol code is licensed  under the GNU Affero General Public License v3.0. Please see LICENSE.txt for details.
