# Introduction
Measures the shape properties of the object boundaries and displays the curvature.

Original author: Dr. Meghan Driscoll
Modified and compacted/concised a complicated codebase by: Preetham Manjunatha

Thanks to Dr. Meghan Driscoll who kindly shared her code for academic purpose.
If you use this code for visualization and other academic/research/any purposes. 

Please cite:

Reference:
Driscoll MK, McCann C, Kopace R, Homan T, Fourkas JT, Parent C, et al. (2012) 
Cell Shape Dynamics: From Waves to Migration. 
PLoS Comput Biol 8(3): e1002392. 
https://doi.org/10.1371/journal.pcbi.1002392

Important note: This code uses parfor to speed up the things. If you do not have 
the Matlab parallel computing toolbox. Please make 'parfor' as 'for' in this
function.

%%%%%%%% This code is way too slow! (curvature should not be in a for loop) %%%%%%%%%
If I have time I will try to improve this. If anyone improves it, please
share the modified version code with me.


# Quick Pipeline Visualization
## Example: Curvature visualization
| Blob |
| ------------- |
| <img width="1396" alt="blobviz" src="https://user-images.githubusercontent.com/28588878/127915873-2641a1c1-01e2-45a5-ad80-2ad0f57cce0a.png"> |

| T-section |
| ------------- |
| <img width="1399" alt="Tsectionviz" src="https://user-images.githubusercontent.com/28588878/127915909-23db293b-5967-4fea-a301-c29999631615.png"> |

# Requirements
MATLAB <br />
MATLAB - Parallel Computing Toolbox

# MATLAB Central


