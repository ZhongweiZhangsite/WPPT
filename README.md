# 'WPPT': A solver for Wave- and Particle-like Phonon Transport (WPPT)

## Note that we have only uploaded a draft of softeware to this page. The well designed sourcecode will be available soon.
## For any request, please contact the authors.

## Authors:

- Zhongwei Zhang <zhongweizhang123@gmail.com> / <zhongwei@tongji.edu.cn>, Institute of Industrial Science, The University of Tokyo
- Sebastian Volz <volz@iis.u-tokyo.ac.jp>, LIMMS, CNRS-IIS UMI 2820, The University of Tokyo
- Jie Chen <jie@tongji.edu.cn>, Tongji University
- Masahiro Nomura <nomura@iis.u-tokyo.ac.jp>, Institute of Industrial Science, The University of Tokyo

## The code contains three modules: (in preparing)

- NMA: Time-dependent normal mode velocity or displacement calculation.
- PWave: Phonon Wavelet transform for the coherence time distribution, coherence length distribution and phonon decay calculation. (upcoming soon for the coherence length calculation part)
- WPpt: 1) particlelike lifetimes, coherence corrected lifetimes and coherence time; 2) phonon coherence from spectroscopy or MD spectral energy; 3) thermal conductivity calculation including wavelike and particlelike phonons. 

## The code should have interface to: (in preparing)

- LAMMPS
- GPUMD
- VASP & CP2K & DFTB

## Manual:

- Please refers to the the [online manual](https://zhongweizhang123.wixsite.com/wppt).
  
## We are gradually improving the code and creating more examples for friendly use!

Before that please refers to our publications or preprint about the used theories.

- [1] Z. Zhang, Y. Guo, M. Bescond, J. Chen, M. Nomura, and S. Volz, Generalized Decay Law for Particlelike and Wavelike Thermal Phonons, Phys. Rev. B 103, 184307 (2021). [Download](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.184307).
- [2] Z. Zhang, Y. Guo, M. Bescond, J. Chen, M. Nomura, and S. Volz, Heat conduction theory including phonon coherence, Phys. Rev. Lett. 128, 015901 (2022). [Download](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.015901).
- [3] Z. Zhang, Y. Guo, M. Bescond, J. Chen, M. Nomura, and S. Volz, How coherence is governing diffuson heat transfer in amorphous solids, npj Comput. Mater. 8, 96 (2022). [Download](https://doi.org/10.1038/s41524-022-00776-w).
- [4] Z. Zhang, Y. Guo, M. Bescond, M. Nomura, and S. Volz, J. Chen, Assessing Phonon Coherence Using Spectroscopy, arXiv:2206.12172 (2022). [Download](https://arxiv.org/abs/2206.12172).

## Acknowledgement:

- This work is partially supported by CREST JST (Grant Nos. JPMJCR19I1 and JPMJCR19Q3).
- This project is also supported in part by the National Natural Science Foundation of China (Grant Nos. 12075168 and 11890703), and Science and Technology Commission of Shanghai Municipality (Grant No.19ZR1478600).
- Zhongwei Zhang also thanks for the support from China Scholarship Council (20018-2020).
