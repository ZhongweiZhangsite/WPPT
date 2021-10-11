# 'WPPT': A solver for Wave- and Particle-like Phonon Transport (WPPT)

## Authors:

- Zhongwei Zhang <zhongweizhang123@gmail.com> / <zhongwei@tongji.edu.cn>, Institute ofIndustrial Science, The University of Tokyo
- Sebastian Volz <volz@iis.u-tokyo.ac.jp>, LIMMS, CNRS-IIS UMI 2820, The University of Tokyo
- Jie Chen <jie@tongji.edu.cn>, Tongji University
- Masahiro Nomura <nomura@iis.u-tokyo.ac.jp>, Institute ofIndustrial Science, The University of Tokyo

## The code contans three modules:

- NMA: Time-dependent normal mode velocity or displacement calculation.
- PWave: Phonon Wavelet transform for the coherence time distribution, coherence length distribution and phonon decay calculation. (upcoming soon for the coherence length calculation part)
- WPpt: Particlelike lifetimes, coherence corrected lifetimes and coherence time. And, thermal conductivity calculation including wavelike and particlelike phonons. (upcoming soon)

## The code should have interface to:

- LAMMPS
- GPUMD
- VASP & CP2K & DFTB

## Manual:

- Please refers to the the online manual (in preparing): 
  https://zhongweizhang123.wixsite.com/wppt
  
## We will gradually improve the code and prepare more examples for friendly use!

Before that please refers to our publications or preprint about the used theories.

- [1] Z. Zhang, Y. Guo, M. Bescond, J. Chen, M. Nomura, and S. Volz, Generalized Decay Law for Particlelike and Wavelike Thermal Phonons, Phys. Rev. B 103, 184307 (2021).

## Acknowledgement:

- This work is partially supported by CREST JST (No. JPMJCR19I1 and JPMJCR19Q3).
- This project is also supported in part by the grants from the National Natural Science Foundation of China (Grant Nos. 12075168 and 11890703), and Science and Technology Commission of Shanghai Municipality (Grant No.19ZR1478600).
- Zhongwei Zhang also thanks the supports from China Scholarship Council (20018-2020).
