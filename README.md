# Active Noise Control for Impulsive Noise

This repository provides MATLAB implementations of various **Feedforward Active Noise Control (ANC)** algorithms designed to suppress **impulsive noise**. These algorithms are evaluated and compared in terms of their Adaptive Noise Reduction (ANR) performance.

![System Overview](ffanc.png)

## Requirements

- MATLAB R2024b
- Signal Processing Toolbox

## Algorithms

All algorithm implementations are located in the `algorithms/` directory.

| Abbreviation            | Description                                                                   | Reference |
|-------------------------|-------------------------------------------------------------------------------|-----------|
| **FxNLMS**              | Standard Filtered-x Normalized LMS                                             | –         |
| **FxgsnLMS**[*1]        | Gaussian Step-size Normalized LMS                                              | [1]       |
| **FxLMM**[*2]           | Hampel-type Least Mean M-estimate Filtered-x LMS                               | [2]       |
| **MFxLMM**[*3]          | Enhanced M-estimate with dynamic σ update                                     | [3]       |
| **Th-FxNLMS**           | Thresholded FxNLMS                                                             | –         |
| **MFxLCH**[*4]          | Median Filtered Hybrid LMS with L-Criterion                                   | [4]       |
| **Fair**[*5]            | Fair-type M-estimator adaptive algorithm                                       | [5]       |
| **FxlogNLMS**[*6]       | Logarithmic-transformed FxNLMS                                                 | [6]       |
| **NS-FxlogLMS**[*7]     | FxlogLMS with Nonlinear Selection                                               | [7]       |
| **FxlogNLMS-plus**[*8]  | Improved FxlogNLMS with stability and speed-up modifications                   | [8]       |
| **NSS-FxlogLMS-plus**[*8]| Noise-shaped step-size version of FxlogNLMS-plus                             | [8]       |

### References

[*1] Y. Zhou, Q. Zhang, and Y. Yin,  
"Active control of impulsive noise with symmetric α-stable distribution based on an improved step-size normalized adaptive algorithm,"  
*Mechanical Systems and Signal Processing*, vol. 56, pp. 320–333, 2015.  

[*2] T. Thanigai, S. M. Kuo, and R. Yenduri,  
"Nonlinear active noise control for infant incubators in neonatal intensive care units,"  
in *Proc. ICASSP*, 2007, pp. 109–112.  

[*3] G. Sun, M. T. Li, and T. C. Lim,  
"Enhanced filtered-x least mean M-estimate algorithm for active impulsive noise control,"  
*Applied Acoustics*, vol. 90, pp. 31–41, 2015.  

[*4] A. Mirza, A. Zeb, M. Y. Umair, et al.,  
"Less complex solutions for active noise control of impulsive noise,"  
*Analog Integrated Circuits and Signal Processing*, vol. 102, no. 3, pp. 507–521, 2020.  

[*5] L. Wu and X. Qiu,  
"An M-estimator based algorithm for active impulse-like noise control,"  
*Applied Acoustics*, vol. 74, pp. 407–412, 2013.  

[*6] L. Wu, H. He, and X. Qiu,  
"An active impulsive noise control algorithm with logarithmic transformation,"  
*IEEE Trans. Audio, Speech, and Language Processing*, vol. 19, no. 4, pp. 1041–1044, 2011.  

[*7] M. Pawelczyk, W. Wierzchowski, L. Wu, and X. Qiu,  
"An extension to the filtered-x LMS algorithm with logarithmic transformation,"  
in *Proc. ISSPIT*, 2015, pp. 454–459.  

[*8] A. Haneda, Y. Sugiura, and T. Shimamura,  
"FxlogLMS+: Modified FxlogLMS Algorithm for Active Impulsive Noise Control,"  
*Lecture Notes in Electrical Engineering*, vol. 1322, Springer, 2025, pp. 342–351.  
https://link.springer.com/chapter/10.1007/978-981-96-1535-3_34

## Usage

1. Downlaad all files

2. Open MATLAB and run the following script:

```matlab
main.m
```

3. The script will:

- Load predefined ANC scenarios with impulsive noise
- Run all methods sequentially
- Plot the Adaptive Noise Reduction (ANR) curves
- Save the result as `comparison.png`

## Results

The figure below shows the comparative performance of each algorithm under impulsive noise (α-stable noise with α=1.65).

![Comparison Result](comparison.png)

---

For questions or collaborations, feel free to contact [Yosuke Sugiura](mailto:ysugi@xxxx.ac.jp).
