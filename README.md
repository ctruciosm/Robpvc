# Robpvc

R implementation of the Robust Principal Volatility Components procedure of Trucíos et al. (2019) [Journal of Empirical Finance]. The RPVC procedure is a robust-to-outliers extension of the procedure proposed by Li et al. (2016) which, by its turn, extends the work of Hu and Tsay (2014). The RPVC procedure should be used in companion with a robust volatility procedure such as, for instance; Trucíos et al. (2017, 2018) among others.

## References

- Hu, Yu-Pin, and Ruey S. Tsay. "Principal volatility component analysis." Journal of Business & Economic Statistics 32.2 (2014): 153-164.
- Li, Weiming, Gao, Jing, Li, Kunpeng, and Yao, Qiwei. "Modeling multivariate volatilities via latent common factors". Journal of Business & Economic Statistics (2016): 34.4 (2016): 564-573.
- Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. "Robust bootstrap forecast densities for GARCH returns and volatilities." Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
- Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. "Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk." Journal of Statistical Computation and Simulation 88.10 (2018): 1976-2000.
- Trucíos, Carlos, Luiz K. Hotta, and Pedro L. Valls Pereira. "On the robustness of the principal volatility components." Journal of Empirical Finance 52 (2019): 201-219.

## Instalation

Robpvc is not available on CRAN yet, but you can install this version using these commands:

install.packages("devtools")

devtools::install_github("ctruciosm/Robpvc")

## Abstract of the RPVC paper

In this paper, we analyse the recent principal volatility components analysis procedure. The procedure overcomes several difficulties in modelling and forecasting the conditional covariance matrix in large dimensions arising from the curse of dimensionality. We show that outliers have a devastating effect on the construction of the principal volatility components and on the forecast of the conditional covariance matrix and consequently in economic and financial applications based on this forecast. We propose a robust procedure and analyse its finite sample properties by means of Monte Carlo experiments and also illustrate it using empirical data. The robust procedure outperforms the classical method in simulated and empirical data. [https://www.sciencedirect.com/science/article/abs/pii/S0927539819300386]
