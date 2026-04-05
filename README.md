# Gc-TSEB (Canopy conductance-based Two-Source Energy Balance (TSEB) model)

<img width="805" height="556" alt="image" src="https://github.com/user-attachments/assets/80c5f391-e53f-4f35-9233-2dc3c5e8a8de" />

🧩 Model Framework

The Gc-TSEB model partitions land surface fluxes into canopy and soil components while introducing canopy conductance (Gc) as a physiological constraint on transpiration. In the canopy module, transpiration is directly regulated by stomatal conductance, linking plant physiology to energy balance. The soil module represents evaporation as a function of available energy and soil drying conditions. The model is solved iteratively by updating canopy and soil temperatures until the simulated land surface temperature (LST) converges. Unlike traditional approaches, LST is treated as a model output rather than an input, allowing Gc-TSEB to provide a physically and biologically consistent representation of evapotranspiration dynamics, particularly under water-stressed conditions. 

Compared to earlier versions of Gc-TSEB (Gan et al., 2019), we introduce several key improvements:
(1) canopy transpiration is constrained using SIF-derived photosynthesis coupled with the Medlyn stomatal conductance model, strengthening the link between carbon and water fluxes;
(2) net radiation is partitioned into canopy and soil components using a radiative transfer scheme that accounts for spectral differences between visible and near-infrared wavelengths (Campbell and Norman, 1998; Kustas and Norman, 1999a), reducing biases in sparsely vegetated regions compared to Beer’s Law-based approaches;
(3) soil evaporation is improved by incorporating a moisture–dependent smoothing factor (Zhang et al., 2010), allowing better representation of evaporation pulses following precipitation events in dryland ecosystems.


<img width="1240" height="573" alt="image" src="https://github.com/user-attachments/assets/ca29558d-f2a4-4d4c-be5d-59d838849b30" />




📚 References：

Bu, J., Gan, G. *, Chen, J., Su, Y., Yuan, M., Gao, Y. *, Domingo, F., Lopez-Ballesteros, A., Migliavacca, M., El-Madany, T.S., Gentine, P., Xiao, J., Garcia, M., 2024. Dryland evapotranspiration from remote sensing solar-induced chlorophyll fluorescence: constraining an optimal stomatal model within a two-source energy balance model. Remote Sensing of Environment, 303, 113999.

Gan, G., Kang, T., Yang, S., Bu, J.*, Feng, Z., Gao, Y.*, 2019. An optimized two source energy balance model based on complementary concept and canopy conductance. Remote Sensing of Environment, 223, 243-256. 

Gan, G., Gao, Y.*, 2015. Estimating time series of land surface energy fluxes using optimized two source energy balance schemes: Model formulation, calibration, and validation, Agricultural and Forest Meteorology, 208, 62-75 

Gao, Y., Gan, G.*, Liu, M. and Wang, J., 2016. Evaluating soil evaporation parameterizations at near-instantaneous scales using surface dryness indices. Journal of Hydrology, 541: 1199-1211.

Bu, J., Gan, G.*, Chen, J., Su, Y., Garcia, M., Gao, Y. *, 2021. Biophysical constraints on evapotranspiration partitioning for a conductance-based two source energy balance model. Journal of Hydrology, 603, 127179.
