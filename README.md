---
title: "Update of Kar et al. (2020) Pipeline"
output: github_document
---

# Update of Kar et al. (2020) Pipeline  

This repository provides an updated version of the pipeline originally described in:  

**Kar, S., Tanaka, R., Korbu, L.B., Kholová, J., Iwata, H., Durbha, S.S., Adinarayana, J., Vadez, V., 2020.**  
*Automated discretization of 'transpiration restriction to increasing VPD' features from outdoors high-throughput phenotyping data.*  
Plant Methods 16, 140. https://doi.org/10.1186/s13007-020-00680-8  

🔗 Original pipeline: [EZTr (GitHub)](https://github.com/KSoumya/EZTr)  

---

## About This Update  

This version is based on data from the **“Exp60” trial** carried out at **ICRISAT, Hyderabad, India** (Sept–Oct 2023):  

- Measurement period: *Sept 27 – Oct 17 (excluding 6 irrigation days)*  
- Plant material: *320 sorghum genotypes (RefSet, ICRISAT GenBank)*  
- Leaf area: measured by PlantEye (2 scans/day) and complemented by planimeter measurements at harvest  

During the analysis, we observed:  

- Presence of weight and evapotranspiration outliers **not properly filtered**  
- Calculation of **aberrant** reference evapotranspiration ETref values
- **Typo or command in comment** for the calculation of transpiration Tr and TR.  


These issues motivated us to:  

- **Revise the filtering process** to better handle outliers and improve data quality.  
- **Enhance data visualization** by plotting loadcell data after each major pipeline step, allowing users to track changes and detect anomalies.  
- **Make the pipeline more user-friendly** by converting it into an R Markdown format, with natural stopping points where users should provide inputs. 

---

## Platform Specificity  

The original script was designed for the **LeasyScan platform at ICRISAT**. While similar platforms are now being developed elsewhere (e.g., IRD, Montpellier, France), some steps remain **highly specific to the Indian setup**:  
For example:
- **Loadcell size** may differ between platforms, which affects the conversion of ETr from grams to mm.  
- The pipeline assumes **15-minute intervals between weighings**; different intervals require adjustment.  
- Certain **environmental assumptions** (e.g., filtering night/day time windows) are tailored to tropical conditions and may not apply elsewhere.  

Users applying this pipeline to other LeasyScan platforms should carefully review these steps and adjust parameters accordingly. Currently, the script allows a simplified choice between **ICRISAT** and **IRD** setups. For other platforms, users will need to modify the code to fit their context.  

## Notes  

- This work will be linked to an upcoming publication (**DOI coming soon**).  
- The dataset from **Exp60** will be made publicly available upon publication.  

---

## System Requirements

This pipeline has been tested on the following systems:

- **Operating System:** Windows 10   
- **R Version:** 4.3  
- **Required R Packages:** `tidyverse`, `lubridate`, `ggplot2`, `readr`, `dplyr`, `data.table` (and dependencies)  

---

## Modifications made to the initial version

**1) Outlier detection in the curateRawLC package**

Initial version: detection of weight outliers by boxplot and negative values. 

However, if a LC contains a majority of outliers, the data of this LC will not be removed. For exemple, data of LC "141:03:01" with a mean at 200 kg on the boxplot, instead of a +/- 70kg expected weight, will not be removed.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/66e27a6f-6767-4f4d-98b3-24f84aea5577)



Current version: you defined an expected weight parameter (avg_wgt, in grams) at the beginning of the pipeline. If the weight is outside a range of +/- 30% of this expected weight, it is considered an outlier.


**2) ETref calculation in calculateETref package**
   
Initial version: ETref values too large for the expected order of magnitude. We identified a conversion problem of the solar radiation value in the initial version. In addition, ETref calculation is adapted to a 15 min interval only.

Current version: the conversion problem is fixed and ETref calculation is suitable for all measurement intervals (defined by the seq_by parameter).


**3) PlantEye :**
   
Initial version: There is no visualization of leaf area data measured by the Plant Eye scanner to check whether the 3D-LA data used to calculate Tr and TRrate are consistent. However, Plant Eye can provide unreliable or missing data when plants have reached a certain phenological stage (barcode reading problems, overlapping, etc.). This is why, when we observe the leaf area profile over time, we see leaf area growth at the start of the experiment, followed by an inflection point, then followed by missing or too low data in many cases.

Current version: we propose here to visualize roughly the 3D-LA data of all LC over time, and you can remove by yourself the measurements days when you consider 3D-LA values as outlier (The aim being to keep only the first linear part of the 3D-LA profile.)

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/c98b24f8-6ec1-4c7b-a694-d0b9054deac3)

Once you have confirmed your choice, you can always correct it by looking at the plots of LAI per day/LAI per Timestamp or LA per Timestamp of each LC later on the process. (example here where Oct 10-11, 2023 were kept at first sight, then deleted when viewing the LAI/LA graphs).

![1LAI_over_time](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/7aab842f-edd6-49f0-bb20-b3d94516cd57)

![1LAI_over_time](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/12fce14e-4a14-4487-96c5-27b209af492f)


**4) Management of missing LAI/LA data:**

Initial version: for missing LAI data (often at the end of the experiment), the pipeline calculates the median LAI on the days measured and assigns this value to all missing data. The problem is that this value is therefore lower than the last days measured, and that this value is a constant for the missing days. It therefore does not represent leaf growth.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/aad9c10c-d4f3-4714-9d10-55f18d8e8632)

Current version: with Plant Eye data that the you have chosen to keep, previously converted to observed LA (+ planimeter data if available), the pipeline calculates LA regression as a function of time, assuming linear regression. The slope value is recovered, and the missing values are assigned the previous day's value + slope value. 

**5) Tr calculation in the calculateTr package:**

Tr calculation follows this formula detailled in Kar et al.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/7e35be2f-27d7-4dbb-9263-985442b32e15)

Initial version: Tr= 1-(1-exp(-0.463*LAI)))*ETr

Current version: Tr= (1-exp(-0.463*LAI)))*ETr

**6) TRrate calculation in the calculateTr package:**

Initial version: This part is in comment and therefore not used.

Current version : 

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/f74c56cb-7525-42c6-a15e-b1d4b419a91f)

as described in Kar et al.2020

## Improvements proposed

The main priority for us was to work on data visualization. Indeed, in the initial version of the pipeline, we have just the final outputs. So, we have no visual verification of the filtering process, no clues as to whether the imputed data are logical, whether the transpiration profiles are consistent, etc.. 

**1) You can plot the results of an LC to see the outlier filtering process, and check that the profiles obtained are correct.**

Ex: Weight profile before and after filtering

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/0c5c45d2-9acc-45f1-8f27-a8fef92c1a9c)
![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/37fe5ffe-144d-4cdf-ad9e-8ab4683d6527)

Ex : Evapotranspiration profile before and after filtering 

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/b6a3443a-60ba-4ae9-b901-f99584fa6b93)
![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/fdc64e8d-420c-4863-a2cc-965927eff225)

**2) Visualise data one LC:**
- Leaf area over time (LA, in m²)
  
![1LA_over_time](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/0f39c90c-e052-4bf9-9fdd-b75304af5147)

- Evapotranpsiration (ETr non normalized, mm.15min-1) and Reference Evapotranspiration (mm.15min-1) *to check if ETr is lower than ETref*
  
![1_ETr_before_norm](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/dc86106e-2281-4202-88d5-ca297d7dda81)

- Transpiration (Tr, mm.15min-1)
  
![1_Transpiration](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/78c329d6-35f5-42f6-b720-119d0ad29ce2)

- Transpiration and ETr (mm.15min-1) *to check if Transpiration is lower than Evapotranspiration*

![ETr_1_ETr_and_Transp](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/fd2deeed-e583-4506-83ac-002a4976b00d)

- Transpiration rate (TRrate, mm.15min-1.m-²)

![1_TRrate](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/a8ede257-13ec-448c-9145-c03d2816e50e)

(all the profiles (weights/ETr and Tr/leaf area) are available in png format in the folder results/profile.

**3) Save/load objects created in the process** (folder /saved_objects/)

**4) Outlier detection :**

Count of negative and max values before and after filtering processes, number of remaining LC, 

**5) In the generateETrRatio package of the initial version**

ETr filtering according to ETref value is based on an hourly system to define the diurnal (6:30am-6:30pm) and nocturnal (6:30pm-6:30am) periods. This works in the case of India, but cannot be exported for other conditions. We have therefore chosen to base this filter on solar radiation.

- If SR<0  night time  ETr= 0
- If SR>0  daytime  -0.01<ETr<ETref 

**6) Filtering ETr/TRrate according to profile**

After filtering, we observed that some high peaks remained on the ETr profile normalized by LA or TRrate. The weather data were similar from one day to the next, that did not allow us to determine the cause of this difference. This does not concern a particular day. It affects many LCs, but not all. We decided to retrieve the peak max for each day of measurement and make a boxplot with the max values. If a value is considered an outlier on the boxplot, then the entire day of measurement is removed. 



## Platform-Specific Notes

The pipeline was originally optimized for the **LeasyScan platform at ICRISAT**. Key points:

- **ET conversion:** From grams to mm, based on ICRISAT LC area (¼ m²) – handled in `convETr`
- **Leaf area calculations:** LAI = 0.26 is used for Tr and TRrate calculations – in `calculateTr`

**Platform selection:**

- Users must select which platform they are using: **ICRISAT** or **IRD** in the update. The platform-dependent steps are platform-specific and will adjust automatically based on this selection
- If using another LeasyScan platform, users should review and adapt these steps manually. (Upcoming improvement: Make the pipeline **adaptable to all LeasyScan platforms**, regardless of the specificity of the platform ? )

**Upcoming improvement?**


The equations to convert 3DLA in observed LA come from the article Vadez et al 2015. These equations are based on three crops: pearl millet, cowpea and peanut. In the initial version of the pipeline, the equation is correct in the case of pearl millet, incorrect in the other cases (cowpea and peanut)
It has been corrected and you can decide which crop (between pearl millet, cow pea and peanut) what corresponds best to yourcase. In addition, it would be interesting to have more crops proposed in order to refine the relationship.







