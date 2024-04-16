# Update of Kar et al 2020 pipeline 

This project concerns an update of the pipeline described in the article :Kar, S., Tanaka, R., Korbu, L.B., Kholová, J., Iwata, H., Durbha, S.S., Adinarayana, J., Vadez, V., 2020. Automated discretization of 'transpiration restriction to increasing VPD' features from outdoors high-throughput phenotyping data. Plant Methods 16, 140. https://doi.org/10.1186/s13007-020-00680-8

Link GITHUB : https://github.com/KSoumya/EZTr

This script is based on the analysis of "Exp60" carried out at ICRISAT, Hyderabad, India, in September-October 2023. The measurement period is from September 27 to October 17, with 6 irrigation days removed. 
Weight measurements (Loadcells measurements) and climate data are recorded on a 15 min interval. The data are coupled with PlantEye data measuring 3D-LA leaf area twice a day, and planimeter data (observed LA) at the end of the trial. 
The plant material is composed of 320 sorghum genotypes of RefSet list defined by the Genbank at ICRISAT.

During the analysis of this experiment, we observed some data outliers that were not cleaned, and aberrants ETref values.We therefore decided to take a closer look at the filtering stages and identified certain steps where modifications would be necessary. We present below the modifications we made to the initial script, followed by some simple suggested improvements. We also noticed that the pipeline was adapted to the specific case of the ICRISAT platform. A number of enhancements are envisaged as possible improvements (not yet present in the script) 


## Modifications made to the initial version

**1) Outlier detection in the curateRawLC package:**

Initial version: detection of weight outliers by boxplot and negative values. However, if the LC contains a majority of upper outliers (e.g. 200 kg), the LC will not be removed.
![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/8387dcdd-0f9f-481d-9167-305f35c00421)

Current version: an expected weight parameter (avg_wgt) is defined by the user at the start of the pipeline. If the weight is outside a range of +/- 30% of this expected weight, it is considered an outlier.


**2) ETref calculation in calculateETref package:**
   
Initial version: ETref values too large for the expected order of magnitude. We identified a conversion problem in the initial version concerning the solar radiation. In addition, ETref calculation is adapted to a 15 min interval only.

Current version: the conversion problem is fixed, and the ETref calculation is suitable for all measurement intervals (because it is defined by the seq_by parameter).


**3) PlantEye:** 

Initial version: There is no visualization of leaf area data measured by the Plant Eye scanner to check whether the 3D-LA data used to calculate Tr and TRrate are consistent. However, Plant Eye can provide aberrant/missing data when plants have reached a certain phenological stage (barcode reading problems, overlapping, etc.). In this case, pearl millet crop )

Current version: we propose here to visualize the 3D-LA data by graph, and the user removes by himself/herself the days he/she considers as outlier.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/c98b24f8-6ec1-4c7b-a694-d0b9054deac3)

Once the user has confirmed his choice, he can always correct it by looking at the plots of LAI per day/LAI per Timestamp or LA per Timestamp of each LC. (example here where Oct 10-11, 2023 were kept at first sight, then deleted when viewing the LAI/LA graphs).

![1LAI_over_time](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/dbe12962-bca7-43f7-9ad1-3ae886a8430d)


4) Management of missing LAI/LA data. 
Initial version: for missing LAI data (often at the end of the experiment), the pipeline calculates the median LAI on the days measured and assigns this value to all missing data. The problem is that this value is therefore lower than the last days measured, and that this value is a constant for the missing days. It therefore does not represent leaf growth.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/aad9c10c-d4f3-4714-9d10-55f18d8e8632)

Current version: with Plant Eye data that the user has chosen to keep, previously converted to observed LA (+ planimeter data if available), the pipeline calculates LA regression as a function of time, assuming linear regression. The slope value is recovered, and the missing values are assigned the previous day's value + slope value. 

5) Tr calculation in the calculateTr package:

Tr calculation follows this formula destailled in Kar et al.

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/7e35be2f-27d7-4dbb-9263-985442b32e15)

Initial version: Tr= 1-(1-exp(-0.463*LAI)))*ETr

Current version: Tr= (1-exp(-0.463*LAI)))*ETr

## IMPROVEMENTS PROPOSED

1) The user can plot the results of an LC to see the outlier filtering process, and check that the profiles obtained are correct.
  Ex :Weight profile before and after filtering

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/0c5c45d2-9acc-45f1-8f27-a8fef92c1a9c)
![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/37fe5ffe-144d-4cdf-ad9e-8ab4683d6527)

Ex : Evapotranspiration profile befor and after filtering 

![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/b6a3443a-60ba-4ae9-b901-f99584fa6b93)
![image](https://github.com/lcmgregoire/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/assets/96241863/fdc64e8d-420c-4863-a2cc-965927eff225)

2) Save/load objects created in the process (folder /saved_objects/)

3) Outlier detection: count of negative and max values before and after filtering processes

4) In the generateETrRatio package of the initial version, ETr filtering according to ETref value is based on an hourly system to define the diurnal (6:30am-6:30pm) and nocturnal (6:30pm-6:30am) periods. This works in the case of India, but cannot be exported for other conditions. We have therefore chosen to base this filter on solar radiation.

- If SR<0  night time  ETr= 0
- If SR>0  daytime  -0.01<ETr<ETref 

5) Filtering according to TRrate profile

After filtering, we observed that some outliers remained on the ETr profile normalized by LA or TRrate on certain days... The measured weather does not allow us to determine the cause of this difference, as the days are similar. This does not concern a particular day. It affects many LCs, but not all. We decided to retrieve the peak max for each day of measurement and make a boxplot with the max values for each day. If a value is considered an outlier on the boxplot, then the entire day of measurement is removed. 

## IMPROVEMENTS TO COME : 

1) ET conversion from grams to mm is based on the ICRISAT LC area, ¼ m² (in the convETr package)
  
2) LAI and LA calculations (and therefore affect the value of Tr and TRrate) are based on the LAI of 0.26, adapted to the ICRISAT LC.
  
3) Conversion calculations from the 3D value estimated by PlantEye and the value observed on the planimeter. The equations come from the article Vadez et al 2015, based on three crops: pearl millet, cowpea and peanut. In al initial version of the pipeline, the equation is correct in the case of pearl millet, incorrect in the other cases (cowpea and peanut). It would be wise to i) correct the equation for both crops, ii) make a loop so that the user can decide which crop corresponds best to his case, iii) In addition, it would be interesting to have more crops proposed in order to refine the relationship. 

