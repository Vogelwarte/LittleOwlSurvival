# LittleOwlSurvival
Analysis of first-year survival of Little Owls in Germany.

This repository contains the data and code to replicate the analysis and results of the following scientific article: (Juvenile survival of little owls decreases with snow cover)[https://doi.org/10.1002/ece3.11379]

Raw and prepared data of little owl telemetry detections are contained in the folder /data, while JAGS models that were used to estimate biweekly survival probability are contained in the folder /models. The R script "LIOW_telemetry_data_prep.r" details the various steps taken to prepare the raw detection data for analysis, and the R script "LIOW_telemetry_survival_reduced.r" performs the actual analysis reported in the article.

The folder /output contains graphical and tabular output from the survival model as well as an R workspace "LIOW_survival_output.RData" that can be used to replicate output graphs and tables as it already contains the model output objects.

To replicate the analysis the user must have a recent R installation (>R3.5.0)[https://mirrors.cicku.me/cran/] and a recent JAGS installation (v4.3.1)[https://sourceforge.net/projects/mcmc-jags/files/]. For any questions please contact the corresponding author of the (associated study)[https://doi.org/10.1002/ece3.11379].

## Summary of the study:
Global environmental changes are associated with warmer average temperatures and more extreme weather events, potentially affecting wildlife population dynamics by altering demographic processes. Extreme weather events can reduce food resources and survival in all seasons of the year. Estimates of season-specific survival probabilities are therefore crucial to understand the moderating effect of extreme events on annual mortality. Here, we analysed survival probabilities of 307 radio-tracked juvenile little owls (Athene noctua) over two-week periods from fledging to their first breeding attempt in the following spring to assess the contribution of extreme weather events. Survival probabilities were typically lowest during the first weeks after fledging in summer, but were moderated by seasonal extremes in winter. The duration of snow cover in winter had a strong negative effect on survival probability, while being food-supplemented during the nestling stage increased survival during the first weeks after fledging in summer and ultimately led to a larger proportion of birds surviving the first year. Overall annual survival probability over the first year varied by 34.3% between 0.117 (95% credible interval 0.052 – 0.223) and 0.178 (0.097 – 0.293) depending on the severity of the winter, and was as high as 0.233 (0.127 – 0.373) for food-supplemented fledglings. In years with mild winters, the season with the lowest survival was the summer post-fledging period (0.508; 0.428 – 0.594), but in years with extensive snow cover the winter was the season with the lowest survival (0.481; 0.337 – 0.626). We therefore show that extreme weather events occurring in a particular season reduced the proportion of first-year survivors. Increasing extreme weather events can moderate seasonal survival probability through altering food supply of juvenile little owls either during the nestling period or in winter, with similarly large effects on annual survival and the viability of populations.


