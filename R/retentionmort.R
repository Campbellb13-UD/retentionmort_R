#' Estimating user-based tagging mortality and tag shedding in field
#' mark-recapture studies
#'
#' This model estimates the percent loss in tagged animals at large for
#' field-based recapture studies based on a linear decrease in survival and tag
#' retention (including lost tags and missidentified tags) for five weeks per
#' tagging cohort based on laboratory retention/survival studies. The
#' `retentionmort()` function can be used following a recapture field study to
#' estimate user-based tag loss in animals at large. The model is changed by
#' linear regression coefficients of weekly tag loss rate, weekly mortality
#' rate, and their respective intercepts. The coefficients used can be selected
#' from the currently included list using the `err` input or be customized.
#' This function is also capable of working with a cofactor with two conditions
#' (e.g. class1 individuals and large individuals) to improve resolution for
#' more specified studies.
#'
#' @param nT           A vector of the number of tagged individuals for each
#'                     tagging effort.
#'
#' @param n_c1         (optional) A vector of the number of tagged individuals
#'                     in one of two categorical variables. If this is not being
#'                     used then n_c1 will be equal to nT.
#'
#' @param TaL          A vector of the cumulative number of tagged individuals
#'                     following each effort.
#'
#' @param c            One value for the total number of tagging efforts.This
#'                     value must be greater than or equal to 6.
#'
#' @param R            A vector of the number of recaptured individuals per
#'                     effort.
#'
#' @param err          A value (between 1 and 26) that represents the weekly
#'                     mortality rate and weekly tag loss rate from a preloaded
#'                     case study listed in the metadata. Alternatively, model
#'                     coefficients can be manually included using a combination
#'                     of the preceding parameters. While the preloaded data are
#'                     based on weekly time stamps, customized model
#'                     coefficients can reflect any time period specified and
#'                     the projection will predict loss at 5 times the time
#'                     interval.
#'
#'         1 = Large (> 61mm TL) and class1 (< 61mm TL) Mummichogs tagged with
#'             VIE in caudal peduncle (avg + 95% CI) - McCutcheon et al. in prep
#'         2 = Large (> 61mm TL) and class1 (< 61mm TL) Mummichogs tagged with
#'             VIE in caudal peduncle (avg) - McCutcheon et al. in prep
#'         3 = Large (> 61mm TL) and class1 (< 61mm TL) Mummichogs tagged with
#'             VIE in caudal peduncle (avg - 95% CI) -  McCutcheon et al. in
#'             prep
#'         4 = American Eel elvers (80 – 149 mm TL) tagged with 2 VIE tags in
#'             anterior, posterior, central of body - Eissenhauer et al. 2024
#'             https://doi.org/10.1002/nafm.11016
#'         5 = Mummichogs (45 - 82 mm TL) tagged with 8mm PIT tags in abdominal
#'             cavity - Kimball & Mace 2020
#'             https://doi.org/10.1007/s12237-019-00657-4
#'         6 = Mummichogs (45 - 82 mm TL) tagged with 12mm PIT tags in abdominal
#'             cavity - Kimball & Mace 2020
#'             https://doi.org/10.1007/s12237-019-00657-4
#'         7 = Pinfish (45 - 82 mm TL) tagged with 8mm or 12mm PIT tags in
#'             abdominal cavity - Kimball & Mace 2020
#'             https://doi.org/10.1007/s12237-019-00657-4
#'         8 = Cichlids (29 - 59 mm TL) tagged with VIE in various locations on
#'             body - Jungwirth et al. 2019
#'             https://doi.org/10.1007/s00265-019-2659-y
#'         9 = River Shiners (36 - 49 mm TL) tagged with VIE using anesthesia in
#'             various locations - Moore & Brewer 2021
#'             https://doi.org/10.1002/nafm.10607
#'        10 = River Shiners (50 - 56 mm TL) tagged with 8 mm PIT using
#'             anesthesia in various locations - Moore & Brewer 2021
#'             https://doi.org/10.1002/nafm.10607
#'        11 = River Shiners (40 - 51 mm TL) tagged with VIE using no anesthesia
#'             in various locations - Moore & Brewer 2021
#'             https://doi.org/10.1002/nafm.10607
#'        12 = River Shiners (50 - 55 mm TL) tagged with 8 mm PIT using no
#'             anesthesia in various locations - Moore & Brewer 2021
#'             https://doi.org/10.1002/nafm.10607
#'        13 = Delta Smelt (> 70 mm FL) tagged with injected acoustic tag -
#'             Wilder et al. 2016
#'             https://doi.org/10.1080/02755947.2016.1198287
#'        14 = Delta Smelt (> 70 mm FL) surgically tagged with acoustic tag -
#'             Wilder et al. 2016
#'             https://doi.org/10.1080/02755947.2016.1198287
#'        15 = Rohu Carp tagged with floy tags under dorsal fin - Hadiuzzaman
#'             et al. 2015
#'             https://www.researchgate.net/publication/289460932_Feasibility_study_of_using_floy_tag_and_visible_implant_fluorescent_elastomer_marker_in_major_carps
#'        16 = Silver Carp tagged with floy tags under dorsal fin - Hadiuzzaman
#'             et al. 2015
#'             https://www.researchgate.net/publication/289460932_Feasibility_study_of_using_floy_tag_and_visible_implant_fluorescent_elastomer_marker_in_major_carps
#'        17 = Black Bullhead (mean TL = 153.3 mm) tagged with VIE near dorsal
#'             fin - Schumann et al. 2013
#'             https://benthamopen.com/contents/pdf/TOFISHSJ/TOFISHSJ-6-41.pdf
#'        18 = Bluegill (mean TL = 75.8 mm) tagged with VIE near dorsal fin -
#'             Schumann et al. 2013
#'             https://benthamopen.com/contents/pdf/TOFISHSJ/TOFISHSJ-6-41.pdf
#'        19 = Channel Catfish (mean TL = 127.9 mm) tagged with VIE near dorsal
#'             fin - Schumann et al. 2013
#'             https://benthamopen.com/contents/pdf/TOFISHSJ/TOFISHSJ-6-41.pdf
#'        20 = Juvenile Burbot (88 - 144 mm TL) tagged with coded wire tag on
#'             snout, periocular region, nape, pectoral fin base, dorsal fin
#'             base, and anal fin base - Ashton et al. 2013
#'             https://doi.org/10.1080/02755947.2014.882458
#'        21 = Delta Smelt adults (45 - 77 mm FL) and juveniles (20 - 40 mm FL)
#'             tagged with calein markers - Castillo et al. 2014
#'             https://doi.org/10.1080/02755947.2013.839970
#'        22 = Juvenile Seabass (mean 173 g) tagged with dummy acoustic
#'             transmitters in external or intraperitoneal cavity -
#'             Begout Anras et al. 2003
#'             https://doi.org/10.1016/S1054-3139(03)00135-8
#'        23 = Juvenile American Eels (113 - 175 mm TL) tagged with
#'             micro-acoustic transmitter in body cavity - Mueller et al. 2017
#'             https://doi.org/10.1016/j.fishres.2017.06.017
#'        24 = Juvenile European Eels (7 - 25 g) tagged with 12mm PIT tags -
#'             Jepsen et al. 2022
#'             https://doi.org/10.1111/jfb.15183
#'        25 = Adult Atlantic Croaker (147 - 380 mm TL) tagged with VIE tags in
#'             caudal fin - Torre et al. 2017
#'             https://doi.org/10.1080/00028487.2017.1360391
#'        26 = Adult Spot (65 - 222 mm FL) tagged with VIE tags in caudal fin -
#'             Torre et al. 2017
#'             https://doi.org/10.1080/00028487.2017.1360391
#'
#' @param m_mort_c1    A value that represents the slope of the mortality rate
#'                     for the class1 individuals. While all the preloaded
#'                     datasets work in weekly time intervals, these can be
#'                     customized to any time interval to match the sampling
#'                     interval. The resulting model will then project mortality
#'                     and tag loss for 5 times the time interval. If this value
#'                     is added, then, at minimum, `b_mort_c1`, `m_ret_c1`, and
#'                     `b_ret_c1` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param b_mort_c1    A value that represents the intercept of the mortality
#'                     rate for the class1 individuals. While all the preloaded
#'                     datasets work in weekly time intervals, these can be
#'                     customized to any time interval to match the sampling
#'                     interval. The resulting model will then project mortality
#'                     and tag loss for 5 times the time interval. If this value
#'                     is added, then, at minimum, `m_mort_c1`, `m_ret_c1`, and
#'                     `b_ret_c1` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param m_ret_c1     A value that represents the slope of the tag loss
#'                     (represented as tag loss and missidentification) rate for
#'                     the class1 individuals. While all the preloaded datasets
#'                     work in weekly time intervals, these can be customized to
#'                     any time interval to match the sampling
#'                     interval. The resulting model will then project mortality
#'                     and tag loss for 5 times the time interval. If this value
#'                     is added, then, at minimum, `m_mort_c1`, `b_mort_c1`, and
#'                     `b_ret_c1` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param b_ret_c1     A value that represents the intercept of the tag loss
#'                     (represented as tag loss and missidentification) rate for
#'                     the class1 individuals. While all the preloaded datasets
#'                     work in weekly time intervals, these can be customized to
#'                     any time interval to match the sampling interval. The
#'                     resulting model will then project mortality and tag loss
#'                     for 5 times the time interval. If this value is added,
#'                     then, at minimum, `m_mort_c1`, `b_mort_c1`, and
#'                     `m_ret_c1` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param m_mort_c2    A value that represents the slope of the mortality rate
#'                     for the class2 individuals. While all the preloaded
#'                     datasets work in weekly time intervals, these can be
#'                     customized to any time interval to match the sampling
#'                     interval. The resulting model will then project mortality
#'                     and tag loss for 5 times the time interval. If this value
#'                     is added, then, at minimum, `b_mort_c2`, `m_ret_c2`, and
#'                     `b_ret_c2` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param b_mort_c2    A value that represents the intercept of the mortality
#'                     rate for the class1 individuals. While all the preloaded
#'                     datasets work in weekly time intervals, these can be
#'                     customized to any time interval to match the sampling
#'                     interval. The resulting model will then project mortality
#'                     and tag loss for 5 times the time interval. If this value
#'                     is added, then, at minimum, `m_mort_c2`, `m_ret_c2`, and
#'                     `b_ret_c2` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param m_ret_c2     A value that represents the slope of the tag loss
#'                     (represented as tag loss and missidentification) rate for
#'                     the class1 individuals. While all the preloaded datasets
#'                     work in weekly time intervals, these can be customized to
#'                     any time interval to match the sampling interval. The
#'                     resulting model will then project mortality and tag loss
#'                     for 5 times the time interval. If this value is added,
#'                     then, at minimum, `m_mort_c2`, `b_mort_c2`, and
#'                     `b_ret_c2` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @param b_ret_c2     A value that represents the intercept of the tag loss
#'                     (represented as tag loss and missidentification) rate for
#'                     the class1 individuals. While all the preloaded datasets
#'                     work in weekly time intervals, these can be customized to
#'                     any time interval to match the sampling interval. The
#'                     resulting model will then project mortality and tag loss
#'                     for 5 times the time interval. If this value is added,
#'                     then, at minimum, `m_mort_c2`, `b_mort_c2`, and
#'                     `m_ret_c2` need to be used. The use of these coefficients
#'                     will override the `err` term.
#'
#' @returns  This returns a dataframe `datacomp` that contains summary
#'           information from each mark-recapture effort, several parameters
#'           used in the calculation of adjusted recaptures, and basic error
#'           values between the expected and observed number of recaptures. The
#'           `datacomp` dataframe can be used in the `retentionmort_figure`
#'           function to generate some preliminary figures that can be used to
#'           assess model performance and factors that influence the error
#'           between expected and observed recaptures.
#'
#'           Values that will be returned include:
#'             week =      The week in the study
#'             nT =        The number of tagged individuals per tagging effort
#'             n_c1 =      The number of tagged individuals in class1 per tagging
#'                         effort
#'             TaL =       The cumulative number of tagged individuals at large
#'             TAsum =     The weekly sum of adjusted number of tags at large
#'             TDF =       Tag depreciation factor
#'             YSs =       Resultant survival rate of class1
#'             YSl =       Resultant survival rate of class2
#'             YMs =       Resultant tag loss rate of class1
#'             YMl =       Resultant tag loss rate of class2
#'             TaLs =      The cumulative number of class1 individuals tagged at
#'                         large
#'             TaLl =      The cumulative number of class2 individuals tagged at
#'                         large
#'             R =         The number of recaptured individuals per effort
#'             Rpercent =  The proportion of recaptured individuals
#'             RA =        The adjusted number of recaptured individuals
#'             PSE =       The percent standard error between observed and
#'                         estimated recaptured individuals
#'
#' @examples
#'
#' #To formulate a dataset for each example
#'   list2env(test_dataset_retentionmort(), envir = .GlobalEnv)
#'
#'
#' #Using a preloaded set of model parameters
#'   datacomp = retentionmort(n_c1=n_c1, nT=nT, TaL=TaL, c=c, R=R, err=2
#'                         )
#'
#' #Using custom model parameters for one class type
#'   datacomp = retentionmort(n_c1=n_c1, nT=nT, TaL=TaL, c=c, R=R,
#'                         m_mort_c1=-0.0625,b_mort_c1=1.06, m_ret_c1=-0.113,
#'                         b_ret_c1=1.05
#'                         )
#'             #or
#'   datacomp = retentionmort(n_c1=n_c1, nT=nT, TaL=TaL, c=c, R=R,
#'                         m_mort_c2=-0.0203, b_mort_c2=1.03, m_ret_c2=-0.0541,
#'                         b_ret_c2=0.993
#'                         )
#'
#' #Using custom model parameters for two class types
#'   datacomp = retentionmort(n_c1=n_c1, nT=nT, TaL=TaL, c=c, R=R,
#'                         m_mort_c1=-0.0625, b_mort_c1=1.06, m_ret_c1=-0.113,
#'                         b_ret_c1=1.05, m_mort_c2=-0.0203, b_mort_c2=1.03,
#'                         m_ret_c2=-0.0541, b_ret_c2=0.993
#'                         )
#' @export
retentionmort<- function(nT,                 # number of tagged individuals per tagging effort
                         n_c1 = nT,          # number of tagged individuals in class1 per effort (default = nT)
                         TaL,                # cumulative number of tagged individuals per effort
                         c,                  # number of sampling cohorts - should be equal to the length of each vector
                         R,                  # number of recaptured individuals per effort
                         err = 2,            # choose error values to use from metadata list (default McCutcheon et al. - average case)
                         m_mort_c1 = NA,     # slope of weekly mortality rate for class1 (default NA)
                         b_mort_c1 = NA,     # intercept of weekly mortality rate for class1 (default NA)
                         m_ret_c1 = NA,      # slope of weekly tag loss rate for class1 (default NA)
                         b_ret_c1 = NA,      # intercept of weekly tag loss rate for class1 (default NA)
                         m_mort_c2 = NA,     # slope of weekly mortality rate for class2 (default NA)
                         b_mort_c2 = NA,     # intercept of weekly mortality rate for class2 (default NA)
                         m_ret_c2 = NA,      # slope of weekly tag loss rate for class2 (default NA)
                         b_ret_c2 = NA       # intercept of weekly tag loss rate for class1 (default NA)
  ){

  #1. Establish known coefficients based on laboratory studies - custom values can be added to the end of the mSs, bSs, mSl, bSl, mMs, bMs, mMl, and bMl vectors based on weekly tag loss rate and mortality rate.
  mSs <- c(-0.313,-0.0625,0, -0.001575,-0.014799, -0.029598, -0.009259, -0.01925, -0.027417, -0.01925, -0.030917, -0.1, -0.125, -0.0105, -0.014, -0.005133, -0.036867, -0.042, -0.004375,-0.01, -0.06,0,0,-0.012895,0,-0.00000001)    #slope of survival - class1
  mSs<- matrix(mSs, ncol=1)
  bSs <- c(1.31, 1.06, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)       #intercept of survival - class1
  bSs<- matrix(bSs, ncol=1)
  mSl <- c(-0.0697, -0.0203, 0, -0.001575, -0.014799, -0.029598, -0.009259, -0.01925, -0.027417, -0.01925, -0.030917, -0.1, -0.125, -0.0105, -0.014, -0.005133, -0.036867, -0.042, -0.004375, -0.01, -0.06,0,0,-0.012895,0,-0.000000001) #slope of survival - class2
  mSl<- matrix(mSl, ncol=1)
  bSl <- c(1.009, 1.03, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)      #intercept of survival - class2 individuals
  bSl<- matrix(bSl, ncol=1)
  mMs <- c(-0.245, -0.113, 0, -0.000875, -0.016279, -0.036176, -0.010101, -0.007, -0.029167, -0.004667, -0.025083, -0.0125, 0, -0.013067, -0.046083, -0.003967, -0.016333, -0.021, -0.004375,-0.01, 0, -0.089362,-0.0133,-0.000921,0,-0.001167)   #slope of tag retention - class1
  mMs<- matrix(mMs, ncol=1)
  bMs <- c(0.906, 1.05, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)      #intercept of tag retention - class1
  bMs<- matrix(bMs, ncol=1)
  mMl <- c(-0.0118, -0.0541, 0, -0.000875, -0.016279, -0.036176, -0.010101, -0.007, -0.029167, -0.004667, -0.025083, -0.0125, 0, -0.013067, -0.046083, -0.003967, -0.016333, -0.021, -0.004375, -0.01, 0, -0.089362,-0.0133,-0.000921,0,-0.001167) #slope of tag retention - class2 individuals
  mMl<- matrix(mMl, ncol=1)
  bMl <- c(0.716, 0.993, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)     #intercept of tag retention - class2 individuals
  bMl<- matrix(bMl, ncol=1)

  datacomp = data.frame() # make datacomp to save all outputs

  # Set up each timestep
  t<- c(0:5)                       #for the 5 week time step
  tadj<- c(0:5, rep(5, c-6))       #used to force no more regression after 5 weeks
  ltadj<- length(tadj)

  #2. Add errors to check inputs
  # check n_c1
  if(!is.vector(n_c1))
    stop("n_c1 must be the number of class1 individuals each tagging effort")
  # check nT
  if(!is.vector(nT))
    stop("nT must be a vector of total number of tags per effort")
  # check TaL
  if(!is.vector(TaL))
    stop("TaL must be a vector of...")
  # check c
  if(!is.numeric(c) | !c >= 6)
    stop("c must be a number of tagging efforts and must be at least 6")
  if(c != length(nT) | c != length(n_c1) | c != length(TaL) | c != length(R))
    stop("c must be equal to the length of nT, n_c1, TaL, and R")
  # check R
  if(!is.vector(R))
    stop("R must be a vector of % recaptures")
  if(!R[1] == 0)
    stop("the first recapture event (R) must be 0")
  #check err
  if((err > 26 | err < 1)&(!is.na(err)))
    stop("err can only be a value between 1 and 26 or be NA if custom coefficients are being used")
  # check m coefficients
  if((m_mort_c1 > 0 & !is.na(m_mort_c1)) | (m_mort_c2 > 0 & !is.na(m_mort_c2)) | (m_ret_c1 > 0 & !is.na(m_ret_c1)) | (m_ret_c2 > 0 & !is.na(m_ret_c2)))
    print("WARNING: m values should be a negative value between -1 and 0 or NA if not used")
  #check b coefficients
  if(((b_mort_c1 > 1.5 | b_mort_c1 < 0.5) & !is.na(m_mort_c1) ) | ((b_mort_c2 > 1.5 | b_mort_c2 < 0.5 )& !is.na(m_mort_c2)) | ((b_ret_c1 > 1.5 | b_ret_c1 < 0.5  )& !is.na(m_ret_c1)) | ((b_ret_c2 > 1.5 | b_ret_c2 < 0.5) & !is.na(m_ret_c2) ))
    print("WARNING: b values should generally fall near 1 or be NA if not used. 1 can often be used for b values if one is not predefined")
  #coefficients
  if(!is.na(m_mort_c1) & !is.na(b_mort_c1) & !is.na(m_mort_c2) & !is.na(b_mort_c2) & !is.na(m_ret_c1) & !is.na(b_ret_c1) & !is.na(m_ret_c2) & !is.na(b_ret_c2) & is.na(err) ) {
    err = 27
    mSs[27] = m_mort_c1
    bSs[27] = b_mort_c1
    mSl[27] = m_mort_c2
    bSl[27] = b_mort_c2
    mMs[27] = m_ret_c1
    bMs[27] = b_ret_c1
    mMl[27] = m_ret_c2
    bMl[27] = b_ret_c2
    print ("Custom coefficients used")
  } else if((!is.na(m_mort_c1) & !is.na(b_mort_c1) & !is.na(m_ret_c1) & !is.na(b_ret_c1) & (is.na(m_mort_c2) & is.na(b_mort_c2) & is.na(m_ret_c2) & is.na(b_ret_c2)) )) {
    err = 27
    mSs[27] = m_mort_c1
    bSs[27] = b_mort_c1
    mMs[27] = m_ret_c1
    bMs[27] = b_ret_c1
    mSl[27] = mSs[27]
    bSl[27] = bSs[27]
    mMl[27] = mMs[27]
    bMl[27] = bMs[27]
    print ("Custom coefficients used")
  } else if(( (!is.na(m_mort_c2) & !is.na(b_mort_c2) & !is.na(m_ret_c2) & !is.na(b_ret_c2) & is.na(m_mort_c1) & is.na(b_mort_c1) & is.na(m_ret_c1) & is.na(b_ret_c1)))) {
    err = 27
    mSl[27] = m_mort_c2
    bSl[27] = b_mort_c2
    mMl[27] = m_ret_c2
    bMl[27] = b_ret_c2
    mSs[27] = mSl[27]
    bSs[27] = bSl[27]
    mMs[27] = mMl[27]
    bMs[27] = bMl[27]
    print ("Custom coefficients used")
  } else if ((is.na(m_mort_c2) & is.na(b_mort_c2) & is.na(m_ret_c2) & is.na(b_ret_c2) & is.na(m_mort_c1) & is.na(b_mort_c1) & is.na(m_ret_c1) & is.na(b_ret_c1))) {
    err=err
    print("prepopulated error coefficients used from case number:")
    print(err)
  } else if((!is.na(m_mort_c2) & !is.na(b_mort_c2) & !is.na(m_ret_c2) & !is.na(b_ret_c2) & !is.na(m_mort_c1) & !is.na(b_mort_c1) & !is.na(m_ret_c1) & !is.na(b_ret_c1) & !is.na(err)) | (!is.na(m_mort_c2) & !is.na(b_mort_c2) & !is.na(m_ret_c2) & !is.na(b_ret_c2) & is.na(m_mort_c1) & is.na(b_mort_c1) & is.na(m_ret_c1) & is.na(b_ret_c1) & !is.na(err)) | (is.na(m_mort_c2) & is.na(b_mort_c2) & is.na(m_ret_c2) & is.na(b_ret_c2) & !is.na(m_mort_c1) & !is.na(b_mort_c1) & !is.na(m_ret_c1) & !is.na(b_ret_c1) & !is.na(err))) {
    err = 27
    mSs[27] = m_mort_c1
    bSs[27] = b_mort_c1
    mSl[27] = m_mort_c2
    bSl[27] = b_mort_c2
    mMs[27] = m_ret_c1
    bMs[27] = b_ret_c1
    mMl[27] = m_ret_c2
    bMl[27] = b_ret_c2
    print("Values for both err and custom coefficients used - this will use custom values for model prediction")
  } else {
    stop("There is not a matching set of customized coefficients. Fill all 8 coefficeint variables or all the 'class1' or 'class2' values. A '1' can be used for b coefficients if previously undefined")
  }

  #3. Run dataset through the equations

    nrow=length(t)
    ncol=c
    TApresum<- numeric(c)
    plc1<- numeric(c)
    TaLl<- numeric(c)
    TaLs<- numeric(c)

    #Survival regression
    YSs <- (mSs[err] * tadj) + bSs[err]       #Eq. 1a - survival rate of class1
    YSs[YSs > 1] <- 1
    YSl <- (mSl[err] * tadj) + bSl[err]       #Eq. 1b - survival rate of class2 individuals
    YSl[YSl > 1] <- 1

    #Tag retention regression
    YMs <- (mMs[err] * tadj) + bMs[err]        #Eq. 2a - tag retention rate of class1
    YMs[YMs > 1] <- 1
    YMl <- (mMl[err] * tadj) + bMl[err]        #Eq. 2b - tag retention rate of class2 individuals
    YMl[YMl > 1] <- 1

    plc1 <- n_c1 / nT                #Eq. 3 - proportion of class1 individuals
    TaLs <- nT * plc1               #Eq. 4a - number of class2 tagged individual at large
    TaLl <- nT - (nT * plc1)         #Eq. 4b - number of class1 individuals at large

    TSs <- matrix(NA,nrow=ltadj+c+1,ncol=ncol)        #Eq. 5a - living number of class1 tagged individuals at large
    for (q in 1:c) {
     # first value in column
      TSs[q,q]<- TaLs[q]
      for (i in 1:(ltadj)) {
        # each following value
        TSs[(q+i),q] <- TSs[q,q] * YSs[i]
      }}

    TSl <- matrix(NA,nrow=ltadj+c+1,ncol=ncol)        #Eq. 5b - living number of class2 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TSl[q,q]<- TaLl[q]
      for (i in 1:(ltadj)) {
        # each following value
        TSl[(q+i),q] <- TSl[q,q] * YSl[i]
      }}

    TAs <- matrix(NA,nrow=ltadj+c+1,ncol=ncol)        #Eq. 6a - adjusted number of class1 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TAs[q,q]<- TSs[q,q]
      for (i in 1:(ltadj)) {
        # each following value
        TAs[(q+i),q] <- TSs[q,q] * YMs[i]
      }}

    TAl <- matrix(NA,nrow=ltadj+c+1,ncol=ncol)       #Eq. 6b - adjusted number of class2 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TAl[q,q]<- TSl[q,q]
      for (i in 1:(ltadj)) {
        # each following value
        TAl[(q+i),q] <- TSl[q,q] * YMl[i]
      }}

    TA <- TAs + TAl                               #Eq. 7 - total adjusted number of tagged individuals at large
    TA <- TA[1:(c), ]

    TAsum <- rowSums(TA, na.rm=T)           #Eq. 8 - weekly sum of adjusted number of tags at large

    TDF <- 1- (TAsum /TaL)                 #Eq. 9 - Tag depreciation factor

    RA <-  (R+0.001) * (1+TDF)       #Eq. 10 - Adjusted number of recaptures
    RA[1]<- 0

    # save all outputs
    datacomp <- data.frame(week = seq_along(TaL), nT=nT, n_c1=n_c1, TaL=TaL, TAsum=TAsum, TDF=TDF, YSs=YSs, YSl=YSl, YMs=YMs, YMl=YMl,TaLs=TaLs, TaLl=TaLl, R=R, Rpercent=R/TaL, RA=RA, PSE=100*((abs(R-RA))/RA))

    return(datacomp) #Produce a dataframe of primary calculations.

}
