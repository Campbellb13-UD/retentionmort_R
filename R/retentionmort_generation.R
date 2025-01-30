#' Estimating user-based tagging mortality and tag shedding error using
#' artificial mark-recapture data
#'
#' This model estimates the percent loss in tagged animals at large for
#' field-based recapture studies based on a linear decrease in survival and tag
#' retention (including lost tags and missidentified tags) that gets projected
#' for five weeks per tagging cohort based on laboratory retention/survival
#' studies. This `retentionmort_generation` differs from `retentionmort` because
#' it generates a mark-recapture dataset instead of relying on field data,
#' making it possible to estimate the expected error associated with an upcoming
#' field effort to provide insight on methods development. The model is changed
#' by linear regression coefficients of weekly tag loss rate, weekly mortality
#' rate, and their respective intercepts. The coefficients used can be selected
#' from the currently included list using the `err` input or be customized. This
#' function is also capable of working with a cofactor with two conditions
#' (e.g. class1 individuals and large individuals) to improve resolution for
#' more specified studies.
#'
#' @param n           The number of iterations (i.e. generated mark-recapture
#'                    datasets) the model will run through (default = 100).
#'
#' @param min_weeks   The minimum number of efforts that each generated
#'                    mark-recapture dataset will operate for (default = 6),
#'                    must  be at least 6.
#'
#' @param max_weeks   The maximum number of efforts that each generated
#'                    mark-recapture dataset will run for (default = 100).
#'
#' @param max_tags    The maximum number of individuals that will be tagged per
#'                    effort (default = 500).
#'
#' @param prop_class1 (Optional) The estimated proportion of tagged individuals
#'                    that belong to the first of two classifications (use 1 or
#'                    0 if none; default = 0).
#'
#' @param max_recap   The maximum proportion of recaptured individuals per
#'                    effort (default = 0.5).
#'
#' @param err         A value (between 1 and 26) that represents the weekly
#'                    mortality rate and weekly tag loss rate from a preloaded
#'                    case study listed in the metadata (default = 2).
#'                    Alternatively, model coefficients can be manually included
#'                    using a combination of the preceding parameters. While the
#'                    preloaded data are based on weekly time stamps, customized
#'                    model coefficients can reflect any time period specified
#'                    and the projection will predict loss at 5 times the time
#'                    interval.
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
#'             week =       The week in the study
#'             c =          The total number of efforts within the dataset
#'             iteration =  The iteration number of the dataset
#'             nT =         The number of tagged individuals per tagging effort
#'             nlc1 =       The number of tagged individuals in class1 per
#'                          tagging effort
#'             TaL =        The cumulative number of tagged individuals at large
#'             TAsum =      The weekly sum of adjusted number of tags at large
#'             TDF =        Tag depreciation factor
#'             YSs =        Resultant survival rate of class1
#'             YSl =        Resultant survival rate of class2
#'             YMs =        Resultant tag loss rate of class1
#'             YMl =        Resultant tag loss rate of class2
#'             TaLs =       The cumulative number of class1 individuals tagged
#'                          at large
#'             TaLl =       The cumulative number of class2 individuals tagged
#'                          at large
#'             R =          The number of recaptured individuals per effort
#'             Rpercent =   The proportion of recaptured individuals
#'             RA =         The adjusted number of recaptured individuals
#'             PSE =        The percent standard error between observed and
#'                          estimated recaptured individuals
#'
#' @examples
#'
#' #Using only default variables
#'   datacomp = retentionmort_generation()
#'
#'
#' #Using custom model parameters for one class type
#'   datacomp = retentionmort_generation(n = 100, max_weeks = 100, prop_class1 =
#'                         0, max_recap = 0.5, err = NA, m_mort_c1=-0.0625,
#'                         b_mort_c1=1.06, m_ret_c1=-0.113, b_ret_c1=1.05
#'                         )
#'            #or
#'   datacomp = retentionmort(n = 100, max_weeks = 100, prop_class1 = 0,
#'                         max_recap = 0.5, err = NA, m_mort_c2=-0.0203,
#'                         b_mort_c2=1.03, m_ret_c2=-0.0541, b_ret_c2=0.993
#'                         )
#'
#'
#' #Using custom model parameters for two class types
#'   datacomp = retentionmort(n = 100, max_weeks = 100, prop_class1 =
#'                         0, max_recap = 0.5, err = NA, m_mort_c1=-0.0625,
#'                         b_mort_c1=1.06, m_ret_c1=-0.113, b_ret_c1=1.05,
#'                         m_mort_c2=-0.0203, b_mort_c2=1.03, m_ret_c2=-0.0541,
#'                         b_ret_c2=0.993
#'                         )
#'           #or
#'   datacomp = retentionmort(n = 100, max_weeks = 100, prop_class1 =
#'                         0, max_recap = 0.5, err = 2, m_mort_c1=-0.0625,
#'                         b_mort_c1=1.06, m_ret_c1=-0.113, b_ret_c1=1.05,
#'                         m_mort_c2=-0.0203, b_mort_c2=1.03, m_ret_c2=-0.0541,
#'                         b_ret_c2=0.993
#'                         ) #err gets overrided by customized coefficients
#'
#'
#' @export
retentionmort_generation = function(n = 100,         # number of datasets generated and ran through the model
                                    min_weeks = 6,   # minimum number of weeks (default 6), must be at least 6
                                    max_weeks = 100, # maximum number of weeks (default 100)
                                    max_tags = 500 , # max number of animals tagged per effort (default 500)
                                    prop_class1 = 0, # estimated proportion of individuals in the first class (default 0)
                                    max_recap = 0.5, # max proportion of recaptured individuals per cohort (default 0.5)
                                    err = 2,         # choose error values to use from metadata list (default McCutcheon et al. - average case)
                                    m_mort_c1 = NA,  # slope of weekly mortality rate for class1 (default NA)
                                    b_mort_c1 = NA,  # intercept of weekly mortality rate for class1 (default NA)
                                    m_ret_c1 = NA,   # slope of weekly tag loss rate for class1 (default NA)
                                    b_ret_c1 = NA,   # intercept of weekly tag loss rate for class1 (default NA)
                                    m_mort_c2 = NA,  # slope of weekly mortality rate for class2 (default NA)
                                    b_mort_c2 = NA,  # intercept of weekly mortality rate for class2 (default NA)
                                    m_ret_c2 = NA,   # slope of weekly tag loss rate for class2 (default NA)
                                    b_ret_c2 = NA    # intercept of weekly tag loss rate for class1 (default NA)
){

  #1. Establish known coefficients based on laboratory studies - custom values can be added to the end of the mSs, bSs, mSl, bSl, mMs, bMs, mMl, and bMl vectors based on weekly tag loss rate and mortality rate.
  mSs <- c(-0.313,-0.0625,0, -0.001575,-0.014799, -0.029598, -0.009259, -0.01925, -0.027417, -0.01925, -0.030917, -0.1, -0.125, -0.0105, -0.014, -0.005133, -0.036867, -0.042, -0.004375,-0.01, -0.06,0,0,-0.012895,0,-0.00000001)    #slope of survival - class 1
  mSs<- matrix(mSs, ncol=1)
  bSs <- c(1.31, 1.06, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)       #intercept of survival - class 1
  bSs<- matrix(bSs, ncol=1)
  mSl <- c(-0.0697, -0.0203, 0, -0.001575, -0.014799, -0.029598, -0.009259, -0.01925, -0.027417, -0.01925, -0.030917, -0.1, -0.125, -0.0105, -0.014, -0.005133, -0.036867, -0.042, -0.004375, -0.01, -0.06,0,0,-0.012895,0,-0.000000001) #slope of survival - class 2
  mSl<- matrix(mSl, ncol=1)
  bSl <- c(1.009, 1.03, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)      #intercept of survival - class 2
  bSl<- matrix(bSl, ncol=1)
  mMs <- c(-0.245, -0.113, 0, -0.000875, -0.016279, -0.036176, -0.010101, -0.007, -0.029167, -0.004667, -0.025083, -0.0125, 0, -0.013067, -0.046083, -0.003967, -0.016333, -0.021, -0.004375,-0.01, 0, -0.089362,-0.0133,-0.000921,0,-0.001167)   #slope of tag retention - class 1
  mMs<- matrix(mMs, ncol=1)
  bMs <- c(0.906, 1.05, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)      #intercept of tag retention - class 1
  bMs<- matrix(bMs, ncol=1)
  mMl <- c(-0.0118, -0.0541, 0, -0.000875, -0.016279, -0.036176, -0.010101, -0.007, -0.029167, -0.004667, -0.025083, -0.0125, 0, -0.013067, -0.046083, -0.003967, -0.016333, -0.021, -0.004375, -0.01, 0, -0.089362,-0.0133,-0.000921,0,-0.001167) #slope of tag retention - class 2
  mMl<- matrix(mMl, ncol=1)
  bMl <- c(0.716, 0.993, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)     #intercept of tag retention - class 2
  bMl<- matrix(bMl, ncol=1)

  datacomp = data.frame() # make empty dataframe to save data into

  #2. Add errors to check inputs
  # check prop_class1
  if((max(prop_class1) > 1) | (min(prop_class1) < 0))
    stop("prop_class1 is the proportion of animals in the first group which is multiplied by each tagging effort and must be a value between 0 and 1. If no distinction is used, either a single 1 or 0 will remove this conditon")
  # check max_tags
  if(!is.vector(max_tags))
    stop("max_tags must be a vector of total number of tags per effort")
  # check min_weeks and max_weeks
  if(!is.numeric(min_weeks) | !is.numeric(max_weeks) | !min_weeks >= 6 | !max_weeks >= 6)
    stop("min_weeks and max_weeks must be a number of tagging efforts and must be greater than or equal to 6")
  if(min_weeks > max_weeks)
    stop("min_weeks needs to be less than max_weeks")
  # check max_recap
  if(max_recap > 1 | max_recap < 0)
    stop("max_recap is the proportion of recaptured individuals per effort and must be a value between 0 and 1.")
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
  if(!is.na(m_mort_c1) & !is.na(b_mort_c1) & !is.na(m_mort_c2) & !is.na(b_mort_c2) & !is.na(m_ret_c1) & !is.na(b_ret_c1) & !is.na(m_ret_c2) & !is.na(b_ret_c2) & is.na(err)) {
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

  # 3. Start the loop to generate datasets and run the model through them
  for (k in 1:n) {  #change the range of k to alter the number of randomly generated model runs. If using custom field data, the for loop can be removed

    # 3b. Adjusting parameters to generate a randomized recapture dataset -- these values can be customized depending on subject species and expected field effort.
    c <- sample(min_weeks : max_weeks, 1, replace = TRUE) #total number of tagging efforts - minimum 6
    nT<- sample(seq(0,max_tags, by = 1), c, replace = TRUE) #number of tagged individuals per tagging effort

    TaL<- length(nT) #cumulative number of tagged individuals per effort
    TaL[1] <- nT[1]
    for (p in 1:c) {
      TaL[p+1]<- TaL[p]+nT[p+1]
    }
    TaL<- TaL[-length(TaL)]

    # number of individuals in class1
    nlc1 = nT*prop_class1
    nlc1 = ceiling(nlc1)

    # number of recaptures per effort
    R<- runif(c, min=0*TaL, max=max_recap*TaL)
    R<- ceiling(R) #rounding up to the nearest whole number
    R[1]<-0 #must make the first effort have a recapture rate of 0

    # Set up each timestep
    t<- c(0:5)                       #for the 5 week time step
    tadj<- c(0:5, rep(5, c-6))       #used to force no more regression after 5 weeks
    ltadj<- length(tadj)

    # 3c. Produce regression equations from known coefficients

    #Survival regression
    YSs <- (mSs[err] * tadj) + bSs[err]       #Eq. 1a - survival rate of class 1
    YSs[YSs > 1] <- 1
    YSl <- (mSl[err] * tadj) + bSl[err]       #Eq. 1b - survival rate of class 2
    YSl[YSl > 1] <- 1

    #Tag retention regression
    YMs <- (mMs[err] * tadj) + bMs[err]        #Eq. 2a - tag retention rate of class 1
    YMs[YMs > 1] <- 1
    YMl <- (mMl[err] * tadj) + bMl[err]        #Eq. 2b - tag retention rate of class 2
    YMl[YMl > 1] <- 1

    # 3d. Perform projections using the artificial dataset

    ncol2=c
    plc1<- numeric(c)
    TaLl<- numeric(c)
    TaLs<- numeric(c)

    plc1 <- nlc1 / nT                #Eq. 3 - proportion of class1 individuals
    plc1<- matrix(plc1, ncol=1)

    TaLs <- nT * plc1               #Eq. 4a - number of class1 tagged individual at large
    TaLl <- nT - (nT * plc1)         #Eq. 4b - number of class2 individuals at large

    TSs <- matrix(NA,nrow=ltadj+c+1,ncol=ncol2)        #Eq. 5a - living number of class1 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TSs[q,q]<- TaLs[q]
      for (i in 1:(ltadj)) {
        # each following value
        TSs[(q+i),q] <- TSs[q,q] * YSs[i]
      }}

    TSl <- matrix(NA,nrow=ltadj+c+1,ncol=ncol2)        #Eq. 5b - living number of class2 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TSl[q,q]<- TaLl[q]
      for (i in 1:(ltadj)) {
        # each following value
        TSl[(q+i),q] <- TSl[q,q] * YSl[i]
      }}

    TAs <- matrix(NA,nrow=ltadj+c+1,ncol=ncol2)        #Eq. 6a - adjusted number of class1 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TAs[q,q]<- TSs[q,q]
      for (i in 1:(ltadj)) {
        # each following value
        TAs[(q+i),q] <- TSs[q,q] * YMs[i]
      }}

    TAl <- matrix(NA,nrow=ltadj+c+1,ncol=ncol2)       #Eq. 6b - adjusted number of class2 tagged individuals at large
    for (q in 1:c) {
      # first value in column
      TAl[q,q]<- TSl[q,q]
      for (i in 1:(ltadj)) {
        # each following value
        TAl[(q+i),q] <- TSl[q,q] * YMl[i]
      }}

    TA <- TAs + TAl                               #Eq. 7 - total adjusted number of tagged individuals at large
    TA <- TA[1:(c), ]

    TAsum <- rowSums(TA, na.rm=T)           #Eq. 8 - weekly sum of total adjusted number of tags at large

    TDF <- 1- (TAsum /TaL)                 #Eq. 9 - Tag depreciation factor

    RA <-  (R+0.001) * (1+TDF)       #Eq. 10 - Adjusted number of recaptures
    RA[1]<- 0

    #4. Produce a dataframe of primary calculations.
    df <- data.frame(week = seq_along(TaL), c=rep(c), iteration=rep(k), nT=nT, nlc1=nlc1, TaL=TaL, TAsum=TAsum, TDF=TDF, YSs=YSs, YSl=YSl, YMs=YMs, YMl=YMl,TaLs=TaLs, TaLl=TaLl, R=R, Rpercent=R/TaL, RA=RA, PSE=100*((abs(R-RA))/RA))

    datacomp<- rbind(datacomp,df) #add model run to the existing dataset from previous iterations

    # save certain variables for each iteration in the loop
    rm(list = setdiff(ls(), c("datacomp", "max_weeks", "max_tags", "prop_class1", "max_recap","min_weeks", "err","mSs","bSs","mMs","bMs","mSl","bSl","mMl","bMl" ))
    )
  }
  return(datacomp) #return the dataset
}



