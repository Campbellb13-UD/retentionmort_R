#' Generate a .html markdown file of preliminary figures pertaining from either
#' the `retenionmort` or `retentionmort_generation` function.
#'
#' By inputing the `datacomp` dataframe, this function will save a markdown file
#' in the working directory named `retentionmort.html`that provides helpful
#' information on the error associated with the number of recaptured individuals
#' compared to the expected number provided by the model. Some figures will be
#' less applicable for the field data application using the `retenionmort`
#' function due to a low sample size, specifically figures 5 and 6.
#'
#' @param datacomp The file generated from either the `retentionmort` or
#'                 `retentionmort_generation` functions.
#'
#' @returns
#'
#' This function will return one markdown file named `retentionmort.html` in your
#' current working directory listing some helpful information for analyzing
#' model data generated from the `retentionmort` or `retentionmort_generation`
#' functions
#'
#' @examples
#'
#' #Using `retentionmort_generation` to produce multiple iterations of data to
#' #run the model through
#'    datacomp = retentionmort_generation()
#'    Rmark = retentionmort_figure(datacomp)
#'
#'
#' #Creating a dataset with `test_dataset_retentionmort` and running the
#' #`retentionmort` function
#'    list2env(test_dataset_retentionmort(), envir = .GlobalEnv)
#'    datacomp = retentionmort(n_c1=n_c1, nT=nT, TaL=TaL, c=c, R=R, err=2)
#' #Creating the markdown on datacomp
#'    Rmark = retentionmort_figure(datacomp)
#'
#' @export
retentionmort_figure <- function(datacomp){

  cat("---
title: \"Preliminary analysis from the output from the *retentionmort* or *retentionmort_generation* functions\"
author: \"Brendan Campbell* (bpc@udel.edu), Jasper McCutcheon, Rileigh E. Hudock, Noah Motz, Madison Windsor, Aaron Carlisle, Edward Hale\"
date: \"V1: 2024-10-07\"
output:
  html_document: default
  pdf_document: default
---

## This markdown will generate a series of preliminary figures using the *datacomp* dataframe that was created from either retentionmort function.

Install packages and load the dataset named *datacomp* from either the *retentionmort_generation* or the *retentionmort* function

- dplyr
- ggplot2
- gridExtra
- patchwork
- readr

\`\`\`{r include=FALSE}
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

if (!require(readr)) install.packages('readr')
library(readr)

if (!require(patchwork)) install.packages('patchwork')
library(patchwork)

if (!require(rmarkdown)) install.packages('rmarkdown')
library(rmarkdown)


\`\`\`

## Preliminary data exploration

Table of relevant data outputs
\`\`\`{r, echo = FALSE}
  #make a table to print
data <- matrix(c(
  round(min(datacomp$TaL, na.rm = TRUE), 2), round(max(datacomp$TaL, na.rm = TRUE), 2), round(mean(datacomp$TaL, na.rm = TRUE), 2), round(sd(datacomp$TaL, na.rm = TRUE), 2),
  round(min(datacomp$TAsum, na.rm = TRUE), 2), round(max(datacomp$TAsum, na.rm = TRUE), 2), round(mean(datacomp$TAsum, na.rm = TRUE), 2), round(sd(datacomp$TAsum, na.rm = TRUE), 2),
  round(min(datacomp$R, na.rm = TRUE), 2), round(max(datacomp$R, na.rm = TRUE), 2), round(mean(datacomp$R, na.rm = TRUE), 2), round(sd(datacomp$R, na.rm = TRUE), 2),
  round(min(datacomp$RA, na.rm = TRUE), 2), round(max(datacomp$RA, na.rm = TRUE), 2), round(mean(datacomp$RA, na.rm = TRUE), 2), round(sd(datacomp$RA, na.rm = TRUE), 2),
  round(min(datacomp$PSE, na.rm = TRUE), 2), round(max(datacomp$PSE, na.rm = TRUE), 2), round(mean(datacomp$PSE, na.rm = TRUE), 2), round(sd(datacomp$PSE, na.rm = TRUE), 2)
), nrow=5, byrow=TRUE)
data <- as.table(data)
colnames(data) <- c('Minimum', 'Maximum', 'Average', 'SD')
rownames(data) <- c('Obs. Tags at Large', 'Exp. Tags at Large', 'Obs. Recaptures', 'Exp. Recaptures', 'PSE of Recaptures (%)')
print(data)

#paired t-tests
ttest_recapture<- with(datacomp, t.test(RA,R, paired = T))
ttest_tagsatlarge<- with(datacomp, t.test(TAsum,TaL, paired = T))

\`\`\`

Paired T-Test of expected vs. observed tags at large
\`\`\`{r, echo = TRUE}
print(ttest_tagsatlarge)
\`\`\`

Paired T-Test of expected vs. observed recaptures
\`\`\`{r, echo = TRUE}
print(ttest_recapture)
\`\`\`


## Exploratory Figures


1. Histogram of percent standard error (PSE) between estimated and observed number of tags at large

\`\`\`{r, echo = FALSE}
hist(datacomp$PSE)
\`\`\`

2. Percent standard error as a function of the number of observed tags at large. A red line is drawn at the 20% standard error, orange at 15%, yellow at 10%, and green at 5% to denote potential error thresholds for data.

\`\`\`{r, echo = FALSE}
ggplot(datacomp, aes(x = TaL, y = PSE)) +
  geom_point() +
  geom_hline(yintercept = 20, color = \"#D55E00\") +
geom_hline(yintercept = 15, color = \"#E69F00\") +
  geom_hline(yintercept = 10, color = \"#F0E442\") +
  geom_hline(yintercept = 5, color = \"#009E73\") +
  ggtitle(\"Error threshold for model runs compared to Total number of tags at large\")+
  ylab(\"Percent Standard Error\")+
  xlab(\"Tags at Large\")+
  theme_bw()
\`\`\`

3. Deviation between predicted and observed number of recaptures. Blue line indicates a perfect 1:1 match while green indicates 5% deviation, yellow 10%, orange 15%, red 20%.

\`\`\`{r, echo = FALSE}
ggplot(datacomp, aes(x = R, y = RA)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, color = \"#56B4E9\") +
geom_abline(slope=1.2, intercept=0, color = \"#D55E00\")+
  geom_abline(slope=0.8, intercept=0, color = \"#D55E00\")+
  geom_abline(slope=1.15, intercept=0, color =\"#E69F00\")+
  geom_abline(slope=0.85, intercept=0, color = \"#E69F00\")+
  geom_abline(slope=1.10, intercept=0, color =\"#F0E442\")+
  geom_abline(slope=0.9, intercept=0, color = \"#F0E442\")+
  geom_abline(slope=1.05, intercept=0, color =\"#009E73\")+
  geom_abline(slope=0.95, intercept=0, color = \"#009E73\")+
  ggtitle(\"Deviation in tag recovery across number of tagged animals at large\")+
  ylab(\"Expected number of recaptures after accounting for mortality and tag loss\")+
  xlab(\"Observed number of recaptures\")+
  theme_bw()
\`\`\`

4. Percent standard error between expected and observed tags at large compared to the number of tagging efforts conducted. Red line indicates 20% PSE, orange 15%, yellow 10%, and green 5% PSE.

\`\`\`{r, echo = FALSE}
ggplot(datacomp, aes(x = week, y = PSE)) +
  geom_point() +
  geom_hline(yintercept = 20, color = \"#D55E00\") +
geom_hline(yintercept = 15, color = \"#E69F00\") +
  geom_hline(yintercept = 10, color = \"#F0E442\") +
  geom_hline(yintercept = 5, color = \"#009E73\") +
  ggtitle(\"Amount of error associated with number of sampling efforts\")+
  ylab(\"Percent Standard Error\")+
  xlab(\"Number of tagging efforts\")+
  theme_bw()
\`\`\`

5. Determine the critical points in the number of tags at large from each model run that results in percent standard errors exceeding 5%, 10%, 15%, and 20% thresholds. A blank figure denotes no critical points that exceed the PSE threshold. This plot will likely only work for generated data using the *retentionmort_generation* function. The table below summarizes the data represented for each error threshold.

\`\`\`{r, echo = FALSE}
#filter out error thresholds
df.20 <- datacomp %>%
  filter(PSE >= 19.5 & PSE <= 20.5)
df.15 <- datacomp %>%
  filter(PSE>= 14.5 & PSE <= 15.5)
df.10 <- datacomp %>%
  filter(PSE>= 9.5 & PSE <= 10.5)
df.5 <- datacomp %>%
  filter(PSE>= 4.5 & PSE <= 5.5)

#Tags at Large
if (nrow(df.20) > 2) {
bptal20<- ggplot(df.20, aes(x = factor(1), y = TaL)) +
  geom_violin(fill = \"white\", color = \"black\") +
  geom_boxplot(width = 0.25, notch = TRUE, fill = \"#D55E00\", color = \"black\") +
xlab(\"\") +
  ylab(\"Range of tags at large to exceed a 20% error threshold\") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
} else {bptal20<- plot_spacer()}

if (nrow(df.15) > 2) {
  bptal15<- ggplot(df.15, aes(x = factor(1), y = TaL)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#E69F00\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of tags at large to exceed a 15% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bptal15<- plot_spacer()}

if (nrow(df.10) > 2) {
  bptal10<- ggplot(df.10, aes(x = factor(1), y = TaL)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#F0E442\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of tags at large to exceed a 10% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bptal10<- plot_spacer()}

if (nrow(df.5) > 2) {
  bptal5<- ggplot(df.5, aes(x = factor(1), y = TaL)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#009E73\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of tags at large to exceed a 5% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bptal5<- plot_spacer()}

grid.arrange(bptal20, bptal15, bptal10, bptal5, nrow = 1)

dataTALthresh <- matrix(c(
  round(min(df.5$TaL, na.rm = TRUE), 2), round(max(df.5$TaL, na.rm = TRUE)), round(median(df.5$TaL, na.rm = TRUE), 2), round(mean(df.5$TaL, na.rm = TRUE), 2), round(sd(df.5$TaL, na.rm = TRUE), 2), as.numeric(length(df.5$TaL)),
  round(min(df.10$TaL, na.rm = TRUE), 2), round(max(df.10$TaL, na.rm = TRUE)), round(median(df.10$TaL, na.rm = TRUE), 2), round(mean(df.10$TaL, na.rm = TRUE), 2), round(sd(df.10$TaL, na.rm = TRUE), 2),as.numeric(length(df.10$TaL)),
  round(min(df.15$TaL, na.rm = TRUE), 2), round(max(df.15$TaL, na.rm = TRUE)), round(median(df.15$TaL, na.rm = TRUE), 2), round(mean(df.15$TaL, na.rm = TRUE), 2), round(sd(df.15$TaL, na.rm = TRUE), 2),as.numeric(length(df.15$TaL)),
  round(min(df.20$TaL, na.rm = TRUE), 2), round(max(df.20$TaL, na.rm = TRUE)), round(median(df.20$TaL, na.rm = TRUE), 2), round(mean(df.20$TaL, na.rm = TRUE), 2), round(sd(df.20$TaL, na.rm = TRUE), 2), as.numeric(length(df.20$TaL))
), nrow=4, byrow=TRUE)
dataTALthresh <- as.table(dataTALthresh)
colnames(dataTALthresh) <- c('Minimum', 'Maximum', 'Median', 'Average', 'SD', 'n')
rownames(dataTALthresh) <- c('4.5 - 5.5% PSE','9.5 - 10.5% PSE', '14.5 - 15.5% PSE', '19.5 - 20.5% PSE')
print(dataTALthresh)

```

6. Determine the critical points in recapture efforts from each model run that results in percent standard errors exceeding 5%, 10%, 15%, and 20% thresholds. A blank figure denotes no critical points that exceed the PSE threshold. This plot will likely only work for generated data using the *retentionmort_generation* function. The table below summarizes the data represented for each error threshold.

```{r, echo = FALSE}
df.20 <- datacomp %>%
  filter(PSE >= 19.5 & PSE <= 20.5)
df.15 <- datacomp %>%
  filter(PSE>= 14.5 & PSE <= 15.5)
df.10 <- datacomp %>%
  filter(PSE>= 9.5 & PSE <= 10.5)
df.5 <- datacomp %>%
  filter(PSE>= 4.5 & PSE <= 5.5)

if (nrow(df.20) > 2) {
bprt20<- ggplot(df.20, aes(x = factor(1), y = R)) +
  geom_violin(fill = \"white\", color = \"black\") +
  geom_boxplot(width = 0.25, notch = TRUE, fill = \"#D55E00\", color = \"black\") +
xlab(\"\") +
  ylab(\"Range of recaptures to exceed a 20% error threshold\") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
} else {bprt20<- plot_spacer()}

if (nrow(df.15) > 2) {
  bprt15<- ggplot(df.15, aes(x = factor(1), y = R)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#E69F00\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of recaptures to exceed a 15% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bprt15<- plot_spacer()}

if (nrow(df.10) > 2) {
  bprt10<- ggplot(df.10, aes(x = factor(1), y = R)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#F0E442\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of recaptures to exceed a 10% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bprt10<- plot_spacer()}

if (nrow(df.5) > 2) {
  bprt5<- ggplot(df.5, aes(x = factor(1), y = R)) +
    geom_violin(fill = \"white\", color = \"black\") +
    geom_boxplot(width = 0.25, notch = TRUE, fill = \"#009E73\", color = \"black\") +
    xlab(\"\") +
    ylab(\"Range of recaptures to exceed a 5% error threshold\") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
} else {bprt5<- plot_spacer()}

grid.arrange(bprt20, bprt15, bprt10, bprt5, nrow = 1)


datarecapthresh <- matrix(c(
  round(min(df.5$R, na.rm = TRUE), 2), round(max(df.5$R, na.rm =TRUE), 2), round(median(df.5$R, na.rm = TRUE), 2), round(mean(df.5$R, na.rm = TRUE), 2), round(sd(df.5$R, na.rm = TRUE), 2), as.numeric(length(df.5$R)),
  round(min(df.10$R, na.rm = TRUE), 2),  round(max(df.10$R, na.rm =TRUE), 2), round(median(df.10$R, na.rm = TRUE), 2), round(mean(df.10$R, na.rm = TRUE), 2), round(sd(df.10$R, na.rm = TRUE), 2),as.numeric(length(df.10$R)),
  round(min(df.15$R, na.rm = TRUE), 2),  round(max(df.15$R, na.rm =TRUE), 2), round(median(df.15$R, na.rm = TRUE), 2), round(mean(df.15$R, na.rm = TRUE), 2), round(sd(df.15$R, na.rm = TRUE), 2),as.numeric(length(df.15$R)),
  round(min(df.20$R, na.rm = TRUE), 2),  round(max(df.20$R, na.rm =TRUE), 2), round(median(df.20$R, na.rm = TRUE), 2), round(mean(df.20$R, na.rm = TRUE), 2), round(sd(df.20$R, na.rm = TRUE), 2), as.numeric(length(df.20$R))
), nrow=4, byrow=TRUE)
datarecapthresh <- as.table(datarecapthresh)
colnames(datarecapthresh) <- c('Minimum', 'Maximum','Median', 'Average', 'SD', 'n')
rownames(datarecapthresh) <- c('4.5 - 5.5% PSE','9.5 - 10.5% PSE', '14.5 - 15.5% PSE', '19.5 - 20.5% PSE')
print(datarecapthresh)
```



",file="retentionmort.Rmd")
  rmarkdown::render("retentionmort.Rmd")
}

