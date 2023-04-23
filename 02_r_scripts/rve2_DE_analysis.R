# code for DE analysis of rve2-2/Col-0 RNA-seq data

library(readxl)
library(dplyr)
library(zoo)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(waffle)
library(readr)
library(RColorBrewer)
library(patchwork)
library(magrittr)
library(ggpubr)
library(UpSetR)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(GOplot)
library(org.At.tair.db)

# 1. RNAseq dataset----
# Temperature & time series RNA-seq experiment
# read-in dataset
transcript_RNAseq <- read_excel('./00_raw_data/Suppl Dataset S4.xlsx', sheet = 2)

# rename 1st column
colnames(transcript_RNAseq)[1] ="Isoform"

# for the transient-cooling day (day 2) calculate the mean for the time-points
# means are as new columns
transcript_RNAseq_means <- transcript_RNAseq %>%
  mutate(Col_T9 = rowMeans(cbind(Col.T9.rep1, Col.T9.rep2, Col.T9.rep3), na.rm = T),
         Col_T10 = rowMeans(cbind(Col.T10.rep1, Col.T10.rep2, Col.T10.rep3), na.rm = T),
         Col_T11 = rowMeans(cbind(Col.T11.rep1, Col.T11.rep2, Col.T11.rep3), na.rm = T),
         Col_T12 = rowMeans(cbind(Col.T12.rep1, Col.T12.rep2, Col.T12.rep3), na.rm = T),
         Col_T13 = rowMeans(cbind(Col.T13.rep1, Col.T13.rep2, Col.T13.rep3), na.rm = T),
         Col_T14 = rowMeans(cbind(Col.T14.rep1, Col.T14.rep2, Col.T14.rep3), na.rm = T),
         Col_T15 = rowMeans(cbind(Col.T15.rep1, Col.T15.rep2, Col.T15.rep3), na.rm = T),
         Col_T16 = rowMeans(cbind(Col.T16.rep1, Col.T16.rep2, Col.T16.rep3), na.rm = T),
         Col_T17 = rowMeans(cbind(Col.T17.rep1, Col.T17.rep2, Col.T17.rep3), na.rm = T),
         Col_T18 = rowMeans(cbind(Col.T18.rep1, Col.T18.rep2, Col.T18.rep3), na.rm = T),
         rve2_T9 = rowMeans(cbind(Rve2.T9.rep1, Rve2.T9.rep2, Rve2.T9.rep3), na.rm = T),
         rve2_T10 = rowMeans(cbind(Rve2.T10.rep1, Rve2.T10.rep2, Rve2.T10.rep3), na.rm = T),
         rve2_T11 = rowMeans(cbind(Rve2.T11.rep1, Rve2.T11.rep2, Rve2.T11.rep3), na.rm = T),
         rve2_T12 = rowMeans(cbind(Rve2.T12.rep1, Rve2.T12.rep2, Rve2.T12.rep3), na.rm = T),
         rve2_T13 = rowMeans(cbind(Rve2.T13.rep1, Rve2.T13.rep2, Rve2.T13.rep3), na.rm = T),
         rve2_T14 = rowMeans(cbind(Rve2.T14.rep1, Rve2.T14.rep2, Rve2.T14.rep3), na.rm = T),
         rve2_T15 = rowMeans(cbind(Rve2.T15.rep1, Rve2.T15.rep2, Rve2.T15.rep3), na.rm = T),
         rve2_T16 = rowMeans(cbind(Rve2.T16.rep1, Rve2.T16.rep2, Rve2.T16.rep3), na.rm = T),
         rve2_T17 = rowMeans(cbind(Rve2.T17.rep1, Rve2.T17.rep2, Rve2.T17.rep3), na.rm = T),
         rve2_T18 = rowMeans(cbind(Rve2.T18.rep1, Rve2.T18.rep2, Rve2.T18.rep3), na.rm = T)) 

# just select the Isoform and means columns
transcript_RNAseq_means_select <- transcript_RNAseq_means %>% 
  dplyr::select(Isoform, Col_T9:Col_T18, rve2_T9:rve2_T18) 

# filter remove rows where the sum of all the day two time-points (rve2-2 and Col-0 samples) are less than TPM 80
transcript_RNAseq_means_select_filtered_low <- transcript_RNAseq_means_select %>% 
  mutate(sum_rows = rowSums(across(where(is.numeric)))) %>% 
  filter(sum_rows >= 80)

# *1.1 time points----
# function to get time-point data for rve2-2 and Col-0

get_hour_simple <- function(df, var1, var2){
  
  hour <- df 
  
  hour_select <- df %>% dplyr::select(Isoform, {{var1}}, {{var2}}) 
  
  return(hour_select)
  
}

h0_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T9, rve2_T9)
h3_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T10, rve2_T10)
h6_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T11, rve2_T11)
h7_5_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T12, rve2_T12)
h9_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T13, rve2_T13)
h12_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T14, rve2_T14)
h15_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T15, rve2_T15)
h18_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T16, rve2_T16)
h21_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T17, rve2_T17)
h24_simple <- get_hour_simple(transcript_RNAseq_means_select_filtered_low, Col_T18, rve2_T18)

# *1.2 DE for time points----

get_consistent_DE <- function(df_time_point, new_col_str, numerator, denominator) {
  
  DE_isoforms <- df_time_point %>%
    mutate({{new_col_str}} := {{numerator}}/{{denominator}}) %>% 
    dplyr::select(Isoform, {{new_col_str}})
  
}

# *1.2.1 repressed----

h0_simple_repressed <- get_consistent_DE(h0_simple, h0_repressed, rve2_T9, Col_T9)
h3_simple_repressed <- get_consistent_DE(h3_simple, h3_repressed, rve2_T10, Col_T10)
h6_simple_repressed <- get_consistent_DE(h6_simple, h6_repressed, rve2_T11, Col_T11)
h7_5_simple_repressed <- get_consistent_DE(h7_5_simple, h7_5_repressed, rve2_T12, Col_T12)
h9_simple_repressed <- get_consistent_DE(h9_simple, h9_repressed, rve2_T13, Col_T13)
h12_simple_repressed <- get_consistent_DE(h12_simple, h12_repressed, rve2_T14, Col_T14)
h15_simple_repressed <- get_consistent_DE(h15_simple, h15_repressed, rve2_T15, Col_T15)
h18_simple_repressed <- get_consistent_DE(h18_simple, h18_repressed, rve2_T16, Col_T16)
h21_simple_repressed <- get_consistent_DE(h21_simple, h21_repressed, rve2_T17, Col_T17)
h24_simple_repressed <- get_consistent_DE(h24_simple, h24_repressed, rve2_T18, Col_T18)

# h0_simple_repressed <- h0_simple %>% 
#   mutate(h0_repressed = rve2_T9/Col_T9) %>% 
#   dplyr::select(Isoform, h0_repressed)

# bind repressed together
repressed_simple <- bind_cols(h0_simple_repressed,
                              h3_simple_repressed,
                              h6_simple_repressed, 
                              h7_5_simple_repressed,
                              h9_simple_repressed,
                              h12_simple_repressed,
                              h15_simple_repressed,
                              h18_simple_repressed,
                              h21_simple_repressed,
                              h24_simple_repressed) %>% 
  dplyr::select(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

RVE1_repressed_simple <- repressed_simple %>% filter(Isoform %in% 'AT5G17300_P1')
RVE7_repressed_simple <- repressed_simple %>% filter(Isoform %in% 'AT1G18330_P1')
CBF1_repressed_simple <- repressed_simple %>% filter(Isoform %in% 'AT4G25490.1')

# rename first column
colnames(repressed_simple)[1] ="Isoform"

# *1.2.2 activated----

h0_simple_activated <- get_consistent_DE(h0_simple, h0_activated, Col_T9, rve2_T9)
h3_simple_activated <- get_consistent_DE(h3_simple, h3_activated, Col_T10, rve2_T10)
h6_simple_activated <- get_consistent_DE(h6_simple, h6_activated, Col_T11, rve2_T11)
h7_5_simple_activated <- get_consistent_DE(h7_5_simple, h7_5_activated, Col_T12, rve2_T12)
h9_simple_activated <- get_consistent_DE(h9_simple, h9_activated, Col_T13, rve2_T13)
h12_simple_activated <- get_consistent_DE(h12_simple, h12_activated, Col_T14, rve2_T14)
h15_simple_activated <- get_consistent_DE(h15_simple, h15_activated, Col_T15, rve2_T15)
h18_simple_activated <- get_consistent_DE(h18_simple, h18_activated, Col_T16, rve2_T16)
h21_simple_activated <- get_consistent_DE(h21_simple, h21_activated, Col_T17, rve2_T17)
h24_simple_activated <- get_consistent_DE(h24_simple, h24_activated, Col_T18, rve2_T18)

# bind activated together
activated_simple <- bind_cols(h0_simple_activated,
                              h3_simple_activated,
                              h6_simple_activated, 
                              h7_5_simple_activated,
                              h9_simple_activated,
                              h12_simple_activated,
                              h15_simple_activated,
                              h18_simple_activated,
                              h21_simple_activated,
                              h24_simple_activated) %>% 
  dplyr::select(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

# rename first column
colnames(activated_simple)[1] ="Isoform"

# 1.3 consistent DE----
# *1.3.1 consistent repressed----
# CBF1 <- repressed_simple_flagged %>% 
#   filter(Isoform %in% 'AT4G25490.1')
repressed_simple_flagged <- repressed_simple %>%
  mutate(consistent_h0_to_h6 = case_when(h0_repressed > 1.2 & h3_repressed > 1.2 & h6_repressed > 1.2 ~ 1, TRUE ~ 0),
         consistent_h3_to_h9 = case_when((h3_repressed > 1.2 | h6_repressed > 1.2) & h7_5_repressed > 1.2 & h9_repressed > 1.2 ~ 1, TRUE ~ 0),
         consistent_h6_to_h12 = case_when((h6_repressed > 1.2 | h7_5_repressed > 1.2) & h9_repressed > 1.2 & h12_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h9_to_h15 = case_when(h9_repressed > 1.2 & h12_repressed > 1.2 & h15_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h12_to_h18 = case_when(h12_repressed > 1.2 & h15_repressed > 1.2 & h18_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h15_to_h21 = case_when(h15_repressed > 1.2 & h18_repressed > 1.2 & h21_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h18_to_h24 = case_when(h18_repressed > 1.2 & h21_repressed > 1.2 & h24_repressed > 1.2  ~ 1, TRUE ~ 0))

repressed_simple_flagged_new <- repressed_simple %>%
  mutate(consistent_h0_to_h6 = case_when(h0_repressed > 1.2 & h3_repressed > 1.2 & h6_repressed > 1.2 ~ 1, TRUE ~ 0),
         consistent_h6_to_h12 = case_when((h6_repressed > 1.2 | h7_5_repressed > 1.2) & h9_repressed > 1.2 & h12_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h12_to_h18 = case_when(h12_repressed > 1.2 & h15_repressed > 1.2 & h18_repressed > 1.2  ~ 1, TRUE ~ 0),
         consistent_h18_to_h24 = case_when(h18_repressed > 1.2 & h21_repressed > 1.2 & h24_repressed > 1.2  ~ 1, TRUE ~ 0))

CBF1_repressed_simple_flagged_new <- repressed_simple_flagged_new %>% filter(Isoform %in% 'AT4G25490.1')
CBF2_repressed_simple_flagged_new <- repressed_simple_flagged_new %>% filter(Isoform %in% 'AT4G25470_P1')
CBF3_repressed_simple_flagged_new <- repressed_simple_flagged_new %>% filter(Isoform %in% 'AT4G25480_JC1')

# *1.3.2 consistent activated----

activated_simple_flagged <- activated_simple %>%
  mutate(consistent_h0_to_h6 = case_when(h0_activated > 1.2 & h3_activated > 1.2 & h6_activated > 1.2 ~ 1, TRUE ~ 0),
         consistent_h3_to_h9 = case_when((h3_activated > 1.2 | h6_activated > 1.2) & h7_5_activated > 1.2 & h9_activated > 1.2 ~ 1, TRUE ~ 0),
         consistent_h6_to_h12 = case_when((h6_activated > 1.2 | h7_5_activated > 1.2) & h9_activated > 1.2 & h12_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h9_to_h15 = case_when(h9_activated > 1.2 & h12_activated > 1.2 & h15_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h12_to_h18 = case_when(h12_activated > 1.2 & h15_activated > 1.2 & h18_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h15_to_h21 = case_when(h15_activated > 1.2 & h18_activated > 1.2 & h21_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h18_to_h24 = case_when(h18_activated > 1.2 & h21_activated > 1.2 & h24_activated > 1.2  ~ 1, TRUE ~ 0))

activated_simple_flagged_new <- activated_simple %>%
  mutate(consistent_h0_to_h6 = case_when(h0_activated > 1.2 & h3_activated > 1.2 & h6_activated > 1.2 ~ 1, TRUE ~ 0),
         consistent_h6_to_h12 = case_when((h6_activated > 1.2 | h7_5_activated > 1.2) & h9_activated > 1.2 & h12_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h12_to_h18 = case_when(h12_activated > 1.2 & h15_activated > 1.2 & h18_activated > 1.2  ~ 1, TRUE ~ 0),
         consistent_h18_to_h24 = case_when(h18_activated > 1.2 & h21_activated > 1.2 & h24_activated > 1.2  ~ 1, TRUE ~ 0))

# *1.4 consistent DE ranked----
#ranking of targets in 6hr windows----

get_consistent_windows <- function(df, filter_col, new_col_str1, new_col_str2, new_col_str3, 
                                   col_to_rank1, col_to_rank2, col_to_rank3,
                                   new_col_str4) {
  
  rank_prep <- df %>% 
    dplyr::filter({{filter_col}} == 1) %>% 
    mutate({{new_col_str1}} := dense_rank(dplyr::desc({{col_to_rank1}})),
           {{new_col_str2}} := dense_rank(dplyr::desc({{col_to_rank2}})),
           {{new_col_str3}} := dense_rank(dplyr::desc({{col_to_rank3}})))
  
  avg_rank <- rank_prep %>% 
    mutate({{new_col_str4}} := rowMeans(cbind({{new_col_str1}}, {{new_col_str2}}, {{new_col_str3}}))) %>% 
    arrange({{new_col_str4}})
  
  rank_arrange <- avg_rank %>%
    dplyr::select(Isoform, {{new_col_str4}})
  
  return(rank_arrange)
}

# *1.4.1 consistent repressed ranked----

consistent_h0_to_h6_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h0_to_h6, 
                                                         h0_rank, h3_rank, h6_rank,
                                                         h0_repressed, h3_repressed, h6_repressed,
                                                         mean_rank_h0_to_h6) %>% 
  write_csv('consistent_h0_to_h6_ranked_rep.csv')

consistent_h0_to_h6_ranked_rep %>% filter(Isoform %in% 'AT4G25490.1')

consistent_h3_to_h9_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h3_to_h9, 
                                                         h3_rank, h6_rank, h9_rank,
                                                         h3_repressed, h6_repressed, h9_repressed,
                                                         mean_rank_h3_to_h9) %>% 
  write_csv('consistent_h3_to_h9_ranked_rep.csv')

#consistent_h3_to_h9_ranked_rep %>% filter(Isoform %in% 'AT4G25490.1')

consistent_h6_to_h12_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h6_to_h12, 
                                                          h6_rank, h9_rank, h12_rank,
                                                          h6_repressed, h9_repressed, h12_repressed,
                                                          mean_rank_h6_to_h12) %>% 
  write_csv('consistent_h6_to_h12_ranked_rep.csv')

consistent_h9_to_h15_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h9_to_h15, 
                                                          h9_rank, h12_rank, h15_rank,
                                                          h9_repressed, h12_repressed, h15_repressed,
                                                          mean_rank_h9_to_h15) %>% 
  write_csv('consistent_h9_to_h15_ranked_rep.csv')

#consistent_h9_to_h15_ranked_rep %>% filter(Isoform %in% 'AT4G25490.1')

consistent_h12_to_h18_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h12_to_h18, 
                                                           h12_rank, h15_rank, h18_rank,
                                                           h12_repressed, h15_repressed, h18_repressed,
                                                           mean_rank_h12_to_h18) %>% 
  write_csv('consistent_h12_to_h18_ranked_rep.csv')

consistent_h15_to_h21_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h15_to_h21, 
                                                           h15_rank, h18_rank, h21_rank,
                                                           h15_repressed, h18_repressed, h21_repressed,
                                                           mean_rank_h15_to_h21) %>% 
  write_csv('consistent_h15_to_h21_ranked_rep.csv')

consistent_h18_to_h24_ranked_rep <- get_consistent_windows(repressed_simple_flagged, consistent_h18_to_h24, 
                                                           h18_rank, h21_rank, h24_rank,
                                                           h18_repressed, h21_repressed, h24_repressed,
                                                           mean_rank_h18_to_h24) %>% 
  write_csv('consistent_h18_to_h24_ranked_rep.csv')

#repressed
#h0_h6
# consistent_h0_to_h6_ranked_rep <- repressed_simple_flagged %>% 
#   filter(consistent_h0_to_h6 == 1)%>%
#   mutate(h0_rank = dense_rank(dplyr::desc(h0_repressed)),
#          h3_rank = dense_rank(dplyr::desc(h3_repressed)),
#          h6_rank = dense_rank(dplyr::desc(h6_repressed)))
# 
# consistent_h0_to_h6_ranked_rep <- consistent_h0_to_h6_ranked_rep %>%
#   mutate(mean_rank_h0_to_h6 = rowMeans(cbind(h0_rank, h3_rank, h6_rank))) %>%
#   arrange(mean_rank_h0_to_h6)
# 
# consistent_h0_to_h6_ranked_rep <- consistent_h0_to_h6_ranked_rep %>%
#   dplyr::select(Isoform, mean_rank_h0_to_h6)

# *1.4.2 consistent activated ranked----

consistent_h0_to_h6_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h0_to_h6, 
                                                         h0_rank, h3_rank, h6_rank,
                                                         h0_activated, h3_activated, h6_activated,
                                                         mean_rank_h0_to_h6) %>% 
  write_csv('consistent_h0_to_h6_ranked_act.csv')

consistent_h3_to_h9_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h3_to_h9, 
                                                         h3_rank, h6_rank, h9_rank,
                                                         h3_activated, h6_activated, h9_activated,
                                                         mean_rank_h3_to_h9) %>% 
  write_csv('consistent_h3_to_h9_ranked_act.csv')

consistent_h6_to_h12_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h6_to_h12, 
                                                          h6_rank, h9_rank, h12_rank,
                                                          h6_activated, h9_activated, h12_activated,
                                                          mean_rank_h6_to_h12) %>% 
  write_csv('consistent_h6_to_h12_ranked_act.csv')

consistent_h9_to_h15_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h9_to_h15, 
                                                          h9_rank, h12_rank, h15_rank,
                                                          h9_activated, h12_activated, h15_activated,
                                                          mean_rank_h9_to_h15) %>% 
  write_csv('consistent_h9_to_h15_ranked_act.csv')

consistent_h12_to_h18_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h12_to_h18, 
                                                           h12_rank, h15_rank, h18_rank,
                                                           h12_activated, h15_activated, h18_activated,
                                                           mean_rank_h12_to_h18) %>% 
  write_csv('consistent_h12_to_h18_ranked_act.csv')

consistent_h15_to_h21_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h15_to_h21, 
                                                           h15_rank, h18_rank, h21_rank,
                                                           h15_activated, h18_activated, h21_activated,
                                                           mean_rank_h15_to_h21) %>% 
  write_csv('consistent_h15_to_h21_ranked_act.csv')

consistent_h18_to_h24_ranked_act <- get_consistent_windows(activated_simple_flagged, consistent_h18_to_h24, 
                                                           h18_rank, h21_rank, h24_rank,
                                                           h18_activated, h21_activated, h24_activated,
                                                           mean_rank_h18_to_h24) %>% 
  write_csv('consistent_h18_to_h24_ranked_act.csv')

#ranking of targets in 6hr windows----
#activated
#h0_h6
# consistent_h0_to_h6_ranked_act <- activated_simple_flagged %>% 
#   filter(consistent_h0_to_h6 == 1)%>%
#   mutate(h0_rank = dense_rank(dplyr::desc(h0_activated)),
#          h3_rank = dense_rank(dplyr::desc(h3_activated)),
#          h6_rank = dense_rank(dplyr::desc(h6_activated)))
# 
# consistent_h0_to_h6_ranked_act <- consistent_h0_to_h6_ranked_act %>%
#   mutate(mean_rank_h0_to_h6 = rowMeans(cbind(h0_rank, h3_rank, h6_rank))) %>%
#   arrange(mean_rank_h0_to_h6)
# 
# consistent_h0_to_h6_ranked_act <- consistent_h0_to_h6_ranked_act %>%
#   dplyr::select(Isoform, mean_rank_h0_to_h6)

# *1.5 split the isoform tag----

split_isoforms <- function(df, new_col_str1, new_col_str2, new_col_str3) {
  
  isoform_col <- df %>% 
    dplyr::select(Isoform)
  
  tags <- isoform_col %>% 
    mutate({{new_col_str1}} := str_split(Isoform, "_", simplify = TRUE)[ ,2],
           {{new_col_str2}} := str_split(Isoform, "[.]", simplify = TRUE)[ ,2]) %>%
    unite({{new_col_str3}}, {{new_col_str1}} : {{new_col_str2}}, sep ="")
  
  return(tags)
  
}

h0_to_h6_tags_rep <- split_isoforms(consistent_h0_to_h6_ranked_rep, tag1, tag2, tags)
h3_to_h9_tags_rep <- split_isoforms(consistent_h3_to_h9_ranked_rep, tag1, tag2, tags)
h6_to_h12_tags_rep <- split_isoforms(consistent_h6_to_h12_ranked_rep, tag1, tag2, tags)
h9_to_h15_tags_rep <- split_isoforms(consistent_h9_to_h15_ranked_rep, tag1, tag2, tags)
h12_to_h18_tags_rep <- split_isoforms(consistent_h12_to_h18_ranked_rep, tag1, tag2, tags)
h15_to_h21_tags_rep <- split_isoforms(consistent_h15_to_h21_ranked_rep, tag1, tag2, tags)
h18_to_h24_tags_rep <- split_isoforms(consistent_h18_to_h24_ranked_rep, tag1, tag2, tags)

h0_to_h6_tags_act <- split_isoforms(consistent_h0_to_h6_ranked_act, tag1, tag2, tags)
h3_to_h9_tags_act <- split_isoforms(consistent_h3_to_h9_ranked_act, tag1, tag2, tags)
h6_to_h12_tags_act <- split_isoforms(consistent_h6_to_h12_ranked_act, tag1, tag2, tags)
h9_to_h15_tags_act <- split_isoforms(consistent_h9_to_h15_ranked_act, tag1, tag2, tags)
h12_to_h18_tags_act <- split_isoforms(consistent_h12_to_h18_ranked_act, tag1, tag2, tags)
h15_to_h21_tags_act <- split_isoforms(consistent_h15_to_h21_ranked_act, tag1, tag2, tags)
h18_to_h24_tags_act <- split_isoforms(consistent_h18_to_h24_ranked_act, tag1, tag2, tags)

# *1.6 count the isoform tag----
count_isoforms <- function(df, col_to_count, new_col_str1, new_col_str2, new_col_str3) {
  
  count <- df %>% 
    count({{col_to_count}}) %>% 
    mutate({{new_col_str1}} := round((n/sum(n)*100), 1))
  
  rank_and_group <- count %>%
    mutate({{new_col_str2}} := rank(-percent),
           {{new_col_str3}} := ifelse(percent >= 5, {{col_to_count}}, 'Others'))
  
  return(rank_and_group)
}

# repressed
h0_to_h6_tags_count_rep <- count_isoforms(h0_to_h6_tags_rep, tags, percent, rank, group)
h3_to_h9_tags_count_rep <- count_isoforms(h3_to_h9_tags_rep, tags, percent, rank, group)
h6_to_h12_tags_count_rep <- count_isoforms(h6_to_h12_tags_rep, tags, percent, rank, group)
h9_to_h15_tags_count_rep <- count_isoforms(h9_to_h15_tags_rep, tags, percent, rank, group)
h12_to_h18_tags_count_rep <- count_isoforms(h12_to_h18_tags_rep, tags, percent, rank, group)
h15_to_h21_tags_count_rep <- count_isoforms(h15_to_h21_tags_rep, tags, percent, rank, group)
h18_to_h24_tags_count_rep <- count_isoforms(h18_to_h24_tags_rep, tags, percent, rank, group)

# activated
h0_to_h6_tags_count_act <- count_isoforms(h0_to_h6_tags_act, tags, percent, rank, group)
h3_to_h9_tags_count_act <- count_isoforms(h3_to_h9_tags_act, tags, percent, rank, group)
h6_to_h12_tags_count_act <- count_isoforms(h6_to_h12_tags_act, tags, percent, rank, group)
h9_to_h15_tags_count_act <- count_isoforms(h9_to_h15_tags_act, tags, percent, rank, group)
h12_to_h18_tags_count_act <- count_isoforms(h12_to_h18_tags_act, tags, percent, rank, group)
h15_to_h21_tags_count_act <- count_isoforms(h15_to_h21_tags_act, tags, percent, rank, group)
h18_to_h24_tags_count_act <- count_isoforms(h18_to_h24_tags_act, tags, percent, rank, group)

# Isoforms_consistent_h0_to_h6_ranked_rep <- consistent_h0_to_h6_ranked_rep %>% 
#   dplyr::select(Isoform)
# 
# h0_to_h6_tags_rep <- Isoforms_consistent_h0_to_h6_ranked_rep %>%
#   mutate(tag1 = str_split(Isoform, "_", simplify = TRUE)[ ,2],
#          tag2 = str_split(Isoform, "[.]", simplify = TRUE)[ ,2]) %>% 
#   unite(tag, c("tag1", "tag2"), sep="") 

# h0_to_h6_tags_count_rep <- h0_to_h6_tags_rep %>%
#   count(tags) %>% 
#   mutate(percent = round((n/sum(n)*100), 1))
# 
# h0_to_h6_tags_count_plotting_data_rep <- h0_to_h6_tags_count_rep %>% 
#   mutate(rank = rank(-percent),
#          group = ifelse(percent >= 5, tags, 'Others'))

# *1.7 Waffle plot prep----

prep_waffle_plot <- function(df, new_col_str1) {
  
  summarise <- df %>% 
    mutate(group = factor(group)) %>% 
    mutate(group = fct_relevel(group, c("Others", "1", "P2", "P1"))) %>%
    arrange(group) %>% 
    group_by(group) %>%
    summarise({{new_col_str1}} := sum(n))
  
  return(summarise)
  
}

# repressed
h0_to_h6_rep_waffle_prep <- prep_waffle_plot(h0_to_h6_tags_count_rep, sum_isoforms)
h3_to_h9_rep_waffle_prep <- prep_waffle_plot(h3_to_h9_tags_count_rep, sum_isoforms)
h6_to_h12_rep_waffle_prep <- prep_waffle_plot(h6_to_h12_tags_count_rep, sum_isoforms)
h9_to_h15_rep_waffle_prep <- prep_waffle_plot(h9_to_h15_tags_count_rep, sum_isoforms)
h12_to_h18_rep_waffle_prep <- prep_waffle_plot(h12_to_h18_tags_count_rep, sum_isoforms)
h15_to_h21_rep_waffle_prep <- prep_waffle_plot(h15_to_h21_tags_count_rep, sum_isoforms)
h18_to_h24_rep_waffle_prep <- prep_waffle_plot(h18_to_h24_tags_count_rep, sum_isoforms)

# activated

h0_to_h6_act_waffle_prep <- prep_waffle_plot(h0_to_h6_tags_count_act, sum_isoforms)
h3_to_h9_act_waffle_prep <- prep_waffle_plot(h3_to_h9_tags_count_act, sum_isoforms)
h6_to_h12_act_waffle_prep <- prep_waffle_plot(h6_to_h12_tags_count_act, sum_isoforms)
h9_to_h15_act_waffle_prep <- prep_waffle_plot(h9_to_h15_tags_count_act, sum_isoforms)
h12_to_h18_act_waffle_prep <- prep_waffle_plot(h12_to_h18_tags_count_act, sum_isoforms)
h15_to_h21_act_waffle_prep <- prep_waffle_plot(h15_to_h21_tags_count_act, sum_isoforms)
h18_to_h24_act_waffle_prep <- prep_waffle_plot(h18_to_h24_tags_count_act, sum_isoforms)

# h0_to_h6_tags_count_plotting_data_waffle_rep <- h0_to_h6_tags_count_plotting_data_rep %>% 
#   mutate(group = factor(group)) %>% 
#   mutate(group = fct_relevel(group, c("Others", "1", "P2", "P1"))) %>%
#   arrange(group) %>% 
#   group_by(group) %>%
#   summarise(sum_isoforms = sum(n))

# *1.8 Waffle plot prep----

plot_waffle <- function(df, title) {
  
  plot <- df %>% 
    ggplot(aes(fill = group, values = sum_isoforms)) +
    geom_waffle(n_rows = 4, size = 0.2, colour = "grey30", flip = TRUE, make_proportional = TRUE) +
    scale_fill_manual(name = NULL,
                      values = c('#bf5b17','#beaed4','#386cb0','#7fc97f')) +
    coord_equal() +
    theme_void() +
    guides(fill = guide_legend(reverse = TRUE, title = "Isoform")) +
    ggtitle(glue('{title}')) +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          legend.position = "left",
          legend.key.size = unit(0.5, 'cm'))
  
  return(plot)
  
}

# repressed
h0_to_h6_rep_waffle <- plot_waffle(h0_to_h6_act_waffle_prep, "h0 to h6")
h3_to_h9_rep_waffle <- plot_waffle(h3_to_h9_act_waffle_prep, "h3 to h9")
h6_to_h12_rep_waffle <- plot_waffle(h6_to_h12_act_waffle_prep, "h6 to h12")
h9_to_h15_rep_waffle <- plot_waffle(h9_to_h15_act_waffle_prep, "h9 to h15")
h12_to_h18_rep_waffle <- plot_waffle(h12_to_h18_act_waffle_prep, "h12 to h18")
h15_to_h21_rep_waffle <- plot_waffle(h15_to_h21_act_waffle_prep, "h15 to h21")
h18_to_h24_rep_waffle <- plot_waffle(h18_to_h24_act_waffle_prep, "h18 to h24")

# activated
h0_to_h6_act_waffle <- plot_waffle(h0_to_h6_rep_waffle_prep, "h0 to h6")
h3_to_h9_act_waffle <- plot_waffle(h3_to_h9_rep_waffle_prep, "h3 to h9")
h6_to_h12_act_waffle <- plot_waffle(h6_to_h12_rep_waffle_prep, "h6 to h12")
h9_to_h15_act_waffle <- plot_waffle(h9_to_h15_rep_waffle_prep, "h9 to h15")
h12_to_h18_act_waffle <- plot_waffle(h12_to_h18_rep_waffle_prep, "h12 to h18")
h15_to_h21_act_waffle <- plot_waffle(h15_to_h21_rep_waffle_prep, "h15 to h21")
h18_to_h24_act_waffle <- plot_waffle(h18_to_h24_rep_waffle_prep, "h18 to h24")

ggplot(h0_to_h6_tags_count_plotting_data_waffle_rep, aes(fill = group, values = sum_isoforms)) +
  geom_waffle(n_rows = 4, size = 0.2, colour = "grey30", flip = TRUE, make_proportional = TRUE) +
  scale_fill_manual(name = NULL,
                    values = c('#bf5b17','#beaed4','#386cb0','#7fc97f')) +
  coord_equal() +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE, title = "Isoform")) +
  ggtitle('h0 to h6', subtitle = '(n = 93)') +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "left",
        legend.key.size = unit(0.5, 'cm'))

ggsave('waffle_h0_to_h6_rep.png', height = 4, width = 1.8, units = 'in')

write_csv(Isoforms_consistent_h0_to_h6_ranked_rep, 'Isoforms_consistent_h0_to_h6_ranked_rep.csv')

# 1.9 Odds ratio tile maps ----

# *1.9.1 repressed h0 to h15----
# combine consistent groups
repressed_h0_to_h15 <- bind_rows(h0_to_h6_tags_rep, h3_to_h9_tags_rep, h6_to_h12_tags_rep, h9_to_h15_tags_rep) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(repressed_h0_to_h15, 'repressed_h0_to_h15_isoforms.csv')

repressed_h0_to_h12_new <- bind_rows(h0_to_h6_tags_rep, h6_to_h12_tags_rep) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(repressed_h0_to_h12_new, 'repressed_h0_to_h12_new_isoforms.csv')

# convert AGI_IsoformTag into AGI - could be multiple AGIs i.e. ATxxxxxxx.1 and ATxxxxxxx_ID2 will aggregate to two times ATxxxxxxx

repressed_h0_to_h15_AGI <- repressed_h0_to_h15 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9))

repressed_h0_to_h12_new_AGI <- repressed_h0_to_h12_new %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) 

# select just the AGI and reduce to distinct AGI

repressed_h0_to_h15_AGI_distinct <- repressed_h0_to_h15_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(repressed_h0_to_h15_AGI_distinct, 'repressed_h0_to_h15_AGI_distinct.csv')

repressed_h0_to_h12_new_AGI_distinct <- repressed_h0_to_h12_new_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(repressed_h0_to_h12_new_AGI_distinct, 'repressed_h0_to_h12_new_AGI_distinct.csv')

# associate with a TFcluster number
inner_join_clusters_h0_h15_distinct_rep <- inner_join(TF_network_clusters, repressed_h0_to_h15_AGI_distinct)

inner_join_clusters_h0_h15_distinct_rep_count <- inner_join_clusters_h0_h15_distinct_rep %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

rep_h0_h15_fishers_prep <- inner_join_clusters_h0_h15_distinct_rep_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

rep_h0_h15_fishers <- rep_h0_h15_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(rep_h0_h15_fishers)[5] <- "Odds_Ratio"

rep_h0_h15_fishers_p <- rep_h0_h15_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(rep_h0_h15_fishers_p)[5] <- "p_value"

rep_h0_h15_fishers_Odds <- rep_h0_h15_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

rep_h0_h15_fishers_pval <- rep_h0_h15_fishers_p %>% 
  dplyr::select(p_value)

rep_inner_join_clusters_h0_h15_distinct_count_Odds_pval <- inner_join_clusters_h0_h15_distinct_rep_count %>% 
  bind_cols(rep_h0_h15_fishers_Odds, rep_h0_h15_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(rep_inner_join_clusters_h0_h15_distinct_count_Odds_pval, 'h0_h15_repressed_odds_p_val.csv')

rep_h0_h15_odds_heatmap1 <- rep_inner_join_clusters_h0_h15_distinct_count_Odds_pval %>%
  mutate(time_points = "early") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

rep_h0_h15_odds_heatmap2 <- rep_h0_h15_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  #ggplot(aes(x = cluster, y= time_points, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="Cluster") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 14),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 14),
        legend.key.height = grid::unit(0.8,"cm"),
        legend.key.width = grid::unit(0.8,"cm"),
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14, vjust = 0.2,colour = 'grey30'),
        axis.title.y = element_text(size = 14, colour = 'grey30'),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank())
#       #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
#       plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 1)) + 
# ggtitle('Repressed in Col-0') 

rep_h0_h15_odds_heatmap2

#ggsave(rep_h0_h15_odds_heatmap2, filename="rep_h0_h15_odds_heatmap2.png",height=9,width=3,units="in",dpi=200)

# rep h0 to h12----

# associate with a TFcluster number
inner_join_clusters_h0_h12_distinct_rep <- inner_join(TF_network_clusters, repressed_h0_to_h12_new_AGI_distinct)

inner_join_clusters_h0_h12_distinct_rep_count <- inner_join_clusters_h0_h12_distinct_rep %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

rep_h0_h12_fishers_prep <- inner_join_clusters_h0_h12_distinct_rep_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

rep_h0_h12_fishers <- rep_h0_h12_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(rep_h0_h12_fishers)[5] <- "Odds_Ratio"

rep_h0_h12_fishers_p <- rep_h0_h12_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(rep_h0_h12_fishers_p)[5] <- "p_value"

rep_h0_h12_fishers_Odds <- rep_h0_h12_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

rep_h0_h12_fishers_pval <- rep_h0_h12_fishers_p %>% 
  dplyr::select(p_value)

rep_inner_join_clusters_h0_h12_distinct_count_Odds_pval <- inner_join_clusters_h0_h12_distinct_rep_count %>% 
  bind_cols(rep_h0_h12_fishers_Odds, rep_h0_h12_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(rep_inner_join_clusters_h0_h12_distinct_count_Odds_pval, 'h0_h12_repressed_odds_p_val.csv')

rep_h0_h12_odds_heatmap1 <- rep_inner_join_clusters_h0_h12_distinct_count_Odds_pval %>%
  mutate(time_points = "h0 to h12") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 2.5, 5, 7.5, max(Odds_Ratio)), labels=c("0-1","1-2.5","2.5-5","5-7.5",">7.5")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

rep_h0_h12_odds_heatmap2 <- rep_h0_h12_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  #ggplot(aes(x = cluster, y= time_points, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="Cluster") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c')), drop = FALSE) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 10),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 10),
        legend.key.height = grid::unit(0.4,"cm"),
        legend.key.width = grid::unit(0.4,"cm"),
        axis.text.x = element_text(size = 10,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, vjust = 0.2,colour = 'grey30'),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 0.5)) + 
  ggtitle('Repressed in Col-0') 

rep_h0_h12_odds_heatmap2

ggsave(rep_h0_h12_odds_heatmap2, filename="rep_h0_h12_odds_heatmap2.png",height=9,width=3,units="in",dpi=200)


# *1.9.2 activated h0 to h15----

activated_h0_to_h15 <- bind_rows(h0_to_h6_tags_act, h3_to_h9_tags_act, h6_to_h12_tags_act, h9_to_h15_tags_act) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(activated_h0_to_h15, 'activated_h0_to_h15_isoforms.csv')

activated_h0_to_h12 <- bind_rows(h0_to_h6_tags_act, h6_to_h12_tags_act) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(activated_h0_to_h12, 'activated_h0_to_h12_isoforms.csv')

# convert AGI_IsoformTag into AGI - could be multiple AGIs i.e. ATxxxxxxx.1 and ATxxxxxxx_ID2 will aggregate to two times ATxxxxxxx
activated_h0_to_h15_AGI <- activated_h0_to_h15 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9))

activated_h0_to_h12_AGI <- activated_h0_to_h12 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) 

# select just the AGI and reduce to distinct AGI

activated_h0_to_h15_AGI_distinct <- activated_h0_to_h15_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(activated_h0_to_h15_AGI_distinct, 'activated_h0_to_h15_AGI_distinct.csv')

activated_h0_to_h12_AGI_distinct <- activated_h0_to_h12_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(activated_h0_to_h12_AGI_distinct, 'activated_h0_to_h12_AGI_distinct.csv')

# associate with a TFcluster number
inner_join_clusters_h0_h15_distinct_act <- inner_join(TF_network_clusters, activated_h0_to_h15_AGI_distinct)

inner_join_clusters_h0_h15_distinct_act_count <- inner_join_clusters_h0_h15_distinct_act %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

act_h0_h15_fishers_prep <- inner_join_clusters_h0_h15_distinct_act_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

act_h0_h15_fishers <- act_h0_h15_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(act_h0_h15_fishers)[5] <- "Odds_Ratio"

act_h0_h15_fishers_p <- act_h0_h15_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(act_h0_h15_fishers_p)[5] <- "p_value"

act_h0_h15_fishers_Odds <- act_h0_h15_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

act_h0_h15_fishers_pval <- act_h0_h15_fishers_p %>% 
  dplyr::select(p_value)

act_inner_join_clusters_h0_h15_distinct_count_Odds_pval <- inner_join_clusters_h0_h15_distinct_act_count %>% 
  bind_cols(act_h0_h15_fishers_Odds, act_h0_h15_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(act_inner_join_clusters_h0_h15_distinct_count_Odds_pval, 'h0_h15_activated_odds_p_val.csv')

act_h0_h15_odds_heatmap1 <- act_inner_join_clusters_h0_h15_distinct_count_Odds_pval %>%
  mutate(time_points = "early") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

act_h0_h15_odds_heatmap2 <- act_h0_h15_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  #labs(x="",y="Cluster") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank())
#plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
#       plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 0.5),
#       axis.title.x = element_blank()) + 
# ggtitle('Activated in Col-0')
# scale_fill_manual(values=rev(brewer.pal(4, "YlGnBu")))
act_h0_h15_odds_heatmap2

#ggsave(act_h0_h15_odds_heatmap2, filename="act_h0_h15_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# act h0 to h12----

# associate with a TFcluster number
inner_join_clusters_h0_h12_distinct_act <- inner_join(TF_network_clusters, activated_h0_to_h12_AGI_distinct)

inner_join_clusters_h0_h12_distinct_act_count <- inner_join_clusters_h0_h12_distinct_act %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

act_h0_h12_fishers_prep <- inner_join_clusters_h0_h12_distinct_act_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

act_h0_h12_fishers <- act_h0_h12_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(act_h0_h12_fishers)[5] <- "Odds_Ratio"

act_h0_h12_fishers_p <- act_h0_h12_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(act_h0_h12_fishers_p)[5] <- "p_value"

act_h0_h12_fishers_Odds <- act_h0_h12_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

act_h0_h12_fishers_pval <- act_h0_h12_fishers_p %>% 
  dplyr::select(p_value)

act_inner_join_clusters_h0_h12_distinct_count_Odds_pval <- inner_join_clusters_h0_h12_distinct_act_count %>% 
  bind_cols(act_h0_h12_fishers_Odds, act_h0_h12_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(act_inner_join_clusters_h0_h12_distinct_count_Odds_pval, 'h0_h12_activated_odds_p_val.csv')

act_h0_h12_odds_heatmap1 <- act_inner_join_clusters_h0_h12_distinct_count_Odds_pval %>%
  mutate(time_points = "h0 to h12") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 2.5, 5, 7.5, max(Odds_Ratio)), labels=c("0-1","1-2.5","2.5-5","5-7.5",">7.5")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

act_h0_h12_odds_heatmap2 <- act_h0_h12_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  #labs(x="",y="Cluster") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c')), drop = FALSE) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 10),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 10),
        legend.key.height = grid::unit(0.4,"cm"),
        legend.key.width = grid::unit(0.4,"cm"),
        axis.text.x = element_text(size = 10,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, vjust = 0.2,colour = 'grey30'),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 0.5),
        axis.title.x = element_blank()) + 
  ggtitle('Activated in Col-0')
# scale_fill_manual(values=rev(brewer.pal(4, "YlGnBu")))
act_h0_h12_odds_heatmap2

ggsave(act_h0_h12_odds_heatmap2, filename="act_h0_h12_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# *1.9.3 repressed h9 to h24----
repressed_h9_to_h24 <- bind_rows(h9_to_h15_tags_rep, h12_to_h18_tags_rep, h15_to_h21_tags_rep, h18_to_h24_tags_rep) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(repressed_h9_to_h24, 'repressed_h9_to_h24_isoforms.csv')

repressed_h12_to_h24 <- bind_rows(h12_to_h18_tags_rep, h18_to_h24_tags_rep) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(repressed_h12_to_h24, 'repressed_h12_to_h24_isoforms.csv')

# convert AGI_IsoformTag into AGI - could be multiple AGIs i.e. ATxxxxxxx.1 and ATxxxxxxx_ID2 will aggregate to two times ATxxxxxxx
repressed_h9_to_h24_AGI <- repressed_h9_to_h24 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) 

repressed_h12_to_h24_AGI <- repressed_h12_to_h24 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) 

# select just the AGI and reduce to distinct AGI

repressed_h9_to_h24_AGI_distinct <- repressed_h9_to_h24_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(repressed_h9_to_h24_AGI_distinct, 'repressed_h9_to_h24_AGI_distinct.csv')

repressed_h12_to_h24_AGI_distinct <- repressed_h12_to_h24_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(repressed_h12_to_h24_AGI_distinct, 'repressed_h12_to_h24_AGI_distinct.csv')

# associate with a TFcluster number
inner_join_clusters_h9_h24_distinct_rep <- inner_join(TF_network_clusters, repressed_h9_to_h24_AGI_distinct)

inner_join_clusters_h9_h24_distinct_rep_count <- inner_join_clusters_h9_h24_distinct_rep %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

rep_h9_h24_fishers_prep <- inner_join_clusters_h9_h24_distinct_rep_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

rep_h9_h24_fishers <- rep_h9_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(rep_h9_h24_fishers)[5] <- "Odds_Ratio"

rep_h9_h24_fishers_p <- rep_h9_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(rep_h9_h24_fishers_p)[5] <- "p_value"

rep_h9_h24_fishers_Odds <- rep_h9_h24_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

rep_h9_h24_fishers_pval <- rep_h9_h24_fishers_p %>% 
  dplyr::select(p_value)

rep_inner_join_clusters_h9_h24_distinct_count_Odds_pval <- inner_join_clusters_h9_h24_distinct_rep_count %>% 
  bind_cols(rep_h9_h24_fishers_Odds, rep_h9_h24_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(rep_inner_join_clusters_h9_h24_distinct_count_Odds_pval, 'h9_h24_repressed_odds_p_val.csv')

filter(rep_inner_join_clusters_h9_h24_distinct_count_Odds_pval, Odds_Ratio >= 5)

rep_h9_h24_odds_heatmap1 <- rep_inner_join_clusters_h9_h24_distinct_count_Odds_pval %>%
  mutate(time_points = "late") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

rep_h9_h24_odds_heatmap2 <- rep_h9_h24_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 14),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 14),
        legend.key.height = grid::unit(0.8,"cm"),
        legend.key.width = grid::unit(0.8,"cm"),
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        #plot.title = element_text(colour = 'grey30', size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

rep_h9_h24_odds_heatmap2

#ggsave(rep_h9_h24_odds_heatmap2, filename="rep_h9_h24_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# rep h12 to h24----

# associate with a TFcluster number
inner_join_clusters_h12_h24_distinct_rep <- inner_join(TF_network_clusters, repressed_h12_to_h24_AGI_distinct)

inner_join_clusters_h12_h24_distinct_rep_count <- inner_join_clusters_h12_h24_distinct_rep %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

rep_h12_h24_fishers_prep <- inner_join_clusters_h12_h24_distinct_rep_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

rep_h12_h24_fishers <- rep_h12_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(rep_h12_h24_fishers)[5] <- "Odds_Ratio"

rep_h12_h24_fishers_p <- rep_h12_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(rep_h12_h24_fishers_p)[5] <- "p_value"

rep_h12_h24_fishers_Odds <- rep_h12_h24_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

rep_h12_h24_fishers_pval <- rep_h12_h24_fishers_p %>% 
  dplyr::select(p_value)

rep_inner_join_clusters_h12_h24_distinct_count_Odds_pval <- inner_join_clusters_h12_h24_distinct_rep_count %>% 
  bind_cols(rep_h12_h24_fishers_Odds, rep_h12_h24_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(rep_inner_join_clusters_h12_h24_distinct_count_Odds_pval, 'h12_h24_repressed_odds_p_val.csv')

filter(rep_inner_join_clusters_h12_h24_distinct_count_Odds_pval, Odds_Ratio >= 5)

rep_h12_h24_odds_heatmap1 <- rep_inner_join_clusters_h12_h24_distinct_count_Odds_pval %>%
  mutate(time_points = "h12 to h24") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 2.5, 5, 7.5, max(Odds_Ratio)), labels=c("0-1","1-2.5","2.5-5","5-7.5",">7.5")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

rep_h12_h24_odds_heatmap2 <- rep_h12_h24_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c')), drop = FALSE) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 10),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 10),
        legend.key.height = grid::unit(0.4,"cm"),
        legend.key.width = grid::unit(0.4,"cm"),
        axis.text.x = element_text(size = 10,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        #plot.title = element_text(colour = 'grey30', size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

rep_h12_h24_odds_heatmap2

ggsave(rep_h12_h24_odds_heatmap2, filename="rep_h12_h24_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# *1.9.4 activated h9 to h24----
activated_h9_to_h24 <- bind_rows(h9_to_h15_tags_act, h12_to_h18_tags_act, h15_to_h21_tags_act, h18_to_h24_tags_act) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(activated_h9_to_h24, 'activated_h9_to_h24_isoforms.csv')

activated_h12_to_h24 <- bind_rows(h12_to_h18_tags_act, h18_to_h24_tags_act) %>%
  dplyr::select(Isoform) %>% 
  distinct()

write_csv(activated_h12_to_h24, 'activated_h12_to_h24_isoforms.csv')

# convert AGI_IsoformTag into AGI - could be multiple AGIs i.e. ATxxxxxxx.1 and ATxxxxxxx_ID2 will aggregate to two times ATxxxxxxx
activated_h9_to_h24_AGI <- activated_h9_to_h24 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9))

activated_h12_to_h24_AGI <- activated_h12_to_h24 %>% 
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) 

# select just the AGI and reduce to distinct AGI

activated_h9_to_h24_AGI_distinct <- activated_h9_to_h24_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(activated_h9_to_h24_AGI_distinct, 'activated_h9_to_h24_AGI_distinct.csv')

activated_h12_to_h24_AGI_distinct <- activated_h12_to_h24_AGI %>%
  dplyr::select(gene_ID) %>% 
  distinct()

write_csv(activated_h12_to_h24_AGI_distinct, 'activated_h12_to_h24_AGI_distinct.csv')

# associate with a TFcluster number
inner_join_clusters_h9_h24_distinct_act <- inner_join(TF_network_clusters, activated_h9_to_h24_AGI_distinct)

inner_join_clusters_h9_h24_distinct_act_count <- inner_join_clusters_h9_h24_distinct_act %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

act_h9_h24_fishers_prep <- inner_join_clusters_h9_h24_distinct_act_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

act_h9_h24_fishers <- act_h9_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(act_h9_h24_fishers)[5] <- "Odds_Ratio"

act_h9_h24_fishers_p <- act_h9_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(act_h9_h24_fishers_p)[5] <- "p_value"

act_h9_h24_fishers_Odds <- act_h9_h24_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

act_h9_h24_fishers_pval <- act_h9_h24_fishers_p %>% 
  dplyr::select(p_value)

act_inner_join_clusters_h9_h24_distinct_count_Odds_pval <- inner_join_clusters_h9_h24_distinct_act_count %>% 
  bind_cols(act_h9_h24_fishers_Odds, act_h9_h24_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(act_inner_join_clusters_h9_h24_distinct_count_Odds_pval, 'h9_h24_activated_odds_p_val.csv')

act_h9_h24_odds_heatmap1 <- act_inner_join_clusters_h9_h24_distinct_count_Odds_pval %>%
  mutate(time_points = "late") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

act_h9_h24_odds_heatmap2 <- act_h9_h24_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30", size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  #labs(x="",y="") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 14),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 14),
        legend.key.height = grid::unit(0.8,"cm"),
        legend.key.width = grid::unit(0.8,"cm"),
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        #plot.title = element_text(colour = 'grey30', size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

act_h9_h24_odds_heatmap2

#ggsave(act_h9_h24_odds_heatmap2, filename="act_h9_h24_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# act h12 to h24----

# associate with a TFcluster number
inner_join_clusters_h12_h24_distinct_act <- inner_join(TF_network_clusters, activated_h12_to_h24_AGI_distinct)

inner_join_clusters_h12_h24_distinct_act_count <- inner_join_clusters_h12_h24_distinct_act %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

act_h12_h24_fishers_prep <- inner_join_clusters_h12_h24_distinct_act_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

act_h12_h24_fishers <- act_h12_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(act_h12_h24_fishers)[5] <- "Odds_Ratio"

act_h12_h24_fishers_p <- act_h12_h24_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(act_h12_h24_fishers_p)[5] <- "p_value"

act_h12_h24_fishers_Odds <- act_h12_h24_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

act_h12_h24_fishers_pval <- act_h12_h24_fishers_p %>% 
  dplyr::select(p_value)

act_inner_join_clusters_h12_h24_distinct_count_Odds_pval <- inner_join_clusters_h12_h24_distinct_act_count %>% 
  bind_cols(act_h12_h24_fishers_Odds, act_h12_h24_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(act_inner_join_clusters_h12_h24_distinct_count_Odds_pval, 'h12_h24_activated_odds_p_val.csv')

act_h12_h24_odds_heatmap1 <- act_inner_join_clusters_h12_h24_distinct_count_Odds_pval %>%
  mutate(time_points = "late") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 2.5, 5, 7.5, max(Odds_Ratio)), labels=c("0-1","1-2.5","2.5-5","5-7.5",">7.5")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

act_h12_h24_odds_heatmap2 <- act_h12_h24_odds_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30", size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  #labs(x="",y="") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c')), drop = FALSE) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 10),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 10),
        legend.key.height = grid::unit(0.4,"cm"),
        legend.key.width = grid::unit(0.4,"cm"),
        axis.text.x = element_text(size = 10,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        #plot.title = element_text(colour = 'grey30', size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

act_h12_h24_odds_heatmap2

ggsave(act_h12_h24_odds_heatmap2, filename="act_h12_h24_odds_heatmap2_new_scale.png",height=9,width=3,units="in",dpi=200)

# 1.10 Patchworked Odds Ratio plots----
# TF cluster sizes plot
TF_cluster_sizes_plot <- act_inner_join_clusters_h9_h24_distinct_count_Odds_pval %>%
  ggplot(aes(x = cluster, y= cluster_number)) +
  geom_bar(stat="identity", colour = 'grey30', fill="grey90") +
  scale_x_continuous(expand=c(0, 0), breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)) +
  scale_y_reverse(expand = c(0,0), breaks = c(0, 100, 200, 300, 400, 500, 600, 700)) +
  coord_flip() +
  theme_minimal() +
  #ylab('Cluster size') +
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.x=element_line(size=0.4),
        plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 0.5)) 
#+ggtitle('Cluster sizes')

TF_cluster_sizes_plot

rep_patched_tiles <- TF_cluster_sizes_plot + (rep_h0_h15_odds_heatmap2 + theme(legend.position = "none")) + rep_h9_h24_odds_heatmap2 +
  plot_layout(guides = 'auto')

# rep patched tiles overlap----
rep_patched_tiles2 <- TF_cluster_sizes_plot + rep_h0_h15_odds_heatmap2 + theme(legend.position = "none") + rep_h9_h24_odds_heatmap2 + theme(legend.position = "none")

rep_patched_tiles2

ggsave(rep_patched_tiles2, filename="rep_patched_tiles_alt_scale2.png", height = 9, width = 4.5, units = "in", dpi = 300)

# rep patched tiles no overlap----

rep_patched_tiles3 <- TF_cluster_sizes_plot + rep_h0_h12_odds_heatmap2 + theme(legend.position = "none") + rep_h12_h24_odds_heatmap2 + theme(legend.position = "none")

rep_patched_tiles3

ggsave(rep_patched_tiles3, filename="rep_patched_tiles_alt_scale3.png", height = 9, width = 4.5, units = "in", dpi = 300)

# act patched tiles overlap----

act_patched_tiles <- TF_cluster_sizes_plot + (act_h0_h15_odds_heatmap2 + theme(legend.position = "none")) + act_h9_h24_odds_heatmap2 +
  plot_layout(guides = 'auto')

act_patched_tiles2 <- act_h0_h15_odds_heatmap2 + theme(legend.position = "none") + act_h9_h24_odds_heatmap2 +
  plot_layout(guides = 'auto')

act_patched_tiles2

ggsave(act_patched_tiles2, filename="act_patched_tiles_alt_scale2.png", height = 9, width = 4.5, units = "in", dpi = 300)

# act patched tiles no overlap----

act_patched_tiles3 <- act_h0_h12_odds_heatmap2 + act_h12_h24_odds_heatmap2 + plot_layout(guides = "collect")

act_patched_tiles3

ggsave(act_patched_tiles3, filename="act_patched_tiles_alt_scale3.png", height = 9, width = 4.5, units = "in", dpi = 300)

patch_plot <- TF_cluster_sizes_plot + (rep_h0_h12_odds_heatmap2 + rep_h12_h24_odds_heatmap2 + act_h0_h12_odds_heatmap2 + act_h12_h24_odds_heatmap2 + plot_layout(ncol = 4) + plot_layout(guides = "collect")) 

ggsave(patch_plot2, filename="patch_plot2.png", height = 9, width = 7, units = "in", dpi = 300)

patch_plot2 <- TF_cluster_sizes_plot + (rep_h0_h15_odds_heatmap2 + rep_h9_h24_odds_heatmap2 + act_h0_h15_odds_heatmap2 + act_h9_h24_odds_heatmap2 + plot_layout(ncol = 4) + plot_layout(guides = "collect"))

patch_plot3

# 1.11 z-score plots----
# *1.11.1 cluster 9----
# h0_15
# get cluster 9 genes
# this is just gene_IDs
cluster9_h0_h15_rep <- inner_join_clusters_h0_h15_distinct_rep %>% 
  filter(cluster %in% '9')

# get the associated isoforms
cluster9_genes_isoforms_rep <- repressed_h0_to_h15_AGI %>% 
  inner_join(cluster9_h0_h15_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster9_genes_isoforms_rep, 'cluster9_genes_isoforms_rep.csv')

# write the just gene_IDs to file
write_csv(cluster9_h0_h15_rep, 'cluster9_h0_h15_rep.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster9_gene_ids_rep <- read_csv('cluster9_h0_h15_rep_with_gene_IDs.csv')

cluster9_h0_h15_rep_with_gene_IDs <- cluster9_gene_ids_rep %>%
  inner_join(cluster9_h0_h15_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster9_h0_h15_rep_with_gene_IDs, 'cluster9_h0_h15_rep_with_gene_IDs.csv')

#cluster9_h0_h15_rep_isoforms
cluster9_h0_h15_rep_isoforms_z_scores <- cluster9_genes_isoforms_rep %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 9', subtitle = '132 isoforms (101 gene loci)')

ggsave('cluster9_h0_h15_rep_isoforms_z_scores2.png', height = 5, width = 4, units = 'in')

cluster9_h0_h15_rep_isoforms_z_scores

# h0_12----
# get cluster 9 genes
# this is just gene_IDs
cluster9_h0_h12_rep <- inner_join_clusters_h0_h12_distinct_rep %>% 
  filter(cluster %in% '9')

# get the associated isoforms
cluster9_genes_isoforms_rep_new <- repressed_h0_to_h12_new_AGI %>% 
  inner_join(cluster9_h0_h12_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster9_genes_isoforms_rep_new, 'cluster9_genes_isoforms_rep_new.csv')

# write the just gene_IDs to file
write_csv(cluster9_h0_h12_rep, 'cluster9_h0_h12_rep.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster9_gene_ids_rep_new <- read_csv('cluster9_h0_h12_rep_with_gene_IDs.csv')

cluster9_h0_h12_rep_with_gene_IDs <- cluster9_gene_ids_rep_new %>%
  inner_join(cluster9_h0_h12_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster9_h0_h12_rep_with_gene_IDs, 'cluster9_h0_h12_rep_with_gene_IDs.csv')

#cluster9_h0_h12_rep_isoforms
cluster9_h0_h12_rep_isoforms_z_scores <- cluster9_genes_isoforms_rep_new %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 9', subtitle = '111 isoforms (85 gene loci)')

ggsave('cluster9_h0_h12_rep_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster9_h0_h12_rep_isoforms_z_scores

# *1.11.2 cluster 20----

#get cluster 20 genes
# this is just gene_IDs
cluster20_h0_h15_rep <- inner_join_clusters_h0_h15_distinct_rep %>% 
  filter(cluster %in% '20')

# get associated isoforms
cluster20_genes_isoforms_rep <- repressed_h0_to_h15_AGI %>% 
  inner_join(cluster20_h0_h15_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster20_genes_isoforms_rep, 'cluster20_genes_isoforms_rep.csv')

# write the just gene_IDs to file
write_csv(cluster20_h0_h15_rep, 'cluster20_h0_h15_rep.csv')

cluster20_gene_ids_rep <- read_csv('cluster20_h0_h15_rep_gene_id.csv')

cluster20_h0_h15_rep_with_gene_IDs <- cluster20_gene_ids_rep %>%
  inner_join(cluster20_h0_h15_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster20_h0_h15_rep_with_gene_IDs, 'cluster20_h0_h15_rep_with_gene_IDs.csv')

cluster20_h0_h15_rep_isoforms_z_scores <- cluster20_genes_isoforms_rep %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 20', subtitle = '70 isoforms (53 gene loci)')

ggsave('cluster20_h0_h15_rep_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster20_h0_h15_rep_isoforms_z_scores

# h0_12----

#get cluster 20 genes
# this is just gene_IDs
cluster20_h0_h12_rep <- inner_join_clusters_h0_h12_distinct_rep %>% 
  filter(cluster %in% '20')

# get associated isoforms
cluster20_genes_isoforms_rep_new <- repressed_h0_to_h12_new_AGI %>% 
  inner_join(cluster20_h0_h12_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster20_genes_isoforms_rep_new, 'cluster20_genes_isoforms_rep_new.csv')

# write the just gene_IDs to file
write_csv(cluster20_h0_h12_rep, 'cluster20_h0_h12_rep.csv')

cluster20_gene_ids_rep_new <- read_csv('cluster20_h0_h12_rep_gene_id.csv')

cluster20_h0_h12_rep_with_gene_IDs <- cluster20_gene_ids_rep_new %>%
  inner_join(cluster20_h0_h12_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster20_h0_h12_rep_with_gene_IDs, 'cluster20_h0_h12_rep_with_gene_IDs.csv')

cluster20_h0_h12_rep_isoforms_z_scores <- cluster20_genes_isoforms_rep_new %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 20', subtitle = '56 isoforms (45 gene loci)')

ggsave('cluster20_h0_h12_rep_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster20_h0_h12_rep_isoforms_z_scores

# *1.11.3 cluster 25----

#get cluster 25 genes
# this is just gene_IDs
cluster25_h0_h15_rep <- inner_join_clusters_h0_h15_distinct_rep %>% 
  filter(cluster %in% '25')

# get associated isoforms
cluster25_genes_isoforms_rep <- repressed_h0_to_h15_AGI %>% 
  inner_join(cluster25_h0_h15_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster25_genes_isoforms_rep, 'cluster25_genes_isoforms_rep.csv')

# write the just gene_IDs to file
write_csv(cluster25_h0_h15_rep, 'cluster25_h0_h15_rep.csv')

cluster25_gene_ids_rep <- read_csv('cluster25_h0_h15_rep_gene_id.csv')

cluster25_h0_h15_rep_with_gene_IDs <- cluster25_gene_ids_rep %>%
  inner_join(cluster25_h0_h15_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster25_h0_h15_rep_with_gene_IDs, 'cluster25_h0_h15_rep_with_gene_IDs.csv')

cluster25_h0_h15_rep_isoforms_z_scores <- cluster25_h0_h15_rep_isoforms %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(1,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 25', subtitle = '56 isoforms (38 gene loci)')

ggsave('cluster25_h0_h15_rep_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster25_h0_h15_rep_isoforms_z_scores

# h0_12----

#get cluster 25 genes
# this is just gene_IDs
cluster25_h0_h12_rep <- inner_join_clusters_h0_h12_distinct_rep %>% 
  filter(cluster %in% '25')

# get associated isoforms
cluster25_genes_isoforms_rep_new <- repressed_h0_to_h12_new_AGI %>% 
  inner_join(cluster25_h0_h12_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster25_genes_isoforms_rep_new, 'cluster25_genes_isoforms_rep_new.csv')

# write the just gene_IDs to file
write_csv(cluster25_h0_h12_rep, 'cluster25_h0_h12_rep.csv')

cluster25_gene_ids_rep_new <- read_csv('cluster25_h0_h12_rep_gene_id.csv')

cluster25_h0_h12_rep_with_gene_IDs <- cluster25_gene_ids_rep_new %>%
  inner_join(cluster25_h0_h12_rep, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster25_h0_h12_rep_with_gene_IDs, 'cluster25_h0_h12_rep_with_gene_IDs.csv')

cluster25_h0_h12_rep_isoforms_z_scores <- cluster25_genes_isoforms_rep_new %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 25', subtitle = '47 isoforms (32 gene loci)')

ggsave('cluster25_h0_h12_rep_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster25_h0_h12_rep_isoforms_z_scores

# *1.11.4 cluster 35----

# get cluster 35 genes
# this is just gene_IDs
cluster35_h0_h15_rep <- inner_join_clusters_h0_h15_distinct_rep %>% 
  filter(cluster %in% '35')

# get associated isoforms
cluster35_genes_isoforms_rep <- repressed_h0_to_h15_AGI %>% 
  inner_join(cluster35_h0_h15_rep, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

#write_csv(cluster35_genes_isoforms_rep, 'cluster35_genes_isoforms_rep.csv')

# write the just gene_IDs to file
write_csv(cluster35_h0_h15_rep, 'cluster35_h0_h15_rep.csv')

cluster35_gene_ids_rep <- read_csv('cluster35_h0_h15_rep_gene_id.csv')

cluster35_h0_h15_rep_with_gene_IDs <- cluster35_gene_ids_rep %>%
  inner_join(cluster35_h0_h15_rep, by= 'gene_ID') %>%
  select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster35_h0_h15_rep_with_gene_IDs, 'cluster35_h0_h15_rep_with_gene_IDs.csv')

cluster35_h0_h15_rep_isoforms_z_scores <- cluster35_genes_isoforms_rep %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 35', subtitle = '22 isoforms (19 gene loci)')

ggsave('cluster35_h0_h15_rep_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster35_h0_h15_rep_isoforms_z_scores

# *1.11.5 cluster 17----

# get cluster 17 genes
# this is just gene_IDs
cluster17_h9_h24_act <- inner_join_clusters_h9_h24_distinct_act %>% 
  filter(cluster %in% '17')

# get associated isoforms
cluster17_genes_isoforms_act <- activated_h9_to_h24_AGI %>% 
  inner_join(cluster17_h9_h24_act, by = 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster17_genes_isoforms_act, 'cluster17_genes_isoforms_act.csv')

# write the just gene_IDs to file
write_csv(cluster17_h9_h24_act, 'cluster17_h9_h24_act.csv')

cluster17_gene_ids_act <- read_csv('cluster17_h9_h24_act_gene_id.csv')

cluster17_h9_h24_act_with_gene_IDs <- cluster17_gene_ids_act %>%
  inner_join(cluster17_h9_h24_act, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster17_h9_h24_act_with_gene_IDs, 'cluster17_h9_h24_act_with_gene_IDs.csv')

cluster17_h9_h24_act_isoforms_z_scores <- cluster17_genes_isoforms_act %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#B3E3A0"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 17', subtitle = '45 isoforms (45 gene loci)')

ggsave('cluster17_h9_h24_act_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster17_h9_h24_act_isoforms_z_scores

# h12_24----

# get cluster 17 genes
# this is just gene_IDs
cluster17_h12_h24_act <- inner_join_clusters_h12_h24_distinct_act %>% 
  filter(cluster %in% '17')

# get associated isoforms
cluster17_genes_isoforms_act_new <- activated_h12_to_h24_AGI %>% 
  inner_join(cluster17_h12_h24_act, by = 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster17_genes_isoforms_act_new, 'cluster17_genes_isoforms_act_new.csv')

# write the just gene_IDs to file
write_csv(cluster17_h12_h24_act, 'cluster17_h12_h24_act.csv')

cluster17_gene_ids_act_new <- read_csv('cluster17_h12_h24_act_gene_id.csv')

cluster17_h12_h24_act_with_gene_IDs <- cluster17_gene_ids_act_new %>%
  inner_join(cluster17_h12_h24_act, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster17_h12_h24_act_with_gene_IDs, 'cluster17_h12_h24_act_with_gene_IDs.csv')

cluster17_h12_h24_act_isoforms_z_scores <- cluster17_genes_isoforms_act_new %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#B3E3A0"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 17', subtitle = '45 isoforms (45 gene loci)')

ggsave('cluster17_h12_h24_act_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster17_h12_h24_act_isoforms_z_scores

# *1.11.6 cluster 50----

# get cluster 50 genes
# this is just gene_IDs

cluster50_h0_h15_act <- inner_join_clusters_h0_h15_distinct_act %>% 
  filter(cluster %in% '50')

# get associated isoforms
cluster50_genes_isoforms_act <- activated_h0_to_h15_AGI %>% 
  inner_join(cluster50_h0_h15_act, by = 'gene_ID') %>% 
  select(-gene_ID)

#write_csv(cluster50_genes_isoforms_act, 'cluster50_genes_isoforms_act.csv')

# write the just gene_IDs to file
write_csv(cluster50_h0_h15_act, 'cluster50_h0_h15_act.csv')

cluster50_gene_ids_act <- read_csv('cluster50_h0_h15_act_gene_id.csv')

cluster50_h0_h15_act_with_gene_IDs <- cluster50_gene_ids_act %>%
  inner_join(cluster50_h0_h15_act, by= 'gene_ID') %>%
  select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster50_h0_h15_act_with_gene_IDs, 'cluster50_h0_h15_act_with_gene_IDs.csv')

cluster50_h0_h15_act_isoforms_z_scores <- cluster50_genes_isoforms_act %>% 
  inner_join(transcript_recovered_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#B3E3A0"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.85),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 50', subtitle = '12 isoforms (10 gene loci)')

ggsave('cluster50_h0_h15_act_isoforms_z_scores.png', height = 4, width = 4, units = 'in')

cluster50_h0_h15_act_isoforms_z_scores

# 1.12 TF network----

# 7302 gene loci grouped in 75 clusters (0-74)
TF_network_clusters <- read_csv('TF Network Cluster Nov2018 long format.csv') %>%
  filter(!row_number() %in% 5772)

# count of gene loci in each TF cluster
TF_network_clusters_count <-  TF_network_clusters %>% 
  group_by(cluster) %>%
  summarise(cluster_number = n())

# *1.12.1 TF network plot----
cluster_total_bar_plot <- rep_inner_join_clusters_h0_h15_distinct_count_Odds_pval %>%
  ggplot(aes(x = cluster, y= cluster_number)) +
  geom_bar(stat="identity", colour = 'grey30', fill="#7fc97f") +
  scale_x_continuous(expand=c(0, 0), breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)) +
  scale_y_reverse(expand = c(0,0), breaks = c(0, 100, 200, 300, 400, 500, 600, 700)) +
  coord_flip() +
  theme_minimal(base_size=12) +
  ylab('Cluster size') +
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.x=element_line(size=0.4))

cluster_total_bar_plot

ggsave(cluster_total_bar_plot,filename="cluster_total_bar_plot.png",height=9,width=3,units="in",dpi=200)

ggarrange(cluster_total_bar_plot, cluster_bar_plot_h3_h15, h3_h15_odds_heatmap2, 
          ncol = 3, nrow = 1) %>% 
  ggexport(filename = "test1.png")

# 1.13 ANOVA----
# *1.13.1 repressed h0 to h15----

h0_h15_extract <- transcript_RNAseq %>% 
  filter(Isoform %in% repressed_h0_to_h15_AGI$Isoform)

h0_h15_extract_pivot <- h0_h15_extract %>%
  dplyr::select(Isoform, Col.T9.rep1:Col.T15.rep1, Col.T9.rep2:Col.T15.rep2, Col.T9.rep3:Col.T15.rep3,
                Rve2.T9.rep1:Rve2.T15.rep1, Rve2.T9.rep2:Rve2.T15.rep2, Rve2.T9.rep3:Rve2.T15.rep3) %>% 
  pivot_longer(cols = c(Col.T9.rep1:Rve2.T15.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

RVE7_extract_pivot <- h0_h15_extract_pivot %>% filter(Isoform %in% 'AT1G18330_P1')

h0_h15_extract_pivot_ANOVAs <- h0_h15_extract_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

h0_h15_extract_pivot_ANOVAs %>% filter(Isoform %in% 'AT5G17300_P1')

repressed_h0_to_h15_ANOVA_prep_clusters <-repressed_h0_to_h15_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h0_h15_extract_pivot_ANOVAs_stats <- repressed_h0_to_h15_ANOVA_prep_clusters %>% 
  left_join(h0_h15_extract_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h0_h15_extract_pivot_ANOVAs_stats %>% filter(Isoform %in% 'AT1G18330_P1')

h0_h15_extract_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("rep_h0_h15_1_star.csv")

table(h0_h15_extract_pivot_ANOVAs_stats$p_flag)

h0_h15_extract_pivot_ANOVAs_stats_summary <- h0_h15_extract_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# rep h0 to h12----

h0_h12_extract <- transcript_RNAseq %>% 
  filter(Isoform %in% repressed_h0_to_h12_new_AGI$Isoform)

h0_h12_extract_pivot <- h0_h12_extract %>%
  dplyr::select(Isoform, Col.T9.rep1:Col.T14.rep1, Col.T9.rep2:Col.T14.rep2, Col.T9.rep3:Col.T14.rep3,
                Rve2.T9.rep1:Rve2.T14.rep1, Rve2.T9.rep2:Rve2.T14.rep2, Rve2.T9.rep3:Rve2.T14.rep3) %>% 
  pivot_longer(cols = c(Col.T9.rep1:Rve2.T14.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

#RVE7_extract_pivot <- h0_h12_extract_pivot %>% filter(Isoform %in% 'AT1G18330_P1')

h0_h12_extract_pivot_ANOVAs <- h0_h12_extract_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

h0_h12_extract_pivot_ANOVAs %>% filter(Isoform %in% 'AT5G17300_P1')

repressed_h0_to_h12_ANOVA_prep_clusters <-repressed_h0_to_h12_new_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h0_h12_extract_pivot_ANOVAs_stats <- repressed_h0_to_h12_ANOVA_prep_clusters %>% 
  left_join(h0_h12_extract_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h0_h12_extract_pivot_ANOVAs_stats %>% filter(Isoform %in% 'AT1G18330_P1')

h0_h12_extract_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("rep_h0_h12_1_star.csv")

#table(h0_h12_extract_pivot_ANOVAs_stats$p_flag)

h0_h12_extract_pivot_ANOVAs_stats_summary <- h0_h12_extract_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# *1.13.2 repressed h9 to h24----
h9_h24_extract <- transcript_RNAseq %>% 
  filter(Isoform %in% repressed_h9_to_h24_AGI$Isoform)

h9_h24_extract_pivot <- h9_h24_extract %>%
  dplyr::select(Isoform, Col.T13.rep1:Col.T18.rep1, Col.T13.rep2:Col.T18.rep2, Col.T13.rep3:Col.T18.rep3,
                Rve2.T13.rep1:Rve2.T18.rep1, Rve2.T13.rep2:Rve2.T18.rep2, Rve2.T13.rep3:Rve2.T18.rep3) %>% 
  pivot_longer(cols = c(Col.T13.rep1:Rve2.T18.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h9_h24_extract_pivot_ANOVAs <- h9_h24_extract_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

repressed_h9_to_h24_ANOVA_prep_clusters <-repressed_h9_to_h24_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h9_h24_extract_pivot_ANOVAs_stats <- repressed_h9_to_h24_ANOVA_prep_clusters %>% 
  left_join(h9_h24_extract_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns')) 

h9_h24_extract_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("rep_h9_h24_1_star.csv")

h9_h24_extract_pivot_ANOVAs_stats_summary <- h9_h24_extract_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# rep h12 to h24----

h12_h24_extract <- transcript_RNAseq %>% 
  filter(Isoform %in% repressed_h12_to_h24_AGI$Isoform)

h12_h24_extract_pivot <- h12_h24_extract %>%
  dplyr::select(Isoform, Col.T14.rep1:Col.T18.rep1, Col.T14.rep2:Col.T18.rep2, Col.T14.rep3:Col.T18.rep3,
                Rve2.T14.rep1:Rve2.T18.rep1, Rve2.T14.rep2:Rve2.T18.rep2, Rve2.T14.rep3:Rve2.T18.rep3) %>% 
  pivot_longer(cols = c(Col.T14.rep1:Rve2.T18.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h12_h24_extract_pivot_ANOVAs <- h12_h24_extract_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

repressed_h12_to_h24_ANOVA_prep_clusters <-repressed_h12_to_h24_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h12_h24_extract_pivot_ANOVAs_stats <- repressed_h12_to_h24_ANOVA_prep_clusters %>% 
  left_join(h12_h24_extract_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns')) 

h12_h24_extract_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("rep_h12_h24_1_star.csv")

h12_h24_extract_pivot_ANOVAs_stats_summary <- h12_h24_extract_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# *1.13.3 activated h0 to h15----

h0_h15_extract_act <- transcript_RNAseq %>% 
  filter(Isoform %in% activated_h0_to_h15_AGI$Isoform)

h0_h15_extract_act_pivot <- h0_h15_extract_act %>%
  dplyr::select(Isoform, Col.T9.rep1:Col.T15.rep1, Col.T9.rep2:Col.T15.rep2, Col.T9.rep3:Col.T15.rep3,
                Rve2.T9.rep1:Rve2.T15.rep1, Rve2.T9.rep2:Rve2.T15.rep2, Rve2.T9.rep3:Rve2.T15.rep3) %>% 
  pivot_longer(cols = c(Col.T9.rep1:Rve2.T15.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h0_h15_extract_act_pivot_ANOVAs <- h0_h15_extract_act_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

activated_h0_to_h15_ANOVA_prep_clusters <-activated_h0_to_h15_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h0_h15_extract_act_pivot_ANOVAs_stats <- activated_h0_to_h15_ANOVA_prep_clusters %>% 
  left_join(h0_h15_extract_act_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h0_h15_extract_act_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("act_h0_h15_1_star.csv")

h0_h15_extract_act_pivot_ANOVAs_stats_summary <- h0_h15_extract_act_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# act h0 to h12----

h0_h12_extract_act <- transcript_RNAseq %>% 
  filter(Isoform %in% activated_h0_to_h12_AGI$Isoform)

h0_h12_extract_act_pivot <- h0_h12_extract_act %>%
  dplyr::select(Isoform, Col.T9.rep1:Col.T14.rep1, Col.T9.rep2:Col.T14.rep2, Col.T9.rep3:Col.T14.rep3,
                Rve2.T9.rep1:Rve2.T14.rep1, Rve2.T9.rep2:Rve2.T14.rep2, Rve2.T9.rep3:Rve2.T14.rep3) %>% 
  pivot_longer(cols = c(Col.T9.rep1:Rve2.T14.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h0_h12_extract_act_pivot_ANOVAs <- h0_h12_extract_act_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

activated_h0_to_h12_ANOVA_prep_clusters <-activated_h0_to_h12_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h0_h12_extract_act_pivot_ANOVAs_stats <- activated_h0_to_h12_ANOVA_prep_clusters %>% 
  left_join(h0_h12_extract_act_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h0_h12_extract_act_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("act_h0_h12_1_star.csv")

h0_h12_extract_act_pivot_ANOVAs_stats_summary <- h0_h12_extract_act_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# *1.13.4 activated h9 to h24----

h9_h24_extract_act <- transcript_RNAseq %>% 
  filter(Isoform %in% activated_h9_to_h24_AGI$Isoform)

h9_h24_extract_act_pivot <- h9_h24_extract_act %>%
  dplyr::select(Isoform, Col.T13.rep1:Col.T18.rep1, Col.T13.rep2:Col.T18.rep2, Col.T13.rep3:Col.T18.rep3,
                Rve2.T13.rep1:Rve2.T18.rep1, Rve2.T13.rep2:Rve2.T18.rep2, Rve2.T13.rep3:Rve2.T18.rep3) %>% 
  pivot_longer(cols = c(Col.T13.rep1:Rve2.T18.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h9_h24_extract_act_pivot_ANOVAs <- h9_h24_extract_act_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

activated_h9_to_h24_ANOVA_prep_clusters <-activated_h9_to_h24_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h9_h24_extract_act_pivot_ANOVAs_stats <- activated_h9_to_h24_ANOVA_prep_clusters %>% 
  left_join(h9_h24_extract_act_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h9_h24_extract_act_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("act_h9_h24_1_star.csv")

h9_h24_extract_act_pivot_ANOVAs_stats_summary <- h9_h24_extract_act_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# act h12 to h24----

h12_h24_extract_act <- transcript_RNAseq %>% 
  filter(Isoform %in% activated_h12_to_h24_AGI$Isoform)

h12_h24_extract_act_pivot <- h12_h24_extract_act %>%
  dplyr::select(Isoform, Col.T14.rep1:Col.T18.rep1, Col.T14.rep2:Col.T18.rep2, Col.T14.rep3:Col.T18.rep3,
                Rve2.T14.rep1:Rve2.T18.rep1, Rve2.T14.rep2:Rve2.T18.rep2, Rve2.T14.rep3:Rve2.T18.rep3) %>% 
  pivot_longer(cols = c(Col.T14.rep1:Rve2.T18.rep3),
               names_to = c('genotype', 'time_point', 'rep'),
               names_sep = "\\.",
               values_to ='TPM')

h12_h24_extract_act_pivot_ANOVAs <- h12_h24_extract_act_pivot %>%
  group_by(Isoform) %>%
  do(broom::tidy(Anova(lm(TPM ~ time_point + genotype, data = .), type = 'II'))) %>%
  ungroup %>% 
  filter(term == 'genotype') %>%
  dplyr::select(c(1,6))

activated_h12_to_h24_ANOVA_prep_clusters <-activated_h12_to_h24_AGI %>% 
  left_join(TF_network_clusters, by = 'gene_ID') %>% 
  mutate_all(as.character) %>% 
  mutate(cluster = replace_na(cluster, 'nd'))

h12_h24_extract_act_pivot_ANOVAs_stats <- activated_h12_to_h24_ANOVA_prep_clusters %>% 
  left_join(h12_h24_extract_act_pivot_ANOVAs, by = 'Isoform') %>% 
  mutate(p_flag = case_when(p.value <= 0.0001 ~ '4_star', 
                            p.value <= 0.001 & p.value > 0.0001 ~ '3_star',
                            p.value <= 0.01 & p.value > 0.001 ~ '2_star',
                            p.value <= 0.05 & p.value > 0.01 ~ '1_star',
                            TRUE ~ 'ns'))

h12_h24_extract_act_pivot_ANOVAs_stats %>% filter(p_flag == '1_star') %>% 
  write_csv("act_h12_h24_1_star.csv")

h12_h24_extract_act_pivot_ANOVAs_stats_summary <- h12_h24_extract_act_pivot_ANOVAs_stats %>% 
  group_by(p_flag) %>% 
  summarise(total = n()) %>% 
  mutate(freq = round(total / sum(total) *100, 1)) %>% 
  arrange(desc(freq))

# 1.14 ANOVA enrichment----
# final plot made in GraphPad Prism

# *1.14.1 h0 to h15 repressed
h0_h15_extract_pivot_ANOVAs_stats_no_ns <- h0_h15_extract_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h0_h15_extract_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h0_h15_extract_pivot_ANOVAs_stats_no_ns$cluster))

h0_h15_extract_pivot_ANOVAs_stats_no_ns<-h0_h15_extract_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("rep_h0_15_ANOVA_by_R.csv")

# rep h0 to h12----

h0_h12_extract_pivot_ANOVAs_stats_no_ns <- h0_h12_extract_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h0_h12_extract_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h0_h12_extract_pivot_ANOVAs_stats_no_ns$cluster))

h0_h12_extract_pivot_ANOVAs_stats_no_ns<-h0_h12_extract_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("rep_h0_12_ANOVA_by_R.csv")

# *1.14.2 h9 to h24 repressed
h9_h24_extract_pivot_ANOVAs_stats_no_ns <- h9_h24_extract_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h9_h24_extract_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h9_h24_extract_pivot_ANOVAs_stats_no_ns$cluster))

h9_h24_extract_pivot_ANOVAs_stats_no_ns <- h9_h24_extract_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("rep_h9_24_ANOVA_by_R.csv")

# rep h12 to h24----

h12_h24_extract_pivot_ANOVAs_stats_no_ns <- h12_h24_extract_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h12_h24_extract_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h12_h24_extract_pivot_ANOVAs_stats_no_ns$cluster))

h12_h24_extract_pivot_ANOVAs_stats_no_ns <- h12_h24_extract_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("rep_h12_24_ANOVA_by_R.csv")

# *1.14.3 h0 to h15 activated
h0_h15_extract_act_pivot_ANOVAs_stats_no_ns <- h0_h15_extract_act_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h0_h15_extract_act_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h0_h15_extract_act_pivot_ANOVAs_stats_no_ns$cluster))

h0_h15_extract_act_pivot_ANOVAs_stats_no_ns <- h0_h15_extract_act_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("act_h0_15_ANOVA_by_R.csv")

# act h0 to h12----

h0_h12_extract_act_pivot_ANOVAs_stats_no_ns <- h0_h12_extract_act_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h0_h12_extract_act_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h0_h12_extract_act_pivot_ANOVAs_stats_no_ns$cluster))

h0_h12_extract_act_pivot_ANOVAs_stats_no_ns <- h0_h12_extract_act_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("act_h0_12_ANOVA_by_R.csv")

# *1.14.4 h9 to h24 activated
h9_h24_extract_act_pivot_ANOVAs_stats_no_ns <- h9_h24_extract_act_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h9_h24_extract_act_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h9_h24_extract_act_pivot_ANOVAs_stats_no_ns$cluster))

h9_h24_extract_act_pivot_ANOVAs_stats_no_ns <- h9_h24_extract_act_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("act_h9_24_ANOVA_by_R.csv")

# act h12 to h24----

h12_h24_extract_act_pivot_ANOVAs_stats_no_ns <- h12_h24_extract_act_pivot_ANOVAs_stats %>% 
  filter(!p_flag == 'ns') %>% 
  group_by(cluster) %>% 
  summarise(total = n()) 

h12_h24_extract_act_pivot_ANOVAs_stats_no_ns$cluster <- as.numeric(as.character(h12_h24_extract_act_pivot_ANOVAs_stats_no_ns$cluster))

h12_h24_extract_act_pivot_ANOVAs_stats_no_ns <- h12_h24_extract_act_pivot_ANOVAs_stats_no_ns %>% 
  dplyr::arrange(cluster) %>% 
  complete(.,cluster = 0:74, fill = list(total = 0)) %>% 
  write_csv("act_h12_24_ANOVA_by_R.csv")

# 1.15 UpSetR overlap plot----

h0_h15_extract_pivot_ANOVAs_filter_out_ns <- h0_h15_extract_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h9_h24_extract_pivot_ANOVAs_filter_out_ns <- h9_h24_extract_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h0_h15_extract_act_pivot_ANOVAs_filter_out_ns <- h0_h15_extract_act_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h9_h24_extract_act_pivot_ANOVAs_filter_out_ns <- h9_h24_extract_act_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

# filter out CCX2 and LBD11 from the h0_h15_extract_pivot_ANOVAs_filter_out_ns dataset due to tp 7.5 issue
h0_h15_extract_pivot_ANOVAs_filter_out_ns <- h0_h15_extract_pivot_ANOVAs_filter_out_ns %>% 
  filter(! Isoform == 'AT5G17850_P1' & ! Isoform == 'AT2G28500_P1')

myGeneSets <- list('repressed early' = h0_h15_extract_pivot_ANOVAs_filter_out_ns$Isoform,
                   'repressed late' = h9_h24_extract_pivot_ANOVAs_filter_out_ns$Isoform,
                   'activated early' = h0_h15_extract_act_pivot_ANOVAs_filter_out_ns$Isoform,
                   'activated late' = h9_h24_extract_act_pivot_ANOVAs_filter_out_ns$Isoform
)

h0_h15_extract_pivot_ANOVAs_filter_out_ns %>% filter(Isoform %in% 'AT5G17850_P1')

# fromList: a function to convert a list of named vectors to a data frame compatible with UpSetR
sets <- fromList(myGeneSets)

UpSet<-upset(sets, nsets=4, keep.order = T, sets = c("repressed early", "repressed late",
                                                     "activated early", "activated late"),
             sets.bar.color= c('#D53E4F', '#D53E4F', '#4DAF4A', '#4DAF4A'), matrix.color= 'grey30', point.size = 2.5, sets.x.label = "Consistent DE groups", mainbar.y.label = "Sizes of intersections")
UpSet
png("Upset_repressed_and_activated_early_late.png", width = 4, height = 3, units = 'in', res = 300)
UpSet

dev.off()

# 12 h no overlap----

h0_h12_extract_pivot_ANOVAs_filter_out_ns <- h0_h12_extract_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h12_h24_extract_pivot_ANOVAs_filter_out_ns <- h12_h24_extract_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h0_h12_extract_act_pivot_ANOVAs_filter_out_ns <- h0_h12_extract_act_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

h12_h24_extract_act_pivot_ANOVAs_filter_out_ns <- h12_h24_extract_act_pivot_ANOVAs %>% 
  filter(! p.value > 0.05)

# filter out CCX2 and LBD11 from the h0_h15_extract_pivot_ANOVAs_filter_out_ns dataset due to tp 7.5 issue
# h0_h15_extract_pivot_ANOVAs_filter_out_ns <- h0_h15_extract_pivot_ANOVAs_filter_out_ns %>% 
#   filter(! Isoform == 'AT5G17850_P1' & ! Isoform == 'AT2G28500_P1')

myGeneSets_new <- list('repressed h0-12' = h0_h12_extract_pivot_ANOVAs_filter_out_ns$Isoform,
                       'repressed h12-24' = h12_h24_extract_pivot_ANOVAs_filter_out_ns$Isoform,
                       'activated h0-12' = h0_h12_extract_act_pivot_ANOVAs_filter_out_ns$Isoform,
                       'activated h12-24' = h12_h24_extract_act_pivot_ANOVAs_filter_out_ns$Isoform
)

#h0_h12_extract_pivot_ANOVAs_filter_out_ns %>% filter(Isoform %in% 'AT5G17850_P1')

# fromList: a function to convert a list of named vectors to a data frame compatible with UpSetR
sets_new <- fromList(myGeneSets_new)

UpSet_new <-upset(sets_new, nsets=4, keep.order = T, sets = c("repressed h0-12", "repressed h12-24",
                                                              "activated h0-12", "activated h12-24"),
                  sets.bar.color= c('#D53E4F', '#D53E4F', '#4DAF4A', '#4DAF4A'), matrix.color= 'grey30', point.size = 2.5, sets.x.label = "Consistent DE groups", mainbar.y.label = "Sizes of intersections")
UpSet_new
png("Upset_repressed_and_activated_new.png", width = 4, height = 3, units = 'in', res = 300)
UpSet_new

dev.off()

# 1.16 Non-redundant isoform lists----
# *1.16.1 repressed----
repressed_distinct <- bind_rows(h0_h15_extract_pivot_ANOVAs_filter_out_ns,
                                h9_h24_extract_pivot_ANOVAs_filter_out_ns) %>% 
  distinct(Isoform) %>%
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) %>% 
  left_join(., TF_network_clusters) %>%
  mutate_all(as.character) %>%
  mutate(cluster = replace_na(cluster, 'nd')) %>% 
  write_csv("repressed_distinct.csv")

repressed_distinct_AGI <- repressed_distinct %>% 
  distinct(gene_ID) %>% 
  write_csv("repressed_distinct_AGI.csv")

Allan_list_AGI <- Allan_list_AGI %>% 
  rename('AGI' = 'gene_ID')

Allan_list_AGI_new <- Allan_list_AGI_new %>% 
  rename('AGI' = 'gene_ID')

overlap_isform_old <- anti_join(Allan_list_AGI_new, repressed_distinct_AGI, by = 'gene_ID')

# 12 h no overlap----

repressed_distinct_new <- bind_rows(h0_h12_extract_pivot_ANOVAs_filter_out_ns,
                                    h12_h24_extract_pivot_ANOVAs_filter_out_ns) %>% 
  distinct(Isoform) %>%
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) %>% 
  left_join(., TF_network_clusters) %>%
  mutate_all(as.character) %>%
  mutate(cluster = replace_na(cluster, 'nd')) %>% 
  write_csv("repressed_distinct_new.csv")

repressed_distinct_AGI_new <- repressed_distinct_new %>% 
  distinct(gene_ID) %>% 
  write_csv("repressed_distinct_AGI_new.csv")

Allan_list_AGI <- Allan_list_AGI %>% 
  rename('AGI' = 'gene_ID')

Allan_list_AGI_new <- Allan_list_AGI_new %>% 
  rename('AGI' = 'gene_ID')

overlap_isform_old <- anti_join(Allan_list_AGI_new, repressed_distinct_AGI_new, by = 'gene_ID')

# *1.16.2 activated----
activated_distinct <- bind_rows(h0_h15_extract_act_pivot_ANOVAs_filter_out_ns,
                                h9_h24_extract_act_pivot_ANOVAs_filter_out_ns) %>% 
  distinct(Isoform) %>%
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) %>% 
  left_join(., TF_network_clusters) %>%
  mutate_all(as.character) %>%
  mutate(cluster = replace_na(cluster, 'nd')) %>% 
  write_csv("activated_distinct.csv")

activated_distinct_AGI <- activated_distinct %>% 
  distinct(gene_ID) %>% 
  write_csv("activated_distinct_AGI.csv")

# 12 h no overlap----

activated_distinct_new <- bind_rows(h0_h12_extract_act_pivot_ANOVAs_filter_out_ns,
                                    h12_h24_extract_act_pivot_ANOVAs_filter_out_ns) %>% 
  distinct(Isoform) %>%
  mutate(gene_ID = substr(Isoform, start = 1, stop = 9)) %>% 
  left_join(., TF_network_clusters) %>%
  mutate_all(as.character) %>%
  mutate(cluster = replace_na(cluster, 'nd')) %>% 
  write_csv("activated_distinct_new.csv")

activated_distinct_AGI_new <- activated_distinct_new %>% 
  distinct(gene_ID) %>% 
  write_csv("activated_distinct_AGI_new.csv")

# 1.17 GO enrich plots----
At_genes_all<-read.csv('Gene TPM.csv')
all_genes<- At_genes_all[,1]

# *1.17.1 GO repressed----
GO_analysis <- enrichGO(gene = repressed_distinct_AGI$gene_ID,
                        universe = all_genes,
                        OrgDb = org.At.tair.db,
                        keyType = "TAIR",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10,
                        readable = TRUE,
                        pool = TRUE)

dotplot(GO_analysis, showCategory=10, font.size=12, label_format = 55)
dp_repressed <- dotplot(GO_analysis, showCategory=10, font.size=12, label_format = 55)

dp_repressed <-  dp_repressed + ggtitle("Repressed DE group")
ggsave('dp_repressed_GO.png', dp_repressed, height=3,width=9.5,units="in",dpi=200)

# 12 h no overlap----

GO_analysis_rep_new <- enrichGO(gene = repressed_distinct_AGI_new$gene_ID,
                                universe = all_genes,
                                OrgDb = org.At.tair.db,
                                keyType = "TAIR",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.10,
                                readable = TRUE,
                                pool = TRUE)

dotplot(GO_analysis_rep_new, showCategory=10, font.size=12, label_format = 55)
dp_repressed_new <- dotplot(GO_analysis_rep_new, showCategory=10, font.size=12, label_format = 55)

dp_repressed_new <-  dp_repressed_new + ggtitle("Repressed DE group")
ggsave('dp_repressed_GO_new.png', dp_repressed_new, height=3,width=9.5,units="in",dpi=200)

# *1.17.2 GO activated----
GO_analysis_act <- enrichGO(gene = activated_distinct_AGI$gene_ID,
                            universe = all_genes,
                            OrgDb = org.At.tair.db,
                            keyType = "TAIR",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.10,
                            readable = TRUE,
                            pool = TRUE)

dotplot(GO_analysis_act, showCategory=10, font.size=12, label_format = 55)
dp_activated <- dotplot(GO_analysis_act, showCategory=10, font.size=12, label_format = 55)

dp_activated <-  dp_activated + ggtitle("Activated DE group")
ggsave('dp_activated_GO.png', dp_activated, height=3,width=9.5,units="in",dpi=200)

# 12 h no overlap----

GO_analysis_act_new <- enrichGO(gene = activated_distinct_AGI_new$gene_ID,
                                universe = all_genes,
                                OrgDb = org.At.tair.db,
                                keyType = "TAIR",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.10,
                                readable = TRUE,
                                pool = TRUE)

dotplot(GO_analysis_act_new, showCategory=10, font.size=12, label_format = 55)
dp_activated_new <- dotplot(GO_analysis_act_new, showCategory=10, font.size=12, label_format = 55)

dp_activated_new <-  dp_activated_new + ggtitle("Activated DE group")
ggsave('dp_activated_GO_new.png', dp_activated, height=3,width=9.5,units="in",dpi=200)



GO_rep_act <- dp_repressed / dp_activated 

GO_rep_act_horizontal <- dp_repressed + dp_activated

GO_rep_act_horizontal_new <- dp_repressed_new + dp_activated_new

ggsave('GO_rep_act_horiz_new.png', GO_rep_act_horizontal_new, height=4,width=18,units="in",dpi=200)

# 1.18 Odd plots on h0-24 repressed and activated sets----
# 1.18.1 repressed----
# associate with a TFcluster number
inner_join_clusters_repressed <- inner_join(TF_network_clusters, repressed_distinct_AGI)

inner_join_clusters_repressed_count <- inner_join_clusters_repressed %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

repressed_fishers_prep <- inner_join_clusters_repressed_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

repressed_fishers <- repressed_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(repressed_fishers)[5] <- "Odds_Ratio"

repressed_fishers_p <- repressed_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(repressed_fishers_p)[5] <- "p_value"

repressed_fishers_Odds <- repressed_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

repressed_fishers_pval <- repressed_fishers_p %>% 
  dplyr::select(p_value)

repressed_fishers_Odds_pval <- inner_join_clusters_repressed_count %>% 
  bind_cols(repressed_fishers_Odds, repressed_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(repressed_fishers_Odds_pval, 'repressed_fishers_Odds_pval.csv')

repressed_fishers_Odds_pval_heatmap1 <- repressed_fishers_Odds_pval %>%
  mutate(time_points = "repressed") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

repressed_fishers_Odds_pval_heatmap2 <- repressed_fishers_Odds_pval_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  #ggplot(aes(x = cluster, y= time_points, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="Cluster") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 14),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 14),
        legend.key.height = grid::unit(0.8,"cm"),
        legend.key.width = grid::unit(0.8,"cm"),
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14, vjust = 0.2,colour = 'grey30'),
        axis.title.y = element_text(size = 14, colour = 'grey30'),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank())
#       #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
#       plot.title = element_text(colour = 'grey30', size = 10, face = "bold", hjust = 1)) + 
# ggtitle('Repressed in Col-0') 

repressed_fishers_Odds_pval_heatmap2

# 1.18.2 activated----

inner_join_clusters_activated <- inner_join(TF_network_clusters, activated_distinct_AGI)

inner_join_clusters_activated_count <- inner_join_clusters_activated %>%
  group_by(cluster) %>%
  count(cluster) %>%
  ungroup() %>% 
  mutate(percent = round((n/sum(n)*100), 1)) %>%
  complete(.,cluster = 0:74, fill = list(n = 0, percent = 0)) %>% 
  inner_join(., TF_network_clusters_count) 

activated_fishers_prep <- inner_join_clusters_activated_count %>% 
  mutate(fishers_col2 = sum(n) - n,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(n, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = n,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

activated_fishers <- activated_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))

colnames(activated_fishers)[5] <- "Odds_Ratio"

activated_fishers_p <- activated_fishers_prep %>% 
  data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$p.value)) 

colnames(activated_fishers_p)[5] <- "p_value"

activated_fishers_Odds <- activated_fishers %>% 
  dplyr::select(Odds_Ratio) %>% 
  round(2)

activated_fishers_pval <- activated_fishers_p %>% 
  dplyr::select(p_value)

activated_fishers_Odds_pval <- inner_join_clusters_activated_count %>% 
  bind_cols(activated_fishers_Odds, activated_fishers_pval) %>% 
  dplyr::rename(number_in_cluster = n)

write_csv(activated_fishers_Odds_pval, 'activated_fishers_Odds_pval.csv')

activated_fishers_Odds_pval_heatmap1 <- activated_fishers_Odds_pval %>%
  mutate(time_points = "activated") %>% 
  dplyr::select(cluster, Odds_Ratio, time_points) %>% 
  mutate(cluster=factor(cluster),
         time_points =factor(time_points),
         Enrichfactor=cut(Odds_Ratio, breaks=c(-1, 1, 3, 5, 7, max(Odds_Ratio)), labels=c("0-1","1-3","3-5","5-7",">7")),
         Enrichfactor=factor(as.character(Enrichfactor),levels=rev(levels(Enrichfactor))))

activated_fishers_Odds_pval_heatmap2 <- activated_fishers_Odds_pval_heatmap1 %>% 
  ggplot(aes(x = time_points, y= cluster, fill= Enrichfactor))+
  geom_tile(colour="grey30",size=0.3, width = 2.5) +
  guides(fill=guide_legend(title="Enrichment")) +
  labs(x="",y="") +
  scale_y_discrete(expand=c(0,0), breaks=c("0","5","10","15","20","25","30","35","40","45","50","55","60","65","70")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=rev(c('#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#810f7c'))) +
  coord_fixed(ratio=0.75) +
  theme_grey() +
  theme(legend.position = "right",legend.direction = "vertical",
        legend.title = element_text(colour = 'grey30', size = 14),
        legend.margin = margin(grid::unit(0,"cm")),
        legend.text = element_text(colour = 'grey30',size = 14),
        legend.key.height = grid::unit(0.8,"cm"),
        legend.key.width = grid::unit(0.8,"cm"),
        axis.text.x = element_text(size = 14,colour = 'grey30',angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        plot.background = element_blank(),
        panel.border = element_blank(),
        #plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        #plot.title = element_text(colour = 'grey30', size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

activated_fishers_Odds_pval_heatmap2

patch_plot3 <- TF_cluster_sizes_plot + (repressed_fishers_Odds_pval_heatmap2 + activated_fishers_Odds_pval_heatmap2 + plot_layout(ncol = 2) + plot_layout(guides = "collect"))
patch_plot3

patch_plot4 <- TF_cluster_sizes_plot + repressed_fishers_Odds_pval_heatmap2 + activated_fishers_Odds_pval_heatmap2 + plot_layout(guides = "collect")
patch_plot4

ggsave(patch_plot4, filename="patch_plot4.png", height = 9, width = 5, units = "in", dpi = 300)

# 1.19 z-score plots on h0-24 repressed and activated sets--------
# *1.19.1 cluster 9----

# get cluster 9 genes
# this is just gene_IDs
cluster9_repressed <- inner_join_clusters_repressed %>% 
  filter(cluster %in% '9')

repressed_distinct2 <- repressed_distinct %>% 
  dplyr::select(-cluster)

# get the associated isoforms
cluster9_genes_isoforms_represssed <- repressed_distinct2 %>% 
  inner_join(cluster9_repressed, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster9_genes_isoforms_represssed, 'cluster9_genes_isoforms_represssed.csv')

# write the just gene_IDs to file
write_csv(cluster9_repressed, 'cluster9_repressed.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster9_gene_ids_repressed <- read_csv('cluster9_repressed_with_gene_IDs.csv')

cluster9_repressed_with_gene_IDs <- cluster9_gene_ids_repressed %>%
  inner_join(cluster9_repressed, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster9_repressed_with_gene_IDs, 'cluster9_repressed_with_gene_IDs.csv')

#cluster9_repressed_isoforms
cluster9_repressed_isoforms_z_scores <- cluster9_genes_isoforms_represssed %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 9', subtitle = '98 isoforms (82 gene loci)')

ggsave('cluster9_repressed_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster9_repressed_isoforms_z_scores


# *1.19.2 cluster 20----

# get cluster 20 genes
# this is just gene_IDs
cluster20_repressed <- inner_join_clusters_repressed %>% 
  filter(cluster %in% '20')

repressed_distinct2 <- repressed_distinct %>% 
  dplyr::select(-cluster)

# get the associated isoforms
cluster20_genes_isoforms_represssed <- repressed_distinct2 %>% 
  inner_join(cluster20_repressed, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster20_genes_isoforms_represssed, 'cluster20_genes_isoforms_represssed.csv')

# write the just gene_IDs to file
write_csv(cluster20_repressed, 'cluster20_repressed.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster20_gene_ids_repressed <- read_csv('cluster20_repressed_with_gene_IDs.csv')

cluster20_repressed_with_gene_IDs <- cluster20_gene_ids_repressed %>%
  inner_join(cluster20_repressed, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster20_repressed_with_gene_IDs, 'cluster20_repressed_with_gene_IDs.csv')

#cluster20_repressed_isoforms
cluster20_repressed_isoforms_z_scores <- cluster20_genes_isoforms_represssed %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 20', subtitle = '53 isoforms (45 gene loci)')

ggsave('cluster20_repressed_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster20_repressed_isoforms_z_scores

# *1.19.3 cluster 25----
# get cluster 25 genes
# this is just gene_IDs
cluster25_repressed <- inner_join_clusters_repressed %>% 
  filter(cluster %in% '25')

repressed_distinct2 <- repressed_distinct %>% 
  dplyr::select(-cluster)

# get the associated isoforms
cluster25_genes_isoforms_represssed <- repressed_distinct2 %>% 
  inner_join(cluster25_repressed, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster25_genes_isoforms_represssed, 'cluster25_genes_isoforms_represssed.csv')

# write the just gene_IDs to file
write_csv(cluster25_repressed, 'cluster25_repressed.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster25_gene_ids_repressed <- read_csv('cluster25_repressed_with_gene_IDs.csv')

cluster25_repressed_with_gene_IDs <- cluster25_gene_ids_repressed %>%
  inner_join(cluster25_repressed, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster25_repressed_with_gene_IDs, 'cluster25_repressed_with_gene_IDs.csv')

#cluster25_repressed_isoforms
cluster25_repressed_isoforms_z_scores <- cluster25_genes_isoforms_represssed %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>% 
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#DA777E"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 25', subtitle = '37 isoforms (32 gene loci)')

ggsave('cluster25_repressed_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster25_repressed_isoforms_z_scores 

# *1.19.4 cluster 17----
# get cluster 17 genes
# this is just gene_IDs
cluster17_activated <- inner_join_clusters_activated %>% 
  filter(cluster %in% '17')

activated_distinct2 <- activated_distinct %>% 
  dplyr::select(-cluster)

# get the associated isoforms
cluster17_genes_isoforms_activated <- activated_distinct2 %>% 
  inner_join(cluster17_activated, by= 'gene_ID') %>% 
  dplyr::select(-gene_ID)

write_csv(cluster17_genes_isoforms_activated, 'cluster17_genes_isoforms_activated.csv')

# write the just gene_IDs to file
write_csv(cluster17_activated, 'cluster17_activated.csv')

# get gene descriptions from TAIR add to the 'genes_ID' file then read back in
cluster17_gene_ids_activated <- read_csv('cluster17_activated_with_gene_IDs.csv')

cluster17_activated_with_gene_IDs <- cluster17_gene_ids_activated %>%
  inner_join(cluster17_activated, by= 'gene_ID') %>%
  dplyr::select(-c(2, 4, 7))

#colnames(cluster9_h0_h15_rep_with_gene_IDs)[2] <- "Isoform"

write_csv(cluster17_activated_with_gene_IDs, 'cluster17_activated_with_gene_IDs.csv')

#cluster17_activated_isoforms

cluster17_activated_isoforms_z_scores <- cluster17_genes_isoforms_activated %>% 
  inner_join(transcript_RNAseq_means_select_filtered_low, by = 'Isoform') %>%
  dplyr::select(-(c(2,23))) %>% 
  pivot_longer(cols= Col_T9:rve2_T18,
               names_to='time_point',
               values_to='TPM') %>% 
  group_by(Isoform) %>% 
  mutate(z_score = scale(TPM)) %>% 
  ungroup() %>%
  group_by(time_point) %>% 
  summarise(mean=mean(z_score), sd=sd(z_score)) %>% 
  mutate(time = case_when(grepl('_T9', time_point) ~ 0,
                          grepl('_T10', time_point) ~ 3,
                          grepl('_T11', time_point) ~ 6,
                          grepl('_T12', time_point) ~ 7.5,
                          grepl('_T13', time_point) ~ 9,
                          grepl('_T14', time_point) ~ 12,
                          grepl('_T15', time_point) ~ 15,
                          grepl('_T16', time_point) ~ 18,
                          grepl('_T17', time_point) ~ 21,
                          TRUE ~ 24),
         genotype = case_when(grepl('Col', time_point) ~ 'Col',
                              TRUE ~ 'rve')) %>% 
  arrange(time) %>%
  ggplot(aes(time, mean)) +          
  geom_line(aes(color = genotype), size = 1.5) + 
  scale_color_manual(values = c("grey30", "#B3E3A0"), labels=c('Col-0', 'rve2-2')) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  theme_light(base_family = 'Arial',
              base_size = 14) +
  xlab('h after cooling (day 2)') +
  ylab('mean z-score') +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_rect(linewidth=0.8, linetype="solid", 
                                         colour ="grey30"),
        legend.title = element_blank()) +
  ggtitle('Cluster 17', subtitle = '31 isoforms (31 gene loci)')

ggsave('cluster17_activated_isoforms_z_scores.png', height = 5, width = 4, units = 'in')

cluster17_activated_isoforms_z_scores