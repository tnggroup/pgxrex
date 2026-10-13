### Load packages
library(data.table)
library(dplyr)
library(tidyr)

### Load data - before loading, set blanks to NA and removed duplicated rows
PRADASim <- fread("Desktop/Joni/Senior_Lecturer/Grants/Successful/PRADA/final_predicted_score_petrushka_arm.scorediff.txt", data.table=F)

### Make unique ID
PRADASim$uid <- paste(PRADASim$pid, PRADASim$gid, sep="_")

### Make relative gscore as difference from max
PRADASim$gscore_rel <- NA
PRADASim$gscore_rel <- (PRADASim$gscore / PRADASim$gscore.max) - 1

### Scale relative gscore to approximate max scale from pheno score
PRADASim$gscore_rel_scaled_07 <- PRADASim$gscore_rel * 0.07

### Make combined scores 
### Note that missing pheno scores are NA, which is fine
### Note that missing geno scores are max, which is up for discussion
PRADASim$combscore <- PRADASim$score + PRADASim$gscore_rel
PRADASim$combscore_scaled_07 <- PRADASim$score + PRADASim$gscore_rel_scaled_07

### Rank within uid
PRADASim <- PRADASim %>%
  group_by(uid) %>%
  mutate(combscore_rank = min_rank(desc(combscore))) %>%
  ungroup()

PRADASim <- PRADASim %>%
  group_by(uid) %>%
  mutate(combscore_scaled_07_rank = min_rank(desc(combscore_scaled_07))) %>%
  ungroup()

### Check if top 3 differ

top3_check <- PRADASim %>%
  group_by(uid) %>%
  summarise(
    top_rank = list(drug_name[order(rank)][1:3]),
    top_comb = list(drug_name[order(combscore_rank)][1:3]),
    same_top3 = setequal(top_rank[[1]], top_comb[[1]]),
    n_overlap = length(intersect(top_rank[[1]], top_comb[[1]])),
    .groups = "drop"
  )

top3_check <- top3_check %>%
  separate(uid, into = c("pid", "gid"), sep = "_")

table(top3_check$same_top3)
## FALSE  TRUE 
## 14644  6914 

top3_check_by_pid <- top3_check %>%
  group_by(pid) %>%
  summarise(
    n = n(),
    n_same = sum(same_top3),
    prop_same = mean(same_top3),
    .groups = "drop"
  )

quantile(top3_check_by_pid$prop_same, probs = seq(0,1,0.1))
# 0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
# 0.0000000 0.1686747 0.2054859 0.2656130 0.2988506 0.3179144 0.3372093 0.3602882 0.3750000 0.4137931 1.0000000 

top3_check_scaled_07 <- PRADASim %>%
  group_by(uid) %>%
  summarise(
    top_rank = list(drug_name[order(rank)][1:3]),
    top_comb = list(drug_name[order(combscore_scaled_07_rank)][1:3]),
    same_top3 = setequal(top_rank[[1]], top_comb[[1]]),
    n_overlap = length(intersect(top_rank[[1]], top_comb[[1]])),
    .groups = "drop"
  )

top3_check_scaled_07 <- top3_check_scaled_07 %>%
  separate(uid, into = c("pid", "gid"), sep = "_")

table(top3_check_scaled_07$same_top3)
## FALSE  TRUE 
## 11907  9651 

top3_check_by_pid_scaled_07 <- top3_check_scaled_07 %>%
  group_by(pid) %>%
  summarise(
    n = n(),
    n_same = sum(same_top3),
    prop_same = mean(same_top3),
    .groups = "drop"
  )

quantile(top3_check_by_pid_scaled_07$prop_same, probs = seq(0,1,0.1))
# 0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
# 0.0000000 0.1741318 0.2195322 0.2910397 0.3318352 0.3662602 0.4100313 0.4818449 0.7704545 0.9185897 1.0000000 

