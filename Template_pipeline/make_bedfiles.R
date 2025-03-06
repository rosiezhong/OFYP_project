library(dplyr)

var <- commandArgs(trailingOnly = TRUE)

for (i in 1:2) {
    depth<-data.frame(read.table(var[i]))
    depth<-depth %>% dplyr::rename(dp = V3)
    depth<-depth %>% dplyr::rename(pos = V2)
    depth<-depth %>% dplyr::rename(chr = V1)

    pos_start<-c()
    pos_end<-c()
    chr_start<-c()
    chr_end<-c()

    if (depth$dp[1]==2) { 
        pos_start<-append(pos_start,depth$pos[1])
        chr_start<-append(chr_start,depth$chr[1])
    }

    for (j in 2:length(depth$dp)) { 
        if (depth$dp[j]==2 & depth$dp[j-1] != 2) {
            pos_start<-append(pos_start,depth$pos[j])
            chr_start<-append(chr_start,depth$chr[j])
        }
    }

    for (j in 1:(length(depth$dp)-1)) { 
        if (depth$dp[j]==2 & depth$dp[j+1] != 2) {
            pos_end<-append(pos_end,depth$pos[j])
            chr_end<-append(chr_end,depth$chr[j])
        }
    }

    if (depth$dp[length(depth$dp)]==2) {
        pos_end<-append(pos_end,depth$pos[length(depth$dp)])
        chr_end<-append(chr_end,depth$chr[length(depth$dp)])
    }
    
    if (identical(chr_start,chr_end)==TRUE) { 
        positions_name <- paste0("positions", i)
        assign(positions_name, data.frame(chr_start, pos_start, pos_end))
    }
}

print("Start and end of 2 strings extracted.")

merged_pos<-rbind(positions1,positions2)
merged_pos<-unique(merged_pos %>% arrange(pos_start))

Counter<-1
merged_pos$counter <- vector("list", nrow(merged_pos))

print("Proceeding to Overlap phase 1.")

for (i in 1:(nrow(merged_pos)-1)) {
  merged_pos$counter[[i]]<-c(merged_pos$counter[[i]],Counter)
  for (j in (i+1):nrow(merged_pos)) {
    if (merged_pos$pos_start[j]<=merged_pos$pos_end[i]) {
      merged_pos$counter[[j]]<-c(merged_pos$counter[[j]],Counter)
    }
  }
  Counter<-Counter + 1
}

print("Overlap phase 1 completed.")

merged_pos$counter2 <- NA
label <- 1 
print("Proceeding to Overlap phase 2.")

counter_index <- list()

for (i in 1:nrow(merged_pos)) {
  for (val in merged_pos$counter[[i]]) {
    if (!val %in% names(counter_index)) {
      counter_index[[as.character(val)]] <- integer(0)
    }
    counter_index[[as.character(val)]] <- c(counter_index[[as.character(val)]], i)
  }
}

for (i in 1:nrow(merged_pos)) {
  if (is.na(merged_pos$counter2[i])) {
    queue <- c(i)
    merged_pos$counter2[i] <- label
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      shared_rows <- unique(unlist(lapply(merged_pos$counter[[current]], function(val) counter_index[[as.character(val)]])))
      for (j in shared_rows) {
        if (i != j && is.na(merged_pos$counter2[j])) {
          merged_pos$counter2[j] <- label
          queue <- c(queue, j)
        }
      }
    }
    label <- label + 1
  }
}

print("Overlap phase 2 completed.")
print("Proceeding to Overlap phase 3(grouping).")

grouped_pos <- merged_pos %>%
  group_by(counter2) %>%
  summarize(
    pos_start = min(pos_start),
    pos_end = max(pos_end),
    chr_start = if(length(unique(chr_start)) > 1) {
      stop("Error: Non-unique chr_pos values within a group")
    } else {
      unique(chr_start)
    },
    .groups = 'drop'
  ) %>%
  arrange(pos_start)

print("Overlap phase 3(grouping) completed.")
print("Proceeding to extract vcfs.")

library(vcfR)
my.vcf<-read.vcfR(var[3])
vcf_fix<-data.frame(my.vcf@fix)
vcf_fix$POS<-as.numeric(vcf_fix$POS)

library(data.table)
vcf_dt <- as.data.table(vcf_fix)
grouped_dt <- as.data.table(grouped_pos)
setkey(grouped_dt, pos_start, pos_end)
vcf_range <- vcf_dt[, .(start = POS, end = POS)]
final_pos <- foverlaps(vcf_range, grouped_dt, by.x = c("start", "end"), type = "within")
final_pos <- final_pos[!is.na(pos_start)]
final_pos <- final_pos %>%
  select(-start, -end)
print("vcf regions extracted.")
final_pos<-unique(final_pos)
final_pos<-final_pos%>%arrange(pos_start)

print("Proceeding to create bed files.")

for (i in 1:nrow(final_pos)) {
  bed_file_name <- paste(var[4], i, ".bed", sep="")
  
  bed_file <- file(bed_file_name, "w")
  
  line <- paste(final_pos$chr_start[i], final_pos$pos_start[i], final_pos$pos_end[i], sep="\t")
  cat(line, file = bed_file, sep="")
  
  close(bed_file)
}
