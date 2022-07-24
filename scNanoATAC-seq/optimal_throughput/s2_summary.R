library(ggplot2)
library(dplyr)
library(RColorBrewer)

result <- list()
for (tp in (2:10) * 500) {
  for (seed in 1:15) {
    try({
      tb <-
        readRDS(paste0(
          './archr/', tp, '.', seed, '/', tp, '.', seed, '.clustering.Rds'
        ))
      if (nrow(tb)==5) {
        precision <- c()
        for (i in c("K562", "HEK293T", "HFF1", "GM12878", "eHAP1")) {
          precision[i] <- max(tb[, i]) / sum(tb[, i])
        }
        result[[paste0(tp, seed)]] <- data.frame(
          precision = precision,
          cell.type = names(precision),
          throughput = tp
      )
      }
    })
  }
}

df <- do.call(rbind,result)
rownames(df) <- NULL
df$throughput <- factor(df$throughput)

df.sd <- df %>%
    group_by(cell.type,throughput) %>%
    summarise_each(list(mean=mean,sd=sd)) %>%
    mutate(
        precision=mean,
        mean=NULL,
        upper=sapply(precision+sd,function(x)min(1,x)),
        lower=sapply(precision-sd,function(x)max(0,x))
    )

pdf('dowmsapled_clustering_precision.pdf',width=6,height=6)

cp=colorRampPalette(rev(brewer.pal(9,'Blues')))(18)[seq_along(unique(df.sd$throughput))]

ggplot(data=df.sd,aes(x=cell.type, y=precision)) +
  geom_pointrange(aes(ymin=lower, ymax=upper,color=throughput),
                 position=position_dodge(.9)) +
    scale_color_manual(values=cp) +
    theme_classic() +
    ylim(c(0,1))


cp=colorRampPalette(rev(brewer.pal(9,'Reds')))(18)[seq_along(unique(df.sd$throughput))]

ggplot(data=df.sd,aes(x=cell.type, y=precision)) +
  geom_pointrange(aes(ymin=lower, ymax=upper,color=throughput),
                 position=position_dodge(.9)) +
    scale_color_manual(values=cp) +
    theme_classic() +
    ylim(c(0,1))

dev.off()

