library(ggplot2)
library(ggrepel)

# 1. 先在完整 object 上算 centrality（不要 subset，避免 subsetCellChat 的 bug）
object.list <- lapply(object.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})

# 2. 提取每个 object 的 outgoing/incoming strength + interaction count
extractSignalingRoleDF <- function(object) {
  centr <- object@netP$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  num.link <- rowSums(object@net$count) + colSums(object@net$count) - diag(object@net$count)
  data.frame(
    celltype = rownames(outgoing),
    x = rowSums(outgoing),
    y = rowSums(incoming),
    Count = num.link
  )
}

# 3. 合并所有 condition 成一个长表，并打上 condition 标签
df.all <- do.call(rbind, lapply(names(object.list), function(nm) {
  df <- extractSignalingRoleDF(object.list[[nm]])
  df$condition <- nm
  df
}))

# 4. 只保留你要的 cell type
celltypes.use <- c("CellTypeA", "CellTypeB", "CellTypeC")
df.sel <- df.all[df.all$celltype %in% celltypes.use, ]

# 5. 点大小范围（用全体数据算，保证跨图/跨条件可比；如果只想按选中的算，用 df.sel 代替 df.all）
weight.MinMax <- c(min(df.all$Count), max(df.all$Count))

# 6. 画在同一张图上，颜色映射 condition
ggplot(df.sel, aes(x = x, y = y, colour = condition)) +
  geom_point(aes(size = Count), alpha = 0.75) +
  scale_size_continuous(range = c(3, 10), limits = weight.MinMax, name = "Strength") +
  geom_text_repel(aes(label = celltype), size = 3, max.overlaps = Inf,
                   show.legend = FALSE) +
  theme_bw() +
  labs(x = "Outgoing interaction strength",
       y = "Incoming interaction strength",
       colour = "Condition") +
  theme(legend.position = "right")
