
getwd()


orig_data <- read.delim("Suppl_Figure6.cdt",
                        header=FALSE,
                        as.is=TRUE)

expressions <- orig_data[4:nrow(orig_data), 5:ncol(orig_data)]

mat <- matrix(as.numeric(unlist(expressions)),
              nrow=nrow(expressions), ncol=ncol(expressions))

dimnames(mat) <- list(orig_data[4:nrow(orig_data),3],
                      orig_data[1, 5:ncol(orig_data)])


mat.thresh <- mat[rep(seq_len(nrow(mat)), each = 2), ]
indices <- seq(2, nrow(mat.thresh), by=2)
rownames(mat.thresh)[indices] <- lapply(rownames(mat.thresh)[indices], function(x) paste0("dup", x))


express.thresh <- 4.0
mat.thresh[mat.thresh <= express.thresh & mat.thresh >= -express.thresh] <- NA

indices <- which(mat.thresh < 0, arr.ind = TRUE)
mat.thresh[indices[seq(1, nrow(indices), by=2), ] ] <- NA

indices <- which(mat.thresh > 0, arr.ind = TRUE)
mat.thresh[indices[seq(2, nrow(indices), by=2), ] ] <- NA


# First add all the gene names from the rows as new columns (on the right)
extra <- matrix(ncol = nrow(mat.thresh), nrow = nrow(mat.thresh))
colnames(extra) <- rownames(mat.thresh)
mat.sq <- cbind(mat.thresh, extra)

# Now add all the array names from the columns as new rows (at the top)
extra <- matrix(ncol = ncol(mat.sq), nrow = ncol(mat.thresh))
rownames(extra) <- colnames(mat.thresh)
mat.sq <- rbind(extra, mat.sq)


mat.sq[is.na(mat.sq)] <- 0
mat.sq <- abs(mat.sq)


require(HiveR)
f6hive <- adj2HPD(M = mat.sq,
                  axis.cols = c("white"),
                  type = "2D",
                  desc = "the array and gene expression data")


str(f6hive)


f6hive$nodes$axis <- c(rep(as.integer(1), ncol(mat)),
                       rep(c(as.integer(2), as.integer(3)), nrow(mat)))


f6hive <- manipAxis(f6hive, "rank")


plotHive(f6hive, bkgnd = "white")


f6hive$edges$weight <- f6hive$edges$weight / max(f6hive$edges$weight)
f6hive$nodes$size <- 0.4
f6hive <- manipAxis(f6hive, "scale", action = c(3.75, 1.0, 1.0))
plotHive(f6hive, bkgnd = "white")


f6hive$edges$color <- ifelse(f6hive$edges$id1 %% 2 == 0, "green", "red")
plotHive(f6hive, bkgnd = "white")


f6hive$nodes$color <- "gray60"
subclass.colors <- read.csv("subtype_colors.csv",
                            header=TRUE, stringsAsFactors = FALSE)
replace <- subclass.colors$color[match(f6hive$nodes$lab,
                                       subclass.colors$array)]
f6hive$nodes$color[!is.na(replace)] <- replace[!is.na(replace)]
plotHive(f6hive, bkgnd = "white")


f6hive$edges$color <- f6hive$nodes$color[match(f6hive$edges$id2, f6hive$nodes$id)]
plotHive(f6hive, bkgnd = "white")
