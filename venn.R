install.packages("VennDiagram")
library(VennDiagram)

# ven diagram
degs.prim <- rownames(filtered.prim)
degs.prim.up <- rownames(dplyr::filter(filtered.prim, filtered.prim$log2FoldChange > 1))
degs.prim.down <- rownames(dplyr::filter(filtered.prim, filtered.prim$log2FoldChange < -1))
degs.recu <- rownames(filtered.recu)
degs.recu.up <- rownames(dplyr::filter(filtered.recu, filtered.recu$log2FoldChange > 1))
degs.recu.down <- rownames(dplyr::filter(filtered.recu, filtered.recu$log2FoldChange < -1))

dev.off()
draw.pairwise.venn(length(degs.prim),
                   length(degs.recu),
                   length( intersect(degs.prim, degs.recu) ),
                   category = c("Primary vs Normal", "Recurrence vs Normal"), scaled = F,
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                   cat.pos = c(0, 0))
#grid.draw(venn.plot)
dev.off()
draw.pairwise.venn(length(degs.prim.up),
                   length(degs.recu.up),
                   length( intersect(degs.prim.up, degs.recu.up)),
                   #category = c("Primary vs Normal", "Recurrence vs Normal"), 
                   scaled = F,
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                   #fontfamily = rep("arial", 15),
                   cat.pos = c(0, 0))

dev.off()
draw.pairwise.venn(length(degs.prim.down),
                   length(degs.recu.down),
                   length( intersect(degs.prim.down, degs.recu.down)),
                   #category = c("Primary vs Normal", "Recurrence vs Normal"), 
                   scaled = T,
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                   #fontfamily = rep("arial", 15),
                   cat.pos = c(0, 0))


degs.prim.spec <- setdiff(filtered.prim$genes, filtered.recu$genes)
degs.recu.spec <- setdiff(filtered.recu$genes, filtered.prim$genes)
degs.prim.spec2 <- setdiff(filtered.prim2$genes, filtered.recu2$genes)
degs.recu.spec2 <- setdiff(filtered.recu2$genes, filtered.prim2$genes)

draw.triple.venn(length(degs.prim),
                 length(degs.recu),
                 length(degs.deseq.recu.prim),
                 length(intersect(degs.prim, degs.recu)),
                 length(intersect(degs.deseq.recu.prim, degs.recu)),
                 length(intersect(degs.prim, degs.deseq.recu.prim)),
                 length(intersect(degs.deseq.recu.prim, intersect(degs.prim, degs.recu))),
                 category = c("Primary vs Normal", "Recurrence vs Normal", "DESeq"), 
                 scaled = T,
                 fill = c("light blue", "pink", "light yellow"), 
                 #alpha = rep(0.5, 2),
                 #fontfamily = rep("arial", 15),
)