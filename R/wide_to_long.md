```R
pseudocount.use = 1
log1pdata.mean.fxn <- function(x) {
  return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
}
scaledata.mean.fxn <- rowMeans
counts.mean.fxn <- function(x) {
  return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
}


degs = cosg(sobj, layer = 'counts', assay = 'RNA', mu = 1.3, remove_lowly_expressed = T, expressed_pct = .1, n_genes_user = 200)

degs.long = cbind(tidyr::pivot_longer(degs$names, cols = everything(), names_to = 'cluster', values_to = 'gene'), tidyr::pivot_longer(degs$scores, cols = everything(), names_to = 'cluster', values_to = 'score')[,-1]) %>% arrange(cluster, desc(score)) %>% split(f = .[['cluster']]) %>% plapply(function(x){
    g = unique(x$cluster)
    gg = unique(x$gene)
    fc = FoldChange(sobj, ident.1 = g, slot = 'counts', features = gg, mean.fxn = scaledata.mean.fxn, base = 2, pseudocount.use = 1) %>% mutate('gene' = rownames(.))
    x = x %>% left_join(fc, by = 'gene')
    x
}, cl = 20) %>% Reduce(f = rbind)
```