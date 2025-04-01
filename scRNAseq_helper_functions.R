#========================================================================================================================
## Correlation plot
##========================================================================================================================
plot.cor <- function(object = final.pop.call.integrated.full.seurat, assay = "RNA", features = "GAPDH", cor.feature = "CLEC4A", group.by = "ident", ncol = 2, file.name = "cor.plot.pdf", width = 10, height = 10){
  # Get average expression for the requested features per group
  avg.expression.gene.set <- as.data.frame(AverageExpression(object   = object, 
                                                             assays   = assay, 
                                                             features = c(features, cor.feature), 
                                                             group.by = group.by)[[assay]])

  # Scale from 0 to 1
  avg.expression.gene.set <- as.data.frame(apply(avg.expression.gene.set, 1, rescale))
  
  
  # Iterate over all comparisons with CLEC4A
  ggscatter(data           = avg.expression.gene.set,
            x              = cor.feature,
            y              = colnames(avg.expression.gene.set)[-length(colnames(avg.expression.gene.set))],
            combine        = T,
            ncol           = ncol,
            add            = "reg.line",
            add.params     = list(color = "red", fill = "grey"),
            conf.int       = TRUE, 
            cor.coef       = TRUE, 
            cor.method     = "spearman",
            xlab           = paste(cor.feature,"(AU)", sep = " "),
            ylab           = "(AU)", 
            cor.coef.coord = c(0,1.1), 
            ylim           = c(0,1.1),
            xlim           = c(0,1),
  ) +
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    theme_pubr(border = T)+
    theme(axis.title = element_text(face = "italic"),
          aspect.ratio = 1)
  ggsave(file.name, width = width, height = height)
}

#=========================================================================================================================
## Plot a series of customised plots for stratified expression
##========================================================================================================================
stratPlots <- function(object = samples.seurat, assay = "RNA", group.by = NULL, features = "GAPDH", cols = NULL, x.label = "Stratified", name = "plot_basename", ncol = 1){
   p1 <- VlnPlot(object  = object,
                 assay   = "RNA",  
                features = features, 
                group.by = group.by, 
                cols     = cols, 
                ncol     = ncol, combine = T) & theme(axis.title.x = element_blank()) 
  p1 <- p1 + plot_annotation(caption = x.label, theme = theme(plot.caption = element_text(hjust = 0.5)))
  
  ggsave(paste(name, " - violin plot.pdf", sep = ""), height = 5*(length(features)/ncol), plot = p1)

  DotPlot(object = object, assay = "RNA", features = features, group.by = group.by, cols = "Spectral", dot.scale = 15) + scale_colour_gsea()
  ggsave(paste(name, " - dot plot.pdf", sep = ""), width = 2*length(features))
  
  FeaturePlot(object = object, assay = "RNA", features = features, pt.size = 2, split.by = group.by, cols = c("grey", "firebrick"), order = T, ncol = ncol)
  ggsave(paste(name, " - feature plot.pdf", sep = ""), height = 5*(length(features)/ncol), width = 5*(length(features)/ncol))
}

#==================================================================================================================================
## Stratify the seurat object based on expresion quantile of a given gene, and Plot a series of custom plots for stratified data.
##=================================================================================================================================
stratifyByExpression <- function(object = all.seur.combined, assay = "RNA", strat.by = "GAPDH", gene.groups = list("Foamy genes" = c("ABCA1", "ABCG1", "OLR1"), "Resident genes" = c("LYVE1", "MRC1", "FOLR2"), "LAM genes" = c("TREM2", "CD9", "GPNMB"), "Inflammatory genes" = c("IL1B", "TNF", "NLRP3", "CASP1")), file.name = "myplot", return.object = T, do.plot = T, verbose = T, onlyUMAP = F, lower.quantile = 0.25, upper.quantile = 0.75){
  # Sanity check gene.groups
  if(!onlyUMAP & type(gene.groups) != "list"){
    cat("ERROR: gene.groups argument has to be a named list of vectors with gene names!\n")
    if(return.object){
      return(object)
    }else{
      return(invisible(NULL))
    }
  }
  
  # Get marker expression levels
  if(!strat.by %in% row.names(object)){
    cat("ERROR:", strat.by, "not present in this Seurat object's RNA assay!\n")
    if(return.object){
      return(object)
    }else{
      return(invisible(NULL))
    }
  }
  
  theMarker <- GetAssayData(object = object, assay = "RNA")[strat.by,]
  
  # Check some stats
  if(verbose){
    cat("Cells not expressing ",          strat.by, " ", sum(theMarker == 0), "\n", sep = "")
    cat("Cells expressing ",              strat.by, " ", sum(theMarker != 0), "\n", sep = "")
    cat("Perentage of cells expressing ", strat.by, " ", round((sum(theMarker != 0) / length(theMarker)) * 100), "%\n", sep = "")
  }
  
  # First subset the cells with expression of theMarker
  theMarker.expressed <- theMarker[theMarker != 0]
  
  # Then stratify the expressing cells
  if(verbose){
    cat("Quantiles of", strat.by, "expressing cells:\n", names(summary(theMarker.expressed)),"\n", summary(theMarker.expressed), "\n")
  }
  
  theMarker.low  <- names(theMarker.expressed)[theMarker.expressed <= quantile(theMarker.expressed, lower.quantile)]
  theMarker.med  <- names(theMarker.expressed)[theMarker.expressed >  quantile(theMarker.expressed, lower.quantile) & theMarker.expressed <= quantile(theMarker.expressed, upper.quantile)]
  theMarker.high <- names(theMarker.expressed)[theMarker.expressed >  quantile(theMarker.expressed, upper.quantile)]

  # Add as metadata to the seurat object
  marker.col <- paste(strat.by, "expr", sep = "_")
  object     <- AddMetaData(object, metadata = rep("Zero", ncol(object)), col.name = marker.col)

  object@meta.data[row.names(object@meta.data) %in% theMarker.low,  marker.col] <- "Low"
  object@meta.data[row.names(object@meta.data) %in% theMarker.med,  marker.col] <- "Medium"
  object@meta.data[row.names(object@meta.data) %in% theMarker.high, marker.col] <- "High"

    object@meta.data[,marker.col] <- factor( object@meta.data[,marker.col], levels = c("Zero", "Low", "Medium", "High"))
  
  # Plot the stratification
  if(do.plot){
    VlnPlot(object   = object, 
            assay    = assay, 
            features = strat.by, 
            group.by = marker.col, 
            cols     = c("grey", "bisque", "coral", "firebrick")) + xlab(paste("Stratified by", strat.by ,"expression", sep = " "))
    ggsave(filename  = paste(file.name, " - violin.pdf"))
    
    customUMAP(object    = object, 
               group.by  = marker.col,
	       order     = c("High", "Medium", "Low", "Zero"),
               pt.size   = 4, 
               cols      = c("grey", "bisque", "coral", "firebrick"), 
               title     = paste(strat.by, "expression", sep = " "),
               file.name = paste(file.name, " - UMAP.pdf"))
    
    # And some genes, if requested
    if(!onlyUMAP){
      for(theGroup in names(gene.groups)){
        stratPlots(object   = object, 
                   group.by = marker.col, 
                   features = gene.groups[[theGroup]], 
                   cols     = c("grey", "bisque", "coral", "firebrick"),
                   x.label  = paste("Stratified by", strat.by ,"expression", sep = " "), 
                   name     = paste(file.name, " - ", theGroup, sep = ""),
                   ncol     = 1)
      }
    }
  }
    
  if(return.object){
    return(object)
  }
}

#=========================================================================================================================
## Plot a customised UMAP (dimplot), keeping standard values as most often used in this project.
##========================================================================================================================
# Wrapper to extract plot limits so we can dynamically alter the position of the 'axes arrows' on the umap plot
get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

customUMAP <- function(object = object, group.by = NULL, order = NULL, pt.size = 4, label = F,label.size = NULL, cols = NULL, title = NULL, font.size = 14, reduction = "umap", shuffle = T, legend.pos = "top", seed = 666, file.name = "myplot.pdf", plot.height = 10, plot.width = 10, sizes.highlight = 1, cells.highlight = NULL, cells = NULL){
  p1 <- DimPlot(object = object, group.by = group.by, order = order, pt.size = pt.size, label = label, label.size = label.size, cols = cols, reduction = reduction, shuffle = shuffle, seed = seed, sizes.highlight = sizes.highlight, cells.highlight = cells.highlight, cells = cells) +
        ggtitle(title) +
        xlab("UMAP 2") +
        ylab("UMAP 1") + 
        theme_pubr(base_size = font.size, legend = legend.pos) +
        theme(axis.line        = element_blank(), 
              axis.text        = element_blank(), 
              axis.ticks       = element_blank(),
              axis.title       = element_text(hjust = 0.025),
              panel.background = element_blank(),
              title            = element_text(size = (font.size + 2), face = "bold"))

  # Add dynamic 'axes arrows'
    p1 <- p1 + coord_cartesian(xlim   = c((floor(get_plot_limits(p1)$xmin) - 0.2), ceiling(get_plot_limits(p1)$xmax)),
                               ylim   = c((floor(get_plot_limits(p1)$ymin) - 0.2), ceiling(get_plot_limits(p1)$ymax)),
                               expand = F, 
                               clip   = "off") +
          annotate(geom = "segment", 
                   x    = floor(get_plot_limits(p1)$xmin), 
                   xend = (floor(get_plot_limits(p1)$xmin) + 1.5), 
                   y    = floor(get_plot_limits(p1)$ymin), 
                   yend = floor(get_plot_limits(p1)$ymin), 
                   lineend = "round", lwd = 2, arrow = grid::arrow(length = unit(15, "pt"), angle = 20)) +
          annotate(geom = "segment", 
                   x    = floor(get_plot_limits(p1)$xmin), 
                   xend = floor(get_plot_limits(p1)$xmin),
                   y    = floor(get_plot_limits(p1)$ymin), 
                   yend = (floor(get_plot_limits(p1)$ymin) + 1.5), 
                   lineend = "round", lwd = 2, arrow = grid::arrow(length = unit(15, "pt"), angle = 20))
ggsave(filename = file.name, height = plot.height, width = plot.width, plot = p1)
}

#=========================================================================================================================
## Plot a customised violinPlot, keeping standard values as most often used in this project.
##========================================================================================================================
customVln <- function(object = samples.seurat, group.by = NULL, idents = NULL, features = "GAPDH", assay = "RNA", draw.names = T, name = "plot_name.pdf", splitPlot = F, ncol = NULL, stack = F, pt.size = 0, width = 15, height = 15, split.by = NULL, cols = NULL){
  if(stack == T){
    # If we stack we also wanna flip
    VlnPlot(pt.size    = pt.size, 
            group.by   = group.by,
            idents     = idents,
            object     = object, 
            features   = features,
            assay      = assay,
            split.by   = split.by, 
            split.plot = splitPlot,
            stack      = T,
            flip       = T) &
      theme_pubr(base_size = 14, x.text.angle = 45) &
      theme(panel.background = element_blank(),
            plot.margin      = unit(c(1,1,1,5), units = "cm"),
             title           = element_text(size = 16, face = "bold")
      ) 
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
  }else{
    if(draw.names == T){
    VlnPlot(pt.size    = pt.size, 
            group.by   = group.by, 
            idents     = idents,
            object     = object, 
            features   = features,
            assay      = assay,
            ncol       = ncol,
            split.by   = split.by, 
            split.plot = splitPlot,
            cols       = cols) & 
      theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") &
      theme(panel.background = element_blank(),
            plot.margin      = unit(c(1,1,1,5), units = "cm"),
            title            = element_text(size = 16, face = "bold")
      )
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
    }else{
      VlnPlot(pt.size    = pt.size, 
              group.by   = group.by, 
              idents     = idents,
              object     = object, 
              features   = features,
              assay      = assay,
              ncol       = ncol,
              split.by   = split.by, 
              split.plot = splitPlot,
              cols       = cols) & 
        theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") &
        theme(panel.background = element_blank(),
              plot.margin      = unit(c(1,1,1,5), units = "cm"),
              title            = element_text(size = 16, face = "bold"),
              axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()
        )
      ggsave(filename  = name, 
             width     = width, 
             height    = height,
             limitsize = F)
    }
  }
}


#=========================================================================================================================
## Plot a customised DotPlot, keeping standard values as most often used in this project.
##========================================================================================================================
customDot <- function(object = samples.seurat, group.by = NULL, idents = NULL, features = "GAPDH", assay = "RNA", cluster.idents = T, name = "plot_name.pdf", split.by = NULL, dot.scale = 7, width = "auto", height = "auto"){
  if(width == "auto"){
    # Base width 5. Add 0.2 for every extra feature, and add the maximum string width of the cluster names.
    pop.width <- max(strwidth(levels(Idents(object)), units = "inches")) * 2.54 # Convert to cm
    if(length(features) <= 5){
      width <- 5 + pop.width
      } else{
        width <- 5 + ((length(features) - 5) * 0.5) + pop.width
      }
  }

  if(height == "auto"){
    # Determine number of categories on the y axis:
    # Number of idents...
    if(is.null(group.by)){
      num.cats <- length(unique(Idents(object)))
    # Or number of whatever we're grouping by...
    }else{
      num.cats <- length(unique(object@meta.data[,group.by]))
    }
    # multiplied by the number of categories we are splitting by
    if(!is.null(split.by)){
      num.cats <- num.cats * length(unique(object@meta.data[,split.by]))
    }
    
    # Base height 5. Add 0.25 for every extra category
    if(num.cats <= 5){
      height <- 5}
    else{
        height <- 5 + ((num.cats - 5) * 0.25)
      }
    }

    DotPlot(object       = object,
          group.by       = group.by, 
          idents         = idents,
          cluster.idents = cluster.idents,
          dot.scale      = dot.scale, 
          split.by       = split.by, 
          cols           = "Spectral", 
          assay          = assay,
          features       = features) + 
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")) +
    scale_color_gsea()
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
}


#=========================================================================================================================
## Plot a customised FeaturePlot, keeping standard values as most often used in this project.
##========================================================================================================================
customFeature <- function(object = samples.seurat, cols = c("grey", "blue"), features = "GAPDH", name = "plot_name.pdf", reduction = "umap", pt.size = 1, order = T, width = 15, height = 15, ncol = NULL){
  FeaturePlot(object    = object,
              features  = features,
              cols      = cols,
              reduction = reduction, 
              pt.size   = pt.size, 
              order     = order, 
              ncol      = ncol) & 
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") &
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold"), aspect.ratio = 1)
  ggsave(filename  = name, 
         width     = width, 
         height    = height,
         limitsize = F)
}

#=========================================================================================================================
## Plot a bunch of custom plots at once, using most commonly ussed variations.
## Don't add .pdf to the name here as we are making compound names based on the settings used!
## Implemented: feature, violins and DotPlots
##========================================================================================================================
bunchOfCustomPlots <- function(object = samples.seurat, idents = NULL, group.by = NULL, features = "GAPDH", assay = "RNA", reduction = "umap", Vln.draw.names = T, name = "plot_name", feature.pt.size = 3, Vln.pt.size = 0, dot.scale = 7, Vln.width = 15, Vln.height = 15, Vln.stack = FALSE, Vln.color = NULL, Dot.width = "auto", Dot.height = "auto", ncol = NULL){
  # Feature plot
  customFeature(object = object, features = features, name = paste(name, " - feature plot.pdf", sep = ""), ncol = ncol, pt.size = feature.pt.size, reduction = reduction)

  # Violin plots
  customVln(object = object, idents = idents, group.by = group.by, features = features, assay = assay, draw.names = Vln.draw.names, name = paste(name, " - violin plot.pdf", sep = ""), width = Vln.width, height = Vln.height, pt.size = Vln.pt.size, ncol = ncol, stack = Vln.stack, cols = Vln.color)
  
  # Dot plots
  customDot(object = object, idents = idents, group.by = group.by, features = features, assay = assay, name = paste(name, " - dot plot.pdf", sep = ""),    width = Dot.width, height = Dot.height, dot.scale = dot.scale)
}