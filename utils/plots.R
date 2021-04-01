library(ggplot2)
library(plotly)

#' Theme for PCA plot
#' @Fsize font size
#' @Fcolor font color
#' @Fangle font angle
#' 
theme_PCA <- function(Fsize = 10, Fcolor = "black", Fangle = 0, Fhjust = 0.5) {
    theme_classic(base_size = Fsize) +
    theme(axis.text.x = element_text(size = Fsize, color = Fcolor,
          angle = Fangle, hjust = Fhjust),
          axis.text.y = element_text(size = Fsize, color = Fcolor),
          axis.title.x = element_text(size = Fsize, color = Fcolor),
          axis.title.y = element_text(size = Fsize, color = Fcolor),
          legend.text = element_text(color = Fcolor, size = Fsize - 1,
                                     margin = margin(l = 0, r = 20)),
          legend.title = element_text(color = Fcolor, size = Fsize,
                                      hjust = 0.5),
          legend.position = "bottom")
}

#' Generate PCA plot with ggplot2
#' 
plot_PCA_gg <- function(df, variance, pc1, pc2, category, category_name="", palette, cat_order, Lnr=2, Lnc=3) {
    ggplot(data = df, aes(x = !!sym(pc1), y = !!sym(pc2),
           color = !!sym(category))) +
        geom_point(size = 1.5) +
        scale_color_manual(breaks = cat_order, labels = cat_order,
                           values = palette[cat_order]) +
        guides(color = guide_legend(title = category_name,
                                    nrow = Lnr, ncol = Lnc)) +
        labs(x = paste0(pc1, " (", variance[pc1], "%)"),
             y = paste0(pc2, " (", variance[pc2], "%)"))
}


#' Generate PCA plot with plotly
#' @df
#' @variance
#' @pc1
#' @pc2
#' @category
#' 
plot_PCA_plotly <- function(df, variance, pc1, pc2, label_category, label_id, palette, Lposition="h") {
    p <- plot_ly(df,
                 x = df[, pc1], y = df[, pc2],
                 text = ~paste("</br> Sample:", df[, label_id],
                               "</br>", df[, label_category]),
                 color = df[, label_category],
                 colors = palette)
    p <- layout(p,
                legend = list(orientation = Lposition),
                xaxis = list(title = paste0(pc1, " (", variance[pc1], "%)"),
                             showgrid = F, zeroline = F),
                yaxis = list(title = paste0(pc2, " (", variance[pc2], "%)"),
                             showgrid = F, zeroline = F))
    p
}

# Theme for ADMIXTURE plot
theme_ADMIX <- function(Fsize=10, Fcolor="black", Fangle=0, Lposition="none", Fhjust=1) {
    theme_minimal() +
        theme(strip.text.x = element_text(angle = Fangle, hjust = Fhjust,
                                         size = Fsize - 1, color = "black"),
            panel.spacing.x = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = Fsize, color = "black"),
            axis.text.y = element_text(size = Fsize, color = "black",
                                       angle = 90),
            panel.grid = element_blank(),
            legend.position = Lposition)
}

# Plot ADMIXTURE
plot_ADMIX_gg <- function(df, formula, label_id="sample", label_Q = "Q", label_K="K", palette, Fangle=0, Lposition="none") {
    uniqK <- df %>% select(label_K) %>% n_distinct()
    ggplot(data = df, aes(x = !!sym(label_id), y = !!sym(label_Q),
                          fill = !!sym(label_K))) +
        geom_col(width = 1) +
        facet_grid(as.formula(formula),
                   switch = "x", scales = "free", space = "free") +
        scale_y_continuous(expand = c(0, 0.05), breaks = c(0, 0.5, 1.0)) +
        scale_x_discrete(expand = expansion(add = 1)) +
        scale_fill_manual(values = pals[1:uniqK],
                          guide = guide_legend(nrow = 1)) +
        labs(x = "Individuals", y = "Ancestry", fill = "K") +
        theme_ADMIX(Lposition = Lposition)



}