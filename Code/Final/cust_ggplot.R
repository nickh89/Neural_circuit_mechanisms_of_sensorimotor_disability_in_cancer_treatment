cust_ggplot <- function (dts,
                  x.var,
                  y.var,
                  fill.var)  {
  
  dts <- dts[, f.var := as.factor(dts[, get(fill.var)])] 
  #NB this will modify your original table
  # setorderv(dts, order.var)
  
  ggplot(dts, 
         aes_string(x = x.var, y = y.var, fill = "f.var")) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values =  c("#808080", "#7F4B9A"), name = "Treatment")+
    theme_bw(base_size = 16)+
    theme_classic()

    # scale_y_continuous(labels = function(val.var) { 
    #   format(val.var, big.mark = " ", scientific = FALSE)
    # })+
    # labs(x = "x.label", y = "y.label", title = "str.title")
  # https://stackoverflow.com/questions/38917301/another-way-to-pass-variable-as-factor-to-fill-argument-in-ggplot2-function
}


