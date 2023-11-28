#Calculates p-values from two sets of columns in a dataframe, returns a vector
ttestFun <- function(df, col1, col2){
  result = NULL
  for(i in 1:(nrow(df))){
    v1 = (as.numeric(df[i, col1]))
    v2 = (as.numeric(df[i, col2]))
    t_test = t.test(v1, v2, var.equal = TRUE)
    pval = t_test$p.value
    result[i] = pval
  }
  result
}


#Calculates Pearson correlation between column from dataframe and vector, returns dataframe with new "cor" column
CorFun <- function(x, xcols, y){
  for(i in 1:nrow(x)){
    v = as.numeric(x[i,xcols])
    c = cor(v, y, method = "pearson")
    x[i,"cor"] = c
  }
  x
}

#Retrieves legend from ggplot
get_legend<-function(myggplot){
  tmp = ggplot_gtable(ggplot_build(myggplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}