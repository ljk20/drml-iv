
library(tidyverse)
library(gridExtra)

rm(list=ls())

element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}


load("./output/scenario1_results_ivstr.RData" )

sims = scenarios
sims = unnest( sims, res )


out = sims %>% group_by( iv.str, method ) %>%
  summarise( rmse = mean( est, trim = .025 ),
             bias = sqrt(mean( clate^2, trim = .05 )), 
             .groups = "drop" )

sqrt( mean((est - LATE)^2))
##### Facet Set up
df = out %>% mutate( method = recode(method, 
                                    DR = "DRML - Nonparametric",
                                    FF = "Frauen & Feuerriegel"))


dflL = df %>% 
  pivot_longer( cols = c( rmse, bias ), names_to = "metric", values_to = "stat" )
  
dflL$metric = factor(dflL$metric,
                   levels = c( "rmse", "bias" ),
                   labels = c( "RMSE", "Bias" )) 
                   

dflL$iv.str = factor(dflL$iv.str,
                   levels = c(1.1, 2, 2.9, 5.5),
                   labels = c( "15", "30", "50", "65"))  

               
data.bias <- dflL %>% filter(metric == "Bias") 

plot.1 <- ggplot(data.bias, aes(x=iv.str, y=stat, group=method)) + 
                 geom_line(aes(color=method)) + 
                 geom_point(aes(shape=method), size=1.5) + 
                 xlab("Share of Compliers") +  
                 ylab( "" ) +
                 theme_minimal() +
                 theme(legend.position="none") +
                 labs(title = "RMSE for CLATE") +
  						theme(
    					plot.title = element_textbox(
      					hjust = 0.5, margin = margin(t = 5, b = 5)
    					)
  					)          
            
plot.1

data.rmse <- dflL %>% filter(metric == "RMSE") 
                    
plot.2 <- ggplot(data.rmse, aes(x=iv.str, y=stat, group=method)) + 
                 geom_line(aes(color=method)) + 
                 geom_point(aes(shape=method), size=1.5) + 
                 xlab("Share of Compliers") +  
                 ylab( "" ) +
                 theme_minimal() +
                 theme(legend.title=element_blank()) +
                 labs(title = "RMSE for ITE Distribution") +
  						theme(
    					plot.title = element_textbox(
      					hjust = 0.5, margin = margin(t = 5, b = 5)
    					)
  					)          
            
plot.2

grid.arrange(plot.1, plot.2, nrow=1, widths = c(.6,.8))




