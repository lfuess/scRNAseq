##this code will create the base version of the graph shown in figure 3a of the manuscript##
##code was generated using data from the output in the loupe cell broswer##

library(ggplot2)
##create a matrix displaying the number of cells in each of our 8 groups for each of our three populations
  ##fill in the population column##
  pop <- c(rep("SAY" , 9), rep("GOS" , 9) , rep("ROB" , 9))
  ##then the type column##
  type <- rep(c("HCs" , "Neutrophils" , "APCs", "B-cells", "Erythrocytes", "Platelets", "Fibroblasts", "NK Cells", "Unclassified") , 3)
  ##finally the count column##
  count <- c(1029, 3829, 1585, 3807, 2054, 43, 29, 160, 267,
             3555, 10673, 1514, 5156, 2422, 83, 35, 318, 718,
             2836, 15075, 2219, 4116, 953, 211, 47, 214, 594)
  ##merge it all together##
  data <- data.frame(pop,type,count)
  ##adjust the levels for later graphing
  data$pop = factor(data$pop, levels = c("SAY", "GOS", "ROB"))
  data$type = factor(data$type, levels = c("HCs" , "Neutrophils" , "APCs", "B-cells", "Erythrocytes", "Platelets", "Fibroblasts", "NK Cells", "Unclassified"))

##then we just make a stacked percentage graph##
ggplot(data, aes(fill=type, y=count, x=pop)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#440154","#46337E","#365C8D","#277F8E","#1FA187","#4AC16D","#9FDA3A",
                               "#FDE725", "grey")) +
  ylab("Proportional Abundance") + xlab("Population") + labs(fill = "Cell Type") +
  scale_y_continuous(expand = expansion(mult = c(0,0))) +
  theme(panel.background = element_blank(), axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none")
