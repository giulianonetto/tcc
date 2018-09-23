library(ggplot2)
data("mtcars")
qplot(mpg, wt, data = mtcars, colour = factor(cyl))
ggplot(mtcars, aes(x=mpg, y=wt))+
  geom_point(aes(col=factor(rownames(mtcars))))+
  geom_smooth(method="lm")+
  ggtitle("Hello world")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_void()


library(stringr)
dat <- read.table('/home/giulianonetto/windows/tcc/checklist/genesigs/papers/Xue, 2014/modules.csv', sep = '\t', header = T, stringsAsFactors = F)

dat_colnames <- colnames(dat)
dat_colnames_new = vector()
for (i in 1:ncol(dat)){
  old <- str_extract(dat_colnames[i], "[0-9]+")
  new = paste0("Module", old)
  dat_colnames_new = c(dat_colnames_new, new)
}
colnames(dat) <- dat_colnames_new

datatable(dat, rownames = FALSE, filter = 'top',
          caption = "Table 1. Modules defined by Xue, 2014.",
          options=list(pageLength = 100, scrollX='400px', scrollY='400px',
                       lengthMenu = c(30, 100, 1000) 
                       
          )
)

colcount <- data.frame()
for(i in 1:ncol(dat)){
  module = colnames(dat)[i]
  Counts =  length(unique(dat[[i]])) - 1
  print(c(module, Counts))
  colcount = rbind(colcount,data.frame(module, Counts))
}
# colcount <- data.frame('numbers'=as.numeric(colcount))
qplot(colcount$colcount, geom = "histogram", 
      bins = length(colcount$colcount)/2)

ggplot(data=colcount, aes(colcount$numbers))+
  geom_density(aes(fill = 1))
ggplot(colcount)+
  geom_col(aes(x=module, y=Counts))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=FALSE, colour=FALSE)+
  ylab("Number of genes")+
  ggtitle("Number of genes per module")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  scale_fill_manual("red")

ggplot(colcount, aes(Counts))+
  geom_dotplot(aes(fill=Counts))
               