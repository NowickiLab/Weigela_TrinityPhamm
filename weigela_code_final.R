## Weigela Dataset analysis

library("poppr")
library("hierfstat")
library("data.table")
library("sjPlot")
library("tibble")
library('phytools')
library("ggplot2")
library("RColorBrewer")


##################### General Setup ####################

setwd("____") #insert working directory here

genbank_start <- 158066
num_of_primers <- 19 ## subtract 1 from the actual number
gb_letter <- "OM" ## letter the genbank numbers start with


# Cultivar names
cultivars <- c('Spilled Wine®', 'Date Night™ Maroon Swoon®', 'Date Night™ Maroon Swoon®', 
               'Electric Love®', 'Electric Love®', 'Tango', 'Crimson Kisses®',
               'Crimson Kisses®', 'Sonic Bloom® Pink', 'Sonic Bloom® Pink', 'Dark Horse',
               'Dark Horse', 'Minuet', 'Minuet', 'Pink Poppet', 'Pink Poppet', 'Rainbow Sensation™',
               'Rainbow Sensation™', 'Red Prince', 'Red Prince', 'Shining Sensation™',
               'Shining Sensation™', 'Suzanne', 'Suzanne', 'Minor Black', 'Minor Black',
               'Towers of Flowers® Apple Blossom', 'Towers of Flowers® Apple Blossom',
               'Towers of Flowers® Cherry', 'Towers of Flowers® Cherry', 'Tuxedo™', 'Tuxedo™',
               'Stunner™', 'Stunner™')


# Read in the data
genalex <- read.genalex(genalex = "weigela_genalex.csv", #the csv file has to be in genalex format
                        ploidy = 2, 
                        genclone = FALSE) #this one assumes 1 population
genalex
summary(genalex)

# Make sure no loci have more than 30% missing data or samples with > 50%
genalex <- missingno(genalex, type = "loci", cutoff = 0.30, quiet = FALSE,
                     freq = FALSE) 
genalex <- missingno(genalex, type = "genotype", cutoff = .5)
info_table(genalex, plot = TRUE) # gives overview of % missing data
setPop(genalex) <- ~Pop 

genalex
summary(genalex)


######################## Diversity Indices ###################

basic.stats(genalex) ## just to get the overal Ho

DI <- locus_table(genalex)
n_alleles <- as.data.frame(round(DI[,1], digits = 0))
colnames(n_alleles) <-  "N"

## adding % missing data to table
missing_table <- info_table(genalex)
missing <- missing_table[1,]*100
missing <- unname(missing)
missing <- round(missing, 2)
datatable <- add_column(n_alleles, missing)

## adding Genbank numbers to table
genbank <- seq(genbank_start, genbank_start+num_of_primers, by=1)
genbank <- as.character(genbank)
gb_letter_vector <- character(length(genbank))
gb_letter_vector[] <- gb_letter
final_genbank <- paste(gb_letter_vector, genbank, sep = "")
final_genbank <- c(final_genbank, 0)
datatable <- add_column(datatable, final_genbank, .before = "N")
colnames(datatable)[1] <- "GenBank #"

## adding allele size ranges to table
size_ranges <- c('135-150',
                 '148-181',
                 '150-202',
                 '159-192',
                 '156-175',
                 '174-201',
                 '237-249',
                 '280-297',
                 '292-325',
                 '201-217',
                 '304-329',
                 '340-362',
                 '308-340',
                 '348-361',
                 '383-432',
                 '388-406',
                 '138-197',
                 '171-221',
                 '266-287',
                 '281-320',
                 '--')

datatable <- add_column(datatable, size_ranges, .before = "N")

tab_df(datatable,
       show.rownames = TRUE,
       title = "Table 1. Characteristics of the 20 gSSRs developed for Weigela",
       file = "diversity_indices.doc")


################## Linkage disequlibrium ##################

## Index of Association
set.seed(1044)
ia(genalex, hist = TRUE, sample = 999)


## pairwise linkage disequilibrium
setPop(genalex) <- ~Pop 
genalex_clone<-clonecorrect(genalex, strata = ~Pop, combine = FALSE, keep = 1) # have to clone correct before doing this

genalex_pair_clone <- pair.ia(genalex_clone)
plotrange <- range(genalex_pair_clone, na.rm = TRUE)
png("pairwise_IA_cc.png", width = 9, height = 9, units = 'in', res = 600)
plot(genalex_pair_clone, limits = plotrange, label=TRUE, low = "#ffffcc", high = "#41b6c4") +
  theme(legend.text = element_text(size = 10), axis.text = element_text(size = 15))
dev.off()



##################### PCoA #########################

## bruvo distance

repeat_length <- c(4, 2, 2, 3, 4, 4, 4, 3, 2, 3, 2, 3, 4, 2, 4, 2, 2, 2, 2, 2)

dist_matrix <- bruvo.dist(genalex, replen = repeat_length)
bruvo_pcoa <- ape::pcoa(dist_matrix)
bruvo_pcoa$values
var_bruvo_pcoa <- round(bruvo_pcoa$values$Relative_eig[1:2],4)

bruvo_df <- data.frame(PC1 = bruvo_pcoa$vectors[, 1], PC2 = bruvo_pcoa$vectors[, 2])
bruvo_df$Cultivar <- cultivars
bruvo_df$Cultivar <- factor(bruvo_df$Cultivar, levels = c("Electric Love®", "Minuet", "Stunner™", 
                                                          "Tango", "Dark Horse", "Towers of Flowers® Cherry",
                                                          "Tuxedo™", "Crimson Kisses®", "Date Night™ Maroon Swoon®",
                                                          "Red Prince", "Minor Black", "Pink Poppet",
                                                          "Rainbow Sensation™", "Shining Sensation™", "Sonic Bloom® Pink",
                                                          "Spilled Wine®", "Suzanne", "Towers of Flowers® Apple Blossom"))



### for the alphabetical legend

uniqueInitials <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
                    "N", "O", "P", "Q", "R")
initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))

zp2 <- ggplot(bruvo_df,aes(x = PC1, y = PC2, shape = Cultivar))
zp2 <- zp2 + geom_point(size = 3.75) + labs(x= paste("PC1 (",var_bruvo_pcoa[1]*100,"%)"), 
                                            y = paste("PC2 (",var_bruvo_pcoa[2]*100,"%)"))
zp2 <- zp2 + theme_bw()  # \/ Manually select letter-shapes
zp2 <- zp2 + scale_shape_manual(values = initialShapes)


png("PCoA_black_letters.png", 
    width = 9, height = 5.5, units = 'in', res = 600)
print(zp2)
dev.off()

## colors to use
ten_col <- c("#08519c", "#3182bd", "#6baed6", "#bdd7e7", #blue
             "#a63603", "#e6550d", "#fd8d3c", #orange
             "#54278f", "#756bb1", "#9e9ac8", rep('black', 8)) #purple



## base PCoA plot
bp = ggplot(bruvo_df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Cultivar), size = 4, alpha = 0.85) + 
  labs(x= paste("PC1 (",var_bruvo_pcoa[1]*100,"%)"), 
       y = paste("PC2 (",var_bruvo_pcoa[2]*100,"%)"))



png("PCoA_bruvo.png", width = 9, height = 5.5, units = 'in', res = 600)
bp + scale_colour_manual(values = ten_col) 
dev.off()




###################### Neighbor-Joining ###################

## cailliez correction results in higher correlation coefficients
D <- cailliez(dist_matrix)
trenj_cal <- nj(D)

trebionj_cal<- bionj(D)

trefastmebal_cal <- fastme.bal(D)

trefastmeols_cal <- fastme.ols(D)

trehclust_cal <- hclust(D)
trehclust_cal <- as.phylo(hclust(D,method="average")) # is the same as UPGMA

#CHECK WHICH ALGORITHM IS BEST
##hclust
a <- as.vector(D)
b <- as.vector(as.dist(cophenetic(trehclust_cal)))
plot(a, b, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is UPMJA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(a~b), col="red")
cor(a,b)^2 #.858

##fastmeols
a <- as.vector(D)
b <- as.vector(as.dist(cophenetic(trefastmeols_cal)))
plot(a, b, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is fastme.ols appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(a~b), col="red")
cor(a,b)^2 #.145

##fastmebal
a <- as.vector(D)
b <- as.vector(as.dist(cophenetic(trefastmebal_cal)))
plot(a, b, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is fastme.bal appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(a~b), col="red")
cor(a,b)^2 #.319

##trenj
a <- as.vector(D)
b <- as.vector(as.dist(cophenetic(trenj_cal)))
plot(a, b, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(a~b), col="red")
cor(a,b)^2 #.900

######### Best #########
##bionj
a <- as.vector(D)
b <- as.vector(as.dist(cophenetic(trebionj_cal)))
plot(a, b, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is bionj appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(a~b), col="red")
cor(a,b)^2 #.901


##### Bootstrapped bionj Plot #####

set.seed(1044)
myBoots <- bruvo.boot(genalex, replen = repeat_length, sample = 1000, cutoff = 50, quiet = FALSE, tree = bionj)


png("bionj_bruvo.png", width = 6, height = 6, units = 'in', res = 600)
plot(myBoots, type="unr", tip.col="black", font=2, show.tip.label = FALSE)
add.scale.bar()
dev.off()

# to figure out labels and branch lengths

png("bionj_bruvo_tb_labels.png", width = 6, height = 6, units = 'in', res = 600)
l.col <- funky(54)

plot(myBoots, type="unr", tip.col="black", font=2, show.tip.label = FALSE)
tiplabels(cultivars, bg=l.col, cex=.5, fg="transparent")
annot <- round(myBoots$edge.length,2)
edgelabels(annot[annot>0.1], which(annot>0.1), frame="n")
add.scale.bar()
dev.off()

# for bootstrap labels
png("bionj_bruvo_bs_labels.png", width = 6, height = 6, units = 'in', res = 600)
plot(myBoots, type="unr", tip.col="black", font=2, show.tip.label = FALSE)
nodelabels(text = myBoots$node.label, frame = 'none', font = 3, adj = c(1.2,-.49))
add.scale.bar()
dev.off()


############# Private alleles ################

pop(genalex) <- cultivars

# number of cultivars
table(cultivars)

# includes all cultivars
pop_private_all <- private_alleles(genalex, form = locus ~., count.alleles = FALSE, level = "population")
# to get a list of the ones that have private alleles
pop_private <- private_alleles(genalex, form = locus ~., count.alleles = FALSE, level = "population",
                               drop = TRUE)

# find the cultivars that don't have private alleles
setdiff(rownames(pop_private_all), rownames(pop_private))

# Totals of private alleles for each cultivar (dosage independent)
no_dose <- rowSums(pop_private_all)

# Totals of prviate alleles for each cultivar (dosage dependent)
pop_private_count_true <- private_alleles(genalex, form = locus ~., level = "population")
dose <- rowSums(pop_private_count_true)

private <- data.frame(no_dose)
colnames(private) <- "Private Alleles"
tab_df(private,
       show.rownames = TRUE,
       title = "Table _. Number of Private Alleles Found in Each Cultivar",
       file = "private_alleles.doc")

