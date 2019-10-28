## missing genome reconstruction

Let's make a reference graph, align reads to it, and use surjection and k-means clustering to attempt to find the reads that map to the 5th genome. Finally, we'll assemble these and see if it matches known information about the fifth genome.

First we'll build and index the graph.

We build the graph using `vg msga` from the input sequences. (Note that the assembly of these 4 genomes will always be circular, why is this?)

    Notes on paramters:
    vg msga:  
      -w, --band-width INT band/chunk width for long read alignment [default: 128]
      -D, --debug print debugging information about construction to stderr
    vg mod:
      -U, --until-normal N iterate normalization until convergence, or at most N times
      -c, --compact-ids should we sort and compact the id space? (default false)
```
vg msga -f REF.fasta -w 1024 -D | vg mod -U 10 - | vg mod -c -  >REF5.vg
```

We then index the graph with xg:

```
vg index -x REF5.xg REF5.vg
```

And we prune and index with GCSA2:
    vg mod:
      -g, --subgraph ID gets the subgraph rooted at node ID, multiple allowed
      -p, --prune-complex remove nodes that are reached by paths of --length which cross more than --edge-max edges
      -l, --length N for pruning complex regions and short subgraphs
```
vg index -g REF5.gcsa -k 16 -p <(vg mod -pl 16 -e 3 REF5.vg)
```

### Surjection based exploration

We then map the first 10k Illumina reads:
 multipath mapping included here as an example as well
```
vg map -d REF5 -f <(cat SRR961514_1.fastq | head -40000) >SRR961514_1.first10k.gam
vg mpmap -x REF5.xg -g REF5.gcsa -S -f<(cat SRR961514_1.fastq | head -40000) > SRR961514_1.first10k_mpmap.gam
```
View the first alignment:
    vg view -a SRR961514_1.first10k.gam | head -1 | vg view -JaG - >first_aln.gam
    vg find -x REF5.xg -G first_aln.gam | vg view -dA first_aln.gam - | dot -Tpdf -o first_aln.pdf

Our goal is to find the un-assembled genome within the Illumina reads. To do so we will use a technique based on obtaining the path-relative identity for each read against each of our 4 [5] reference genomes. We hypothesize that the reads from the 5th genome will tend to have a unique pattern of relationships between the known genomes, and should show up in clustering analysis.

First we surject the reads into each one of the references and build a table of the output:

NOTE: Need to sort and index the GAM file, just like manipulating a BAM file

```
vg gamsort SRR961514_1.first10k.gam > SRR961514_1.first10k.sorted.gam
vg index -l SRR961514_1.first10k.sorted.gam
```

See the different alignments of the first read for the two paths.
    vg surject -x REF5.xg -p HXB2 SRR961514_1.first10k.sorted.gam - > test.gam
    vg view -a test.gam | head -1 | vg view -JaG - >test_HXB2_aln.gam
    vg find -x REF5.xg -G test_HXB2_aln.gam | vg view -dA test_HXB2_aln.gam - | dot -Tpdf -o test_HXB2_aln.pdf

    vg surject -x REF5.xg -p YU2 SRR961514_1.first10k.sorted.gam - > test.gam
    vg view -a test.gam | head -1 | vg view -JaG - >test_YU2_aln.gam
    vg find -x REF5.xg -G test_YU2_aln.gam | vg view -dA test_YU2_aln.gam - | dot -Tpdf -o test_YU2_aln.pdf

###--
Make file with all of the path names:
    vg paths -L -v REF5.vg > REF5_path_names
    IDs: 'vg paths -L -v REF5.vg'
      896
      HXB2
      JRCSF
      NL43
      YU2

Surject the reads to non-overlapping path names listed within the above file
    vg surject -x REF5.xg -F REF5_path_names SRR961514_1.first10k.gam > SRR961514_1.first10k.surjected.gam

```bash
for ref in $(vg paths -L -v REF5.vg)
    do (vg view -a SRR961514_1.first10k.gam - | jq -cr '[.name, "'$ref'", .identity]' | sed s/null/0/g | jq -cr @tsv ) | gzip >first10k.surj.$ref.tsv.gz
done
```
###### this just gives identical results... can't seem to get the path-specific identity scores for each read

####### augment the graph with all variation from the GAM except that implied by soft clips, saving to SRR961514_1.first10k.aug.vg.  SRR961514_1.first10k.aug.gam contains the same reads as SRR961514_1.first10k.gam but mapped to SRR961514_1.first10k.aug.vg
    vg augment REF5.vg SRR961514_1.first10k.gam -C -A SRR961514_1.first10k.aug.gam > SRR961514_1.first10k.aug.vg


##### from original PANgenomics github: (non-functional)
```bash
for ref in $(vg paths -L -v REF5.vg)
    do ( vg surject -x REF5.xg -p $ref SRR961514_1.first10k.sorted.gam - | vg view -a - \
        | jq -cr '[.name, "'$ref'", .identity]' | sed s/null/0/g | jq -cr @tsv ) | gzip >first10k.surj.$ref.tsv.gz
done
```
Now we join up the results to build one table: (Or we would, had we managed to create different identities for each read relative to each reference)
```
( echo read.name id.896 id.HXB2 id.JRCSF id.NL43 id.YU2
join <(zcat first10k.surj.896.tsv.gz | cut -f 1,3 | sort -V) <(zcat first10k.surj.HXB2.tsv.gz | cut -f 1,3 | sort -V) \
| join  <(zcat first10k.surj.JRCSF.tsv.gz | cut -f 1,3| sort -V) - \
| join  <(zcat first10k.surj.NL43.tsv.gz | cut -f 1,3| sort -V) - \
| join  <(zcat first10k.surj.YU2.tsv.gz | cut -f 1,3| sort -V) - ) \
| tr ' ' '\t' | gzip >first10k.surj.join.tsv.gz
```

In R we can run k-means clustering PCA to evaluate if there is a subset that doesn't match any of the references.
We take a number of steps to ensure that the PCA makes sense. First, we remove cases where we don't have a relatively good match for all of the references, as the number that don't are relatively few we consider these outliers.
Then, we normalize each row so that the identity sum is 1.
Finally, we run the PCA and k-means clustering with a center count of 5 (as we are assuming there are 5 samples mixed together).

module load R
```R
require(tidyverse)
require(broom)
require(ggplot2)
first10k <- (read.delim("first10k.surj.join.tsv.gz") %>% mutate(id.sum=id.896+id.HXB2+id.JRCSF+id.NL43) %>% subset(id.sum > 0.95*4))
first10k.pca <- prcomp(first10k[,2:5]/first10k$id.sum)
first10k.pca.df <- as.data.frame(first10k.pca$x)
first10k.pca.df$read.name <- first10k$read.name
first10k.kmeans <- kmeans(first10k[,2:5]/first10k$id.sum, 5)
first10k$cluster <- as.factor(first10k.kmeans$cluster)
first10k.pca.df$cluster <- as.factor(first10k.kmeans$cluster)
ggplot(first10k.pca.df, aes(x=PC1, y=PC2, color=cluster)) + geom_point()
ggplot(first10k.pca, alpha=0.2, groups=first10k$cluster)
```

### Vectorize based analysis

The pacbio sequencing data from the experiment includes reads that are approximately 1300bp on average.
We can use vg vectorize to explore the relationships between the given reference genomes and the read set.
This generates a vector representation of the reads over the graph, where each alignment is represented by a vector in the order of the graph.
For each node that the alignment touches we record a 1, otherwise a 0.
(Different encodings are possible, relating to the goodness of alignment of each Mapping in the Alignment's Path to the graph.)

First we can run an alignment of part of the read set.

```
vg map -d REF5 -f <(cat SRR961669_1.fastq | head -4000) >SRR961669.first1k.gam
```

Then we extract the paths themselves as an alignment:

```
vg paths -v REF5.vg -X >REF5.paths.gam
```

Finally we can make a matrix representation of the read set and the graph paths.

```
vg vectorize -f -x REF5.xg <(cat REF5.paths.gam SRR961669.first1k.gam) | gzip >SRR961669.vectorized.1k.tsv.gz
```

In R we can explore the relationship by building a distance metric between each read and each reference path.
This is similar to surjection but avoids the need to actually run the local realignment.
Then we can explore the relationship between the sequences using PCA.

```R
pacbio <- read.delim('SRR961669.vectorized.1k.tsv.gz')
pacbio.m <- pacbio[2:ncol(pacbio)]
pacbio.dist <- data.frame(
    aln.name = pacbio$aln.name,
    node.count = rowSums(pacbio.m),
    id.896 = apply(pacbio.m, 1, function(x) dist(rbind(pacbio.m[1,], x))),
    id.HXB2 = apply(pacbio.m, 1, function(x) dist(rbind(pacbio.m[2,], x))),
    id.JRCSF = apply(pacbio.m, 1, function(x) dist(rbind(pacbio.m[3,], x))),
    id.NL43 = apply(pacbio.m, 1, function(x) dist(rbind(pacbio.m[4,], x))))
pacbio.dist.pca <- prcomp(-pacbio.dist[5:nrow(pacbio.dist),3:6])
ggbiplot(pacbio.dist.pca) + geom_point(aes(color=pacbio.dist[5:nrow(pacbio.dist),]$node.count+1)) + scale_color_continuous("node count") + theme_bw()
ggbiplot(pacbio.dist.pca) + geom_point(aes(color=pacbio.dist[5:nrow(pacbio.dist),]$node.count+1)) + scale_color_continuous("node count") + theme_bw()
```

Creates the image below:
![HIV PCA and Vectorization](https://raw.githubusercontent.com/dduchen/PANGenomics/master/images/HIV_Exercise2_PacBio_PCA.PNG)
