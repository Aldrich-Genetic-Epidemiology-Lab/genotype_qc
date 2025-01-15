#Written by Peter Straub


#--workingdir  Path to the output directory
#--pcaresultspath  Path to the the FlashPCA2 PCA results, run with three 1000genomes reference populations spiked in, that are to be parse
#--persongroupfilepath  Path to a two-column file with group assignments for samples
#--genomes1ksamplepoppath  Path to the 1000Genomes sample population file
#--genomes1krefgroups  a comma-delimited list of the three 1000genomes population used in the pcaresultspath flashPCA2 set   
#--plottinggroups  A "/"-deimited list of comma-delimited lists of groups selected for plotting

#--groupassignmentfilepath  Path to the groups assignment file generate from a previous run of this script where you specified the PC-axis boundaries for defining onr or more selected populations
#--rescaledpcaresultspath  Path to the rescaled PCA results from a previous run of this script
#--assignedgroups  A comma-delimited list of rows in the groupassignmentfilepath file that you want processed, based on the values in the assigned.eth column


library("getopt")



spec=matrix(c("help", "h", 0, "logical",
              "workingdir", "d", 1, "character",
              "pcaresultspath", "f", 1, "character",
              "persongroupfilepath", "t", 1, "character",
              "genomes1ksamplepoppath", "w", 1, "character",
              "genomes1krefgroups", "r", 1, "character",
              "plottinggroups", "n", 1, "character",
              "groupassignmentfilepath", "a", 1, "character",
              "rescaledpcaresultspath", "p", 1, "character",
              "assignedgroups", "g", 1, "character"
), byrow=TRUE, ncol=4)


opt=getopt(spec)

#if help, just exit
if (!is.null(opt$help))
{
  cat(getopt(spec, usage=TRUE))
  cat("Required parameters:
--workingdir  Path to the output directory
--pcaresultspath  Path to the the FlashPCA2 PCA results, run with three 1000genomes reference populations spiked in, that are to be parse
--persongroupfilepath  Path to a two-column file with group assignments for non-reference samples.  The first column should be the sample ID that matches the 'IID' second column found in the flashPCA2 file specified in pcaresultspath  
--genomes1ksamplepoppath  Path to the 1000Genomes sample population file
--genomes1krefgroups  a comma-delimited list of the three 1000genomes population used in the pcaresultspath flashPCA2 set.  Reference populations should follow the 'Super pop'-'pop' format (e.g. if CEU, YRI, and CHB are your reference populations, use 'EUR-CEU,AFR-YRI,EAS-CHB' for this parameter )

Optional parameters:      
--groupassignmentfilepath  Path to the groups assignment file generate from a previous run of this script where you specified the PC-axis boundaries for defining onr or more selected populations
--rescaledpcaresultspath  Path to the rescaled PCA results from a previous run of this script
--plottinggroups  A '/'-deimited list of comma-delimited lists of groups selected for plotting.  All found groups will be plotted in a single set if this parameter is not provided.  Use 'selected' to display the group of samples selected by your filter thresholds.  e.g. if you have already rescaled you PCs and are filtering for a subset you might choose 'W,B,A/W/B/A/W,B,A,selected/selected' for a sample with 'W', 'B', and 'A' as the sample groups
--assignedgroups  A comma-delimited list of rows in the groupassignmentfilepath file that you want processed, based on the values in the assigned.eth column

General instructions:
This script is set up to be run on sets of FlashPCA2 results files for sets that have already had three 1000genomes reference populations spiked in.  The first time the script is run on the set it will generate a rescaled PC set and plots color coded based on the persongroupfilepath data, along with a *_PC1_vs_PC2_assigned_group_boundaries.txt file.  This file can be used to define filter boundaries relative to the triangle formed by the three centers of mass for the three reference 1000genomes populations.  The three sides of triangle define three new axes for rescaling and filtering.  The *__PC1_vs_PC2_assigned_group_boundaries.txt can be used to assign refilter and rescale along these three reference axes.  When specifying a plot of the samples selected from your filter, use 'selected' in your plottinggroups lists.  

For example:  Let's assumed you have a flashPCA2 result set with the CEU, YRI, and CHB 1000genomes reference populations spiked in located at /my/data/targetplusref.pcs and a sample groups file located at /my/data/targetgroups.txt with groups 'W', 'B', 'A'.  You also have your Phase3 1000genomes sample population file located at /shared/1000genomes/sources/1000GP_Phase3/1000GP_Phase3.sample.  You would first want to generate a rescaled version of this set with the following command:
Rscript rescale_and_filter_pca.R --workingdir=/my/workingdir --pcaresultspath=/my/data/targetplusref.pcs --persongroupfilepath=/my/data/targetgroups.txt --genomes1ksamplepoppath=/shared/1000genomes/sources/1000GP_Phase3/1000GP_Phase3.sample --genomes1krefgroups='EUR-CEU,AFR-YRI,EAS-CHB' --plottinggroups='W,B,A/W/B/A'

After running the script, a targetplusref.scaled.pcs.txt file will be generated in /my/workingdir, along with a targetplusref_PC1_vs_PC2_assigned_group_boundaries.txt file and four plots (a rescaled plot of W,B,A together and plots of W, B and A separately as specified with our 'W,B,A/W/B/A' plottinggroups parameter)   

If you, for example, want to define an African ancestry subset of people relative to the AFR-YRI reference (the AFR-YRI point in the boundary triangle), you can modify the targetplusref_PC1_vs_PC2_assigned_group_boundaries.txt file to define those boundaries and rerun the script using the previously-generated rescaled PCs to generate new plots and plink 'keep' files of the selected samples that met your filter.

Let's say you want to only grab the samples that are within '20%' of the AFR-YRI cluster along the EUR-CEU:AFR-YRI axis and within '40%' along the AFR-YRI:EAS-CHB axis, you would modify the targetplusref_PC1_vs_PC2_assigned_group_boundaries.txt file as follows:

For the first group row in the file ('default'), set the value for the 'EUR-CEU_vs_AFR-YRI.AFR-YRI' column to 20 and 'AFR-YRI_vs_EAS-CHB.AFR-YRI' to 10.  Set the rest of the axis thresholds to NA.  Give this group a name in the assigned.eth column (e.g. 'AFR.subset1').  The assigned.ref.pop.for.refilter column is optionally available for assigning a new set of three 1000Genomes reference populations to your filtered subset.  A new plink keep file with your filtered subset and the three selected reference populations will be generated.  For example, if you wanted to regenerate a plink set of your AFR-filtered subset with the reference populations AFR-YRI, EAS-CHB, and SAS-GIH, you could set this column to 'AFR-YRI,EAS-CHB,SAS-GIH'. 
Example:
assignment.group.name   first.PC        second.PC       EUR-CEU_vs_AFR-YRI.EUR-CEU      EUR-CEU_vs_AFR-YRI.AFR-YRI      EUR-CEU_vs_EAS-CHB.EUR-CEU      EUR-CEU_vs_EAS-CHB.EAS-CHB      AFR-YRI_vs_EAS-CHB.AFR-YRI      AFR-YRI_vs_EAS-CHB.EAS-CHB      assigned.eth    assigned.ref.pop.for.refilter
default        PC1     PC2     NA      20      NA      40      NA      NA      AFR.subset1     AFR-YRI,EAS-CHB,SAS-GIH



You can rename your modified *_PC1_vs_PC2_assigned_group_boundaries.txt file for the refiltering step.  For example, targetplusref_PC1_vs_PC2_assigned_group_boundaries_AFRsubset1.txt

The following command would generate the refiltered subset focused around the AFR-YRI reference 'point' in the triangle:

Rscript rescale_and_filter_pca.R --workingdir=/my/workingdir --pcaresultspath=/my/data/targetplusref.pcs --persongroupfilepath=/my/data/targetgroups.txt --genomes1ksamplepoppath=/shared/1000genomes/sources/1000GP_Phase3/1000GP_Phase3.sample --genomes1krefgroups='EUR-CEU,AFR-YRI,EAS-CHB' --plottinggroups='W,B,A,selected/W/B/A/selected' --assignedgroups=AFR.subset1 --groupassignmentfilepath=/my/workingdir/targetplusref_PC1_vs_PC2_assigned_group_boundaries_AFRsubset1.txt --rescaledpcaresultspath=/my/workingdir/targetplusref.scaled.pcs.txt   

In addition to the newly generate plots, there will also be two plink keep files generated:  /my/data/targetplusref_AFR.subset1.refilter.no.refs.keep.txt will contain just the non-refence samples that passed your filter.  The /my/data/targetplusref_AFR.subset1.refilter.keep.txt keep will include the reference samples for three selected reference population in the assigned.ref.pop.for.refilter column of the groupassignmentfilepath file (AFR-YRI, EAS-CHB, and SAS-GIH in the example)   


            ")
  q(status=1)
}

print(opt)


spec=matrix(c("help", "h", 0, "logical",
              "workingdir", "p", 1, "character",
              "pcaresultspath", "f", 1, "character",
              "persongroupfilepath", "t", 1, "character",
              "genomes1ksamplepoppath", "w", 1, "character",
              "genomes1krefgroups", "r", 1, "character",
              "plottinggroups", "n", 1, "character",
              "groupassignmentfilepath", "a", 1, "character",
              "rescaledpcaresultspath", "r", 1, "character",
              "assignedgroups", "g", 1, "character"
), byrow=TRUE, ncol=4)



if (is.null(opt$workingdir)){
  print("you must provide a --workingdir value")
  stop()
}

if (is.null(opt$pcaresultspath)){
  print("you must provide a --pcaresultspath value")
  stop()
}

if (is.null(opt$persongroupfilepath)){
  print("you must provide a --pcaresultspath value")
  stop()
}

if (is.null(opt$genomes1ksamplepoppath)){
  print("you must provide a --genomes1ksamplepoppath value")
  stop()
}

if (is.null(opt$genomes1krefgroups)){
  print("you must provide a --genomes1krefgroups value")
  stop()
}


if (is.null(opt$plottinggroups)){
  opt$plottinggroups<-NA
}


if (is.null(opt$groupassignmentfilepath)){
  opt$groupassignmentfilepath<-NA
}

if (is.null(opt$rescaledpcaresultspath)){
  opt$rescaledpcaresultspath<-NA
}

if (is.null(opt$assignedgroups)){
  opt$assignedgroups<-NA
}



working.dir<-opt$workingdir


eth.assignment.file.path<-opt$groupassignmentfilepath


target.ref.added.bed.pcs.path<-opt$pcaresultspath

target.race.file.path<-opt$persongroupfilepath

scaled.pcs.path<-opt$rescaledpcaresultspath

genomes1k.sample.pop.path<-opt$genomes1ksamplepoppath
if (!is.na(opt$assignedgroups))
  assigned.eth.to.make.subset<-unlist(strsplit(opt$assignedgroups, ","))

target.ref.groups<-unlist(strsplit(opt$genomes1krefgroups, ","))


target.race.group.dat<-read.table(target.race.file.path, header=TRUE, stringsAsFactors = FALSE)
colnames(target.race.group.dat)<-c("GRID", "RACE")

dir.create(working.dir)


if (!is.na(opt$plottinggroups)){
  plotting.sample.groups.subsets<-strsplit(opt$plottinggroups, "/")
  plotting.sample.groups.subsets.lst<-list()
  
  for (l in 1:length(plotting.sample.groups.subsets))
  {
    this.group<-plotting.sample.groups.subsets[[l]]
    print(this.group)
    this.group.vec<-unlist(strsplit(this.group, ","))
    
    plotting.sample.groups.subsets.lst<-append(plotting.sample.groups.subsets.lst, this.group.vec)
  }
}else{
  #let's just create a single plotting group
  plotting.sample.groups.subsets.lst<-append(plotting.sample.groups.subsets.lst, unique(target.race.group.dat$RACE))
}

output.file.prefix.path<-file.path(working.dir, paste0(basename(target.ref.added.bed.pcs.path), ".rescaled_pca.resubsetted") )



remap_flashPCA_to_ref<-function(pcs.dat, sample.groups.dat, ref.groups, output.dir, plot.name.prefix, target.PCs, plotting.sample.groups.subsets.lst, plotting.sample.groups.subsets.labels=NULL, legend.title="race", rescaled.PCs.path=NULL)
{
  
  GGpcs.dat<<-pcs.dat
  GGsample.groups.dat<<-sample.groups.dat
  GGref.groups<<-ref.groups
  GGoutput.dir<<-output.dir
  GGplot.name.prefix<<-plot.name.prefix
  GGtarget.PCs<<-target.PCs
  GGplotting.sample.groups.subsets.lst<<-plotting.sample.groups.subsets.lst
  GGplotting.sample.groups.subsets.labels<<-plotting.sample.groups.subsets.labels
  
  pcs.dat<-GGpcs.dat
  sample.groups.dat<-GGsample.groups.dat
  ref.groups<-GGref.groups
  output.dir<-GGoutput.dir
  plot.name.prefix<-GGplot.name.prefix
  target.PCs<-GGtarget.PCs
  plotting.sample.groups.subsets.lst<-GGplotting.sample.groups.subsets.lst
  plotting.sample.groups.subsets.labels<-GGplotting.sample.groups.subsets.labels
  
  if (1==2)
  {
    
    pcs.dat<-flash.pca.pcs.dat
    # pcs.dat<-flash.pca.pcs.test.dat
    
    sample.groups.dat<-demo.dat
    ref.groups<-c("CEU", "YRI", "CHB")    
    output.dir<-"C:\\Users\\straubp\\Documents\\Projects\\MEGA\\mega_recalled_imputation_results\\PCA\\all_sets\\rescaled\\"
    plot.name.prefix<-"test.pca_new"
    target.PCs<-c(1,2)
    plotting.sample.groups.subsets.lst<-list(c("W", "B", "U", "N", "A", "I"), c("W"), c("B"), c("A") )
  }
  
  if (!file.exists(output.dir))
  {
    dir.create(output.dir, recursive=TRUE)
  }
  
  
  colnames(sample.groups.dat)<-c("id", "group")
  colnames(pcs.dat)<-c("id", paste0("PC", 2:ncol(pcs.dat)-1))
  
  ref.sample.ids <- sample.groups.dat[sample.groups.dat[,2] %in% ref.groups, 1]
  test.sample.ids <- sample.groups.dat[!(sample.groups.dat[,2] %in% ref.groups) & (sample.groups.dat[,1] %in% pcs.dat[,1]), 1]
  
  
  test.sample.dat<-sample.groups.dat[sample.groups.dat[,1] %in% test.sample.ids, ]
  ref.sample.dat<-sample.groups.dat[sample.groups.dat[,1] %in% ref.sample.ids, ]
  
  #because we need our plots to have an aspect ratio of 1:1, let's rescale out PCA values so each PC is on the same scale
  pc.scale.dat<-data.frame(PC=c(), min.val=c(), max.val=c())
  pcs.scaled.dat<-data.frame(id=pcs.dat$id)
  
  
  c=1
  #for (c in 1:(ncol(pcs.dat) - 1))
  for (c in 1:length(target.PCs))
  {
    this.PC<-target.PCs[c]
    #this.min<-min(pcs.dat[, this.PC+1])
    #this.max<-max(pcs.dat[, this.PC+1])
    this.min<-min(pcs.dat[, paste0("PC", this.PC)])
    this.max<-max(pcs.dat[, paste0("PC", this.PC)])
    pc.scale.dat<-rbind(pc.scale.dat, data.frame(PC=this.PC, min.val=this.min, max.val=this.max))
    
    
    
    this.pcs.dat<-pcs.dat[,c(1,this.PC+1)]
    this.range<-this.max-this.min
    this.pcs.dat$dist.from.min<-this.pcs.dat[,2] - this.min
    this.pcs.dat$scaled.val<-this.pcs.dat$dist.from.min/this.range
    
    pcs.scaled.dat$scaled.pc<-this.pcs.dat$scaled.val
    colnames(pcs.scaled.dat)[ncol(pcs.scaled.dat)]<-paste0("PC", this.PC)
    
  }
  
  #return(pcs.scaled.dat)
  #break;
  
  
  ref.pcs.dat<-pcs.scaled.dat[pcs.scaled.dat[,1] %in% ref.sample.ids,]
  test.pcs.dat<-pcs.scaled.dat[pcs.scaled.dat[,1] %in% test.sample.ids,]
  
  
  test.sample.pcs.dat<-merge(test.sample.dat, test.pcs.dat, by="id")
  
  
  ref.PC.cntOfMass<-data.frame(pop=c(), PC=c(), avg=c())
  
  #now find the center of mass for each of the reference pops along each of the PC dimensions
  p=1  
  for (p in 1:length(ref.groups))
  {
    
    this.ref.pop.ids<-sample.groups.dat[sample.groups.dat[,2]== ref.groups[p],1]
    this.pop.dat<-pcs.scaled.dat[pcs.scaled.dat[,1] %in% this.ref.pop.ids,]
    
    
    #for (c in 2:11)
    for (c in 1:length(target.PCs))
    {
      this.PC<-target.PCs[c]
      avg.vec=mean(this.pop.dat[, paste0("PC", this.PC)])
      
      ref.PC.cntOfMass<-rbind(ref.PC.cntOfMass,data.frame(pop=ref.groups[p], PC=this.PC, avg=avg.vec))
      
    }
  }
  
  
  #if the rescaled PCs haven't yet been generated, generate them
  if (is.null(rescaled.PCs.path))  
  {
    dist.scaled=data.frame(test.id=c(), ref.pop1=c(), ref.pop2=c(), PC1=c(), PC2=c(), dist.ref1=c(), dist.ref2=c())
    
    line.dat=data.frame(ref.pop1=c(), ref.pop2=c(), PC1=c(), PC2=c(), ref.pop1.cntOfMass.PC1=c(), ref.pop1.cntOfMass.PC2=c(), ref.pop2.cntOfMass.PC1=c(), ref.pop2.cntOfMass.PC2=c(), slope=c(), intcpt=c())
    
    print("finding distances between ref population pairs")
    p=1
    #cycle through each pair of ref populations
    for (p in 1:(length(ref.groups)-1))
    {
      print(paste0("ref pop 1 = ", ref.groups[p]))
      
      p2=2  
      for (p2 in (p+1):length(ref.groups))
      {
        print(paste0("ref pop 2 = ", ref.groups[p2]))
        #cycle through each pair of PCs
        c=1
        #for (c in 1:2)
        for (c in 1:(length(target.PCs) - 1))
        {
          this.PC<-target.PCs[c]
          print(paste0("PC1 = ", this.PC))
          
          
          ref.pop1.PC1.cntOfMass=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p]  & ref.PC.cntOfMass$PC==this.PC,]
          ref.pop2.PC1.cntOfMass=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p2]  & ref.PC.cntOfMass$PC==this.PC,]
          
          c2=2
          for (c2 in (c+1):length(target.PCs))
          {
            this.PC2<-target.PCs[c2]
            print(paste0("PC2 = ", this.PC2))
            
            ref.pop1.PC2.cntOfMass=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p]  & ref.PC.cntOfMass$PC==this.PC2,]
            ref.pop2.PC2.cntOfMass=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p2]  & ref.PC.cntOfMass$PC==this.PC2,]
            
            #find slope
            xdiff=(ref.pop2.PC1.cntOfMass$avg - ref.pop1.PC1.cntOfMass$avg)
            ydiff=(ref.pop2.PC2.cntOfMass$avg - ref.pop1.PC2.cntOfMass$avg)
            slope=ydiff/xdiff
            
            #find intercept
            y1=ref.pop1.PC2.cntOfMass$avg
            x1=ref.pop1.PC1.cntOfMass$avg
            
            y1.2=-y1*xdiff
            x1.2=-x1*ydiff
            
            int1=x1.2 - y1.2
            intcpt=int1/xdiff
            
            line.dat=rbind(line.dat, data.frame(ref.pop1=p, ref.pop2=p2, PC1=this.PC, PC2=this.PC2, ref.pop1.cntOfMass.PC1=ref.pop1.PC1.cntOfMass$avg, ref.pop1.cntOfMass.PC2=ref.pop1.PC2.cntOfMass$avg, ref.pop2.cntOfMass.PC1=ref.pop2.PC1.cntOfMass$avg, ref.pop2.cntOfMass.PC2=ref.pop2.PC2.cntOfMass$avg, slope=slope, intcpt=intcpt))
            
            
            dist.refs=((ref.pop2.PC2.cntOfMass$avg - ref.pop1.PC2.cntOfMass$avg)^2 + (ref.pop2.PC1.cntOfMass$avg - ref.pop1.PC1.cntOfMass$avg)^2)^0.5
            
            
            #now, for each non-ref individual, find the perpendicular that connects that individual to the the line we just created 
            i=1
            
            this.dist.scaled<- do.call(rbind,lapply(1:nrow(test.pcs.dat), function(i, test.pcs.dat, pc, pc2, slope, xdiff, ydiff, intcpt, dist.refs, ref.pop1.PC1.cntOfMass, ref.pop1.PC2.cntOfMass, ref.pop2.PC1.cntOfMass, ref.pop2.PC2.cntOfMass, ref.groups, p, p2)
            {
              
              #ind.x=test.pcs.dat[i, (c + 1)]
              #ind.y=test.pcs.dat[i, (c2 + 1)]
              
              ind.x=test.pcs.dat[i, pc]
              ind.y=test.pcs.dat[i, pc2]
              
              yi.2=-ind.y*-ydiff
              xi.2=-ind.x*xdiff
              
              int.i=xi.2 - yi.2
              intcpt.i=int.i/-ydiff
              
              #now let's find where the lines intersect
              x.cross=(intcpt.i-intcpt)/(slope + 1/slope )
              y.cross=x.cross*slope + intcpt
              
              
              #finally, let's describe the intesection point as a percentage of the distance between our two reference points
              
              dist.ref1=((y.cross - ref.pop1.PC2.cntOfMass$avg)^2 + (x.cross - ref.pop1.PC1.cntOfMass$avg)^2)^0.5
              dist.ref2=((y.cross - ref.pop2.PC2.cntOfMass$avg)^2 + (x.cross - ref.pop2.PC1.cntOfMass$avg)^2)^0.5
              
              dist.ref1.scaled=dist.ref1/dist.refs
              dist.ref2.scaled=dist.ref2/dist.refs
              
              
              #check to see if the intersection point is NOT between the two reference points and correct the sign
              if ((dist.ref1.scaled > 1) | (dist.ref2.scaled > 1))
              {
                if (dist.ref1.scaled > dist.ref2.scaled)
                {
                  dist.ref2.scaled= -dist.ref2.scaled
                }else
                {
                  dist.ref1.scaled= -dist.ref1.scaled
                }
                
              }
              
              
              return(data.frame(test.id=test.pcs.dat[i,1], ref.pop1=ref.groups[p], ref.pop2=ref.groups[p2], PC1=pc, PC2=pc2, dist.ref1=dist.ref1.scaled, dist.ref2=dist.ref2.scaled, stringsAsFactors = FALSE))
              
            }, test.pcs.dat, paste0("PC",this.PC), paste0("PC", this.PC2), slope, xdiff, ydiff, intcpt, dist.refs, ref.pop1.PC1.cntOfMass, ref.pop1.PC2.cntOfMass, ref.pop2.PC1.cntOfMass, ref.pop2.PC2.cntOfMass, ref.groups, p, p2))
            
            dist.scaled<-rbind(dist.scaled, this.dist.scaled)
            
          }
        }
      }
    }
    
    write.table(dist.scaled, file=file.path(output.dir, paste0(plot.name.prefix, ".scaled.pcs.txt")), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
    write.table(line.dat, file=file.path(output.dir, paste0(plot.name.prefix, ".scaled.reference.lines.txt")), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
  }else{
    dist.scaled=read.table(rescaled.PCs.path, header = TRUE, sep="\t", stringsAsFactors = FALSE)
    line.dat=read.table(gsub("pcs.txt$","reference.lines.txt",rescaled.PCs.path), header = TRUE, sep="\t", stringsAsFactors = FALSE)
    
  }
  
  
  c=1
  for (c in 1:(length(target.PCs) - 1))
  {
    this.PC <- target.PCs[c]
    print(paste0("PC1 = ", this.PC))
    c2=2
    for (c2 in (c+1):length(target.PCs))
    {
      this.PC2<-target.PCs[c2]
      print(paste0("PC2 = ", this.PC2))
      
      xmin=min(pcs.scaled.dat[,paste0("PC", this.PC)])
      ymin=min(pcs.scaled.dat[,paste0("PC", this.PC2)])
      
      xmax=max(pcs.scaled.dat[,paste0("PC", this.PC)])
      ymax=max(pcs.scaled.dat[,paste0("PC", this.PC2)])
      
      
      #now, for each of the plotting.sample.groups.subsets.lst groups, make a plot
      g=1
      for (g in 1:length(plotting.sample.groups.subsets.lst))
      {
        this.subset.group<-plotting.sample.groups.subsets.lst[[g]]
        this.subset.group.label<-paste(this.subset.group, collapse=".")
        print(paste0("making plot for ", this.subset.group.label, " group"))
        this.subset.group.plotting.sample.ids<-sample.groups.dat[sample.groups.dat$group %in% this.subset.group, "id"]
        
        plotting.test.sample.pcs.dat<-test.sample.pcs.dat[test.sample.pcs.dat$id %in% this.subset.group.plotting.sample.ids, ]
        
        
        if (length(this.subset.group.plotting.sample.ids) == 0)
        {
          print(paste0("no group data found for ", this.subset.group.label, " group"))
          next;
        }else if (nrow(plotting.test.sample.pcs.dat) == 0)
        {
          print(paste0("no people found in ", this.subset.group.label, " group"))
          next;
          
        }
        
        
        
        if (is.null(plotting.sample.groups.subsets.labels))
        {
          png.path<-file.path(output.dir, paste0(plot.name.prefix, "_scaled_PC", this.PC, "_vs_PC", this.PC2, ".", this.subset.group.label, ".png"))
        }else
        {
          png.path<-file.path(output.dir, paste0(plot.name.prefix, "_scaled_PC", this.PC, "_vs_PC", this.PC2, ".", plotting.sample.groups.subsets.labels[g], ".png"))
        }
        
        png(file=png.path)
        tick.cols=c("red", "green", "blue")
        tick.col.counter=1
        
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        
        #now create a plot with the reference lines and tick marks
        plot(c(xmin,xmax), c(ymin,ymax), type = "n", xlab=paste0("PC",this.PC), ylab=paste0("PC",this.PC2), asp=1)
        
        p=1
        for (p in 1:(length(ref.groups)-1))
        {
          print(paste0("ref pop 1 = ", ref.groups[p]))
          
          p2=2  
          for (p2 in (p+1):length(ref.groups))
          {
            #create the center of mass points for the ref populations and draw the lines connecting them 
            this.line.dat=line.dat[line.dat$ref.pop1==p & line.dat$ref.pop2==p2 & line.dat$PC1==this.PC & line.dat$PC2==this.PC2,]
            
            points(this.line.dat$ref.pop1.cntOfMass.PC1, this.line.dat$ref.pop1.cntOfMass.PC2)
            points(this.line.dat$ref.pop2.cntOfMass.PC1, this.line.dat$ref.pop2.cntOfMass.PC2)
            
            abline(this.line.dat$intcpt, this.line.dat$slope)
            
            #now, draw dotted perpenticular lines at each 10% tick between the ref points
            
            dist.ref=((this.line.dat$ref.pop2.cntOfMass.PC2 - this.line.dat$ref.pop1.cntOfMass.PC2)^2 + (this.line.dat$ref.pop2.cntOfMass.PC1 - this.line.dat$ref.pop1.cntOfMass.PC1)^2)^0.5
            
            ticks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
            t=1
            for (t in 1:length(ticks))
            {
              PC1.tick=((1 - ticks[t])* this.line.dat$ref.pop1.cntOfMass.PC1 + ticks[t]*this.line.dat$ref.pop2.cntOfMass.PC1)
              PC2.tick=((1 - ticks[t])* this.line.dat$ref.pop1.cntOfMass.PC2 + ticks[t]*this.line.dat$ref.pop2.cntOfMass.PC2)
              
              #now find the formula for the perpendicular line that passes through this tick mark
              #get the intercept
              perp.slope=-1/this.line.dat$slope
              
              perp.intcpt=PC2.tick - perp.slope*PC1.tick
              
              lty.type="dotted"
              if (ticks[t]==0.5)
              {
                lty.type="solid"
              }
              
              abline(perp.intcpt, perp.slope,lty=lty.type, col=tick.cols[tick.col.counter])
              
            }
            tick.col.counter=tick.col.counter+1
            
          }
        }
        
        #now create dots for all the non-reference individuals in this set
        
        
        num.race=unique(plotting.test.sample.pcs.dat$group)
        color.pal=cbind(num.race, rainbow(length(num.race)))
        colnames(color.pal)=c("race", "col")
        t=1
        for (t in 1:nrow(plotting.test.sample.pcs.dat))
        {
          this.col=color.pal[color.pal[,"race"]==plotting.test.sample.pcs.dat[t,"group"],"col"]
          
          #points(plotting.test.sample.pcs.dat[t, this.PC+2], plotting.test.sample.pcs.dat[t, this.PC2+2], col=this.col, pch=4)
          points(plotting.test.sample.pcs.dat[t, paste0("PC", this.PC)], plotting.test.sample.pcs.dat[t, paste0("PC", this.PC2)], col=this.col, pch=4)
          
        }
        
        legend("topright", inset=c(-0.2,0), legend=color.pal[,"race"], col=color.pal[,"col"], title=legend.title, pch =4 )
        
        #now add the labels for the reference populations
        p=1
        for (p in 1:length(ref.groups))
        {
          cntofMass1=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p] & ref.PC.cntOfMass$PC==this.PC, "avg"]
          cntofMass2=ref.PC.cntOfMass[ref.PC.cntOfMass$pop==ref.groups[p] & ref.PC.cntOfMass$PC==this.PC2, "avg"]
          
          text(cntofMass1, cntofMass2, ref.groups[p])
          
          
        }
        
        
        dev.off() 
      }
      ref.pop1.pop2.pop1=paste0(ref.groups[1], "_vs_", ref.groups[2],".",ref.groups[1])
      ref.pop1.pop2.pop2=paste0(ref.groups[1], "_vs_", ref.groups[2],".",ref.groups[2])
      
      ref.pop1.pop3.pop1=paste0(ref.groups[1], "_vs_", ref.groups[3],".",ref.groups[1])
      ref.pop1.pop3.pop3=paste0(ref.groups[1], "_vs_", ref.groups[3],".",ref.groups[3])
      
      ref.pop2.pop3.pop2=paste0(ref.groups[2], "_vs_", ref.groups[3],".",ref.groups[2])
      ref.pop2.pop3.pop3=paste0(ref.groups[2], "_vs_", ref.groups[3],".",ref.groups[3])
      
      
      
      #now we need to create the file that will be manually filled out for users to specify cutoffs along our three dimensions to define
      man.fill.dat=data.frame(assignment.group.name="default", first.PC=1, second.PC=2, ref.pop1.pop2_pop1=50, ref.pop1.pop2_pop2=50, ref.pop1.pop3_pop1=50, ref.pop1.pop3_pop3=50, ref.pop2.pop3_pop2=50, ref.pop2.pop3_pop3=50, assigned.eth=NA, assigned.ref.pop.for.refilter=NA, stringsAsFactors = FALSE)
      colnames(man.fill.dat)=c("assignment.group.name", "first.PC", "second.PC", ref.pop1.pop2.pop1, ref.pop1.pop2.pop2, ref.pop1.pop3.pop1, ref.pop1.pop3.pop3, ref.pop2.pop3.pop2, ref.pop2.pop3.pop3, "assigned.eth", "assigned.ref.pop.for.refilter")
      
      assignment.group.path<-file.path(output.dir, paste0(plot.name.prefix, "_scaled_PC", this.PC, "_vs_PC", this.PC2, "_assigned_group_boundaries.txt"))	
      print("writing file to fill out with boundaries: ")
      
      newpath.found<-F
      newpathcntr<-2
      while(!newpath.found)
      {
        if (file.exists(assignment.group.path))
        {
          if (newpathcntr==2)
          {
            assignment.group.path<-gsub(".txt$", ".2.txt", assignment.group.path)  
          }else{
            assignment.group.path<-gsub(paste0(".", newpathcntr-1,".txt$"), paste0(".", newpathcntr,".txt$"), assignment.group.path)
          }
          newpathcntr<-newpathcntr+1
        }else{
          newpath.found<-T
        }
      }
      write.table(man.fill.dat, file=assignment.group.path, sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)
      
    }
  }
  
  
  
  
}








genomes1k.sample.pop.dat<-read.table(genomes1k.sample.pop.path, header=TRUE, sep="", stringsAsFactors = FALSE)


genomes1k.flash.race.dat<-genomes1k.sample.pop.dat[, c("ID","GROUP", "POP")]
genomes1k.flash.race.dat$POP.GROUP<-paste0(genomes1k.flash.race.dat$GROUP,"-",genomes1k.flash.race.dat$POP)
colnames(genomes1k.flash.race.dat)<-c("GRID", "GROUP", "POP","RACE")



#target.race.group.dat<-mega.race.eth.dat
target.genomes1k.race.group.dat<-rbind(target.race.group.dat[,c("GRID", "RACE")], genomes1k.flash.race.dat[,c("GRID", "RACE")])


flash.pca.pcs.dat<-read.table(target.ref.added.bed.pcs.path, header=TRUE, sep="", stringsAsFactors = FALSE)
flash.pca.pcs.dat<-flash.pca.pcs.dat[,-1]




#if there's an existing group assignment file, we need to read that in
if (!is.na(eth.assignment.file.path))
{
  
  if (!is.na(scaled.pcs.path))
  {
    #read in the scaled pca values for the target samples
    scaled.pcs.dat<-read.table(scaled.pcs.path, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    
  }else{
    print("You provided a group assignment file from a previous run but no rescaled PCs")
    stop()
    
  }
  
  
  #read in the ethnicity assignments
  eth.assignment.dat<-read.table(eth.assignment.file.path, header=TRUE, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
  
  
  #now let's limit it to just hte se
  
  #now go through each filter and get the people that meet that filter criteria
  
  r=1
  for (r in 1:nrow(eth.assignment.dat))
  {
    #r=1
    this.filt<- eth.assignment.dat[r,]
    print(paste0("processing ", this.filt))  
    
    
    #go through the different reference population axes and apply each filter
    c=4
    this.filtered.ids<-unique(scaled.pcs.dat$test.id)
    
    for (c in 4:9)
    {
      this.filt.val<-as.numeric(this.filt[c])
      print(c)
      
      #if the cutoff filter is NA or otherwise not numeric, skip it
      if (is.na(this.filt.val) |  !is.numeric(this.filt.val))
        next;
      
      
      this.axis<-colnames(eth.assignment.dat)[c]
      
      this.axis.pops<-unlist(strsplit(this.axis, "_vs_|\\."))
      this.ref.pop1<-this.axis.pops[1]
      this.ref.pop2<-this.axis.pops[2]
      this.target.ref.pop<-this.axis.pops[3]
      this.filt.cutoff<-this.filt.val/100
      
      #now get the scaled PC data for these reference pops
      this.scaled.pc.dat<-subset(scaled.pcs.dat, (ref.pop1 %in% c(this.ref.pop1, this.ref.pop2)) &  (ref.pop2 %in% c(this.ref.pop1, this.ref.pop2) & (PC1 %in% c(this.filt$first.PC, this.filt$second.PC)) & (PC2 %in% c(this.filt$first.PC, this.filt$second.PC)))   )
      this.scaled.ref.pop1<-unique(this.scaled.pc.dat$ref.pop1)
      this.scaled.ref.pop2<-unique(this.scaled.pc.dat$ref.pop2)
      
      #get the orientation of this filter (i.e. if we should use the dist.ref1 or dist.ref2 columns from the scaled data)
      if (this.target.ref.pop==this.scaled.ref.pop1)
      {
        this.filtered.scaled.pc.dat<-subset(this.scaled.pc.dat, dist.ref1 <= this.filt.cutoff)
      }else
      {
        this.filtered.scaled.pc.dat<-subset(this.scaled.pc.dat, dist.ref2 <= this.filt.cutoff)
      }
      
      #now filter out filtered list down to just those ids that met this filter
      this.filtered.ids<-this.filtered.ids[this.filtered.ids %in% this.filtered.scaled.pc.dat$test.id]
      
    }
    
    this.filtered.ids.dat<-data.frame(FID=0, IID=this.filtered.ids, stringsAsFactors = FALSE)
    
    #now get the selected reference pop.  Might be a comma-delimited list
    this.selected.ref.pops<-this.filt$assigned.ref.pop.for.refilter
    this.selected.ref.pops<-unlist(strsplit(this.selected.ref.pops, ","))
    
    this.ref.pop.ids<-subset(genomes1k.flash.race.dat, RACE %in% this.selected.ref.pops)
    this.filtered.ids.and.refs.dat<-rbind(this.filtered.ids.dat, data.frame(FID=this.ref.pop.ids$GRID, IID=this.ref.pop.ids$GRID, stringsAsFactors = FALSE))
    
    output.file.path<-paste0(output.file.prefix.path, ".", this.filt$assigned.eth, ".refilter.keep.txt")  
    write.table(this.filtered.ids.and.refs.dat, file=output.file.path, col.names = FALSE, row.names = FALSE, quote=FALSE, sep="\t")  
    
    #let's also make a keep file without the reference populations that we'll use for running PCA on the selected target popuation alone
    output.file.no.ref.path<-paste0(output.file.prefix.path, ".", this.filt$assigned.eth, ".refilter.no.refs.keep.txt")  
    write.table(this.filtered.ids.dat, file=output.file.no.ref.path, col.names = FALSE, row.names = FALSE, quote=FALSE, sep="\t")  
    
    
  }
  
  
  
  e=1
  for (e in 1:length(assigned.eth.to.make.subset))
  {
    #e=1
    this.filt.dat<-subset(eth.assignment.dat, assigned.eth==assigned.eth.to.make.subset[e])
    
    this.eth.ids.path<-file.path(paste0(output.file.prefix.path, ".", assigned.eth.to.make.subset[e], ".refilter.keep.txt")  )
    this.eth.ids.dat<-read.table(this.eth.ids.path, header=FALSE, sep="", stringsAsFactors = FALSE)
    
    #first, let's make the rescaled PCA with CEU/YRI/CHB but showing the selected group
    this.eth.ids.non.ref.dat<-this.eth.ids.dat[!(this.eth.ids.dat$V2 %in% genomes1k.sample.pop.dat[,1]),]
    target.genome1k.race.selected.dat<-target.genomes1k.race.group.dat
    target.genome1k.race.selected.dat[target.genome1k.race.selected.dat[,1] %in% this.eth.ids.non.ref.dat[,2], "RACE"]<-"selected"
    remap_flashPCA_to_ref(pcs.dat=flash.pca.pcs.dat, sample.groups.dat=target.genome1k.race.selected.dat[,c("GRID", "RACE")], ref.groups=target.ref.groups, output.dir=working.dir, plot.name.prefix=paste0(basename(target.ref.added.bed.pcs.path), ".",assigned.eth.to.make.subset[e], ".show.selected"), target.PCs=c(1,2), plotting.sample.groups.subsets.lst=plotting.sample.groups.subsets.lst, legend.title="population", rescaled.PCs.path = scaled.pcs.path)
    
    
    
  }
  
  
  
  
}else #otherwise, just rescale this PCA set as is without any "selected" set
{
  remap_flashPCA_to_ref(pcs.dat=flash.pca.pcs.dat, sample.groups.dat=target.genomes1k.race.group.dat[,c("GRID", "RACE")], ref.groups=target.ref.groups, output.dir=working.dir, plot.name.prefix=paste0(basename(target.ref.added.bed.pcs.path), ".nonselected"), target.PCs=c(1,2), plotting.sample.groups.subsets.lst=plotting.sample.groups.subsets.lst, legend.title="population")
  
}









