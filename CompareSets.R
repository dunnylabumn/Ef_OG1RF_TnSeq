#!!!!! Monte Carlo Simulation is drawing on 24 CORES !!!!!!!!!!
MonteCarlo = TRUE
Monte_Limit = 10
CPUcount = 24

#6700_biofilm_1_GATCAG_L002_R1_001.fastq.all.tab.annot   4
#6700_biofilm_2_GTTTCG_L002_R1_001.fastq.all.tab.annot   4
#6700_input_1_ACTTGA_L002_R1_001.fastq.all.tab.annot     5
#6700_input_2_GTGGCC_L002_R1_001.fastq.all.tab.annot     5
#6700_planktonic_1_TAGCTT_L002_R1_001.fastq.all.tab.annot        6
#6700_planktonic_2_CGTACG_L002_R1_001.fastq.all.tab.annot        6

#Rscript --vanilla sample_def_176k.tab results_176k.tab
args = commandArgs(trailingOnly=TRUE)
# test if there is at least two argument: if not, return an error
if (length(args)<1) {
  stop("At least one argument (comparision file) must be supplied (comparision file), [output file]", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "comparison_data.tab"
}


# "EggNOG"        "locus_tag"     "old_locus_tag" "COG"          
#[5] "minPos"        "maxPos"        "Orient"        "gene_id"      
#[9] "NCBI_desc"     "protein_id"    "gi_id"         "pfam_id"      
#[13] "interpro_id"   "EC_num"        "unknown"       "COG_desc"     
#[17] "unknown2"      "NOG_desc"      "Position"      "TA_site"      
#[21] "Integrations" 

#Read in filename
read.specific = function(z) tryCatch(
  read.table(as.character(z), quote="", header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char=''),
  error=function(e) {print(paste("Problem with file",z,e));stop()} )

#Input_list is a two column file with filepaths in the first column and groupIDs in the second column.
input_list = read.table(args[1], quote="", header=FALSE, sep="\t")
#print(input_list)
filesets=split(input_list[,1], input_list[,2])
input_data = lapply(input_list[,1], read.specific)
input_data = lapply(input_data, function(x){ x[which(!is.na(x[,"Position"])),]})

#Group info for comparisons
group_divisions = input_list[,2]
group_sets = unique(group_divisions)
paired_combin = combn(group_sets,2)


#Helper Functions
converge_dat = function(x){
  total_length = length(x)
  allnames <- unique(unlist(lapply(x, names)))
  base_dat = do.call(rbind, lapply(x, function(df) {
  	return(df[,-21])
  }))
  base_dat =unique(base_dat)
  
  names(base_dat) = names(x[[1]])[-21]
  return(base_dat)
}

cov_fun = function(y){
	return(sum(y>0))
}

#Get unique list of all positions. Merge all input datasets against this unique list.
base_data = converge_dat(input_data)
f_col_name = paste("Integrations", input_list[1,1],sep="")
for(i in 1:length(input_data)){
	base_data = tryCatch(merge(x = base_data, y = input_data[[i]][,c("Position","Integrations")], by = "Position", all.x = TRUE,suffixes=c("", as.character(input_list[i,1]))),error=function(e) { print(i); print(base_data[1:10,]); print(input_data[[i]][1:10,c("Position","Integrations")]); stop(e)})
}
base_data = unique(base_data)
colnames(base_data)[which(colnames(base_data) == "Integrations")] = f_col_name

#Number of non-TA sites with an integration
off_rate = base_data[which(is.na(base_data[,"TA_site"])),]
off_rate = off_rate[,21:ncol(off_rate)]
off_rate[is.na(off_rate)] = 0
off_rate = apply(off_rate,2,cov_fun)

int_sites = base_data[which(!is.na(base_data[,"TA_site"])),]

total_integrations = colSums(int_sites[21:ncol(int_sites)],na.rm=TRUE)

#number of potential non-TA sites within annotated regions
off_site = sum(base_data[,"maxPos"] - base_data[,"minPos"]) - sum(!is.na(base_data[,"TA_site"]))

#percentage of non-TA sites with Integrations (off targetting rate)
off_rate = off_rate/off_site

#Split data be annotated regions
spl_dat = split(base_data, base_data[,"old_locus_tag"])


#Computes odds for all samples. This is a test of each dataset against the NULL distribution of random background.
# gene_summarize = function(x){
# 	offsite = x[which(is.na(x[,"TA_site"])),]
# 	offsite = offsite[,21:ncol(offsite)]
# 	offsite[is.na(offsite)] = 0
# 
# 	onsite = x[which(!is.na(x[,"TA_site"])),]
# 	onsite = onsite[,21:ncol(onsite)]
# 	onsite[is.na(onsite)] = 0
# 
# 	gene_size = x[1,7] - x[1,6] + 1
# 
# 	site_cov = apply(onsite,2, cov_fun)
# 	back_cov = apply(offsite,2,cov_fun)
# 
# 	TA_count = nrow(onsite)
# 	bg_count = nrow(offsite)
# 	
# 	all_cov = mapply(c, site_cov, off_rate, SIMPLIFY=FALSE)
# 
# 	paired_phyper = function(x,TAs, gene_len){
# 		#ret_res = phyper(x[1], (gene_len-TAs)*(1-x[2]), (gene_len-TAs)*x[2], TAs, lower.tail=FALSE)
# 		ret_res = pbinom(x[1], size=TAs, prob=x[2], lower.tail=FALSE)
# 		return(ret_res)
# 	}
# 
# 	prob_res = lapply(all_cov, paired_phyper, TAs=TA_count, gene_len=gene_size)	
#Rscript --vanilla /home/support/jballer/Dunny/transposon_pipeline/CompareSets.R $basedir/sample_def_6700.tab $basedir/results_6700.tab
#Rscript --vanilla /home/support/jballer/Dunny/transposon_pipeline/CompareSets.R $basedir/sample_def_6700.tab $basedir/results_6700.tab
#Rscript --vanilla /home/support/jballer/Dunny/transposon_pipeline/CompareSets.R $basedir/sample_def_6700.tab $basedir/results_6700.tab
#Rscript --vanilla /home/support/jballer/Dunny/transposon_pipeline/CompareSets.R $basedir/sample_def_6700.tab $basedir/results_6700.tab
#   	prob_res = as.numeric(unlist(prob_res))
# 	return(c(x[1,1:20], site_cov, TA_sites=TA_count, prob_res*total_tests))  
# }

#out_dat = lapply(spl_dat, gene_summarize) # Calculates per gene odds


#Monte Carlo simulation to compare each data group against all other data groups.
monte_carlo_sim = function(dist_list_a, dist_list_b, reps=1000, tails="lower"){
  #tails options upper, lower, either
  #dist_list Structure:
  #list(list(dist_function_1,par1=par1_val, par2=par2_val, ...) ,list(dist_function_2,par1=par1_val, par2=par2_val, ...))
  #The first parameter to dist_function_1/2 should be the number of reps 
	inner_reps = 10000000
	proc_iter = function(x){
		proc_dist = function(y){
			y.minus = y
			y.minus[[1]] = NULL
			return(do.call(y[[1]], c(inner_reps, y.minus)))
		}
		prob_res = lapply(x, proc_dist)
		prob_res = Reduce('+', prob_res)/length(prob_res)
		return(prob_res)
	}
	inner_iter = function(){
		a_mean = proc_iter(dist_list_a)
		b_mean = proc_iter(dist_list_b)

		if(tails=="lower"){
			in_res=a_mean<b_mean
			in_res=sum(!in_res)
		}else if(tails=="upper"){
			in_res=a_mean>b_mean
			in_res=sum(!in_res)
		}else if(tails == "either"){
			if(sum(a_mean<b_mean,na.rm=TRUE) > sum(a_mean>b_mean,na.rm=TRUE)){
				in_res =  a_mean<b_mean
			}else{
				in_res =  a_mean>b_mean
			}
  			in_res=sum(!in_res)
		}
		return(in_res)
	}
	
	bin_res=sum(replicate(reps, inner_iter()),na.rm=TRUE)	
	if(bin_res  == 0){
		bin_res = 1
	}
	if(tails == "either"){
		return((bin_res/(reps*inner_reps))*2)
	}else{
		return(bin_res/(reps*inner_reps))
	} 
}

#Pass base_data as x 
comp_bp = function(x, num_reps=NULL){
  onsite = x[which(!is.na(x[,"TA_site"])),]

  per_bp_proc = function(bp_line){
    on_target = as.numeric(bp_line[21:length(bp_line)])
    on_target[is.na(on_target)] = 0
    
    #Make array for output of all comparisions of a single annotated regions
    final_stat = array(0, dim=ncol(paired_combin))
    names(final_stat) = apply(paired_combin,2, function(x) return(paste(group_sets[x],collapse="_vs_")))  
    
    final_stat2 = array(0, dim=ncol(paired_combin))
    names(final_stat2) = apply(paired_combin,2, function(x) return(paste(group_sets[x],collapse="_vs_")))
    
    #Iterate across all comparisions
    for(i in 1:ncol(paired_combin)){
      #browser()
      
      iA = paired_combin[1,i]
      iB = paired_combin[2,i]
      
      sel_A = which(group_divisions == iA)
      sel_B = which(group_divisions == iB)
      
      set_A = on_target[sel_A]
      set_B = on_target[sel_B]
      
      tot_A = total_integrations[sel_A]
      tot_B = total_integrations[sel_B]
      
      all_A = data.frame(rbind(set_A, tot_A))
      all_B = data.frame(rbind(set_B, tot_B))
     
 
      #Pass to monte carlo sim two distributions to be compared.
      if((sum(all_A)+sum(all_B) >=Monte_Limit) && MonteCarlo ){
        final_stat[i] = monte_carlo_sim(lapply(all_A, function(x) return(list('rbeta', shape1=x[1]+1, shape2=(x[2]-x[1])+1))), lapply(all_B, function(x) return(list('rbeta', shape1=x[1]+1, shape2=(x[2]-x[1])+1))), reps=10, tails="either")
      }else{
        final_stat[i] = NA
      }
      final_stat2[i] = chisq.test(matrix(c(sum(set_A),sum(tot_A)-sum(set_A), sum(set_B), sum(tot_B)-sum(set_B)),nrow=2))$p.value
    }
   
    val_arr = as.numeric(bp_line[21:length(bp_line)]) 
    names(val_arr) = names(bp_line[21:length(bp_line)])
    return(c(bp_line[c(1,3:4,6:11,13)], val_arr/as.numeric(total_integrations) ,final_stat*ncol(paired_combin)*num_reps,final_stat2*ncol(paired_combin)*num_reps))
    #return(c(bp_line[c(1,3:4,6:11,13)], bp_line[21:length(bp_line)] ,as.numeric(bp_line[21:length(bp_line)])/as.numeric(total_integrations) ,final_stat2*ncol(paired_combin)*num_reps))
  }  
  out_dat = as.data.frame(t(sfApply(onsite, 1, per_bp_proc)))
  return(out_dat)
}
 
comp_genes = function(full_set, num_reps=NULL){
  onsite = full_set[which(!is.na(full_set[,"TA_site"])),]

  onsite_spl = split(onsite, paste(onsite[,"locus_tag"],onsite[,"old_locus_tag"]))
  test_count = length(onsite_spl)
  
  
  per_gene_proc = function(gene_set){
  	on_target = gene_set[,21:ncol(gene_set)]
   	on_target[is.na(on_target)] = 0
    
  	#Observed Sites
  	site_cov = apply(on_target,2, sum)
  	gene_size = gene_set[1,7] - gene_set[1,6] +1
  	
  	#Max Sites
  	TA_count = nrow(gene_set) 
  
  	#Make array for output of all comparisions of a single annotated regions
  	final_stat = array(0, dim=ncol(paired_combin))
  	names(final_stat) = apply(paired_combin,2, function(x) return(paste(group_sets[x],collapse="_vs_")))  
    
  	final_stat2 = array(0, dim=ncol(paired_combin))
  	names(final_stat2) = apply(paired_combin,2, function(x) return(paste(group_sets[x],collapse="_vs_")))
  	
  	#Iterate across all comparisions
  	for(i in 1:ncol(paired_combin)){
  		iA = paired_combin[1,i]
  		iB = paired_combin[2,i]
     
  		sel_A = which(group_divisions == iA)
  		sel_B = which(group_divisions == iB)
  		
  		tot_A = total_integrations[sel_A]
  		tot_B = total_integrations[sel_B]
      
  		set_A = site_cov[sel_A]
  		set_B = site_cov[sel_B]
  	  
  		#Pass to monte carlo sim two distributions to be compared.
  		if(MonteCarlo && (sum(set_A) + sum(set_B)) >= Monte_Limit){
			final_stat[i] = monte_carlo_sim(lapply(set_A, function(x) return(list('rbeta', shape1=x+1, shape2=(TA_count-x)+1))), lapply(set_B, function(x) return(list('rbeta', shape1=x+1, shape2=(TA_count-x)+1))), reps=10, tails="either")
		}else{
			final_stat[i] = NA
		}
  		final_stat2[i] = tryCatch(chisq.test(matrix(c(sum(set_A),sum(tot_A)-sum(set_A), sum(set_B), sum(tot_B)-sum(set_B)),nrow=2))$p.value, error=function(e){
  		  print(matrix(c(sum(set_A),TA_count*length(set_A)-sum(set_A), sum(set_B), TA_count*length(set_B)-sum(set_B)),nrow=2));
  		  print(gene_set);
  		  stop(e)
  		})
  			
  	}  
  	#return(c(x[1,c(2:15,17,19)],TA_count=TA_count,site_cov,final_stat*ncol(paired_combin)*num_reps,final_stat2*ncol(paired_combin)*num_reps))
  	midpoint = (gene_set[1,7] + gene_set[1,6])/2
	front_score = colMeans(on_target[which(gene_set[,"Position"] < midpoint),])
	back_score = colMeans(on_target[which(gene_set[,"Position"] >= midpoint),])
  	return(c(gene_set[1,c(2:15,17,19)],TA_count=TA_count,site_cov, site_cov/total_integrations, final_stat*ncol(paired_combin)*num_reps,final_stat2*ncol(paired_combin)*test_count, int_mean_5=front_score,int_mean_3=back_score))
  }
  out_dat = sfLapply(onsite_spl, per_gene_proc)
  return(do.call("rbind", out_dat))
}
 
library('snowfall')
if(!MonteCarlo) sfInit( parallel=FALSE) else sfInit(parallel=TRUE, cpus=CPUcount)
sfExport('comp_genes')
sfExport('comp_bp')
sfExport('monte_carlo_sim') 
sfExport('cov_fun')
sfExport('paired_combin')
sfExport('group_divisions')
sfExport('group_sets')
sfExport('total_integrations')
sfExport('MonteCarlo','Monte_Limit','base_data')

#Multiprocess across annotated regions 
#out_dat= comp_bp(base_data, num_reps=sum(!is.na(base_data[,"TA_site"])))
out_dat= comp_bp(base_data, num_reps=round(sum(!is.na(base_data[,"TA_site"])/1000)))


######out_dat = do.call("rbind", out_dat)
write.table(out_dat, paste(args[2],".per_bp.tab", sep=""), row.names=FALSE, sep="\t")



#out_dat= comp_genes(base_data, num_reps=sum(!is.na(base_data[,"TA_site"])))
out_dat= comp_genes(base_data, num_reps=round(sum(!is.na(base_data[,"TA_site"])/1000)))


write.table(out_dat, paste(args[2],".per_gene.tab", sep=""), row.names=FALSE, sep="\t")

sfStop()
