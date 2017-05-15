#----------------------------------------------------------
#----------------------------------------------------------
# Notes or References:
# 
#----------------------------------------------------------
########
# define filepaths and create directory structure
########
proj_path<-"/PATH/TO/DIR/WWTP_Impact_on_Stream/"
data_path<-file.path(proj_path,"data")
figs_path<-file.path(proj_path,"figs")
if(!file_test("-d",figs_path)) dir.create(figs_path)
R_path<-file.path(proj_path,"R")
if(!file_test("-d",R_path)) dir.create(R_path)
output_path<-file.path(proj_path,"output")
if(!file_test("-d",output_path)) dir.create(output_path)
rds_path<-file.path(output_path,"RDS")
if(!file_test("-d",rds_path)) dir.create(rds_path)
#----------------------------------------------------------
########
# Set Seed and Options
########
set.seed(20141004)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########