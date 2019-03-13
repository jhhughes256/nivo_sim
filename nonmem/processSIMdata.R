  
dir.extension <- ".nm7"
processSIMdata <- function(model.path.file.ext)
   {
    #Debug model.path.file.ext <- "NM14/base_CRCLCL_min_VPC.ctl"
    
    #Work out the nonmem model dir and model file name
    model.path.file <- gsub(".ctl", "", model.path.file.ext)
    model.path <- dirname(model.path.file)
    model.file.ext <- basename(model.path.file.ext)
    model.file <- gsub(".ctl","",model.file.ext)
    
    model.dir <- paste(master.dir,"/",model.path, sep="")
    
    #Work out the nonmem output dir
    output.dir <- paste(master.dir,"/",model.path.file,dir.extension, sep="")
    setwd(output.dir)
    
    #Work out the SIM fit filename
    SIM.file.ext <- paste(model.file,".fit",sep="")
    #Work out the SIM fit filename
    SIM.file.ext.out <- paste(model.file,".fit.csv",sep="")
    
    
    #Need to remove lines with "TABLE NO. 1 in them" and the header lines for each new subject
    #Strip the unwanted lines
    indata <- readLines(SIM.file.ext)
    tablelines <- grep("TABLE NO.  1",indata)  #May be installation specific
    headerlines <- grep(" ID",indata)          #May be installation specific 
    
    #Extract the information in the header line for column names
    header <- indata[headerlines[1]]
    header <- scan(textConnection(header), what="char")
    colnum <- length(header)
    
    #Strip out the unwanted lines
    badlines <- c(tablelines,headerlines)
    indata <- indata[-badlines]
    
    #replace white space with commas
    for (i in 1:length(indata))
    {
      indata[i] <- gsub('[[:space:]]+',',',indata[i])
    }
    
    #write to a file
    writeLines(indata, "SIMtemp.txt")
    
    #read again as csv file
    SIMdata <- read.csv("SIMtemp.txt", header=F)
    SIMdata <- SIMdata[,-1]     #delete first blank column
    names(SIMdata) <- header
    
    #The NONMEM output does not make a new ID number for each simulation
    #Therefore make a list of SIM numbers
    nsims <- length(tablelines)
    numtimepoints <- length(SIMdata$ID)
    numtimespersubject <- numtimepoints/nsims
    SIMdata$SIM <- rep(1:nsims, each=numtimespersubject)
    
    #Write the SIMdata to a file
    write.csv(SIMdata,SIM.file.ext.out, row.names=F)
    
    #tidy up
    #consider deleting the original fit file?
    file.remove("SIMtemp.txt")
    setwd(master.dir)
    
    #SIMdata
   }