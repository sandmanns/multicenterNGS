library(shiny)
library(BSgenome)
library(DT)
load(file=("./ofInterest.RData"))

shinyServer(function(input, output) {


    #Inital results:
    #coverage
    progress1 <- shiny::Progress$new()
    progress1$set(message = "Analyzing Coverage", value = 0)
    all_means<-data.frame()
    cov_threshold<-100
    for(i in 1:11){
        progress1$inc(1/11)
        binary<-input_data[[1]][[i]]
        binary[input_data[[1]][[i]]<cov_threshold]<-0
        binary[input_data[[1]][[i]]>=cov_threshold]<-1
        binary[is.na(input_data[[1]][[i]])]<-0
        binary_sum<-rowSums(binary,na.rm=T)
        if(i==1){
            all_means<-binary_sum
        }
        if(i>1){
            all_means<-rbind(all_means,binary_sum)
        }
    }
    target_cov<-input_data[[5]]
    genes<-unique(target_cov[,4])
    info<-data.frame(genes,start=NA,end=NA,length=NA,stringsAsFactors=F)
    for(i in 1:length(genes)){
        relevant<-target_cov[target_cov[,4]==genes[i],]
        if(i==1){
            start<-1
            end<-sum(relevant[,3]-relevant[,2]+1)
        }
        if(i>1){
            start<-end+1
            end<-start+sum(relevant[,3]-relevant[,2]+1)-1
        }
        info[i,2]<-start
        info[i,3]<-end
        info[i,4]<-end-start+1
    }
    
    coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                            V8=NA,V9=NA,V10=NA,V11=NA)
    
    for(i in 1:length(coverage_quality[,1])){
        coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
    }

    score_cov<-coverage_quality
    overview<-data.frame(Gene=info[,1],Length=info[,4],Coverage_Score=1-rowSums(score_cov[,5:15],na.rm = T)/831)
    progress1$close()
    
    #MQ
    progress2 <- shiny::Progress$new()
    progress2$set(message = "Analyzing Mapping Quality", value = 0)
    all_means<-data.frame()
    mq_threshold<-58
    for(i in 1:11){
        progress2$inc(1/11)
        binary<-input_data[[4]][[i]]
        binary[input_data[[4]][[i]]<mq_threshold]<-0
        binary[input_data[[4]][[i]]>=mq_threshold]<-1
        binary[is.na(input_data[[4]][[i]])]<-0
        binary_sum<-rowSums(binary,na.rm = T)
        if(i==1){
            all_means<-binary_sum
        }
        if(i>1){
            all_means<-rbind(all_means,binary_sum)
        }
    }
    
    coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                            V8=NA,V9=NA,V10=NA,V11=NA)
    
    for(i in 1:length(coverage_quality[,1])){
        coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
    }

    score_mq<-coverage_quality
    overview<-cbind(overview,MQ_Score=1-rowSums(score_mq[,5:15],na.rm = T)/831)
    progress2$close()
    
    
    #BQ
    progress3 <- shiny::Progress$new()
    progress3$set(message = "Analyzing Base Quality", value = 0)
    all_means<-data.frame()
    bq_threshold<-30
    for(i in 1:11){
        progress3$inc(1/11)
        binary<-input_data[[3]][[i]]
        binary[input_data[[3]][[i]]<bq_threshold]<-0
        binary[input_data[[3]][[i]]>=bq_threshold]<-1
        binary[is.na(input_data[[3]][[i]])]<-0
        binary_sum<-rowSums(binary,na.rm = T)
        if(i==1){
            all_means<-binary_sum
        }
        if(i>1){
            all_means<-rbind(all_means,binary_sum)
        }
    }
    
    coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                            V8=NA,V9=NA,V10=NA,V11=NA)
    
    for(i in 1:length(coverage_quality[,1])){
        coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
    }

    score_bq<-coverage_quality
    overview<-cbind(overview,BQ_Score=1-rowSums(score_bq[,5:15],na.rm = T)/831)
    progress3$close()
    
    #Noise
    progress4 <- shiny::Progress$new()
    progress4$set(message = "Analyzing Background Noise", value = 0)
    all_means<-data.frame()
    bn_threshold<-1
    for(i in 1:11){
        progress4$inc(1/11)
        binary<-input_data[[2]][[i]]
        binary[input_data[[2]][[i]]<(100-bn_threshold)/100]<-0
        binary[input_data[[2]][[i]]>=(100-bn_threshold)/100]<-1
        binary[is.na(input_data[[2]][[i]])]<-0
        binary_sum<-rowSums(binary,na.rm = T)
        if(i==1){
            all_means<-binary_sum
        }
        if(i>1){
            all_means<-rbind(all_means,binary_sum)
        }
    }
    
    coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                            V8=NA,V9=NA,V10=NA,V11=NA)
    
    for(i in 1:length(coverage_quality[,1])){
        coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
    }

    score_bn<-coverage_quality
    overview<-cbind(overview,Noise_Score=1-rowSums(score_bn[,5:15],na.rm = T)/831)
    progress4$close()
    
    
    score_bn<-as.numeric(overview[,6])
    score_bq<-as.numeric(overview[,5])
    score_mq<-as.numeric(overview[,4])
    score_cov<-as.numeric(overview[,3])
    
    qual_score<-rowSums(cbind(score_bn,score_bq,score_mq),na.rm = T)/3
    mean_difficulty<-(qual_score+score_cov)/2
    
    seq_simp<-mean_difficulty
    qual_simp<-qual_score
    cov_simp<-score_cov

    overview<-cbind(overview,Qual_Score=round(qual_score,digits = 4),Seq_Score=round(seq_simp,digits = 4))
    
    clin_impact<-c()
    clin_impact[1]<-14.56
    clin_impact[2]<-2.36
    clin_impact[3]<-2.47
    clin_impact[4]<-1.15
    clin_impact[5]<-9.54
    clin_impact[6]<-1.16
    clin_impact[7]<-5.19
    clin_impact[8]<-1.94
    clin_impact[9]<-1.16
    clin_impact[10]<-2.5
    clin_impact[11]<-3.98
    clin_impact[12]<-9.31
    clin_impact[13]<-1.58
    clin_impact[14]<-1.56
    clin_impact[15]<-4.14
    clin_impact[16]<-1.42
    clin_impact[17]<-1.68
    clin_impact[18]<-1.37
    clin_impact[19]<-5.71
    clin_impact[20]<-1.1
    clin_impact[21]<-8.92
    clin_impact[22]<-23.42
    clin_impact[23]<-1.37
    clin_impact[24]<-1.37
    clin_impact[25]<-7.43
    clin_impact[26]<-3.91
    clin_impact[27]<-18.22
    clin_impact[28]<-8.48
    clin_impact[29]<-7.18
    clin_impact[30]<-0.61
    clin_impact[31]<-5.1
    overview<-cbind(overview,Mutation_Frequency=clin_impact)
    
    output$table <- renderDataTable(datatable(overview))
    output$table1 <- renderDataTable(datatable(overview))
    
    results<-cbind(overview[,1:8],cov_simp,qual_simp,overview[,9])
    #results<-results[order(results[,11],results[,8]),]
    results[,9]<-as.numeric(results[,9])
    results[,10]<-as.numeric(results[,10])
    
    result_backup<-results
    
    color_inside<-c()
    color_outside<-c()
    for(i in 1:length(results[,1])){
        if(results[i,9]<0.8){
            color_inside[i]<-"red"
        }
        if(results[i,9]>=0.8){
            color_inside[i]<-"black"
        }
        if(results[i,10]<0.8){
            color_outside[i]<-"lightblue"
        }
        if(results[i,10]>=0.8){
            color_outside[i]<-"black"
        }
    }

    output$plot<-output$plot1<-renderPlot({
        par(mar=c(4,4,0.5,0.5))
        plot(results[,8],results[,11],xlab="Expected rate of bases with minimum requirements [%]",
             ylab="Mutation frequency [%]",xaxt="n",
             xlim=c((round(min(result_backup[,8]),digits = 2)-0.05),(round(max(result_backup[,8]),digits = 2)+0.05)),
             ylim=c(0,min(100,ceiling(max(results[,11])))),
             cex.lab=1.5,pch=21,col=color_outside,bg=color_inside,lwd=3.5,cex=2)
        legend(x=max(result_backup[,8])-0.32,y=max(result_backup[,11])-0.1,
               legend = c("Problematic coverage and quality","Problematic coverage",
                          "Problematic quality","Unproblematic coverage and quality"),
               pch=21,bty="n",pt.bg=c("red","red","black","black"),pt.cex = 1.7,pt.lwd = 3.5,
               col = c("lightblue","black","lightblue","black"),cex=1.3)
        axis(1,at=seq(0,1,0.1),labels=seq(0,100,10),cex.axis=1.3)
        text(results[,8],results[,11],results[,1],
             pos=c(3,3,2,1,3,1,3,3,3,3,
                   3,3,2,1,3,3,3,3,3,3,
                   3,3,1,3,3,3,3,3,3,1,
                   3),cex=1.2,offset = 0.7)
    },width=1000,height=1000)
    
    output$covUI1<-renderUI({
        conditionalPanel(
            condition="input.cov!='Default'&&input.evaluationMethod=='Absolute'",
            numericInput('cov_threshold','Threshold [0x;20,000x]',min=0,max=20000,value=100)
        )})
    output$covUI1r<-renderUI({
        conditionalPanel(
            condition="input.cov!='Default'&&input.evaluationMethod=='Relative'",
            sliderInput('cov_threshold_r','Threshold (percentile)',min=0,max=100,value=5)
        )})
    
    output$mqUI1<-renderUI({
        conditionalPanel(
            condition="input.mq!='Default'&&input.evaluationMethod=='Absolute'",
            numericInput('mq_threshold','Threshold [0;70]',min=0,max=70,value=58)
        )})
    output$mqUI1r<-renderUI({
        conditionalPanel(
            condition="input.mq!='Default'&&input.evaluationMethod=='Relative'",
            sliderInput('mq_threshold_r','Threshold (percentile)',min=0,max=100,value=5)
        )})
    
    output$bqUI1<-renderUI({
        conditionalPanel(
            condition="input.bq!='Default'&&input.evaluationMethod=='Absolute'",
            numericInput('bq_threshold','Threshold [0;50]',min=0,max=50,value=30)
        )})
    output$bqUI1r<-renderUI({
        conditionalPanel(
            condition="input.bq!='Default'&&input.evaluationMethod=='Relative'",
            sliderInput('bq_threshold_r','Threshold (percentile)',min=0,max=100,value=5)
        )})
    
    output$bnUI1<-renderUI({
        conditionalPanel(
            condition="input.bn!='Default'&&input.evaluationMethod=='Absolute'",
            numericInput('bn_threshold','Threshold [0;100]',min=0,max=100,value=1.0,step = 0.1)
        )})
    output$bnUI1r<-renderUI({
        conditionalPanel(
            condition="input.bn!='Default'&&input.evaluationMethod=='Relative'",
            sliderInput('bn_threshold_r','Threshold (percentile)',min=0,max=100,value=5)
        )})

    output$weightUI1<-renderUI({
        conditionalPanel(
            condition="input.weight!='Default'",
            numericInput("weight1","Coverage",min=0,max=10,value=3)
        )})
    output$weightUI2<-renderUI({
        conditionalPanel(
            condition="input.weight!='Default'",
            numericInput("weight2","Mapping quality",min=0,max=10,value=1)
        )})
    output$weightUI3<-renderUI({
        conditionalPanel(
            condition="input.weight!='Default'",
            numericInput("weight3","Base quality",min=0,max=10,value=1)
        )})
    output$weightUI4<-renderUI({
        conditionalPanel(
            condition="input.weight!='Default'",
            numericInput("weight4","Background noise",min=0,max=10,value=1)
        )})
    
    output$datasetUI1<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset1","Set 1",value=T)
        )})
    output$datasetUI2<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset2","Set 2",value=T)
        )})
    output$datasetUI3<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset3","Set 3",value=T)
        )})
    output$datasetUI4<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset4","Set 4",value=T)
        )})
    output$datasetUI5<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset5","Set 5",value=T)
        )})
    output$datasetUI6<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset6","Set 6",value=T)
        )})
    output$datasetUI7<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset7","Set 7",value=T)
        )})
    output$datasetUI8<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset8","Set 8",value=T)
        )})
    output$datasetUI9<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset9","Set 9",value=T)
        )})
    output$datasetUI10<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset10","Set 10",value=T)
        )})
    output$datasetUI11<-renderUI({
        conditionalPanel(
            condition="input.datasets!='Default'",
            checkboxInput("dataset11","Set 11",value=T)
        )})
    
    output$colorUI1<-renderUI({
        conditionalPanel(
            condition="input.color!='Default'",
            sliderInput("color1","Threshold for marking coverage as problematic [%]",min=0,max=100,value=80)
        )})
    output$colorUI2<-renderUI({
        conditionalPanel(
            condition="input.color!='Default'",
            sliderInput("color2","Threshold for marking quality as problematic [%]",min=0,max=100,value=80)
        )})
    
    output$asxl1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('asxl1','ASXL1',min=0,max=100,value=14.56)
        )})
    output$bcorUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('bcor','BCOR',min=0,max=100,value=2.36)
        )})
    output$cblUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('cbl','CBL',min=0,max=100,value=2.47)
        )})
    output$cebpaUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('cebpa','CEBPA',min=0,max=100,value=1.15)
        )})
    output$dnmt3aUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('dnmt3a','DNMT3A',min=0,max=100,value=9.54)
        )})
    output$etv6UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('etv6','ETV6',min=0,max=100,value=1.16)
        )})
    output$ezh2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('ezh2','EZH2',min=0,max=100,value=5.19)
        )})
    output$flt3UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('flt3','FLT3',min=0,max=100,value=1.94)
        )})
    output$gata2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('gata2','GATA2',min=0,max=100,value=1.16)
        )})
    output$idh1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('idh1','IDH1',min=0,max=100,value=2.50)
        )})
    output$idh2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('idh2','IDH2',min=0,max=100,value=3.98)
        )})
    output$jak2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('jak2','JAK2',min=0,max=100,value=9.31)
        )})
    output$kdm6aUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('kdm6a','KDM6A',min=0,max=100,value=1.58)
        )})
    output$kitUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('kit','KIT',min=0,max=100,value=1.56)
        )})
    output$kmt2aUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('kmt2a','KMT2A',min=0,max=100,value=4.14)
        )})
    output$krasUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('kras','KRAS',min=0,max=100,value=1.42)
        )})
    output$mplUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('mpl','MPL',min=0,max=100,value=1.68)
        )})
    output$npm1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('npm1','NPM1',min=0,max=100,value=1.37)
        )})
    output$nrasUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('nras','NRAS',min=0,max=100,value=5.71)
        )})
    output$rad21UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('rad21','RAD21',min=0,max=100,value=1.10)
        )})
    output$runx1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('runx1','RUNX1',min=0,max=100,value=8.92)
        )})
    output$sf3b1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('sf3b1','SF3B1',min=0,max=100,value=23.42)
        )})
    output$smc1aUI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('smc1a','SMC1A',min=0,max=100,value=1.37)
        )})
    output$smc3UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('smc3','SMC3',min=0,max=100,value=1.37)
        )})
    output$srsf2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('srsf2','SRSF2',min=0,max=100,value=7.43)
        )})
    output$stag2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('stag2','STAG2',min=0,max=100,value=3.91)
        )})
    output$tet2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('tet2','TET2',min=0,max=100,value=18.22)
        )})
    output$tp53UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('tp53','TP53',min=0,max=100,value=8.48)
        )})
    output$u2af1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('u2af1','U2AF1',min=0,max=100,value=7.18)
        )})
    output$wt1UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('wt1','WT1',min=0,max=100,value=0.61)
        )})
    output$zrsr2UI<-renderUI({
        conditionalPanel(
            condition="input.change_impact=='Customize'",
            numericInput('zrsr2','ZRSR2',min=0,max=100,value=5.10)
        )})

    observeEvent(input$do1,{
        #select target
        if(input$targetRegion=="Region of Interest"){
            load(file=("./ofInterest.RData"))
        }
        if(input$targetRegion=="Intersecting Region"){
            load(file=("./intersecting.RData"))
        }
        
        #select data sets
        dataset_vector<-c()
        if(input$dataset1==T){
            dataset_vector<-c(dataset_vector,1)
        }
        if(input$dataset2==T){
            dataset_vector<-c(dataset_vector,2)
        }
        if(input$dataset3==T){
            dataset_vector<-c(dataset_vector,3)
        }
        if(input$dataset4==T){
            dataset_vector<-c(dataset_vector,4)
        }
        if(input$dataset5==T){
            dataset_vector<-c(dataset_vector,5)
        }
        if(input$dataset6==T){
            dataset_vector<-c(dataset_vector,6)
        }
        if(input$dataset7==T){
            dataset_vector<-c(dataset_vector,7)
        }
        if(input$dataset8==T){
            dataset_vector<-c(dataset_vector,8)
        }
        if(input$dataset9==T){
            dataset_vector<-c(dataset_vector,9)
        }
        if(input$dataset10==T){
            dataset_vector<-c(dataset_vector,10)
        }
        if(input$dataset11==T){
            dataset_vector<-c(dataset_vector,11)
        }
        all_samples<-c(111,200,23,22,114,129,46,91,66,13,16)
        
        #coverage
        progress1 <- shiny::Progress$new()
        progress1$set(message = "Updating Coverage", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$cov=="Default"){
                cov_threshold<-100
            }
            if(input$cov!="Default"){
                cov_threshold<-input$cov_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$cov=="Default"){
                cov_threshold<-5
            }
            if(input$cov!="Default"){
                cov_threshold<-input$cov_threshold_r
            }
        }

        for(i in 1:11){
            progress1$inc(1/11)
            binary<-input_data[[1]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[1]][[i]]<cov_threshold]<-0
                binary[input_data[[1]][[i]]>=cov_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[1]][[i]]),cov_threshold/100,na.rm=T)
                binary[input_data[[1]][[i]]<grenze]<-0
                binary[input_data[[1]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[1]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm=T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        target_cov<-input_data[[5]]
        genes<-unique(target_cov[,4])
        info<-data.frame(genes,start=NA,end=NA,length=NA,stringsAsFactors=F)
        for(i in 1:length(genes)){
            relevant<-target_cov[target_cov[,4]==genes[i],]
            if(i==1){
                start<-1
                end<-sum(relevant[,3]-relevant[,2]+1)
            }
            if(i>1){
                start<-end+1
                end<-start+sum(relevant[,3]-relevant[,2]+1)-1
            }
            info[i,2]<-start
            info[i,3]<-end
            info[i,4]<-end-start+1
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)

        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }

        score_cov<-coverage_quality
        overview<-data.frame(Gene=info[,1],Length=info[,4],Coverage_Score=1-rowSums(score_cov[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress1$close()
        
        #MQ
        progress2 <- shiny::Progress$new()
        progress2$set(message = "Updating Mapping Quality", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$mq=="Default"){
                mq_threshold<-58
            }
            if(input$mq!="Default"){
                mq_threshold<-input$mq_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$mq=="Default"){
                mq_threshold<-5
            }
            if(input$mq!="Default"){
                mq_threshold<-input$mq_threshold_r
            }
        }
        for(i in 1:11){
            progress2$inc(1/11)
            binary<-input_data[[4]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[4]][[i]]<mq_threshold]<-0
                binary[input_data[[4]][[i]]>=mq_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[4]][[i]]),mq_threshold/100,na.rm=T)
                binary[input_data[[4]][[i]]<grenze]<-0
                binary[input_data[[4]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[4]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }

        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }

        score_mq<-coverage_quality
        overview<-cbind(overview,MQ_Score=1-rowSums(score_mq[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress2$close()
        
        
        #BQ
        progress3 <- shiny::Progress$new()
        progress3$set(message = "Updating Base Quality", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$bq=="Default"){
                bq_threshold<-30
            }
            if(input$bq!="Default"){
                bq_threshold<-input$bq_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$bq=="Default"){
                bq_threshold<-5
            }
            if(input$bq!="Default"){
                bq_threshold<-input$bq_threshold_r
            }
        }
        for(i in 1:11){
            progress3$inc(1/11)
            binary<-input_data[[3]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[3]][[i]]<bq_threshold]<-0
                binary[input_data[[3]][[i]]>=bq_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[3]][[i]]),bq_threshold/100,na.rm=T)
                binary[input_data[[3]][[i]]<grenze]<-0
                binary[input_data[[3]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[3]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }

        score_bq<-coverage_quality
        overview<-cbind(overview,BQ_Score=1-rowSums(score_bq[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress3$close()
        
        #Noise
        progress4 <- shiny::Progress$new()
        progress4$set(message = "Updating Background Noise", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$bn=="Default"){
                bn_threshold<-1
            }
            if(input$bn!="Default"){
                bn_threshold<-input$bn_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$bn=="Default"){
                bn_threshold<-5
            }
            if(input$bn!="Default"){
                bn_threshold<-input$bn_threshold_r
            }
        }
        for(i in 1:11){
            progress4$inc(1/11)
            binary<-input_data[[2]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[2]][[i]]<(100-bn_threshold)/100]<-0
                binary[input_data[[2]][[i]]>=(100-bn_threshold)/100]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[2]][[i]]),bn_threshold/100,na.rm=T)
                binary[input_data[[2]][[i]]<grenze]<-0
                binary[input_data[[2]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[2]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }
 
        score_bn<-coverage_quality
        overview<-cbind(overview,Noise_Score=1-rowSums(score_bn[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress4$close()
        
        if(input$targetRegion=="Intersecting Region"){
            overview<-overview[order(overview[,1]),]
        }
        
        score_bn<-as.numeric(overview[,6])
        score_bq<-as.numeric(overview[,5])
        score_mq<-as.numeric(overview[,4])
        score_cov<-as.numeric(overview[,3])
        
        if(input$weight=="Default"){
            qual_score<-rowSums(cbind(score_bn,score_bq,score_mq),na.rm = T)/3
            mean_difficulty<-(qual_score+score_cov)/2
            
            seq_simp<-mean_difficulty
            qual_simp<-qual_score
            cov_simp<-score_cov
        }
        if(input$weight!="Default"){
            qual_score<-rowSums(cbind(input$weight4*score_bn,input$weight3*score_bq,input$weight2*score_mq),na.rm = T)/(input$weight2+input$weight3+input$weight4)
            mean_difficulty<-(input$weight1*qual_score+(input$weight2+input$weight3+input$weight4)*score_cov)/(input$weight1+input$weight2+input$weight3+input$weight4)
            
            seq_simp<-mean_difficulty
            qual_simp<-qual_score
            cov_simp<-score_cov
        }

        overview<-cbind(overview,Qual_Score=round(qual_score,digits = 4),Seq_Score=round(seq_simp,digits = 4))
        
        clin_impact<-c()
        if(input$change_impact=="Cosmic"){
            clin_impact[1]<-14.56
            clin_impact[2]<-2.36
            clin_impact[3]<-2.47
            clin_impact[4]<-1.15
            clin_impact[5]<-9.54
            clin_impact[6]<-1.16
            clin_impact[7]<-5.19
            clin_impact[8]<-1.94
            clin_impact[9]<-1.16
            clin_impact[10]<-2.5
            clin_impact[11]<-3.98
            clin_impact[12]<-9.31
            clin_impact[13]<-1.58
            clin_impact[14]<-1.56
            clin_impact[15]<-4.14
            clin_impact[16]<-1.42
            clin_impact[17]<-1.68
            clin_impact[18]<-1.37
            clin_impact[19]<-5.71
            clin_impact[20]<-1.1
            clin_impact[21]<-8.92
            clin_impact[22]<-23.42
            clin_impact[23]<-1.37
            clin_impact[24]<-1.37
            clin_impact[25]<-7.43
            clin_impact[26]<-3.91
            clin_impact[27]<-18.22
            clin_impact[28]<-8.48
            clin_impact[29]<-7.18
            clin_impact[30]<-0.61
            clin_impact[31]<-5.1
        }
        if(input$change_impact=="Haferlach et al. 2014"){
            clin_impact[1]<-23.41
            clin_impact[2]<-4.03
            clin_impact[3]<-5.08
            clin_impact[4]<-0.95
            clin_impact[5]<-13.14
            clin_impact[6]<-2.33
            clin_impact[7]<-5.51
            clin_impact[8]<-1.17
            clin_impact[9]<-0.74
            clin_impact[10]<-2.54
            clin_impact[11]<-3.92
            clin_impact[12]<-4.77
            clin_impact[13]<-0.00
            clin_impact[14]<-0.85
            clin_impact[15]<-0.00
            clin_impact[16]<-2.54
            clin_impact[17]<-2.97
            clin_impact[18]<-0.95
            clin_impact[19]<-3.81
            clin_impact[20]<-1.38
            clin_impact[21]<-10.59
            clin_impact[22]<-32.94
            clin_impact[23]<-1.06
            clin_impact[24]<-1.69
            clin_impact[25]<-17.48
            clin_impact[26]<-7.52
            clin_impact[27]<-33.26
            clin_impact[28]<-6.36
            clin_impact[29]<-7.73
            clin_impact[30]<-0.53
            clin_impact[31]<-7.63
        }
        if(input$change_impact=="Papaemmanuil et al. 2013"){
            clin_impact[1]<-14.26
            clin_impact[2]<-2.91
            clin_impact[3]<-3.68
            clin_impact[4]<-0.31
            clin_impact[5]<-10.58
            clin_impact[6]<-0.46
            clin_impact[7]<-4.91
            clin_impact[8]<-0.31
            clin_impact[9]<-1.07
            clin_impact[10]<-2.15
            clin_impact[11]<-3.99
            clin_impact[12]<-2.76
            clin_impact[13]<-0.46
            clin_impact[14]<-0.77
            clin_impact[15]<-0
            clin_impact[16]<-1.38
            clin_impact[17]<-0.77
            clin_impact[18]<-1.23
            clin_impact[19]<-3.07
            clin_impact[20]<-0.61
            clin_impact[21]<-8.9
            clin_impact[22]<-27.3
            clin_impact[23]<-0
            clin_impact[24]<-0
            clin_impact[25]<-16.72
            clin_impact[26]<-4.14
            clin_impact[27]<-26.99
            clin_impact[28]<-5.37
            clin_impact[29]<-6.75
            clin_impact[30]<-0.77
            clin_impact[31]<-3.99
        }
        if(input$change_impact=="Customize"){
            clin_impact[1]<-input$asxl1
            clin_impact[2]<-input$bcor
            clin_impact[3]<-input$cbl
            clin_impact[4]<-input$cebpa
            clin_impact[5]<-input$dnmt3a
            clin_impact[6]<-input$etv6
            clin_impact[7]<-input$ezh2
            clin_impact[8]<-input$flt3
            clin_impact[9]<-input$gata2
            clin_impact[10]<-input$idh1
            clin_impact[11]<-input$idh2
            clin_impact[12]<-input$jak2
            clin_impact[13]<-input$kdm6a
            clin_impact[14]<-input$kit
            clin_impact[15]<-input$kmt2a
            clin_impact[16]<-input$kras
            clin_impact[17]<-input$mpl
            clin_impact[18]<-input$npm1
            clin_impact[19]<-input$nras
            clin_impact[20]<-input$rad21
            clin_impact[21]<-input$runx1
            clin_impact[22]<-input$sf3b1
            clin_impact[23]<-input$smc1a
            clin_impact[24]<-input$smc3
            clin_impact[25]<-input$srsf2
            clin_impact[26]<-input$stag2
            clin_impact[27]<-input$tet2
            clin_impact[28]<-input$tp53
            clin_impact[29]<-input$u2af1
            clin_impact[30]<-input$wt1
            clin_impact[31]<-input$zrsr2
        }
        
        if(input$targetRegion=="Intersecting Region"){
            clin_impact<-clin_impact[c(1:24,26:31)]
            overview<-overview[order(overview[,1]),]
        }
        overview<-cbind(overview,Mutation_Frequency=clin_impact)
        
        output$table <- renderDataTable(datatable(overview))
        output$table1 <- renderDataTable(datatable(overview))
        
        results<-cbind(overview[,1:8],cov_simp,qual_simp,overview[,9])
        results[,9]<-as.numeric(results[,9])
        results[,10]<-as.numeric(results[,10])
        
        if(input$change_impact!="Cosmic"){
            result_backup<-results
            results<-results[results[,11]!=0,]
        }

        color_inside<-c()
        color_outside<-c()
        for(i in 1:length(results[,1])){
            if(input$color=="Default"){
                #message(results[i,1]," i=",i," cov=",results[i,9]," qual=",results[i,10]," insg=",results[i,8])
                if(results[i,9]<0.8){
                    color_inside[i]<-"red"
                }
                if(results[i,9]>=0.8){
                    color_inside[i]<-"black"
                }
                if(results[i,10]<0.8){
                    color_outside[i]<-"lightblue"
                }
                if(results[i,10]>=0.8){
                    color_outside[i]<-"black"
                }
            }
            if(input$color!="Default"){
                if(results[i,9]<(input$color1/100)){
                    color_inside[i]<-"red"
                }
                if(results[i,9]>=(input$color1/100)){
                    color_inside[i]<-"black"
                }
                if(results[i,10]<(input$color2/100)){
                    color_outside[i]<-"lightblue"
                }
                if(results[i,10]>=(input$color2/100)){
                    color_outside[i]<-"black"
                }
            }
        }
        
        if(input$change_impact=="Cosmic"){
            output$plot<-output$plot1<-renderPlot({
                par(mar=c(4,4,0.5,0.5))
                plot(results[,8],results[,11],xlab="Expected rate of bases with minimum requirements [%]",
                     ylab="Mutation frequency [%]",xaxt="n",
                     xlim=c((round(min(result_backup[,8]),digits = 2)-0.05),(round(max(result_backup[,8]),digits = 2)+0.05)),
                     ylim=c(0,min(100,ceiling(max(results[,11])))),
                     cex.lab=1.5,pch=21,col=color_outside,bg=color_inside,lwd=3.5,cex=2)
                legend(x=max(result_backup[,8])-0.32,y=max(result_backup[,11])-0.1,
                       legend = c("Problematic coverage and quality","Problematic coverage",
                                  "Problematic quality","Unproblematic coverage and quality"),
                       pch=21,bty="n",pt.bg=c("red","red","black","black"),pt.cex = 1.7,pt.lwd = 3.5,
                       col = c("lightblue","black","lightblue","black"),cex=1.3)
                axis(1,at=seq(0,1,0.1),labels=seq(0,100,10),cex.axis=1.3)
                if(input$evaluationMethod=="Absolute"&&input$targetRegion=="Region of Interest"){
                    text(results[,8],results[,11],results[,1],
                         pos=c(3,3,2,1,3,1,3,3,3,3,
                               3,3,2,1,3,3,3,3,3,3,
                               3,3,1,3,3,3,3,3,3,1,
                               3),cex=1.2,offset = 0.7)
                }
                if(input$evaluationMethod=="Relative"&&input$targetRegion=="Intersecting Region"){
                    text(results[,8],results[,11],results[,1],
                         pos=c(3,3,2,3,3,1,3,2,3,4,
                               4,4,3,4,3,3,3,3,3,1,
                               3,3,4,2,3,3,3,3,1,3),cex=1.2,offset = 0.7)
                }
                if((input$evaluationMethod=="Relative"&&input$targetRegion=="Region of Interest")||
                   (input$evaluationMethod=="Absolute"&&input$targetRegion=="Intersecting Region")){
                    text(results[,8],results[,11],results[,1],
                         pos=c(3,1,3,1,3,1,3,1,3,1,
                               3,1,3,1,3,1,3,1,3,1,
                               3,1,3,1,3,1,3,1,3,1),cex=1.2,offset = 0.7)
                }
            },width=1000,height=1000)
        }
        if(input$change_impact!="Cosmic"){
            if(input$change_impact=="Haferlach et al. 2014"||input$change_impact=="Papaemmanuil et al. 2013"){
                max_cat<-max(results[,11])
                min_cat<-min(results[,11])
            }
            if(input$change_impact=="Customize"){
                max_cat<-max(input$asxl1,input$bcor,input$cbl,input$cebpa,input$dnmt3a,
                             input$etv6,input$ezh2,input$flt3,input$gata2,input$idh1,
                             input$idh2,input$jak2,input$kdm6a,input$kit,input$kmt2a,
                             input$kras,input$mpl,input$npm1,input$nras,input$rad21,
                             input$runx1,input$sf3b1,input$smc1a,input$smc3,input$srsf2,
                             input$stag2,input$tet2,input$tp53,input$u2af1,input$wt1,input$zrsr2)
                min_cat<-min(input$asxl1,input$bcor,input$cbl,input$cebpa,input$dnmt3a,
                             input$etv6,input$ezh2,input$flt3,input$gata2,input$idh1,
                             input$idh2,input$jak2,input$kdm6a,input$kit,input$kmt2a,
                             input$kras,input$mpl,input$npm1,input$nras,input$rad21,
                             input$runx1,input$sf3b1,input$smc1a,input$smc3,input$srsf2,
                             input$stag2,input$tet2,input$tp53,input$u2af1,input$wt1,input$zrsr2)
            }
            results2<-results[results[,11]!=0,]
            
            output$plot<-output$plot1<-renderPlot({
                par(mar=c(4,4,0.5,0.5))
                plot(results2[,8],results2[,11],xlab="Expected rate of unproblematic bases",
                     ylab="Mutation frequency [%]",xaxt="n",
                     xlim=c((round(min(result_backup[,8]),digits = 2)-0.05),(round(max(result_backup[,8]),digits = 2)+0.05)),
                     ylim=c(0,(max_cat+0.5)),
                     cex.lab=1.5,pch=21,col=color_outside,bg=color_inside,lwd=3.5,cex=2)
                legend(x=max(results2[,8])-0.15,y=max(results2[,11]),
                       legend = c("Problematic coverage and quality","Problematic coverage",
                                  "Problematic quality","Unproblematic coverage and quality"),
                       pch=21,bty="n",pt.bg=c("red","red","black","black"),pt.cex = 1.7,pt.lwd = 3.5,
                       col = c("lightblue","black","lightblue","black"),cex=1.3)
                
                axis(1,at=seq(0,1,0.1),labels=seq(0,100,10),cex.axis=1.3)
                
                if(input$change_impact=="Haferlach et al. 2014"&&input$evaluationMethod=="Absolute"&&input$targetRegion=="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,3,2,3,3,1,3,4,3,4,
                               3,3,1,1,3,3,4,1,2,3,
                               3,3,2,2,3,1,3,3,3,1,
                               3),cex=1,offset = 0.7)
                }
                if(input$change_impact=="Papaemmanuil et al. 2013"&&input$evaluationMethod=="Absolute"&&input$targetRegion=="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,3,3,3,3,1,3,4,3,3,
                               3,3,2,2,3,3,4,3,3,1,
                               3,3,2,1,3,2,3,3,3,3,
                               3),cex=1,offset = 0.7)
                }
                if((input$change_impact!="Haferlach et al. 2014"&&input$change_impact!="Papaemmanuil et al. 2013")||
                   input$evaluationMethod!="Absolute"||input$targetRegion!="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,1),cex=1,offset = 0.7)
                }

            },width=1000,height=1000)
        }
    })
  
    observeEvent(input$do2,{
        #select target
        if(input$targetRegion=="Region of Interest"){
            load(file=("./ofInterest.RData"))
        }
        if(input$targetRegion=="Intersecting Region"){
            load(file=("./intersecting.RData"))
        }
        
        #select data sets
        dataset_vector<-c()
        if(input$dataset1==T){
            dataset_vector<-c(dataset_vector,1)
        }
        if(input$dataset2==T){
            dataset_vector<-c(dataset_vector,2)
        }
        if(input$dataset3==T){
            dataset_vector<-c(dataset_vector,3)
        }
        if(input$dataset4==T){
            dataset_vector<-c(dataset_vector,4)
        }
        if(input$dataset5==T){
            dataset_vector<-c(dataset_vector,5)
        }
        if(input$dataset6==T){
            dataset_vector<-c(dataset_vector,6)
        }
        if(input$dataset7==T){
            dataset_vector<-c(dataset_vector,7)
        }
        if(input$dataset8==T){
            dataset_vector<-c(dataset_vector,8)
        }
        if(input$dataset9==T){
            dataset_vector<-c(dataset_vector,9)
        }
        if(input$dataset10==T){
            dataset_vector<-c(dataset_vector,10)
        }
        if(input$dataset11==T){
            dataset_vector<-c(dataset_vector,11)
        }
        all_samples<-c(111,200,23,22,114,129,46,91,66,13,16)
        
        #coverage
        progress1 <- shiny::Progress$new()
        progress1$set(message = "Updating Coverage", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$cov=="Default"){
                cov_threshold<-100
            }
            if(input$cov!="Default"){
                cov_threshold<-input$cov_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$cov=="Default"){
                cov_threshold<-5
            }
            if(input$cov!="Default"){
                cov_threshold<-input$cov_threshold_r
            }
        }
        
        for(i in 1:11){
            progress1$inc(1/11)
            binary<-input_data[[1]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[1]][[i]]<cov_threshold]<-0
                binary[input_data[[1]][[i]]>=cov_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[1]][[i]]),cov_threshold/100,na.rm=T)
                binary[input_data[[1]][[i]]<grenze]<-0
                binary[input_data[[1]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[1]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm=T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        target_cov<-input_data[[5]]
        genes<-unique(target_cov[,4])
        info<-data.frame(genes,start=NA,end=NA,length=NA,stringsAsFactors=F)
        for(i in 1:length(genes)){
            relevant<-target_cov[target_cov[,4]==genes[i],]
            if(i==1){
                start<-1
                end<-sum(relevant[,3]-relevant[,2]+1)
            }
            if(i>1){
                start<-end+1
                end<-start+sum(relevant[,3]-relevant[,2]+1)-1
            }
            info[i,2]<-start
            info[i,3]<-end
            info[i,4]<-end-start+1
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }
        
        score_cov<-coverage_quality
        overview<-data.frame(Gene=info[,1],Length=info[,4],Coverage_Score=1-rowSums(score_cov[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress1$close()
        
        #MQ
        progress2 <- shiny::Progress$new()
        progress2$set(message = "Updating Mapping Quality", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$mq=="Default"){
                mq_threshold<-58
            }
            if(input$mq!="Default"){
                mq_threshold<-input$mq_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$mq=="Default"){
                mq_threshold<-5
            }
            if(input$mq!="Default"){
                mq_threshold<-input$mq_threshold_r
            }
        }
        for(i in 1:11){
            progress2$inc(1/11)
            binary<-input_data[[4]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[4]][[i]]<mq_threshold]<-0
                binary[input_data[[4]][[i]]>=mq_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[4]][[i]]),mq_threshold/100,na.rm=T)
                binary[input_data[[4]][[i]]<grenze]<-0
                binary[input_data[[4]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[4]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }
        
        score_mq<-coverage_quality
        overview<-cbind(overview,MQ_Score=1-rowSums(score_mq[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress2$close()
        
        
        #BQ
        progress3 <- shiny::Progress$new()
        progress3$set(message = "Updating Base Quality", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$bq=="Default"){
                bq_threshold<-30
            }
            if(input$bq!="Default"){
                bq_threshold<-input$bq_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$bq=="Default"){
                bq_threshold<-5
            }
            if(input$bq!="Default"){
                bq_threshold<-input$bq_threshold_r
            }
        }
        for(i in 1:11){
            progress3$inc(1/11)
            binary<-input_data[[3]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[3]][[i]]<bq_threshold]<-0
                binary[input_data[[3]][[i]]>=bq_threshold]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[3]][[i]]),bq_threshold/100,na.rm=T)
                binary[input_data[[3]][[i]]<grenze]<-0
                binary[input_data[[3]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[3]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }
        
        score_bq<-coverage_quality
        overview<-cbind(overview,BQ_Score=1-rowSums(score_bq[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress3$close()
        
        #Noise
        progress4 <- shiny::Progress$new()
        progress4$set(message = "Updating Background Noise", value = 0)
        all_means<-data.frame()
        if(input$evaluationMethod=="Absolute"){
            if(input$bn=="Default"){
                bn_threshold<-1
            }
            if(input$bn!="Default"){
                bn_threshold<-input$bn_threshold
            }
        }
        if(input$evaluationMethod=="Relative"){
            if(input$bn=="Default"){
                bn_threshold<-5
            }
            if(input$bn!="Default"){
                bn_threshold<-input$bn_threshold_r
            }
        }
        for(i in 1:11){
            progress4$inc(1/11)
            binary<-input_data[[2]][[i]]
            if(input$evaluationMethod=="Absolute"){
                binary[input_data[[2]][[i]]<(100-bn_threshold)/100]<-0
                binary[input_data[[2]][[i]]>=(100-bn_threshold)/100]<-1
            }
            if(input$evaluationMethod=="Relative"){
                grenze<-quantile(unlist(input_data[[2]][[i]]),bn_threshold/100,na.rm=T)
                binary[input_data[[2]][[i]]<grenze]<-0
                binary[input_data[[2]][[i]]>=grenze]<-1
            }
            binary[is.na(input_data[[2]][[i]])]<-0
            binary_sum<-rowSums(binary,na.rm = T)
            if(i==1){
                all_means<-binary_sum
            }
            if(i>1){
                all_means<-rbind(all_means,binary_sum)
            }
        }
        
        coverage_quality<-cbind(info,V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,
                                V8=NA,V9=NA,V10=NA,V11=NA)
        
        for(i in 1:length(coverage_quality[,1])){
            coverage_quality[i,5]<-111-sum(all_means[1,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,6]<-200-sum(all_means[2,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,7]<-23-sum(all_means[3,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,8]<-22-sum(all_means[4,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,9]<-114-sum(all_means[5,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,10]<-129-sum(all_means[6,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,11]<-46-sum(all_means[7,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,12]<-91-sum(all_means[8,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,13]<-66-sum(all_means[9,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,14]<-13-sum(all_means[10,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
            coverage_quality[i,15]<-16-sum(all_means[11,coverage_quality[i,2]:coverage_quality[i,3]],na.rm=T)/(coverage_quality[i,4])
        }
        
        score_bn<-coverage_quality
        overview<-cbind(overview,Noise_Score=1-rowSums(score_bn[,(4+dataset_vector)],na.rm = T)/sum(all_samples[dataset_vector]))
        progress4$close()
        
        if(input$targetRegion=="Intersecting Region"){
            overview<-overview[order(overview[,1]),]
        }
        
        score_bn<-as.numeric(overview[,6])
        score_bq<-as.numeric(overview[,5])
        score_mq<-as.numeric(overview[,4])
        score_cov<-as.numeric(overview[,3])
        
        if(input$weight=="Default"){
            qual_score<-rowSums(cbind(score_bn,score_bq,score_mq),na.rm = T)/3
            mean_difficulty<-(qual_score+score_cov)/2
            
            seq_simp<-mean_difficulty
            qual_simp<-qual_score
            cov_simp<-score_cov
        }
        if(input$weight!="Default"){
            qual_score<-rowSums(cbind(input$weight4*score_bn,input$weight3*score_bq,input$weight2*score_mq),na.rm = T)/(input$weight2+input$weight3+input$weight4)
            mean_difficulty<-(input$weight1*qual_score+(input$weight2+input$weight3+input$weight4)*score_cov)/(input$weight1+input$weight2+input$weight3+input$weight4)
            
            seq_simp<-mean_difficulty
            qual_simp<-qual_score
            cov_simp<-score_cov
        }
        
        overview<-cbind(overview,Qual_Score=round(qual_score,digits = 4),Seq_Score=round(seq_simp,digits = 4))
        
        clin_impact<-c()
        if(input$change_impact=="Cosmic"){
            clin_impact[1]<-14.56
            clin_impact[2]<-2.36
            clin_impact[3]<-2.47
            clin_impact[4]<-1.15
            clin_impact[5]<-9.54
            clin_impact[6]<-1.16
            clin_impact[7]<-5.19
            clin_impact[8]<-1.94
            clin_impact[9]<-1.16
            clin_impact[10]<-2.5
            clin_impact[11]<-3.98
            clin_impact[12]<-9.31
            clin_impact[13]<-1.58
            clin_impact[14]<-1.56
            clin_impact[15]<-4.14
            clin_impact[16]<-1.42
            clin_impact[17]<-1.68
            clin_impact[18]<-1.37
            clin_impact[19]<-5.71
            clin_impact[20]<-1.1
            clin_impact[21]<-8.92
            clin_impact[22]<-23.42
            clin_impact[23]<-1.37
            clin_impact[24]<-1.37
            clin_impact[25]<-7.43
            clin_impact[26]<-3.91
            clin_impact[27]<-18.22
            clin_impact[28]<-8.48
            clin_impact[29]<-7.18
            clin_impact[30]<-0.61
            clin_impact[31]<-5.1
        }
        if(input$change_impact=="Haferlach et al. 2014"){
            clin_impact[1]<-23.41
            clin_impact[2]<-4.03
            clin_impact[3]<-5.08
            clin_impact[4]<-0.95
            clin_impact[5]<-13.14
            clin_impact[6]<-2.33
            clin_impact[7]<-5.51
            clin_impact[8]<-1.17
            clin_impact[9]<-0.74
            clin_impact[10]<-2.54
            clin_impact[11]<-3.92
            clin_impact[12]<-4.77
            clin_impact[13]<-0.00
            clin_impact[14]<-0.85
            clin_impact[15]<-0.00
            clin_impact[16]<-2.54
            clin_impact[17]<-2.97
            clin_impact[18]<-0.95
            clin_impact[19]<-3.81
            clin_impact[20]<-1.38
            clin_impact[21]<-10.59
            clin_impact[22]<-32.94
            clin_impact[23]<-1.06
            clin_impact[24]<-1.69
            clin_impact[25]<-17.48
            clin_impact[26]<-7.52
            clin_impact[27]<-33.26
            clin_impact[28]<-6.36
            clin_impact[29]<-7.73
            clin_impact[30]<-0.53
            clin_impact[31]<-7.63
        }
        if(input$change_impact=="Papaemmanuil et al. 2013"){
            clin_impact[1]<-14.26
            clin_impact[2]<-2.91
            clin_impact[3]<-3.68
            clin_impact[4]<-0.31
            clin_impact[5]<-10.58
            clin_impact[6]<-0.46
            clin_impact[7]<-4.91
            clin_impact[8]<-0.31
            clin_impact[9]<-1.07
            clin_impact[10]<-2.15
            clin_impact[11]<-3.99
            clin_impact[12]<-2.76
            clin_impact[13]<-0.46
            clin_impact[14]<-0.77
            clin_impact[15]<-0
            clin_impact[16]<-1.38
            clin_impact[17]<-0.77
            clin_impact[18]<-1.23
            clin_impact[19]<-3.07
            clin_impact[20]<-0.61
            clin_impact[21]<-8.9
            clin_impact[22]<-27.3
            clin_impact[23]<-0
            clin_impact[24]<-0
            clin_impact[25]<-16.72
            clin_impact[26]<-4.14
            clin_impact[27]<-26.99
            clin_impact[28]<-5.37
            clin_impact[29]<-6.75
            clin_impact[30]<-0.77
            clin_impact[31]<-3.99
        }
        if(input$change_impact=="Customize"){
            clin_impact[1]<-input$asxl1
            clin_impact[2]<-input$bcor
            clin_impact[3]<-input$cbl
            clin_impact[4]<-input$cebpa
            clin_impact[5]<-input$dnmt3a
            clin_impact[6]<-input$etv6
            clin_impact[7]<-input$ezh2
            clin_impact[8]<-input$flt3
            clin_impact[9]<-input$gata2
            clin_impact[10]<-input$idh1
            clin_impact[11]<-input$idh2
            clin_impact[12]<-input$jak2
            clin_impact[13]<-input$kdm6a
            clin_impact[14]<-input$kit
            clin_impact[15]<-input$kmt2a
            clin_impact[16]<-input$kras
            clin_impact[17]<-input$mpl
            clin_impact[18]<-input$npm1
            clin_impact[19]<-input$nras
            clin_impact[20]<-input$rad21
            clin_impact[21]<-input$runx1
            clin_impact[22]<-input$sf3b1
            clin_impact[23]<-input$smc1a
            clin_impact[24]<-input$smc3
            clin_impact[25]<-input$srsf2
            clin_impact[26]<-input$stag2
            clin_impact[27]<-input$tet2
            clin_impact[28]<-input$tp53
            clin_impact[29]<-input$u2af1
            clin_impact[30]<-input$wt1
            clin_impact[31]<-input$zrsr2
        }
        
        if(input$targetRegion=="Intersecting Region"){
            clin_impact<-clin_impact[c(1:24,26:31)]
            overview<-overview[order(overview[,1]),]
        }
        
        overview<-cbind(overview,Mutation_Frequency=clin_impact)
        
        output$table <- renderDataTable(datatable(overview))
        output$table1 <- renderDataTable(datatable(overview))
        
        results<-cbind(overview[,1:8],cov_simp,qual_simp,overview[,9])
        results[,9]<-as.numeric(results[,9])
        results[,10]<-as.numeric(results[,10])
        
        if(input$change_impact!="Cosmic"){
            result_backup<-results
            results<-results[results[,11]!=0,]
        }
        
        color_inside<-c()
        color_outside<-c()
        for(i in 1:length(results[,1])){
            if(input$color=="Default"){
                #message(results[i,1]," i=",i," cov=",results[i,9]," qual=",results[i,10]," insg=",results[i,8])
                if(results[i,9]<0.8){
                    color_inside[i]<-"red"
                }
                if(results[i,9]>=0.8){
                    color_inside[i]<-"black"
                }
                if(results[i,10]<0.8){
                    color_outside[i]<-"lightblue"
                }
                if(results[i,10]>=0.8){
                    color_outside[i]<-"black"
                }
            }
            if(input$color!="Default"){
                if(results[i,9]<(input$color1/100)){
                    color_inside[i]<-"red"
                }
                if(results[i,9]>=(input$color1/100)){
                    color_inside[i]<-"black"
                }
                if(results[i,10]<(input$color2/100)){
                    color_outside[i]<-"lightblue"
                }
                if(results[i,10]>=(input$color2/100)){
                    color_outside[i]<-"black"
                }
            }
        }
        
        if(input$change_impact=="Cosmic"){
            output$plot<-output$plot1<-renderPlot({
                par(mar=c(4,4,0.5,0.5))
                plot(results[,8],results[,11],xlab="Expected rate of bases with minimum requirements [%]",
                     ylab="Mutation frequency [%]",xaxt="n",
                     xlim=c((round(min(result_backup[,8]),digits = 2)-0.05),(round(max(result_backup[,8]),digits = 2)+0.05)),
                     ylim=c(0,min(100,ceiling(max(results[,11])))),
                     cex.lab=1.5,pch=21,col=color_outside,bg=color_inside,lwd=3.5,cex=2)
                legend(x=max(result_backup[,8])-0.32,y=max(result_backup[,11])-0.1,
                       legend = c("Problematic coverage and quality","Problematic coverage",
                                  "Problematic quality","Unproblematic coverage and quality"),
                       pch=21,bty="n",pt.bg=c("red","red","black","black"),pt.cex = 1.7,pt.lwd = 3.5,
                       col = c("lightblue","black","lightblue","black"),cex=1.3)
                axis(1,at=seq(0,1,0.1),labels=seq(0,100,10),cex.axis=1.3)
                text(results[,8],results[,11],results[,1],
                     pos=c(3,3,2,1,3,1,3,3,3,3,
                           3,3,2,1,3,3,3,3,3,3,
                           3,3,1,3,3,3,3,3,3,1,
                           3),cex=1.2,offset = 0.7)
            },width=1000,height=1000)
        }
        if(input$change_impact!="Cosmic"){
            if(input$change_impact=="Haferlach et al. 2014"||input$change_impact=="Papaemmanuil et al. 2013"){
                max_cat<-max(results[,11])
                min_cat<-min(results[,11])
            }
            if(input$change_impact=="Customize"){
                max_cat<-max(input$asxl1,input$bcor,input$cbl,input$cebpa,input$dnmt3a,
                             input$etv6,input$ezh2,input$flt3,input$gata2,input$idh1,
                             input$idh2,input$jak2,input$kdm6a,input$kit,input$kmt2a,
                             input$kras,input$mpl,input$npm1,input$nras,input$rad21,
                             input$runx1,input$sf3b1,input$smc1a,input$smc3,input$srsf2,
                             input$stag2,input$tet2,input$tp53,input$u2af1,input$wt1,input$zrsr2)
                min_cat<-min(input$asxl1,input$bcor,input$cbl,input$cebpa,input$dnmt3a,
                             input$etv6,input$ezh2,input$flt3,input$gata2,input$idh1,
                             input$idh2,input$jak2,input$kdm6a,input$kit,input$kmt2a,
                             input$kras,input$mpl,input$npm1,input$nras,input$rad21,
                             input$runx1,input$sf3b1,input$smc1a,input$smc3,input$srsf2,
                             input$stag2,input$tet2,input$tp53,input$u2af1,input$wt1,input$zrsr2)
            }
            results2<-results[results[,11]!=0,]
            
            output$plot<-output$plot1<-renderPlot({
                par(mar=c(4,4,0.5,0.5))
                plot(results2[,8],results2[,11],xlab="Expected rate of unproblematic bases",
                     ylab="Mutation frequency [%]",xaxt="n",
                     xlim=c((round(min(result_backup[,8]),digits = 2)-0.05),(round(max(result_backup[,8]),digits = 2)+0.05)),
                     ylim=c(0,(max_cat+0.5)),
                     cex.lab=1.5,pch=21,col=color_outside,bg=color_inside,lwd=3.5,cex=2)
                legend(x=max(results2[,8])-0.15,y=max(results2[,11]),
                       legend = c("Problematic coverage and quality","Problematic coverage",
                                  "Problematic quality","Unproblematic coverage and quality"),
                       pch=21,bty="n",pt.bg=c("red","red","black","black"),pt.cex = 1.7,pt.lwd = 3.5,
                       col = c("lightblue","black","lightblue","black"),cex=1.3)
                
                axis(1,at=seq(0,1,0.1),labels=seq(0,100,10),cex.axis=1.3)
                
                if(input$change_impact=="Haferlach et al. 2014"&&input$evaluationMethod=="Absolute"&&input$targetRegion=="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,3,2,3,3,1,3,4,3,4,
                               3,3,1,1,3,3,4,1,2,3,
                               3,3,2,2,3,1,3,3,3,1,
                               3),cex=1,offset = 0.7)
                }
                if(input$change_impact=="Papaemmanuil et al. 2013"&&input$evaluationMethod=="Absolute"&&input$targetRegion=="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,3,3,3,3,1,3,4,3,3,
                               3,3,2,2,3,3,4,3,3,1,
                               3,3,2,1,3,2,3,3,3,3,
                               3),cex=1,offset = 0.7)
                }
                if((input$change_impact!="Haferlach et al. 2014"&&input$change_impact!="Papaemmanuil et al. 2013")||
                   input$evaluationMethod!="Absolute"||input$targetRegion!="Region of Interest"){
                    text(results2[,8],results2[,11],results2[,1],
                         pos=c(3,1),cex=1,offset = 0.7)
                }
                
            },width=1000,height=1000)
        }
    })
})




