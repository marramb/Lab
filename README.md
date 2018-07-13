# Analysis of data generated in a qPCR assay

From the qPCR software, an Excel file carrying the run information can be exported. To start the analysis it is important to modify the file in Excel to the following structure:

![structure_data](https://github.com/jmuribes/images/blob/master/structure_data.png)

Just as a reminder, the **Ct** value corresponds to the threshold cycle, the number of cycles it took to detect a real signal from your samples ([source](https://bitesizebio.com/24581/what-is-a-ct-value/)).

![Ct image](https://bitesizebio.com/wp-content/uploads/2015/06/PCR.png)

Now, the qPCR data analysis is based on determining the average Ct value per primer pair, calculating the difference between that average Ct value and the Ct value (dCt) from a reference primer pair from a housekeeping gene, and also if there is a difference between conditions (i.ex. het/homozygous against the wild-type).

To start, we have to save the Excel sheet as a .txt or .csv and upload it to R Studio.

```
qPCR_raw_data <- read.delim("2018.05.30_wdfy3_test2_setup.txt", head=T)
str(qPCR_raw_data) # Confirm that "SAMPLE" and "PRIMER" are factors, and "Ct" is numeric.

# If necessary, use these commands to change your data to the needed format:
qPCR_raw_data$SAMPLE <- as.factor(qPCR_raw_data$SAMPLE)
qPCR_raw_data$PRIMER <- as.factor(qPCR_raw_data$PRIMER)
qPCR_raw_data$Ct <- as.numeric(as.character(qPCR_raw_data$Ct))

```

We can remove the "NTC" and "NRT" samples, since their use is as a general control and will not be taken into account in the formal analysis of the data.

```
qPCR_raw_data <- qPCR_raw_data[qPCR_raw_data$SAMPLE!="NTC" & qPCR_raw_data$PRIMER!="NRT",]
```

Now, we can process the analysis step by step:

1. **Obtain average Ct values per primer** 

We have to obtain the average (mean) and standard deviation (sd) per primer per sample:

```
avg_Ct <- aggregate(Ct ~ PRIMER + SAMPLE, data=qPCR_raw_data, FUN="mean")
sd_Ct <- aggregate(Ct ~ PRIMER + SAMPLE, data=qPCR_raw_data, FUN="sd")

```

We can create a plot to quickly observe the behavior of the Ct values across our samples:

```
library(ggplot2)  
ggplot(Ct_estimates,aes(x=SAMPLE,y=avg_Ct, fill=PRIMER))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(aes(ymin=avg_Ct-sd_Ct,ymax=avg_Ct+sd_Ct), width=.2, position= position_dodge(.9))+xlab("Sample")+ylab("Ct")+theme_minimal()

```

![Ct_plot](https://github.com/jmuribes/images/blob/master/Ct_plot.png)


2. **Comparisson to a housekeeping gene**

First, by looking at the previous graph, you can decide which housekeeping gene to use in your further analysis. In this case, I will use "beta-actin".

We can split our original raw qPCR data into subsets of data per **sample**:
```
qPCR_data_split <- split(qPCR_raw_data,qPCR_raw_data$SAMPLE,drop=T)
```
Now, we can generate a for-loop to substract to each Ct value in each subset, the average of the desired housekeeping gene:

```
dCt_results <- NULL
for(i in 1:length(qPCR_data_split)){
  temporal_data <- x[[i]]
  temporal_data <- temporal_data$Ct - Ct_estimates$avg_Ct[Ct_estimates$SAMPLE==temporal_data[1,1]&Ct_estimates$PRIMER=="beta-actin"]
  dCt_results<-c(dCt_results,temporal_data)
}
qPCR_data_ordered <- qPCR_raw_data[order(qPCR_raw_data$SAMPLE),]
qPCR_data_ordered$dCt <- dCt_results
```
At this point, we should have a new dataset called "qPCR_data_ordered" that includes a new column called "dCt". Now, we can estimate the average (mean) and standard deviation (sd) of the dCt as calculated in the previous step.

```
avg_dCt <- aggregate(dCt ~ PRIMER + SAMPLE, data=qPCR_data_ordered, FUN="mean")
sd_dCt <- aggregate(dCt ~ PRIMER + SAMPLE, data=qPCR_data_ordered, FUN="sd")
dCt_estimates <- merge(avg_dCt,sd_dCt,by=c("SAMPLE","PRIMER"))
colnames(dCt_estimates) <-  c("SAMPLE","PRIMER","avg_dCt","sd_dCt")
```

We can also plot this information:

```
ggplot(dCt_estimates,aes(x=SAMPLE,y=avg_dCt, fill=PRIMER))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(aes(ymin=avg_dCt-sd_dCt,ymax=avg_dCt+sd_dCt), width=.2, position= position_dodge(.9))+xlab("Sample")+ylab("dCt")+theme_minimal()
```

![dCt_plot](https://github.com/jmuribes/images/blob/master/dCt_plot.png)


3. **Comparisson to one specific sample**

In most cases, we want to determine if there is a significant difference between sample types (i.ex. heterozyous/homozygous against the wild-type).

In this case, I will use as reference the sample named "wt-1_A6".

One strategy would be to generate a new column that includes the information from the reference sample repeatedly, so we can further substract it to the other samples.

In this specific example, we have the samples in **triplicates** and the housekeeping genes in **duplicates**. So, we have to generate the repetitions according to this information.

First, let's be sure the data is ordered, so it further fits after the split.

```
qPCR_data_ordered_ref <- qPCR_data_ordered[order(qPCR_data_ordered$SAMPLE,qPCR_data_ordered$PRIMER),]
head(qPCR_data_ordered_ref)
```
![triplicates_duplicates](https://github.com/jmuribes/images/blob/master/triplicates_duplicates.png)

Second, let's create a vector that has the information from the chosen reference sample (i.ex. wt-2_A6), where the dCt values are triplicated for the sample primers and duplicated for the housekeeping primers.

```
ref_sample <- dCt_estimates[dCt_estimates$SAMPLE=="wt-1_A6",] #isolate the information from the reference sample
ref_sample$repeats <- c(2,2,3,3) #a new column "repeats" determines the amount of times that row will be repeated.

```
We can confirm the number of samples in our analysis before going forward:

```
levels(qPCR_data_ordered$SAMPLE)
```

In this example, we have 7 samples: "NTC", "het-1_A1", "het-2_A12", "hom-del-1_A9", "hom-del-2_A10", "wt-1_A6", "wt-2_A7".

So, we create a vector that repeats the dCt values from the chosen reference to substract it to each sample. Then, we add it to the general ordered data that has the dCt values:

```
ref_repeated_vector <- rep(ref_sample_expanded$avg_dCt,6) #6 because we only need 7 repeats in this case.
qPCR_data_ordered_ref$ref_dCt <- ref_repeated_vector
```

Now, to each dCt we have to substract the reference dCt to calculate (ddCt).

```
qPCR_data_ordered_ref$ddCt <- qPCR_data_ordered_ref$dCt - qPCR_data_ordered_ref$ref_dCt

```

Once the ddCt has been calculated, we can calculate the **fold change**, which is given by the formula: **2^-ddCt**.

```
qPCR_data_ordered_ref$fold <- 2^(-qPCR_data_ordered_ref$ddCt)

```

Then, we can obtain the average (mean) and standard deviation (sd) of this fold change as previously done:

```
avg_fold <- aggregate(fold ~ PRIMER + SAMPLE, data=qPCR_data_ordered_ref, FUN="mean")
sd_fold <- aggregate(fold ~ PRIMER + SAMPLE, data=qPCR_data_ordered_ref, FUN="sd")
fold_estimates <- merge(avg_fold,sd_fold,by=c("SAMPLE","PRIMER"))
colnames(fold_estimates) <-  c("SAMPLE","PRIMER","avg_fold","sd_fold")

```

Finally, we can generate a quick plot to observe the results:

```
ggplot(fold_estimates,aes(x=SAMPLE,y=avg_fold, fill=PRIMER))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(aes(ymin=avg_fold-sd_fold,ymax=avg_fold+sd_fold), width=.2, position= position_dodge(.9))+xlab("Sample")+ylab("Fold change with respect to wt-1_A6")+theme_minimal()

## Note: in "ylab" remember to put the name of the sample you chose as reference!

```

![fold_plot](https://github.com/jmuribes/images/blob/master/fold_plot.png)
